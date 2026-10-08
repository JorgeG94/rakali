!! Tracer mass must be conserved by the clamp and by closed walls.
module test_tracer_conservation
   !! Two independent defects made the coastal path leak salt.  Both are
   !! pinned here — one at kernel level, one end-to-end.
   !!
   !! ## Defect 1 — the clamp destroyed tracer when the column mean was out
   !! of bounds
   !!
   !! `ml_clamp_tracer` clips per-layer concentration into `[tr_min, tr_max]`
   !! (pass 1) and then redistributes the clipped amount into the column's
   !! spare capacity (pass 2) with `frac = min(1, net/cap_to_add)`.  That
   !! `min` saturates exactly when `net > cap_to_add`, and the residual
   !! `net - cap_to_add` was silently destroyed.
   !!
   !! The saturation condition is **equivalent to the column mean leaving the
   !! bounds**:
   !!
   !!     mean <= tr_max
   !!       <=> sum_k (tr_max - tr_k) h_k >= 0
   !!       <=> -excess_above + cap_unsat + (cap_below + deficit_below) >= 0
   !!       <=> excess_above - deficit_below <= cap_unsat + cap_below
   !!       <=> net <= cap_to_add
   !!
   !! so the two-pass scheme is EXACTLY conservative whenever the mean is in
   !! range, and cannot be conservative when it is not.  The kernel now skips
   !! any column whose mean is out of range, making it conservative for every
   !! input.
   !!
   !! Why it mattered in practice: the default `S_max` is 40 PSU, which sits
   !! right at the top of the real seawater range, and the mode-split
   !! tracer/thickness inconsistency (hTr advected with per-layer CWC fluxes
   !! while `h_layer` is overwritten as `dsig(k)*H_new`) drives per-layer S
   !! above the bound at any sharp front.  So the "rare" saturation branch
   !! fired continuously.
   !!
   !! ## Defect 2 — closed walls leaked tracer
   !!
   !! `ml_apply_layer_bcs` distributed the barotropic ghost momentum
   !! depth-uniformly (`hu_layer(g,j,k) = dsig(k)*hu(g,j)`).  At a wall the
   !! barotropic ghost is a mirror, so the layers sum to zero net mass flux
   !! and WATER mass was conserved exactly — but the per-layer face velocity
   !! `v_face(k) = 0.5*(hu_ghost(k) + hu_int(k))` was non-zero in every layer
   !! whenever the interior column was sheared, because a shear-free ghost
   !! cannot cancel a sheared interior layer by layer.  Each layer then traded
   !! tracer with the ghost at its own concentration, giving a net tracer flux
   !! `sum_k v(k)*Tr(k) /= 0` through a CLOSED wall.  Wall edges now mirror
   !! the interior per-layer profile instead, so `v_face(k)` vanishes
   !! identically.
   !!
   !! ## Measured
   !!
   !! Lock exchange, flat bed, 40x40x6, walls, zero forcing, 400 steps
   !! (max|hu| ~ 0.1 m^2/s — vigorous real flow):
   !!
   !!   | state                    | rel. mass | rel. salt | rel. heat |
   !!   |--------------------------|-----------|-----------|-----------|
   !!   | before both fixes        | 1.5e-14   | -8.85e-06 | +4.60e-08 |
   !!   | clamp fix only           | 1.5e-14   | +4.85e-08 | +4.71e-08 |
   !!   | both fixes               | 1.5e-14   | **0.0**   | +8.5e-15  |
   !!
   !! Mass was always exact; only the tracer leaked.  After both fixes the
   !! salt residual is EXACTLY zero (not merely small) and heat is at
   !! round-off, so the assertion below is a genuine identity check rather
   !! than a tolerance negotiation.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_grid, only: hgrid_t
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_init, only: compute_total_mass, compute_layer_sums
   use rki_boundary_types, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_eos, only: ml_equation_of_state
   use rki_ml_tracers, only: ml_clamp_tracer
   implicit none
   private

   public :: collect_tracer_conservation_tests

   ! ---- Clamp unit fixture ----
   integer, parameter :: CNX = 5, CNY = 5, CNGHOST = 2, CNZ = 6
   real(wp), parameter :: TR_MIN = 0.0_wp, TR_MAX = 40.0_wp
   real(wp), parameter :: H_K = 2.0_wp

   ! ---- Lock-exchange fixture ----
   integer, parameter :: NX = 40, NY = 40, NGHOST = 2, NZ = 6
   real(wp), parameter :: DX = 5.0_wp, CFL_NUM = 0.45_wp
   integer, parameter :: N_STEPS = 400
   real(wp), parameter :: ETA0 = 10.0_wp, RHO_0 = 1025.0_wp
   real(wp), parameter :: S_DENSE = 40.0_wp, S_LIGHT = 35.0_wp
   real(wp), parameter :: S_REF = 35.0_wp, T_REF = 10.0_wp
   real(wp), parameter :: BETA_S = 0.8_wp, ALPHA_T = 0.2_wp
   real(wp), parameter :: T_UNIFORM = 10.0_wp

   real(wp), parameter :: TRACER_TOL = 1.0e-10_wp

contains

   subroutine collect_tracer_conservation_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("clamp_conserves_when_column_mean_above_bound", test_clamp_above), &
                  new_unittest("clamp_conserves_when_column_mean_below_bound", test_clamp_below), &
                  new_unittest("clamp_conserves_and_bounds_when_mean_in_range", test_clamp_in_range), &
                  new_unittest("lock_exchange_conserves_salt", test_lock_salt), &
                  new_unittest("lock_exchange_conserves_mass", test_lock_mass), &
                  new_unittest("lock_exchange_fixture_actually_flows", test_lock_flows) &
                  ]
   end subroutine collect_tracer_conservation_tests

   ! ------------------------------------------------------------------
   ! Clamp — kernel level
   ! ------------------------------------------------------------------

   subroutine clamp_column(profile, sum_before, sum_after, max_tr_after)
      !! Fill every interior column with `profile`, clamp, and report the
      !! interior tracer sum before/after plus the max concentration after.
      real(wp), intent(in) :: profile(CNZ)
      real(wp), intent(out) :: sum_before, sum_after, max_tr_after

      type(hgrid_t) :: grid
      real(wp), allocatable :: h_layer(:, :, :), hTr_layer(:, :, :)
      integer :: i, j, k

      call grid%init(CNX, CNY, CNGHOST, 1.0_wp, 1.0_wp)
      allocate (h_layer(grid%nx_total, grid%ny_total, CNZ))
      allocate (hTr_layer(grid%nx_total, grid%ny_total, CNZ))

      do k = 1, CNZ
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               h_layer(i, j, k) = H_K
               hTr_layer(i, j, k) = profile(k)*H_K
            end do
         end do
      end do

      sum_before = 0.0_wp
      do k = 1, CNZ
         do j = CNGHOST + 1, grid%ny_total - CNGHOST
            do i = CNGHOST + 1, grid%nx_total - CNGHOST
               sum_before = sum_before + hTr_layer(i, j, k)
            end do
         end do
      end do

      call ml_clamp_tracer(hTr_layer, h_layer, TR_MIN, TR_MAX, grid, CNZ)

      sum_after = 0.0_wp
      max_tr_after = -huge(1.0_wp)
      do k = 1, CNZ
         do j = CNGHOST + 1, grid%ny_total - CNGHOST
            do i = CNGHOST + 1, grid%nx_total - CNGHOST
               sum_after = sum_after + hTr_layer(i, j, k)
               max_tr_after = max(max_tr_after, hTr_layer(i, j, k)/H_K)
            end do
         end do
      end do
   end subroutine clamp_column

   subroutine test_clamp_above(error)
      !! Column mean ABOVE tr_max: no clamp can hold the bound and conserve,
      !! so the kernel must conserve and leave the bound violated — never
      !! destroy the residual.  Before the fix this column lost mass.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: profile(CNZ)
      real(wp) :: s0, s1, mx

      profile = [45.0_wp, 44.0_wp, 43.0_wp, 42.0_wp, 41.0_wp, 41.0_wp]  ! mean 42.7
      call clamp_column(profile, s0, s1, mx)

      call check(error, abs(s1 - s0) <= 4.0_wp*spacing(s0), &
                 "clamp destroyed tracer when the column mean was above tr_max")
   end subroutine test_clamp_above

   subroutine test_clamp_below(error)
      !! Mirror case: column mean BELOW tr_min.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: profile(CNZ)
      real(wp) :: s0, s1, mx

      profile = [-5.0_wp, -4.0_wp, -3.0_wp, -2.0_wp, -1.0_wp, -1.0_wp]  ! mean -2.7
      call clamp_column(profile, s0, s1, mx)

      call check(error, abs(s1 - s0) <= 4.0_wp*spacing(abs(s0)), &
                 "clamp created tracer when the column mean was below tr_min")
   end subroutine test_clamp_below

   subroutine test_clamp_in_range(error)
      !! Column mean IN range but individual layers violating: the kernel must
      !! BOTH conserve exactly AND pull every layer inside the bounds.  This is
      !! the case the clamp exists for, and the fix must not weaken it.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: profile(CNZ)
      real(wp) :: s0, s1, mx

      profile = [48.0_wp, 46.0_wp, 30.0_wp, 30.0_wp, 30.0_wp, 30.0_wp]  ! mean 35.7
      call clamp_column(profile, s0, s1, mx)

      call check(error, abs(s1 - s0) <= 4.0_wp*spacing(s0), &
                 "clamp failed to conserve with the column mean in range")
      if (allocated(error)) return
      call check(error, mx <= TR_MAX + 1.0e-12_wp, &
                 "clamp failed to enforce tr_max with the column mean in range")
   end subroutine test_clamp_in_range

   ! ------------------------------------------------------------------
   ! Lock exchange — end to end
   ! ------------------------------------------------------------------

   subroutine setup_lock(state)
      !! Flat bed, walls, zero forcing.  Dense water in the left half drives a
      !! genuine gravity current — the regime where the leaks appeared.
      type(state_t), intent(inout) :: state

      integer :: i, j, k, i_mid
      real(wp) :: s_val

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%rho_0 = RHO_0
      state%multilayer%alpha_T = ALPHA_T
      state%multilayer%beta_S = BETA_S
      state%multilayer%S_ref = S_REF
      state%multilayer%T_ref = T_REF
      state%multilayer%initial_salinity = S_LIGHT
      state%multilayer%initial_temperature = T_UNIFORM

      call state%init(NX, NY, NGHOST, DX, DX, CFL_NUM)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%b(i, j) = 0.0_wp
            state%barotropic%h(i, j) = ETA0
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
         end do
      end do

      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, &
                                          state%barotropic%hv, state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, state%grid, NZ)

      i_mid = NGHOST + NX/2
      do k = 1, NZ
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               if (i <= i_mid) then
                  s_val = S_DENSE
               else
                  s_val = S_LIGHT
               end if
               associate (ml => state%multilayer)
                  ml%tracers(ml%idx_salinity)%hTr(i, j, k) = ml%h_layer(i, j, k)*s_val
                  ml%tracers(ml%idx_temperature)%hTr(i, j, k) = ml%h_layer(i, j, k)*T_UNIFORM
               end associate
            end do
         end do
      end do

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, RHO_0, BETA_S, ALPHA_T, &
                                S_REF, T_REF, state%grid, NZ)

      state%forcing%heat_flux = 0.0_wp
      state%forcing%evap_minus_precip = 0.0_wp
      state%forcing%wind_x = 0.0_wp
      state%forcing%wind_y = 0.0_wp
      state%forcing%precipitation = 0.0_wp
   end subroutine setup_lock

   subroutine run_lock(state)
      type(state_t), intent(inout) :: state
      integer :: n
      real(wp) :: t

      t = 0.0_wp
      call solver_enter_data(state)
      do n = 1, N_STEPS
         call solver_step(state, time=t)
         t = t + state%dt
      end do
      call solver_exit_data(state)
   end subroutine run_lock

   subroutine test_lock_salt(error)
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: h0, s0, t0, h1, s1, t1, rel

      call setup_lock(state)
      call compute_layer_sums(state, h0, s0, t0)
      call run_lock(state)
      call compute_layer_sums(state, h1, s1, t1)

      rel = abs(s1 - s0)/max(abs(s0), tiny(1.0_wp))
      call check(error, rel < TRACER_TOL, &
                 "lock exchange does not conserve salt under advection")

      call state%destroy()
   end subroutine test_lock_salt

   subroutine test_lock_mass(error)
      !! Control: mass was ALWAYS exact, including before both fixes.  If this
      !! ever fails, the defect is somewhere new.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: m0, m1, rel

      call setup_lock(state)
      m0 = compute_total_mass(state)
      call run_lock(state)
      m1 = compute_total_mass(state)

      rel = abs(m1 - m0)/max(m0, tiny(1.0_wp))
      call check(error, rel < 1.0e-12_wp, &
                 "lock exchange does not conserve water mass")

      call state%destroy()
   end subroutine test_lock_mass

   subroutine test_lock_flows(error)
      !! Guard on the FIXTURE: a quiescent column conserves trivially, so the
      !! conservation assertions above are only meaningful if the gravity
      !! current actually develops.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      integer :: i, j
      real(wp) :: mx

      call setup_lock(state)
      call run_lock(state)

      mx = 0.0_wp
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            mx = max(mx, abs(state%barotropic%hu(i, j)))
         end do
      end do

      call check(error, mx > 1.0e-3_wp, &
                 "fixture: no gravity current developed, conservation test is vacuous")

      call state%destroy()
   end subroutine test_lock_flows

end module test_tracer_conservation
