!! Evaporation as a surface water-mass sink (coastal multilayer, structured).
!!
!! Exercises the `evap_mass_sink` path in `ml_apply_surface_tracer_fluxes`:
!! net evaporation removes REAL water from the column (column-proportional
!! scaling, salt held) rather than injecting a virtual salt flux.  The
!! physical contract this test pins:
!!
!!   (a) DEPTH DRAWDOWN — a closed basin under a spatially-uniform, steady
!!       evaporation rate E loses water linearly: h(t) = h0 - E*t.
!!   (b) SALT MASS CONSTANT — evaporation removes WATER not salt.  The salt
!!       mass hS (PSU*m) is left untouched, so the column total salt mass
!!       is conserved to round-off while the CONCENTRATION S = hS/h rises
!!       to the analytic h0*S0/h(t).
!!   (c) TEMPERATURE CONCENTRATION CONSTANT — the evaporated water carries
!!       its own temperature away (hT scales with h), so the temperature
!!       concentration T = hT/h stays constant: no spurious warming.  (v1
!!       does NOT model latent-heat cooling — documented follow-up.)
!!   (d) THIN-FLAT DRYING — a shallow flat under evaporation must dry
!!       gracefully: no negative depth, no NaN.  The sink clamps to the
!!       water available above ML_MIN_DEPTH and the unsatisfied fraction
!!       is dropped honestly (the water was not there to evaporate; no
!!       compensating source is invented).  Salt accounting at the dry
!!       margin: hS is held while h floors, so the concentration spikes
!!       and the S_max clamp then destroys the excess — the irreducible
!!       cost of a hard bound, stated openly.
!!
!! Bit-identity of the default path (knob off / E=0) is covered by the
!! existing ML regression suite (surface_flux_strat, lock_exchange,
!! river_plume) — this suite only exercises the new mass-sink path.
module test_coastal_evap
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use rki_constants, only: wp, CP_WATER, RHO_WATER, LATENT_HEAT_VAPORIZATION
   use rki_ml_common, only: ml_set_min_depth, ML_MIN_DEPTH
   use rki_state, only: state_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_eos, only: ml_equation_of_state
   use rki_ml_tracers, only: ml_apply_surface_tracer_fluxes
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_coastal_evap_tests

   ! Domain — small closed basin
   integer, parameter :: NX_PHYS = 6
   integer, parameter :: NY_PHYS = 6
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4

   real(wp), parameter :: DX_VAL = 100.0_wp

   ! Deep basin so evaporation over the test window stays well above the
   ! dry floor and the drawdown is cleanly linear.
   real(wp), parameter :: H0 = 20.0_wp
   real(wp), parameter :: S_INIT = 35.0_wp
   real(wp), parameter :: T_INIT = 15.0_wp

   ! EOS — referenced at the IC so the initial density field is uniform;
   ! the only thing that should move the fields is the evaporation sink.
   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: BETA_S = 0.78_wp
   real(wp), parameter :: ALPHA_T = 0.17_wp
   real(wp), parameter :: T_REF = 15.0_wp
   real(wp), parameter :: S_REF = 35.0_wp

   real(wp), parameter :: ML_MIN_DEPTH_DEFAULT = 0.5_wp
      !! Restored after the thin-flat test (module-global protected var).

contains

   subroutine collect_coastal_evap_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("depth_drawdown_linear", test_depth_drawdown), &
                  new_unittest("salt_mass_constant_concentration_rises", &
                               test_salt_mass_constant), &
                  new_unittest("temperature_concentration_constant", &
                               test_temperature_concentration), &
                  new_unittest("thin_flat_dries_gracefully", &
                               test_thin_flat_dries), &
                  new_unittest("fused_nz1_mass_sink", test_fused_nz1), &
                  new_unittest("latent_cooling_kernel_exact", &
                               test_latent_cooling_kernel), &
                  new_unittest("latent_cooling_solver_cools_column", &
                               test_latent_cooling_solver) &
                  ]
   end subroutine collect_coastal_evap_tests

   subroutine setup_evap_basin(state, h0_val, evap_rate)
      !! Closed basin, uniform S/T, uniform prescribed evaporation rate E
      !! (m/s), evaporation-as-mass-sink knob ON.  A generous S_max keeps
      !! the conservative clamp from eating the concentration rise.
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: h0_val, evap_rate

      integer :: i, j, k

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_salinity = S_INIT
      state%multilayer%initial_temperature = T_INIT
      state%multilayer%alpha_T = ALPHA_T
      state%multilayer%beta_S = BETA_S
      state%multilayer%rho_0 = RHO_0
      state%multilayer%T_ref = T_REF
      state%multilayer%S_ref = S_REF
      ! Wide salinity bound so evaporative concentration is not clamped.
      state%multilayer%S_max = 500.0_wp
      state%multilayer%S_min = 0.0_wp
      ! The new mass-sink path (default off elsewhere ⇒ bit-identical).
      state%multilayer%evap_mass_sink = .true.

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 1.0_wp)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      state%barotropic%h = h0_val
      state%barotropic%b = 0.0_wp
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp

      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, &
                                          state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      do k = 1, NZ
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = &
                  state%multilayer%h_layer(i, j, k)*S_INIT
               state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = &
                  state%multilayer%h_layer(i, j, k)*T_INIT
            end do
         end do
      end do

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, &
                                RHO_0, BETA_S, S_REF, ALPHA_T, T_REF, &
                                state%grid, NZ)

      ! Uniform net evaporation (positive) everywhere; no precip, no heat,
      ! no wind.
      state%forcing%heat_flux = 0.0_wp
      state%forcing%evap_minus_precip = evap_rate
      state%forcing%wind_x = 0.0_wp
      state%forcing%wind_y = 0.0_wp
      state%forcing%precipitation = 0.0_wp
   end subroutine setup_evap_basin

   subroutine run_to(state, t_final, t_actual)
      !! Advance to at least t_final; return the exact accumulated
      !! simulated time (dt is CFL-adaptive, so callers must use this for
      !! analytic comparisons).
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: t_final
      real(wp), intent(out) :: t_actual

      call solver_enter_data(state)
      t_actual = 0.0_wp
      do while (t_actual < t_final)
         call solver_step(state, time=t_actual)
         t_actual = t_actual + state%dt
      end do
      call solver_exit_data(state)
   end subroutine run_to

   pure function column_depth(state, i, j) result(h_col)
      type(state_t), intent(in) :: state
      integer, intent(in) :: i, j
      real(wp) :: h_col
      integer :: k
      h_col = 0.0_wp
      do k = 1, NZ
         h_col = h_col + state%multilayer%h_layer(i, j, k)
      end do
   end function column_depth

   pure function column_salt_mass(state, i, j) result(hs_col)
      !! Column-integrated salt mass (PSU*m).
      type(state_t), intent(in) :: state
      integer, intent(in) :: i, j
      real(wp) :: hs_col
      integer :: k
      hs_col = 0.0_wp
      do k = 1, NZ
         hs_col = hs_col + state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k)
      end do
   end function column_salt_mass

   pure function domain_salt_mass(state) result(hs_tot)
      !! Domain-total salt mass over the physical interior (PSU*m*cells).
      type(state_t), intent(in) :: state
      real(wp) :: hs_tot
      integer :: i, j
      hs_tot = 0.0_wp
      do j = NGHOST + 1, NGHOST + NY_PHYS
         do i = NGHOST + 1, NGHOST + NX_PHYS
            hs_tot = hs_tot + column_salt_mass(state, i, j)
         end do
      end do
   end function domain_salt_mass

   pure function domain_depth(state) result(h_tot)
      !! Domain-total water column depth over the physical interior.
      type(state_t), intent(in) :: state
      real(wp) :: h_tot
      integer :: i, j
      h_tot = 0.0_wp
      do j = NGHOST + 1, NGHOST + NY_PHYS
         do i = NGHOST + 1, NGHOST + NX_PHYS
            h_tot = h_tot + column_depth(state, i, j)
         end do
      end do
   end function domain_depth

   pure function column_heat_mass(state, i, j) result(ht_col)
      !! Column-integrated heat content (degC*m).
      type(state_t), intent(in) :: state
      integer, intent(in) :: i, j
      real(wp) :: ht_col
      integer :: k
      ht_col = 0.0_wp
      do k = 1, NZ
         ht_col = ht_col + state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k)
      end do
   end function column_heat_mass

   subroutine test_depth_drawdown(error)
      !! Uniform evaporation E ⇒ h(t) = h0 - E*t.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: E_RATE = 1.0e-3_wp   !! m/s (strong, for SNR)
      real(wp), parameter :: T_END = 100.0_wp     !! s
      integer :: ic, jc
      real(wp) :: h_expected, h_final, t

      call setup_evap_basin(state, H0, E_RATE)
      ic = NGHOST + NX_PHYS/2
      jc = NGHOST + NY_PHYS/2

      call run_to(state, T_END, t)

      h_expected = H0 - E_RATE*t

      ! Domain-mean drawdown is exact to round-off (the sink is applied
      ! identically everywhere; wall-adjacent columns exchange O(1e-8) of
      ! water with their neighbours but the domain total is conserved).
      h_final = domain_depth(state)/real(NX_PHYS*NY_PHYS, wp)
      call check(error, abs(h_final - h_expected) < 1.0e-10_wp*H0, &
                 "domain-mean evaporative drawdown must match h0 - E*t")
      if (allocated(error)) then
         call state%destroy(); return
      end if

      ! Per-column drawdown: same analytic value, loose enough to admit
      ! the round-off-scale horizontal redistribution.
      h_final = column_depth(state, ic, jc)
      call check(error, abs(h_final - h_expected) < 1.0e-7_wp*H0, &
                 "per-column evaporative drawdown must match h0 - E*t")

      call state%destroy()
   end subroutine test_depth_drawdown

   subroutine test_salt_mass_constant(error)
      !! Salt mass is conserved to round-off; concentration rises to the
      !! analytic h0*S0/h(t).
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: E_RATE = 1.0e-3_wp
      real(wp), parameter :: T_END = 100.0_wp
      integer :: ic, jc
      real(wp) :: hs_initial, hs_final
      real(wp) :: s_mean_final, s_expected, t

      call setup_evap_basin(state, H0, E_RATE)
      ic = NGHOST + NX_PHYS/2
      jc = NGHOST + NY_PHYS/2

      hs_initial = domain_salt_mass(state)

      call run_to(state, T_END, t)

      hs_final = domain_salt_mass(state)

      ! (b1) TOTAL salt mass exactly constant to round-off (evaporation
      ! removes water, not salt).  Domain total is the physically
      ! meaningful conservation statement: individual columns exchange
      ! O(1e-9) relative salt with neighbours through round-off-scale
      ! currents, but the domain integral must hold exactly.
      call check(error, abs(hs_final - hs_initial) < 1.0e-12_wp*hs_initial, &
                 "domain-total salt mass must be conserved under evaporation")
      if (allocated(error)) then
         call state%destroy(); return
      end if

      s_mean_final = column_salt_mass(state, ic, jc)/column_depth(state, ic, jc)

      ! (b2) Column-mean concentration must have RISEN above the initial S.
      call check(error, s_mean_final > S_INIT, &
                 "evaporation must raise the column-mean salinity concentration")
      if (allocated(error)) then
         call state%destroy(); return
      end if

      ! (b3) Column-mean salinity follows h0*S0/h(t): salt mass fixed,
      ! depth drops linearly ⇒ S_mean(t) = hS0/h(t) = H0*S0/(H0 - E*t).
      ! (Concentration is robust to the redistribution: water and salt
      ! move together at uniform S.  Measured error ~1e-9 relative.)
      s_expected = (H0*S_INIT)/(H0 - E_RATE*t)
      call check(error, abs(s_mean_final - s_expected) < 1.0e-7_wp*s_expected, &
                 "column-mean salinity must follow the analytic concentration rise")

      call state%destroy()
   end subroutine test_salt_mass_constant

   subroutine test_temperature_concentration(error)
      !! Heat leaves with the evaporated water ⇒ the temperature
      !! CONCENTRATION is unchanged (no spurious warming/cooling in v1).
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: E_RATE = 1.0e-3_wp
      real(wp), parameter :: T_END = 100.0_wp
      integer :: ic, jc
      real(wp) :: t_mean_final, h_final, ht_final, t

      call setup_evap_basin(state, H0, E_RATE)
      ic = NGHOST + NX_PHYS/2
      jc = NGHOST + NY_PHYS/2

      call run_to(state, T_END, t)

      ht_final = column_heat_mass(state, ic, jc)
      h_final = column_depth(state, ic, jc)
      t_mean_final = ht_final/h_final

      ! Column-mean temperature concentration must be UNCHANGED: the
      ! evaporated water leaves at its own temperature (hT scales with h),
      ! so removing water does not heat or cool the remaining column.
      ! (v1 does NOT apply latent-heat cooling — documented follow-up.)
      call check(error, abs(t_mean_final - T_INIT) < 1.0e-9_wp*T_INIT, &
                 "column-mean temperature concentration must be preserved")

      call state%destroy()
   end subroutine test_temperature_concentration

   subroutine test_thin_flat_dries(error)
      !! A thin flat under evaporation must dry without negative h or NaN.
      !! ML_MIN_DEPTH is lowered so the sink actually engages on the thin
      !! column (the default 0.5 m floor would gate a 0.05 m flat off) and
      !! restored afterwards.  The sink clamps to the water available above
      !! the floor; the unsatisfied fraction is dropped honestly.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: H_THIN = 0.05_wp    !! shallow flat (m)
      real(wp), parameter :: E_RATE = 1.0e-3_wp
      real(wp), parameter :: T_END = 400.0_wp    !! long enough to over-draw
      real(wp), parameter :: FLOOR_THIN = 1.0e-3_wp
      integer :: i, j, k, ic, jc
      real(wp) :: h_min, h_val, hs_val, t, h_final
      logical :: any_nan

      call ml_set_min_depth(FLOOR_THIN)

      call setup_evap_basin(state, H_THIN, E_RATE)
      ic = NGHOST + NX_PHYS/2
      jc = NGHOST + NY_PHYS/2

      call run_to(state, T_END, t)

      ! No layer thickness may be negative; no NaN anywhere in h or hS.
      h_min = huge(1.0_wp)
      any_nan = .false.
      do k = 1, NZ
         do j = NGHOST + 1, NGHOST + NY_PHYS
            do i = NGHOST + 1, NGHOST + NX_PHYS
               h_val = state%multilayer%h_layer(i, j, k)
               hs_val = state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k)
               h_min = min(h_min, h_val)
               if (ieee_is_nan(h_val) .or. ieee_is_nan(hs_val)) any_nan = .true.
            end do
         end do
      end do

      call check(error,.not. any_nan, "thin-flat evaporation must not produce NaN")
      if (allocated(error)) then
         call ml_set_min_depth(ML_MIN_DEPTH_DEFAULT)
         call state%destroy(); return
      end if

      call check(error, h_min >= 0.0_wp, &
                 "thin-flat evaporation must not drive layer depth negative")
      if (allocated(error)) then
         call ml_set_min_depth(ML_MIN_DEPTH_DEFAULT)
         call state%destroy(); return
      end if

      ! The demand (E*t = 0.4 m) exceeds the water available (0.05 m), so
      ! the column must have drawn down toward the floor and STOPPED there —
      ! the sink clamps rather than over-drafting.
      h_final = column_depth(state, ic, jc)
      call check(error, h_final < H_THIN .and. h_final >= 0.5_wp*FLOOR_THIN, &
                 "thin flat must draw down and hold at the dry floor")

      call ml_set_min_depth(ML_MIN_DEPTH_DEFAULT)
      call state%destroy()
   end subroutine test_thin_flat_dries

   subroutine test_fused_nz1(error)
      !! nz=1 exercises the FUSED tracer-post kernel (ml_fused_tracer_post)
      !! — the solver dispatches to the fused fast path at a single layer,
      !! so the mass-sink branch there needs its own coverage.  Same
      !! invariants as the nz=4 tests: analytic drawdown, salt mass
      !! constant (concentration rises), temperature invariant.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: E_RATE = 1.0e-3_wp
      real(wp), parameter :: T_END = 100.0_wp
      integer :: i, j, ic, jc
      real(wp) :: t, hs_tot0, hs_tot, h_dom, h_expected
      real(wp) :: s_cell, t_cell

      state%use_multilayer = .true.
      state%multilayer%nz_ml = 1
      state%multilayer%initial_salinity = S_INIT
      state%multilayer%initial_temperature = T_INIT
      state%multilayer%alpha_T = ALPHA_T
      state%multilayer%beta_S = BETA_S
      state%multilayer%rho_0 = RHO_0
      state%multilayer%T_ref = T_REF
      state%multilayer%S_ref = S_REF
      state%multilayer%S_max = 500.0_wp
      state%multilayer%S_min = 0.0_wp
      state%multilayer%evap_mass_sink = .true.

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 1.0_wp)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      state%barotropic%h = H0
      state%barotropic%b = 0.0_wp
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp

      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, &
                                          state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, 1)

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, 1) = &
               state%multilayer%h_layer(i, j, 1)*S_INIT
            state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, 1) = &
               state%multilayer%h_layer(i, j, 1)*T_INIT
         end do
      end do

      state%forcing%heat_flux = 0.0_wp
      state%forcing%evap_minus_precip = E_RATE
      state%forcing%wind_x = 0.0_wp
      state%forcing%wind_y = 0.0_wp
      state%forcing%precipitation = 0.0_wp

      hs_tot0 = 0.0_wp
      do j = NGHOST + 1, NGHOST + NY_PHYS
         do i = NGHOST + 1, NGHOST + NX_PHYS
            hs_tot0 = hs_tot0 + state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, 1)
         end do
      end do

      call run_to(state, T_END, t)

      ic = NGHOST + NX_PHYS/2
      jc = NGHOST + NY_PHYS/2

      hs_tot = 0.0_wp
      h_dom = 0.0_wp
      do j = NGHOST + 1, NGHOST + NY_PHYS
         do i = NGHOST + 1, NGHOST + NX_PHYS
            hs_tot = hs_tot + state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, 1)
            h_dom = h_dom + state%multilayer%h_layer(i, j, 1)
         end do
      end do

      h_expected = H0 - E_RATE*t
      call check(error, abs(h_dom/real(NX_PHYS*NY_PHYS, wp) - h_expected) < 1.0e-10_wp*H0, &
                 "nz=1 fused: domain-mean drawdown must match h0 - E*t")
      if (allocated(error)) then
         call state%destroy(); return
      end if

      call check(error, abs(hs_tot - hs_tot0) < 1.0e-12_wp*hs_tot0, &
                 "nz=1 fused: domain-total salt mass must be conserved")
      if (allocated(error)) then
         call state%destroy(); return
      end if

      s_cell = state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(ic, jc, 1) &
               /state%multilayer%h_layer(ic, jc, 1)
      t_cell = state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ic, jc, 1) &
               /state%multilayer%h_layer(ic, jc, 1)

      call check(error, abs(s_cell - (H0*S_INIT)/(H0 - E_RATE*t)) < 1.0e-7_wp*S_INIT, &
                 "nz=1 fused: salinity concentration must follow h0*S0/h(t)")
      if (allocated(error)) then
         call state%destroy(); return
      end if

      call check(error, abs(t_cell - T_INIT) < 1.0e-9_wp*T_INIT, &
                 "nz=1 fused: temperature concentration must be preserved")

      call state%destroy()
   end subroutine test_fused_nz1

   subroutine test_latent_cooling_kernel(error)
      !! EXACT bookkeeping + energy-closure check on the latent-heat sink.
      !! A single call to `ml_apply_surface_tracer_fluxes` with
      !! evap_latent_cooling ON must:
      !!   (i)  scale every layer's hT by fac = 1 - dh_eff/h_col (sensible
      !!        heat leaves with the removed mass), AND
      !!   (ii) additionally cool the TOP layer (k=nz) by
      !!        latent_gain*dh_eff, latent_gain = rho_water*L_v/(rho_0*cp).
      !! Energy accounting: the latent heat withdrawn, converted back to
      !! J/m^2 (× rho_0*cp), must equal rho_water*L_v*dh_eff — the phase-
      !! change energy of the evaporated mass — closing to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: E_RATE = 1.0e-3_wp
      real(wp), parameter :: DT_K = 10.0_wp     !! fixed dt for the single call
      integer :: i, j, k, ic, jc
      real(wp), allocatable :: hT_pre(:, :, :), h_pre(:, :, :)
      real(wp) :: h_col, dh_eff, fac, latent_gain
      real(wp) :: hT_top_expected, hT_k_expected, latent_energy, evap_mass_energy

      call setup_evap_basin(state, H0, E_RATE)
      state%multilayer%evap_latent_cooling = .true.   ! v2 latent cooling ON
      ic = NGHOST + NX_PHYS/2
      jc = NGHOST + NY_PHYS/2

      call solver_enter_data(state)

      ! Registry lock WIRING (the guard itself is unit-tested in
      ! test_multilayer): once the solver has mapped tracers(:)
      ! element-by-element, growing the registry would strand every
      ! existing tracer's device mapping, so enter_data must leave the
      ! registry locked.  Asserted here rather than in the unit test
      ! because only a real enter_data exercises the wiring.
      call check(error, state%multilayer%registry_locked, &
                 "solver_enter_data must lock the tracer registry")
      if (allocated(error)) return

      ! Host snapshot of the pre-call layer fields (host == device right
      ! after enter_data's copyin).
      allocate (hT_pre, source=state%multilayer%tracers(state%multilayer%idx_temperature)%hTr)
      allocate (h_pre, source=state%multilayer%h_layer)

      ! One surface-flux pass at a fixed dt (heat_flux == 0 in the setup, so
      ! the only thermal change is the latent sink).
      call ml_apply_surface_tracer_fluxes( &
         state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
         state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
         state%multilayer%h_layer, &
         state%multilayer%hu_layer, state%multilayer%hv_layer, &
         state%forcing%heat_flux, state%forcing%evap_minus_precip, &
         state%multilayer%evap_mass_sink, state%multilayer%evap_latent_cooling, &
         state%multilayer%rho_0, CP_WATER, state%grid, NZ, DT_K)

      !$acc update self(state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
      !$acc&            state%multilayer%h_layer)
      call solver_exit_data(state)

      call check(error,.not. state%multilayer%registry_locked, &
                 "solver_exit_data must release the tracer-registry lock")
      if (allocated(error)) return

      h_col = 0.0_wp
      do k = 1, NZ
         h_col = h_col + h_pre(ic, jc, k)
      end do
      dh_eff = min(DT_K*E_RATE, max(h_col - ML_MIN_DEPTH, 0.0_wp))
      fac = 1.0_wp - dh_eff/h_col
      latent_gain = RHO_WATER*LATENT_HEAT_VAPORIZATION/(state%multilayer%rho_0*CP_WATER)

      ! (i) deeper layers: pure sensible scaling by fac.
      do k = 1, NZ - 1
         hT_k_expected = hT_pre(ic, jc, k)*fac
         call check(error, abs(state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ic, jc, k) &
                               - hT_k_expected) < 1.0e-10_wp*abs(hT_k_expected), &
                    "sub-surface hT must scale by fac (sensible heat leaves with water)")
         if (allocated(error)) then
            call state%destroy(); return
         end if
      end do

      ! (ii) top layer: sensible scaling PLUS latent withdrawal.
      hT_top_expected = hT_pre(ic, jc, NZ)*fac - latent_gain*dh_eff
      call check(error, abs(state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ic, jc, NZ) &
                            - hT_top_expected) < 1.0e-10_wp*abs(hT_top_expected), &
                 "top-layer hT must carry the additional latent-heat sink")
      if (allocated(error)) then
         call state%destroy(); return
      end if

      ! Energy accounting: measure the latent heat ACTUALLY removed from the
      ! post-kernel top-layer hT — the decrement beyond the sensible fac
      ! scaling, (hT_pre(nz)*fac - hT_post(nz)) — converted to J/m^2
      ! (x rho_0*cp).  It must equal the phase-change energy of the
      ! evaporated mass, rho_water * dh_eff * L_v.  Reading the real field
      ! makes this a genuine kernel check, not a host-side identity.
      latent_energy = (hT_pre(ic, jc, NZ)*fac &
                       - state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ic, jc, NZ)) &
                      *state%multilayer%rho_0*CP_WATER
      evap_mass_energy = RHO_WATER*dh_eff*LATENT_HEAT_VAPORIZATION
      call check(error, abs(latent_energy - evap_mass_energy) < 1.0e-9_wp*evap_mass_energy, &
                 "latent energy removed from post-kernel hT must equal rho_water*L_v*dh_eff (budget closes)")

      call state%destroy()
   end subroutine test_latent_cooling_kernel

   subroutine test_latent_cooling_solver(error)
      !! Integration check: with latent cooling wired through the full RK2
      !! solver, a closed basin under steady evaporation must COOL.  The
      !! column-mean temperature follows the continuum drawdown-cooling law
      !!   T_mean(t) = T_INIT - L*ln(H0/h(t)),   L = rho_water*L_v/(rho_0*cp),
      !! obtained by integrating dT/dt = -L*E/h(t) with h(t) = H0 - E*t.
      !! (Discrete RK2 introduces an O(E*dt/H0) quadrature error, hence a
      !! finite tolerance.)  A latent-OFF twin must NOT cool (guards against
      !! the knob leaking when disabled).
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: E_RATE = 1.0e-3_wp
      real(wp), parameter :: T_END = 100.0_wp
      integer :: ic, jc
      real(wp) :: t, h_final, t_mean, latent_L, drop_expected, drop_actual

      ic = NGHOST + NX_PHYS/2
      jc = NGHOST + NY_PHYS/2

      ! ---- latent-OFF twin: temperature concentration must stay put ----
      call setup_evap_basin(state, H0, E_RATE)
      state%multilayer%evap_latent_cooling = .false.
      call run_to(state, T_END, t)
      t_mean = column_heat_mass(state, ic, jc)/column_depth(state, ic, jc)
      call check(error, abs(t_mean - T_INIT) < 1.0e-9_wp*T_INIT, &
                 "latent OFF: temperature concentration must be unchanged")
      call state%destroy()
      if (allocated(error)) return

      ! ---- latent-ON: column cools along the continuum law ----
      call setup_evap_basin(state, H0, E_RATE)
      state%multilayer%evap_latent_cooling = .true.
      call run_to(state, T_END, t)

      h_final = column_depth(state, ic, jc)
      t_mean = column_heat_mass(state, ic, jc)/h_final

      latent_L = RHO_WATER*LATENT_HEAT_VAPORIZATION/(state%multilayer%rho_0*CP_WATER)
      drop_expected = latent_L*log(H0/(H0 - E_RATE*t))
      drop_actual = T_INIT - t_mean

      ! Must cool by a clearly non-trivial amount, and match the analytic
      ! drop to a few percent (RK2 quadrature error on the 1/h integral).
      call check(error, drop_actual > 0.5_wp, &
                 "latent ON: column must cool measurably")
      if (allocated(error)) then
         call state%destroy(); return
      end if
      call check(error, abs(drop_actual - drop_expected) < 3.0e-2_wp*drop_expected, &
                 "latent ON: cooling must follow T_INIT - L*ln(H0/h(t))")

      call state%destroy()
   end subroutine test_latent_cooling_solver

end module test_coastal_evap
