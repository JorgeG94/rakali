!! Sediment transport tests (structured): constant-velocity settling
!! conserves total (suspended + bed) sediment and sinks the suspended load
!! into the bed reservoir.  Runs a quiescent flat closed basin (no flow), so
!! the only process acting on the sediment tracer is settling — making this a
!! clean conservation + monotonicity check on `ml_settling_tracer` and the
!! once-per-step solver wiring.  Covers both the nz>1 path and the nz=1 fused
!! fast path.
module test_sediment
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, VCOORD_SIGMA
   use rki_state, only: state_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_eos, only: ml_equation_of_state
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_sediment_tests

   integer, parameter :: NX_PHYS = 8, NY_PHYS = 8, NGHOST = 2
   real(wp), parameter :: DX = 100.0_wp
   real(wp), parameter :: H0 = 20.0_wp
   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: C_INIT = 2.0_wp    !! initial suspended conc (kg/m^3)
   real(wp), parameter :: W_S = 0.02_wp      !! settling velocity (m/s)
   real(wp), parameter :: BED_INIT = 10.0_wp  !! initial bed mass (kg/m^2)
   real(wp), parameter :: TAU_CE = 0.1_wp    !! critical erosion stress (Pa)
   real(wp), parameter :: E_RATE = 1.0e-3_wp  !! erosion rate E0 (kg/m^2/s)

contains

   subroutine collect_sediment_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("sediment_settling_conserves_nz5", test_conserve_nz5), &
                  new_unittest("sediment_settling_conserves_nz1", test_conserve_nz1), &
                  new_unittest("sediment_ppm_multilayer_descent", test_ppm_multilayer_descent), &
                  new_unittest("sediment_erosion_fires", test_erosion_fires), &
                  new_unittest("sediment_erosion_quiescent_none", test_erosion_quiescent_none) &
                  ]
   end subroutine collect_sediment_tests

   subroutine setup_sediment_state(state, nz)
      !! Quiescent flat closed basin (no flow), nz sigma layers, sediment on.
      !! Uniform suspended C = C_INIT, empty bed.  At rest + no mixing so only
      !! settling moves the sediment tracer.
      type(state_t), intent(inout) :: state
      integer, intent(in) :: nz
      integer :: idx_sed

      state%use_multilayer = .true.
      state%multilayer%nz_ml = nz
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%rho_0 = RHO_0
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%kappa_S_bg = 0.0_wp
      state%multilayer%kappa_T_bg = 0.0_wp
      state%multilayer%hdiff_kappa = 0.0_wp
      state%multilayer%smagorinsky_cs = 0.0_wp
      state%multilayer%vcoord%coord_type = VCOORD_SIGMA
      ! Enable sediment BEFORE state%init so bed_sediment allocates.
      state%multilayer%use_sediment = .true.
      state%multilayer%sediment_w_s = W_S
      state%multilayer%sediment_bed_init = 0.0_wp

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX, DX, 0.45_wp)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      state%barotropic%coriolis_f = 0.0_wp
      state%barotropic%manning_n = 0.0_wp
      state%barotropic%b = -H0
      state%barotropic%h = H0
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp

      call ml_init_layers_from_barotropic(state%barotropic%h, &
                                          state%barotropic%hu, &
                                          state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, &
                                          state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, nz)

      ! state%init registers S/T only; register the sediment tracer here (the
      ! config path auto-registers it, but this test builds state directly).
      call state%multilayer%register_passive_tracer(state%grid, "sediment", &
                                                    C_INIT, 0.0_wp, huge(1.0_wp), &
                                                    0.0_wp, 0.0_wp)
      idx_sed = size(state%multilayer%tracers)
      state%multilayer%idx_sediment = idx_sed
      state%multilayer%tracers(idx_sed)%hTr = state%multilayer%h_layer*C_INIT

      ! Passive S/T (alpha_T = beta_S = 0) seeded to sane values.
      state%multilayer%tracers(state%multilayer%idx_salinity)%hTr = &
         state%multilayer%h_layer*35.0_wp
      state%multilayer%tracers(state%multilayer%idx_temperature)%hTr = &
         state%multilayer%h_layer*15.0_wp

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, &
                                RHO_0, 0.0_wp, 0.0_wp, 0.0_wp, 15.0_wp, &
                                state%grid, nz)
   end subroutine setup_sediment_state

   function total_suspended(state) result(total)
      type(state_t), intent(in) :: state
      real(wp) :: total
      integer :: i, j, k, nz
      total = 0.0_wp
      nz = state%multilayer%nz_ml
      do k = 1, nz
         do j = NGHOST + 1, state%grid%ny_total - NGHOST
            do i = NGHOST + 1, state%grid%nx_total - NGHOST
               total = total + state%multilayer%tracers(state%multilayer%idx_sediment)%hTr(i, j, k)
            end do
         end do
      end do
   end function total_suspended

   function total_bed(state) result(total)
      type(state_t), intent(in) :: state
      real(wp) :: total
      integer :: i, j
      total = 0.0_wp
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            total = total + state%multilayer%bed_sediment(i, j)
         end do
      end do
   end function total_bed

   subroutine run_to(state, t_final)
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: t_final
      real(wp) :: t
      call solver_enter_data(state)
      t = 0.0_wp
      do while (t < t_final)
         call solver_step(state, time=t)
         t = t + state%dt
      end do
      call solver_exit_data(state)
   end subroutine run_to

   subroutine check_settling(error, nz)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: nz
      type(state_t) :: state
      real(wp), parameter :: T_END = 60.0_wp
      real(wp), parameter :: REL_TOL = 1.0e-10_wp
      real(wp) :: susp0, bed0, total0, suspT, bedT

      checks: block
         call setup_sediment_state(state, nz)
         susp0 = total_suspended(state)
         bed0 = total_bed(state)
         total0 = susp0 + bed0

         call run_to(state, T_END)

         suspT = total_suspended(state)
         bedT = total_bed(state)

         ! Total (suspended + bed) sediment conserved to round-off.
         call check(error, abs((suspT + bedT) - total0)/total0 < REL_TOL, &
                    "sediment: total (suspended + bed) conserved under settling")
         if (allocated(error)) exit checks
         ! Settling actually transported mass: bed grew, suspended shrank.
         call check(error, bedT > bed0 .and. suspT < susp0, &
                    "sediment: settling deposits suspended load into the bed")
      end block checks
      call state%destroy()
   end subroutine check_settling

   subroutine test_ppm_multilayer_descent(error)
      !! PPM remap settling moves sediment more than one layer per step when
      !! w_s*dt spans several layers — a property the old donor-cell flux
      !! (each layer's downward flux capped at its own mass) provably cannot
      !! reproduce.  Release all sediment in the surface layer, take one step
      !! with w_s*dt ~ 1.6 layer thicknesses, and check meaningful mass landed
      !! >= 2 layers below the release while total (suspended + bed) is
      !! conserved.  Donor-cell would leave every layer below the next-down
      !! one exactly empty after a single step.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      integer, parameter :: NZ = 5
      real(wp), parameter :: REL_TOL = 1.0e-10_wp
      real(wp) :: dz_layer, total0, mass_lo, suspT, bedT
      integer :: i, j, k

      checks: block
         call setup_sediment_state(state, NZ)

         ! Re-seed: all sediment in the surface layer (k=NZ), bed empty.
         state%multilayer%tracers(state%multilayer%idx_sediment)%hTr = 0.0_wp
         state%multilayer%tracers(state%multilayer%idx_sediment)%hTr(:, :, NZ) = &
            state%multilayer%h_layer(:, :, NZ)*C_INIT

         ! Settle ~1.6 layer thicknesses in a single step (donor-cell caps at
         ! one layer regardless of w_s*dt).  Pin dt to a known CFL-stable value
         ! — state%dt is only computed inside the first solver_step, so it is
         ! still 0 here.
         state%dt_fixed = 3.0_wp
         dz_layer = H0/real(NZ, wp)
         state%multilayer%sediment_w_s = 1.6_wp*dz_layer/state%dt_fixed

         total0 = total_suspended(state)   ! bed starts empty

         call solver_enter_data(state)
         call solver_step(state, time=0.0_wp)
         call solver_exit_data(state)

         suspT = total_suspended(state)
         bedT = total_bed(state)

         ! Total (suspended + bed) sediment conserved to round-off.
         call check(error, abs((suspT + bedT) - total0)/total0 < REL_TOL, &
                    "ppm settling: total (suspended + bed) conserved in one step")
         if (allocated(error)) exit checks

         ! Mass reached layers k <= NZ-2 (>= 2 layers below the release).
         ! These are exactly zero under one donor-cell step.
         mass_lo = 0.0_wp
         do k = 1, NZ - 2
            do j = NGHOST + 1, state%grid%ny_total - NGHOST
               do i = NGHOST + 1, state%grid%nx_total - NGHOST
                  mass_lo = mass_lo + &
                            state%multilayer%tracers(state%multilayer%idx_sediment)%hTr(i, j, k)
               end do
            end do
         end do
         call check(error, mass_lo > 1.0e-6_wp*total0, &
                    "ppm settling: mass descended more than one layer in a single step")
      end block checks
      call state%destroy()
   end subroutine test_ppm_multilayer_descent

   subroutine test_conserve_nz5(error)
      type(error_type), allocatable, intent(out) :: error
      call check_settling(error, 5)
   end subroutine test_conserve_nz5

   subroutine test_conserve_nz1(error)
      type(error_type), allocatable, intent(out) :: error
      call check_settling(error, 1)
   end subroutine test_conserve_nz1

   subroutine setup_erosion_state(state, nz, u_init)
      !! Flat closed basin, nz sigma layers, sediment ON with a LOADED bed and
      !! NO suspended load.  Settling off (w_s=0), erosion on, Manning drag on.
      !! A uniform initial velocity u_init drives a bed stress that — for
      !! u_init big enough — exceeds tau_ce and erodes the bed.  With settling
      !! off, erosion is the only process that moves sediment.
      type(state_t), intent(inout) :: state
      integer, intent(in) :: nz
      real(wp), intent(in) :: u_init
      integer :: idx_sed

      state%use_multilayer = .true.
      state%multilayer%nz_ml = nz
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%rho_0 = RHO_0
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%kappa_S_bg = 0.0_wp
      state%multilayer%kappa_T_bg = 0.0_wp
      state%multilayer%hdiff_kappa = 0.0_wp
      state%multilayer%smagorinsky_cs = 0.0_wp
      state%multilayer%vcoord%coord_type = VCOORD_SIGMA
      state%multilayer%use_sediment = .true.
      state%multilayer%sediment_w_s = 0.0_wp     ! settling off — isolate erosion
      state%multilayer%sediment_bed_init = BED_INIT
      state%multilayer%sediment_tau_crit_erosion = TAU_CE
      state%multilayer%sediment_erosion_rate = E_RATE
      state%multilayer%sediment_porosity = 0.4_wp

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX, DX, 0.45_wp)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      state%barotropic%coriolis_f = 0.0_wp
      state%barotropic%manning_n = 0.025_wp      ! drag on -> bed stress > 0
      state%barotropic%b = -H0
      state%barotropic%h = H0
      state%barotropic%hu = H0*u_init
      state%barotropic%hv = 0.0_wp

      call ml_init_layers_from_barotropic(state%barotropic%h, &
                                          state%barotropic%hu, &
                                          state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, &
                                          state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, nz)

      call state%multilayer%register_passive_tracer(state%grid, "sediment", &
                                                    0.0_wp, 0.0_wp, huge(1.0_wp), &
                                                    0.0_wp, 0.0_wp)
      idx_sed = size(state%multilayer%tracers)
      state%multilayer%idx_sediment = idx_sed
      state%multilayer%tracers(idx_sed)%hTr = 0.0_wp   ! no initial suspended load

      state%multilayer%tracers(state%multilayer%idx_salinity)%hTr = &
         state%multilayer%h_layer*35.0_wp
      state%multilayer%tracers(state%multilayer%idx_temperature)%hTr = &
         state%multilayer%h_layer*15.0_wp

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, &
                                RHO_0, 0.0_wp, 0.0_wp, 0.0_wp, 15.0_wp, &
                                state%grid, nz)
   end subroutine setup_erosion_state

   subroutine test_erosion_fires(error)
      !! Driven flow over a loaded bed: the Manning-consistent bed stress
      !! exceeds tau_ce, so bed sediment erodes into the bottom layer.  Bed
      !! mass drops, the suspended load grows from zero, and total (suspended +
      !! bed) is conserved.  Settling is off, so erosion is the only process.
      !! A single step keeps the seeded uniform velocity from sloshing into a
      !! wall-reflected transient (the IC is unphysical for sustained flow).
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: REL_TOL = 1.0e-10_wp
      real(wp) :: bed0, susp0, total0, bedT, suspT

      checks: block
         call setup_erosion_state(state, 5, 0.5_wp)
         bed0 = total_bed(state)
         susp0 = total_suspended(state)
         total0 = bed0 + susp0

         call solver_enter_data(state)
         call solver_step(state, time=0.0_wp)
         call solver_exit_data(state)

         bedT = total_bed(state)
         suspT = total_suspended(state)

         call check(error, abs((bedT + suspT) - total0)/total0 < REL_TOL, &
                    "erosion: total (suspended + bed) conserved")
         if (allocated(error)) exit checks
         call check(error, suspT > susp0 .and. bedT < bed0, &
                    "erosion: bed stress > tau_ce lifts bed sediment into suspension")
      end block checks
      call state%destroy()
   end subroutine test_erosion_fires

   subroutine test_erosion_quiescent_none(error)
      !! No flow -> zero bed stress < tau_ce -> no erosion.  Settling is also
      !! off and there is no initial suspended load, so the loaded bed is left
      !! untouched and the water stays clear.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: REL_TOL = 1.0e-12_wp
      real(wp) :: bed0, bedT, suspT

      checks: block
         call setup_erosion_state(state, 5, 0.0_wp)   ! quiescent
         bed0 = total_bed(state)

         call run_to(state, 10.0_wp)

         bedT = total_bed(state)
         suspT = total_suspended(state)

         call check(error, abs(bedT - bed0) < REL_TOL*bed0 .and. suspT < REL_TOL, &
                    "erosion: quiescent bed (tau_b = 0 < tau_ce) does not erode")
      end block checks
      call state%destroy()
   end subroutine test_erosion_quiescent_none

end module test_sediment
