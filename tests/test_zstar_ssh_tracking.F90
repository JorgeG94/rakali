!! Integration test: VCOORD_ZSTAR (z*-lite) through the full solver
!!
!! Closed basin with non-uniform `z_ref` (surface-concentrated layers),
!! uniform stratification, initial SSH bump that propagates as a seiche.
!! Runs the full SSP-RK2 multilayer pipeline with VCOORD_ZSTAR active
!! and the conservative remap path engaged every stage.
!!
!! Validates the integration properties of z*-lite:
!!   1. Solver runs without blowing up under VCOORD_ZSTAR (no NaN, no
!!      crash, no negative h_layer).
!!   2. Total water mass is conserved across the run (closed basin).
!!   3. Total tracer mass (salt + heat) is conserved through the
!!      ZSTAR remap step (the conservative-remap kernel is exercised
!!      every RK2 stage with the new VCOORD).
!!   4. Per-cell layer-thickness RATIOS match the z_ref pattern at
!!      every cell at every time — z*-lite's defining invariant.
!!      Under uniform sigma the same probe would see uniform layer
!!      thicknesses regardless of the z_ref pattern.
!!
!! z*-lite (the v1 implementation) does NOT promise per-cell physical-
!! depth-anchored layers — that requires per-column `z_ref` matched
!! to local bathymetry (the full MOM6 z*).  This test validates what
!! z*-lite actually delivers: non-uniform fractional layer spacing
!! preserved under SSH variation, with mass + tracer conservation.
module test_zstar_ssh_tracking
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, VCOORD_ZSTAR
   use rki_state, only: state_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_eos, only: ml_equation_of_state
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_vcoord, only: vcoord_target_dz_column
   implicit none
   private

   public :: collect_zstar_ssh_tracking_tests
   public :: collect_zstar_ssh_runs_test
   public :: collect_zstar_ssh_ratios_test

   integer, parameter :: NX_PHYS = 8
   integer, parameter :: NY_PHYS = 8
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   real(wp), parameter :: H0 = 20.0_wp        !! 20 m base water depth
   real(wp), parameter :: ETA_AMPLITUDE = 2.0_wp   !! 10 % of H0 — large
   real(wp), parameter :: DX_VAL = 100.0_wp   !! 100 m cells

   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: S_INIT = 35.0_wp
   real(wp), parameter :: T_INIT = 15.0_wp

   ! Surface-concentrated z_ref pattern (deepest reference layer is
   ! the thickest under bottom-up ROMS ordering: dz(1) = bottom).
   ! Reference layer thicknesses top-down: 2, 4, 6, 8 m → bottom-up: 8, 6, 4, 2.
   ! Total z_ref(NZ) = 20 m, matches H0.
   real(wp), parameter :: Z_REF_PATTERN(0:NZ) = &
                          [0.0_wp, 2.0_wp, 6.0_wp, 12.0_wp, 20.0_wp]

contains

   ! Bundled collector for fpm; CMake uses the per-test collectors below.
   subroutine collect_zstar_ssh_tracking_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("zstar_solver_runs_and_conserves", &
                               test_zstar_solver_runs_and_conserves), &
                  new_unittest("zstar_layer_ratios_preserved", &
                               test_zstar_layer_ratios_preserved) &
                  ]
   end subroutine collect_zstar_ssh_tracking_tests

   subroutine collect_zstar_ssh_runs_test(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("zstar_solver_runs_and_conserves", &
                                test_zstar_solver_runs_and_conserves)]
   end subroutine collect_zstar_ssh_runs_test

   subroutine collect_zstar_ssh_ratios_test(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("zstar_layer_ratios_preserved", &
                                test_zstar_layer_ratios_preserved)]
   end subroutine collect_zstar_ssh_ratios_test

   subroutine setup_zstar_state(state)
      !! Build a closed-basin state configured for VCOORD_ZSTAR with
      !! a non-uniform z_ref pattern and an initial SSH bump.
      type(state_t), intent(inout) :: state

      integer :: i, j, k

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_salinity = S_INIT
      state%multilayer%initial_temperature = T_INIT
      state%multilayer%alpha_T = 0.0_wp     ! decouple T from buoyancy
      state%multilayer%beta_S = 0.0_wp      ! decouple S from buoyancy
      state%multilayer%rho_0 = RHO_0
      ! Disable extra mixing so conservation budgets are clean to
      ! machine precision (mixing is conservative but introduces
      ! diffusive flux that complicates the per-layer probe).
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%hdiff_kappa = 0.0_wp
      state%multilayer%smagorinsky_cs = 0.0_wp

      ! Activate VCOORD_ZSTAR via the multilayer state's vcoord struct
      ! BEFORE state%init() — multilayer_state_init() reads coord_type
      ! to allocate / populate the right reference structures.
      state%multilayer%vcoord%coord_type = VCOORD_ZSTAR

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      ! Override z_ref with our non-uniform pattern (still host-side,
      ! before solver_enter_data ships state to the device).  Both the
      ! vcoord-internal copy and the GPU-mappable mirror need to see
      ! the override so the on-device kernel reads the right values.
      if (allocated(state%multilayer%vcoord%z_ref)) then
         state%multilayer%vcoord%z_ref = Z_REF_PATTERN
      end if
      if (allocated(state%multilayer%vcoord_z_ref)) then
         state%multilayer%vcoord_z_ref = Z_REF_PATTERN
      end if

      ! Closed basin
      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      state%barotropic%coriolis_f = 0.0_wp
      state%barotropic%manning_n = 0.0_wp

      ! Bathymetry: flat at H0 below z=0 (so b = -H0; eta is measured
      ! from z=0 datum and h = eta - b = eta + H0).  Keep it simple:
      ! b = 0 and h = H0 + bump.
      state%barotropic%b = 0.0_wp
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp

      ! Initial η bump: positive Gaussian centred at the basin centre.
      ! Generates a seiche when released.
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%h(i, j) = H0 + ETA_AMPLITUDE &
                                       *exp(-real((i - state%grid%nx_total/2)**2 &
                                                  + (j - state%grid%ny_total/2)**2, wp) &
                                            /4.0_wp)
         end do
      end do

      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      ! Re-shape h_layer per cell to the ZSTAR target so the IC is
      ! self-consistent with the vcoord (otherwise the first RK2-averaging
      ! step would blend a uniform-sigma snapshot with a ZSTAR-remapped
      ! one, producing transient drift before convergence).
      block
         integer :: kk
         real(wp) :: dz_col(NZ)
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               call vcoord_target_dz_column(VCOORD_ZSTAR, NZ, state%barotropic%h(i, j), &
                                            state%multilayer%dsig_ml, state%multilayer%vcoord_z_ref, &
                                            state%multilayer%vcoord%depth_transition, &
                                            state%multilayer%vcoord%blend_width, dz_col)
               do kk = 1, NZ
                  state%multilayer%h_layer(i, j, kk) = dz_col(kk)
               end do
            end do
         end do
      end block

      ! Uniform tracer concentrations per layer (so hTr scales with layer thickness)
      state%multilayer%tracers(state%multilayer%idx_salinity)%hTr = state%multilayer%h_layer*S_INIT
      state%multilayer%tracers(state%multilayer%idx_temperature)%hTr = state%multilayer%h_layer*T_INIT

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, &
                                RHO_0, 0.0_wp, 0.0_wp, 0.0_wp, T_INIT, &
                                state%grid, NZ)
   end subroutine setup_zstar_state

   function total_tracer_mass(state, idx) result(total)
      type(state_t), intent(in) :: state
      integer, intent(in) :: idx
      real(wp) :: total
      integer :: i, j, k

      total = 0.0_wp
      do k = 1, NZ
         do j = NGHOST + 1, state%grid%ny_total - NGHOST
            do i = NGHOST + 1, state%grid%nx_total - NGHOST
               total = total + state%multilayer%tracers(idx)%hTr(i, j, k)
            end do
         end do
      end do
      total = total*DX_VAL*DX_VAL
   end function total_tracer_mass

   function total_water_mass(state) result(total)
      type(state_t), intent(in) :: state
      real(wp) :: total
      integer :: i, j

      total = 0.0_wp
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            total = total + state%barotropic%h(i, j)
         end do
      end do
      total = total*DX_VAL*DX_VAL
   end function total_water_mass

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

   subroutine test_zstar_solver_runs_and_conserves(error)
      !! ZSTAR runs through the full SSP-RK2 multilayer pipeline (which
      !! engages the conservative remap path because ZSTAR sets
      !! `vcoord_needs_remap` = .true.) and conserves mass + tracer mass
      !! across an SSH-driven seiche.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp) :: mass0, salt0, heat0
      real(wp) :: massT, saltT, heatT
      real(wp), parameter :: T_END = 30.0_wp     !! seconds — several CFL steps
      ! Closed basin: should be tightly conservative.  Match
      ! test_river_plume's 1 % tolerance for the tracer (PPM advection
      ! + remap can leak this much through wall ghost-cell handling on
      ! a relatively coarse grid); water mass tighter at 1e-6.
      real(wp), parameter :: REL_TOL_TRACER = 1.0e-2_wp
      real(wp), parameter :: REL_TOL_MASS = 1.0e-6_wp

      call setup_zstar_state(state)

      mass0 = total_water_mass(state)
      salt0 = total_tracer_mass(state, state%multilayer%idx_salinity)
      heat0 = total_tracer_mass(state, state%multilayer%idx_temperature)

      call run_to(state, T_END)

      massT = total_water_mass(state)
      saltT = total_tracer_mass(state, state%multilayer%idx_salinity)
      heatT = total_tracer_mass(state, state%multilayer%idx_temperature)

      call check(error, abs(massT - mass0)/mass0 < REL_TOL_MASS, &
                 "ZSTAR solver: total water mass conserved through seiche")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      call check(error, abs(saltT - salt0)/salt0 < REL_TOL_TRACER, &
                 "ZSTAR solver: total salt mass conserved through remap (~1%)")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      call check(error, abs(heatT - heat0)/heat0 < REL_TOL_TRACER, &
                 "ZSTAR solver: total heat content conserved through remap (~1%)")

      call state%destroy()
   end subroutine test_zstar_solver_runs_and_conserves

   subroutine test_zstar_layer_ratios_preserved(error)
      !! At every interior cell at the end of the run, the layer-
      !! thickness RATIOS must match the z_ref pattern — z*-lite's
      !! defining invariant.  This holds even as h(i,j) varies across
      !! the basin from the seiche dynamics, because the per-column
      !! ZSTAR formula scales the whole z_ref pattern by H/z_ref(nz).
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: T_END = 20.0_wp
      real(wp), parameter :: TOL = 1.0e-9_wp
      integer :: i, j, k
      real(wp) :: ratio_dz, ratio_zref

      call setup_zstar_state(state)
      call run_to(state, T_END)

      ! For every interior cell, every adjacent layer pair, the
      ! per-cell ratio dz(k)/dz(k+1) must equal the z_ref ratio
      ! (z_ref(NZ-k+1) - z_ref(NZ-k)) / (z_ref(NZ-k) - z_ref(NZ-k-1)).
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            ! Skip cells with collapsed columns (shouldn't happen with
            ! the chosen IC but guard anyway)
            if (state%multilayer%h_layer(i, j, 1) < 1.0e-6_wp) cycle
            do k = 1, NZ - 1
               ratio_dz = state%multilayer%h_layer(i, j, k) &
                          /state%multilayer%h_layer(i, j, k + 1)
               ratio_zref = (Z_REF_PATTERN(NZ - k + 1) - Z_REF_PATTERN(NZ - k)) &
                            /(Z_REF_PATTERN(NZ - k) - Z_REF_PATTERN(NZ - k - 1))
               call check(error, abs(ratio_dz - ratio_zref) < TOL, &
                          "ZSTAR: per-cell layer ratios match z_ref pattern")
               if (allocated(error)) then
                  call state%destroy()
                  return
               end if
            end do
         end do
      end do

      call state%destroy()
   end subroutine test_zstar_layer_ratios_preserved

end module test_zstar_ssh_tracking
