!! Regression test for the VCOORD_ZSTAR_FULL tracer-mass leak observed
!! in the Lakes Entrance 4-hour run (2026-05-11): the structured solver
!! lost ~7% of total salt mass over 4 hours of simulated tidal forcing.
!!
!! The existing `test_zstar_full_shelf_break::total_tracer_mass_conserved_*`
!! sets up a static walled basin with zero IC velocity and `beta_S = 0`,
!! which means nothing moves — the remap path is exercised but on
!! identical `dz_old` / `dz_new`, so no real remap work happens.
!!
!! This test drives a barotropic seiche by perturbing the initial η
!! (Gaussian SSH bump at the basin centre) so the barotropic step
!! actually does work and `ml_distribute_with_remap_full` is invoked
!! with non-trivial `dz_old != dz_new` pairs.
module test_zstar_full_tracer_conservation
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, VCOORD_ZSTAR_FULL, VCOORD_SIGMA
   use rki_state, only: state_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic, &
                              ml_build_vcoord_zstar_full
   use rki_ml_eos, only: ml_equation_of_state
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_vcoord, only: vcoord_target_dz_column_zstar_full, STRETCH_UNIFORM
   implicit none
   private

   public :: collect_zstar_full_tracer_conservation_tests
   public :: collect_zstar_full_seiche_test
   public :: collect_zstar_full_vdiff_test
   public :: collect_zstar_full_hdiff_test
   public :: collect_zstar_full_step_basin_test
   public :: collect_zstar_full_passive_tracer_test
   public :: collect_nz1_passive_tracer_test

   integer, parameter :: NX_PHYS = 8
   integer, parameter :: NY_PHYS = 8
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 5
   real(wp), parameter :: H0 = 20.0_wp
   real(wp), parameter :: ETA_AMPLITUDE = 2.0_wp
   real(wp), parameter :: DX_VAL = 100.0_wp
   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: S_INIT = 35.0_wp
   real(wp), parameter :: T_INIT = 15.0_wp

   real(wp), parameter :: H_SURF_TARGET = 1.0_wp
   real(wp), parameter :: H_MIN = 1.0e-4_wp
   integer, parameter :: N_SURF = 2

   ! Deep-shallow step-basin test (regression for the LE 4h leak):
   ! sharp bath contrast plus initial SSH tilt to drive sloshing across
   ! the interface.  Picked to make the bug bite hard within ~200s.
   real(wp), parameter :: H_DEEP = 50.0_wp
   real(wp), parameter :: H_SHALLOW = 2.0_wp
   real(wp), parameter :: ETA_TILT = 1.0_wp

contains

   ! Bundled collector for fpm (single-process). CMake uses the per-test
   ! collectors below — see the rki_api split for the same rationale
   ! (module-scope kernel workspaces + OpenMP-target present table leak
   ! between sequential tests in one process).
   subroutine collect_zstar_full_tracer_conservation_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("zstar_full_tracer_mass_conserved_in_seiche", &
                               test_tracer_mass_seiche), &
                  new_unittest("zstar_full_vertical_diffusion_conserves_tracer", &
                               test_vertical_diffusion_conserves), &
                  new_unittest("zstar_full_horizontal_diffusion_conserves_tracer", &
                               test_horizontal_diffusion_conserves), &
                  new_unittest("zstar_full_deep_shallow_step_basin_conserves_tracer", &
                               test_deep_shallow_step_basin), &
                  new_unittest("zstar_full_passive_tracer_tracks_salinity", &
                               test_passive_tracer_tracks_salinity), &
                  new_unittest("nz1_passive_tracer_tracks_salinity", &
                               test_nz1_passive_tracer_tracks_salinity) &
                  ]
   end subroutine collect_zstar_full_tracer_conservation_tests

   subroutine collect_zstar_full_seiche_test(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("zstar_full_tracer_mass_conserved_in_seiche", &
                                test_tracer_mass_seiche)]
   end subroutine collect_zstar_full_seiche_test

   subroutine collect_zstar_full_vdiff_test(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("zstar_full_vertical_diffusion_conserves_tracer", &
                                test_vertical_diffusion_conserves)]
   end subroutine collect_zstar_full_vdiff_test

   subroutine collect_zstar_full_hdiff_test(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("zstar_full_horizontal_diffusion_conserves_tracer", &
                                test_horizontal_diffusion_conserves)]
   end subroutine collect_zstar_full_hdiff_test

   subroutine collect_zstar_full_step_basin_test(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("zstar_full_deep_shallow_step_basin_conserves_tracer", &
                                test_deep_shallow_step_basin)]
   end subroutine collect_zstar_full_step_basin_test

   subroutine collect_zstar_full_passive_tracer_test(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("zstar_full_passive_tracer_tracks_salinity", &
                                test_passive_tracer_tracks_salinity)]
   end subroutine collect_zstar_full_passive_tracer_test

   subroutine collect_nz1_passive_tracer_test(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("nz1_passive_tracer_tracks_salinity", &
                                test_nz1_passive_tracer_tracks_salinity)]
   end subroutine collect_nz1_passive_tracer_test

   subroutine setup_zstar_full_state(state)
      !! Closed-basin state configured for VCOORD_ZSTAR_FULL with an
      !! initial SSH bump (drives barotropic seiche).
      type(state_t), intent(inout) :: state
      integer :: i, j, k
      real(wp) :: dz_col(NZ)

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_salinity = S_INIT
      state%multilayer%initial_temperature = T_INIT
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%rho_0 = RHO_0
      ! Disable mixing so any drift we observe is from the remap itself,
      ! not from PP81 / Smag / hdiff fluxes.
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%kappa_S_bg = 0.0_wp
      state%multilayer%kappa_T_bg = 0.0_wp
      state%multilayer%hdiff_kappa = 0.0_wp
      state%multilayer%smagorinsky_cs = 0.0_wp

      state%multilayer%vcoord%coord_type = VCOORD_ZSTAR_FULL
      state%multilayer%vcoord%zstar_h_surf_target = H_SURF_TARGET
      state%multilayer%vcoord%zstar_h_min = H_MIN
      state%multilayer%vcoord%zstar_stretching = STRETCH_UNIFORM
      state%multilayer%vcoord%zstar_n_surf = N_SURF

      ! DEBUG: widen tracer bounds so the clamp never fires.  If this
      ! makes the salt drift disappear, the leak is in the clamp.
      state%multilayer%S_min = -1000.0_wp
      state%multilayer%S_max = 1000.0_wp
      state%multilayer%T_min = -1000.0_wp
      state%multilayer%T_max = 1000.0_wp

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      state%barotropic%coriolis_f = 0.0_wp
      state%barotropic%manning_n = 0.0_wp

      ! Flat bed at H0 below z=0 (bed_cv = -H0).  The structured side
      ! stores b directly (positive elev = above z=0); for our setup
      ! bed elevation = -H0, so b = -H0 everywhere.
      state%barotropic%b = -H0
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%h(i, j) = H0 + ETA_AMPLITUDE &
                                       *exp(-real((i - state%grid%nx_total/2)**2 &
                                                  + (j - state%grid%ny_total/2)**2, wp) &
                                            /4.0_wp)
         end do
      end do

      ! Build the per-cell vcoord_z_ref_col from bathymetry, then lay
      ! down the initial sigma h_layer (the first solver step will
      ! remap it to the ZSTAR_FULL target).
      call ml_build_vcoord_zstar_full(state%barotropic%b, &
                                      state%multilayer%vcoord_z_ref_col, &
                                      state%grid, NZ, H_SURF_TARGET, &
                                      N_SURF, STRETCH_UNIFORM)

      call ml_init_layers_from_barotropic(state%barotropic%h, &
                                          state%barotropic%hu, &
                                          state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, &
                                          state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      ! Re-shape h_layer to ZSTAR_FULL target so IC is self-consistent.
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            call vcoord_target_dz_column_zstar_full(NZ, state%barotropic%h(i, j), &
                                                    state%multilayer%vcoord_z_ref_col(:, i, j), &
                                                    H_MIN, dz_col)
            do k = 1, NZ
               state%multilayer%h_layer(i, j, k) = dz_col(k)
            end do
         end do
      end do

      ! Uniform tracer concentration (hTr scales with h_layer)
      state%multilayer%tracers(state%multilayer%idx_salinity)%hTr = &
         state%multilayer%h_layer*S_INIT
      state%multilayer%tracers(state%multilayer%idx_temperature)%hTr = &
         state%multilayer%h_layer*T_INIT

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, &
                                RHO_0, 0.0_wp, 0.0_wp, 0.0_wp, T_INIT, &
                                state%grid, NZ)
   end subroutine setup_zstar_full_state

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

   subroutine test_tracer_mass_seiche(error)
      !! Closed walled basin with a Gaussian SSH bump → seiche → exercises
      !! `ml_distribute_with_remap_full` repeatedly with non-trivial
      !! dz_old/dz_new pairs.  No physical flux through boundaries, so
      !! domain-integrated water mass must be conserved to machine
      !! precision; tracer mass to remap-roundoff (~1e-3, matches the
      !! existing `test_zstar_ssh_tracking` lite-path tolerance of 1e-2).
      !!
      !! REGRESSION for the structured-builder sign bug found 2026-05-11:
      !! `ml_build_vcoord_zstar_full` was passing `b` raw to the column
      !! builder, which expects POSITIVE bed depth.  For a typical
      !! rakali setup where `b < 0` for submerged cells, this hit the
      !! `h_bed <= 0` degenerate branch and produced a `z_ref_col` of
      !! zeros, causing catastrophic mass blowup.  This test would have
      !! caught it.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: T_END = 30.0_wp
      real(wp), parameter :: REL_TOL_MASS = 1.0e-10_wp
      real(wp), parameter :: REL_TOL_TRACER = 1.0e-3_wp
      real(wp) :: mass0, massT, salt0, saltT, heat0, heatT

      call setup_zstar_full_state(state)

      mass0 = total_water_mass(state)
      salt0 = total_tracer_mass(state, state%multilayer%idx_salinity)
      heat0 = total_tracer_mass(state, state%multilayer%idx_temperature)

      call run_to(state, T_END)

      massT = total_water_mass(state)
      saltT = total_tracer_mass(state, state%multilayer%idx_salinity)
      heatT = total_tracer_mass(state, state%multilayer%idx_temperature)

      call check(error, abs(massT - mass0)/mass0 < REL_TOL_MASS, &
                 "ZSTAR_FULL closed basin: total water mass conserved through seiche")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      call check(error, abs(saltT - salt0)/salt0 < REL_TOL_TRACER, &
                 "ZSTAR_FULL closed basin: total salt mass conserved through remap")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      call check(error, abs(heatT - heat0)/heat0 < REL_TOL_TRACER, &
                 "ZSTAR_FULL closed basin: total heat content conserved through remap")

      call state%destroy()
   end subroutine test_tracer_mass_seiche

   subroutine setup_diffusion_state(state)
      !! Closed walled basin at rest with VCOORD_ZSTAR_FULL layers
      !! (non-uniform dz: H_SURF_TARGET=1 m at the top, thicker below)
      !! and a vertical salinity gradient.  PP81 is disabled (nu_max=0)
      !! but kappa_S_bg is enabled so the vertical-diffusion tracer
      !! kernel does real work.  beta_S = 0 keeps rho uniform so PP81
      !! does not introduce shear-dependent terms.
      type(state_t), intent(inout) :: state
      integer :: i, j, k
      real(wp) :: dz_col(NZ), S_k

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_salinity = S_INIT
      state%multilayer%initial_temperature = T_INIT
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%rho_0 = RHO_0
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%kappa_S_bg = 1.0e-3_wp
      state%multilayer%kappa_T_bg = 0.0_wp
      state%multilayer%hdiff_kappa = 0.0_wp
      state%multilayer%smagorinsky_cs = 0.0_wp

      state%multilayer%vcoord%coord_type = VCOORD_ZSTAR_FULL
      state%multilayer%vcoord%zstar_h_surf_target = H_SURF_TARGET
      state%multilayer%vcoord%zstar_h_min = H_MIN
      state%multilayer%vcoord%zstar_stretching = STRETCH_UNIFORM
      state%multilayer%vcoord%zstar_n_surf = N_SURF

      state%multilayer%S_min = -1000.0_wp
      state%multilayer%S_max = 1000.0_wp
      state%multilayer%T_min = -1000.0_wp
      state%multilayer%T_max = 1000.0_wp

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

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

      call ml_build_vcoord_zstar_full(state%barotropic%b, &
                                      state%multilayer%vcoord_z_ref_col, &
                                      state%grid, NZ, H_SURF_TARGET, &
                                      N_SURF, STRETCH_UNIFORM)

      call ml_init_layers_from_barotropic(state%barotropic%h, &
                                          state%barotropic%hu, &
                                          state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, &
                                          state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            call vcoord_target_dz_column_zstar_full(NZ, state%barotropic%h(i, j), &
                                                    state%multilayer%vcoord_z_ref_col(:, i, j), &
                                                    H_MIN, dz_col)
            do k = 1, NZ
               state%multilayer%h_layer(i, j, k) = dz_col(k)
            end do
         end do
      end do

      ! Vertical salinity gradient: 30 PSU at bed (k=1), 40 PSU at surface (k=NZ).
      ! This guarantees dtr != 0 across every interior interface, so the
      ! vertical-diffusion kernel actually moves mass — non-conservation
      ! under non-uniform dz would be immediately visible.
      do k = 1, NZ
         S_k = 30.0_wp + 10.0_wp*real(k - 1, wp)/real(NZ - 1, wp)
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = &
                  state%multilayer%h_layer(i, j, k)*S_k
               state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = &
                  state%multilayer%h_layer(i, j, k)*T_INIT
            end do
         end do
      end do

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, &
                                RHO_0, 0.0_wp, 0.0_wp, 0.0_wp, T_INIT, &
                                state%grid, NZ)
   end subroutine setup_diffusion_state

   subroutine test_vertical_diffusion_conserves(error)
      !! Regression for the dimensionally-wrong + non-conservative
      !! vertical-diffusion tracer flux found 2026-05-11.  The kernel
      !! used to multiply each side of the interface flux by its own
      !! h_k, so when h_k != h_{k+1} (typical for VCOORD_ZSTAR_FULL)
      !! the layer above gained more than the layer below lost (or
      !! vice-versa).  Symptom on LE: salt-per-water ratio 2x the
      !! ambient concentration at the open boundary.
      !!
      !! This test runs a closed walled basin at rest, with PP81
      !! disabled but kappa_S_bg active and a S(z) gradient that
      !! guarantees the kernel does work every step.  Total salt mass
      !! must be conserved to machine precision (no boundary flux,
      !! no remap because the basin is at rest, only diffusion runs).
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: T_END = 600.0_wp
      real(wp), parameter :: REL_TOL = 1.0e-12_wp
      real(wp) :: salt0, saltT

      call setup_diffusion_state(state)
      salt0 = total_tracer_mass(state, state%multilayer%idx_salinity)
      call run_to(state, T_END)
      saltT = total_tracer_mass(state, state%multilayer%idx_salinity)

      call check(error, abs(saltT - salt0)/salt0 < REL_TOL, &
                 "ZSTAR_FULL vertical diffusion: total salt mass conserved under non-uniform dz")

      call state%destroy()
   end subroutine test_vertical_diffusion_conserves

   subroutine setup_hdiff_state(state)
      !! Closed walled basin at rest with VCOORD_ZSTAR_FULL layers
      !! (non-uniform dz) and a HORIZONTAL salinity gradient.  PP81
      !! disabled, vertical diffusion disabled, only `hdiff_kappa`
      !! active so the horizontal-diffusion tracer kernel does real
      !! work.  beta_S = 0 keeps rho uniform.
      type(state_t), intent(inout) :: state
      integer :: i, j, k
      real(wp) :: dz_col(NZ), S_x

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_salinity = S_INIT
      state%multilayer%initial_temperature = T_INIT
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%rho_0 = RHO_0
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%kappa_S_bg = 0.0_wp
      state%multilayer%kappa_T_bg = 0.0_wp
      ! Active horizontal diffusion (the kernel under test).
      state%multilayer%hdiff_kappa = 1.0_wp
      state%multilayer%smagorinsky_cs = 0.0_wp

      state%multilayer%vcoord%coord_type = VCOORD_ZSTAR_FULL
      state%multilayer%vcoord%zstar_h_surf_target = H_SURF_TARGET
      state%multilayer%vcoord%zstar_h_min = H_MIN
      state%multilayer%vcoord%zstar_stretching = STRETCH_UNIFORM
      state%multilayer%vcoord%zstar_n_surf = N_SURF

      state%multilayer%S_min = -1000.0_wp
      state%multilayer%S_max = 1000.0_wp
      state%multilayer%T_min = -1000.0_wp
      state%multilayer%T_max = 1000.0_wp

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

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

      call ml_build_vcoord_zstar_full(state%barotropic%b, &
                                      state%multilayer%vcoord_z_ref_col, &
                                      state%grid, NZ, H_SURF_TARGET, &
                                      N_SURF, STRETCH_UNIFORM)

      call ml_init_layers_from_barotropic(state%barotropic%h, &
                                          state%barotropic%hu, &
                                          state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, &
                                          state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            call vcoord_target_dz_column_zstar_full(NZ, state%barotropic%h(i, j), &
                                                    state%multilayer%vcoord_z_ref_col(:, i, j), &
                                                    H_MIN, dz_col)
            do k = 1, NZ
               state%multilayer%h_layer(i, j, k) = dz_col(k)
            end do
         end do
      end do

      ! Horizontal salinity gradient: 30 PSU at the west edge, 40 PSU
      ! at the east edge.  Guarantees nabla S != 0 across every
      ! interior face, so the horizontal-diffusion kernel does work
      ! every step.  Non-conservative kernel + non-uniform dz would
      ! show up immediately as drift.
      do i = 1, state%grid%nx_total
         S_x = 30.0_wp + 10.0_wp*real(i - NGHOST - 1, wp)/real(NX_PHYS - 1, wp)
         do j = 1, state%grid%ny_total
            do k = 1, NZ
               state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = &
                  state%multilayer%h_layer(i, j, k)*S_x
               state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = &
                  state%multilayer%h_layer(i, j, k)*T_INIT
            end do
         end do
      end do

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, &
                                RHO_0, 0.0_wp, 0.0_wp, 0.0_wp, T_INIT, &
                                state%grid, NZ)
   end subroutine setup_hdiff_state

   subroutine test_horizontal_diffusion_conserves(error)
      !! Regression for the dimensionally-wrong + non-conservative
      !! horizontal-diffusion tracer flux.  The kernel previously
      !! computed `dtr = kappa·dt·∇²Tr` and updated `hTr_new = h_self
      !! · (Tr_self + dtr)` — each cell used its OWN h on the RHS, so
      !! adjacent cells with different h_layer(k) (zsigma / zstar at
      !! bathymetric transitions) saw different mass changes through
      !! the same shared face, leaking total mass.
      !!
      !! Closed walled basin at rest with VCOORD_ZSTAR_FULL
      !! (non-uniform dz) and a horizontal S gradient.  Only
      !! horizontal diffusion is active.  Total salt mass must be
      !! conserved to machine precision (no boundary flux, no remap,
      !! no advection — only the hdiff kernel runs).
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: T_END = 600.0_wp
      real(wp), parameter :: REL_TOL = 1.0e-12_wp
      real(wp) :: salt0, saltT

      call setup_hdiff_state(state)
      salt0 = total_tracer_mass(state, state%multilayer%idx_salinity)
      call run_to(state, T_END)
      saltT = total_tracer_mass(state, state%multilayer%idx_salinity)

      call check(error, abs(saltT - salt0)/salt0 < REL_TOL, &
                 "ZSTAR_FULL horizontal diffusion: total salt mass conserved under non-uniform dz")

      call state%destroy()
   end subroutine test_horizontal_diffusion_conserves

   subroutine setup_deep_shallow_step_state(state)
      !! Closed walled basin with a STEP bathymetry: left half deep
      !! (H_DEEP), right half shallow (H_SHALLOW).  Initial SSH tilt
      !! across the basin drives a sloshing barotropic flow back and
      !! forth through the deep-shallow interface — exactly the
      !! configuration that activates the LE leak under VCOORD_ZSTAR_FULL.
      !!
      !! Production clamp bounds [S_min=20, S_max=40] are LEFT AT
      !! DEFAULTS (not widened).  Initial S = 34 PSU uniform — well
      !! within bounds.  Any per-layer S drift induced by the centered-
      !! CWC over-allocation at the interface eventually pushes S out of
      !! the [20, 40] window in some layers; the clamp clips, destroying
      !! mass.  The accumulated loss is what this test detects.
      !!
      !! All other physics (PP81, hdiff, Smagorinsky, EOS coupling) are
      !! disabled so the test isolates the horizontal-advection +
      !! distribute + clamp pipeline.
      type(state_t), intent(inout) :: state
      integer :: i, j, k, mid
      real(wp) :: dz_col(NZ)
      real(wp) :: bath, x_frac

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_salinity = 34.0_wp
      state%multilayer%initial_temperature = T_INIT
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%rho_0 = RHO_0
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%kappa_S_bg = 0.0_wp
      state%multilayer%kappa_T_bg = 0.0_wp
      state%multilayer%hdiff_kappa = 0.0_wp
      state%multilayer%smagorinsky_cs = 0.0_wp

      state%multilayer%vcoord%coord_type = VCOORD_ZSTAR_FULL
      state%multilayer%vcoord%zstar_h_surf_target = H_SURF_TARGET
      state%multilayer%vcoord%zstar_h_min = H_MIN
      state%multilayer%vcoord%zstar_stretching = STRETCH_UNIFORM
      state%multilayer%vcoord%zstar_n_surf = N_SURF

      ! Production clamp window — DO NOT widen.  We need the clamp to
      ! fire on the centered-CWC induced per-layer S drift so that the
      ! test detects the LE-style mass destruction.
      state%multilayer%S_min = 20.0_wp
      state%multilayer%S_max = 40.0_wp
      state%multilayer%T_min = -10.0_wp
      state%multilayer%T_max = 50.0_wp

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      state%barotropic%coriolis_f = 0.0_wp
      state%barotropic%manning_n = 0.0_wp

      ! Step bathymetry: deep on the left, shallow on the right.  The
      ! interface sits between physical columns mid and mid+1, where
      ! mid = NX_PHYS/2 (ghost cells on the left at indices 1..NGHOST,
      ! physical columns at NGHOST+1..NGHOST+NX_PHYS).
      mid = NGHOST + NX_PHYS/2
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            if (i <= mid) then
               bath = -H_DEEP
            else
               bath = -H_SHALLOW
            end if
            state%barotropic%b(i, j) = bath
         end do
      end do

      ! IC: tilt SSH across the basin so eta varies from -ETA_TILT on
      ! the left edge to +ETA_TILT on the right edge.  This produces a
      ! barotropic pressure gradient that drives a slosh from right to
      ! left through the deep-shallow interface, then back, then again
      ! as the basin oscillates.  Many cycles → many opportunities for
      ! the per-layer flux misallocation to leak.
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            x_frac = real(i - NGHOST, wp)/real(NX_PHYS, wp) - 0.5_wp
            state%barotropic%h(i, j) = -state%barotropic%b(i, j) + ETA_TILT*x_frac
            if (state%barotropic%h(i, j) < H_MIN) state%barotropic%h(i, j) = H_MIN
         end do
      end do
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp

      call ml_build_vcoord_zstar_full(state%barotropic%b, &
                                      state%multilayer%vcoord_z_ref_col, &
                                      state%grid, NZ, H_SURF_TARGET, &
                                      N_SURF, STRETCH_UNIFORM)

      call ml_init_layers_from_barotropic(state%barotropic%h, &
                                          state%barotropic%hu, &
                                          state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, &
                                          state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      ! Reshape h_layer to ZSTAR_FULL target per cell so the IC is
      ! self-consistent with the per-cell vcoord.  k=1 is the bed,
      ! k=nz is the surface (bottom-up indexing — see Rakali style guide).
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            call vcoord_target_dz_column_zstar_full(NZ, state%barotropic%h(i, j), &
                                                    state%multilayer%vcoord_z_ref_col(:, i, j), &
                                                    H_MIN, dz_col)
            do k = 1, NZ
               state%multilayer%h_layer(i, j, k) = dz_col(k)
            end do
         end do
      end do

      ! Uniform S = 34 PSU everywhere — well inside [20, 40].  Any
      ! per-layer drift introduced by the advection+remap pipeline that
      ! pushes S out of bounds will be clipped by the clamp.
      state%multilayer%tracers(state%multilayer%idx_salinity)%hTr = &
         state%multilayer%h_layer*34.0_wp
      state%multilayer%tracers(state%multilayer%idx_temperature)%hTr = &
         state%multilayer%h_layer*T_INIT

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, &
                                RHO_0, 0.0_wp, 0.0_wp, 0.0_wp, T_INIT, &
                                state%grid, NZ)
   end subroutine setup_deep_shallow_step_state

   subroutine test_deep_shallow_step_basin(error)
      !! REGRESSION for the Lakes Entrance 4-hour salt-mass leak under
      !! VCOORD_ZSTAR_FULL (2026-05-11).  Mechanism:
      !!
      !!  1. ZSTAR_FULL has vanishing bed layers (h ≈ h_min = 1e-4 m) in
      !!     shallow cells.
      !!  2. CWC tracer advection uses centred per-layer velocity
      !!         v_xL(k) = 0.5 * (hu_layer(deep, k) + hu_layer(shal, k))
      !!     which averages incompatible layer geometries when adjacent
      !!     cells have very different h_layer distributions.
      !!  3. The shallow side's vanishing layers receive flux scaled by
      !!     ~deep/2 per layer, while they only hold h_min worth of
      !!     water.
      !!  4. After ml_distribute_with_remap_full, S_old = hS / dz_old_col
      !!     blows up to thousands of PSU on the affected layers.
      !!  5. The production [S_min=20, S_max=40] clamp clips, destroying
      !!     ~99% of the affected layer's mass per event.
      !!
      !! This test stages exactly that configuration in miniature:
      !! closed walled basin with deep-shallow step bathymetry plus an
      !! initial SSH tilt to drive sloshing flow across the interface.
      !! Production clamp bounds are active.  Uniform initial S = 34 PSU.
      !!
      !! No physical flux through boundaries → total salt mass MUST be
      !! conserved.  Centred-CWC fails by losing ~%-level mass over a
      !! few hundred steps.  A fix that respects per-cell layer budgets
      !! (e.g. shear-aware self-frac:
      !!   f_layer(k) = hu_self(k) + frac_self(k) * (mass_flux - HU_self)
      !! ) should keep the drift below the tolerance.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: T_END = 200.0_wp
      real(wp), parameter :: REL_TOL_TRACER = 1.0e-4_wp
      real(wp), parameter :: REL_TOL_MASS = 1.0e-10_wp
      real(wp) :: mass0, massT, salt0, saltT

      call setup_deep_shallow_step_state(state)

      mass0 = total_water_mass(state)
      salt0 = total_tracer_mass(state, state%multilayer%idx_salinity)

      call run_to(state, T_END)

      massT = total_water_mass(state)
      saltT = total_tracer_mass(state, state%multilayer%idx_salinity)

      call check(error, abs(massT - mass0)/mass0 < REL_TOL_MASS, &
                 "ZSTAR_FULL deep-shallow step basin: water mass conserved")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      call check(error, abs(saltT - salt0)/salt0 < REL_TOL_TRACER, &
                 "ZSTAR_FULL deep-shallow step basin: salt mass conserved "// &
                 "(regression for LE 4h leak; centered-CWC fails this, "// &
                 "per-cell-frac variants should pass)")

      call state%destroy()
   end subroutine test_deep_shallow_step_basin

   subroutine test_passive_tracer_tracks_salinity(error)
      !! Tracer-registry generalization: register a 3rd PASSIVE tracer
      !! with salinity's IC + transport config and run the seiche.  With
      !! beta_S = 0 salinity is effectively passive (EOS reads but never
      !! writes it), and the dye goes through the identical generic
      !! transport (BCs, advection, ALE remap, vertical exchange, clamp,
      !! RK2 save/average), so it must track salinity bit-for-bit.  Any
      !! kernel that silently skips the 3rd tracer (e.g. an
      !! un-generalized remap) makes the two diverge — strictly stronger
      !! than a bare sum-conservation check, which a skipped tracer would
      !! falsely pass.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: T_END = 30.0_wp
      real(wp), parameter :: TOL = 1.0e-9_wp
      integer :: idx_dye, idx_s
      real(wp) :: max_diff

      call setup_zstar_full_state(state)

      ! Append a passive dye with salinity's widened bounds + zero
      ! diffusivity (matching the seiche setup), seeded identically to S.
      call state%multilayer%register_passive_tracer(state%grid, "dye", S_INIT, &
                                                    -1000.0_wp, 1000.0_wp, &
                                                    0.0_wp, 0.0_wp)
      idx_dye = size(state%multilayer%tracers)
      idx_s = state%multilayer%idx_salinity
      state%multilayer%tracers(idx_dye)%hTr = state%multilayer%h_layer*S_INIT

      call run_to(state, T_END)

      max_diff = maxval(abs(state%multilayer%tracers(idx_dye)%hTr &
                            - state%multilayer%tracers(idx_s)%hTr))
      call check(error, max_diff < TOL, &
                 "passive tracer must track salinity through the generic "// &
                 "transport pipeline (advection + remap + exchange + clamp + RK2)")

      call state%destroy()
   end subroutine test_passive_tracer_tracks_salinity

   subroutine setup_nz1_state(state)
      !! Single-layer (nz_ml=1) closed walled basin, sigma coordinate, with
      !! a Gaussian SSH bump that drives a seiche -- exercises the nz=1
      !! fused fast path (solver_multilayer_rk_stage_nz1) under active
      !! barotropic flow.  beta_S = 0 (salinity passive), no mixing, widened
      !! tracer bounds so the clamp never fires.
      type(state_t), intent(inout) :: state
      integer :: i, j

      state%use_multilayer = .true.
      state%multilayer%nz_ml = 1
      state%multilayer%initial_salinity = S_INIT
      state%multilayer%initial_temperature = T_INIT
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

      state%multilayer%S_min = -1000.0_wp
      state%multilayer%S_max = 1000.0_wp
      state%multilayer%T_min = -1000.0_wp
      state%multilayer%T_max = 1000.0_wp

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      state%barotropic%coriolis_f = 0.0_wp
      state%barotropic%manning_n = 0.0_wp
      state%barotropic%b = -H0
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%h(i, j) = H0 + ETA_AMPLITUDE &
                                       *exp(-real((i - state%grid%nx_total/2)**2 &
                                                  + (j - state%grid%ny_total/2)**2, wp) &
                                            /4.0_wp)
         end do
      end do

      call ml_init_layers_from_barotropic(state%barotropic%h, &
                                          state%barotropic%hu, &
                                          state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, &
                                          state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, 1)

      state%multilayer%tracers(state%multilayer%idx_salinity)%hTr = &
         state%multilayer%h_layer*S_INIT
      state%multilayer%tracers(state%multilayer%idx_temperature)%hTr = &
         state%multilayer%h_layer*T_INIT

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, &
                                RHO_0, 0.0_wp, 0.0_wp, 0.0_wp, T_INIT, &
                                state%grid, 1)
   end subroutine setup_nz1_state

   subroutine test_nz1_passive_tracer_tracks_salinity(error)
      !! Tracer-registry generalization on the nz=1 fused fast path: register
      !! a 3rd PASSIVE tracer seeded like salinity and run a single-layer
      !! seiche.  With beta_S = 0 salinity is passive, so the dye must track
      !! it bit-for-bit through the fused-path transport (tracer BCs +
      !! advection + clamp; RK2 average via ml_fused_rk2_finalise).  Before
      !! the nz=1 path was generalized it transported S/T only, leaving the
      !! dye frozen -- this test catches that.  Companion to
      !! test_passive_tracer_tracks_salinity (which covers the nz>1 path).
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: T_END = 30.0_wp
      real(wp), parameter :: TOL = 1.0e-9_wp
      integer :: idx_dye, idx_s
      real(wp) :: max_diff

      call setup_nz1_state(state)

      ! Append a passive dye with salinity's IC + transport config, seeded
      ! identically to S.
      call state%multilayer%register_passive_tracer(state%grid, "dye", S_INIT, &
                                                    -1000.0_wp, 1000.0_wp, &
                                                    0.0_wp, 0.0_wp)
      idx_dye = size(state%multilayer%tracers)
      idx_s = state%multilayer%idx_salinity
      state%multilayer%tracers(idx_dye)%hTr = state%multilayer%h_layer*S_INIT

      call run_to(state, T_END)

      max_diff = maxval(abs(state%multilayer%tracers(idx_dye)%hTr &
                            - state%multilayer%tracers(idx_s)%hTr))
      call check(error, max_diff < TOL, &
                 "nz=1 fused path: passive tracer must track salinity "// &
                 "through BCs + advection + clamp + RK2")

      call state%destroy()
   end subroutine test_nz1_passive_tracer_tracks_salinity

end module test_zstar_full_tracer_conservation
