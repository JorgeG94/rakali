!! Unit tests for coupled hydrostatic multi-layer kernels
module test_multilayer
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, CP_WATER
   use rki_ml_dynamics, only: ml_distribute_barotropic, &
                              ml_compute_total_depth, &
                              ml_compute_barotropic_momentum, &
                              ml_compute_vertical_velocity, &
                              ml_compute_vertical_exchange, &
                              ml_rk2_average_layers, &
                              ml_init_layers_from_barotropic, &
                              ml_save_layer_state
   use rki_ml_tracers, only: ml_init_tracer, &
                             ml_init_tracer_stratified, &
                             ml_distribute_tracer, &
                             ml_vertical_exchange_tracer, &
                             ml_vertical_diffusion_tracer, &
                             ml_clamp_tracer, &
                             ml_apply_surface_tracer_fluxes, &
                             ml_advect_tracer, &
                             ml_advect_tracer_ppm, &
                             ml_horizontal_diffusion_tracer, &
                             tracer_scratch_t, &
                             tracer_scratch_ensure, &
                             tracer_scratch_destroy
   use rki_ml_eos, only: ml_equation_of_state
   use rki_ml_pressure, only: ml_baroclinic_pressure_gradient, &
                              ml_baroclinic_pressure_gradient_jacobian
   use rki_ml_vertical_viscosity, only: ml_vertical_viscosity
   use rki_ml_horizontal_viscosity, only: ml_horizontal_viscosity_smag
   use rki_ml_friction, only: ml_apply_bottom_friction
   use rki_ml_truncate, only: ml_truncate_velocity
   use rki_multilayer_state, only: multilayer_state_t
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: collect_multilayer_tests

   integer, parameter :: NX = 12, NY = 12, NG = 2, NZ = 2
   real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp

   type(hgrid_t) :: test_grid
   logical :: test_grid_initialised = .false.

   ! Shared caller-owned scratch for the `pure` tracer kernels; each test
   ! ensures the groups it needs then destroys them (module-level to avoid
   ! threading a local through every tracer test).
   type(tracer_scratch_t) :: scratch

   ! Dummy reference-density placeholders for tests that don't exercise
   ! the Phase B2 anomaly (rho_ref_dz=0 makes the kernel fall back to
   ! rho_0). Tests that probe B2 should build a real profile via
   ! ml_compute_reference_density.
   real(wp), parameter :: TEST_RHO_REF(0:0) = [0.0_wp]
   real(wp), parameter :: TEST_RHO_REF_DZ = 0.0_wp
   integer, parameter :: TEST_NZ_REF = 0

contains

   subroutine ensure_grid()
      if (.not. test_grid_initialised) then
         call test_grid%init(NX - 2*NG, NY - 2*NG, NG, DX, DY)
         test_grid_initialised = .true.
      end if
   end subroutine ensure_grid

   subroutine collect_multilayer_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("total_depth_from_layers", test_total_depth), &
                  new_unittest("barotropic_momentum", test_barotropic_momentum), &
                  new_unittest("init_layers_from_barotropic", test_init_layers), &
                  new_unittest("save_layer_state", test_save_layers), &
                  new_unittest("rk2_average_layers", test_rk2_avg), &
                  new_unittest("distribute_uniform", test_distribute_uniform), &
                  new_unittest("distribute_preserves_shear", test_distribute_shear), &
                  new_unittest("vertical_velocity_uniform", test_w_uniform), &
                  new_unittest("vertical_exchange_no_w", test_vexch_no_w), &
                  new_unittest("vertical_exchange_transfers", test_vexch_transfers), &
                  new_unittest("init_salinity", test_init_salinity), &
                  new_unittest("init_tracer_temperature", test_init_tracer_temperature), &
                  new_unittest("init_tracer_stratified_linear_profile", test_init_tracer_stratified), &
                  new_unittest("eos_linear", test_eos_linear), &
                  new_unittest("eos_two_tracer", test_eos_two_tracer), &
                  new_unittest("eos_alpha_T_zero_regression", test_eos_alpha_T_zero_regression), &
                  new_unittest("eos_S_ref_nonzero", test_eos_S_ref_nonzero), &
                  new_unittest("clamp_tracer_both_tracers", test_clamp_tracer_both_tracers), &
                  new_unittest("distribute_salinity_preserves_S", test_distribute_salinity), &
                  new_unittest("vertical_exchange_salinity", test_vexch_salinity), &
                  new_unittest("bpg_zero_uniform_density", test_bpg_uniform), &
                  new_unittest("bpg_lock_exchange", test_bpg_lock_exchange), &
                  new_unittest("vmix_zero_shear", test_vmix_zero_shear), &
                  new_unittest("vmix_unstratified_shear", test_vmix_unstratified_shear), &
                  new_unittest("vmix_strong_stratification", test_vmix_strong_strat), &
                  new_unittest("vmix_salinity_mixes", test_vmix_salinity), &
                  new_unittest("jacobian_zero_uniform", test_jacobian_uniform), &
                  new_unittest("jacobian_zero_varying_depth", test_jacobian_varying_depth), &
                  new_unittest("jacobian_lock_exchange", test_jacobian_lock_exchange), &
                  new_unittest("smag_uniform_no_change", test_smag_uniform), &
                  new_unittest("smag_shear_damps", test_smag_shear), &
                  new_unittest("bottom_friction_damps", test_bottom_friction), &
                  new_unittest("surface_heat_flux_warms_top", test_surface_heat_flux), &
                  new_unittest("surface_salt_flux_salts_top", test_surface_salt_flux), &
                  new_unittest("advect_zero_velocity_no_change", &
                               test_advect_zero_velocity_no_change), &
                  new_unittest("advect_uniform_tracer_no_change", &
                               test_advect_uniform_tracer_no_change), &
                  new_unittest("advect_upwind_direction", test_advect_upwind_direction), &
                  new_unittest("hdiff_zero_kappa_no_change", test_hdiff_zero_kappa_no_change), &
                  new_unittest("hdiff_uniform_no_change", test_hdiff_uniform_no_change), &
                  new_unittest("hdiff_step_diffuses_symmetrically", &
                               test_hdiff_step_diffuses_symmetrically), &
                  new_unittest("truncate_below_threshold_untouched", &
                               test_truncate_below_threshold_untouched), &
                  new_unittest("truncate_above_threshold_clipped", &
                               test_truncate_above_threshold_clipped), &
                  new_unittest("truncate_per_cell_independence", &
                               test_truncate_per_cell_independence), &
                  new_unittest("ppm_zero_velocity_no_change", &
                               test_ppm_zero_velocity_no_change), &
                  new_unittest("ppm_uniform_tracer_no_change", &
                               test_ppm_uniform_tracer_no_change), &
                  new_unittest("ppm_upwind_direction", &
                               test_ppm_upwind_direction), &
                  new_unittest("ppm_less_diffusive_than_upwind", &
                               test_ppm_less_diffusive_than_upwind), &
                  new_unittest("tracer_registry_lock_refuses", &
                               test_tracer_registry_lock_refuses) &
                  ]
   end subroutine collect_multilayer_tests

   subroutine test_tracer_registry_lock_refuses(error)
      !! Tracer-registry lock guard: registering BEFORE the device map
      !! grows the registry by 1; once `registry_locked` is set (as
      !! `solver_tracers_enter_data` does after the element-by-element
      !! device map) registration is REFUSED — size unchanged and the
      !! existing tracer names intact.  Host-only: the lock is set
      !! directly, no GPU or solver run required.
      type(error_type), allocatable, intent(out) :: error

      type(multilayer_state_t) :: ms
      integer :: n_before

      call ensure_grid()

      checks: block
         ms%nz_ml = NZ
         call ms%init(test_grid)
         n_before = size(ms%tracers) ! salinity + temperature

         ! Before the lock: registration succeeds and grows the registry.
         call ms%register_passive_tracer(test_grid, "dye", 1.0_wp, 0.0_wp, &
                                         10.0_wp, 0.0_wp, 0.0_wp)
         call check(error, size(ms%tracers) == n_before + 1, &
                    "register before the lock must grow the registry by 1")
         if (allocated(error)) exit checks

         ! Locked: registration must be refused without touching the registry.
         ms%registry_locked = .true.
         call ms%register_passive_tracer(test_grid, "dye2", 1.0_wp, 0.0_wp, &
                                         10.0_wp, 0.0_wp, 0.0_wp)
         call check(error, size(ms%tracers) == n_before + 1, &
                    "register after the lock must leave the registry size unchanged")
         if (allocated(error)) exit checks
         call check(error, trim(ms%tracers(ms%idx_salinity)%name) == "salinity" .and. &
                    trim(ms%tracers(ms%idx_temperature)%name) == "temperature" .and. &
                    trim(ms%tracers(n_before + 1)%name) == "dye", &
                    "a refused registration must leave existing tracer names intact")
      end block checks

      call ms%destroy()
      if (.not. allocated(error)) then
         call check(error,.not. ms%registry_locked, &
                    "destroy must clear registry_locked")
      end if
   end subroutine test_tracer_registry_lock_refuses

   subroutine test_total_depth(error)
      !! H = sum_k h_k
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), h(NX, NY)

      h_layer(:, :, 1) = 1.5_wp
      h_layer(:, :, 2) = 0.5_wp

      !$acc data copyin( h_layer) copy( h)
      call ensure_grid()
      call ml_compute_total_depth(h_layer, h, test_grid, NZ)
      !$acc end data

      call check(error, abs(h(5, 5) - 2.0_wp) < 1.0e-12_wp, &
                 "H = 1.5 + 0.5 = 2.0")
   end subroutine test_total_depth

   subroutine test_barotropic_momentum(error)
      !! hu_total = sum_k hu_k
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: hu(NX, NY), hv(NX, NY)

      h_layer(:, :, 1) = 1.0_wp
      h_layer(:, :, 2) = 1.0_wp
      hu_layer(:, :, 1) = 2.0_wp
      hu_layer(:, :, 2) = 4.0_wp
      hv_layer(:, :, 1) = -1.0_wp
      hv_layer(:, :, 2) = 1.0_wp

      !$acc data copyin( h_layer, hu_layer, hv_layer) copy( hu, hv)
      call ensure_grid()
      call ml_compute_barotropic_momentum(h_layer, hu_layer, hv_layer, &
                                          hu, hv, test_grid, NZ)
      !$acc end data

      call check(error, abs(hu(5, 5) - 6.0_wp) < 1.0e-12_wp, &
                 "hu = 2 + 4 = 6")
      if (allocated(error)) return
      call check(error, abs(hv(5, 5) - 0.0_wp) < 1.0e-12_wp, &
                 "hv = -1 + 1 = 0")
   end subroutine test_barotropic_momentum

   subroutine test_init_layers(error)
      !! Initialise layers from barotropic: h_k = dsig_k * H
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), dsig(NZ)
      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)

      h = 4.0_wp
      hu = 2.0_wp
      hv = -1.0_wp
      dsig = 1.0_wp/real(NZ, wp)

      !$acc data copyin( h, hu, hv, dsig) copy( h_layer, hu_layer, hv_layer)
      call ensure_grid()
      call ml_init_layers_from_barotropic(h, hu, hv, dsig, &
                                          h_layer, hu_layer, hv_layer, &
                                          test_grid, NZ)
      !$acc end data

      call check(error, abs(h_layer(5, 5, 1) - 2.0_wp) < 1.0e-12_wp, &
                 "h_layer(k=1) = 0.5 * 4 = 2")
      if (allocated(error)) return
      call check(error, abs(hu_layer(5, 5, 1) - 1.0_wp) < 1.0e-12_wp, &
                 "hu_layer(k=1) = 0.5 * 2 = 1")
   end subroutine test_init_layers

   subroutine test_save_layers(error)
      !! Save layer state copies to saved arrays
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: h0(NX, NY, NZ), hu0(NX, NY, NZ), hv0(NX, NY, NZ)

      h_layer = 3.0_wp
      hu_layer = 1.0_wp
      hv_layer = -2.0_wp
      h0 = 0.0_wp
      hu0 = 0.0_wp
      hv0 = 0.0_wp

      !$acc data copyin( h_layer, hu_layer, hv_layer) &
      !$acc&            copy( h0, hu0, hv0)
      call ensure_grid()
      call ml_save_layer_state(h_layer, hu_layer, hv_layer, &
                               h0, hu0, hv0, test_grid, NZ)
      !$acc end data

      call check(error, abs(h0(5, 5, 2) - 3.0_wp) < 1.0e-12_wp, &
                 "saved h_layer")
      if (allocated(error)) return
      call check(error, abs(hv0(5, 5, 1) - (-2.0_wp)) < 1.0e-12_wp, &
                 "saved hv_layer")
   end subroutine test_save_layers

   subroutine test_rk2_avg(error)
      !! RK2 average: U = 0.5*(U0 + U)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: h0(NX, NY, NZ), hu0(NX, NY, NZ), hv0(NX, NY, NZ)

      h_layer = 4.0_wp
      hu_layer = 2.0_wp
      hv_layer = -2.0_wp
      h0 = 2.0_wp
      hu0 = 0.0_wp
      hv0 = 0.0_wp

      !$acc data copy( h_layer, hu_layer, hv_layer) &
      !$acc&            copyin( h0, hu0, hv0)
      call ensure_grid()
      call ml_rk2_average_layers(h_layer, hu_layer, hv_layer, &
                                 h0, hu0, hv0, test_grid, NZ)
      !$acc end data

      call check(error, abs(h_layer(5, 5, 1) - 3.0_wp) < 1.0e-12_wp, &
                 "h = 0.5*(2+4) = 3")
      if (allocated(error)) return
      call check(error, abs(hu_layer(5, 5, 1) - 1.0_wp) < 1.0e-12_wp, &
                 "hu = 0.5*(0+2) = 1")
   end subroutine test_rk2_avg

   subroutine test_distribute_uniform(error)
      !! Uniform velocity: after barotropic update, all layers get
      !! the new barotropic momentum proportionally.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: hu_pre(NX, NY), hv_pre(NX, NY)
      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: dsig(NZ)

      dsig = 1.0_wp/real(NZ, wp)

      ! Before barotropic update: H=2, HU=4 (U=2), uniform layers
      hu_pre = 4.0_wp
      hv_pre = 0.0_wp
      h_layer(:, :, 1) = 1.0_wp   ! dsig * H_old = 0.5 * 2
      h_layer(:, :, 2) = 1.0_wp
      hu_layer(:, :, 1) = 2.0_wp  ! h_k * u = 1.0 * 2.0
      hu_layer(:, :, 2) = 2.0_wp
      hv_layer = 0.0_wp

      ! After barotropic update: H=2, HU=6 (U=3, accelerated by 1)
      h = 2.0_wp
      hu = 6.0_wp
      hv = 0.0_wp

      !$acc data copyin( h, hu, hv, hu_pre, hv_pre, dsig) &
      !$acc&            copy( h_layer, hu_layer, hv_layer)
      call ensure_grid()
      call ml_distribute_barotropic(h, hu, hv, hu_pre, hv_pre, dsig, &
                                    h_layer, hu_layer, hv_layer, &
                                    test_grid, NZ)
      !$acc end data

      ! h_k = dsig * H_new = 0.5 * 2 = 1.0
      call check(error, abs(h_layer(5, 5, 1) - 1.0_wp) < 1.0e-12_wp, &
                 "h_layer = dsig * H = 1.0")
      if (allocated(error)) return
      ! u_k_new = u_k_old + du = 2.0 + (6-4)/2 = 3.0
      ! hu_k_new = h_k * u_k_new = 1.0 * 3.0 = 3.0
      call check(error, abs(hu_layer(5, 5, 1) - 3.0_wp) < 1.0e-12_wp, &
                 "hu_layer = h_k * (u_old + du) = 1.0 * 3.0")
      if (allocated(error)) return
      call check(error, abs(hu_layer(5, 5, 2) - 3.0_wp) < 1.0e-12_wp, &
                 "hu_layer(k=2) same as k=1 for uniform velocity")
   end subroutine test_distribute_uniform

   subroutine test_distribute_shear(error)
      !! Sheared velocity: layers preserve their velocity anomaly
      !! relative to the depth average after barotropic acceleration.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: hu_pre(NX, NY), hv_pre(NX, NY)
      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: dsig(NZ)
      real(wp) :: u1_new, u2_new

      dsig = 1.0_wp/real(NZ, wp)

      ! Before: H=2, layers have shear: u1=1, u2=3 => U_bar=2, HU=4
      hu_pre = 4.0_wp
      hv_pre = 0.0_wp
      h_layer(:, :, 1) = 1.0_wp
      h_layer(:, :, 2) = 1.0_wp
      hu_layer(:, :, 1) = 1.0_wp  ! h_k * u_k = 1.0 * 1.0
      hu_layer(:, :, 2) = 3.0_wp  ! h_k * u_k = 1.0 * 3.0
      hv_layer = 0.0_wp

      ! After barotropic update: HU accelerated to 8 (U=4, du=2)
      h = 2.0_wp
      hu = 8.0_wp
      hv = 0.0_wp

      !$acc data copyin( h, hu, hv, hu_pre, hv_pre, dsig) &
      !$acc&            copy( h_layer, hu_layer, hv_layer)
      call ensure_grid()
      call ml_distribute_barotropic(h, hu, hv, hu_pre, hv_pre, dsig, &
                                    h_layer, hu_layer, hv_layer, &
                                    test_grid, NZ)
      !$acc end data

      ! du = (8 - 4) / 2 = 2
      ! u1_new = 1.0 + 2.0 = 3.0, hu1_new = 1.0 * 3.0 = 3.0
      ! u2_new = 3.0 + 2.0 = 5.0, hu2_new = 1.0 * 5.0 = 5.0
      u1_new = 3.0_wp
      u2_new = 5.0_wp
      call check(error, abs(hu_layer(5, 5, 1) - u1_new) < 1.0e-12_wp, &
                 "bottom layer: u_old=1 + du=2 => hu_k = 3")
      if (allocated(error)) return
      call check(error, abs(hu_layer(5, 5, 2) - u2_new) < 1.0e-12_wp, &
                 "top layer: u_old=3 + du=2 => hu_k = 5")
      if (allocated(error)) return
      ! Shear preserved: u2 - u1 = 5 - 3 = 2 (same as before: 3 - 1 = 2)
      call check(error, abs((hu_layer(5, 5, 2) - hu_layer(5, 5, 1)) - 2.0_wp) < 1.0e-12_wp, &
                 "velocity shear preserved: du2 - du1 = 2")
   end subroutine test_distribute_shear

   subroutine test_w_uniform(error)
      !! Uniform horizontal velocity with no divergence => w = 0 everywhere
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), dsig(NZ)

      h_layer = 1.0_wp
      hu_layer = 2.0_wp
      hv_layer = 0.0_wp
      w = 999.0_wp
      dsig = 1.0_wp/real(NZ, wp)

      !$acc data copyin( h_layer, hu_layer, hv_layer, dsig) copy( w)
      call ensure_grid()
      call ml_compute_vertical_velocity(h_layer, hu_layer, hv_layer, &
                                        dsig, w, test_grid, NZ)
      !$acc end data

      call check(error, abs(w(5, 5, 1)) < 1.0e-12_wp, "w(bottom) = 0")
      if (allocated(error)) return
      call check(error, abs(w(5, 5, 2)) < 1.0e-12_wp, &
                 "w(interface) = 0 for uniform flow")
      if (allocated(error)) return
      call check(error, abs(w(5, 5, 3)) < 1.0e-12_wp, &
                 "w(top) = 0 for uniform flow")
   end subroutine test_w_uniform

   subroutine test_vexch_no_w(error)
      !! Zero vertical velocity => no momentum exchange
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1)
      real(wp) :: hu_before

      h_layer = 1.0_wp
      hu_layer(:, :, 1) = 2.0_wp
      hu_layer(:, :, 2) = 4.0_wp
      hv_layer = 0.0_wp
      w = 0.0_wp

      hu_before = hu_layer(5, 5, 1)

      !$acc data copyin( h_layer, w) copy( hu_layer, hv_layer)
      call ensure_grid()
      call ml_compute_vertical_exchange(h_layer, hu_layer, hv_layer, &
                                        w, test_grid, NZ, 0.1_wp)
      !$acc end data

      call check(error, abs(hu_layer(5, 5, 1) - hu_before) < 1.0e-12_wp, &
                 "no exchange when w = 0")
   end subroutine test_vexch_no_w

   subroutine test_vexch_transfers(error)
      !! Positive w at interface transfers momentum upward (upwind)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1)
      real(wp) :: dt, flux_expected
      real(wp) :: hu1_before, hu2_before

      dt = 0.1_wp
      h_layer = 1.0_wp
      hu_layer(:, :, 1) = 2.0_wp  ! u1 = 2
      hu_layer(:, :, 2) = 4.0_wp  ! u2 = 4
      hv_layer = 0.0_wp
      w = 0.0_wp
      w(:, :, 2) = 0.5_wp  ! upward at interface between layer 1 and 2

      hu1_before = hu_layer(5, 5, 1)
      hu2_before = hu_layer(5, 5, 2)

      ! w > 0 => upwind from below: flux = w * u1 * dt = 0.5 * 2.0 * 0.1 = 0.1
      flux_expected = 0.5_wp*2.0_wp*dt

      !$acc data copyin( h_layer, w) copy( hu_layer, hv_layer)
      call ensure_grid()
      call ml_compute_vertical_exchange(h_layer, hu_layer, hv_layer, &
                                        w, test_grid, NZ, dt)
      !$acc end data

      ! Bottom layer loses momentum
      call check(error, abs(hu_layer(5, 5, 1) - (hu1_before - flux_expected)) < 1.0e-12_wp, &
                 "bottom layer loses momentum upward")
      if (allocated(error)) return
      ! Top layer gains momentum
      call check(error, abs(hu_layer(5, 5, 2) - (hu2_before + flux_expected)) < 1.0e-12_wp, &
                 "top layer gains momentum from below")
   end subroutine test_vexch_transfers

   ! ================================================================
   ! Baroclinic pressure gradient tests
   ! ================================================================

   subroutine test_bpg_uniform(error)
      !! Uniform density: BPG should be zero (no density-driven flow)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), rho_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: rho_0, dt

      rho_0 = 1000.0_wp
      dt = 0.1_wp
      h_layer = 1.0_wp
      rho_layer = 1025.0_wp  ! uniform density everywhere
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp

      !$acc data copyin( h_layer, rho_layer) copy( hu_layer, hv_layer)
      call ensure_grid()
      call ml_baroclinic_pressure_gradient(h_layer, rho_layer, &
                                           hu_layer, hv_layer, rho_0, &
                                           TEST_RHO_REF, TEST_RHO_REF_DZ, TEST_NZ_REF, &
                                           test_grid, NZ, dt)
      !$acc end data

      ! No density gradient => no BPG => momentum unchanged
      call check(error, abs(hu_layer(5, 5, 1)) < 1.0e-12_wp, &
                 "BPG zero for uniform density (hu)")
      if (allocated(error)) return
      call check(error, abs(hv_layer(5, 5, 1)) < 1.0e-12_wp, &
                 "BPG zero for uniform density (hv)")
   end subroutine test_bpg_uniform

   subroutine test_bpg_lock_exchange(error)
      !! Bottom-layer density gradient with uniform top layer.
      !! Dense bottom water on the left should spread right:
      !! the internal pressure at the bottom is higher on the left
      !! (heavier bottom), so the BPG pushes the bottom layer right.
      !! The top layer has uniform density so sees zero BPG.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), rho_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: rho_0, dt
      integer :: i, j

      rho_0 = 1000.0_wp
      dt = 0.01_wp
      h_layer = 1.0_wp
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp

      ! Top layer: uniform density everywhere
      rho_layer(:, :, 2) = 1000.0_wp

      ! Bottom layer: heavy on left, light on right
      do j = 1, NY
         do i = 1, NX
            if (i <= NX/2) then
               rho_layer(i, j, 1) = 1030.0_wp
            else
               rho_layer(i, j, 1) = 1020.0_wp
            end if
         end do
      end do

      !$acc data copyin( h_layer, rho_layer) copy( hu_layer, hv_layer)
      call ensure_grid()
      call ml_baroclinic_pressure_gradient(h_layer, rho_layer, &
                                           hu_layer, hv_layer, rho_0, &
                                           TEST_RHO_REF, TEST_RHO_REF_DZ, TEST_NZ_REF, &
                                           test_grid, NZ, dt)
      !$acc end data

      ! Internal pressure at bottom layer centre:
      !   p_1 = rho_2*h_2 + 0.5*rho_1*h_1
      !   Left:  1000 + 0.5*1030 = 1515
      !   Right: 1000 + 0.5*1020 = 1510
      ! dp/dx < 0 (higher pressure on left) => force pushes RIGHT
      ! Dense bottom water spreads to the right
      call check(error, hu_layer(6, 5, 1) > 0.0_wp, &
                 "dense bottom water pushed right (toward lighter side)")
      if (allocated(error)) return

      ! Top layer: uniform density and uniform h => dp_top/dx = 0
      ! No force on the top layer
      call check(error, abs(hu_layer(6, 5, 2)) < 1.0e-12_wp, &
                 "uniform top layer has zero BPG")
   end subroutine test_bpg_lock_exchange

   ! ================================================================
   ! Salinity and equation of state tests
   ! ================================================================

   subroutine test_init_salinity(error)
      !! hS = h_k * S_init
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hS_layer(NX, NY, NZ)
      real(wp) :: S_init

      h_layer = 1.5_wp
      S_init = 35.0_wp

      !$acc data copyin( h_layer) copy( hS_layer)
      call ensure_grid()
      call ml_init_tracer(h_layer, hS_layer, S_init, test_grid, NZ)
      !$acc end data

      call check(error, abs(hS_layer(5, 5, 1) - 1.5_wp*35.0_wp) < 1.0e-12_wp, &
                 "hS = h * S_init = 1.5 * 35 = 52.5")
   end subroutine test_init_salinity

   subroutine test_eos_linear(error)
      !! rho = rho_0 + beta_S * S
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hS_layer(NX, NY, NZ), hT_layer(NX, NY, NZ)
      real(wp) :: rho_layer(NX, NY, NZ)
      real(wp) :: rho_0, beta_S, S_ref, alpha_T, T_ref, expected

      rho_0 = 1000.0_wp
      beta_S = 0.78_wp
      S_ref = 0.0_wp
      alpha_T = 0.0_wp          ! decouple T: passive tracer, legacy form
      T_ref = 15.0_wp
      h_layer = 1.0_wp
      hS_layer(:, :, 1) = 1.0_wp*35.0_wp  ! S = 35 PSU (seawater)
      hS_layer(:, :, 2) = 1.0_wp*0.0_wp   ! S = 0 PSU (fresh water)
      hT_layer = 1.0_wp*T_ref             ! T = T_ref so T contribution = 0

      !$acc data copyin( h_layer, hS_layer, hT_layer) copy( rho_layer)
      call ensure_grid()
      call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_layer, &
                                rho_0, beta_S, S_ref, alpha_T, T_ref, &
                                test_grid, NZ)
      !$acc end data

      expected = 1000.0_wp + 0.78_wp*35.0_wp  ! = 1027.3
      call check(error, abs(rho_layer(5, 5, 1) - expected) < 1.0e-10_wp, &
                 "rho(S=35) = 1000 + 0.78*35 = 1027.3")
      if (allocated(error)) return
      call check(error, abs(rho_layer(5, 5, 2) - 1000.0_wp) < 1.0e-10_wp, &
                 "rho(S=0) = 1000.0 (fresh water)")
   end subroutine test_eos_linear

   subroutine test_init_tracer_temperature(error)
      !! ml_init_tracer also initialises temperature: hT = h * T_init
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hT_layer(NX, NY, NZ)
      real(wp), parameter :: T_INIT = 10.0_wp

      h_layer = 1.5_wp

      !$acc data copyin( h_layer) copyout( hT_layer)
      call ensure_grid()
      call ml_init_tracer(h_layer, hT_layer, T_INIT, test_grid, NZ)
      !$acc end data

      call check(error, abs(hT_layer(5, 5, 1) - 1.5_wp*T_INIT) < 1.0e-12_wp, &
                 "hT = h * T_init = 1.5 * 10 = 15")
      if (allocated(error)) return
      call check(error, abs(hT_layer(3, 3, 2) - 1.5_wp*T_INIT) < 1.0e-12_wp, &
                 "hT = h * T_init layer 2")
   end subroutine test_init_tracer_temperature

   subroutine test_init_tracer_stratified(error)
      !! ml_init_tracer_stratified seeds a linear T(z): T_BOT at the bed
      !! (k=1) -> T_SURF at the surface (k=nz), hT_k = h_k * T(k).  Mirrors
      !! the unstructured `tracer_init_stratified_linear_profile` analytical
      !! test for matched-pair coastal parity.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hT_layer(NX, NY, NZ)
      real(wp), parameter :: T_BOT = 12.0_wp, T_SURF = 22.0_wp
      real(wp) :: t_expect
      integer :: k

      h_layer = 1.5_wp

      !$acc data copyin( h_layer) copyout( hT_layer)
      call ensure_grid()
      call ml_init_tracer_stratified(h_layer, hT_layer, T_BOT, T_SURF, test_grid, NZ)
      !$acc end data

      do k = 1, NZ
         t_expect = T_BOT + (T_SURF - T_BOT)*real(k - 1, wp)/real(NZ - 1, wp)
         call check(error, abs(hT_layer(5, 5, k) - 1.5_wp*t_expect) < 1.0e-12_wp, &
                    "hT = h * linear T(z)")
         if (allocated(error)) return
      end do
   end subroutine test_init_tracer_stratified

   subroutine test_eos_two_tracer(error)
      !! Full two-tracer EOS: rho = rho_0 + beta_S*(S-S_ref) - alpha_T*(T-T_ref)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hS_layer(NX, NY, NZ), hT_layer(NX, NY, NZ)
      real(wp) :: rho_layer(NX, NY, NZ)
      real(wp), parameter :: rho_0 = 1000.0_wp
      real(wp), parameter :: beta_S = 0.78_wp
      real(wp), parameter :: S_ref = 0.0_wp
      real(wp), parameter :: alpha_T = 0.2_wp
      real(wp), parameter :: T_ref = 15.0_wp
      real(wp) :: expected

      h_layer = 1.0_wp
      hS_layer(:, :, 1) = 35.0_wp        ! S = 35, warm: 20
      hT_layer(:, :, 1) = 20.0_wp        !   rho = 1000 + 0.78*35 - 0.2*(20-15) = 1026.3
      hS_layer(:, :, 2) = 0.0_wp        ! S = 0,  cold: 5
      hT_layer(:, :, 2) = 5.0_wp        !   rho = 1000 + 0        - 0.2*(5-15)  = 1002.0

      !$acc data copyin( h_layer, hS_layer, hT_layer) copy( rho_layer)
      call ensure_grid()
      call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_layer, &
                                rho_0, beta_S, S_ref, alpha_T, T_ref, &
                                test_grid, NZ)
      !$acc end data

      expected = rho_0 + beta_S*(35.0_wp - S_ref) - alpha_T*(20.0_wp - T_ref)
      call check(error, abs(rho_layer(5, 5, 1) - expected) < 1.0e-10_wp, &
                 "rho(S=35,T=20) matches two-tracer formula")
      if (allocated(error)) return
      expected = rho_0 + beta_S*(0.0_wp - S_ref) - alpha_T*(5.0_wp - T_ref)
      call check(error, abs(rho_layer(5, 5, 2) - expected) < 1.0e-10_wp, &
                 "rho(S=0,T=5) matches two-tracer formula")
   end subroutine test_eos_two_tracer

   subroutine test_eos_alpha_T_zero_regression(error)
      !! With alpha_T=0 and S_ref=0, the two-tracer EOS reduces to the
      !! legacy single-tracer form: rho = rho_0 + beta_S * S.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hS_layer(NX, NY, NZ), hT_layer(NX, NY, NZ)
      real(wp) :: rho_layer(NX, NY, NZ)
      real(wp), parameter :: rho_0 = 1000.0_wp
      real(wp), parameter :: beta_S = 0.78_wp

      h_layer = 1.0_wp
      hS_layer(:, :, 1) = 35.0_wp
      hS_layer(:, :, 2) = 0.0_wp
      hT_layer = 1.0_wp*99.0_wp   ! arbitrary T — must NOT affect rho when alpha_T=0

      !$acc data copyin( h_layer, hS_layer, hT_layer) copy( rho_layer)
      call ensure_grid()
      call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_layer, &
                                rho_0, beta_S, 0.0_wp, 0.0_wp, 15.0_wp, &
                                test_grid, NZ)
      !$acc end data

      call check(error, abs(rho_layer(5, 5, 1) - (rho_0 + beta_S*35.0_wp)) < 1.0e-10_wp, &
                 "alpha_T=0 regression: rho = rho_0 + beta_S*S")
      if (allocated(error)) return
      call check(error, abs(rho_layer(5, 5, 2) - rho_0) < 1.0e-10_wp, &
                 "alpha_T=0 regression: rho = rho_0 at S=0 (T ignored)")
   end subroutine test_eos_alpha_T_zero_regression

   subroutine test_eos_S_ref_nonzero(error)
      !! With S=S_ref, density anomaly from salinity is exactly zero.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hS_layer(NX, NY, NZ), hT_layer(NX, NY, NZ)
      real(wp) :: rho_layer(NX, NY, NZ)

      h_layer = 1.0_wp
      hS_layer = 1.0_wp*35.0_wp          ! S = 35
      hT_layer = 1.0_wp*15.0_wp          ! T = T_ref

      !$acc data copyin( h_layer, hS_layer, hT_layer) copy( rho_layer)
      call ensure_grid()
      call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_layer, &
                                1000.0_wp, 0.78_wp, 35.0_wp, &
                                0.17_wp, 15.0_wp, &
                                test_grid, NZ)
      !$acc end data

      call check(error, abs(rho_layer(5, 5, 1) - 1000.0_wp) < 1.0e-10_wp, &
                 "S=S_ref and T=T_ref -> rho = rho_0 exactly")
   end subroutine test_eos_S_ref_nonzero

   subroutine test_clamp_tracer_both_tracers(error)
      !! Same ml_clamp_tracer kernel clamps S to [0,40] and T to [-2,40].
      !! Verifies the conservative limiter:
      !!   (a) post-clamp S/T values stay within [tr_min, tr_max] in
      !!       every layer, and
      !!   (b) per-cell column total tracer is preserved exactly when
      !!       the column-mean concentration is itself in bounds (the
      !!       redistribute step pushes the rejected excess back into
      !!       spare capacity in saturated layers).
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hS_layer(NX, NY, NZ), hT_layer(NX, NY, NZ)
      real(wp) :: hS_total_pre, hS_total_post, hT_total_pre, hT_total_post
      real(wp) :: S1, S2, T1, T2

      h_layer = 1.0_wp
      ! Salinity: overshoot high (50) and undershoot negative (-5)
      ! Column-mean S = (50 + -5) / 2 = 22.5, well inside [0, 40] →
      ! column-conservative redistribute should preserve the total.
      hS_layer(:, :, 1) = 50.0_wp
      hS_layer(:, :, 2) = -5.0_wp
      ! Temperature: overshoot high (45) and undershoot freezing (-10)
      ! Column-mean T = (45 + -10) / 2 = 17.5, inside [-2, 40].
      hT_layer(:, :, 1) = 45.0_wp
      hT_layer(:, :, 2) = -10.0_wp

      hS_total_pre = sum(hS_layer(5, 5, :))
      hT_total_pre = sum(hT_layer(5, 5, :))

      !$acc data copyin( h_layer) copy( hS_layer, hT_layer)
      call ensure_grid()
      call ml_clamp_tracer(hS_layer, h_layer, 0.0_wp, 40.0_wp, test_grid, NZ)
      call ml_clamp_tracer(hT_layer, h_layer, -2.0_wp, 40.0_wp, test_grid, NZ)
      !$acc end data

      S1 = hS_layer(5, 5, 1)/h_layer(5, 5, 1)
      S2 = hS_layer(5, 5, 2)/h_layer(5, 5, 2)
      T1 = hT_layer(5, 5, 1)/h_layer(5, 5, 1)
      T2 = hT_layer(5, 5, 2)/h_layer(5, 5, 2)
      hS_total_post = sum(hS_layer(5, 5, :))
      hT_total_post = sum(hT_layer(5, 5, :))

      ! Bounds: every layer must end up within [tr_min, tr_max].
      call check(error, S1 <= 40.0_wp + 1.0e-12_wp .and. S1 >= 0.0_wp - 1.0e-12_wp, &
                 "S layer 1 within [0, 40] after clamp")
      if (allocated(error)) return
      call check(error, S2 <= 40.0_wp + 1.0e-12_wp .and. S2 >= 0.0_wp - 1.0e-12_wp, &
                 "S layer 2 within [0, 40] after clamp")
      if (allocated(error)) return
      call check(error, T1 <= 40.0_wp + 1.0e-12_wp .and. T1 >= -2.0_wp - 1.0e-12_wp, &
                 "T layer 1 within [-2, 40] after clamp")
      if (allocated(error)) return
      call check(error, T2 <= 40.0_wp + 1.0e-12_wp .and. T2 >= -2.0_wp - 1.0e-12_wp, &
                 "T layer 2 within [-2, 40] after clamp")
      if (allocated(error)) return

      ! Column-conservative: total tracer preserved when column-mean in bounds.
      call check(error, abs(hS_total_post - hS_total_pre) < 1.0e-12_wp, &
                 "S column total preserved by conservative limiter")
      if (allocated(error)) return
      call check(error, abs(hT_total_post - hT_total_pre) < 1.0e-12_wp, &
                 "T column total preserved by conservative limiter")
   end subroutine test_clamp_tracer_both_tracers

   subroutine test_distribute_salinity(error)
      !! Distributing salinity through barotropic step preserves S value
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hS_layer(NX, NY, NZ)
      real(wp) :: dsig(NZ), h_new(NX, NY)
      real(wp) :: S_before, S_after

      dsig = 1.0_wp/real(NZ, wp)

      ! Before: h_k = 1.0, S = 35.0 => hS = 35.0
      h_layer = 1.0_wp
      hS_layer = 35.0_wp  ! hS = h * S = 1.0 * 35.0

      ! After barotropic update: total depth changed from 2.0 to 3.0
      h_new = 3.0_wp

      S_before = hS_layer(5, 5, 1)/h_layer(5, 5, 1)  ! = 35.0

      !$acc data copyin( h_layer, dsig, h_new) copy( hS_layer)
      call ensure_grid()
      call ml_distribute_tracer(h_layer, hS_layer, dsig, h_new, &
                                test_grid, NZ)
      !$acc end data

      ! h_k_new = dsig * H_new = 0.5 * 3.0 = 1.5
      ! S should be preserved: hS_new / h_k_new = S_before
      S_after = hS_layer(5, 5, 1)/(dsig(1)*h_new(5, 5))
      call check(error, abs(S_after - S_before) < 1.0e-12_wp, &
                 "salinity preserved through barotropic step")
      if (allocated(error)) return
      ! hS_new = 1.5 * 35.0 = 52.5
      call check(error, abs(hS_layer(5, 5, 1) - 52.5_wp) < 1.0e-10_wp, &
                 "hS = h_new * S = 1.5 * 35 = 52.5")
   end subroutine test_distribute_salinity

   subroutine test_vexch_salinity(error)
      !! Vertical exchange transfers salinity upward with positive w
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hS_layer(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1)
      real(wp) :: dt, hS1_before, hS2_before, flux_expected

      dt = 0.1_wp
      h_layer = 1.0_wp
      hS_layer(:, :, 1) = 35.0_wp   ! S1 = 35 (salty bottom)
      hS_layer(:, :, 2) = 0.0_wp    ! S2 = 0 (fresh top)
      w = 0.0_wp
      w(:, :, 2) = 0.5_wp  ! upward at interface

      hS1_before = hS_layer(5, 5, 1)
      hS2_before = hS_layer(5, 5, 2)
      ! w > 0 => upwind from below: flux = w * S1 * dt = 0.5 * 35.0 * 0.1
      flux_expected = 0.5_wp*35.0_wp*dt

      !$acc data copyin( h_layer, w) copy( hS_layer)
      call ensure_grid()
      call ml_vertical_exchange_tracer(h_layer, hS_layer, w, &
                                       test_grid, NZ, dt)
      !$acc end data

      call check(error, abs(hS_layer(5, 5, 1) - (hS1_before - flux_expected)) < 1.0e-10_wp, &
                 "bottom layer loses salinity upward")
      if (allocated(error)) return
      call check(error, abs(hS_layer(5, 5, 2) - (hS2_before + flux_expected)) < 1.0e-10_wp, &
                 "top layer gains salinity from below")
   end subroutine test_vexch_salinity

   ! ================================================================
   ! Ri-dependent vertical mixing tests
   ! ================================================================

   subroutine test_vmix_zero_shear(error)
      !! Uniform velocity across layers: no shear => zero diffusive flux
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: rho_layer(NX, NY, NZ)
      real(wp) :: rho_0, nu_max, nu_bg, ri_alpha, dt
      real(wp) :: hu_before

      rho_0 = 1000.0_wp
      nu_max = 0.01_wp
      nu_bg = 1.0e-4_wp
      ri_alpha = 5.0_wp
      dt = 0.1_wp

      h_layer = 1.0_wp
      hu_layer = 2.0_wp   ! uniform velocity u=2 in both layers
      hv_layer = 0.0_wp
      rho_layer = 1025.0_wp

      hu_before = hu_layer(5, 5, 1)

      !$acc data copyin( h_layer, rho_layer) copy( hu_layer, hv_layer)
      call ensure_grid()
      call ml_vertical_viscosity(h_layer, hu_layer, hv_layer, rho_layer, &
                                 rho_0, nu_max, nu_bg, ri_alpha, &
                                 test_grid, NZ, dt)
      !$acc end data

      ! No shear => du=0 => zero flux regardless of Ri
      call check(error, abs(hu_layer(5, 5, 1) - hu_before) < 1.0e-14_wp, &
                 "zero shear gives zero mixing flux")
   end subroutine test_vmix_zero_shear

   subroutine test_vmix_unstratified_shear(error)
      !! Uniform density + velocity shear: Ri=0 => nu_t = nu_max + nu_bg.
      !! Analytic flux between layers k and k+1 (with h_k = h_{k+1} = 1):
      !!   flux_u = nu_t * du/h_avg * dt
      !!   hu(k)    += h_k    * flux_u
      !!   hu(k+1)  -= h_{k+1}* flux_u
      !! so the update magnitudes must match to machine precision, and
      !! momentum transfer is exactly conserved between the two layers.
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: RHO0_VAL = 1000.0_wp, NU_MAX = 0.01_wp
      real(wp), parameter :: NU_BG = 1.0e-4_wp, RI_ALPHA = 5.0_wp
      real(wp), parameter :: DT_VAL = 0.1_wp, H_VAL = 1.0_wp
      real(wp), parameter :: U_TOP = 4.0_wp, U_BOT = 0.0_wp
      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: rho_layer(NX, NY, NZ)
      real(wp) :: nu_t, du, flux_u, hu_bot_expected, hu_top_expected

      h_layer = H_VAL
      hu_layer(:, :, 1) = U_BOT
      hu_layer(:, :, 2) = U_TOP
      hv_layer = 0.0_wp
      rho_layer = RHO0_VAL    ! uniform => Ri = 0

      !$acc data copyin( h_layer, rho_layer) copy( hu_layer, hv_layer)
      call ensure_grid()
      call ml_vertical_viscosity(h_layer, hu_layer, hv_layer, rho_layer, &
                                 RHO0_VAL, NU_MAX, NU_BG, RI_ALPHA, &
                                 test_grid, NZ, DT_VAL)
      !$acc end data

      du = U_TOP - U_BOT
      nu_t = NU_MAX + NU_BG            ! Ri = 0 => denom = 1
      flux_u = nu_t*du/H_VAL*DT_VAL    ! h_avg = H_VAL

      hu_bot_expected = U_BOT + H_VAL*flux_u
      hu_top_expected = U_TOP - H_VAL*flux_u

      call check(error, abs(hu_layer(5, 5, 1) - hu_bot_expected) < 1.0e-12_wp, &
                 "bottom layer hu must match analytic PP flux update")
      if (allocated(error)) return
      call check(error, abs(hu_layer(5, 5, 2) - hu_top_expected) < 1.0e-12_wp, &
                 "top layer hu must match analytic PP flux update")
      if (allocated(error)) return
      ! Column-integrated momentum conserved (flux-form diffusion).
      call check(error, abs((hu_layer(5, 5, 1) + hu_layer(5, 5, 2)) - (U_BOT + U_TOP)) &
                 < 1.0e-12_wp, &
                 "PP vertical mixing conserves column-integrated hu")
   end subroutine test_vmix_unstratified_shear

   subroutine test_vmix_strong_strat(error)
      !! Strong stable stratification + weak shear: Ri >> 1 => nu_t ~ nu_bg.
      !! Compare flux magnitude against the unstratified case — should be
      !! much smaller (suppressed by Ri).
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: hu_unstrat(NX, NY, NZ), hu_strat(NX, NY, NZ)
      real(wp) :: rho_layer(NX, NY, NZ)
      real(wp) :: rho_0, nu_max, nu_bg, ri_alpha, dt
      real(wp) :: change_unstrat, change_strat

      rho_0 = 1000.0_wp
      nu_max = 0.01_wp
      nu_bg = 1.0e-4_wp
      ri_alpha = 5.0_wp
      dt = 0.1_wp

      h_layer = 1.0_wp
      hu_layer(:, :, 1) = 0.0_wp
      hu_layer(:, :, 2) = 1.0_wp  ! mild shear
      hv_layer = 0.0_wp

      ! Unstratified case
      rho_layer = 1000.0_wp

      !$acc data copyin( h_layer, rho_layer) copy( hu_layer, hv_layer)
      call ensure_grid()
      call ml_vertical_viscosity(h_layer, hu_layer, hv_layer, rho_layer, &
                                 rho_0, nu_max, nu_bg, ri_alpha, &
                                 test_grid, NZ, dt)
      !$acc end data
      hu_unstrat = hu_layer

      ! Reset and do strongly stratified case
      hu_layer(:, :, 1) = 0.0_wp
      hu_layer(:, :, 2) = 1.0_wp
      hv_layer = 0.0_wp

      ! Stable stratification: dense at bottom, light at top
      rho_layer(:, :, 1) = 1030.0_wp
      rho_layer(:, :, 2) = 1000.0_wp

      !$acc data copyin( h_layer, rho_layer) copy( hu_layer, hv_layer)
      call ml_vertical_viscosity(h_layer, hu_layer, hv_layer, rho_layer, &
                                 rho_0, nu_max, nu_bg, ri_alpha, &
                                 test_grid, NZ, dt)
      !$acc end data
      hu_strat = hu_layer

      change_unstrat = abs(hu_unstrat(5, 5, 1))
      change_strat = abs(hu_strat(5, 5, 1))

      ! Stable stratification with large Ri should mix much less
      call check(error, change_strat < 0.5_wp*change_unstrat, &
                 "strong stratification suppresses mixing")
   end subroutine test_vmix_strong_strat

   subroutine test_vmix_salinity(error)
      !! Salinity gradient with low Ri (unstratified): salinity difference decreases
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: hS_layer(NX, NY, NZ), rho_layer(NX, NY, NZ)
      real(wp) :: rho_0, nu_max, kappa_bg, ri_alpha, dt
      real(wp) :: S1_before, S2_before, S1_after, S2_after

      rho_0 = 1000.0_wp
      nu_max = 0.01_wp
      kappa_bg = 1.0e-5_wp
      ri_alpha = 5.0_wp
      dt = 0.1_wp

      h_layer = 1.0_wp
      hu_layer(:, :, 1) = 0.0_wp
      hu_layer(:, :, 2) = 1.0_wp  ! some shear to keep Ri low
      hv_layer = 0.0_wp
      rho_layer = 1000.0_wp  ! uniform density => Ri = 0
      ! Bottom salty, top fresh
      hS_layer(:, :, 1) = 1.0_wp*35.0_wp   ! h*S = 35
      hS_layer(:, :, 2) = 1.0_wp*0.0_wp    ! h*S = 0

      S1_before = hS_layer(5, 5, 1)/h_layer(5, 5, 1)
      S2_before = hS_layer(5, 5, 2)/h_layer(5, 5, 2)

      !$acc data copyin( h_layer, hu_layer, hv_layer, rho_layer) &
      !$acc&            copy( hS_layer)
      call ensure_grid()
      call ml_vertical_diffusion_tracer(h_layer, hu_layer, hv_layer, &
                                        hS_layer, rho_layer, rho_0, &
                                        nu_max, kappa_bg, ri_alpha, &
                                        test_grid, NZ, dt)
      !$acc end data

      S1_after = hS_layer(5, 5, 1)/h_layer(5, 5, 1)
      S2_after = hS_layer(5, 5, 2)/h_layer(5, 5, 2)

      ! Salt should mix upward: bottom loses salt, top gains salt
      call check(error, S1_after < S1_before, &
                 "bottom layer loses salinity via mixing")
      if (allocated(error)) return
      call check(error, S2_after > S2_before, &
                 "top layer gains salinity via mixing")
   end subroutine test_vmix_salinity

   ! ================================================================
   ! Pressure Jacobian BPG tests
   ! ================================================================

   subroutine test_jacobian_uniform(error)
      !! Uniform density: Jacobian BPG should be zero (same as sigma method)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), rho_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: rho_0, dt

      rho_0 = 1000.0_wp
      dt = 0.1_wp
      h_layer = 1.0_wp
      rho_layer = 1025.0_wp
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp

      !$acc data copyin( h_layer, rho_layer) copy( hu_layer, hv_layer)
      call ensure_grid()
      call ml_baroclinic_pressure_gradient_jacobian(h_layer, rho_layer, &
                                                    hu_layer, hv_layer, rho_0, &
                                                    TEST_RHO_REF, TEST_RHO_REF_DZ, TEST_NZ_REF, &
                                                    test_grid, NZ, dt)
      !$acc end data

      call check(error, abs(hu_layer(5, 5, 1)) < 1.0e-10_wp, &
                 "Jacobian BPG zero for uniform density (hu)")
      if (allocated(error)) return
      call check(error, abs(hv_layer(5, 5, 1)) < 1.0e-10_wp, &
                 "Jacobian BPG zero for uniform density (hv)")
   end subroutine test_jacobian_uniform

   subroutine test_jacobian_varying_depth(error)
      !! KEY TEST: uniform density but varying total depth.
      !! The sigma method has a PGE here (dp/dx != 0 at depth change).
      !! The Jacobian method should give exactly zero BPG because
      !! pressure at the same physical depth z is identical in all columns.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), rho_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: rho_0, dt
      integer :: i, j

      rho_0 = 1000.0_wp
      dt = 0.1_wp
      rho_layer = 1025.0_wp  ! uniform density
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp

      ! Varying depth: left half deep (h=2), right half shallow (h=1)
      do j = 1, NY
         do i = 1, NX
            if (i <= NX/2) then
               h_layer(i, j, :) = 2.0_wp/NZ  ! H = 2m
            else
               h_layer(i, j, :) = 1.0_wp/NZ  ! H = 1m
            end if
         end do
      end do

      !$acc data copyin( h_layer, rho_layer) copy( hu_layer, hv_layer)
      call ensure_grid()
      call ml_baroclinic_pressure_gradient_jacobian(h_layer, rho_layer, &
                                                    hu_layer, hv_layer, rho_0, &
                                                    TEST_RHO_REF, TEST_RHO_REF_DZ, TEST_NZ_REF, &
                                                    test_grid, NZ, dt)
      !$acc end data

      ! At a cell away from the depth transition, BPG should be zero
      ! (uniform density on both sides, same depth)
      call check(error, abs(hu_layer(4, 5, 1)) < 1.0e-10_wp, &
                 "Jacobian BPG zero in uniform-depth region (left)")
      if (allocated(error)) return

      ! At the depth transition, Jacobian BPG should STILL be zero
      ! because pressure at same z is the same with uniform rho.
      ! Cell 7 (shallow side) looks across at deep neighbour at i=6:
      ! its z_mid is well within the neighbour's column.
      call check(error, abs(hu_layer(7, 5, 1)) < 1.0e-10_wp, &
                 "Jacobian BPG zero at depth transition with uniform density")
      if (allocated(error)) return

      call check(error, abs(hv_layer(7, 5, 1)) < 1.0e-10_wp, &
                 "Jacobian BPG zero in y at depth transition")
      if (allocated(error)) return

      ! Cell 6 (deep side) looks across at shallow neighbour at i=7:
      ! its bottom-layer z_mid = 1.5m lies BELOW the neighbour's bed at 1m.
      ! With the shelf-break wall fallback (`p_nbr = p_here`) BPG is zero
      ! here too. Without the fix, the kernel clamps to the neighbour's
      ! bed pressure and produces a spurious gradient.
      call check(error, abs(hu_layer(6, 5, 1)) < 1.0e-10_wp, &
                 "Jacobian BPG zero on deep flank of shelf break (bottom layer)")
      if (allocated(error)) return

      call check(error, abs(hu_layer(6, 5, 2)) < 1.0e-10_wp, &
                 "Jacobian BPG zero on deep flank of shelf break (top layer)")
   end subroutine test_jacobian_varying_depth

   subroutine test_jacobian_lock_exchange(error)
      !! Lock exchange: Jacobian should give same qualitative result as sigma
      !! method — dense bottom water pushed rightward.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), rho_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: rho_0, dt
      integer :: i, j

      rho_0 = 1000.0_wp
      dt = 0.01_wp
      h_layer = 1.0_wp
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp

      ! Top layer: uniform
      rho_layer(:, :, 2) = 1000.0_wp
      ! Bottom layer: heavy left, light right
      do j = 1, NY
         do i = 1, NX
            if (i <= NX/2) then
               rho_layer(i, j, 1) = 1030.0_wp
            else
               rho_layer(i, j, 1) = 1020.0_wp
            end if
         end do
      end do

      !$acc data copyin( h_layer, rho_layer) copy( hu_layer, hv_layer)
      call ensure_grid()
      call ml_baroclinic_pressure_gradient_jacobian(h_layer, rho_layer, &
                                                    hu_layer, hv_layer, rho_0, &
                                                    TEST_RHO_REF, TEST_RHO_REF_DZ, TEST_NZ_REF, &
                                                    test_grid, NZ, dt)
      !$acc end data

      ! Dense bottom water should be pushed right (toward lighter side)
      call check(error, hu_layer(6, 5, 1) > 0.0_wp, &
                 "Jacobian: dense bottom water pushed right")
   end subroutine test_jacobian_lock_exchange

   ! ================================================================
   ! Smagorinsky horizontal viscosity tests
   ! ================================================================

   subroutine test_smag_uniform(error)
      !! Uniform velocity field: zero strain rate => zero viscosity => no change.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: cs, dt
      real(wp) :: hu_before

      cs = 0.15_wp
      dt = 0.1_wp
      h_layer = 1.0_wp
      hu_layer = 2.0_wp   ! uniform u = 2 m/s everywhere
      hv_layer = 0.5_wp   ! uniform v = 0.5 m/s everywhere

      hu_before = hu_layer(5, 5, 1)

      !$acc data copyin( h_layer) copy( hu_layer, hv_layer)
      call ensure_grid()
      call ml_horizontal_viscosity_smag(h_layer, hu_layer, hv_layer, &
                                        cs, test_grid, NZ, dt)
      !$acc end data

      call check(error, abs(hu_layer(5, 5, 1) - hu_before) < 1.0e-12_wp, &
                 "uniform velocity: no change from Smagorinsky")
   end subroutine test_smag_uniform

   subroutine test_smag_shear(error)
      !! Sharp velocity step u=0 -> u=U at i=NX/2. At the peak cell
      !! (i=NX/2+1) with h=1 and dx=dy=1:
      !!   du/dx = U/2, dv/dy = du/dy = dv/dx = 0
      !!   |S|   = sqrt((du/dx)^2 + 0) = U/2
      !!   nu_h  = (Cs*dx)^2 * |S| = Cs^2 * U / 2
      !!   lap u = (U - 2U + 0)/dx^2 + (U - 2U + U)/dy^2 = -U
      !! so the update is hu_new = U + dt * h * nu_h * (-U) = U * (1 - dt*Cs^2*U/2).
      !! The flat-region cell (i=NX/2+3) has no velocity difference and
      !! must be exactly untouched to machine precision.
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: CS_VAL = 0.2_wp, DT_VAL = 0.1_wp, U_STEP = 4.0_wp
      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: strain_expected, nu_h_expected, lap_u_expected, hu_peak_expected
      integer :: i, j

      h_layer = 1.0_wp
      hv_layer = 0.0_wp

      do j = 1, NY
         do i = 1, NX
            if (i <= NX/2) then
               hu_layer(i, j, :) = 0.0_wp
            else
               hu_layer(i, j, :) = U_STEP
            end if
         end do
      end do

      !$acc data copyin( h_layer) copy( hu_layer, hv_layer)
      call ensure_grid()
      call ml_horizontal_viscosity_smag(h_layer, hu_layer, hv_layer, &
                                        CS_VAL, test_grid, NZ, DT_VAL)
      !$acc end data

      strain_expected = U_STEP/(2.0_wp*DX)              ! |S|
      nu_h_expected = (CS_VAL*sqrt(DX*DY))**2*strain_expected
      lap_u_expected = -U_STEP/(DX*DX)
      hu_peak_expected = U_STEP + DT_VAL*1.0_wp*nu_h_expected*lap_u_expected

      call check(error, abs(hu_layer(NX/2 + 1, 5, 1) - hu_peak_expected) < 1.0e-12_wp, &
                 "Smagorinsky matches analytic update: U + dt*h*nu_h*lap(u)")
      if (allocated(error)) return

      ! Flat region: zero strain, zero laplacian, zero update.
      call check(error, abs(hu_layer(NX/2 + 3, 5, 1) - U_STEP) < 1.0e-12_wp, &
                 "Smagorinsky leaves smooth region exactly unchanged")
   end subroutine test_smag_shear

   ! ================================================================
   ! Per-layer bottom friction test
   ! ================================================================

   subroutine test_bottom_friction(error)
      !! Bottom friction damps the bottom layer (k=1), leaves the top
      !! layer (k=NZ) unchanged.  k=1=bottom under bottom-up ROMS indexing.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: manning_n, dt
      real(wp) :: hu_top_before, hu_bot_before

      manning_n = 0.025_wp
      dt = 0.1_wp
      h_layer = 1.0_wp
      ! Both layers moving at u = 2 m/s
      hu_layer = 2.0_wp
      hv_layer = 0.0_wp

      hu_top_before = hu_layer(5, 5, NZ)
      hu_bot_before = hu_layer(5, 5, 1)

      !$acc data copyin( h_layer) copy( hu_layer, hv_layer)
      call ensure_grid()
      call ml_apply_bottom_friction(h_layer, hu_layer, hv_layer, &
                                    manning_n, .false., 0.0_wp, test_grid, NZ, dt)
      !$acc end data

      ! Top layer (k=NZ) should be untouched
      call check(error, abs(hu_layer(5, 5, NZ) - hu_top_before) < 1.0e-12_wp, &
                 "bottom friction leaves top layer (k=NZ) unchanged")
      if (allocated(error)) return

      ! Bottom layer (k=1) should be decelerated
      call check(error, hu_layer(5, 5, 1) < hu_bot_before, &
                 "bottom friction decelerates bottom layer (k=1)")
      if (allocated(error)) return

      ! Bottom layer should still be positive (not over-damped)
      call check(error, hu_layer(5, 5, 1) > 0.0_wp, &
                 "bottom friction does not reverse flow")
   end subroutine test_bottom_friction

   subroutine test_surface_heat_flux(error)
      !! ml_apply_surface_tracer_fluxes warms hT(top) by dt * Q / (rho_0 * cp).
      !! With Q = +400 W/m^2, rho_0 = 1000, cp = 4000, dt = 1, h_top = 1:
      !!   d(hT) = 1 * 400 / (1000 * 4000) = 1e-4 [m * K]
      !! Top layer = k=NZ under the structured solver's bottom-up indexing
      !! (k=1=bottom, k=nz=surface), matching BPG / rho_ref / NH init.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: hS_layer(NX, NY, NZ), hT_layer(NX, NY, NZ)
      real(wp) :: heat_flux_2d(NX, NY), evap_minus_precip_2d(NX, NY)
      real(wp), parameter :: rho_0 = 1000.0_wp
      real(wp), parameter :: Q = 400.0_wp
      real(wp), parameter :: dt = 1.0_wp
      real(wp) :: hT_before, expected_delta

      h_layer = 1.0_wp
      hS_layer = 35.0_wp
      hT_layer = 15.0_wp*h_layer   ! hT = h * 15 = 15
      heat_flux_2d = Q
      evap_minus_precip_2d = 0.0_wp    ! isolate heat from salt
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp

      hT_before = hT_layer(5, 5, NZ)
      expected_delta = dt*Q/(rho_0*CP_WATER)

      !$acc data copyin( h_layer, hu_layer, hv_layer, heat_flux_2d, evap_minus_precip_2d) &
      !$acc&     copy( hS_layer, hT_layer)
      call ensure_grid()
      call ml_apply_surface_tracer_fluxes(hS_layer, hT_layer, h_layer, &
                                          hu_layer, hv_layer, &
                                          heat_flux_2d, evap_minus_precip_2d, &
                                          .false., .false., &
                                          rho_0, CP_WATER, test_grid, NZ, dt)
      !$acc end data

      call check(error, abs((hT_layer(5, 5, NZ) - hT_before) - expected_delta) < 1.0e-12_wp, &
                 "surface heat flux warms top layer (k=NZ) by dt*Q/(rho_0*cp)")
      if (allocated(error)) return
      ! Bottom layer (k=1) unchanged
      call check(error, abs(hT_layer(5, 5, 1) - 15.0_wp) < 1.0e-12_wp, &
                 "surface heat flux does not touch the bottom layer (k=1)")
      if (allocated(error)) return
      ! Salinity untouched when E-P = 0
      call check(error, abs(hS_layer(5, 5, NZ) - 35.0_wp) < 1.0e-12_wp, &
                 "zero E-P leaves hS unchanged")
   end subroutine test_surface_heat_flux

   subroutine test_surface_salt_flux(error)
      !! Virtual salt flux: hS(NZ) += dt * (E-P) * S_top.
      !! With E-P = 1e-6 m/s, S_top = 35, dt = 1:  d(hS) = 3.5e-5 [m*PSU]
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: hS_layer(NX, NY, NZ), hT_layer(NX, NY, NZ)
      real(wp) :: heat_flux_2d(NX, NY), evap_minus_precip_2d(NX, NY)
      real(wp), parameter :: rho_0 = 1000.0_wp
      real(wp), parameter :: EP = 1.0e-6_wp
      real(wp), parameter :: dt = 1.0_wp
      real(wp) :: hS_before, S_top, expected_delta

      h_layer = 1.0_wp
      hS_layer = 35.0_wp*h_layer    ! hS = h*35 = 35
      hT_layer = 15.0_wp*h_layer
      heat_flux_2d = 0.0_wp             ! isolate salt from heat
      evap_minus_precip_2d = EP
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp

      S_top = hS_layer(5, 5, NZ)/h_layer(5, 5, NZ)
      hS_before = hS_layer(5, 5, NZ)
      expected_delta = dt*EP*S_top      ! = 1e-6 * 35 = 3.5e-5

      !$acc data copyin( h_layer, hu_layer, hv_layer, heat_flux_2d, evap_minus_precip_2d) &
      !$acc&     copy( hS_layer, hT_layer)
      call ensure_grid()
      call ml_apply_surface_tracer_fluxes(hS_layer, hT_layer, h_layer, &
                                          hu_layer, hv_layer, &
                                          heat_flux_2d, evap_minus_precip_2d, &
                                          .false., .false., &
                                          rho_0, CP_WATER, test_grid, NZ, dt)
      !$acc end data

      call check(error, abs((hS_layer(5, 5, NZ) - hS_before) - expected_delta) < 1.0e-12_wp, &
                 "virtual salt flux concentrates top layer (k=NZ): hS += dt*(E-P)*S_top")
      if (allocated(error)) return
      call check(error, abs(hS_layer(5, 5, 1) - 35.0_wp) < 1.0e-12_wp, &
                 "virtual salt flux does not touch the bottom layer (k=1)")
      if (allocated(error)) return
      ! Temperature untouched when Q = 0
      call check(error, abs(hT_layer(5, 5, NZ) - 15.0_wp) < 1.0e-12_wp, &
                 "zero heat flux leaves hT unchanged")
   end subroutine test_surface_salt_flux

   ! ================================================================
   ! Horizontal advection (ml_advect_tracer)
   ! ================================================================

   subroutine test_advect_zero_velocity_no_change(error)
      !! With hu = hv = 0, every upwind flux is zero and hTr must be
      !! bit-identical to the input on every cell and layer.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      real(wp) :: hTr_layer(NX, NY, NZ), hTr_before(NX, NY, NZ)
      real(wp), parameter :: dt = 0.05_wp
      integer :: i, j, k

      h_layer = 1.0_wp
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp
      mass_flux_x = 0.0_wp
      mass_flux_y = 0.0_wp

      ! Non-trivial tracer pattern so any kernel that accidentally
      ! wrote a constant would be caught.
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               hTr_layer(i, j, k) = real(i, wp) + 2.0_wp*real(j, wp) &
                                    + 0.5_wp*real(k, wp)
            end do
         end do
      end do
      hTr_before = hTr_layer

      call tracer_scratch_ensure(scratch, NX, NY, NZ, &
                                 want_ppm=.false., want_upwind=.true., want_hdiff=.false.)
      !$acc data copyin(h_layer, hu_layer, hv_layer, mass_flux_x, mass_flux_y) &
      !$acc&     copy(hTr_layer)
      call ensure_grid()
      call ml_advect_tracer(h_layer, hu_layer, hv_layer, &
                            mass_flux_x, mass_flux_y, &
                            hTr_layer, scratch%adv_hTr_snap, test_grid, NZ, dt)
      !$acc end data
      call tracer_scratch_destroy(scratch)

      do k = 1, NZ
         do j = NG + 1, NY - NG
            do i = NG + 1, NX - NG
               call check(error, &
                          abs(hTr_layer(i, j, k) - hTr_before(i, j, k)) < 1.0e-14_wp, &
                          "zero-velocity advection must leave hTr unchanged")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_advect_zero_velocity_no_change

   subroutine test_advect_uniform_tracer_no_change(error)
      !! A spatially uniform tracer has zero Laplacian and zero gradient,
      !! so the upwind flux difference at every cell is zero regardless
      !! of the velocity field. Concentration Tr = hTr/h must stay
      !! constant up to the round-off from (u_right - u_left) * Tr * dt
      !! cancellation, which here is exactly zero because u is uniform.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      real(wp) :: hTr_layer(NX, NY, NZ)
      real(wp), parameter :: dt = 0.05_wp, U0 = 0.3_wp, V0 = -0.2_wp
      real(wp), parameter :: TR0 = 35.0_wp
      integer :: i, j, k

      h_layer = 1.0_wp
      hu_layer = U0*h_layer  ! uniform hu = u * h
      hv_layer = V0*h_layer
      hTr_layer = TR0*h_layer  ! uniform concentration
      ! Barotropic mass flux equals sum_k hu_layer at the face. With
      ! uniform velocity and uniform thickness this is just NZ*U0.
      mass_flux_x = real(NZ, wp)*U0
      mass_flux_y = real(NZ, wp)*V0

      call tracer_scratch_ensure(scratch, NX, NY, NZ, &
                                 want_ppm=.false., want_upwind=.true., want_hdiff=.false.)
      !$acc data copyin(h_layer, hu_layer, hv_layer, mass_flux_x, mass_flux_y) &
      !$acc&     copy(hTr_layer)
      call ensure_grid()
      call ml_advect_tracer(h_layer, hu_layer, hv_layer, &
                            mass_flux_x, mass_flux_y, &
                            hTr_layer, scratch%adv_hTr_snap, test_grid, NZ, dt)
      !$acc end data
      call tracer_scratch_destroy(scratch)

      do k = 1, NZ
         do j = NG + 1, NY - NG
            do i = NG + 1, NX - NG
               call check(error, &
                          abs(hTr_layer(i, j, k) - TR0*h_layer(i, j, k)) < 1.0e-12_wp, &
                          "uniform tracer with uniform velocity must stay unchanged")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_advect_uniform_tracer_no_change

   subroutine test_advect_upwind_direction(error)
      !! Drop a single tracer patch at (ic, jc) with uniform hu > 0 and
      !! hv = 0. After one step the upwind scheme must move mass
      !! RIGHTWARD: the source cell must lose content to the cell
      !! immediately to its right, and nothing to its left.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      real(wp) :: hTr_layer(NX, NY, NZ)
      real(wp), parameter :: dt = 0.05_wp
      real(wp), parameter :: U0 = 0.5_wp, TR_VALUE = 10.0_wp
      real(wp) :: hTr_src_before, hTr_left_before
      integer :: ic, jc, kc

      h_layer = 1.0_wp
      hu_layer = U0*h_layer
      hv_layer = 0.0_wp
      hTr_layer = 0.0_wp
      mass_flux_x = real(NZ, wp)*U0
      mass_flux_y = 0.0_wp

      ic = NX/2
      jc = NY/2
      kc = 1
      hTr_layer(ic, jc, kc) = TR_VALUE*h_layer(ic, jc, kc)
      hTr_src_before = hTr_layer(ic, jc, kc)
      hTr_left_before = hTr_layer(ic - 1, jc, kc)

      call tracer_scratch_ensure(scratch, NX, NY, NZ, &
                                 want_ppm=.false., want_upwind=.true., want_hdiff=.false.)
      !$acc data copyin(h_layer, hu_layer, hv_layer, mass_flux_x, mass_flux_y) &
      !$acc&     copy(hTr_layer)
      call ensure_grid()
      call ml_advect_tracer(h_layer, hu_layer, hv_layer, &
                            mass_flux_x, mass_flux_y, &
                            hTr_layer, scratch%adv_hTr_snap, test_grid, NZ, dt)
      !$acc end data
      call tracer_scratch_destroy(scratch)

      call check(error, hTr_layer(ic, jc, kc) < hTr_src_before, &
                 "upwind advect: source cell must lose tracer with hu > 0")
      if (allocated(error)) return
      call check(error, hTr_layer(ic + 1, jc, kc) > 0.0_wp, &
                 "upwind advect: downstream cell must gain tracer with hu > 0")
      if (allocated(error)) return
      call check(error, abs(hTr_layer(ic - 1, jc, kc) - hTr_left_before) < 1.0e-12_wp, &
                 "upwind advect: upstream cell must not receive tracer with hu > 0")
   end subroutine test_advect_upwind_direction

   ! ================================================================
   ! Horizontal tracer diffusion (ml_horizontal_diffusion_tracer)
   ! ================================================================

   subroutine test_hdiff_zero_kappa_no_change(error)
      !! kappa = 0 is the explicit early-return in the kernel. hTr
      !! must be untouched bit-for-bit even if the tracer field has
      !! huge Laplacians.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hTr_layer(NX, NY, NZ)
      real(wp) :: hTr_before(NX, NY, NZ)
      real(wp), parameter :: dt = 0.05_wp
      integer :: i, j, k

      h_layer = 1.0_wp
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               hTr_layer(i, j, k) = 0.01_wp*real(i*i + j*j, wp)
            end do
         end do
      end do
      hTr_before = hTr_layer

      call tracer_scratch_ensure(scratch, NX, NY, NZ, &
                                 want_ppm=.false., want_upwind=.false., want_hdiff=.true.)
      !$acc data copyin(h_layer) copy(hTr_layer)
      call ensure_grid()
      call ml_horizontal_diffusion_tracer(h_layer, hTr_layer, 0.0_wp, &
                                          scratch%hdiff_hTr_pre, &
                                          test_grid, NZ, dt)
      !$acc end data
      call tracer_scratch_destroy(scratch)

      do k = 1, NZ
         do j = NG + 1, NY - NG
            do i = NG + 1, NX - NG
               call check(error, &
                          abs(hTr_layer(i, j, k) - hTr_before(i, j, k)) < 1.0e-14_wp, &
                          "kappa = 0 must be a bit-exact no-op")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_hdiff_zero_kappa_no_change

   subroutine test_hdiff_uniform_no_change(error)
      !! A spatially uniform tracer has zero Laplacian, so even with a
      !! large kappa the kernel must leave hTr untouched.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hTr_layer(NX, NY, NZ)
      real(wp), parameter :: dt = 0.05_wp, KAPPA = 1.0_wp, TR0 = 20.0_wp
      integer :: i, j, k

      h_layer = 1.0_wp
      hTr_layer = TR0*h_layer

      call tracer_scratch_ensure(scratch, NX, NY, NZ, &
                                 want_ppm=.false., want_upwind=.false., want_hdiff=.true.)
      !$acc data copyin(h_layer) copy(hTr_layer)
      call ensure_grid()
      call ml_horizontal_diffusion_tracer(h_layer, hTr_layer, KAPPA, &
                                          scratch%hdiff_hTr_pre, &
                                          test_grid, NZ, dt)
      !$acc end data
      call tracer_scratch_destroy(scratch)

      do k = 1, NZ
         do j = NG + 1, NY - NG
            do i = NG + 1, NX - NG
               call check(error, &
                          abs(hTr_layer(i, j, k) - TR0*h_layer(i, j, k)) < 1.0e-12_wp, &
                          "uniform tracer with kappa > 0 must stay unchanged")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_hdiff_uniform_no_change

   subroutine test_hdiff_step_diffuses_symmetrically(error)
      !! Start with an x-step in tracer (Tr = 30 for i <= NX/2, Tr = 40
      !! otherwise) and apply one diffusion step. The signed change
      !! either side of the step must be equal and opposite (1D
      !! Laplacian of a step is a dipole +/-; symmetric diffusion
      !! leaks mass from the high side into the low side in equal
      !! amounts per unit time). Total tracer mass across the two
      !! columns straddling the step must be unchanged.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hTr_layer(NX, NY, NZ)
      real(wp), parameter :: dt = 0.05_wp, KAPPA = 0.5_wp
      real(wp), parameter :: TR_LO = 30.0_wp, TR_HI = 40.0_wp
      real(wp) :: delta_lo, delta_hi
      integer :: i, j, mid, kc

      h_layer = 1.0_wp
      mid = NX/2
      do j = 1, NY
         do i = 1, NX
            if (i <= mid) then
               hTr_layer(i, j, :) = TR_LO*h_layer(i, j, :)
            else
               hTr_layer(i, j, :) = TR_HI*h_layer(i, j, :)
            end if
         end do
      end do

      call tracer_scratch_ensure(scratch, NX, NY, NZ, &
                                 want_ppm=.false., want_upwind=.false., want_hdiff=.true.)
      !$acc data copyin(h_layer) copy(hTr_layer)
      call ensure_grid()
      call ml_horizontal_diffusion_tracer(h_layer, hTr_layer, KAPPA, &
                                          scratch%hdiff_hTr_pre, &
                                          test_grid, NZ, dt)
      !$acc end data
      call tracer_scratch_destroy(scratch)

      kc = 1
      ! Low-side column at i = mid: diffusion must push it UP.
      delta_lo = hTr_layer(mid, NY/2, kc)/h_layer(mid, NY/2, kc) - TR_LO
      ! High-side column at i = mid + 1: diffusion must pull it DOWN.
      delta_hi = hTr_layer(mid + 1, NY/2, kc)/h_layer(mid + 1, NY/2, kc) - TR_HI

      call check(error, delta_lo > 0.0_wp, &
                 "diffusion: low-side column of a step must warm up")
      if (allocated(error)) return
      call check(error, delta_hi < 0.0_wp, &
                 "diffusion: high-side column of a step must cool down")
      if (allocated(error)) return
      call check(error, abs(delta_lo + delta_hi) < 1.0e-12_wp, &
                 "diffusion across a step is antisymmetric (mass-preserving)")
   end subroutine test_hdiff_step_diffuses_symmetrically

   ! ================================================================
   ! Velocity truncation (ml_truncate_velocity)
   ! ================================================================

   subroutine test_truncate_below_threshold_untouched(error)
      !! Velocities below 0.9 * cfl * dx/dt must be left bit-identical.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: hu_before(NX, NY, NZ), hv_before(NX, NY, NZ)
      real(wp), parameter :: dt = 0.1_wp, cfl_trunc = 0.9_wp
      real(wp) :: u_max, u_below, v_below
      integer :: i, j, k

      u_max = 0.9_wp*cfl_trunc*DX/dt   ! = 0.729 m/s with DX=1, dt=0.1, cfl=0.9
      u_below = 0.5_wp*u_max           ! well below the truncation bound
      v_below = 0.3_wp*u_max

      h_layer = 1.0_wp
      hu_layer = u_below*h_layer
      hv_layer = v_below*h_layer
      hu_before = hu_layer
      hv_before = hv_layer

      !$acc data copyin(h_layer) copy(hu_layer, hv_layer)
      call ensure_grid()
      call ml_truncate_velocity(h_layer, hu_layer, hv_layer, &
                                test_grid, NZ, dt, cfl_trunc)
      !$acc end data

      do k = 1, NZ
         do j = NG + 1, NY - NG
            do i = NG + 1, NX - NG
               call check(error, &
                          abs(hu_layer(i, j, k) - hu_before(i, j, k)) < 1.0e-14_wp, &
                          "truncate must be a no-op on in-bounds hu")
               if (allocated(error)) return
               call check(error, &
                          abs(hv_layer(i, j, k) - hv_before(i, j, k)) < 1.0e-14_wp, &
                          "truncate must be a no-op on in-bounds hv")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_truncate_below_threshold_untouched

   subroutine test_truncate_above_threshold_clipped(error)
      !! Velocities set intentionally above the truncation bound must
      !! be clipped EXACTLY to 0.9 * cfl * dx/dt (in absolute value).
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp), parameter :: dt = 0.1_wp, cfl_trunc = 0.9_wp
      real(wp) :: u_max, u_above, v_above, u_observed, v_observed
      integer :: ic, jc, kc

      u_max = 0.9_wp*cfl_trunc*DX/dt   ! = 0.729 m/s
      u_above = 5.0_wp*u_max           ! grossly above the bound
      v_above = -7.0_wp*u_max          ! negative to also check sign preservation

      h_layer = 1.0_wp
      hu_layer = u_above*h_layer
      hv_layer = v_above*h_layer

      !$acc data copyin(h_layer) copy(hu_layer, hv_layer)
      call ensure_grid()
      call ml_truncate_velocity(h_layer, hu_layer, hv_layer, &
                                test_grid, NZ, dt, cfl_trunc)
      !$acc end data

      ic = NX/2; jc = NY/2; kc = 1
      u_observed = hu_layer(ic, jc, kc)/h_layer(ic, jc, kc)
      v_observed = hv_layer(ic, jc, kc)/h_layer(ic, jc, kc)

      call check(error, abs(abs(u_observed) - u_max) < 1.0e-12_wp, &
                 "truncate: |u| must be clipped to 0.9 * cfl_trunc * dx/dt")
      if (allocated(error)) return
      call check(error, abs(abs(v_observed) - u_max) < 1.0e-12_wp, &
                 "truncate: |v| must be clipped to 0.9 * cfl_trunc * dy/dt")
      if (allocated(error)) return
      call check(error, u_observed > 0.0_wp, &
                 "truncate must preserve the sign of u")
      if (allocated(error)) return
      call check(error, v_observed < 0.0_wp, &
                 "truncate must preserve the sign of v")
   end subroutine test_truncate_above_threshold_clipped

   subroutine test_truncate_per_cell_independence(error)
      !! Truncate one column in the interior to an above-threshold
      !! value, leave the rest in-bounds. The kernel must clip the
      !! offending column and leave every other column bit-identical —
      !! no spurious cross-talk between cells.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: hu_before(NX, NY, NZ), hv_before(NX, NY, NZ)
      real(wp), parameter :: dt = 0.1_wp, cfl_trunc = 0.9_wp
      real(wp) :: u_max, u_bad, u_good
      integer :: i, j, k, ic, jc

      u_max = 0.9_wp*cfl_trunc*DX/dt
      u_good = 0.2_wp*u_max
      u_bad = 4.0_wp*u_max

      h_layer = 1.0_wp
      hu_layer = u_good*h_layer
      hv_layer = 0.0_wp

      ic = NX/2; jc = NY/2
      hu_layer(ic, jc, :) = u_bad*h_layer(ic, jc, :)

      hu_before = hu_layer
      hv_before = hv_layer

      !$acc data copyin(h_layer) copy(hu_layer, hv_layer)
      call ensure_grid()
      call ml_truncate_velocity(h_layer, hu_layer, hv_layer, &
                                test_grid, NZ, dt, cfl_trunc)
      !$acc end data

      ! All other columns must be unchanged.
      do k = 1, NZ
         do j = NG + 1, NY - NG
            do i = NG + 1, NX - NG
               if (i == ic .and. j == jc) cycle
               call check(error, &
                          abs(hu_layer(i, j, k) - hu_before(i, j, k)) < 1.0e-14_wp, &
                          "truncate must not touch in-bounds neighbour columns")
               if (allocated(error)) return
            end do
         end do
      end do

      ! The flagged column got clipped.
      do k = 1, NZ
         call check(error, &
                    abs(hu_layer(ic, jc, k)/h_layer(ic, jc, k) - u_max) < 1.0e-12_wp, &
                    "truncate must clip the offending column to u_max")
         if (allocated(error)) return
      end do
   end subroutine test_truncate_per_cell_independence

   ! ================================================================
   ! PPM tracer advection (ml_advect_tracer_ppm)
   ! ================================================================

   subroutine test_ppm_zero_velocity_no_change(error)
      !! With hu = hv = 0, the PPM flux is zero on every face and hTr
      !! must be bit-identical to the input on every cell and layer.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      real(wp) :: hTr_layer(NX, NY, NZ), hTr_before(NX, NY, NZ)
      real(wp), parameter :: dt = 0.05_wp
      integer :: i, j, k

      h_layer = 1.0_wp
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp
      mass_flux_x = 0.0_wp
      mass_flux_y = 0.0_wp

      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               hTr_layer(i, j, k) = real(i, wp) + 2.0_wp*real(j, wp) &
                                    + 0.5_wp*real(k, wp)
            end do
         end do
      end do
      hTr_before = hTr_layer

      call tracer_scratch_ensure(scratch, NX, NY, NZ, &
                                 want_ppm=.true., want_upwind=.false., want_hdiff=.false.)
      !$acc data copyin(h_layer, hu_layer, hv_layer, mass_flux_x, mass_flux_y) &
      !$acc&     copy(hTr_layer)
      call ensure_grid()
      call ml_advect_tracer_ppm(h_layer, hu_layer, hv_layer, &
                                mass_flux_x, mass_flux_y, &
                                hTr_layer, &
                                scratch%ppm_tr_snap, scratch%ppm_ql_x, scratch%ppm_qr_x, &
                                scratch%ppm_ql_y, scratch%ppm_qr_y, &
                                test_grid, NZ, dt)
      !$acc end data
      call tracer_scratch_destroy(scratch)

      do k = 1, NZ
         do j = NG + 1, NY - NG
            do i = NG + 1, NX - NG
               call check(error, &
                          abs(hTr_layer(i, j, k) - hTr_before(i, j, k)) < 1.0e-14_wp, &
                          "PPM zero-velocity advection must leave hTr unchanged")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_ppm_zero_velocity_no_change

   subroutine test_ppm_uniform_tracer_no_change(error)
      !! A spatially uniform tracer under PPM reconstruction produces
      !! constant edge values and zero parabolic correction, so the
      !! time-averaged flux at every face is TR0 * hu_face. The flux
      !! divergence vanishes for uniform hu, leaving hTr unchanged.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      real(wp) :: hTr_layer(NX, NY, NZ)
      real(wp), parameter :: dt = 0.05_wp, U0 = 0.3_wp, V0 = -0.2_wp
      real(wp), parameter :: TR0 = 35.0_wp
      integer :: i, j, k

      h_layer = 1.0_wp
      hu_layer = U0*h_layer
      hv_layer = V0*h_layer
      hTr_layer = TR0*h_layer
      mass_flux_x = real(NZ, wp)*U0
      mass_flux_y = real(NZ, wp)*V0

      call tracer_scratch_ensure(scratch, NX, NY, NZ, &
                                 want_ppm=.true., want_upwind=.false., want_hdiff=.false.)
      !$acc data copyin(h_layer, hu_layer, hv_layer, mass_flux_x, mass_flux_y) &
      !$acc&     copy(hTr_layer)
      call ensure_grid()
      call ml_advect_tracer_ppm(h_layer, hu_layer, hv_layer, &
                                mass_flux_x, mass_flux_y, &
                                hTr_layer, &
                                scratch%ppm_tr_snap, scratch%ppm_ql_x, scratch%ppm_qr_x, &
                                scratch%ppm_ql_y, scratch%ppm_qr_y, &
                                test_grid, NZ, dt)
      !$acc end data
      call tracer_scratch_destroy(scratch)

      do k = 1, NZ
         do j = NG + 1, NY - NG
            do i = NG + 1, NX - NG
               call check(error, &
                          abs(hTr_layer(i, j, k) - TR0*h_layer(i, j, k)) < 1.0e-12_wp, &
                          "PPM uniform tracer with uniform velocity must stay unchanged")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_ppm_uniform_tracer_no_change

   subroutine test_ppm_upwind_direction(error)
      !! A single tracer patch at (ic, jc) with uniform hu > 0 and
      !! hv = 0. After one PPM step the patch must move rightward:
      !! source cell loses mass, and the downstream cell gains more
      !! tracer than the upstream cell. (Unlike first-order upwind, PPM's
      !! wider stencil allows small amounts of tracer to appear upstream
      !! due to the 4th-order edge interpolation — this is expected.)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      real(wp) :: hTr_layer(NX, NY, NZ)
      real(wp), parameter :: dt = 0.05_wp
      real(wp), parameter :: U0 = 0.5_wp, TR_VALUE = 10.0_wp
      real(wp) :: hTr_src_before
      integer :: ic, jc, kc

      h_layer = 1.0_wp
      hu_layer = U0*h_layer
      hv_layer = 0.0_wp
      hTr_layer = 0.0_wp
      mass_flux_x = real(NZ, wp)*U0
      mass_flux_y = 0.0_wp

      ic = NX/2
      jc = NY/2
      kc = 1
      hTr_layer(ic, jc, kc) = TR_VALUE*h_layer(ic, jc, kc)
      hTr_src_before = hTr_layer(ic, jc, kc)

      call tracer_scratch_ensure(scratch, NX, NY, NZ, &
                                 want_ppm=.true., want_upwind=.false., want_hdiff=.false.)
      !$acc data copyin(h_layer, hu_layer, hv_layer, mass_flux_x, mass_flux_y) &
      !$acc&     copy(hTr_layer)
      call ensure_grid()
      call ml_advect_tracer_ppm(h_layer, hu_layer, hv_layer, &
                                mass_flux_x, mass_flux_y, &
                                hTr_layer, &
                                scratch%ppm_tr_snap, scratch%ppm_ql_x, scratch%ppm_qr_x, &
                                scratch%ppm_ql_y, scratch%ppm_qr_y, &
                                test_grid, NZ, dt)
      !$acc end data
      call tracer_scratch_destroy(scratch)

      call check(error, hTr_layer(ic, jc, kc) < hTr_src_before, &
                 "PPM advect: source cell must lose tracer with hu > 0")
      if (allocated(error)) return
      call check(error, hTr_layer(ic + 1, jc, kc) > 0.0_wp, &
                 "PPM advect: downstream cell must gain tracer with hu > 0")
      if (allocated(error)) return
      call check(error, &
                 hTr_layer(ic + 1, jc, kc) > abs(hTr_layer(ic - 1, jc, kc)), &
                 "PPM advect: downstream gain must exceed any upstream leakage")
   end subroutine test_ppm_upwind_direction

   subroutine test_ppm_less_diffusive_than_upwind(error)
      !! Advect a step function in the x-direction with both first-order
      !! upwind and PPM. After one step, PPM should produce a sharper
      !! front (less numerical diffusion), measured by the L2 error
      !! against the exact advected step.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      real(wp) :: hTr_upwind(NX, NY, NZ), hTr_ppm(NX, NY, NZ)
      real(wp), parameter :: dt = 0.1_wp
      real(wp), parameter :: U0 = 0.5_wp
      real(wp), parameter :: TR_LO = 0.0_wp, TR_HI = 10.0_wp
      real(wp) :: l2_upwind, l2_ppm, tr_exact, diff
      integer :: i, j, k

      h_layer = 1.0_wp
      hu_layer = U0*h_layer
      hv_layer = 0.0_wp
      mass_flux_x = real(NZ, wp)*U0
      mass_flux_y = 0.0_wp

      ! Step function: Tr = TR_HI for i <= NX/2, TR_LO otherwise
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               if (i <= NX/2) then
                  hTr_upwind(i, j, k) = TR_HI*h_layer(i, j, k)
               else
                  hTr_upwind(i, j, k) = TR_LO*h_layer(i, j, k)
               end if
            end do
         end do
      end do
      hTr_ppm = hTr_upwind

      call tracer_scratch_ensure(scratch, NX, NY, NZ, &
                                 want_ppm=.true., want_upwind=.true., want_hdiff=.false.)

      ! Advect with first-order upwind
      !$acc data copyin(h_layer, hu_layer, hv_layer, mass_flux_x, mass_flux_y) &
      !$acc&     copy(hTr_upwind)
      call ensure_grid()
      call ml_advect_tracer(h_layer, hu_layer, hv_layer, &
                            mass_flux_x, mass_flux_y, &
                            hTr_upwind, scratch%adv_hTr_snap, test_grid, NZ, dt)
      !$acc end data

      ! Advect with PPM
      !$acc data copyin(h_layer, hu_layer, hv_layer, mass_flux_x, mass_flux_y) &
      !$acc&     copy(hTr_ppm)
      call ml_advect_tracer_ppm(h_layer, hu_layer, hv_layer, &
                                mass_flux_x, mass_flux_y, &
                                hTr_ppm, &
                                scratch%ppm_tr_snap, scratch%ppm_ql_x, scratch%ppm_qr_x, &
                                scratch%ppm_ql_y, scratch%ppm_qr_y, &
                                test_grid, NZ, dt)
      !$acc end data
      call tracer_scratch_destroy(scratch)

      ! Compute L2 error against the exact shifted step.
      ! The step moves by U0*dt/DX cells = 0.05 cells rightward (sub-cell).
      ! The exact solution is the same step (shift < 1 cell).
      l2_upwind = 0.0_wp
      l2_ppm = 0.0_wp
      k = 1
      do j = NG + 1, NY - NG
         do i = NG + 1, NX - NG
            if (i <= NX/2) then
               tr_exact = TR_HI
            else
               tr_exact = TR_LO
            end if
            diff = hTr_upwind(i, j, k)/h_layer(i, j, k) - tr_exact
            l2_upwind = l2_upwind + diff*diff
            diff = hTr_ppm(i, j, k)/h_layer(i, j, k) - tr_exact
            l2_ppm = l2_ppm + diff*diff
         end do
      end do

      call check(error, l2_ppm < l2_upwind, &
                 "PPM must produce less numerical diffusion than first-order upwind")
   end subroutine test_ppm_less_diffusive_than_upwind

end module test_multilayer
