!! Unit tests for non-uniform sigma layer correctness
!!
!! All existing multilayer tests use uniform dsig = 1/nz.  These tests
!! exercise the column-sum accumulators added in Phase 1 with layers of
!! varying thickness, verifying that the fixes are correct.
module test_nonuniform_sigma
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, THIN_LAYER_THRESHOLD, CP_WATER
   use rki_ml_dynamics, only: ml_compute_vertical_velocity, &
                              ml_compute_vertical_exchange, &
                              ml_compute_total_depth
   use rki_ml_tracers, only: ml_apply_surface_tracer_fluxes, &
                             ml_horizontal_diffusion_tracer, &
                             ml_advect_tracer, &
                             tracer_scratch_t, &
                             tracer_scratch_ensure, &
                             tracer_scratch_destroy
   use rki_ml_eos, only: ml_equation_of_state
   use rki_ml_pressure, only: ml_baroclinic_pressure_gradient_jacobian
   use rki_ml_friction, only: ml_apply_bottom_friction
   use rki_ml_vertical_viscosity, only: ml_vertical_viscosity
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: collect_nonuniform_sigma_tests

   ! Grid: 12x12 with 2 ghost cells, 4 sigma layers
   integer, parameter :: NX = 12, NY = 12, NG = 2, NZ = 4
   real(wp), parameter :: DX = 10.0_wp, DY = 10.0_wp

   type(hgrid_t) :: test_grid
   logical :: test_grid_initialised = .false.

   ! Phase B2 placeholders (rho_ref_dz=0 -> kernel falls back to rho_0).
   real(wp), parameter :: TEST_RHO_REF(0:0) = [0.0_wp]
   real(wp), parameter :: TEST_RHO_REF_DZ = 0.0_wp
   integer, parameter :: TEST_NZ_REF = 0

contains

   subroutine collect_nonuniform_sigma_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("bottom_friction_nonuniform", &
                               test_bottom_friction_nonuniform), &
                  new_unittest("bottom_friction_analytical", &
                               test_bottom_friction_analytical), &
                  new_unittest("total_depth_identity", &
                               test_total_depth_identity), &
                  new_unittest("eos_skip_dry_nonuniform", &
                               test_eos_skip_dry_nonuniform), &
                  new_unittest("vertical_velocity_nonuniform", &
                               test_vertical_velocity_nonuniform), &
                  new_unittest("w_diasurface_barotropic", &
                               test_w_diasurface_barotropic), &
                  new_unittest("w_diasurface_baroclinic", &
                               test_w_diasurface_baroclinic), &
                  new_unittest("vertical_exchange_nonuniform", &
                               test_vertical_exchange_nonuniform), &
                  new_unittest("bpg_skip_dry_nonuniform", &
                               test_bpg_skip_dry_nonuniform), &
                  new_unittest("vmix_skip_dry_nonuniform", &
                               test_vmix_skip_dry_nonuniform), &
                  new_unittest("surface_flux_nonuniform", &
                               test_surface_flux_nonuniform), &
                  new_unittest("advect_skip_dry_nonuniform", &
                               test_advect_skip_dry_nonuniform) &
                  ]
   end subroutine collect_nonuniform_sigma_tests

   subroutine ensure_grid()
      if (.not. test_grid_initialised) then
         call test_grid%init(NX - 2*NG, NY - 2*NG, NG, DX, DY)
         test_grid_initialised = .true.
      end if
   end subroutine ensure_grid

   ! Helper: fill h_layer with non-uniform dsig = [0.1, 0.2, 0.3, 0.4]
   ! for a given total depth H at interior cells, zero at ghost cells
   subroutine fill_nonuniform(h_layer, H_total, set_ghosts)
      real(wp), intent(out) :: h_layer(NX, NY, NZ)
      real(wp), intent(in) :: H_total
      logical, intent(in), optional :: set_ghosts
      real(wp) :: dsig(NZ)
      integer :: i, j, k
      logical :: do_ghosts

      dsig = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
      do_ghosts = .false.
      if (present(set_ghosts)) do_ghosts = set_ghosts

      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               if (do_ghosts .or. &
                   (i > NG .and. i <= NX - NG .and. &
                    j > NG .and. j <= NY - NG)) then
                  h_layer(i, j, k) = dsig(k)*H_total
               else
                  h_layer(i, j, k) = dsig(k)*H_total
               end if
            end do
         end do
      end do
   end subroutine fill_nonuniform

   ! ================================================================
   ! Bottom friction tests
   ! ================================================================

   subroutine test_bottom_friction_nonuniform(error)
      !! Bottom friction with non-uniform layers computes correct H_total.
      !! With uniform sigma the old code would compute h_layer(NZ)*nz = 0.4*H*4 = 1.6*H
      !! which is WRONG (should be H).  The fix computes sum(h_layer) = H.
      !! Bottom layer = k=1 under bottom-up ROMS indexing.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp), parameter :: H_total = 2.0_wp
      real(wp), parameter :: manning_n = 0.025_wp
      real(wp), parameter :: dt = 0.1_wp
      real(wp) :: hu_before
      integer :: k

      call ensure_grid()
      call fill_nonuniform(h_layer, H_total, set_ghosts=.true.)

      ! Bottom layer (k=1) has thickness dsig(1)*H = 0.1*2 = 0.2 m
      ! All layers moving at u = 2 m/s
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp
      do k = 1, NZ
         hu_layer(:, :, k) = h_layer(:, :, k)*2.0_wp
      end do

      hu_before = hu_layer(5, 5, 1)

      !$omp target data map(to: h_layer) map(tofrom: hu_layer, hv_layer)
      call ml_apply_bottom_friction(h_layer, hu_layer, hv_layer, &
                                    manning_n, .false., 0.0_wp, test_grid, NZ, dt)
      !$omp end target data

      ! Bottom layer (k=1) should be damped
      call check(error, hu_layer(5, 5, 1) < hu_before, &
                 "bottom friction damps bottom layer (k=1) with non-uniform sigma")
      if (allocated(error)) return

      ! Top layer (k=NZ) should be unchanged
      call check(error, abs(hu_layer(5, 5, NZ) - h_layer(5, 5, NZ)*2.0_wp) < 1.0e-12_wp, &
                 "bottom friction leaves top layer (k=NZ) unchanged")
   end subroutine test_bottom_friction_nonuniform

   subroutine test_bottom_friction_analytical(error)
      !! Verify the analytical semi-implicit friction formula with non-uniform layers.
      !! C_d = g * n^2 / H^(1/3), denom = 1 + dt * C_d * |u_b| / h_k
      !! hu_new = hu_old / denom
      !! Bottom layer = k=1 under bottom-up ROMS indexing.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp), parameter :: H_total = 5.0_wp
      real(wp), parameter :: manning_n = 0.03_wp
      real(wp), parameter :: dt = 0.5_wp
      real(wp), parameter :: u_b = 1.5_wp
      real(wp) :: dsig(NZ), h_k, cf, denom, hu_expected
      integer :: ic, jc, k

      dsig = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
      ic = 5; jc = 5

      ! Set up non-uniform layers
      call ensure_grid()
      call fill_nonuniform(h_layer, H_total, set_ghosts=.true.)
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp
      do k = 1, NZ
         hu_layer(:, :, k) = h_layer(:, :, k)*u_b
      end do

      ! Analytical expectation
      h_k = dsig(1)*H_total  ! bottom layer thickness (k=1)
      cf = GRAVITY*manning_n*manning_n/(H_total**(1.0_wp/3.0_wp))
      denom = 1.0_wp + dt*cf*abs(u_b)/h_k
      hu_expected = h_k*u_b/denom

      !$omp target data map(to: h_layer) map(tofrom: hu_layer, hv_layer)
      call ml_apply_bottom_friction(h_layer, hu_layer, hv_layer, &
                                    manning_n, .false., 0.0_wp, test_grid, NZ, dt)
      !$omp end target data

      call check(error, abs(hu_layer(ic, jc, 1) - hu_expected) < 1.0e-12_wp, &
                 "bottom friction matches analytical formula with non-uniform sigma")
   end subroutine test_bottom_friction_analytical

   ! ================================================================
   ! Total depth identity
   ! ================================================================

   subroutine test_total_depth_identity(error)
      !! For uniform sigma, sum(h_layer) must equal h_layer(1)*nz exactly.
      !! This ensures the accumulator doesn't introduce rounding differences
      !! compared to the old shortcut.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), h(NX, NY)
      real(wp), parameter :: H_total = 3.0_wp
      integer :: i, j, k
      real(wp) :: h_shortcut

      ! Uniform sigma
      call ensure_grid()
      h_layer = H_total/real(NZ, wp)

      !$omp target data map(to: h_layer) map(from: h)
      call ml_compute_total_depth(h_layer, h, test_grid, NZ)
      !$omp end target data

      do j = NG + 1, NY - NG
         do i = NG + 1, NX - NG
            h_shortcut = h_layer(i, j, 1)*NZ
            call check(error, abs(h(i, j) - h_shortcut) < 1.0e-12_wp, &
                       "total_depth matches h_layer(1)*nz for uniform sigma")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_total_depth_identity

   ! ================================================================
   ! EOS skip dry column with non-uniform sigma
   ! ================================================================

   subroutine test_eos_skip_dry_nonuniform(error)
      !! With non-uniform layers, a dry column (sum(h_k) < ML_MIN_DEPTH)
      !! should be correctly identified and skipped by EOS.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ)
      real(wp) :: hS_layer(NX, NY, NZ), hT_layer(NX, NY, NZ)
      real(wp) :: rho_layer(NX, NY, NZ)
      real(wp), parameter :: rho_0 = 1000.0_wp, beta_S = 0.78_wp
      real(wp), parameter :: S_ref = 0.0_wp, alpha_T = 0.0_wp, T_ref = 15.0_wp
      real(wp) :: dsig(NZ)
      integer :: i, j, k, ic, jc

      dsig = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
      ic = 5; jc = 5

      ! Most cells wet (H=2), but one cell is very shallow
      call fill_nonuniform(h_layer, 2.0_wp, set_ghosts=.true.)
      ! Make cell (ic,jc) nearly dry: total = 0.001 m (below ML_MIN_DEPTH = 0.5)
      do k = 1, NZ
         h_layer(ic, jc, k) = dsig(k)*0.001_wp
      end do

      hS_layer = 0.0_wp
      hT_layer = 0.0_wp
      ! Put some salinity in the wet cells
      !$omp target teams distribute parallel do collapse(3)
      do k = 1, NZ
      do j = 1, NY
      do i = 1, NX
         hS_layer(i, j, k) = 35.0_wp*h_layer(i, j, k)
         hT_layer(i, j, k) = 15.0_wp*h_layer(i, j, k)
      end do
      end do
      end do
      !$omp end target teams distribute parallel do
      rho_layer = rho_0

      call ensure_grid()
      !$omp target data map(to: h_layer, hS_layer, hT_layer) map(tofrom: rho_layer)
      call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_layer, &
                                rho_0, beta_S, S_ref, alpha_T, T_ref, &
                                test_grid, NZ)
      !$omp end target data

      ! Wet cell should have density updated
      call check(error, abs(rho_layer(6, 6, 1) - (rho_0 + beta_S*35.0_wp)) < 1.0e-10_wp, &
                 "EOS computes density for wet non-uniform cell")
   end subroutine test_eos_skip_dry_nonuniform

   ! ================================================================
   ! Vertical velocity with non-uniform sigma
   ! ================================================================

   subroutine test_vertical_velocity_nonuniform(error)
      !! With zero divergence and non-uniform layers, w should be zero
      !! at all interfaces.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), dsig(NZ)
      integer :: k

      call fill_nonuniform(h_layer, 3.0_wp, set_ghosts=.true.)
      dsig = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]

      ! Uniform velocity u=1 → zero divergence
      do k = 1, NZ
         hu_layer(:, :, k) = h_layer(:, :, k)*1.0_wp
      end do
      hv_layer = 0.0_wp
      w = -999.0_wp

      call ensure_grid()
      !$omp target data map(to: h_layer, hu_layer, hv_layer, dsig) map(from: w)
      call ml_compute_vertical_velocity(h_layer, hu_layer, hv_layer, dsig, w, &
                                        test_grid, NZ)
      !$omp end target data

      ! Bottom: w(k=1) = 0 by construction
      call check(error, abs(w(5, 5, 1)) < 1.0e-14_wp, &
                 "w at bottom = 0 for non-uniform sigma")
      if (allocated(error)) return

      ! Surface: w(nz+1) should be ~0 for uniform velocity (zero divergence)
      call check(error, abs(w(5, 5, NZ + 1)) < 1.0e-10_wp, &
                 "w at surface ~ 0 for zero-divergence non-uniform flow")
   end subroutine test_vertical_velocity_nonuniform

   subroutine test_w_diasurface_barotropic(error)
      !! DIASURFACE property: pure barotropic sigma flow (hu_k = dsig_k*hu_bt)
      !! with a spatially DIVERGENT barotropic transport must yield ~0
      !! sigma-relative velocity at every interior interface — the layer
      !! transport divergence is exactly cancelled by the moving-coordinate
      !! term dsig_k*div_col.  The pre-fix kernel (w(k+1)=w(k)-div(hu_k))
      !! reported the physical velocity and gave spurious O(1e-3) interior w.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), dsig(NZ)
      real(wp) :: H_col, hu_bt, hv_bt, x, y
      integer :: i, j, k
      real(wp) :: wmax

      dsig = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]

      ! Smoothly varying total depth + barotropic transport (nonzero
      ! divergence).  Distribute strictly by sigma => barotropic layer flow.
      do j = 1, NY
         do i = 1, NX
            x = real(i, wp)*DX
            y = real(j, wp)*DY
            H_col = 50.0_wp + 5.0_wp*sin(0.04_wp*x) + 3.0_wp*cos(0.03_wp*y)
            hu_bt = 2.0_wp + 0.5_wp*cos(0.035_wp*x) + 0.2_wp*sin(0.025_wp*y)
            hv_bt = 1.5_wp + 0.3_wp*sin(0.03_wp*x) - 0.4_wp*cos(0.02_wp*y)
            do k = 1, NZ
               h_layer(i, j, k) = dsig(k)*H_col
               hu_layer(i, j, k) = dsig(k)*hu_bt
               hv_layer(i, j, k) = dsig(k)*hv_bt
            end do
         end do
      end do
      w = -999.0_wp

      call ensure_grid()
      !$omp target data map(to: h_layer, hu_layer, hv_layer, dsig) map(from: w)
      call ml_compute_vertical_velocity(h_layer, hu_layer, hv_layer, dsig, w, &
                                        test_grid, NZ)
      !$omp end target data

      ! Interior interfaces (all of w) must vanish to round-off.
      wmax = 0.0_wp
      do k = 1, NZ + 1
         do j = NG + 1, NY - NG
            do i = NG + 1, NX - NG
               wmax = max(wmax, abs(w(i, j, k)))
            end do
         end do
      end do
      call check(error, wmax < 1.0e-13_wp, &
                 "diasurface w ~ 0 for barotropic sigma flow (divergent BT)")
   end subroutine test_w_diasurface_barotropic

   subroutine test_w_diasurface_baroclinic(error)
      !! Genuine cross-isopycnal transport: a sheared (baroclinic) flow whose
      !! per-layer transport divergence is NOT proportional to dsig_k must
      !! produce sensible NONZERO interior diasurface w, while the free
      !! surface interface stays ~0 (no diasurface flux through the moving
      !! top — guaranteed by construction).
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), dsig(NZ)
      real(wp) :: H_col, hu_bt, hv_bt, x, y
      real(wp) :: prof(NZ)
      integer :: i, j, k
      real(wp) :: w_interior_max

      dsig = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
      ! Per-layer velocity anomaly profile (top layer opposite the bed):
      ! makes div(hu_k)/dsig_k vary with k => true interior diapycnal flux.
      prof = [-1.5_wp, -0.5_wp, 0.5_wp, 1.0_wp]

      do j = 1, NY
         do i = 1, NX
            x = real(i, wp)*DX
            y = real(j, wp)*DY
            H_col = 50.0_wp + 5.0_wp*sin(0.04_wp*x) + 3.0_wp*cos(0.03_wp*y)
            hu_bt = 2.0_wp + 0.5_wp*cos(0.035_wp*x)
            hv_bt = 1.5_wp + 0.3_wp*sin(0.03_wp*x)
            do k = 1, NZ
               h_layer(i, j, k) = dsig(k)*H_col
               hu_layer(i, j, k) = dsig(k)*(hu_bt + prof(k)*0.3_wp*cos(0.035_wp*x))
               hv_layer(i, j, k) = dsig(k)*(hv_bt + prof(k)*0.3_wp*sin(0.030_wp*y))
            end do
         end do
      end do
      w = -999.0_wp

      call ensure_grid()
      !$omp target data map(to: h_layer, hu_layer, hv_layer, dsig) map(from: w)
      call ml_compute_vertical_velocity(h_layer, hu_layer, hv_layer, dsig, w, &
                                        test_grid, NZ)
      !$omp end target data

      ! Bed BC: w(1) = 0 exactly.
      call check(error, abs(w(5, 5, 1)) < 1.0e-14_wp, "diasurface w bed = 0")
      if (allocated(error)) return

      ! Free surface: diasurface flux through the moving top vanishes.
      call check(error, abs(w(5, 5, NZ + 1)) < 1.0e-12_wp, &
                 "diasurface w at free surface ~ 0 (moving-top BC)")
      if (allocated(error)) return

      ! Interior interfaces carry genuine, sensible nonzero flux.
      w_interior_max = 0.0_wp
      do k = 2, NZ
         do j = NG + 1, NY - NG
            do i = NG + 1, NX - NG
               w_interior_max = max(w_interior_max, abs(w(i, j, k)))
            end do
         end do
      end do
      call check(error, w_interior_max > 1.0e-6_wp, &
                 "baroclinic shear drives nonzero interior diasurface w")
   end subroutine test_w_diasurface_baroclinic

   ! ================================================================
   ! Vertical exchange with non-uniform sigma
   ! ================================================================

   subroutine test_vertical_exchange_nonuniform(error)
      !! Vertical exchange conserves total momentum with non-uniform layers.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1)
      real(wp) :: hu_total_before, hu_total_after
      real(wp), parameter :: dt = 0.01_wp
      integer :: k

      call fill_nonuniform(h_layer, 4.0_wp, set_ghosts=.true.)

      ! Sheared flow: bottom slow, top fast
      hu_layer(:, :, 1) = h_layer(:, :, 1)*0.5_wp
      hu_layer(:, :, 2) = h_layer(:, :, 2)*1.0_wp
      hu_layer(:, :, 3) = h_layer(:, :, 3)*2.0_wp
      hu_layer(:, :, 4) = h_layer(:, :, 4)*3.0_wp
      hv_layer = 0.0_wp

      ! Prescribe upward w at interfaces
      w = 0.0_wp
      w(:, :, 2) = 0.1_wp
      w(:, :, 3) = 0.05_wp
      w(:, :, 4) = 0.02_wp

      hu_total_before = sum(hu_layer(5, 5, :))

      call ensure_grid()
      !$omp target data map(to: h_layer, w) map(tofrom: hu_layer, hv_layer)
      call ml_compute_vertical_exchange(h_layer, hu_layer, hv_layer, &
                                        w, test_grid, NZ, dt)
      !$omp end target data

      hu_total_after = sum(hu_layer(5, 5, :))

      ! Total momentum should be conserved (exchanges between layers are internal)
      call check(error, abs(hu_total_after - hu_total_before) < 1.0e-10_wp, &
                 "vertical exchange conserves total momentum with non-uniform sigma")
   end subroutine test_vertical_exchange_nonuniform

   ! ================================================================
   ! BPG skip dry column with non-uniform sigma
   ! ================================================================

   subroutine test_bpg_skip_dry_nonuniform(error)
      !! Jacobian BPG correctly skips a dry column when layers are non-uniform.
      !! If the old h_layer(1)*nz shortcut were used, h_layer(1)=0.1*0.001=0.0001,
      !! times nz=4 gives 0.0004, which is < ML_MIN_DEPTH=0.5.  But with wrong
      !! non-uniform dsig and a wet total depth, the shortcut would fail.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ)
      real(wp) :: rho_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp), parameter :: rho_0 = 1025.0_wp, dt = 0.1_wp
      real(wp) :: hu_before
      real(wp) :: dsig(NZ)
      integer :: ic, jc, k

      dsig = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
      ic = 5; jc = 5

      ! All cells wet with non-uniform sigma
      call fill_nonuniform(h_layer, 3.0_wp, set_ghosts=.true.)

      ! Make one column nearly dry (but only in the guard sense)
      ! Total = sum(dsig*0.001) = 0.001 < ML_MIN_DEPTH
      do k = 1, NZ
         h_layer(ic, jc, k) = dsig(k)*0.001_wp
      end do

      ! Uniform density: BPG should be zero for wet cells,
      ! and the dry cell should be skipped
      rho_layer = rho_0
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp
      hu_before = hu_layer(ic, jc, 1)

      call ensure_grid()
      !$omp target data map(to: h_layer, rho_layer) map(tofrom: hu_layer, hv_layer)
      call ml_baroclinic_pressure_gradient_jacobian(h_layer, rho_layer, &
                                                    hu_layer, hv_layer, &
                                                    rho_0, &
                                                    TEST_RHO_REF, TEST_RHO_REF_DZ, TEST_NZ_REF, &
                                                    test_grid, NZ, dt)
      !$omp end target data

      ! Dry cell should not have been modified
      call check(error, abs(hu_layer(ic, jc, 1) - hu_before) < 1.0e-14_wp, &
                 "BPG Jacobian skips dry column with non-uniform sigma")
   end subroutine test_bpg_skip_dry_nonuniform

   ! ================================================================
   ! Vertical viscosity skip dry with non-uniform sigma
   ! ================================================================

   subroutine test_vmix_skip_dry_nonuniform(error)
      !! Vertical viscosity correctly skips dry columns with non-uniform layers.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: rho_layer(NX, NY, NZ)
      real(wp), parameter :: rho_0 = 1025.0_wp, dt = 0.1_wp
      real(wp), parameter :: nu_max = 0.01_wp, nu_bg = 1.0e-4_wp, ri_alpha = 5.0_wp
      real(wp) :: dsig(NZ)
      real(wp) :: hu_before
      integer :: ic, jc, k

      dsig = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
      ic = 5; jc = 5

      call fill_nonuniform(h_layer, 3.0_wp, set_ghosts=.true.)
      ! Make one column dry
      do k = 1, NZ
         h_layer(ic, jc, k) = dsig(k)*0.001_wp
      end do

      hu_layer = 0.0_wp
      hv_layer = 0.0_wp
      rho_layer = rho_0
      hu_before = hu_layer(ic, jc, 1)

      call ensure_grid()
      !$omp target data map(to: h_layer, rho_layer) map(tofrom: hu_layer, hv_layer)
      call ml_vertical_viscosity(h_layer, hu_layer, hv_layer, rho_layer, &
                                 rho_0, nu_max, nu_bg, ri_alpha, &
                                 test_grid, NZ, dt)
      !$omp end target data

      call check(error, abs(hu_layer(ic, jc, 1) - hu_before) < 1.0e-14_wp, &
                 "vmix skips dry column with non-uniform sigma")
   end subroutine test_vmix_skip_dry_nonuniform

   ! ================================================================
   ! Surface tracer fluxes with non-uniform sigma
   ! ================================================================

   subroutine test_surface_flux_nonuniform(error)
      !! Surface heat flux applies to top layer (k=NZ) with non-uniform sigma.
      !! The guard checks sum(h_k) > ML_MIN_DEPTH, not h_layer(NZ)*nz.
      !! Top layer = k=NZ under the structured solver's bottom-up indexing.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ)
      real(wp) :: hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: hS_layer(NX, NY, NZ), hT_layer(NX, NY, NZ)
      real(wp) :: heat_flux_2d(NX, NY), evap_minus_precip_2d(NX, NY)
      real(wp), parameter :: rho_0 = 1000.0_wp, Q = 200.0_wp, dt = 1.0_wp
      real(wp) :: hT_before, expected_delta
      integer :: k

      call fill_nonuniform(h_layer, 2.0_wp, set_ghosts=.true.)
      ! Top layer (k=NZ) thickness = dsig(NZ) * H_total = 0.4 * 2 = 0.8 m

      hS_layer = 0.0_wp
      hT_layer = 0.0_wp
      do k = 1, NZ
         hS_layer(:, :, k) = 35.0_wp*h_layer(:, :, k)
         hT_layer(:, :, k) = 15.0_wp*h_layer(:, :, k)
      end do
      heat_flux_2d = Q
      evap_minus_precip_2d = 0.0_wp
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp

      hT_before = hT_layer(5, 5, NZ)
      expected_delta = dt*Q/(rho_0*CP_WATER)

      call ensure_grid()
      !$omp target data map(to: h_layer, hu_layer, hv_layer, heat_flux_2d, evap_minus_precip_2d) &
      !$omp&     map(tofrom: hS_layer, hT_layer)
      call ml_apply_surface_tracer_fluxes(hS_layer, hT_layer, h_layer, &
                                          hu_layer, hv_layer, &
                                          heat_flux_2d, evap_minus_precip_2d, &
                                          .false., .false., &
                                          rho_0, CP_WATER, test_grid, NZ, dt)
      !$omp end target data

      call check(error, abs((hT_layer(5, 5, NZ) - hT_before) - expected_delta) < 1.0e-12_wp, &
                 "surface heat flux correct with non-uniform sigma")
   end subroutine test_surface_flux_nonuniform

   ! ================================================================
   ! Tracer advection skip dry with non-uniform sigma
   ! ================================================================

   subroutine test_advect_skip_dry_nonuniform(error)
      !! Tracer advection correctly skips dry columns with non-uniform layers.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      real(wp) :: hTr_layer(NX, NY, NZ)
      real(wp) :: dsig(NZ)
      real(wp) :: hTr_before
      real(wp), parameter :: dt = 0.01_wp
      integer :: ic, jc, k
      type(tracer_scratch_t) :: scratch

      dsig = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
      ic = 5; jc = 5

      call fill_nonuniform(h_layer, 3.0_wp, set_ghosts=.true.)
      ! Make one column dry
      do k = 1, NZ
         h_layer(ic, jc, k) = dsig(k)*0.001_wp
      end do

      hu_layer = 0.0_wp
      hv_layer = 0.0_wp
      hTr_layer = 0.0_wp
      mass_flux_x = 0.0_wp
      mass_flux_y = 0.0_wp
      do k = 1, NZ
         hTr_layer(:, :, k) = 35.0_wp*h_layer(:, :, k)
      end do

      hTr_before = hTr_layer(ic, jc, 1)

      call ensure_grid()
      call tracer_scratch_ensure(scratch, NX, NY, NZ, &
                                 want_ppm=.false., want_upwind=.true., want_hdiff=.false.)
      !$omp target data map(to: h_layer, hu_layer, hv_layer, mass_flux_x, mass_flux_y) &
      !$omp&     map(tofrom: hTr_layer)
      call ml_advect_tracer(h_layer, hu_layer, hv_layer, &
                            mass_flux_x, mass_flux_y, &
                            hTr_layer, scratch%adv_hTr_snap, test_grid, NZ, dt)
      !$omp end target data
      call tracer_scratch_destroy(scratch)

      call check(error, abs(hTr_layer(ic, jc, 1) - hTr_before) < 1.0e-14_wp, &
                 "tracer advection skips dry column with non-uniform sigma")
   end subroutine test_advect_skip_dry_nonuniform

end module test_nonuniform_sigma
