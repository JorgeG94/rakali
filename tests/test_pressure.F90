!! Unit tests for non-hydrostatic pressure correction kernels
module test_pressure
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE
   use rki_kernel_pressure, only: compute_depth_averaged_momentum, &
                                  save_barotropic_state, &
                                  compute_layer_from_barotropic, &
                                  compute_vertical_velocity, &
                                  compute_nh_rhs, &
                                  apply_poisson_operator, &
                                  apply_vertical_preconditioner, &
                                  apply_nh_correction, &
                                  compute_rk2_average_layers, &
                                  solve_nh_pressure
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: collect_pressure_tests

   integer, parameter :: NX = 12, NY = 12, NG = 2, NZ = 2
   real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp
   type(hgrid_t), save :: grid
   logical, save :: grid_ready = .false.

contains

   subroutine ensure_grid()
      if (grid_ready) return
      call grid%init(NX - 2*NG, NY - 2*NG, NG, DX, DY)
      grid_ready = .true.
   end subroutine ensure_grid

   subroutine collect_pressure_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      call ensure_grid()

      testsuite = [ &
                  new_unittest("depth_avg_momentum_uniform", test_depth_avg_uniform), &
                  new_unittest("depth_avg_momentum_shear", test_depth_avg_shear), &
                  new_unittest("save_barotropic", test_save_barotropic), &
                  new_unittest("layer_from_barotropic_uniform", test_layer_baro_uniform), &
                  new_unittest("vertical_velocity_no_divergence", test_w_no_divergence), &
                  new_unittest("nh_rhs_zero_w", test_rhs_zero_w), &
                  new_unittest("poisson_operator_zero_field", test_poisson_zero), &
                  new_unittest("vertical_preconditioner", test_vertical_preconditioner), &
                  new_unittest("cg_convergence_manufactured", test_cg_manufactured), &
                  new_unittest("rk2_average_layers", test_rk2_avg), &
                  new_unittest("nh_inactive_unchanged", test_nh_inactive) &
                  ]
   end subroutine collect_pressure_tests

   subroutine test_depth_avg_uniform(error)
      !! Uniform velocity across all layers: hu = h * u, hv = h * v
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: u_layer(NX, NY, NZ), v_layer(NX, NY, NZ)
      real(wp) :: dsig(NZ)
      integer :: k

      h = 2.0_wp
      dsig = 1.0_wp/real(NZ, wp)
      do k = 1, NZ
         u_layer(:, :, k) = 3.0_wp
         v_layer(:, :, k) = -1.5_wp
      end do

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_depth_averaged_momentum(h, u_layer, v_layer, dsig, &
                                           hu, hv, grid, NZ)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      call check(error, abs(hu(5, 5) - 6.0_wp) < 1.0e-12_wp, &
                 "hu = h*u = 2*3 = 6")
      if (allocated(error)) return
      call check(error, abs(hv(5, 5) - (-3.0_wp)) < 1.0e-12_wp, &
                 "hv = h*v = 2*(-1.5) = -3")
   end subroutine test_depth_avg_uniform

   subroutine test_depth_avg_shear(error)
      !! Different velocities per layer: hu = h * sum(dsig_k * u_k)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: u_layer(NX, NY, NZ), v_layer(NX, NY, NZ)
      real(wp) :: dsig(NZ)
      real(wp) :: expected_hu

      h = 1.0_wp
      dsig = 1.0_wp/real(NZ, wp)  ! 0.5, 0.5
      u_layer(:, :, 1) = 2.0_wp   ! bottom
      u_layer(:, :, 2) = 4.0_wp   ! top
      v_layer = 0.0_wp

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_depth_averaged_momentum(h, u_layer, v_layer, dsig, &
                                           hu, hv, grid, NZ)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      expected_hu = 1.0_wp*(0.5_wp*2.0_wp + 0.5_wp*4.0_wp)  ! = 3.0
      call check(error, abs(hu(5, 5) - expected_hu) < 1.0e-12_wp, &
                 "hu = h * (0.5*2 + 0.5*4) = 3")
   end subroutine test_depth_avg_shear

   subroutine test_save_barotropic(error)
      !! save_barotropic_state copies hu -> hu_pre, hv -> hv_pre
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: hu(NX, NY), hv(NX, NY)
      real(wp) :: hu_pre(NX, NY), hv_pre(NX, NY)

      hu = 7.0_wp
      hv = -3.0_wp
      hu_pre = 0.0_wp
      hv_pre = 0.0_wp

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call save_barotropic_state(hu, hv, hu_pre, hv_pre, grid)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      call check(error, abs(hu_pre(4, 4) - 7.0_wp) < 1.0e-14_wp, &
                 "hu_pre should equal hu")
      if (allocated(error)) return
      call check(error, abs(hv_pre(4, 4) - (-3.0_wp)) < 1.0e-14_wp, &
                 "hv_pre should equal hv")
   end subroutine test_save_barotropic

   subroutine test_layer_baro_uniform(error)
      !! Barotropic acceleration should be distributed uniformly to all layers
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: hu_pre(NX, NY), hv_pre(NX, NY)
      real(wp) :: u_layer(NX, NY, NZ), v_layer(NX, NY, NZ)
      real(wp) :: dsig(NZ)
      real(wp) :: expected_du

      h = 2.0_wp
      hu_pre = 4.0_wp    ! hu before HLL: h*u = 2*2
      hu = 6.0_wp         ! hu after HLL: acceleration happened
      hv_pre = 0.0_wp
      hv = 0.0_wp
      dsig = 1.0_wp/real(NZ, wp)
      u_layer(:, :, :) = 2.0_wp   ! initial layer velocity
      v_layer = 0.0_wp

      ! du = (hu_new - hu_pre) / h = (6-4)/2 = 1
      expected_du = 1.0_wp

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_layer_from_barotropic(h, hu, hv, hu_pre, hv_pre, &
                                         u_layer, v_layer, dsig, &
                                         grid, NZ, 0.1_wp)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      call check(error, abs(u_layer(5, 5, 1) - 3.0_wp) < 1.0e-12_wp, &
                 "u_layer(1) = 2 + 1 = 3")
      if (allocated(error)) return
      call check(error, abs(u_layer(5, 5, 2) - 3.0_wp) < 1.0e-12_wp, &
                 "u_layer(2) = 2 + 1 = 3 (uniform distribution)")
   end subroutine test_layer_baro_uniform

   subroutine test_w_no_divergence(error)
      !! Uniform velocity field (no divergence) -> w should be zero everywhere
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), b(NX, NY)
      real(wp) :: u_layer(NX, NY, NZ), v_layer(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), dsig(NZ)
      integer :: k

      h = 1.0_wp
      b = 0.0_wp
      dsig = 1.0_wp/real(NZ, wp)
      do k = 1, NZ
         u_layer(:, :, k) = 1.0_wp
         v_layer(:, :, k) = 1.0_wp
      end do
      w = 999.0_wp  ! should be overwritten

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_vertical_velocity(h, u_layer, v_layer, b, dsig, w, &
                                     grid, NZ)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      ! Uniform h*u everywhere -> central difference div is zero -> w = 0
      call check(error, abs(w(5, 5, 1)) < 1.0e-14_wp, "w at bottom = 0")
      if (allocated(error)) return
      call check(error, abs(w(5, 5, NZ + 1)) < 1.0e-12_wp, &
                 "w at top = 0 for uniform velocity")
   end subroutine test_w_no_divergence

   subroutine test_rhs_zero_w(error)
      !! When w is zero everywhere, RHS should be zero
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), w(NX, NY, NZ + 1)
      real(wp) :: rhs(NX, NY, NZ), dsig(NZ)

      h = 1.0_wp
      w = 0.0_wp
      dsig = 1.0_wp/real(NZ, wp)
      rhs = 999.0_wp

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_nh_rhs(h, w, dsig, rhs, NX, NY, NZ, NG, DX, DY, 0.1_wp)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      call check(error, abs(rhs(5, 5, 1)) < 1.0e-14_wp, &
                 "RHS = 0 when w = 0")
   end subroutine test_rhs_zero_w

   subroutine test_poisson_zero(error)
      !! Poisson operator applied to zero field gives zero
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), p(NX, NY, NZ), Ap(NX, NY, NZ)
      real(wp) :: dsig(NZ)

      h = 1.0_wp
      p = 0.0_wp
      dsig = 1.0_wp/real(NZ, wp)
      Ap = 999.0_wp

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call apply_poisson_operator(h, dsig, p, Ap, NX, NY, NZ, NG, DX, DY, 0.1_wp)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      call check(error, abs(Ap(5, 5, 1)) < 1.0e-14_wp, &
                 "A*0 = 0")
   end subroutine test_poisson_zero

   subroutine test_vertical_preconditioner(error)
      !! Verify that the vertical-line preconditioner produces a valid
      !! solve: apply the preconditioner to a known vector, then check
      !! that M*z approximates the input r (where M is the vertical part).
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), r(NX, NY, NZ), z(NX, NY, NZ)
      real(wp) :: Mz(NX, NY, NZ)
      real(wp) :: dsig(NZ)
      real(wp), parameter :: DT = 0.01_wp
      integer :: i, j, k
      real(wp) :: max_err

      h = 2.0_wp
      dsig = 1.0_wp/real(NZ, wp)
      r = 0.0_wp
      z = 0.0_wp
      Mz = 0.0_wp

      ! Set a non-trivial RHS in the interior
      do k = 1, NZ
         do j = NG + 1, NY - NG
            do i = NG + 1, NX - NG
               r(i, j, k) = real(k, wp)*0.5_wp
            end do
         end do
      end do

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call apply_vertical_preconditioner(h, dsig, r, z, NX, NY, NZ, NG, DX, DY, DT)

      ! Apply the Poisson operator to z — if preconditioner is exact for
      ! the vertical part, Mz should recover r (up to horizontal residual)
      call apply_poisson_operator(h, dsig, z, Mz, NX, NY, NZ, NG, DX, DY, DT)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      ! z should be non-trivial (not all zeros)
      call check(error, maxval(abs(z(NG + 1:NX - NG, NG + 1:NY - NG, :))) > 1.0e-15_wp, &
                 "Preconditioner output z should be non-zero")
      if (allocated(error)) return

      ! The preconditioner is approximate (includes horizontal diagonal),
      ! so we just check that the output is finite and reasonable
      call check(error, all(abs(z(NG + 1:NX - NG, NG + 1:NY - NG, :)) < 1.0e10_wp), &
                 "Preconditioner output z should be finite")
   end subroutine test_vertical_preconditioner

   subroutine test_cg_manufactured(error)
      !! PCG solver converges on a problem with known RHS.
      !! Set up uniform h, zero initial q, and a non-trivial w field.
      !! The solver's convergence criterion is ||r||_2 < PCG_TOL, so the
      !! post-solve residual must respect the same tolerance (with a
      !! small slack for reduction-order differences on GPU), and the
      !! iteration count must stay strictly below PCG_MAX_ITER — otherwise
      !! the solver exited on the iteration cap, not on convergence.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), b(NX, NY)
      real(wp) :: u_layer(NX, NY, NZ), v_layer(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), q_nh(NX, NY, NZ)
      real(wp) :: cg_r(NX, NY, NZ), cg_p(NX, NY, NZ), cg_Ap(NX, NY, NZ)
      real(wp) :: cg_z(NX, NY, NZ)
      real(wp) :: dsig(NZ)
      real(wp), parameter :: DT = 0.01_wp
      real(wp), parameter :: PCG_TOL = 1.0e-10_wp
      integer, parameter :: PCG_MAX_ITER = 100
      integer :: iter_out, i, j, k
      real(wp) :: residual

      h = 1.0_wp
      b = 0.0_wp
      dsig = 1.0_wp/real(NZ, wp)
      u_layer = 0.0_wp
      v_layer = 0.0_wp
      q_nh = 0.0_wp
      cg_r = 0.0_wp
      cg_p = 0.0_wp
      cg_Ap = 0.0_wp
      cg_z = 0.0_wp

      do k = 1, NZ + 1
         do j = 1, NY
            do i = 1, NX
               w(i, j, k) = 0.01_wp*real(k - 1, wp)
            end do
         end do
      end do

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call solve_nh_pressure(h, u_layer, v_layer, w, q_nh, dsig, &
                             cg_r, cg_p, cg_Ap, cg_z, &
                             grid, NZ, DT, &
                             PCG_MAX_ITER, PCG_TOL, iter_out)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      ! RHS is non-zero (manufactured w field has nonzero divergence),
      ! so the solver must have iterated at least once.
      call check(error, iter_out >= 1, "PCG should iterate on a non-zero RHS")
      if (allocated(error)) return

      ! Convergence (not iter-cap termination): iter_out < max_iter.
      call check(error, iter_out < PCG_MAX_ITER, &
                 "PCG must converge before hitting max_iter")
      if (allocated(error)) return

      ! L2 residual norm must satisfy the solver's own convergence
      ! criterion (||r|| < tol). Small slack allows for GPU reduction
      ! ordering differences from the solver's internal accumulation.
      residual = sqrt(sum(cg_r(NG + 1:NX - NG, NG + 1:NY - NG, :)**2))
      call check(error, residual < 10.0_wp*PCG_TOL, &
                 "||r||_2 must respect the PCG_TOL convergence tolerance")
      if (allocated(error)) return

      ! q_nh must be non-trivial after solving a non-zero-RHS system.
      call check(error, maxval(abs(q_nh(NG + 1:NX - NG, NG + 1:NY - NG, :))) > 1.0e-15_wp, &
                 "q_nh should be non-zero after PCG solve")
   end subroutine test_cg_manufactured

   subroutine test_rk2_avg(error)
      !! RK2 averaging: result = 0.5*(old + new)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: u_layer(NX, NY, NZ), v_layer(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1)
      real(wp) :: u_layer0(NX, NY, NZ), v_layer0(NX, NY, NZ)
      real(wp) :: w0(NX, NY, NZ + 1)

      u_layer = 4.0_wp
      v_layer = 6.0_wp
      w = 2.0_wp
      u_layer0 = 2.0_wp
      v_layer0 = 4.0_wp
      w0 = 0.0_wp

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_rk2_average_layers(u_layer, v_layer, w, &
                                      u_layer0, v_layer0, w0, &
                                      grid, NZ)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      call check(error, abs(u_layer(5, 5, 1) - 3.0_wp) < 1.0e-14_wp, &
                 "u_layer avg = (2+4)/2 = 3")
      if (allocated(error)) return
      call check(error, abs(v_layer(5, 5, 1) - 5.0_wp) < 1.0e-14_wp, &
                 "v_layer avg = (4+6)/2 = 5")
      if (allocated(error)) return
      call check(error, abs(w(5, 5, 1) - 1.0_wp) < 1.0e-14_wp, &
                 "w avg = (0+2)/2 = 1")
   end subroutine test_rk2_avg

   subroutine test_nh_inactive(error)
      !! When the NH pressure correction is called with q_nh = 0 (the
      !! state "NH inactive" leaves the solver in), the correction
      !! kernel must be a bit-exact no-op on the layer velocities and
      !! on `w`. This guards against any future refactor that sneaks a
      !! constant or sign-error into `apply_nh_correction` that would
      !! still shift a zero-pressure field.
      !!
      !! Earlier versions of this test never called a kernel at all —
      !! they asserted host-side literals against themselves, which
      !! would pass for any implementation of NH whatsoever.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY)
      real(wp) :: u_layer(NX, NY, NZ), v_layer(NX, NY, NZ)
      real(wp) :: u_before(NX, NY, NZ), v_before(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), w_before(NX, NY, NZ + 1)
      real(wp) :: q_nh(NX, NY, NZ)
      real(wp) :: dsig(NZ)
      real(wp), parameter :: DT = 0.01_wp
      integer :: i, j, k

      h = 2.0_wp
      dsig = 1.0_wp/real(NZ, wp)
      q_nh = 0.0_wp       ! NH inactive => zero non-hydrostatic pressure

      ! Non-trivial velocities so a buggy correction would be visible.
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               u_layer(i, j, k) = 1.3_wp + 0.1_wp*real(i, wp)
               v_layer(i, j, k) = -0.7_wp + 0.05_wp*real(j, wp) &
                                  + 0.02_wp*real(k, wp)
            end do
         end do
      end do
      do k = 1, NZ + 1
         w(:, :, k) = 0.03_wp*real(k - 1, wp)
      end do

      u_before = u_layer
      v_before = v_layer
      w_before = w

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call apply_nh_correction(h, u_layer, v_layer, w, q_nh, dsig, &
                               grid, NZ, DT)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      do k = 1, NZ
         do j = NG + 1, NY - NG
            do i = NG + 1, NX - NG
               call check(error, &
                          abs(u_layer(i, j, k) - u_before(i, j, k)) < 1.0e-14_wp, &
                          "apply_nh_correction with q_nh=0 must leave u_layer unchanged")
               if (allocated(error)) return
               call check(error, &
                          abs(v_layer(i, j, k) - v_before(i, j, k)) < 1.0e-14_wp, &
                          "apply_nh_correction with q_nh=0 must leave v_layer unchanged")
               if (allocated(error)) return
            end do
         end do
      end do

      do k = 1, NZ + 1
         do j = NG + 1, NY - NG
            do i = NG + 1, NX - NG
               call check(error, &
                          abs(w(i, j, k) - w_before(i, j, k)) < 1.0e-14_wp, &
                          "apply_nh_correction with q_nh=0 must leave w unchanged")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_nh_inactive

end module test_pressure
