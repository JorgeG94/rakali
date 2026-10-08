!! Unit tests for the state update and RK2 averaging kernels
module test_update
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, DRY_TOLERANCE
   use rki_kernel_update, only: compute_update, compute_rk2_average
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: collect_update_tests

contains

   subroutine collect_update_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("euler_step_with_constant_rhs", test_euler_constant_rhs), &
                  new_unittest("euler_step_zero_rhs", test_euler_zero_rhs), &
                  new_unittest("rk2_average_identical_states", test_rk2_identical), &
                  new_unittest("rk2_average_known_values", test_rk2_known_values), &
                  new_unittest("rk2_average_clamps_dry", test_rk2_clamps_dry), &
                  new_unittest("update_only_interior", test_update_only_interior) &
                  ]
   end subroutine collect_update_tests

   subroutine test_euler_constant_rhs(error)
      !! h = h + dt * flux_h with constant RHS
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: flux_h(NX, NY), flux_hu(NX, NY), flux_hv(NX, NY)
      real(wp), parameter :: DT = 0.1_wp
      type(hgrid_t) :: grid

      h = 1.0_wp
      hu = 0.5_wp
      hv = 0.3_wp
      flux_h = 2.0_wp
      flux_hu = -1.0_wp
      flux_hv = 0.5_wp

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$omp target data map(tofrom:  h, hu, hv) map(to:  flux_h, flux_hu, flux_hv)
      call compute_update(h, hu, hv, flux_h, flux_hu, flux_hv, grid, DT)
      !$omp end target data

      ! Interior cell (4,4): h = 1.0 + 0.1*2.0 = 1.2
      call check(error, abs(h(4, 4) - 1.2_wp) < 1.0e-14_wp, &
                 "h should be 1.0 + dt*flux_h = 1.2")
      if (allocated(error)) return
      ! hu = 0.5 + 0.1*(-1.0) = 0.4
      call check(error, abs(hu(4, 4) - 0.4_wp) < 1.0e-14_wp, &
                 "hu should be 0.5 + dt*flux_hu = 0.4")
      if (allocated(error)) return
      ! hv = 0.3 + 0.1*0.5 = 0.35
      call check(error, abs(hv(4, 4) - 0.35_wp) < 1.0e-14_wp, &
                 "hv should be 0.3 + dt*flux_hv = 0.35")
   end subroutine test_euler_constant_rhs

   subroutine test_euler_zero_rhs(error)
      !! Zero RHS should leave state unchanged
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: flux_h(NX, NY), flux_hu(NX, NY), flux_hv(NX, NY)
      type(hgrid_t) :: grid

      h = 2.5_wp
      hu = 1.0_wp
      hv = -0.5_wp
      flux_h = 0.0_wp
      flux_hu = 0.0_wp
      flux_hv = 0.0_wp

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$omp target data map(tofrom:  h, hu, hv) map(to:  flux_h, flux_hu, flux_hv)
      call compute_update(h, hu, hv, flux_h, flux_hu, flux_hv, grid, 0.1_wp)
      !$omp end target data

      call check(error, abs(h(4, 4) - 2.5_wp) < 1.0e-14_wp, &
                 "h should be unchanged with zero flux")
      if (allocated(error)) return
      call check(error, abs(hu(4, 4) - 1.0_wp) < 1.0e-14_wp, &
                 "hu should be unchanged with zero flux")
   end subroutine test_euler_zero_rhs

   subroutine test_rk2_identical(error)
      !! RK2 average of identical states should give the same state
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: h0(NX, NY), hu0(NX, NY), hv0(NX, NY)
      type(hgrid_t) :: grid

      h = 1.5_wp
      hu = 0.8_wp
      hv = -0.3_wp
      h0 = h
      hu0 = hu
      hv0 = hv

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$omp target data map(tofrom:  h, hu, hv) map(to:  h0, hu0, hv0)
      call compute_rk2_average(h, hu, hv, h0, hu0, hv0, grid)
      !$omp end target data

      call check(error, abs(h(4, 4) - 1.5_wp) < 1.0e-14_wp, &
                 "Average of identical states should be unchanged")
      if (allocated(error)) return
      call check(error, abs(hu(4, 4) - 0.8_wp) < 1.0e-14_wp, &
                 "hu average of identical states should be unchanged")
   end subroutine test_rk2_identical

   subroutine test_rk2_known_values(error)
      !! RK2 average: result = 0.5 * (h0 + h)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: h0(NX, NY), hu0(NX, NY), hv0(NX, NY)
      type(hgrid_t) :: grid

      h0 = 1.0_wp
      hu0 = 2.0_wp
      hv0 = 0.0_wp
      h = 3.0_wp
      hu = 4.0_wp
      hv = 2.0_wp

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$omp target data map(tofrom:  h, hu, hv) map(to:  h0, hu0, hv0)
      call compute_rk2_average(h, hu, hv, h0, hu0, hv0, grid)
      !$omp end target data

      ! h = 0.5*(1.0 + 3.0) = 2.0
      call check(error, abs(h(4, 4) - 2.0_wp) < 1.0e-14_wp, &
                 "h should be 0.5*(1.0 + 3.0) = 2.0")
      if (allocated(error)) return
      ! hu = 0.5*(2.0 + 4.0) = 3.0
      call check(error, abs(hu(4, 4) - 3.0_wp) < 1.0e-14_wp, &
                 "hu should be 0.5*(2.0 + 4.0) = 3.0")
      if (allocated(error)) return
      ! hv = 0.5*(0.0 + 2.0) = 1.0
      call check(error, abs(hv(4, 4) - 1.0_wp) < 1.0e-14_wp, &
                 "hv should be 0.5*(0.0 + 2.0) = 1.0")
   end subroutine test_rk2_known_values

   subroutine test_rk2_clamps_dry(error)
      !! RK2 average: small positive h is preserved for mass conservation;
      !! momentum is zeroed since h < DRY_TOLERANCE.
      !! Only negative undershoots get clamped to zero.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: h0(NX, NY), hu0(NX, NY), hv0(NX, NY)
      real(wp) :: h_expected
      type(hgrid_t) :: grid

      ! Average of these gives h = 0.5*(0 + DRY_TOL/2) = DRY_TOL/4
      h0 = 0.0_wp
      hu0 = 0.0_wp
      hv0 = 0.0_wp
      h = 0.5_wp*DRY_TOLERANCE
      hu = 10.0_wp
      hv = 5.0_wp
      h_expected = 0.25_wp*DRY_TOLERANCE

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$omp target data map(tofrom:  h, hu, hv) map(to:  h0, hu0, hv0)
      call compute_rk2_average(h, hu, hv, h0, hu0, hv0, grid)
      !$omp end target data

      call check(error, abs(h(4, 4) - h_expected) < 1.0e-15_wp, &
                 "Small positive h preserved for mass conservation")
      if (allocated(error)) return
      call check(error, abs(hu(4, 4)) < 1.0e-15_wp, &
                 "hu should be zeroed in thin-layer cell after RK2 average")
   end subroutine test_rk2_clamps_dry

   subroutine test_update_only_interior(error)
      !! Update should only modify interior cells, not ghost cells
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: flux_h(NX, NY), flux_hu(NX, NY), flux_hv(NX, NY)
      real(wp) :: h_ghost_before
      type(hgrid_t) :: grid

      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      flux_h = 1.0_wp
      flux_hu = 0.0_wp
      flux_hv = 0.0_wp

      h_ghost_before = h(1, 1)

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$omp target data map(tofrom:  h, hu, hv) map(to:  flux_h, flux_hu, flux_hv)
      call compute_update(h, hu, hv, flux_h, flux_hu, flux_hv, grid, 0.1_wp)
      !$omp end target data

      ! Ghost cell at (1,1) should be untouched
      call check(error, abs(h(1, 1) - h_ghost_before) < 1.0e-15_wp, &
                 "Ghost cell h should not be modified by update")
      if (allocated(error)) return
      ! But interior cell should be updated
      call check(error, abs(h(4, 4) - 1.1_wp) < 1.0e-14_wp, &
                 "Interior cell should be updated: 1.0 + 0.1*1.0 = 1.1")
   end subroutine test_update_only_interior

end module test_update
