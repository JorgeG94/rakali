!! Unit tests for the biharmonic (del^4) lateral viscosity kernel.
module test_biharmonic
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_ml_horizontal_viscosity, only: ml_horizontal_viscosity_biharmonic, &
                                          ml_horizontal_viscosity_cleanup
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: collect_biharmonic_tests

   integer, parameter :: NX = 8, NY = 8, NG = 2, NZ = 2
   real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp

   type(hgrid_t), save :: grid
   logical, save :: grid_ready = .false.

contains

   subroutine ensure_grid()
      if (grid_ready) return
      call grid%init(NX - 2*NG, NY - 2*NG, NG, DX, DY)
      grid_ready = .true.
   end subroutine ensure_grid

   subroutine collect_biharmonic_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      call ensure_grid()
      testsuite = [ &
                  new_unittest("biharmonic_off_unchanged", test_off_unchanged), &
                  new_unittest("biharmonic_linear_field_untouched", test_linear_untouched), &
                  new_unittest("biharmonic_grid_mode_damped", test_grid_mode_damped) &
                  ]
   end subroutine collect_biharmonic_tests

   subroutine test_off_unchanged(error)
      !! nu_4 = 0 must be a no-op (default off = bit-identical).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NX, NY, NZ), hu(NX, NY, NZ), hv(NX, NY, NZ)

      h = 1.0_wp
      hu = 3.0_wp
      hv = -2.0_wp

      !$omp target data map(to:  h) map(tofrom:  hu, hv)
      call ml_horizontal_viscosity_biharmonic(h, hu, hv, 0.0_wp, grid, NZ, 0.1_wp)
      !$omp end target data

      call check(error, abs(hu(4, 4, 1) - 3.0_wp) < 1.0e-14_wp, &
                 "nu_4 = 0 leaves hu unchanged")
   end subroutine test_off_unchanged

   subroutine test_linear_untouched(error)
      !! A linear velocity field has del^4 u = 0, so the biharmonic operator
      !! must leave it untouched (exactly, in the interior).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NX, NY, NZ), hu(NX, NY, NZ), hv(NX, NY, NZ)
      real(wp), parameter :: SLOPE = 0.5_wp
      integer :: i, j

      h = 1.0_wp
      hv = 0.0_wp
      do j = 1, NY
         do i = 1, NX
            hu(i, j, :) = SLOPE*real(i, wp)   ! u = SLOPE*i (linear -> del^4 = 0)
         end do
      end do

      !$omp target data map(to:  h) map(tofrom:  hu, hv)
      call ml_horizontal_viscosity_biharmonic(h, hu, hv, 0.04_wp, grid, NZ, 0.1_wp)
      !$omp end target data

      call check(error, abs(hu(4, 4, 1) - SLOPE*4.0_wp) < 1.0e-12_wp, &
                 "biharmonic leaves a linear (del^4 = 0) field untouched")
   end subroutine test_linear_untouched

   subroutine test_grid_mode_damped(error)
      !! A grid-scale checkerboard is the biharmonic operator's stiffest mode.
      !! For (-1)^(i+j): del^2 u = -4(1/dx^2 + 1/dy^2) u, so
      !! del^4 u = 16(1/dx^2 + 1/dy^2)^2 u, and one forward step damps the
      !! amplitude by the exact factor 1 - dt*nu_4*16*(1/dx^2 + 1/dy^2)^2.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NX, NY, NZ), hu(NX, NY, NZ), hv(NX, NY, NZ)
      real(wp), parameter :: AMP = 1.0_wp, NU4 = 0.04_wp, DT = 0.1_wp
      real(wp) :: lambda, factor, hu0
      integer :: i, j

      h = 1.0_wp
      hv = 0.0_wp
      do j = 1, NY
         do i = 1, NX
            hu(i, j, :) = AMP*real((-1)**(i + j), wp)   ! checkerboard
         end do
      end do
      hu0 = hu(4, 4, 1)

      !$omp target data map(to:  h) map(tofrom:  hu, hv)
      call ml_horizontal_viscosity_biharmonic(h, hu, hv, NU4, grid, NZ, DT)
      !$omp end target data

      lambda = 16.0_wp*(1.0_wp/(DX*DX) + 1.0_wp/(DY*DY))**2
      factor = 1.0_wp - DT*NU4*lambda
      call check(error, abs(hu(4, 4, 1) - hu0*factor) < 1.0e-12_wp, &
                 "checkerboard damped by the exact biharmonic factor")
      if (allocated(error)) return
      call check(error, abs(hu(4, 4, 1)) < abs(hu0), &
                 "grid-scale mode amplitude decreases")

      call ml_horizontal_viscosity_cleanup()
   end subroutine test_grid_mode_damped

end module test_biharmonic
