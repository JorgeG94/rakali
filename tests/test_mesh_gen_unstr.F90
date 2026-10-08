!! Unit tests for the native unstructured mesh generator.
module test_mesh_gen_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_mesh_gen_unstr, only: gen_regular_triangulation
   implicit none
   private

   public :: collect_mesh_gen_unstr_tests

contains

   subroutine collect_mesh_gen_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("mesh_gen_regular_counts_and_ids", test_counts), &
                  new_unittest("mesh_gen_regular_tiles_domain_ccw", test_tiles) &
                  ]
   end subroutine collect_mesh_gen_unstr_tests

   subroutine test_counts(error)
      !! Correct vertex/triangle counts, valid 1-based ids, corners on the box.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 4, NY = 3
      real(wp), parameter :: LX = 8.0_wp, LY = 6.0_wp, TOL = 1.0e-12_wp
      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :)
      integer :: nverts, ntri

      call gen_regular_triangulation(NX, NY, LX, LY, vx, vy, tv)
      nverts = (NX + 1)*(NY + 1)
      ntri = 2*NX*NY

      call check(error, size(vx) == nverts .and. size(vy) == nverts, &
                 "vertex arrays sized (nx+1)(ny+1)")
      if (allocated(error)) return
      call check(error, size(tv, 1) == 3 .and. size(tv, 2) == ntri, &
                 "tri_v sized (3, 2*nx*ny)")
      if (allocated(error)) return
      call check(error, minval(tv) >= 1 .and. maxval(tv) <= nverts, &
                 "all triangle vertex ids in [1, nverts]")
      if (allocated(error)) return
      call check(error, abs(vx(1)) < TOL .and. abs(vy(1)) < TOL .and. &
                 abs(vx(nverts) - LX) < TOL .and. abs(vy(nverts) - LY) < TOL, &
                 "corner vertices at (0,0) and (LX,LY)")
   end subroutine test_counts

   subroutine test_tiles(error)
      !! Every triangle positively oriented (CCW) and the triangulation tiles
      !! the rectangle exactly: sum of signed areas = LX*LY.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 5, NY = 4
      real(wp), parameter :: LX = 10.0_wp, LY = 7.0_wp
      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :)
      integer :: t
      real(wp) :: a, total, min_area
      real(wp) :: x1, y1, x2, y2, x3, y3

      call gen_regular_triangulation(NX, NY, LX, LY, vx, vy, tv)

      total = 0.0_wp
      min_area = huge(1.0_wp)
      do t = 1, size(tv, 2)
         x1 = vx(tv(1, t))
         y1 = vy(tv(1, t))
         x2 = vx(tv(2, t))
         y2 = vy(tv(2, t))
         x3 = vx(tv(3, t))
         y3 = vy(tv(3, t))
         a = 0.5_wp*((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))
         total = total + a
         min_area = min(min_area, a)
      end do

      call check(error, min_area > 0.0_wp, &
                 "every triangle is positively oriented (CCW)")
      if (allocated(error)) return
      call check(error, abs(total - LX*LY) < 1.0e-9_wp*LX*LY, &
                 "triangulation tiles the rectangle exactly (sum area = LX*LY)")
   end subroutine test_tiles

end module test_mesh_gen_unstr
