!! Unit tests for the unstructured mesh type
module test_mesh_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t, &
                             mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy
   implicit none
   private

   public :: collect_mesh_unstr_tests

   real(wp), parameter :: TOL = 1.0e-12_wp

contains

   subroutine collect_mesh_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("single_triangle_geometry", test_single_triangle_geometry), &
                  new_unittest("single_triangle_normals_outward", test_single_triangle_normals_outward), &
                  new_unittest("single_triangle_cfl_radius", test_single_triangle_cfl_radius), &
                  new_unittest("single_triangle_all_boundary", test_single_triangle_all_boundary), &
                  new_unittest("two_triangles_share_edge", test_two_triangles_share_edge), &
                  new_unittest("two_triangles_neighbour_symmetry", test_two_triangles_neighbour_symmetry), &
                  new_unittest("square_four_triangles_counts", test_square_four_triangles_counts), &
                  new_unittest("square_four_triangles_neighbours", test_square_four_triangles_neighbours), &
                  new_unittest("square_four_triangles_normals_outward", test_square_normals_outward), &
                  new_unittest("surrogate_self_on_boundary", test_surrogate_self_on_boundary), &
                  new_unittest("destroy_deallocates", test_destroy_deallocates) &
                  ]
   end subroutine collect_mesh_unstr_tests

   ! ----------------------------------------------------------------------
   ! Helper builders
   ! ----------------------------------------------------------------------

   subroutine build_unit_right_triangle(mesh)
      !! Right triangle with vertices (0,0), (1,0), (0,1) -- CCW.
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp) :: vx(3), vy(3)
      integer :: tv(3, 1)

      vx = [0.0_wp, 1.0_wp, 0.0_wp]
      vy = [0.0_wp, 0.0_wp, 1.0_wp]
      tv(:, 1) = [1, 2, 3]

      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_unit_right_triangle

   subroutine build_unit_square_two_triangles(mesh)
      !! Unit square split by the v1-v3 diagonal into two CCW triangles.
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp) :: vx(4), vy(4)
      integer :: tv(3, 2)

      vx = [0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp]
      vy = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp]
      tv(:, 1) = [1, 2, 3]
      tv(:, 2) = [1, 3, 4]

      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_unit_square_two_triangles

   subroutine build_unit_square_four_triangles(mesh)
      !! Unit square with a central vertex, split into 4 CCW triangles.
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp) :: vx(5), vy(5)
      integer :: tv(3, 4)

      vx = [0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.5_wp]
      vy = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.5_wp]
      ! Each triangle: outer edge + center vertex, walked CCW
      tv(:, 1) = [1, 2, 5]   ! bottom
      tv(:, 2) = [2, 3, 5]   ! right
      tv(:, 3) = [3, 4, 5]   ! top
      tv(:, 4) = [4, 1, 5]   ! left

      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_unit_square_four_triangles

   ! ----------------------------------------------------------------------
   ! Single-triangle tests
   ! ----------------------------------------------------------------------

   subroutine test_single_triangle_geometry(error)
      !! Area, centroid, and edge lengths of the unit right triangle.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      real(wp), parameter :: SQRT2 = sqrt(2.0_wp)

      call build_unit_right_triangle(mesh)

      call check(error, abs(mesh%area(1) - 0.5_wp) < TOL, "area should be 0.5")
      if (allocated(error)) return
      call check(error, abs(mesh%centroid_x(1) - 1.0_wp/3.0_wp) < TOL, &
                 "centroid x should be 1/3")
      if (allocated(error)) return
      call check(error, abs(mesh%centroid_y(1) - 1.0_wp/3.0_wp) < TOL, &
                 "centroid y should be 1/3")
      if (allocated(error)) return

      ! Edge 1 is opposite vertex 1 -- the hypotenuse, length sqrt(2)
      call check(error, abs(mesh%edge_len(1, 1) - SQRT2) < TOL, &
                 "edge 1 (hypotenuse) should have length sqrt(2)")
      if (allocated(error)) return
      ! Edge 2 is opposite vertex 2 -- the left side, length 1
      call check(error, abs(mesh%edge_len(2, 1) - 1.0_wp) < TOL, &
                 "edge 2 (left) should have length 1")
      if (allocated(error)) return
      ! Edge 3 is opposite vertex 3 -- the bottom, length 1
      call check(error, abs(mesh%edge_len(3, 1) - 1.0_wp) < TOL, &
                 "edge 3 (bottom) should have length 1")

      call mesh_unstr_destroy(mesh)
   end subroutine test_single_triangle_geometry

   subroutine test_single_triangle_normals_outward(error)
      !! For every edge, the vector from centroid to edge midpoint should
      !! have a positive dot product with the outward normal. This is a
      !! coordinate-free correctness invariant for any CCW triangle.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      integer :: i
      real(wp) :: dx, dy, dot, nlen

      call build_unit_right_triangle(mesh)

      do i = 1, 3
         dx = mesh%edge_mid_x(i, 1) - mesh%centroid_x(1)
         dy = mesh%edge_mid_y(i, 1) - mesh%centroid_y(1)
         dot = dx*mesh%edge_n1(i, 1) + dy*mesh%edge_n2(i, 1)
         call check(error, dot > 0.0_wp, &
                    "normal should point outward from centroid")
         if (allocated(error)) return

         nlen = sqrt(mesh%edge_n1(i, 1)**2 + mesh%edge_n2(i, 1)**2)
         call check(error, abs(nlen - 1.0_wp) < TOL, &
                    "normal should be unit length")
         if (allocated(error)) return
      end do

      call mesh_unstr_destroy(mesh)
   end subroutine test_single_triangle_normals_outward

   subroutine test_single_triangle_cfl_radius(error)
      !! CFL radius = min centroid-to-edge-midpoint distance (ANUGA `radii`).
      !! Unit right triangle (0,0),(1,0),(0,1): centroid = (1/3, 1/3). The
      !! nearest edge midpoint is the hypotenuse's at (1/2, 1/2), distance
      !! sqrt(2)/6 (= 1/3 of the hypotenuse median); the leg midpoints are
      !! farther at sqrt(5)/6. So r = sqrt(2)/6.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      real(wp) :: expected

      call build_unit_right_triangle(mesh)

      expected = sqrt(2.0_wp)/6.0_wp
      call check(error, abs(mesh%radius(1) - expected) < TOL, &
                 "CFL radius should be sqrt(2)/6")

      call mesh_unstr_destroy(mesh)
   end subroutine test_single_triangle_cfl_radius

   subroutine test_single_triangle_all_boundary(error)
      !! An isolated triangle has all three edges on the boundary.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      integer :: i

      call build_unit_right_triangle(mesh)

      call check(error, mesh%n_boundaries(1) == 3, &
                 "isolated triangle should have 3 boundary edges")
      if (allocated(error)) return
      call check(error, mesh%n_boundary_edges == 3, &
                 "mesh should report 3 boundary edges total")
      if (allocated(error)) return

      do i = 1, 3
         call check(error, mesh%neighbours(i, 1) == -1, &
                    "boundary neighbour should be -1")
         if (allocated(error)) return
         call check(error, mesh%neighbour_edges(i, 1) == -1, &
                    "boundary neighbour_edge should be -1")
         if (allocated(error)) return
         call check(error, mesh%surrogate_neighbours(i, 1) == 1, &
                    "surrogate of boundary edge should be self (1)")
         if (allocated(error)) return
      end do

      call mesh_unstr_destroy(mesh)
   end subroutine test_single_triangle_all_boundary

   ! ----------------------------------------------------------------------
   ! Two-triangle tests (square split by a diagonal)
   ! ----------------------------------------------------------------------

   subroutine test_two_triangles_share_edge(error)
      !! Triangle 1 = (1,2,3), triangle 2 = (1,3,4) share the diagonal v1-v3.
      !! For triangle 1 the diagonal is opposite local vertex 2 (the v2 corner)
      !! so it is local edge 2. For triangle 2 the diagonal is opposite local
      !! vertex 3 (the v4 corner) so it is local edge 3.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh

      call build_unit_square_two_triangles(mesh)

      call check(error, mesh%ncells_owned == 2, "should have 2 cells")
      if (allocated(error)) return

      call check(error, mesh%neighbours(2, 1) == 2, &
                 "triangle 1 edge 2 (diagonal) should neighbour triangle 2")
      if (allocated(error)) return
      call check(error, mesh%neighbour_edges(2, 1) == 3, &
                 "the matching edge of triangle 2 should be local edge 3")
      if (allocated(error)) return

      call check(error, mesh%neighbours(3, 2) == 1, &
                 "triangle 2 edge 3 (diagonal) should neighbour triangle 1")
      if (allocated(error)) return
      call check(error, mesh%neighbour_edges(3, 2) == 2, &
                 "the matching edge of triangle 1 should be local edge 2")
      if (allocated(error)) return

      call check(error, mesh%n_boundaries(1) == 2, &
                 "triangle 1 should have 2 boundary edges")
      if (allocated(error)) return
      call check(error, mesh%n_boundaries(2) == 2, &
                 "triangle 2 should have 2 boundary edges")
      if (allocated(error)) return
      call check(error, mesh%n_boundary_edges == 4, &
                 "shared-edge square should have 4 boundary edges total")

      call mesh_unstr_destroy(mesh)
   end subroutine test_two_triangles_share_edge

   subroutine test_two_triangles_neighbour_symmetry(error)
      !! For every linked edge, the reverse link must match.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      integer :: k, i, k2, e2

      call build_unit_square_two_triangles(mesh)

      do k = 1, mesh%ncells_owned
         do i = 1, 3
            k2 = mesh%neighbours(i, k)
            if (k2 < 0) cycle
            e2 = mesh%neighbour_edges(i, k)
            call check(error, mesh%neighbours(e2, k2) == k, &
                       "neighbour link should be symmetric (cell)")
            if (allocated(error)) return
            call check(error, mesh%neighbour_edges(e2, k2) == i, &
                       "neighbour link should be symmetric (edge index)")
            if (allocated(error)) return
         end do
      end do

      call mesh_unstr_destroy(mesh)
   end subroutine test_two_triangles_neighbour_symmetry

   ! ----------------------------------------------------------------------
   ! 4-triangle square tests
   ! ----------------------------------------------------------------------

   subroutine test_square_four_triangles_counts(error)
      !! 4 triangles meeting at a central vertex. Each triangle has one
      !! boundary edge (the outer side of the unit square). Total area
      !! must equal 1.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      integer :: k
      real(wp) :: total_area

      call build_unit_square_four_triangles(mesh)

      call check(error, mesh%ncells_owned == 4, "should have 4 cells")
      if (allocated(error)) return
      call check(error, mesh%nverts == 5, "should have 5 vertices")
      if (allocated(error)) return

      total_area = 0.0_wp
      do k = 1, 4
         total_area = total_area + mesh%area(k)
         call check(error, mesh%area(k) > 0.0_wp, &
                    "all cells should have positive (CCW) area")
         if (allocated(error)) return
         call check(error, mesh%n_boundaries(k) == 1, &
                    "each cell should have exactly 1 boundary edge")
         if (allocated(error)) return
      end do
      call check(error, abs(total_area - 1.0_wp) < TOL, &
                 "areas should sum to 1.0")
      if (allocated(error)) return
      call check(error, mesh%n_boundary_edges == 4, &
                 "should have 4 boundary edges total")

      call mesh_unstr_destroy(mesh)
   end subroutine test_square_four_triangles_counts

   subroutine test_square_four_triangles_neighbours(error)
      !! Every interior edge link must be symmetric.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      integer :: k, i, k2, e2, n_links

      call build_unit_square_four_triangles(mesh)

      n_links = 0
      do k = 1, mesh%ncells_owned
         do i = 1, 3
            k2 = mesh%neighbours(i, k)
            if (k2 < 0) cycle
            n_links = n_links + 1
            e2 = mesh%neighbour_edges(i, k)
            call check(error, k2 >= 1 .and. k2 <= mesh%ncells_owned, &
                       "neighbour cell index in range")
            if (allocated(error)) return
            call check(error, e2 >= 1 .and. e2 <= 3, &
                       "neighbour edge index in 1..3")
            if (allocated(error)) return
            call check(error, mesh%neighbours(e2, k2) == k, &
                       "reverse link cell should match")
            if (allocated(error)) return
            call check(error, mesh%neighbour_edges(e2, k2) == i, &
                       "reverse link edge should match")
            if (allocated(error)) return
         end do
      end do

      ! 4 interior edges, each visited twice (once per side) -> 8 links.
      call check(error, n_links == 8, &
                 "expected 8 directed neighbour links (4 interior edges * 2)")

      call mesh_unstr_destroy(mesh)
   end subroutine test_square_four_triangles_neighbours

   subroutine test_square_normals_outward(error)
      !! Coordinate-free outward-normal invariant on every cell+edge.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      integer :: k, i
      real(wp) :: dx, dy, dot

      call build_unit_square_four_triangles(mesh)

      do k = 1, mesh%ncells_owned
         do i = 1, 3
            dx = mesh%edge_mid_x(i, k) - mesh%centroid_x(k)
            dy = mesh%edge_mid_y(i, k) - mesh%centroid_y(k)
            dot = dx*mesh%edge_n1(i, k) + dy*mesh%edge_n2(i, k)
            call check(error, dot > 0.0_wp, &
                       "all normals should point outward from centroid")
            if (allocated(error)) return
         end do
      end do

      call mesh_unstr_destroy(mesh)
   end subroutine test_square_normals_outward

   ! ----------------------------------------------------------------------
   ! Surrogate neighbours and resource hygiene
   ! ----------------------------------------------------------------------

   subroutine test_surrogate_self_on_boundary(error)
      !! Wherever neighbours == -1, surrogate must equal the cell index.
      !! Wherever neighbours > 0, surrogate must equal neighbours.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      integer :: k, i

      call build_unit_square_four_triangles(mesh)

      do k = 1, mesh%ncells_owned
         do i = 1, 3
            if (mesh%neighbours(i, k) < 0) then
               call check(error, mesh%surrogate_neighbours(i, k) == k, &
                          "surrogate of boundary edge should be self")
            else
               call check(error, &
                          mesh%surrogate_neighbours(i, k) == mesh%neighbours(i, k), &
                          "surrogate of interior edge should equal neighbour")
            end if
            if (allocated(error)) return
         end do
      end do

      call mesh_unstr_destroy(mesh)
   end subroutine test_surrogate_self_on_boundary

   subroutine test_destroy_deallocates(error)
      !! All arrays should be released after destroy.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh

      call build_unit_square_four_triangles(mesh)
      call mesh_unstr_destroy(mesh)

      call check(error,.not. allocated(mesh%vert_x), "vert_x deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(mesh%tri_v), "tri_v deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(mesh%area), "area deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(mesh%edge_n1), "edge_n1 deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(mesh%neighbours), "neighbours deallocated")
      if (allocated(error)) return
      call check(error, mesh%ncells_owned == 0, "ncells_owned reset")
      if (allocated(error)) return
      call check(error, mesh%ncells_total == 0, "ncells_total reset")
   end subroutine test_destroy_deallocates

end module test_mesh_unstr
