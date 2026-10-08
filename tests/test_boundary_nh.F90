!! Unit tests for the non-hydrostatic boundary-condition kernel
!!
!! `test_boundary` covers `apply_boundary_conditions` (the depth-averaged
!! routine).  This module hits `apply_nh_boundary_conditions` directly and
!! also covers `bc_type_from_string` for every supported keyword.
module test_boundary_nh
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_grid, only: hgrid_t
   use rki_boundary, only: apply_nh_boundary_conditions, bc_params_t, &
                           BC_WALL, BC_OPEN, BC_TIDAL, BC_NESTED, &
                           BC_INFLOW, BC_DISCHARGE, BC_CLAMPED, &
                           BC_SPONGE, BC_CHAPMAN
   use rki_boundary_types, only: bc_type_from_string
   implicit none
   private

   public :: collect_boundary_nh_tests

   integer, parameter :: NX_PHYS = 4, NY_PHYS = 3, NG = 2, NZ = 2
   integer, parameter :: NX = NX_PHYS + 2*NG, NY = NY_PHYS + 2*NG

contains

   subroutine collect_boundary_nh_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("nh_wall_west_negates_u", test_nh_wall_west), &
                  new_unittest("nh_wall_east_negates_u", test_nh_wall_east), &
                  new_unittest("nh_wall_south_negates_v", test_nh_wall_south), &
                  new_unittest("nh_wall_north_negates_v", test_nh_wall_north), &
                  new_unittest("nh_open_extrapolates_velocity", test_nh_open_west), &
                  new_unittest("nh_q_zero_gradient", test_nh_q_neumann), &
                  new_unittest("nh_w_zero_gradient", test_nh_w_neumann), &
                  new_unittest("bc_type_from_string_all_keywords", test_bc_string_dispatch), &
                  new_unittest("bc_type_from_string_unknown_defaults_wall", test_bc_string_default) &
                  ]
   end subroutine collect_boundary_nh_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(NX_PHYS, NY_PHYS, NG, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine init_arrays(q_nh, u_layer, v_layer, w, h, b)
      real(wp), intent(out) :: q_nh(NX, NY, NZ)
      real(wp), intent(out) :: u_layer(NX, NY, NZ), v_layer(NX, NY, NZ)
      real(wp), intent(out) :: w(NX, NY, NZ + 1)
      real(wp), intent(out) :: h(NX, NY), b(NX, NY)
      integer :: i, j, k

      h = 1.0_wp
      b = 0.0_wp
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               q_nh(i, j, k) = real(i*100 + j, wp)
               u_layer(i, j, k) = 0.5_wp + 0.1_wp*real(i, wp)
               v_layer(i, j, k) = -0.3_wp - 0.05_wp*real(j, wp)
            end do
         end do
      end do
      do k = 1, NZ + 1
         do j = 1, NY
            do i = 1, NX
               w(i, j, k) = 0.01_wp*real(i + j + k, wp)
            end do
         end do
      end do
   end subroutine init_arrays

   ! ----------------------------------------------------------------------
   ! Wall BC: normal velocity is negated, tangential mirrored, q neumann.
   ! ----------------------------------------------------------------------

   subroutine test_nh_wall_west(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(bc_params_t) :: bc
      real(wp) :: q(NX, NY, NZ), u(NX, NY, NZ), v(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), h(NX, NY), b(NX, NY)
      integer :: g, j, k, i_first

      call make_grid(grid)
      call init_arrays(q, u, v, w, h, b)
      bc%west = BC_WALL; bc%east = BC_WALL
      bc%south = BC_WALL; bc%north = BC_WALL

      call apply_nh_boundary_conditions(q, u, v, w, h, b, grid, NZ, bc)

      i_first = NG + 1
      do k = 1, NZ
         do j = 1, NY
            do g = 1, NG
               call check(error, abs(u(i_first - g, j, k) + u(i_first + g - 1, j, k)) < 1.0e-12_wp, &
                          "wall: west u_layer mirrored-and-negated")
               if (allocated(error)) return
               call check(error, abs(v(i_first - g, j, k) - v(i_first + g - 1, j, k)) < 1.0e-12_wp, &
                          "wall: west v_layer mirrored")
               if (allocated(error)) return
               call check(error, abs(q(i_first - g, j, k) - q(i_first + g - 1, j, k)) < 1.0e-12_wp, &
                          "wall: q neumann at west")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_nh_wall_west

   subroutine test_nh_wall_east(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(bc_params_t) :: bc
      real(wp) :: q(NX, NY, NZ), u(NX, NY, NZ), v(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), h(NX, NY), b(NX, NY)
      integer :: g, j, k, i_last

      call make_grid(grid)
      call init_arrays(q, u, v, w, h, b)
      bc%west = BC_WALL; bc%east = BC_WALL
      bc%south = BC_WALL; bc%north = BC_WALL

      call apply_nh_boundary_conditions(q, u, v, w, h, b, grid, NZ, bc)

      i_last = NX - NG
      do k = 1, NZ
         do j = 1, NY
            do g = 1, NG
               call check(error, abs(u(i_last + g, j, k) + u(i_last - g + 1, j, k)) < 1.0e-12_wp, &
                          "wall: east u_layer negated")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_nh_wall_east

   subroutine test_nh_wall_south(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(bc_params_t) :: bc
      real(wp) :: q(NX, NY, NZ), u(NX, NY, NZ), v(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), h(NX, NY), b(NX, NY)
      integer :: g, i, k, j_first

      call make_grid(grid)
      call init_arrays(q, u, v, w, h, b)
      bc%west = BC_WALL; bc%east = BC_WALL
      bc%south = BC_WALL; bc%north = BC_WALL

      call apply_nh_boundary_conditions(q, u, v, w, h, b, grid, NZ, bc)

      j_first = NG + 1
      do k = 1, NZ
         do i = 1, NX
            do g = 1, NG
               call check(error, abs(v(i, j_first - g, k) + v(i, j_first + g - 1, k)) < 1.0e-12_wp, &
                          "wall: south v_layer negated")
               if (allocated(error)) return
               call check(error, abs(u(i, j_first - g, k) - u(i, j_first + g - 1, k)) < 1.0e-12_wp, &
                          "wall: south u_layer mirrored")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_nh_wall_south

   subroutine test_nh_wall_north(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(bc_params_t) :: bc
      real(wp) :: q(NX, NY, NZ), u(NX, NY, NZ), v(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), h(NX, NY), b(NX, NY)
      integer :: g, i, k, j_last

      call make_grid(grid)
      call init_arrays(q, u, v, w, h, b)
      bc%west = BC_WALL; bc%east = BC_WALL
      bc%south = BC_WALL; bc%north = BC_WALL

      call apply_nh_boundary_conditions(q, u, v, w, h, b, grid, NZ, bc)

      j_last = NY - NG
      do k = 1, NZ
         do i = 1, NX
            do g = 1, NG
               call check(error, abs(v(i, j_last + g, k) + v(i, j_last - g + 1, k)) < 1.0e-12_wp, &
                          "wall: north v_layer negated")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_nh_wall_north

   ! ----------------------------------------------------------------------
   ! Open / non-wall: velocities are mirrored without the sign flip.
   ! ----------------------------------------------------------------------

   subroutine test_nh_open_west(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(bc_params_t) :: bc
      real(wp) :: q(NX, NY, NZ), u(NX, NY, NZ), v(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), h(NX, NY), b(NX, NY)
      integer :: g, j, k, i_first

      call make_grid(grid)
      call init_arrays(q, u, v, w, h, b)
      bc%west = BC_OPEN; bc%east = BC_WALL
      bc%south = BC_WALL; bc%north = BC_WALL

      call apply_nh_boundary_conditions(q, u, v, w, h, b, grid, NZ, bc)

      i_first = NG + 1
      do k = 1, NZ
         do j = 1, NY
            do g = 1, NG
               call check(error, abs(u(i_first - g, j, k) - u(i_first + g - 1, j, k)) < 1.0e-12_wp, &
                          "open: u mirrored without negation")
               if (allocated(error)) return
               call check(error, abs(v(i_first - g, j, k) - v(i_first + g - 1, j, k)) < 1.0e-12_wp, &
                          "open: v mirrored")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_nh_open_west

   subroutine test_nh_q_neumann(error)
      !! q_nh has Neumann (zero-gradient) at all boundaries, regardless of bc type
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(bc_params_t) :: bc
      real(wp) :: q(NX, NY, NZ), u(NX, NY, NZ), v(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), h(NX, NY), b(NX, NY)
      integer :: g, j, k, i_first

      call make_grid(grid)
      call init_arrays(q, u, v, w, h, b)
      bc%west = BC_OPEN; bc%east = BC_TIDAL
      bc%south = BC_INFLOW; bc%north = BC_CLAMPED

      call apply_nh_boundary_conditions(q, u, v, w, h, b, grid, NZ, bc)

      i_first = NG + 1
      do k = 1, NZ
         do j = 1, NY
            do g = 1, NG
               call check(error, abs(q(i_first - g, j, k) - q(i_first + g - 1, j, k)) < 1.0e-12_wp, &
                          "q west neumann under non-wall BC")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_nh_q_neumann

   subroutine test_nh_w_neumann(error)
      !! w array has Neumann at all boundaries
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(bc_params_t) :: bc
      real(wp) :: q(NX, NY, NZ), u(NX, NY, NZ), v(NX, NY, NZ)
      real(wp) :: w(NX, NY, NZ + 1), h(NX, NY), b(NX, NY)
      integer :: g, j, k, i_last

      call make_grid(grid)
      call init_arrays(q, u, v, w, h, b)
      bc%west = BC_WALL; bc%east = BC_OPEN
      bc%south = BC_WALL; bc%north = BC_WALL

      call apply_nh_boundary_conditions(q, u, v, w, h, b, grid, NZ, bc)

      i_last = NX - NG
      do k = 1, NZ + 1
         do j = 1, NY
            do g = 1, NG
               call check(error, abs(w(i_last + g, j, k) - w(i_last - g + 1, j, k)) < 1.0e-12_wp, &
                          "w east neumann")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_nh_w_neumann

   ! ----------------------------------------------------------------------
   ! bc_type_from_string
   ! ----------------------------------------------------------------------

   subroutine test_bc_string_dispatch(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, bc_type_from_string("wall") == BC_WALL, "wall")
      if (allocated(error)) return
      call check(error, bc_type_from_string("open") == BC_OPEN, "open")
      if (allocated(error)) return
      call check(error, bc_type_from_string("tidal") == BC_TIDAL, "tidal")
      if (allocated(error)) return
      call check(error, bc_type_from_string("nested") == BC_NESTED, "nested")
      if (allocated(error)) return
      call check(error, bc_type_from_string("inflow") == BC_INFLOW, "inflow")
      if (allocated(error)) return
      call check(error, bc_type_from_string("discharge") == BC_DISCHARGE, "discharge")
      if (allocated(error)) return
      call check(error, bc_type_from_string("clamped") == BC_CLAMPED, "clamped")
      if (allocated(error)) return
      call check(error, bc_type_from_string("sponge") == BC_SPONGE, "sponge")
      if (allocated(error)) return
      call check(error, bc_type_from_string("chapman") == BC_CHAPMAN, "chapman")
   end subroutine test_bc_string_dispatch

   subroutine test_bc_string_default(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, bc_type_from_string("not_a_bc") == BC_WALL, &
                 "unknown keyword should default to wall")
      if (allocated(error)) return
      call check(error, bc_type_from_string("") == BC_WALL, &
                 "empty string should default to wall")
   end subroutine test_bc_string_default

end module test_boundary_nh
