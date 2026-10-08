!! Unit tests for sponge layer relaxation kernel
module test_sponge
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_kernel_sponge, only: apply_sponge
   use rki_boundary_types, only: bc_params_t, BC_SPONGE, BC_WALL
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: collect_sponge_tests

contains

   subroutine collect_sponge_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("sponge_no_op_when_width_zero", test_noop_zero_width), &
                  new_unittest("sponge_no_op_when_strength_zero", test_noop_zero_strength), &
                  new_unittest("sponge_relaxes_toward_ref", test_relaxation), &
                  new_unittest("sponge_interior_unaffected", test_interior_unaffected), &
                  new_unittest("sponge_only_on_sponge_edges", test_only_sponge_edges) &
                  ]
   end subroutine collect_sponge_tests

   subroutine call_sponge(h, hu, hv, nx, ny, ng, dt, bc)
      !! Wrapper for backend-specific data management
      integer, intent(in) :: nx, ny, ng
      real(wp), intent(inout) :: h(nx, ny), hu(nx, ny), hv(nx, ny)
      real(wp), intent(in) :: dt
      type(bc_params_t), intent(in) :: bc

      type(hgrid_t) :: grid

      call grid%init(nx - 2*ng, ny - 2*ng, ng, 1.0_wp, 1.0_wp)

      !$omp target data map(tofrom:  h, hu, hv)
      call apply_sponge(h, hu, hv, grid, dt, bc)
      !$omp end target data
   end subroutine call_sponge

   subroutine test_noop_zero_width(error)
      !! Sponge with width=0 should not modify any fields
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 10, NY = 10, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: h_orig(NX, NY)
      type(bc_params_t) :: bc

      h = 2.0_wp
      hu = 1.0_wp
      hv = 0.5_wp
      h_orig = h

      bc%west = BC_SPONGE
      bc%sponge_width = 0
      bc%sponge_strength = 10.0_wp
      bc%sponge_h_ref = 1.0_wp

      call call_sponge(h, hu, hv, NX, NY, NG, 1.0_wp, bc)

      call check(error, maxval(abs(h - h_orig)) < 1.0e-14_wp, &
                 "Sponge with width=0 should not modify h")
   end subroutine test_noop_zero_width

   subroutine test_noop_zero_strength(error)
      !! Sponge with strength=0 should not modify any fields
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 10, NY = 10, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: h_orig(NX, NY)
      type(bc_params_t) :: bc

      h = 2.0_wp
      hu = 1.0_wp
      hv = 0.5_wp
      h_orig = h

      bc%west = BC_SPONGE
      bc%sponge_width = 5
      bc%sponge_strength = 0.0_wp
      bc%sponge_h_ref = 1.0_wp

      call call_sponge(h, hu, hv, NX, NY, NG, 1.0_wp, bc)

      call check(error, maxval(abs(h - h_orig)) < 1.0e-14_wp, &
                 "Sponge with strength=0 should not modify h")
   end subroutine test_noop_zero_strength

   subroutine test_relaxation(error)
      !! Sponge should relax h toward h_ref near the boundary
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 12, NY = 12, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      type(bc_params_t) :: bc
      real(wp) :: h_before

      h = 3.0_wp
      hu = 2.0_wp
      hv = 1.0_wp

      bc%west = BC_SPONGE
      bc%east = BC_WALL
      bc%south = BC_WALL
      bc%north = BC_WALL
      bc%sponge_width = 4
      bc%sponge_strength = 1.0_wp
      bc%sponge_h_ref = 1.0_wp

      h_before = h(NG + 1, NG + 1)

      call call_sponge(h, hu, hv, NX, NY, NG, 0.1_wp, bc)

      ! h should have moved toward h_ref = 1.0 (was 3.0)
      call check(error, h(NG + 1, NG + 1) < h_before, &
                 "Sponge should relax h toward h_ref (decrease from 3 toward 1)")
      if (allocated(error)) return

      ! hu should have moved toward 0
      call check(error, abs(hu(NG + 1, NG + 1)) < 2.0_wp, &
                 "Sponge should relax hu toward zero")
   end subroutine test_relaxation

   subroutine test_interior_unaffected(error)
      !! Cells outside sponge width should not be modified
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 20, NY = 20, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      type(bc_params_t) :: bc
      integer :: i_mid, j_mid

      h = 3.0_wp
      hu = 2.0_wp
      hv = 1.0_wp

      bc%west = BC_SPONGE
      bc%east = BC_WALL
      bc%south = BC_WALL
      bc%north = BC_WALL
      bc%sponge_width = 3
      bc%sponge_strength = 10.0_wp
      bc%sponge_h_ref = 1.0_wp

      i_mid = NX/2
      j_mid = NY/2

      call call_sponge(h, hu, hv, NX, NY, NG, 0.1_wp, bc)

      call check(error, abs(h(i_mid, j_mid) - 3.0_wp) < 1.0e-14_wp, &
                 "Interior cells outside sponge should be unmodified")
   end subroutine test_interior_unaffected

   subroutine test_only_sponge_edges(error)
      !! Sponge should only apply near edges with BC_SPONGE, not BC_WALL edges
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 20, NY = 20, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      type(bc_params_t) :: bc

      h = 3.0_wp
      hu = 2.0_wp
      hv = 1.0_wp

      bc%west = BC_SPONGE
      bc%east = BC_WALL
      bc%south = BC_WALL
      bc%north = BC_WALL
      bc%sponge_width = 3
      bc%sponge_strength = 10.0_wp
      bc%sponge_h_ref = 1.0_wp

      call call_sponge(h, hu, hv, NX, NY, NG, 0.1_wp, bc)

      call check(error, abs(h(NX - NG, NY/2) - 3.0_wp) < 1.0e-14_wp, &
                 "Cell near east BC_WALL edge should not be sponge-damped")
   end subroutine test_only_sponge_edges

end module test_sponge
