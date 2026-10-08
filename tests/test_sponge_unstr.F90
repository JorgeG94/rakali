!! Unit tests for the unstructured sponge relaxation kernel
module test_sponge_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t, &
                             mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, &
                             mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, &
                              state_unstr_init, &
                              state_unstr_destroy, &
                              state_unstr_enter_data, &
                              state_unstr_exit_data
   use rki_boundary_unstr, only: BC_WALL, BC_SPONGE
   use rki_kernel_sponge_unstr, only: apply_sponge_unstr
   implicit none
   private

   public :: collect_sponge_unstr_tests

   real(wp), parameter :: TOL = 1.0e-12_wp

contains

   subroutine collect_sponge_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("no_sponge_group_no_op", test_no_sponge_group), &
                  new_unittest("zero_strength_no_op", test_zero_strength), &
                  new_unittest("sponge_relaxes_toward_href", test_sponge_relaxes), &
                  new_unittest("non_sponge_cell_untouched", test_non_sponge_cell_untouched) &
                  ]
   end subroutine collect_sponge_unstr_tests

   ! ----------------------------------------------------------------------
   ! Helper: two-triangle unit square
   !   verts: (0,0)=1, (1,0)=2, (1,1)=3, (0,1)=4
   !   tri 1: [1,2,3]  (lower-right triangle)
   !   tri 2: [1,3,4]  (upper-left triangle)
   !   shared (interior) edge: 1--3
   !   boundary edges: cell-1 edge-1 (2--3), cell-1 edge-3 (1--2)
   !                   cell-2 edge-1 (3--4), cell-2 edge-2 (1--4)
   ! ----------------------------------------------------------------------

   subroutine build_two_triangle_square(mesh)
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp) :: vx(4), vy(4)
      integer :: tv(3, 2)

      vx = [0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp]
      vy = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp]
      tv(:, 1) = [1, 2, 3]
      tv(:, 2) = [1, 3, 4]
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_two_triangle_square

   subroutine run_sponge(state, mesh, dt)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: dt

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call apply_sponge_unstr(state, mesh, dt)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_sponge

   ! ----------------------------------------------------------------------
   ! Tests
   ! ----------------------------------------------------------------------

   subroutine test_no_sponge_group(error)
      !! All bc_type default to BC_WALL (state default) — the sponge
      !! kernel must be a no-op: stage/xmom/ymom unchanged.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = [2.0_wp, 3.0_wp]
      state%barotropic%xmom_cv = [1.0_wp, -1.5_wp]
      state%barotropic%ymom_cv = [0.5_wp, 0.75_wp]
      ! bc_type defaults to BC_WALL for all groups — no sponge active.

      call run_sponge(state, mesh, 0.5_wp)

      call check(error, abs(state%barotropic%stage_cv(1) - 2.0_wp) < TOL, "stage(1) untouched")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%stage_cv(2) - 3.0_wp) < TOL, "stage(2) untouched")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(1) - 1.0_wp) < TOL, "xmom(1) untouched")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(2) + 1.5_wp) < TOL, "xmom(2) untouched")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(1) - 0.5_wp) < TOL, "ymom(1) untouched")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(2) - 0.75_wp) < TOL, "ymom(2) untouched")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_no_sponge_group

   subroutine test_zero_strength(error)
      !! BC_SPONGE group with strength = 0 must be a no-op regardless of
      !! h_ref or mesh adjacency.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Tag boundary edges of cell 1 with group 1.
      mesh%bc_group(1, 1) = 1   ! cell 1, edge 1 (boundary)
      mesh%bc_group(3, 1) = 1   ! cell 1, edge 3 (boundary)
      state%bc_type(1) = BC_SPONGE
      state%bc_sponge_strength(1) = 0.0_wp   ! zero strength => no-op
      state%bc_sponge_h_ref(1) = 1.0_wp

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = [3.0_wp, 2.0_wp]
      state%barotropic%xmom_cv = [2.0_wp, 0.5_wp]
      state%barotropic%ymom_cv = [-1.0_wp, 0.3_wp]

      call run_sponge(state, mesh, 0.5_wp)

      call check(error, abs(state%barotropic%stage_cv(1) - 3.0_wp) < TOL, "zero str: stage(1) unchanged")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(1) - 2.0_wp) < TOL, "zero str: xmom(1) unchanged")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(1) + 1.0_wp) < TOL, "zero str: ymom(1) unchanged")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_zero_strength

   subroutine test_sponge_relaxes(error)
      !! Analytic check: cell 1 is adjacent to a BC_SPONGE edge.
      !!   alpha = 0.1, h_ref = 1.0, dt = 0.5
      !!   initial: bed = 0, stage = 2 -> h = 2, xmom = 3, ymom = 4
      !! Expected:
      !!   h_new = 2 - 0.5*0.1*(2 - 1) = 2 - 0.05 = 1.95 -> stage_new = 1.95
      !!   xmom_new = 3 * (1 - 0.5*0.1) = 3 * 0.95 = 2.85
      !!   ymom_new = 4 * 0.95 = 3.80
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp), parameter :: ALPHA = 0.1_wp, H_REF = 1.0_wp, DT_VAL = 0.5_wp
      real(wp), parameter :: H0 = 2.0_wp, XMOM0 = 3.0_wp, YMOM0 = 4.0_wp
      real(wp) :: h_new, xmom_new, ymom_new

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Tag one boundary edge of cell 1 with sponge group 1.
      mesh%bc_group(1, 1) = 1   ! cell 1, edge 1 (boundary: vertices 2--3)
      state%bc_type(1) = BC_SPONGE
      state%bc_sponge_strength(1) = ALPHA
      state%bc_sponge_h_ref(1) = H_REF

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = [H0, H0]
      state%barotropic%xmom_cv = [XMOM0, XMOM0]
      state%barotropic%ymom_cv = [YMOM0, YMOM0]

      call run_sponge(state, mesh, DT_VAL)

      ! Hand-computed expected values for cell 1.
      h_new = H0 - DT_VAL*ALPHA*(H0 - H_REF)     ! = 1.95
      xmom_new = XMOM0*(1.0_wp - DT_VAL*ALPHA)    ! = 2.85
      ymom_new = YMOM0*(1.0_wp - DT_VAL*ALPHA)    ! = 3.80

      call check(error, abs(state%barotropic%stage_cv(1) - h_new) < TOL, &
                 "cell1 stage relaxes toward h_ref")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(1) - xmom_new) < TOL, &
                 "cell1 xmom decays toward 0")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(1) - ymom_new) < TOL, &
                 "cell1 ymom decays toward 0")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_sponge_relaxes

   subroutine test_non_sponge_cell_untouched(error)
      !! Cell 1 is adjacent to BC_SPONGE boundary edges; cell 2 is NOT
      !! adjacent to any sponge edge (its boundary edges are assigned to
      !! group 2 which stays BC_WALL).  Cell 2 must be exactly unchanged.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp), parameter :: ALPHA = 0.2_wp, H_REF = 0.5_wp, DT_VAL = 0.1_wp
      real(wp), parameter :: STAGE2 = 5.0_wp, XMOM2 = 7.0_wp, YMOM2 = -3.0_wp

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Cell 1: boundary edges tagged with sponge group 1.
      mesh%bc_group(1, 1) = 1
      mesh%bc_group(3, 1) = 1
      state%bc_type(1) = BC_SPONGE
      state%bc_sponge_strength(1) = ALPHA
      state%bc_sponge_h_ref(1) = H_REF

      ! Cell 2: boundary edges tagged with group 2 which stays BC_WALL.
      mesh%bc_group(1, 2) = 2
      mesh%bc_group(2, 2) = 2
      ! state%bc_type(2) defaults to BC_WALL — no sponge for cell 2.

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = [2.0_wp, STAGE2]
      state%barotropic%xmom_cv = [1.0_wp, XMOM2]
      state%barotropic%ymom_cv = [0.5_wp, YMOM2]

      call run_sponge(state, mesh, DT_VAL)

      ! Cell 1 must have changed (sponge active).
      call check(error, abs(state%barotropic%stage_cv(1) - 2.0_wp) > 1.0e-15_wp, &
                 "cell1 stage must change (sponge active)")
      if (allocated(error)) return

      ! Cell 2 must be bit-identical to its initial values.
      call check(error, abs(state%barotropic%stage_cv(2) - STAGE2) < TOL, &
                 "cell2 stage untouched (no sponge edge)")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(2) - XMOM2) < TOL, &
                 "cell2 xmom untouched")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(2) - YMOM2) < TOL, &
                 "cell2 ymom untouched")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_non_sponge_cell_untouched

end module test_sponge_unstr
