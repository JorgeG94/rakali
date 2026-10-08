!! Smoke tests for the BC_NESTED unstructured boundary condition
module test_nested_unstr
   !! BC_NESTED is a smoke-level port: at v1 it is identical to BC_TIDAL
   !! (Flather radiation driven by per-group prescribed elevation
   !! `bc_eta_ext`). These tests verify:
   !!   * The solver step completes without NaN when a BC_NESTED edge is
   !!     active with a nonzero external stage.
   !!   * The BC_NESTED cell produces a different response from the
   !!     all-walls baseline, confirming that the case arm is reachable.
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
   use rki_boundary_unstr, only: BC_WALL, BC_NESTED
   use rki_solver_unstr, only: solver_unstr_step
   implicit none
   private

   public :: collect_nested_unstr_tests

contains

   subroutine collect_nested_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("nested_step_no_nan", test_nested_step_no_nan), &
                  new_unittest("nested_differs_from_wall", test_nested_differs_from_wall) &
                  ]
   end subroutine collect_nested_unstr_tests

   ! ----------------------------------------------------------------------
   ! Helpers
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

   ! ----------------------------------------------------------------------
   ! Tests
   ! ----------------------------------------------------------------------

   subroutine test_nested_step_no_nan(error)
      !! One solver step with a BC_NESTED boundary must produce finite
      !! stage, xmom, ymom at every cell.  Regression guard against the
      !! unimplemented fallthrough that would silently wall-bounce or
      !! NaN out.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: ETA_EXT = 1.2_wp
      integer :: c

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Tag one boundary edge of cell 1 as BC_NESTED group 1.
      mesh%bc_group(1, 1) = 1
      state%bc_type(1) = BC_NESTED
      state%bc_eta_ext(1) = ETA_EXT

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 1.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      state%dt = DT

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call solver_unstr_step(state, mesh, DT)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      do c = 1, state%ncells_owned
         call check(error, state%barotropic%stage_cv(c) == state%barotropic%stage_cv(c), &
                    "stage finite (no NaN) at cell")
         if (allocated(error)) return
         call check(error, state%barotropic%xmom_cv(c) == state%barotropic%xmom_cv(c), &
                    "xmom finite (no NaN) at cell")
         if (allocated(error)) return
         call check(error, state%barotropic%ymom_cv(c) == state%barotropic%ymom_cv(c), &
                    "ymom finite (no NaN) at cell")
         if (allocated(error)) return
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_nested_step_no_nan

   subroutine test_nested_differs_from_wall(error)
      !! BC_NESTED with eta_ext > initial stage must drive the adjacent
      !! cell to a different stage than the all-walls baseline after one
      !! step. Confirms the case arm is live (not silently falling through
      !! to BC_WALL).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh_wall, mesh_nested
      type(unstr_state_t) :: state_wall, state_nested
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: ETA_EXT = 2.0_wp  ! higher than initial 1.0

      ! ----- wall baseline -----
      call build_two_triangle_square(mesh_wall)
      call state_unstr_init(state_wall, mesh_wall)
      state_wall%barotropic%bed_cv = 0.0_wp
      state_wall%barotropic%stage_cv = 1.0_wp
      state_wall%barotropic%xmom_cv = 0.0_wp
      state_wall%barotropic%ymom_cv = 0.0_wp
      state_wall%dt = DT
      ! All bc_type default BC_WALL (no group tags on any edge).

      call mesh_unstr_enter_data(mesh_wall)
      call state_unstr_enter_data(state_wall)
      call solver_unstr_step(state_wall, mesh_wall, DT)
      call state_unstr_exit_data(state_wall)
      call mesh_unstr_exit_data(mesh_wall)

      ! ----- nested boundary -----
      call build_two_triangle_square(mesh_nested)
      call state_unstr_init(state_nested, mesh_nested)
      ! Tag cell-1 edge-1 as BC_NESTED group 1 with a raised external stage.
      mesh_nested%bc_group(1, 1) = 1
      state_nested%bc_type(1) = BC_NESTED
      state_nested%bc_eta_ext(1) = ETA_EXT

      state_nested%barotropic%bed_cv = 0.0_wp
      state_nested%barotropic%stage_cv = 1.0_wp
      state_nested%barotropic%xmom_cv = 0.0_wp
      state_nested%barotropic%ymom_cv = 0.0_wp
      state_nested%dt = DT

      call mesh_unstr_enter_data(mesh_nested)
      call state_unstr_enter_data(state_nested)
      call solver_unstr_step(state_nested, mesh_nested, DT)
      call state_unstr_exit_data(state_nested)
      call mesh_unstr_exit_data(mesh_nested)

      ! Cell 1 must differ between wall and nested because BC_NESTED
      ! injects a nonzero Flather flux whereas BC_WALL reflects.
      call check(error, &
                 abs(state_nested%barotropic%stage_cv(1) - state_wall%barotropic%stage_cv(1)) > 1.0e-10_wp, &
                 "BC_NESTED cell differs from BC_WALL baseline after one step")

      call state_unstr_destroy(state_wall)
      call mesh_unstr_destroy(mesh_wall)
      call state_unstr_destroy(state_nested)
      call mesh_unstr_destroy(mesh_nested)
   end subroutine test_nested_differs_from_wall

end module test_nested_unstr
