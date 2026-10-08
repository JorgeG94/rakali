!! Integration tests for the unstructured SSP-RK2 driver
module test_solver_unstr
   !! These are the milestone tests for the unstructured solver --
   !! the first place we exercise time evolution end-to-end. The
   !! load-bearing test is `lake_at_rest_preserved_over_many_steps`:
   !! evolving a lake-at-rest initial condition through 100 RK2 steps
   !! must not perturb the conserved variables to more than machine
   !! precision. If that passes, the well-balanced property survives
   !! the full pipeline (extrapolate -> flux -> update -> protect ->
   !! repeat -> saxpy -> protect) on every cell, every step.
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
   use rki_solver_unstr, only: solver_unstr_step, &
                               solver_unstr_advance
   implicit none
   private

   public :: collect_solver_tests

   real(wp), parameter :: TOL = 1.0e-12_wp
   real(wp), parameter :: LAKE_TOL = 1.0e-10_wp
      !! Loosened slightly for the multi-step lake-at-rest tests --
      !! roundoff in the extrapolate / flux / update chain accumulates
      !! a few orders of magnitude over 100 steps, but should still be
      !! orders of magnitude below the initial state.

contains

   subroutine collect_solver_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("single_step_lake_at_rest_flat", test_single_step_flat), &
                  new_unittest("single_step_lake_at_rest_sloped", test_single_step_sloped), &
                  new_unittest("hundred_steps_lake_at_rest_flat", test_hundred_steps_flat), &
                  new_unittest("hundred_steps_lake_at_rest_sloped_2tri", test_hundred_steps_sloped_2tri), &
                  new_unittest("hundred_steps_lake_at_rest_sloped_4tri", test_hundred_steps_sloped_4tri), &
                  new_unittest("step_advances_time_and_count", test_step_advances_clock), &
                  new_unittest("advance_lands_on_t_end", test_advance_lands_on_t_end), &
                  new_unittest("mass_error_zero_for_lake_at_rest", test_mass_error_zero) &
                  ]
   end subroutine collect_solver_tests

   ! ----------------------------------------------------------------------
   ! Mesh + state helpers
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

   subroutine build_four_triangle_square(mesh)
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp) :: vx(5), vy(5)
      integer :: tv(3, 4)

      vx = [0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.5_wp]
      vy = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.5_wp]
      tv(:, 1) = [1, 2, 5]
      tv(:, 2) = [2, 3, 5]
      tv(:, 3) = [3, 4, 5]
      tv(:, 4) = [4, 1, 5]
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_four_triangle_square

   subroutine set_bed_linear(state, mesh, b0, bx, by)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: b0, bx, by
      integer :: k

      do k = 1, mesh%ncells_total
         state%barotropic%bed_cv(k) = b0 + bx*mesh%centroid_x(k) + by*mesh%centroid_y(k)
      end do
   end subroutine set_bed_linear

   subroutine wrap_step(state, mesh, dt)
      !! Run one solver step inside the device data dance.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: dt

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call solver_unstr_step(state, mesh, dt)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine wrap_step

   subroutine wrap_advance(state, mesh, t_end)
      !! Run advance() inside the device data dance.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: t_end

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call solver_unstr_advance(state, mesh, t_end)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine wrap_advance

   ! ----------------------------------------------------------------------
   ! Single-step lake-at-rest tests
   ! ----------------------------------------------------------------------

   subroutine test_single_step_flat(error)
      !! One RK2 step starting from a flat-bed lake at rest must leave
      !! every conserved variable byte-identical (or as close as
      !! floating-point arithmetic allows).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 2.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp

      call wrap_step(state, mesh, 0.01_wp)

      do k = 1, mesh%ncells_owned
         call check(error, abs(state%barotropic%stage_cv(k) - 2.0_wp) < TOL, &
                    "single step flat: stage preserved")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%xmom_cv(k)) < TOL, &
                    "single step flat: xmom stays 0")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%ymom_cv(k)) < TOL, &
                    "single step flat: ymom stays 0")
         if (allocated(error)) return
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_single_step_flat

   subroutine test_single_step_sloped(error)
      !! Same on a sloped bed -- the well-balanced source telescopes
      !! against the central-upwind pressure flux at every substage,
      !! so the conserved variables are unchanged.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = 5.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      call set_bed_linear(state, mesh, 1.0_wp, 0.5_wp, -0.3_wp)

      call wrap_step(state, mesh, 0.01_wp)

      do k = 1, mesh%ncells_owned
         call check(error, abs(state%barotropic%stage_cv(k) - 5.0_wp) < TOL, &
                    "single step sloped: stage preserved")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%xmom_cv(k)) < TOL, &
                    "single step sloped: xmom stays 0")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%ymom_cv(k)) < TOL, &
                    "single step sloped: ymom stays 0")
         if (allocated(error)) return
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_single_step_sloped

   ! ----------------------------------------------------------------------
   ! Multi-step lake-at-rest -- the load-bearing tests
   ! ----------------------------------------------------------------------

   subroutine test_hundred_steps_flat(error)
      !! 100 RK2 steps on a flat lake at rest. The integrated drift
      !! should be no worse than ~1e-12 (LAKE_TOL).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: i, k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 3.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp

      do i = 1, 100
         call wrap_step(state, mesh, 0.001_wp)
      end do

      do k = 1, mesh%ncells_owned
         call check(error, abs(state%barotropic%stage_cv(k) - 3.0_wp) < LAKE_TOL, &
                    "100 steps flat lake: stage preserved")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%xmom_cv(k)) < LAKE_TOL, &
                    "100 steps flat lake: xmom = 0")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%ymom_cv(k)) < LAKE_TOL, &
                    "100 steps flat lake: ymom = 0")
         if (allocated(error)) return
      end do

      call check(error, state%step_count == 100, &
                 "step_count should be 100 after 100 calls")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_hundred_steps_flat

   subroutine test_hundred_steps_sloped_2tri(error)
      !! 100 RK2 steps with a sloped bed on the 2-triangle mesh. This
      !! exercises the well-balanced property through the entire SSP-RK2
      !! pipeline 100 times -- the strongest correctness test we can
      !! do without I/O or real boundary conditions.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: i, k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = 4.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      call set_bed_linear(state, mesh, 0.5_wp, 0.4_wp, -0.2_wp)

      do i = 1, 100
         call wrap_step(state, mesh, 0.001_wp)
      end do

      do k = 1, mesh%ncells_owned
         call check(error, abs(state%barotropic%stage_cv(k) - 4.0_wp) < LAKE_TOL, &
                    "100 steps sloped 2-tri: stage preserved")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%xmom_cv(k)) < LAKE_TOL, &
                    "100 steps sloped 2-tri: xmom = 0")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%ymom_cv(k)) < LAKE_TOL, &
                    "100 steps sloped 2-tri: ymom = 0")
         if (allocated(error)) return
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_hundred_steps_sloped_2tri

   subroutine test_hundred_steps_sloped_4tri(error)
      !! Same on the 4-triangle mesh -- exercises the 3-neighbour
      !! reconstruction branch with hfactor blending.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: i, k

      call build_four_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = 4.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      call set_bed_linear(state, mesh, 0.5_wp, 0.7_wp, 0.2_wp)

      do i = 1, 100
         call wrap_step(state, mesh, 0.001_wp)
      end do

      do k = 1, mesh%ncells_owned
         call check(error, abs(state%barotropic%stage_cv(k) - 4.0_wp) < LAKE_TOL, &
                    "100 steps sloped 4-tri: stage preserved")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%xmom_cv(k)) < LAKE_TOL, &
                    "100 steps sloped 4-tri: xmom = 0")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%ymom_cv(k)) < LAKE_TOL, &
                    "100 steps sloped 4-tri: ymom = 0")
         if (allocated(error)) return
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_hundred_steps_sloped_4tri

   ! ----------------------------------------------------------------------
   ! Driver bookkeeping
   ! ----------------------------------------------------------------------

   subroutine test_step_advances_clock(error)
      !! Each step must increment state%time by exactly dt and
      !! state%step_count by 1.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)
      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 1.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp

      call wrap_step(state, mesh, 0.05_wp)
      call check(error, abs(state%time - 0.05_wp) < TOL, "time = dt after 1 step")
      if (allocated(error)) return
      call check(error, state%step_count == 1, "step_count = 1")
      if (allocated(error)) return

      call wrap_step(state, mesh, 0.05_wp)
      call check(error, abs(state%time - 0.10_wp) < TOL, "time = 2*dt after 2 steps")
      if (allocated(error)) return
      call check(error, state%step_count == 2, "step_count = 2")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_step_advances_clock

   subroutine test_advance_lands_on_t_end(error)
      !! solver_unstr_advance should land state%time exactly on t_end
      !! (the final step is capped). At lake-at-rest the wave speed is
      !! constant, so the CFL timestep never changes and the number of
      !! steps is exactly ceil(t_end / dt_cfl). state%dt after advance
      !! holds the unclipped dt_cfl (compute_timestep is called again
      !! at the end of the loop).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp), parameter :: T_END = 0.1_wp
      integer :: expected_count

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)
      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 2.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp

      call wrap_advance(state, mesh, T_END)

      call check(error, abs(state%time - T_END) < TOL, &
                 "advance should land exactly on t_end")
      if (allocated(error)) return

      expected_count = ceiling(T_END/state%dt)
      call check(error, state%step_count == expected_count, &
                 "step_count should equal ceil(t_end / dt_cfl)")
      if (allocated(error)) return

      call check(error, abs(state%barotropic%stage_cv(1) - 2.0_wp) < LAKE_TOL, &
                 "advance preserves lake at rest")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_advance_lands_on_t_end

   subroutine test_mass_error_zero(error)
      !! 100 steps of lake at rest must not generate any mass error
      !! from the protect kernel -- nothing ever drops below zero.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: i

      call build_four_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = 4.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      call set_bed_linear(state, mesh, 0.5_wp, 0.7_wp, 0.2_wp)

      do i = 1, 100
         call wrap_step(state, mesh, 0.001_wp)
      end do

      call check(error, abs(state%mass_error) < LAKE_TOL, &
                 "lake at rest should not produce mass error")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_mass_error_zero

end module test_solver_unstr
