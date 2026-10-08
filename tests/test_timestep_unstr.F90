!! Unit tests for the CFL timestep kernel
module test_timestep_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY
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
   use rki_kernel_extrapolate_unstr, only: extrapolate_second_order_edge_unstr
   use rki_kernel_flux_unstr, only: compute_fluxes_central_unstr
   use rki_kernel_timestep_unstr, only: compute_timestep_unstr, TIMESTEP_LARGE
   implicit none
   private

   public :: collect_timestep_tests

   real(wp), parameter :: TOL = 1.0e-12_wp

contains

   subroutine collect_timestep_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("lake_at_rest_dt_matches_cfl_formula", test_lake_at_rest_dt), &
                  new_unittest("dt_scales_linearly_with_cfl", test_dt_scales_with_cfl), &
                  new_unittest("zero_max_speed_does_not_constrain", test_zero_max_speed), &
                  new_unittest("deeper_water_smaller_dt", test_deeper_water_smaller_dt) &
                  ]
   end subroutine collect_timestep_tests

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

   subroutine run_flux_then_timestep(state, mesh)
      !! Populate max_speed via the extrapolation+flux pipeline, then
      !! compute the CFL timestep. Wraps the data dance.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)

      call extrapolate_second_order_edge_unstr(state, mesh)
      call compute_fluxes_central_unstr(state, mesh)
      call compute_timestep_unstr(state, mesh)

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_flux_then_timestep

   ! ----------------------------------------------------------------------
   ! Tests
   ! ----------------------------------------------------------------------

   subroutine test_lake_at_rest_dt(error)
      !! Lake at rest with depth h. Each cell sees max_speed = sqrt(g*h)
      !! (we already verified this in test_flux_unstr). The CFL timestep
      !! is then dt = cfl * min(radius)/sqrt(g*h). With identical
      !! triangles and identical depth, all cells produce the same
      !! candidate so dt = cfl * radius / sqrt(g*h).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: expected_dt, c

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 3.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      state%cfl = 0.5_wp

      call run_flux_then_timestep(state, mesh)

      ! CFL radius = min centroid-to-edge-midpoint distance (ANUGA `radii`).
      ! For each unit right-isoceles triangle (legs 1) the nearest edge
      ! midpoint is the hypotenuse's, at distance sqrt(2)/6 (= 1/3 of the
      ! hypotenuse median). Both triangles are identical so min == that.
      c = sqrt(GRAVITY*3.0_wp)
      expected_dt = state%cfl*(sqrt(2.0_wp)/6.0_wp)/c

      call check(error, abs(state%dt - expected_dt) < 1.0e-12_wp, &
                 "dt should match cfl * radius / sqrt(g*h)")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_lake_at_rest_dt

   subroutine test_dt_scales_with_cfl(error)
      !! Doubling the CFL number should exactly double the timestep.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: dt_half, dt_full

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 1.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp

      state%cfl = 0.25_wp
      call run_flux_then_timestep(state, mesh)
      dt_half = state%dt

      state%cfl = 0.5_wp
      call run_flux_then_timestep(state, mesh)
      dt_full = state%dt

      call check(error, abs(dt_full - 2.0_wp*dt_half) < 1.0e-12_wp, &
                 "dt should double when cfl doubles")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_dt_scales_with_cfl

   subroutine test_zero_max_speed(error)
      !! If we manually zero max_speed (skipping the flux pass), the
      !! reduction should leave dt at its initial "infinity". This is
      !! a contract check on the epsilon guard.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! No flux call -- max_speed is zero everywhere from init.
      ! Run timestep directly.
      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call compute_timestep_unstr(state, mesh)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      ! With no flux pass, every cell has max_speed = 0 and the reduction
      ! never overwrites the TIMESTEP_LARGE initializer, so dt lands on
      ! exactly cfl * TIMESTEP_LARGE.
      call check(error, abs(state%dt - state%cfl*TIMESTEP_LARGE) < &
                 1.0e-80_wp*TIMESTEP_LARGE, &
                 "no characteristics -> dt = cfl * TIMESTEP_LARGE exactly")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_zero_max_speed

   subroutine test_deeper_water_smaller_dt(error)
      !! Deeper water -> faster gravity waves -> smaller CFL timestep.
      !! The ratio should be exactly sqrt(h_shallow / h_deep).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: dt_shallow, dt_deep, expected_ratio, actual_ratio

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)
      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      state%cfl = 0.5_wp

      state%barotropic%stage_cv = 1.0_wp
      call run_flux_then_timestep(state, mesh)
      dt_shallow = state%dt

      state%barotropic%stage_cv = 4.0_wp
      call run_flux_then_timestep(state, mesh)
      dt_deep = state%dt

      ! dt ~ 1/sqrt(g*h), so dt_deep / dt_shallow = sqrt(1/4) = 0.5
      expected_ratio = sqrt(1.0_wp/4.0_wp)
      actual_ratio = dt_deep/dt_shallow
      call check(error, abs(actual_ratio - expected_ratio) < 1.0e-12_wp, &
                 "dt ratio should be sqrt(h_shallow/h_deep)")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_deeper_water_smaller_dt

end module test_timestep_unstr
