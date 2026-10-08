!! Unit tests for the Manning friction kernel
module test_friction_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, RHO_WATER
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
   use rki_kernel_friction_unstr, only: apply_manning_friction_unstr
   implicit none
   private

   public :: collect_friction_tests

   real(wp), parameter :: TOL = 1.0e-12_wp

contains

   subroutine collect_friction_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("zero_manning_no_op", test_zero_manning_no_op), &
                  new_unittest("zero_dt_no_op", test_zero_dt_no_op), &
                  new_unittest("zero_momentum_stays_zero", test_zero_momentum_stays_zero), &
                  new_unittest("dry_cell_untouched", test_dry_cell_untouched), &
                  new_unittest("hand_computed_semi_implicit_factor", test_hand_computed), &
                  new_unittest("friction_reduces_momentum", test_reduces_momentum), &
                  new_unittest("monotonic_decay_over_many_calls", test_monotonic_decay), &
                  new_unittest("wind_stress_hand_computed", test_wind_hand_computed), &
                  new_unittest("wind_stress_dry_cell_skipped", test_wind_dry_skipped), &
                  new_unittest("wind_zero_no_op", test_wind_zero_no_op) &
                  ]
   end subroutine collect_friction_tests

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

   subroutine run_friction(state, mesh, dt)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: dt

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call apply_manning_friction_unstr(state, mesh, dt)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_friction

   ! ----------------------------------------------------------------------
   ! Tests
   ! ----------------------------------------------------------------------

   subroutine test_zero_manning_no_op(error)
      !! Default manning_n = 0 -> friction kernel is a no-op even with
      !! finite momentum and a finite dt.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 1.0_wp
      state%barotropic%xmom_cv = [0.5_wp, -0.5_wp]
      state%barotropic%ymom_cv = [0.3_wp, 0.4_wp]
      ! manning_n stays at default 0.

      call run_friction(state, mesh, 0.1_wp)

      call check(error, abs(state%barotropic%xmom_cv(1) - 0.5_wp) < TOL, "xmom(1) untouched")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(2) + 0.5_wp) < TOL, "xmom(2) untouched")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(1) - 0.3_wp) < TOL, "ymom(1) untouched")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(2) - 0.4_wp) < TOL, "ymom(2) untouched")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_zero_manning_no_op

   subroutine test_zero_dt_no_op(error)
      !! dt = 0 -> kernel must not divide by anything, even with
      !! a real Manning coefficient.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 1.0_wp
      state%barotropic%xmom_cv = [0.5_wp, -0.5_wp]
      state%barotropic%ymom_cv = [0.3_wp, 0.4_wp]
      state%barotropic%manning_n = 0.025_wp

      call run_friction(state, mesh, 0.0_wp)

      call check(error, abs(state%barotropic%xmom_cv(1) - 0.5_wp) < TOL, "dt=0: xmom(1) untouched")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(2) - 0.4_wp) < TOL, "dt=0: ymom(2) untouched")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_zero_dt_no_op

   subroutine test_zero_momentum_stays_zero(error)
      !! Lake-at-rest case with friction enabled. |u| = 0 -> denom = 1
      !! -> momentum stays exactly at zero. Important sanity check that
      !! turning friction on doesn't somehow generate spurious motion.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 2.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      state%barotropic%manning_n = 0.03_wp

      call run_friction(state, mesh, 0.5_wp)

      call check(error, abs(state%barotropic%xmom_cv(1)) < TOL, "lake-at-rest xmom stays 0")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(2)) < TOL, "lake-at-rest xmom stays 0")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(1)) < TOL, "lake-at-rest ymom stays 0")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(2)) < TOL, "lake-at-rest ymom stays 0")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_zero_momentum_stays_zero

   subroutine test_dry_cell_untouched(error)
      !! Cell with depth below min_h must be skipped entirely (the 1/h
      !! singularity would otherwise nuke the run). Wet cell on the
      !! same mesh: friction must follow the analytic semi-implicit
      !! decay uh^{n+1} = uh^n / (1 + dt*c_f*|u|/h).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp), parameter :: N_MAN = 0.025_wp, DT_VAL = 0.1_wp
      real(wp), parameter :: H_WET = 1.0_wp, UH_WET = 0.5_wp
      real(wp) :: cf, denom, expected_wet

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Cell 1 dry (depth = 0), cell 2 wet (depth = 1).
      state%barotropic%bed_cv = [0.0_wp, 0.0_wp]
      state%barotropic%stage_cv = [0.0_wp, H_WET]
      state%barotropic%xmom_cv = [0.5_wp, UH_WET]   ! cell 1 garbage, should be untouched
      state%barotropic%ymom_cv = [0.0_wp, 0.0_wp]
      state%barotropic%manning_n = N_MAN

      call run_friction(state, mesh, DT_VAL)

      cf = GRAVITY*N_MAN*N_MAN/(H_WET**(1.0_wp/3.0_wp))
      denom = 1.0_wp + DT_VAL*cf*(UH_WET/H_WET)/H_WET    ! |u| = UH_WET/H_WET
      expected_wet = UH_WET/denom

      call check(error, abs(state%barotropic%xmom_cv(1) - 0.5_wp) < TOL, &
                 "dry cell xmom untouched (skipped)")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(2) - expected_wet) < 1.0e-13_wp, &
                 "wet cell xmom must match analytic Manning decay")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_dry_cell_untouched

   subroutine test_hand_computed(error)
      !! One wet cell, one momentum component, exact arithmetic.
      !!
      !! Setup: h = 1, n = 0.025, uh = 1, vh = 0, dt = 0.1.
      !!   u     = uh / h          = 1
      !!   |u|   = sqrt(1 + 0)     = 1
      !!   c_f   = g * n^2 / h^(1/3) = 9.80665 * 0.000625 / 1 = 0.006129...
      !!   denom = 1 + dt * c_f * |u| / h = 1 + 0.1 * 0.006129 * 1 / 1
      !!         = 1.000613
      !!   uh_new = 1 / 1.000613... ~= 0.999388
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: cf, denom, expected, n

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      n = 0.025_wp
      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = [1.0_wp, 1.0_wp]
      state%barotropic%xmom_cv = [1.0_wp, 0.0_wp]   ! cell 1 carries momentum, cell 2 still
      state%barotropic%ymom_cv = [0.0_wp, 0.0_wp]
      state%barotropic%manning_n = n

      call run_friction(state, mesh, 0.1_wp)

      ! Recompute the expected value here so the test is self-checking
      ! against the same algebra the kernel uses.
      cf = GRAVITY*n*n/(1.0_wp**(1.0_wp/3.0_wp))
      denom = 1.0_wp + 0.1_wp*cf*1.0_wp/1.0_wp     ! |u| = 1, h = 1
      expected = 1.0_wp/denom

      call check(error, abs(state%barotropic%xmom_cv(1) - expected) < 1.0e-13_wp, &
                 "semi-implicit factor: uh_new = uh / (1 + dt*cf*|u|/h)")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(2)) < TOL, &
                 "still cell stays still")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_hand_computed

   subroutine test_reduces_momentum(error)
      !! Strict-inequality contract: friction must reduce |momentum|,
      !! never increase it. With finite manning_n and finite dt and
      !! nonzero starting momentum, the result must be strictly smaller.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: initial_xmom, initial_ymom

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 2.0_wp
      state%barotropic%xmom_cv = 0.7_wp
      state%barotropic%ymom_cv = -0.3_wp
      state%barotropic%manning_n = 0.04_wp
      initial_xmom = 0.7_wp
      initial_ymom = -0.3_wp

      call run_friction(state, mesh, 0.2_wp)

      call check(error, state%barotropic%xmom_cv(1) < initial_xmom, &
                 "positive xmom should decrease")
      if (allocated(error)) return
      call check(error, state%barotropic%xmom_cv(1) > 0.0_wp, &
                 "positive xmom should not flip sign")
      if (allocated(error)) return
      call check(error, state%barotropic%ymom_cv(1) > initial_ymom, &
                 "negative ymom should increase toward 0")
      if (allocated(error)) return
      call check(error, state%barotropic%ymom_cv(1) < 0.0_wp, &
                 "negative ymom should not flip sign")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_reduces_momentum

   subroutine test_monotonic_decay(error)
      !! Repeated friction calls should drive momentum monotonically
      !! toward zero -- a spinning down lake.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: prev_xmom
      integer :: i

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 1.0_wp
      state%barotropic%xmom_cv = 1.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      state%barotropic%manning_n = 0.05_wp

      prev_xmom = state%barotropic%xmom_cv(1)
      do i = 1, 10
         call run_friction(state, mesh, 0.1_wp)
         call check(error, state%barotropic%xmom_cv(1) < prev_xmom, &
                    "xmom should decrease monotonically each call")
         if (allocated(error)) return
         call check(error, state%barotropic%xmom_cv(1) > 0.0_wp, &
                    "semi-implicit form should never overshoot zero")
         if (allocated(error)) return
         prev_xmom = state%barotropic%xmom_cv(1)
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_monotonic_decay

   subroutine test_wind_hand_computed(error)
      !! Wet cell, zero friction, nonzero wind. The kernel must add
      !! `dt * tau / RHO_WATER` to each momentum component on every
      !! wet cell. Hand-computed expected value.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: tau_x, tau_y, dt, expected_x, expected_y

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      tau_x = 0.2_wp
      tau_y = -0.1_wp
      dt = 0.5_wp

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 1.0_wp
      state%barotropic%xmom_cv = [0.1_wp, 0.0_wp]
      state%barotropic%ymom_cv = [0.0_wp, 0.05_wp]
      state%barotropic%manning_n = 0.0_wp
      state%barotropic%wind_tau_x = tau_x
      state%barotropic%wind_tau_y = tau_y

      call run_friction(state, mesh, dt)

      expected_x = 0.1_wp + dt*tau_x/RHO_WATER
      expected_y = 0.0_wp + dt*tau_y/RHO_WATER

      call check(error, abs(state%barotropic%xmom_cv(1) - expected_x) < TOL, &
                 "wind: xmom += dt*tau_x/rho")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(1) - expected_y) < TOL, &
                 "wind: ymom += dt*tau_y/rho")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(2) - dt*tau_x/RHO_WATER) < TOL, &
                 "wind: cell 2 picks up wind on initially-still water")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(2) - (0.05_wp + dt*tau_y/RHO_WATER)) < TOL, &
                 "wind: cell 2 ymom updated")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_wind_hand_computed

   subroutine test_wind_dry_skipped(error)
      !! Dry cell must not pick up wind momentum -- otherwise the next
      !! step would have nonzero hu/hv on a zero-depth cell, which
      !! creates infinite velocity.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Cell 1 dry, cell 2 wet.
      state%barotropic%bed_cv = [0.0_wp, 0.0_wp]
      state%barotropic%stage_cv = [0.0_wp, 1.0_wp]
      state%barotropic%xmom_cv = [0.0_wp, 0.0_wp]
      state%barotropic%ymom_cv = [0.0_wp, 0.0_wp]
      state%barotropic%wind_tau_x = 0.5_wp
      state%barotropic%wind_tau_y = 0.5_wp

      call run_friction(state, mesh, 0.1_wp)

      call check(error, abs(state%barotropic%xmom_cv(1)) < TOL, &
                 "dry cell xmom stays zero")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(1)) < TOL, &
                 "dry cell ymom stays zero")
      if (allocated(error)) return
      call check(error, state%barotropic%xmom_cv(2) > 0.0_wp, &
                 "wet cell receives wind kick")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_wind_dry_skipped

   subroutine test_wind_zero_no_op(error)
      !! tau_x = tau_y = 0 and manning_n = 0 -> kernel returns immediately
      !! and does not touch any state.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 1.0_wp
      state%barotropic%xmom_cv = [0.5_wp, -0.5_wp]
      state%barotropic%ymom_cv = [0.3_wp, 0.4_wp]
      ! manning_n, wind_tau_* all default to zero.

      call run_friction(state, mesh, 0.1_wp)

      call check(error, abs(state%barotropic%xmom_cv(1) - 0.5_wp) < TOL, "no-op xmom(1)")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(2) - 0.4_wp) < TOL, "no-op ymom(2)")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_wind_zero_no_op

end module test_friction_unstr
