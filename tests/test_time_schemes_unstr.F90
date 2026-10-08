!! Unit tests for the selectable explicit time integrator on the
!! unstructured (KNP central-upwind) 2D barotropic path.
module test_time_schemes_unstr
   !! Mirrors test_time_schemes for the unstructured solver:
   !!
   !!   * well-balanced: a lake at rest stays at rest under every scheme
   !!     (each stage is well-balanced; convex blends of equal states are
   !!     fixed points) — the strongest cheap correctness gate.
   !!   * temporal-order self-convergence (Richardson three-level) on a
   !!     smooth standing wave in a regular-triangulated channel — recovers
   !!     1 / 2 / 3 / 3 for euler / rk2 / rk3 / rk43.
   !!   * wet/dry positivity: each scheme keeps h >= 0 under a stage bump.
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
   use rki_solver_unstr, only: solver_unstr_step
   use rki_mesh_gen_unstr, only: gen_regular_triangulation
   use rki_timestepper, only: TIME_EULER, TIME_RK2, TIME_RK3, TIME_RK43
   implicit none
   private

   public :: collect_time_schemes_unstr_tests

   ! Standing-wave channel geometry.
   integer, parameter :: NX = 96, NY = 2
   real(wp), parameter :: LX = 19200.0_wp, LY = 400.0_wp
   real(wp), parameter :: H0 = 10.0_wp          !! mean depth (m)
   real(wp), parameter :: AMP = 0.01_wp         !! free-surface amplitude (m)
   real(wp), parameter :: LAKE_TOL = 1.0e-10_wp

contains

   subroutine collect_time_schemes_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("lake_at_rest_all_schemes", test_lake_all_schemes), &
                  new_unittest("euler_first_order", test_order_euler), &
                  new_unittest("rk2_second_order", test_order_rk2), &
                  new_unittest("rk3_third_order", test_order_rk3), &
                  new_unittest("rk43_third_order", test_order_rk43), &
                  new_unittest("positivity_all_schemes", test_positivity) &
                  ]
   end subroutine collect_time_schemes_unstr_tests

   subroutine build_channel(mesh)
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :)

      call gen_regular_triangulation(NX, NY, LX, LY, vx, vy, tv)
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_channel

   subroutine run_standing_wave(scheme, dt, n_steps, h_out)
      !! Integrate the smooth standing wave with a fixed timestep; return
      !! the owned-cell free-surface field (enter/exit data once).
      integer, intent(in) :: scheme
      real(wp), intent(in) :: dt
      integer, intent(in) :: n_steps
      real(wp), allocatable, intent(out) :: h_out(:)

      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k, n
      real(wp) :: kwave

      call build_channel(mesh)
      call state_unstr_init(state, mesh)
      state%time_scheme = scheme
      ! First-order (piecewise-constant) reconstruction: the KNP flux of
      ! cell averages is a smooth function of the state, so the temporal
      ! truncation is the ONLY error that varies with dt — no limiter
      ! noise floor at the wave crests.  Spatial order is irrelevant to a
      ! temporal-order measurement; the production (2nd-order) path is
      ! exercised by the lake-at-rest + positivity tests.
      state%extrapolate_velocity_second_order = .false.

      kwave = acos(-1.0_wp)/LX
      state%barotropic%bed_cv = -H0
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      do k = 1, mesh%ncells_total
         state%barotropic%stage_cv(k) = AMP*cos(kwave*mesh%centroid_x(k))
      end do

      allocate (h_out(mesh%ncells_owned))

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      do n = 1, n_steps
         call solver_unstr_step(state, mesh, dt)
      end do
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      do k = 1, mesh%ncells_owned
         h_out(k) = state%barotropic%stage_cv(k)
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine run_standing_wave

   function l2_diff(a, b) result(d)
      real(wp), intent(in) :: a(:), b(:)
      real(wp) :: d
      integer :: k

      d = 0.0_wp
      do k = 1, size(a)
         d = d + (a(k) - b(k))**2
      end do
      d = sqrt(d/real(size(a), wp))
   end function l2_diff

   subroutine check_order(error, scheme, expected, tol, dt1, n1)
      !! Three-level Richardson order estimate for one scheme.  dt1/n1 are
      !! chosen per scheme so the temporal truncation stays above the
      !! KNP-reconstruction limiter floor (same rationale as the
      !! structured test).
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: scheme
      real(wp), intent(in) :: expected, tol, dt1
      integer, intent(in) :: n1

      real(wp), allocatable :: h1(:), h2(:), h4(:)
      real(wp) :: d12, d24, order

      call run_standing_wave(scheme, dt1, n1, h1)
      call run_standing_wave(scheme, 0.5_wp*dt1, 2*n1, h2)
      call run_standing_wave(scheme, 0.25_wp*dt1, 4*n1, h4)

      d12 = l2_diff(h1, h2)
      d24 = l2_diff(h2, h4)

      call check(error, d24 > 0.0_wp, "convergence differences must be nonzero")
      if (allocated(error)) return

      order = log(d12/d24)/log(2.0_wp)
      write (*, '(a,i0,3(a,es12.4))') "  [unstr scheme ", scheme, "] d12=", d12, &
         " d24=", d24, " order=", order
      call check(error, abs(order - expected) < tol, &
                 "observed temporal order off target")
   end subroutine check_order

   subroutine test_order_euler(error)
      type(error_type), allocatable, intent(out) :: error
      call check_order(error, TIME_EULER, 1.0_wp, 0.25_wp, 1.0_wp, 24)
   end subroutine test_order_euler

   subroutine test_order_rk2(error)
      type(error_type), allocatable, intent(out) :: error
      call check_order(error, TIME_RK2, 2.0_wp, 0.3_wp, 1.0_wp, 24)
   end subroutine test_order_rk2

   subroutine test_order_rk3(error)
      type(error_type), allocatable, intent(out) :: error
      call check_order(error, TIME_RK3, 3.0_wp, 0.5_wp, 1.5_wp, 16)
   end subroutine test_order_rk3

   subroutine test_order_rk43(error)
      type(error_type), allocatable, intent(out) :: error
      call check_order(error, TIME_RK43, 3.0_wp, 0.5_wp, 6.0_wp, 6)
   end subroutine test_order_rk43

   subroutine test_lake_all_schemes(error)
      !! Lake at rest on a sloped bed must stay at rest under every
      !! scheme (well-balanced through the generic N-stage path).
      type(error_type), allocatable, intent(out) :: error
      integer :: schemes(4), si, k, n
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      schemes = [TIME_EULER, TIME_RK2, TIME_RK3, TIME_RK43]

      do si = 1, 4
         call build_channel(mesh)
         call state_unstr_init(state, mesh)
         state%time_scheme = schemes(si)

         state%barotropic%xmom_cv = 0.0_wp
         state%barotropic%ymom_cv = 0.0_wp
         do k = 1, mesh%ncells_total
            state%barotropic%bed_cv(k) = -H0 + 1.0e-4_wp*mesh%centroid_x(k)/LX
            state%barotropic%stage_cv(k) = 1.0_wp
         end do

         call mesh_unstr_enter_data(mesh)
         call state_unstr_enter_data(state)
         do n = 1, 50
            call solver_unstr_step(state, mesh, 0.2_wp)
         end do
         call state_unstr_exit_data(state)
         call mesh_unstr_exit_data(mesh)

         do k = 1, mesh%ncells_owned
            call check(error, abs(state%barotropic%stage_cv(k) - 1.0_wp) < LAKE_TOL, &
                       "lake at rest: stage drifted under a scheme")
            if (allocated(error)) exit
            call check(error, abs(state%barotropic%xmom_cv(k)) < LAKE_TOL, &
                       "lake at rest: xmom nonzero under a scheme")
            if (allocated(error)) exit
            call check(error, abs(state%barotropic%ymom_cv(k)) < LAKE_TOL, &
                       "lake at rest: ymom nonzero under a scheme")
            if (allocated(error)) exit
         end do

         call state_unstr_destroy(state)
         call mesh_unstr_destroy(mesh)
         if (allocated(error)) return
      end do
   end subroutine test_lake_all_schemes

   subroutine test_positivity(error)
      !! A localized free-surface bump draining over the channel must keep
      !! h >= 0 (no negative mass error) under every scheme.
      type(error_type), allocatable, intent(out) :: error
      integer :: schemes(4), si, k, n
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      schemes = [TIME_EULER, TIME_RK2, TIME_RK3, TIME_RK43]

      do si = 1, 4
         call build_channel(mesh)
         call state_unstr_init(state, mesh)
         state%time_scheme = schemes(si)

         state%barotropic%bed_cv = -1.0_wp
         state%barotropic%xmom_cv = 0.0_wp
         state%barotropic%ymom_cv = 0.0_wp
         do k = 1, mesh%ncells_total
            if (mesh%centroid_x(k) < 0.5_wp*LX) then
               state%barotropic%stage_cv(k) = 1.0_wp
            else
               state%barotropic%stage_cv(k) = 0.0_wp
            end if
         end do

         call mesh_unstr_enter_data(mesh)
         call state_unstr_enter_data(state)
         do n = 1, 40
            call solver_unstr_step(state, mesh, 1.0_wp)
         end do
         call state_unstr_exit_data(state)
         call mesh_unstr_exit_data(mesh)

         do k = 1, mesh%ncells_owned
            call check(error, &
                       state%barotropic%stage_cv(k) - state%barotropic%bed_cv(k) >= -LAKE_TOL, &
                       "scheme produced negative depth at a wet/dry front")
            if (allocated(error)) exit
         end do

         call state_unstr_destroy(state)
         call mesh_unstr_destroy(mesh)
         if (allocated(error)) return
      end do
   end subroutine test_positivity

end module test_time_schemes_unstr
