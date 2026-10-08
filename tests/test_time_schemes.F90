!! Unit tests for the selectable explicit time integrator (time_scheme).
module test_time_schemes
   !! Verifies the euler/rk2/rk3/rk43 menu wired into the 2D barotropic
   !! coastal solver:
   !!
   !!   * temporal-order self-convergence (Richardson three-level): a
   !!     smooth small-amplitude standing wave in a closed flat-bottom
   !!     basin, integrated at dt / dt/2 / dt/4 to a fixed final time;
   !!     the order p = log2(||u(dt)-u(dt/2)|| / ||u(dt/2)-u(dt/4)||)
   !!     must recover 1 / 2 / 3 / 3.  The clamp/Froude limiter never
   !!     fires in this gentle regime so the truncation order is clean.
   !!   * wet/dry positivity: each scheme keeps h >= 0 driving a
   !!     dam-break into a dry bed (the SSP guarantee).
   !!   * ssp_coeff: rk43 takes ~2x the CFL-limited timestep of rk2 for
   !!     the same quiescent state (the C=2 stability budget).
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_timestepper, only: TIME_EULER, TIME_RK2, TIME_RK3, TIME_RK43
   implicit none
   private

   public :: collect_time_schemes_tests

   ! Standing-wave basin geometry (shared across the convergence runs).
   integer, parameter :: NX = 64, NY = 4, NGHOST = 2
   real(wp), parameter :: DX = 100.0_wp, DY = 100.0_wp
   real(wp), parameter :: H0 = 10.0_wp          !! mean depth (m)
   real(wp), parameter :: AMP = 0.01_wp         !! free-surface amplitude (m)

contains

   subroutine collect_time_schemes_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("euler_first_order", test_order_euler), &
                  new_unittest("rk2_second_order", test_order_rk2), &
                  new_unittest("rk3_third_order", test_order_rk3), &
                  new_unittest("rk43_third_order", test_order_rk43), &
                  new_unittest("positivity_all_schemes", test_positivity), &
                  new_unittest("rk43_doubles_cfl_dt", test_rk43_ssp_dt) &
                  ]
   end subroutine collect_time_schemes_tests

   subroutine run_standing_wave(scheme, dt, n_steps, h_save)
      !! Integrate the smooth standing wave with a fixed timestep and
      !! return the interior free-surface field.
      integer, intent(in) :: scheme
      real(wp), intent(in) :: dt
      integer, intent(in) :: n_steps
      real(wp), intent(out) :: h_save(NX, NY)

      type(state_t) :: s
      integer :: i, j, n
      real(wp) :: xc, kwave

      call s%init(NX, NY, NGHOST, DX, DY, 0.45_wp)
      s%time_scheme = scheme
      s%dt_fixed = dt

      ! Fundamental cosine mode, fits the wall (reflecting) basin.
      kwave = acos(-1.0_wp)/(real(NX, wp)*DX)
      s%barotropic%b = 0.0_wp
      s%barotropic%hu = 0.0_wp
      s%barotropic%hv = 0.0_wp
      do j = 1, s%grid%ny_total
         do i = 1, s%grid%nx_total
            xc = (real(i - s%grid%nghost, wp) - 0.5_wp)*DX
            s%barotropic%h(i, j) = H0 + AMP*cos(kwave*xc)
         end do
      end do

      call solver_enter_data(s)
      do n = 1, n_steps
         call solver_step(s)
      end do
      call solver_exit_data(s)

      do j = 1, NY
         do i = 1, NX
            h_save(i, j) = s%barotropic%h(i + NGHOST, j + NGHOST)
         end do
      end do

      call s%destroy()
   end subroutine run_standing_wave

   function l2_diff(a, b) result(d)
      !! Discrete L2 norm of the difference of two interior fields.
      real(wp), intent(in) :: a(NX, NY), b(NX, NY)
      real(wp) :: d
      integer :: i, j

      d = 0.0_wp
      do j = 1, NY
         do i = 1, NX
            d = d + (a(i, j) - b(i, j))**2
         end do
      end do
      d = sqrt(d/real(NX*NY, wp))
   end function l2_diff

   subroutine check_order(error, scheme, expected, tol, dt1, n1)
      !! Three-level Richardson order estimate for one scheme.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: scheme
      real(wp), intent(in) :: expected, tol
      real(wp), intent(in) :: dt1
         !! Coarsest timestep.  Chosen per scheme so the temporal
         !! truncation error stays well above the minmod-limiter noise
         !! floor at the wave crests (where the limiter clips smooth
         !! extrema to first order) — otherwise the observed slope of the
         !! higher-order schemes is contaminated by that fixed floor.
      integer, intent(in) :: n1
         !! Step count at the coarsest dt (final time T = n1*dt1).

      real(wp) :: h1(NX, NY), h2(NX, NY), h4(NX, NY)
      real(wp) :: d12, d24, order

      ! Same final time T = n1*dt1 at all three resolutions.
      call run_standing_wave(scheme, dt1, n1, h1)
      call run_standing_wave(scheme, 0.5_wp*dt1, 2*n1, h2)
      call run_standing_wave(scheme, 0.25_wp*dt1, 4*n1, h4)

      d12 = l2_diff(h1, h2)
      d24 = l2_diff(h2, h4)

      call check(error, d24 > 0.0_wp, "convergence differences must be nonzero")
      if (allocated(error)) return

      order = log(d12/d24)/log(2.0_wp)
      write (*, '(a,i0,3(a,es12.4))') "  [scheme ", scheme, "] d12=", d12, &
         " d24=", d24, " order=", order
      call check(error, abs(order - expected) < tol, &
                 "observed temporal order off target")
   end subroutine check_order

   subroutine test_order_euler(error)
      type(error_type), allocatable, intent(out) :: error
      call check_order(error, TIME_EULER, 1.0_wp, 0.25_wp, 1.5_wp, 30)
   end subroutine test_order_euler

   subroutine test_order_rk2(error)
      type(error_type), allocatable, intent(out) :: error
      call check_order(error, TIME_RK2, 2.0_wp, 0.35_wp, 1.5_wp, 30)
   end subroutine test_order_rk2

   subroutine test_order_rk3(error)
      type(error_type), allocatable, intent(out) :: error
      call check_order(error, TIME_RK3, 3.0_wp, 0.4_wp, 4.0_wp, 12)
   end subroutine test_order_rk3

   subroutine test_order_rk43(error)
      type(error_type), allocatable, intent(out) :: error
      call check_order(error, TIME_RK43, 3.0_wp, 0.4_wp, 5.0_wp, 10)
   end subroutine test_order_rk43

   subroutine test_positivity(error)
      !! Every scheme keeps h >= 0 driving a dam break into a dry bed.
      type(error_type), allocatable, intent(out) :: error
      integer :: schemes(4), si, i, j, n
      type(state_t) :: s
      real(wp) :: min_h

      schemes = [TIME_EULER, TIME_RK2, TIME_RK3, TIME_RK43]

      do si = 1, 4
         call s%init(40, 8, NGHOST, 1.0_wp, 1.0_wp, 0.30_wp)
         s%time_scheme = schemes(si)

         s%barotropic%b = 0.0_wp
         s%barotropic%hu = 0.0_wp
         s%barotropic%hv = 0.0_wp
         s%barotropic%h = 0.0_wp
         do j = 1, s%grid%ny_total
            do i = 1, s%grid%nghost + s%grid%nx_phys/2
               s%barotropic%h(i, j) = 2.0_wp
            end do
         end do

         call solver_enter_data(s)
         do n = 1, 30
            call solver_step(s)
         end do
         call solver_exit_data(s)

         min_h = 0.0_wp
         do j = s%grid%nghost + 1, s%grid%ny_total - s%grid%nghost
            do i = s%grid%nghost + 1, s%grid%nx_total - s%grid%nghost
               min_h = min(min_h, s%barotropic%h(i, j))
            end do
         end do

         call s%destroy()

         call check(error, min_h >= 0.0_wp, &
                    "scheme produced negative depth at a wet/dry front")
         if (allocated(error)) return
      end do
   end subroutine test_positivity

   subroutine test_rk43_ssp_dt(error)
      !! rk43 (C=2) should select ~2x the CFL-limited dt that rk2 picks
      !! for the same quiescent state.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s_rk2, s_rk43
      real(wp) :: dt_rk2, dt_rk43

      call s_rk2%init(20, 20, NGHOST, 50.0_wp, 50.0_wp, 0.40_wp)
      s_rk2%time_scheme = TIME_RK2
      s_rk2%barotropic%h = 5.0_wp
      s_rk2%barotropic%hu = 0.0_wp
      s_rk2%barotropic%hv = 0.0_wp
      s_rk2%barotropic%b = 0.0_wp
      call solver_enter_data(s_rk2)
      call solver_step(s_rk2)
      call solver_exit_data(s_rk2)
      dt_rk2 = s_rk2%dt
      call s_rk2%destroy()

      call s_rk43%init(20, 20, NGHOST, 50.0_wp, 50.0_wp, 0.40_wp)
      s_rk43%time_scheme = TIME_RK43
      s_rk43%barotropic%h = 5.0_wp
      s_rk43%barotropic%hu = 0.0_wp
      s_rk43%barotropic%hv = 0.0_wp
      s_rk43%barotropic%b = 0.0_wp
      call solver_enter_data(s_rk43)
      call solver_step(s_rk43)
      call solver_exit_data(s_rk43)
      dt_rk43 = s_rk43%dt
      call s_rk43%destroy()

      call check(error, abs(dt_rk43 - 2.0_wp*dt_rk2) < 1.0e-9_wp*dt_rk2, &
                 "rk43 dt should be 2x the rk2 CFL dt (ssp_coeff=2)")
   end subroutine test_rk43_ssp_dt

end module test_time_schemes
