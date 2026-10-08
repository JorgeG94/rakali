!! Grid convergence study for the dam break problem
module test_convergence
   !! Runs the dam break at 4 resolutions (200, 400, 800, 1600 cells)
   !! and verifies that the L2 error decreases with resolution.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use ritter_dambreak, only: exact_riemann_dambreak, exact_riemann_sample
   implicit none
   private

   public :: collect_convergence_tests

   integer, parameter :: NY = 4
   integer, parameter :: NGHOST = 2
   real(wp), parameter :: CFL_NUM = 0.45_wp
   real(wp), parameter :: T_END = 1.0_wp
   real(wp), parameter :: H_LEFT = 2.0_wp
   real(wp), parameter :: H_RIGHT = 1.0_wp
   real(wp), parameter :: DOMAIN_LENGTH = 200.0_wp
   real(wp), parameter :: MIN_CONVERGENCE_RATE = 0.3_wp
      !! Minimum acceptable L2 convergence rate per level pair.
      !! Per-pair L2 rate on dam break is shock-limited and can dip
      !! below 0.5 at the fine end (800 → 1600), so a per-pair floor
      !! of 0.3 is the right "nothing has stalled" check. Stronger
      !! asymptotic behaviour is enforced separately via the overall
      !! 200 → 1600 rate (see test_rate below).

contains

   subroutine collect_convergence_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("convergence_error_decreases", test_error_decreases), &
                  new_unittest("convergence_rate_acceptable", test_rate) &
                  ]
   end subroutine collect_convergence_tests

   subroutine run_dam_break_l2(nx, dx, h_star, u_star, l2_err)
      !! Run a dam break and return L2 error
      integer, intent(in) :: nx
      real(wp), intent(in) :: dx, h_star, u_star
      real(wp), intent(out) :: l2_err

      type(state_t) :: state
      real(wp) :: dy, x_dam, t_current
      real(wp) :: x, h_exact, u_exact, err_val
      integer :: ii, jj, j_mid, mid_x

      dy = dx
      x_dam = DOMAIN_LENGTH/2.0_wp

      call state%init(nx, NY, NGHOST, dx, dy, CFL_NUM)

      mid_x = NGHOST + nx/2
      do jj = 1, state%grid%ny_total
         do ii = 1, state%grid%nx_total
            state%barotropic%b(ii, jj) = 0.0_wp
            state%barotropic%hu(ii, jj) = 0.0_wp
            state%barotropic%hv(ii, jj) = 0.0_wp
            if (ii <= mid_x) then
               state%barotropic%h(ii, jj) = H_LEFT
            else
               state%barotropic%h(ii, jj) = H_RIGHT
            end if
         end do
      end do

      call solver_enter_data(state)
      t_current = 0.0_wp
      do while (t_current < T_END)
         call solver_step(state, dt_max=(T_END - t_current))
         t_current = t_current + state%dt
      end do
      call solver_exit_data(state)

      j_mid = NGHOST + NY/2 + 1
      l2_err = 0.0_wp
      do ii = NGHOST + 1, state%grid%nx_total - NGHOST
         x = (real(ii - NGHOST, wp) - 0.5_wp)*dx
         call exact_riemann_sample(H_LEFT, H_RIGHT, h_star, u_star, &
                                   x_dam, t_current, x, h_exact, u_exact)
         err_val = state%barotropic%h(ii, j_mid) - h_exact
         l2_err = l2_err + err_val*err_val
      end do
      l2_err = sqrt(l2_err/real(nx, wp))

      call state%destroy()
   end subroutine run_dam_break_l2

   subroutine test_error_decreases(error)
      !! Finest grid error should be smaller than coarsest
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h_star, u_star, l2_coarse, l2_fine
      real(wp) :: dx_coarse, dx_fine

      call exact_riemann_dambreak(H_LEFT, H_RIGHT, h_star, u_star)

      dx_coarse = DOMAIN_LENGTH/200.0_wp
      dx_fine = DOMAIN_LENGTH/1600.0_wp

      call run_dam_break_l2(200, dx_coarse, h_star, u_star, l2_coarse)
      call run_dam_break_l2(1600, dx_fine, h_star, u_star, l2_fine)

      call check(error, l2_fine < l2_coarse, &
                 "Error should decrease with finer grid")
   end subroutine test_error_decreases

   subroutine test_rate(error)
      !! Convergence-rate sanity on the dam-break L2 series.
      !!
      !! Per-pair L2 rate is shock-limited and can dip below 0.5 at
      !! the fine end, but:
      !!   * no pair should stall below MIN_CONVERGENCE_RATE (= 0.3),
      !!   * errors must decrease monotonically with refinement, and
      !!   * the overall 200 → 1600 rate must be >= 0.5 (asymptotic).
      !!
      !! Earlier versions of this test checked only the per-pair rate
      !! at 0.3, which admits an isolated good pair masking a stalled
      !! neighbour. The three checks together pin asymptotic
      !! behaviour without being brittle to single-level noise.
      use, intrinsic :: iso_fortran_env, only: error_unit
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N_LEVELS = 4
      integer :: nx_list(N_LEVELS)
      real(wp) :: dx_list(N_LEVELS), l2_list(N_LEVELS)
      real(wp) :: h_star, u_star, rate, rate_total
      integer :: level

      nx_list = [200, 400, 800, 1600]

      call exact_riemann_dambreak(H_LEFT, H_RIGHT, h_star, u_star)

      do level = 1, N_LEVELS
         dx_list(level) = DOMAIN_LENGTH/real(nx_list(level), wp)
         call run_dam_break_l2(nx_list(level), dx_list(level), &
                               h_star, u_star, l2_list(level))
         write (error_unit, '(a,i0,a,es12.4)') &
            '  level ', level, ': L2 = ', l2_list(level)
      end do

      ! Per-pair: nothing should stall below MIN_CONVERGENCE_RATE.
      do level = 2, N_LEVELS
         rate = log(l2_list(level - 1)/l2_list(level))/ &
                log(dx_list(level - 1)/dx_list(level))
         write (error_unit, '(a,i0,a,i0,a,f6.3)') &
            '  pair ', nx_list(level - 1), ' -> ', nx_list(level), &
            ': rate = ', rate
         call check(error, rate > MIN_CONVERGENCE_RATE, &
                    "Per-pair convergence rate below MIN_CONVERGENCE_RATE")
         if (allocated(error)) return
      end do

      ! Errors must monotonically decrease.
      do level = 2, N_LEVELS
         call check(error, l2_list(level) < l2_list(level - 1), &
                    "L2 error must decrease with refinement")
         if (allocated(error)) return
      end do

      ! Overall rate 200 -> 1600 must be >= 0.5 (asymptotic expectation).
      rate_total = log(l2_list(1)/l2_list(N_LEVELS)) &
                   /log(dx_list(1)/dx_list(N_LEVELS))
      write (error_unit, '(a,f6.3)') '  total 200 -> 1600 rate = ', rate_total
      call check(error, rate_total >= 0.5_wp, &
                 "Overall L2 rate 200 -> 1600 should be >= 0.5")
   end subroutine test_rate

end module test_convergence
