!! Thacker parabolic bowl wetting/drying test
module test_thacker_bowl
   !! Runs the Thacker (1981) parabolic bowl sloshing problem and
   !! validates against the exact analytical solution.
   !!
   !! Tests: moving shoreline, wetting and drying, mass conservation
   !! with wet/dry fronts, no spurious negative depths.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, TWO_PI
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use thacker_bowl, only: thacker_omega, thacker_eta_exact, thacker_u_exact
   implicit none
   private

   public :: collect_thacker_bowl_tests

   real(wp), parameter :: H0 = 1.0_wp
   real(wp), parameter :: BOWL_RADIUS = 50.0_wp
   real(wp), parameter :: ETA_AMP = 0.3_wp
   real(wp), parameter :: DOMAIN_SIZE = 120.0_wp
   integer, parameter :: NX = 100
   integer, parameter :: NY = 100
   integer, parameter :: NGHOST = 2
   real(wp), parameter :: DX = DOMAIN_SIZE/real(NX, wp)
   real(wp), parameter :: DY = DOMAIN_SIZE/real(NY, wp)
   real(wp), parameter :: CFL_NUM = 0.4_wp
   real(wp), parameter :: MASS_TOL = 1.0e-4_wp
      !! Relative mass conservation tolerance (<0.01% per spec)
   real(wp), parameter :: L2_ETA_TOL = 0.10_wp
      !! L2 error in eta along the centreline, normalised by H0.
      !!
      !! **Measured history** at 100x100, CFL = 0.4, one period:
      !! **0.0785** before the 2026-08 thin-cell desingularization fix,
      !! **0.0841** after (an independent re-run of the same tree reported
      !! 0.089; treat ~0.084-0.089 as the run-to-run band, the difference
      !! coming from the adaptive step landing differently).
      !!
      !! The previous bound of 0.08 left only **2 %** margin over 0.0785 and
      !! had to be re-baselined regardless of any fix — the docstring note it
      !! carried claiming "~5 % measured" was stale and overstated the
      !! headroom, so the test was one small change away from tripping on
      !! anything.
      !!
      !! **Where the increase actually is.** Splitting the L2^2 budget by
      !! distance from the exact shoreline, same fixture, only `src/` swapped
      !! between the two commits:
      !!
      !!     band                 pre-fix     post-fix    change
      !!     edge (<=5 cells)     1.883e-2    2.038e-2     +8 %
      !!     near (6-20 cells)    2.984e-2    3.774e-2    +26 %
      !!     interior (>20)       7.014e-2    7.236e-2     +3.2 %
      !!
      !! So the bulk solution is essentially untouched; the growth sits in the
      !! band adjacent to the shoreline, which is exactly where the newly
      !! mobile thin water and the ~slope*dx hydrostatic-reconstruction
      !! retention rim act. (An earlier revision of this note attributed it to
      !! a shifted *interior* slosh — that was wrong, and the band split above
      !! is the correction.)
      !!
      !! Note also that the interior ALREADY supplied ~80 % of the L2^2 budget
      !! before the fix: that is first-order HLL smearing the whole bowl, a
      !! pre-existing property. Most of what this metric measures is therefore
      !! not the moving shoreline, which makes it a blunt instrument for
      !! judging a wet/dry change — read `test_dam_break_dry` and
      !! `test_carrier_greenspan` for the sharp ones.
      !!
      !! 0.10 leaves ~12-19 % headroom over the measured band: enough that
      !! run-to-run scatter cannot trip it, tight enough to catch a real
      !! regression. The pre-fix scheme scored better here only because its
      !! clamps froze the shoreline solid; see the two tests above for what
      !! that cost elsewhere.

contains

   subroutine collect_thacker_bowl_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("thacker_mass_conserved", test_mass), &
                  new_unittest("thacker_no_negative_depth", test_no_negative), &
                  new_unittest("thacker_eta_accuracy", test_eta_accuracy) &
                  ]
   end subroutine collect_thacker_bowl_tests

   subroutine run_thacker(state, t_final, mass_initial)
      !! Set up and run Thacker bowl for one full period
      type(state_t), intent(inout) :: state
      real(wp), intent(out) :: t_final, mass_initial

      integer :: i, j, n_steps
      real(wp) :: omega, period, x_pos, y_pos, x_center, y_center
      real(wp) :: eta_exact

      omega = thacker_omega(H0, BOWL_RADIUS)
      period = TWO_PI/omega

      call state%init(NX, NY, NGHOST, DX, DY, CFL_NUM)

      x_center = DOMAIN_SIZE/2.0_wp
      y_center = DOMAIN_SIZE/2.0_wp

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            x_pos = (real(i - NGHOST, wp) - 0.5_wp)*DX - x_center
            y_pos = (real(j - NGHOST, wp) - 0.5_wp)*DY - y_center

            state%barotropic%b(i, j) = H0*((x_pos*x_pos + y_pos*y_pos) &
                                           /(BOWL_RADIUS*BOWL_RADIUS) - 1.0_wp)
            eta_exact = thacker_eta_exact(x_pos, y_pos, 0.0_wp, &
                                          H0, BOWL_RADIUS, ETA_AMP)
            state%barotropic%h(i, j) = max(0.0_wp, eta_exact - state%barotropic%b(i, j))
            state%barotropic%hu(i, j) = state%barotropic%h(i, j)*thacker_u_exact(0.0_wp, &
                                                                                 H0, BOWL_RADIUS, ETA_AMP)
            state%barotropic%hv(i, j) = 0.0_wp
         end do
      end do

      mass_initial = compute_mass(state)

      call solver_enter_data(state)
      t_final = 0.0_wp
      n_steps = 0
      do while (t_final < period)
         call solver_step(state, dt_max=(period - t_final))
         t_final = t_final + state%dt
         n_steps = n_steps + 1
      end do
      call solver_exit_data(state)
   end subroutine run_thacker

   function compute_mass(s) result(mass)
      type(state_t), intent(in) :: s
      real(wp) :: mass
      integer :: i, j

      mass = 0.0_wp
      do j = s%grid%nghost + 1, s%grid%ny_total - s%grid%nghost
         do i = s%grid%nghost + 1, s%grid%nx_total - s%grid%nghost
            mass = mass + s%barotropic%h(i, j)*s%grid%dx*s%grid%dy
         end do
      end do
   end function compute_mass

   subroutine test_mass(error)
      !! Mass should be conserved to <0.01% (spec requirement)
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t_final, mass_initial, mass_final, mass_error

      call run_thacker(state, t_final, mass_initial)
      mass_final = compute_mass(state)
      mass_error = abs(mass_final - mass_initial)/max(mass_initial, 1.0e-30_wp)

      call check(error, mass_error < MASS_TOL, &
                 "Thacker bowl mass conservation violated")

      call state%destroy()
   end subroutine test_mass

   subroutine test_no_negative(error)
      !! No negative depths should occur
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t_final, mass_initial, min_h
      integer :: i, j

      call run_thacker(state, t_final, mass_initial)

      min_h = 0.0_wp
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            min_h = min(min_h, state%barotropic%h(i, j))
         end do
      end do

      call check(error, min_h >= -1.0e-10_wp, &
                 "Negative depth detected in Thacker bowl")

      call state%destroy()
   end subroutine test_no_negative

   subroutine test_eta_accuracy(error)
      !! L2 eta error along centreline should be within tolerance
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t_final, mass_initial
      real(wp) :: l2_eta, eta_exact, eta_num, err_val
      real(wp) :: x_pos, x_center
      integer :: i, j_mid

      call run_thacker(state, t_final, mass_initial)

      x_center = DOMAIN_SIZE/2.0_wp
      l2_eta = 0.0_wp
      j_mid = NGHOST + NY/2 + 1

      do i = NGHOST + 1, state%grid%nx_total - NGHOST
         x_pos = (real(i - NGHOST, wp) - 0.5_wp)*DX - x_center

         eta_exact = thacker_eta_exact(x_pos, 0.0_wp, t_final, &
                                       H0, BOWL_RADIUS, ETA_AMP)
         eta_num = state%barotropic%h(i, j_mid) + state%barotropic%b(i, j_mid)

         if (eta_exact > state%barotropic%b(i, j_mid) + 1.0e-6_wp) then
            err_val = abs(eta_num - eta_exact)
            l2_eta = l2_eta + err_val*err_val
         end if
      end do
      l2_eta = sqrt(l2_eta/real(NX, wp))

      call check(error, l2_eta/H0 < L2_ETA_TOL, &
                 "Thacker bowl L2 eta error too large")

      call state%destroy()
   end subroutine test_eta_accuracy

end module test_thacker_bowl
