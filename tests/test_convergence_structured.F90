!! Observed convergence order of the structured HLL path on a smooth solution.
module test_convergence_structured
   !! The structured backend had **no refinement test at all** — the only two
   !! in the suite are `test_convergence` (structured, but dam-break and so
   !! shock-limited, asserting only rate >= 0.3) and `test_dam_break_unstr`
   !! (unstructured KNP). Neither verifies that the interior scheme reaches its
   !! design order, because neither uses a smooth solution.
   !!
   !! This runs the structured HLL + minmod-limited MUSCL + SSP-RK2 path on a
   !! smooth analytic solution and reads the observed order off the Phase 2
   !! harness.
   !!
   !! ## The oracle: fundamental standing mode in a closed basin
   !!
   !! Linearised shallow water about depth `H0` in `[0, L]` with walls:
   !!
   !!     eta(x,t) = A cos(kx) cos(wt),   u(x,t) = (A w / (k H0)) sin(kx) sin(wt)
   !!     k = pi / L,   w = k sqrt(g H0)
   !!
   !! Two properties make this the right case rather than merely a convenient
   !! one:
   !!
   !!   * `d(eta)/dx = 0` at both walls, so the mode is *compatible* with the
   !!     wall boundary condition and no boundary-layer error contaminates the
   !!     interior;
   !!   * `cos(pi x / L)` is **monotone** on `[0, L]` — its extrema sit exactly
   !!     ON the two walls, so the interior contains no smooth extremum for
   !!     minmod to clip. That matters: the campaign plan warns that a case
   !!     with an interior crest caps the observed rate near 1 and produces a
   !!     "defensible-but-misleading" number. A monotone profile avoids the
   !!     failure mode instead of tolerating it.
   !!
   !! ## Why the amplitude is so small
   !!
   !! The oracle is the LINEARISED solution, so the model's nonlinearity is an
   !! error source in its own right, of relative size O(A/H0). For the
   !! discretisation error to be what we actually measure, we need
   !! `A/H0 << C h^2`. At the finest level here `h = 0.25 m` and the observed
   !! relative error is ~1e-5, so `A/H0 = 1e-7` leaves roughly two orders of
   !! headroom. Push the amplitude up and the series flattens against the
   !! nonlinear floor — which the harness would expose as a collapsing final
   !! pairwise rate, exactly the plateau signature its self-check pins.
   !!
   !! Round-off is not a concern at that amplitude: `eta = H0 +- 1e-6` resolves
   !! the perturbation to ~1e-9 relative in double precision, four orders below
   !! the discretisation error being measured.
   !!
   !! ## Why dt is set explicitly rather than left to the CFL
   !!
   !! SSP-RK2 is 2nd order in time, so the temporal error must shrink at the
   !! same rate as the spatial one or the measured order is a blend of the two.
   !! Each level therefore takes exactly `2x` the steps of the previous one to
   !! the SAME final time, making `dt` proportional to `h` by construction and
   !! landing every level on an identical `t` for comparison. Letting the
   !! adaptive CFL choose would give slightly different final times per level
   !! and a contaminated series.
   !!
   !! ## Measured
   !!
   !!     nx     h        L1          L2         Linf
   !!     50   2.00   6.22e-08   6.56e-09   8.50e-10
   !!    100   1.00   1.63e-08   1.69e-09   2.21e-10
   !!    200   0.50   4.16e-09   4.31e-10   5.63e-11
   !!    400   0.25   1.05e-09   1.09e-10   1.42e-11
   !!
   !!     norm   pairwise rates          overall
   !!     L1     1.936  1.966  1.981     1.961
   !!     L2     1.953  1.974  1.985     1.971
   !!     Linf   1.942  1.974  1.984     1.967
   !!
   !! The structured HLL path reaches its design order. Two details are worth
   !! more than the headline number:
   !!
   !!   * the pairwise rates **increase** monotonically toward 2
   !!     (1.95 -> 1.97 -> 1.99), which is asymptotic convergence rather than a
   !!     lucky average — the opposite of the plateau signature;
   !!   * **Linf converges at 1.97 too**, not at 1 as it would if minmod were
   !!     clipping a smooth extremum. That is the direct evidence the monotone
   !!     profile did its job; a case with an interior crest would show L1 and
   !!     L2 near 2 while Linf collapsed toward 1, and reporting only L2 would
   !!     have hidden it.
   !!
   !! The 1.8 floor was fixed by the campaign plan BEFORE these numbers were
   !! measured, and is left where it is. It is not a threshold fitted to the
   !! result.
   !!
   !! ## Not duplicated on the unstructured backend
   !!
   !! `test_dam_break_unstr` already carries a KNP refinement test. It is
   !! shock-limited rather than smooth, so it verifies something weaker, but
   !! adding a second unstructured refinement case is not what this phase is
   !! for — the gap being closed here is the structured one.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, PI
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_boundary_types, only: bc_type_from_string
   use convergence_harness, only: error_norms_uniform, convergence_rates, &
                                  spacing_from_counts
   implicit none
   private

   public :: collect_convergence_structured_tests

   integer, parameter :: N_LEVELS = 4
   integer, parameter :: NX_LIST(N_LEVELS) = [50, 100, 200, 400]
   integer, parameter :: NY = 4
   integer, parameter :: NGHOST = 2
   real(wp), parameter :: LX = 100.0_wp
   real(wp), parameter :: H0 = 10.0_wp
   real(wp), parameter :: AMP = 1.0e-6_wp
      !! A/H0 = 1e-7 — see the amplitude note above.
   integer, parameter :: N_STEPS_BASE = 40
      !! Steps at the coarsest level; doubles per level so dt ~ h.

   real(wp), parameter :: MIN_RATE_OVERALL = 1.8_wp
      !! Design order is 2 for the interior scheme; 0.2 of slack absorbs the
      !! limiter/boundary contribution from the two wall-adjacent cells at each
      !! end, whose share of the norm shrinks under refinement.
   real(wp), parameter :: MIN_RATE_PAIR = 1.5_wp

contains

   subroutine collect_convergence_structured_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("structured_hll_reaches_design_order_l2", test_order_l2), &
                  new_unittest("structured_hll_reaches_design_order_l1", test_order_l1), &
                  new_unittest("structured_hll_error_decreases_monotonically", test_monotone), &
                  new_unittest("structured_convergence_fixture_is_smooth", test_fixture) &
                  ]
   end subroutine collect_convergence_structured_tests

   pure function omega() result(w)
      real(wp) :: w
      w = (PI/LX)*sqrt(GRAVITY*H0)
   end function omega

   pure function eta_exact(x, t) result(eta)
      !! Linearised fundamental standing mode, measured from the still level.
      real(wp), intent(in) :: x, t
      real(wp) :: eta
      eta = AMP*cos(PI*x/LX)*cos(omega()*t)
   end function eta_exact

   subroutine run_level(level, l1, l2, linf)
      !! Run one refinement level and return the error norms of the free
      !! surface against the analytic mode at the common final time.
      integer, intent(in) :: level
      real(wp), intent(out) :: l1, l2, linf

      type(state_t) :: state
      integer :: nx, n_steps, i, j, n, j_mid
      real(wp) :: dx, dt, t_end, t, x
      real(wp), allocatable :: numeric(:), exact(:)

      nx = NX_LIST(level)
      dx = LX/real(nx, wp)
      n_steps = N_STEPS_BASE*2**(level - 1)

      ! One eighth of a period: long enough that the solution has genuinely
      ! evolved (cos(wt) ~ 0.92, u well away from zero), short enough to stay
      ! comfortably linear.
      t_end = 0.125_wp*(2.0_wp*PI/omega())
      dt = t_end/real(n_steps, wp)

      call state%init(nx, NY, NGHOST, dx, dx, 0.45_wp)
      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            x = (real(i - NGHOST, wp) - 0.5_wp)*dx
            state%barotropic%b(i, j) = 0.0_wp
            state%barotropic%h(i, j) = H0 + eta_exact(x, 0.0_wp)
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
         end do
      end do

      call solver_enter_data(state)
      t = 0.0_wp
      do n = 1, n_steps
         ! dt_max caps the CFL-derived step. The CFL step here is ~0.09*h/2
         ! at every level, comfortably above `dt`, so `dt` is what is taken —
         ! and it scales exactly with h.
         call solver_step(state, dt_max=dt, time=t)
         t = t + state%dt
      end do
      call solver_exit_data(state)

      ! Sample the mid-row: the solution is 1-D in x, and the interior rows
      ! are identical, so one row carries the whole error.
      j_mid = NGHOST + NY/2
      allocate (numeric(nx), exact(nx))
      do i = 1, nx
         x = (real(i, wp) - 0.5_wp)*dx
         numeric(i) = state%barotropic%h(i + NGHOST, j_mid) - H0
         exact(i) = eta_exact(x, t)
      end do

      call error_norms_uniform(numeric, exact, dx, l1, l2, linf)

      deallocate (numeric, exact)
      call state%destroy()
   end subroutine run_level

   subroutine gather_series(h, l1, l2, linf)
      real(wp), intent(out) :: h(N_LEVELS), l1(N_LEVELS), l2(N_LEVELS), linf(N_LEVELS)
      integer :: level

      call spacing_from_counts(LX, NX_LIST, h)
      do level = 1, N_LEVELS
         call run_level(level, l1(level), l2(level), linf(level))
      end do
   end subroutine gather_series

   subroutine test_order_l2(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(N_LEVELS), l1(N_LEVELS), l2(N_LEVELS), linf(N_LEVELS)
      real(wp) :: p_pair(N_LEVELS - 1)
      real(wp) :: p_overall
      logical :: ok
      integer :: i

      call gather_series(h, l1, l2, linf)
      call convergence_rates(h, l2, p_pair, p_overall, ok)

      call check(error, ok, "harness rejected the L2 series")
      if (allocated(error)) return
      call check(error, p_overall >= MIN_RATE_OVERALL, &
                 "structured HLL L2 convergence below design order")
      if (allocated(error)) return
      do i = 1, N_LEVELS - 1
         call check(error, p_pair(i) >= MIN_RATE_PAIR, &
                    "a pairwise L2 rate fell below the floor — check for an error plateau")
         if (allocated(error)) return
      end do
   end subroutine test_order_l2

   subroutine test_order_l1(error)
      !! L1 as well as L2. minmod degrades the two wall-adjacent cells at each
      !! end, and the two norms weight that differently, so agreeing rates are
      !! evidence the order is a property of the interior scheme rather than an
      !! artefact of one norm.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(N_LEVELS), l1(N_LEVELS), l2(N_LEVELS), linf(N_LEVELS)
      real(wp) :: p_pair(N_LEVELS - 1)
      real(wp) :: p_overall
      logical :: ok

      call gather_series(h, l1, l2, linf)
      call convergence_rates(h, l1, p_pair, p_overall, ok)

      call check(error, ok, "harness rejected the L1 series")
      if (allocated(error)) return
      call check(error, p_overall >= MIN_RATE_OVERALL, &
                 "structured HLL L1 convergence below design order")
   end subroutine test_order_l1

   subroutine test_monotone(error)
      !! The error must actually shrink at every level. A series that stops
      !! improving has hit a floor, and a healthy OVERALL rate would still
      !! report success — the plateau failure mode the harness self-check pins.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(N_LEVELS), l1(N_LEVELS), l2(N_LEVELS), linf(N_LEVELS)
      integer :: i

      call gather_series(h, l1, l2, linf)

      do i = 1, N_LEVELS - 1
         call check(error, l2(i + 1) < l2(i), &
                    "L2 error did not decrease under refinement")
         if (allocated(error)) return
      end do
   end subroutine test_monotone

   subroutine test_fixture(error)
      !! Guard the premises the order claim rests on: the initial profile must
      !! be monotone across the interior (no extremum for minmod to clip), and
      !! the amplitude must be small enough that the linearised oracle is the
      !! right one.
      type(error_type), allocatable, intent(out) :: error

      integer :: i, nx
      real(wp) :: dx, x, prev, cur
      logical :: monotone

      nx = NX_LIST(1)
      dx = LX/real(nx, wp)

      monotone = .true.
      prev = eta_exact(0.5_wp*dx, 0.0_wp)
      do i = 2, nx
         x = (real(i, wp) - 0.5_wp)*dx
         cur = eta_exact(x, 0.0_wp)
         if (cur >= prev) monotone = .false.
         prev = cur
      end do

      call check(error, monotone, &
                 "fixture: initial profile is not monotone — minmod will clip an interior extremum")
      if (allocated(error)) return
      call check(error, AMP/H0 < 1.0e-6_wp, &
                 "fixture: amplitude too large for the linearised oracle to be the dominant error")
   end subroutine test_fixture

end module test_convergence_structured
