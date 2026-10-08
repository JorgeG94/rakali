!! Self-check for the convergence harness.
module test_convergence_harness
   !! The harness is the instrument every later convergence phase reads its
   !! answer off, so it is calibrated here against sequences whose order is
   !! known by construction — no solver involved. If `test_convergence_structured`
   !! (Phase 7) ever reports a suspicious rate, this test says whether to
   !! suspect the scheme or the ruler.
   !!
   !! Coverage:
   !!   * an exact O(h^p) series must return exactly p, for p = 1, 2 and 1.5,
   !!     to round-off — this is an algebraic identity, not an approximation;
   !!   * non-integer refinement ratios and non-uniform ratios within one
   !!     series must work, since not every study refines by 2;
   !!   * a series that PLATEAUS at the fine end (the error-floor signature)
   !!     must show a healthy overall rate but a collapsing final pairwise
   !!     rate — the exact failure mode that motivates reporting both;
   !!   * degenerate input must set `ok = .false.` rather than return a
   !!     plausible number;
   !!   * the norms must reproduce hand-computable values, including the
   !!     weighted (unstructured) form.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use convergence_harness, only: error_norms_uniform, error_norms_weighted, &
                                  convergence_rates, spacing_from_counts
   implicit none
   private

   public :: collect_convergence_harness_tests

   real(wp), parameter :: EXACT_TOL = 1.0e-12_wp
      !! The rate of an exact power-law series is an identity; only round-off
      !! in log/exp separates it from the true exponent.

contains

   subroutine collect_convergence_harness_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("harness_recovers_second_order", test_second_order), &
                  new_unittest("harness_recovers_first_order", test_first_order), &
                  new_unittest("harness_recovers_fractional_order", test_fractional_order), &
                  new_unittest("harness_handles_nonuniform_refinement", test_nonuniform_ratio), &
                  new_unittest("harness_exposes_error_floor_plateau", test_plateau), &
                  new_unittest("harness_rejects_degenerate_input", test_degenerate), &
                  new_unittest("harness_norms_match_hand_computation", test_norms), &
                  new_unittest("harness_weighted_norms_match", test_weighted_norms), &
                  new_unittest("harness_spacing_from_counts", test_spacing) &
                  ]
   end subroutine collect_convergence_harness_tests

   ! ------------------------------------------------------------------
   ! Exact power-law series
   ! ------------------------------------------------------------------

   subroutine power_law_series(p, h, err)
      !! err = C * h**p with an arbitrary C, on a 2x-refined series.
      real(wp), intent(in) :: p
      real(wp), intent(out) :: h(4), err(4)
      real(wp), parameter :: C = 0.375_wp
      integer :: i

      do i = 1, 4
         h(i) = 1.0_wp/real(2**(i - 1), wp)
         err(i) = C*h(i)**p
      end do
   end subroutine power_law_series

   subroutine check_exact_order(error, p_want)
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: p_want

      real(wp) :: h(4), err(4), p_pair(3)
      real(wp) :: p_overall
      logical :: ok
      integer :: i

      call power_law_series(p_want, h, err)
      call convergence_rates(h, err, p_pair, p_overall, ok)

      call check(error, ok, "harness rejected a well-formed series")
      if (allocated(error)) return
      call check(error, abs(p_overall - p_want) < EXACT_TOL, &
                 "overall rate does not match the constructed order")
      if (allocated(error)) return
      do i = 1, 3
         call check(error, abs(p_pair(i) - p_want) < EXACT_TOL, &
                    "a pairwise rate does not match the constructed order")
         if (allocated(error)) return
      end do
   end subroutine check_exact_order

   subroutine test_second_order(error)
      type(error_type), allocatable, intent(out) :: error
      call check_exact_order(error, 2.0_wp)
   end subroutine test_second_order

   subroutine test_first_order(error)
      type(error_type), allocatable, intent(out) :: error
      call check_exact_order(error, 1.0_wp)
   end subroutine test_first_order

   subroutine test_fractional_order(error)
      !! 1.5 is the interesting one: a scheme that is formally 2nd order but
      !! limiter-clipped at extrema lands between 1 and 2, and the harness
      !! must report that honestly rather than rounding toward a design value.
      type(error_type), allocatable, intent(out) :: error
      call check_exact_order(error, 1.5_wp)
   end subroutine test_fractional_order

   subroutine test_nonuniform_ratio(error)
      !! Refinement by 3, then 1.5, then 2 — an uneven series. The rate is
      !! defined by the h ratio, not by any assumption of doubling.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(4), err(4), p_pair(3)
      real(wp) :: p_overall
      real(wp), parameter :: P_WANT = 2.0_wp
      logical :: ok
      integer :: i

      h = [1.0_wp, 1.0_wp/3.0_wp, 1.0_wp/4.5_wp, 1.0_wp/9.0_wp]
      do i = 1, 4
         err(i) = 0.2_wp*h(i)**P_WANT
      end do

      call convergence_rates(h, err, p_pair, p_overall, ok)

      call check(error, ok, "harness rejected a non-uniformly refined series")
      if (allocated(error)) return
      call check(error, abs(p_overall - P_WANT) < EXACT_TOL, &
                 "non-uniform refinement broke the overall rate")
      if (allocated(error)) return
      do i = 1, 3
         call check(error, abs(p_pair(i) - P_WANT) < EXACT_TOL, &
                    "non-uniform refinement broke a pairwise rate")
         if (allocated(error)) return
      end do
   end subroutine test_nonuniform_ratio

   ! ------------------------------------------------------------------
   ! The failure mode that motivates reporting both rates
   ! ------------------------------------------------------------------

   subroutine test_plateau(error)
      !! A clean 2nd-order series for the first three levels, then the error
      !! PLATEAUS — the signature of hitting a floor (round-off, an
      !! unconverged solver tolerance, a fixed time-step error) rather than
      !! continuing to converge.
      !!
      !! The overall rate stays comfortably above 1.8 and would pass a
      !! naive assertion. The final pairwise rate collapses to ~0. Asserting
      !! only the overall number would certify a scheme that has stopped
      !! converging, which is exactly why both are returned.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(4), err(4), p_pair(3)
      real(wp) :: p_overall
      logical :: ok
      integer :: i

      do i = 1, 4
         h(i) = 1.0_wp/real(2**(i - 1), wp)
      end do
      err(1) = 1.0e-2_wp
      err(2) = 2.5e-3_wp          ! clean 2nd order
      err(3) = 6.25e-4_wp         ! clean 2nd order
      err(4) = 6.24e-4_wp         ! plateau: error floor reached

      call convergence_rates(h, err, p_pair, p_overall, ok)

      call check(error, ok, "harness rejected the plateau series")
      if (allocated(error)) return
      call check(error, p_overall > 1.3_wp, &
                 "plateau series should still show a healthy OVERALL rate")
      if (allocated(error)) return
      call check(error, abs(p_pair(1) - 2.0_wp) < 1.0e-9_wp .and. &
                 abs(p_pair(2) - 2.0_wp) < 1.0e-9_wp, &
                 "the converging pairs should read 2nd order")
      if (allocated(error)) return
      call check(error, p_pair(3) < 0.05_wp, &
                 "the final pairwise rate must expose the plateau")
   end subroutine test_plateau

   ! ------------------------------------------------------------------
   ! Degenerate input
   ! ------------------------------------------------------------------

   subroutine test_degenerate(error)
      !! Every one of these must set ok = .false. — never return a number.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(4), err(4), p_pair(3)
      real(wp) :: h2(1), err2(1), p_pair2(1)
      real(wp) :: p_overall
      logical :: ok

      call power_law_series(2.0_wp, h, err)

      ! (a) an exactly-zero error: the scheme reproduced the oracle, so the
      !     order is UNDEFINED, not infinite.
      err(3) = 0.0_wp
      call convergence_rates(h, err, p_pair, p_overall, ok)
      call check(error,.not. ok, "zero error should be rejected, not rated")
      if (allocated(error)) return

      ! (b) a negative error norm — impossible, so it signals a broken caller.
      call power_law_series(2.0_wp, h, err)
      err(2) = -1.0_wp
      call convergence_rates(h, err, p_pair, p_overall, ok)
      call check(error,.not. ok, "negative error should be rejected")
      if (allocated(error)) return

      ! (c) h not strictly decreasing (levels supplied fine-first, or a
      !     duplicated resolution).
      call power_law_series(2.0_wp, h, err)
      h(3) = h(2)
      call convergence_rates(h, err, p_pair, p_overall, ok)
      call check(error,.not. ok, "non-decreasing h should be rejected")
      if (allocated(error)) return

      ! (d) a single level cannot define a rate.
      h2 = [1.0_wp]
      err2 = [1.0e-3_wp]
      call convergence_rates(h2, err2, p_pair2, p_overall, ok)
      call check(error,.not. ok, "a one-level series should be rejected")
   end subroutine test_degenerate

   ! ------------------------------------------------------------------
   ! Norms
   ! ------------------------------------------------------------------

   subroutine test_norms(error)
      !! Errors 1, -2, 3, -4 with dcell = 0.5:
      !!   L1   = (1+2+3+4)*0.5              = 5
      !!   L2   = sqrt((1+4+9+16)*0.5)       = sqrt(15)
      !!   Linf = 4
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: numeric(4), exact(4)
      real(wp) :: l1, l2, linf

      exact = [10.0_wp, 20.0_wp, 30.0_wp, 40.0_wp]
      numeric = [11.0_wp, 18.0_wp, 33.0_wp, 36.0_wp]

      call error_norms_uniform(numeric, exact, 0.5_wp, l1, l2, linf)

      call check(error, abs(l1 - 5.0_wp) < EXACT_TOL, "L1 norm incorrect")
      if (allocated(error)) return
      call check(error, abs(l2 - sqrt(15.0_wp)) < EXACT_TOL, "L2 norm incorrect")
      if (allocated(error)) return
      call check(error, abs(linf - 4.0_wp) < EXACT_TOL, "Linf norm incorrect")
   end subroutine test_norms

   subroutine test_weighted_norms(error)
      !! Same errors, per-cell weights 1, 2, 3, 4:
      !!   L1   = 1*1 + 2*2 + 3*3 + 4*4      = 30
      !!   L2   = sqrt(1*1 + 4*2 + 9*3 + 16*4) = sqrt(100) = 10
      !!   Linf = 4  (unweighted by definition)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: numeric(4), exact(4), weight(4)
      real(wp) :: l1, l2, linf

      exact = [10.0_wp, 20.0_wp, 30.0_wp, 40.0_wp]
      numeric = [11.0_wp, 18.0_wp, 33.0_wp, 36.0_wp]
      weight = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]

      call error_norms_weighted(numeric, exact, weight, l1, l2, linf)

      call check(error, abs(l1 - 30.0_wp) < EXACT_TOL, "weighted L1 incorrect")
      if (allocated(error)) return
      call check(error, abs(l2 - 10.0_wp) < EXACT_TOL, "weighted L2 incorrect")
      if (allocated(error)) return
      call check(error, abs(linf - 4.0_wp) < EXACT_TOL, "weighted Linf incorrect")
   end subroutine test_weighted_norms

   subroutine test_spacing(error)
      !! And that a counts-based series feeds straight in and reproduces the
      !! constructed order — the path a real consumer takes.
      type(error_type), allocatable, intent(out) :: error

      integer :: counts(4)
      real(wp) :: h(4), err(4), p_pair(3)
      real(wp) :: p_overall
      logical :: ok
      integer :: i

      counts = [50, 100, 200, 400]
      call spacing_from_counts(100.0_wp, counts, h)

      call check(error, abs(h(1) - 2.0_wp) < EXACT_TOL, "spacing_from_counts wrong")
      if (allocated(error)) return
      call check(error, abs(h(4) - 0.25_wp) < EXACT_TOL, "spacing_from_counts wrong")
      if (allocated(error)) return

      do i = 1, 4
         err(i) = 3.0_wp*h(i)**2
      end do
      call convergence_rates(h, err, p_pair, p_overall, ok)

      call check(error, ok .and. abs(p_overall - 2.0_wp) < EXACT_TOL, &
                 "counts-derived spacings did not reproduce the order")
   end subroutine test_spacing

end module test_convergence_harness
