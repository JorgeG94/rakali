!! Unit tests for the deterministic polynomial implementations in
!! `rki_safe_math`.  Tests call the `*_polynomial` entry points
!! directly so we validate the safe-mode path regardless of how
!! the module is compiled.
!!
!! Coverage:
!!   * `safe_exp_polynomial` — reference values, accuracy sweep,
!!     monotonicity, idempotence, overflow / underflow edges.
!!   * `safe_log_polynomial` — reference values, accuracy sweep,
!!     `log(exp(x)) == x` round-trip, x ≤ 0 sentinel.
!!   * `safe_sin_polynomial` / `safe_cos_polynomial` — reference
!!     values at quadrant boundaries, accuracy sweep, the
!!     Pythagorean identity `sin² + cos² == 1`.
!!   * `safe_pow_polynomial` — reference values, equivalence with
!!     `safe_exp(b * safe_log(a))`, zero-base behaviour.
module test_safe_math
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_safe_math, only: safe_exp_polynomial, safe_log_polynomial, &
                            safe_sin_polynomial, safe_cos_polynomial, &
                            safe_pow_polynomial
   implicit none
   private

   public :: collect_safe_math_tests

contains

   subroutine collect_safe_math_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("safe_exp_reference_values", test_exp_reference), &
                  new_unittest("safe_exp_accuracy_vs_intrinsic", test_exp_accuracy), &
                  new_unittest("safe_exp_monotonic", test_exp_monotonic), &
                  new_unittest("safe_exp_idempotent", test_exp_idempotent), &
                  new_unittest("safe_exp_edge_cases", test_exp_edges), &
                  new_unittest("safe_log_reference_values", test_log_reference), &
                  new_unittest("safe_log_accuracy_vs_intrinsic", test_log_accuracy), &
                  new_unittest("safe_log_inverse_of_exp", test_log_inverse_exp), &
                  new_unittest("safe_log_domain_guard", test_log_domain), &
                  new_unittest("safe_sin_reference_values", test_sin_reference), &
                  new_unittest("safe_cos_reference_values", test_cos_reference), &
                  new_unittest("safe_sin_cos_accuracy_vs_intrinsic", test_trig_accuracy), &
                  new_unittest("safe_sin_cos_pythagorean", test_pythagorean), &
                  new_unittest("safe_pow_reference_values", test_pow_reference), &
                  new_unittest("safe_pow_via_exp_log", test_pow_via_exp_log), &
                  new_unittest("safe_pow_zero_base", test_pow_zero_base) &
                  ]
   end subroutine collect_safe_math_tests

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_exp_reference(error)
      !! Landmark inputs.  `safe_exp_polynomial(0)` must be EXACTLY
      !! 1.0 — the algorithm short-circuits to `scale(1.0_wp, 0)`
      !! when `x = 0`.  Other landmarks match within a few ULPs of
      !! `exp(intrinsic)`.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: y, e_ref

      y = safe_exp_polynomial(0.0_wp)
      call check(error, y == 1.0_wp, &
                 "safe_exp(0) must be exactly 1.0")
      if (allocated(error)) return

      e_ref = exp(1.0_wp)
      y = safe_exp_polynomial(1.0_wp)
      call check(error, abs(y - e_ref)/e_ref < 1.0e-14_wp, &
                 "safe_exp(1) drifted from e")
      if (allocated(error)) return

      e_ref = exp(-1.0_wp)
      y = safe_exp_polynomial(-1.0_wp)
      call check(error, abs(y - e_ref)/e_ref < 1.0e-14_wp, &
                 "safe_exp(-1) drifted from 1/e")
      if (allocated(error)) return

      e_ref = exp(10.0_wp)
      y = safe_exp_polynomial(10.0_wp)
      call check(error, abs(y - e_ref)/e_ref < 1.0e-14_wp, &
                 "safe_exp(10) off intrinsic")
   end subroutine test_exp_reference

   subroutine test_exp_accuracy(error)
      !! Sweep [-20, 20] in 0.01 steps; the maximum relative error
      !! vs `exp(intrinsic)` must stay below 1e-14.  This bound is
      !! a few ULPs — accounts for the 14-term Taylor truncation
      !! plus libm's own rounding.
      type(error_type), allocatable, intent(out) :: error
      integer :: i, n
      real(wp) :: x, y, y_ref, err, max_err

      max_err = 0.0_wp
      n = 4001  ! -20 to 20 in 0.01 steps
      do i = 1, n
         x = -20.0_wp + real(i - 1, wp)*0.01_wp
         y = safe_exp_polynomial(x)
         y_ref = exp(x)
         err = abs(y - y_ref)/y_ref
         max_err = max(max_err, err)
      end do
      call check(error, max_err < 1.0e-14_wp, &
                 "safe_exp polynomial off intrinsic by > 1e-14")
   end subroutine test_exp_accuracy

   subroutine test_exp_monotonic(error)
      !! Dense sweep across [-50, 50]; `safe_exp_polynomial` must
      !! be strictly increasing.  Guards against polynomial-overflow
      !! glitches near the range-reduction boundary.
      type(error_type), allocatable, intent(out) :: error
      integer :: i, n
      real(wp) :: x_prev, x, y_prev, y
      logical :: ok

      ok = .true.
      n = 2000
      x_prev = -50.0_wp
      y_prev = safe_exp_polynomial(x_prev)
      do i = 2, n
         x = -50.0_wp + real(i - 1, wp)*0.05_wp
         y = safe_exp_polynomial(x)
         if (.not. (y > y_prev)) then
            ok = .false.
            exit
         end if
         x_prev = x
         y_prev = y
      end do
      call check(error, ok, "safe_exp not monotone-increasing")
   end subroutine test_exp_monotonic

   subroutine test_exp_idempotent(error)
      !! Calling the function twice on the same input must give
      !! byte-identical bits.  Self-determinism check — catches an
      !! accidental introduction of a non-deterministic intrinsic
      !! or a random-state read.
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      real(wp), parameter :: x_samples(7) = [-15.0_wp, -3.7_wp, -1.0_wp, 0.0_wp, &
                                             0.5_wp, 4.2_wp, 30.0_wp]
      real(wp) :: y1, y2
      logical :: ok

      ok = .true.
      do i = 1, size(x_samples)
         y1 = safe_exp_polynomial(x_samples(i))
         y2 = safe_exp_polynomial(x_samples(i))
         if (y1 /= y2) then
            ok = .false.
            exit
         end if
      end do
      call check(error, ok, "safe_exp not idempotent on the same input")
   end subroutine test_exp_idempotent

   subroutine test_exp_edges(error)
      !! Overflow / underflow guards.  Beyond the thresholds the
      !! function returns `huge` / `0` rather than letting the
      !! polynomial yield `inf` / `subnormal`.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: y_big, y_tiny, y_overflow

      ! Largest in-range input: 90% of the overflow threshold, so the
      ! bounds are precision-portable.  y_big = huge**0.9, always above
      ! sqrt(huge) and below huge in either single or double precision.
      y_big = safe_exp_polynomial(0.9_wp*log(huge(1.0_wp)))
      call check(error, y_big > sqrt(huge(y_big)) .and. y_big < huge(y_big), &
                 "safe_exp(0.9*log(huge)) outside the large-but-finite range")
      if (allocated(error)) return

      ! Past the overflow / underflow thresholds in any precision.
      y_overflow = safe_exp_polynomial(2.0_wp*log(huge(1.0_wp)))
      call check(error, y_overflow == huge(y_overflow), &
                 "safe_exp overflow guard didn't clamp to huge()")
      if (allocated(error)) return

      y_tiny = safe_exp_polynomial(2.0_wp*log(tiny(1.0_wp)))
      call check(error, y_tiny == 0.0_wp, &
                 "safe_exp underflow guard didn't clamp to 0")
   end subroutine test_exp_edges

   ! -----------------------------------------------------------------
   ! safe_log
   ! -----------------------------------------------------------------

   subroutine test_log_reference(error)
      !! Landmark inputs.  `safe_log(1) == 0` exact (the polynomial
      !! collapses to zero after the √2 split + atanh substitution).
      !! Other landmarks match within a few ULPs of `log`.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: y, ref

      y = safe_log_polynomial(1.0_wp)
      call check(error, y == 0.0_wp, "safe_log(1) must be exactly 0.0")
      if (allocated(error)) return

      ref = log(2.0_wp)
      y = safe_log_polynomial(2.0_wp)
      call check(error, abs(y - ref) < 1.0e-14_wp, &
                 "safe_log(2) off intrinsic")
      if (allocated(error)) return

      ref = log(exp(1.0_wp))
      y = safe_log_polynomial(exp(1.0_wp))
      call check(error, abs(y - ref) < 1.0e-14_wp, &
                 "safe_log(e) off 1.0")
      if (allocated(error)) return

      ref = log(1000.0_wp)
      y = safe_log_polynomial(1000.0_wp)
      call check(error, abs(y - ref)/abs(ref) < 1.0e-14_wp, &
                 "safe_log(1000) off intrinsic")
   end subroutine test_log_reference

   subroutine test_log_accuracy(error)
      !! Logarithmic sweep across [1e-10, 1e10].  Max relative error
      !! must stay below 1e-14.
      type(error_type), allocatable, intent(out) :: error
      integer :: i, n
      real(wp) :: x, y, y_ref, err, max_err

      max_err = 0.0_wp
      n = 4000
      do i = 1, n
         x = 10.0_wp**(-10.0_wp + 20.0_wp*real(i - 1, wp)/real(n - 1, wp))
         y = safe_log_polynomial(x)
         y_ref = log(x)
         if (abs(y_ref) > 1.0e-30_wp) then
            err = abs(y - y_ref)/abs(y_ref)
            max_err = max(max_err, err)
         end if
      end do
      call check(error, max_err < 1.0e-14_wp, &
                 "safe_log polynomial off intrinsic by > 1e-14")
   end subroutine test_log_accuracy

   subroutine test_log_inverse_exp(error)
      !! Round-trip identity `safe_log(safe_exp(x)) ≈ x`.  Allows
      !! ~5 ULP slack since two polynomial passes compose.
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      real(wp), parameter :: x_samples(8) = [-10.0_wp, -3.5_wp, -1.0_wp, -0.1_wp, &
                                             0.1_wp, 1.0_wp, 3.5_wp, 10.0_wp]
      real(wp) :: y, x_obs, max_err

      max_err = 0.0_wp
      do i = 1, size(x_samples)
         y = safe_exp_polynomial(x_samples(i))
         x_obs = safe_log_polynomial(y)
         max_err = max(max_err, abs(x_obs - x_samples(i)))
      end do
      call check(error, max_err < 1.0e-13_wp, &
                 "safe_log(safe_exp(x)) round-trip > 1e-13")
   end subroutine test_log_inverse_exp

   subroutine test_log_domain(error)
      !! `x <= 0` returns the `-huge()` sentinel.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: y

      y = safe_log_polynomial(0.0_wp)
      call check(error, y == -huge(y), "safe_log(0) didn't return -huge")
      if (allocated(error)) return

      y = safe_log_polynomial(-5.0_wp)
      call check(error, y == -huge(y), "safe_log(-5) didn't return -huge")
   end subroutine test_log_domain

   ! -----------------------------------------------------------------
   ! safe_sin / safe_cos
   ! -----------------------------------------------------------------

   subroutine test_sin_reference(error)
      !! `sin(0) == 0` exact; quadrant-boundary landmarks match the
      !! intrinsic.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: y
      real(wp), parameter :: PI = acos(-1.0_wp)

      y = safe_sin_polynomial(0.0_wp)
      call check(error, y == 0.0_wp, "safe_sin(0) must be exactly 0.0")
      if (allocated(error)) return

      y = safe_sin_polynomial(PI/2.0_wp)
      call check(error, abs(y - 1.0_wp) < 1.0e-14_wp, &
                 "safe_sin(π/2) off 1.0")
      if (allocated(error)) return

      y = safe_sin_polynomial(PI)
      call check(error, abs(y) < 1.0e-14_wp, &
                 "safe_sin(π) off 0.0")
      if (allocated(error)) return

      y = safe_sin_polynomial(3.0_wp*PI/2.0_wp)
      call check(error, abs(y + 1.0_wp) < 1.0e-14_wp, &
                 "safe_sin(3π/2) off -1.0")
      if (allocated(error)) return

      y = safe_sin_polynomial(2.0_wp*PI)
      call check(error, abs(y) < 1.0e-14_wp, &
                 "safe_sin(2π) off 0.0")
   end subroutine test_sin_reference

   subroutine test_cos_reference(error)
      !! `cos(0) == 1` exact; quadrant boundaries match the intrinsic.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: y
      real(wp), parameter :: PI = acos(-1.0_wp)

      y = safe_cos_polynomial(0.0_wp)
      call check(error, y == 1.0_wp, "safe_cos(0) must be exactly 1.0")
      if (allocated(error)) return

      y = safe_cos_polynomial(PI/2.0_wp)
      call check(error, abs(y) < 1.0e-14_wp, &
                 "safe_cos(π/2) off 0.0")
      if (allocated(error)) return

      y = safe_cos_polynomial(PI)
      call check(error, abs(y + 1.0_wp) < 1.0e-14_wp, &
                 "safe_cos(π) off -1.0")
      if (allocated(error)) return

      y = safe_cos_polynomial(2.0_wp*PI)
      call check(error, abs(y - 1.0_wp) < 1.0e-14_wp, &
                 "safe_cos(2π) off 1.0")
   end subroutine test_cos_reference

   subroutine test_trig_accuracy(error)
      !! Dense sweep [-100, 100].  Hits ~32 full periods so all four
      !! quadrants are exercised.  Max absolute error must stay
      !! below 1e-13 — the bound is wider than safe_exp's because
      !! Cody-Waite reduction loses ~1 ULP per |x|/2π revolutions.
      type(error_type), allocatable, intent(out) :: error
      integer :: i, n
      real(wp) :: x, max_err_sin, max_err_cos

      max_err_sin = 0.0_wp
      max_err_cos = 0.0_wp
      n = 8001
      do i = 1, n
         x = -100.0_wp + 200.0_wp*real(i - 1, wp)/real(n - 1, wp)
         max_err_sin = max(max_err_sin, abs(safe_sin_polynomial(x) - sin(x)))
         max_err_cos = max(max_err_cos, abs(safe_cos_polynomial(x) - cos(x)))
      end do
      call check(error, max_err_sin < 1.0e-13_wp, &
                 "safe_sin off intrinsic by > 1e-13")
      if (allocated(error)) return
      call check(error, max_err_cos < 1.0e-13_wp, &
                 "safe_cos off intrinsic by > 1e-13")
   end subroutine test_trig_accuracy

   subroutine test_pythagorean(error)
      !! `safe_sin²(x) + safe_cos²(x) == 1` to a few ULPs across a
      !! sweep that crosses multiple periods.  Catches sign errors
      !! in the quadrant-dispatch table.
      type(error_type), allocatable, intent(out) :: error
      integer :: i, n
      real(wp) :: x, s, c, residual, max_residual

      max_residual = 0.0_wp
      n = 2000
      do i = 1, n
         x = -20.0_wp + 40.0_wp*real(i - 1, wp)/real(n - 1, wp)
         s = safe_sin_polynomial(x)
         c = safe_cos_polynomial(x)
         residual = abs(s*s + c*c - 1.0_wp)
         max_residual = max(max_residual, residual)
      end do
      call check(error, max_residual < 1.0e-13_wp, &
                 "safe_sin² + safe_cos² off 1.0 by > 1e-13")
   end subroutine test_pythagorean

   ! -----------------------------------------------------------------
   ! safe_pow
   ! -----------------------------------------------------------------

   subroutine test_pow_reference(error)
      !! Landmark inputs: `pow(a, 0) == 1`, `pow(a, 1) == a`,
      !! `pow(2, 10) == 1024`, `pow(e, 1) == e`.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: y, ref

      y = safe_pow_polynomial(2.5_wp, 0.0_wp)
      call check(error, abs(y - 1.0_wp) < 1.0e-14_wp, &
                 "safe_pow(a, 0) must be 1.0")
      if (allocated(error)) return

      y = safe_pow_polynomial(7.0_wp, 1.0_wp)
      call check(error, abs(y - 7.0_wp)/7.0_wp < 1.0e-13_wp, &
                 "safe_pow(7, 1) off 7.0")
      if (allocated(error)) return

      y = safe_pow_polynomial(2.0_wp, 10.0_wp)
      call check(error, abs(y - 1024.0_wp)/1024.0_wp < 1.0e-13_wp, &
                 "safe_pow(2, 10) off 1024")
      if (allocated(error)) return

      ref = exp(1.0_wp)
      y = safe_pow_polynomial(ref, 1.0_wp)
      call check(error, abs(y - ref)/ref < 1.0e-13_wp, &
                 "safe_pow(e, 1) off e")
   end subroutine test_pow_reference

   subroutine test_pow_via_exp_log(error)
      !! `safe_pow(a, b)` must agree with `safe_exp(b*safe_log(a))`
      !! by construction.  This is the structural identity the
      !! implementation rests on — make sure no edit drifts it.
      type(error_type), allocatable, intent(out) :: error
      integer :: i, j
      real(wp), parameter :: a_samples(4) = [0.5_wp, 1.5_wp, 10.0_wp, 100.0_wp]
      real(wp), parameter :: b_samples(4) = [-2.0_wp, 0.5_wp, 1.7_wp, 3.0_wp]
      real(wp) :: y, ref, max_err

      max_err = 0.0_wp
      do i = 1, size(a_samples)
         do j = 1, size(b_samples)
            y = safe_pow_polynomial(a_samples(i), b_samples(j))
            ref = safe_exp_polynomial(b_samples(j)*safe_log_polynomial(a_samples(i)))
            if (abs(ref) > 1.0e-30_wp) then
               max_err = max(max_err, abs(y - ref)/abs(ref))
            end if
         end do
      end do
      call check(error, max_err < 1.0e-14_wp, &
                 "safe_pow not equal to exp(b*log(a)) by construction")
   end subroutine test_pow_via_exp_log

   subroutine test_pow_zero_base(error)
      !! `pow(0, +) = 0`, `pow(0, 0) = 1` (IEEE convention),
      !! `pow(0, -) = huge`.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: y

      y = safe_pow_polynomial(0.0_wp, 2.0_wp)
      call check(error, y == 0.0_wp, "safe_pow(0, +) must be 0")
      if (allocated(error)) return

      y = safe_pow_polynomial(0.0_wp, 0.0_wp)
      call check(error, y == 1.0_wp, "safe_pow(0, 0) must be 1 (IEEE rule)")
      if (allocated(error)) return

      y = safe_pow_polynomial(0.0_wp, -1.0_wp)
      call check(error, y == huge(y), "safe_pow(0, -) must clamp to huge")
   end subroutine test_pow_zero_base

end module test_safe_math
