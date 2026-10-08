!! Reproducibility-safe wrappers around transcendentals and `**` — a
!! finished library with no production consumer today.
module rki_safe_math
   !! Single point of control the math functions that break bitwise
   !! reproducibility across compilers, GPU vendors, and optimization
   !! levels would route through IF adopted: `exp`, `log`, `sin`,
   !! `cos`, `sqrt`, `pow`.
   !!
   !! **Status (PR-8):** the `safe_*` wrappers (`safe_exp`, `safe_log`,
   !! `safe_sin`, `safe_cos`, `safe_sqrt`, `safe_pow`) are plain
   !! elemental inlines around the bare intrinsic — zero overhead,
   !! codegen identical to writing `exp(x)` directly.  This module
   !! used to carry a second "safe" mode, selected by the
   !! `RAKALI_BITWISE_REPRO` build option, that dispatched the
   !! wrappers to the `*_polynomial` implementations below instead.
   !! That option was **removed in PR-8**: grep across `src/`/`app/`
   !! found zero call sites of any `safe_*` wrapper from a production
   !! kernel (one stale comment, no calls) — the option changed
   !! nothing at all, so it was a lie about a capability nothing used.
   !!
   !! What SURVIVES this PR: every `*_polynomial` function and its
   !! coefficients are correct, unit-tested (`tests/test_safe_math.F90`,
   !! unconditionally compiled — no `#ifdef`), and are the raw material
   !! for a future bitwise-reproducibility PR.  This module is a
   !! **library with no production consumer today** — adopting `safe_*`
   !! in kernels (which would first require re-adding a build-time or
   !! run-time dispatch) is future work, not a promise this module
   !! currently keeps.
   !!
   !! Would-be convention (not currently enforced by any build path):
   !! a kernel wanting repro-safety would call the `safe_*` wrapper
   !! instead of the bare intrinsic.  For integer exponents, `a**2` /
   !! `a**3` etc. compile to multiplications and are already
   !! deterministic — don't wrap those.
   !!
   !! Coverage: `safe_exp`, `safe_log`, `safe_sin`, `safe_cos` each
   !! have real polynomial paths (Cody-Waite reduction + Horner
   !! Taylor).  `safe_sqrt` is IEEE-correctly-rounded by mandate so
   !! the intrinsic path is already deterministic.  `safe_pow`'s
   !! polynomial path is unused by `safe_pow` itself (which is the
   !! bare intrinsic `**`) — it exists for a future
   !! `safe_exp(b * safe_log(a))` decomposition.
   use, intrinsic :: iso_fortran_env, only: real64
   use rki_constants, only: wp
   implicit none
   private

   public :: safe_exp, safe_log, safe_sin, safe_cos, safe_sqrt, safe_pow
   public :: safe_exp_polynomial, safe_log_polynomial
   public :: safe_sin_polynomial, safe_cos_polynomial
   public :: safe_pow_polynomial

   ! =================================================================
   ! Shared constants
   ! =================================================================

   ! ---- Cody-Waite ln(2) split (also used by safe_log) ----
   ! First ~33 bits of ln(2) in LN2_HI; remainder in LN2_LO.  Both
   ! double-precision representable; splitting keeps Cody-Waite
   ! `x - k*LN2_HI - k*LN2_LO` accurate to full precision.
   real(wp), parameter :: LN2_HI = 0.6931471805599452862_wp
   real(wp), parameter :: LN2_LO = 2.319046813846299558e-17_wp
   real(wp), parameter :: INV_LN2 = 1.4426950408889634074_wp

   ! ---- Cody-Waite π/2 split (for safe_sin / safe_cos) ----
   ! Standard fdlibm `__kernel_sin` constants.
   real(wp), parameter :: PI_OVER_2_HI = 1.5707963267948965580_wp
   real(wp), parameter :: PI_OVER_2_LO = 6.123233995736766036e-17_wp
   real(wp), parameter :: INV_PI_OVER_2 = 0.6366197723675813431_wp

   ! ---- Taylor coefficients ----
   ! `EXP_C(k) = 1/k!` for the exp Taylor series.  14 terms gives
   ! r¹⁵/15! ≈ 2e-20 at |r| = ln(2)/2 — well below double-precision
   ! epsilon.
   real(wp), parameter :: EXP_C0 = 1.0_wp
   real(wp), parameter :: EXP_C1 = 1.0_wp
   real(wp), parameter :: EXP_C2 = 1.0_wp/2.0_wp
   real(wp), parameter :: EXP_C3 = 1.0_wp/6.0_wp
   real(wp), parameter :: EXP_C4 = 1.0_wp/24.0_wp
   real(wp), parameter :: EXP_C5 = 1.0_wp/120.0_wp
   real(wp), parameter :: EXP_C6 = 1.0_wp/720.0_wp
   real(wp), parameter :: EXP_C7 = 1.0_wp/5040.0_wp
   real(wp), parameter :: EXP_C8 = 1.0_wp/40320.0_wp
   real(wp), parameter :: EXP_C9 = 1.0_wp/362880.0_wp
   real(wp), parameter :: EXP_C10 = 1.0_wp/3628800.0_wp
   real(wp), parameter :: EXP_C11 = 1.0_wp/39916800.0_wp
   real(wp), parameter :: EXP_C12 = 1.0_wp/479001600.0_wp
   real(wp), parameter :: EXP_C13 = 1.0_wp/6227020800.0_wp

   ! `LOG_Cn = 1/(2n+1)` for the atanh-form log series:
   ! `log((1+f)/(1-f)) = 2 * (f + f³/3 + f⁵/5 + ... )`.  Substitute
   ! `u = f²` and Horner over `u` with these coefficients.  10 terms
   ! on |f| ≤ (√2-1)/(√2+1) ≈ 0.172 gives < 1e-17 error.
   real(wp), parameter :: LOG_C0 = 1.0_wp
   real(wp), parameter :: LOG_C1 = 1.0_wp/3.0_wp
   real(wp), parameter :: LOG_C2 = 1.0_wp/5.0_wp
   real(wp), parameter :: LOG_C3 = 1.0_wp/7.0_wp
   real(wp), parameter :: LOG_C4 = 1.0_wp/9.0_wp
   real(wp), parameter :: LOG_C5 = 1.0_wp/11.0_wp
   real(wp), parameter :: LOG_C6 = 1.0_wp/13.0_wp
   real(wp), parameter :: LOG_C7 = 1.0_wp/15.0_wp
   real(wp), parameter :: LOG_C8 = 1.0_wp/17.0_wp
   real(wp), parameter :: LOG_C9 = 1.0_wp/19.0_wp

   ! Sin Taylor coefficients (odd powers, alternating sign).  9 terms
   ! gives r¹⁹/19! ≈ 8e-19 at |r| = π/4.
   real(wp), parameter :: SIN_C0 = 1.0_wp                       !  r
   real(wp), parameter :: SIN_C1 = -1.0_wp/6.0_wp                ! -r³/3!
   real(wp), parameter :: SIN_C2 = 1.0_wp/120.0_wp              !  r⁵/5!
   real(wp), parameter :: SIN_C3 = -1.0_wp/5040.0_wp             ! -r⁷/7!
   real(wp), parameter :: SIN_C4 = 1.0_wp/362880.0_wp           !  r⁹/9!
   real(wp), parameter :: SIN_C5 = -1.0_wp/39916800.0_wp         ! -r¹¹/11!
   real(wp), parameter :: SIN_C6 = 1.0_wp/6227020800.0_wp       !  r¹³/13!
   real(wp), parameter :: SIN_C7 = -1.0_wp/1307674368000.0_wp    ! -r¹⁵/15!
   real(wp), parameter :: SIN_C8 = 1.0_wp/355687428096000.0_wp  !  r¹⁷/17!

   ! Cos Taylor coefficients (even powers).  9 terms gives
   ! r¹⁸/18! ≈ 1.4e-19 at |r| = π/4.
   real(wp), parameter :: COS_C0 = 1.0_wp                       !  1
   real(wp), parameter :: COS_C1 = -1.0_wp/2.0_wp                ! -r²/2!
   real(wp), parameter :: COS_C2 = 1.0_wp/24.0_wp               !  r⁴/4!
   real(wp), parameter :: COS_C3 = -1.0_wp/720.0_wp              ! -r⁶/6!
   real(wp), parameter :: COS_C4 = 1.0_wp/40320.0_wp            !  r⁸/8!
   real(wp), parameter :: COS_C5 = -1.0_wp/3628800.0_wp          ! -r¹⁰/10!
   real(wp), parameter :: COS_C6 = 1.0_wp/479001600.0_wp        !  r¹²/12!
   real(wp), parameter :: COS_C7 = -1.0_wp/87178291200.0_wp      ! -r¹⁴/14!
   real(wp), parameter :: COS_C8 = 1.0_wp/20922789888000.0_wp   !  r¹⁶/16!

   ! Overflow / underflow thresholds, derived from the working precision
   ! so the guards stay correct in a single-precision build: exp(x)
   ! overflows once x exceeds log(huge(wp)) (~88.7 in real32, ~709.8 in
   ! real64) and rounds to 0 below log(tiny(wp)).  Hardcoding the real64
   ! values (709.78 / -745.13) let real32 exp() reach +Inf for inputs in
   ! (88.7, 709.8), breaking the documented "overflow -> huge" contract.
   real(wp), parameter :: EXP_OVERFLOW = log(huge(1.0_wp))
   real(wp), parameter :: EXP_UNDERFLOW = log(tiny(1.0_wp))

   ! √(1/2) for log's √2-symmetric range reduction.
   real(wp), parameter :: SQRT_HALF = 0.7071067811865475244_wp

contains

   ! =================================================================
   ! safe_exp
   ! =================================================================

   elemental function safe_exp(x) result(y)
      !! `exp(x)`.  Plain elemental inline around the intrinsic — zero
      !! overhead, codegen identical to writing `exp(x)` directly.  No
      !! production module calls this today (PR-8 removed the
      !! `RAKALI_BITWISE_REPRO` dispatch that once routed it to
      !! `safe_exp_polynomial`, which enabled nothing); see the module
      !! docstring.
      real(wp), intent(in) :: x
      real(wp) :: y
      y = exp(x)
   end function safe_exp

   elemental function safe_exp_polynomial(x) result(y)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Cody-Waite range reduction `x = k*ln(2) + r`, |r| ≤ ln(2)/2,
      !! followed by a 14-term Horner Taylor on `exp(r)` and an
      !! IEEE-exact `2^k` multiply via `scale()`.  Max relative
      !! error ≈ 2-3 ULP across [-709, 709].
      !!
      !! Edge guards: overflow → `huge(y)`, underflow → `0`,
      !! `exp(0) == 1.0` exact.
      real(wp), intent(in) :: x
      real(wp) :: y
      real(wp) :: k_real, r, p
      integer :: k_int

      if (x >= EXP_OVERFLOW) then
         y = huge(y)
         return
      end if
      if (x <= EXP_UNDERFLOW) then
         y = 0.0_wp
         return
      end if

      k_real = anint(x*INV_LN2)
      k_int = nint(k_real)
      r = (x - k_real*LN2_HI) - k_real*LN2_LO

      p = EXP_C13
      p = p*r + EXP_C12
      p = p*r + EXP_C11
      p = p*r + EXP_C10
      p = p*r + EXP_C9
      p = p*r + EXP_C8
      p = p*r + EXP_C7
      p = p*r + EXP_C6
      p = p*r + EXP_C5
      p = p*r + EXP_C4
      p = p*r + EXP_C3
      p = p*r + EXP_C2
      p = p*r + EXP_C1
      p = p*r + EXP_C0

#ifdef LFORTRAN_PASSING
      ! LFortran 0.64 runtime SCALE(x,i) returns 0 for negative i (integer
      ! 2**i); use the exact real power-of-two multiply instead.  Bit-identical
      ! to scale() for powers of two on conforming compilers, which keep scale().
      y = p*2.0_wp**k_int
#else
      y = scale(p, k_int)
#endif
   end function safe_exp_polynomial

   ! =================================================================
   ! safe_log
   ! =================================================================

   elemental function safe_log(x) result(y)
      !! `log(x)` (natural log).  Plain elemental inline around the
      !! intrinsic; see `safe_exp` and the module docstring.
      real(wp), intent(in) :: x
      real(wp) :: y
      y = log(x)
   end function safe_log

   elemental function safe_log_polynomial(x) result(y)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! `log(x)` via IEEE exponent extraction + atanh-form series.
      !!
      !!   1. Write `x = m * 2^k` using `fraction()` / `exponent()`.
      !!      `m` initially lives in [0.5, 1).
      !!   2. √2-symmetric split: if `m < √(1/2)`, halve k and
      !!      double m, putting `m` in [√(1/2), √2) — keeps
      !!      |m - 1| ≤ √2 - 1 ≈ 0.414.
      !!   3. atanh substitution `f = (m-1)/(m+1)` (|f| ≤ 0.172),
      !!      `log(m) = 2*f*(1 + f²/3 + f⁴/5 + ... + f¹⁸/19)`.
      !!   4. `log(x) = k*ln(2) + log(m)`, with `k*ln(2)` evaluated
      !!      via Cody-Waite `LN2_HI / LN2_LO`.
      !!
      !! Accuracy: < 1e-15 relative error across positive doubles.
      !! `log(1) == 0` exact (the polynomial collapses to 0 at f=0
      !! and `k = 0` after the √2 split).
      !!
      !! Domain: `x <= 0` returns `-huge()` (sentinel; callers that
      !! care about IEEE NaN semantics should guard upstream).
      real(wp), intent(in) :: x
      real(wp) :: y
      real(wp) :: m, f, u, p, k_real
      integer :: k_int

      if (x <= 0.0_wp) then
         y = -huge(y)
         return
      end if

      k_int = exponent(x)
      m = fraction(x)

      ! √2-symmetric range: shift [0.5, √(1/2)) up by factor 2.
      if (m < SQRT_HALF) then
         m = 2.0_wp*m
         k_int = k_int - 1
      end if
      ! m now in [√(1/2), √2);  k_int is `exponent(x)` (or one less).

      f = (m - 1.0_wp)/(m + 1.0_wp)
      u = f*f

      ! Horner over u = f²: 1 + u/3 + u²/5 + ... + u⁹/19.
      p = LOG_C9
      p = p*u + LOG_C8
      p = p*u + LOG_C7
      p = p*u + LOG_C6
      p = p*u + LOG_C5
      p = p*u + LOG_C4
      p = p*u + LOG_C3
      p = p*u + LOG_C2
      p = p*u + LOG_C1
      p = p*u + LOG_C0

      k_real = real(k_int, wp)
      y = (k_real*LN2_HI + 2.0_wp*f*p) + k_real*LN2_LO
   end function safe_log_polynomial

   ! =================================================================
   ! safe_sin / safe_cos
   ! =================================================================

   elemental function safe_sin(x) result(y)
      !! `sin(x)`.  Plain elemental inline around the intrinsic; see
      !! `safe_exp` and the module docstring.
      real(wp), intent(in) :: x
      real(wp) :: y
      y = sin(x)
   end function safe_sin

   elemental function safe_cos(x) result(y)
      !! `cos(x)`.  Plain elemental inline around the intrinsic; see
      !! `safe_exp` and the module docstring.
      real(wp), intent(in) :: x
      real(wp) :: y
      y = cos(x)
   end function safe_cos

   elemental function safe_sin_polynomial(x) result(y)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! `sin(x)` via Cody-Waite reduction mod π/2 and a 9-term
      !! Horner Taylor on the reduced argument.  The quadrant
      !! `k mod 4` picks sin vs cos and a sign:
      !!
      !!   k mod 4 = 0:  sin(x) =  sin(r)
      !!   k mod 4 = 1:  sin(x) =  cos(r)
      !!   k mod 4 = 2:  sin(x) = -sin(r)
      !!   k mod 4 = 3:  sin(x) = -cos(r)
      !!
      !! Cody-Waite splitting of π/2 keeps the reduction accurate
      !! for |x| ≲ 2¹⁷.  Beyond that argument-precision starts to
      !! erode; callers reducing huge arguments (e.g. multi-year
      !! tidal phases) should pre-reduce.
      !!
      !! `sin(0) == 0` exact (the polynomial collapses at r=0 to
      !! r·SIN_C0 = 0).
      real(wp), intent(in) :: x
      real(wp) :: y
      real(wp) :: k_real, r
      integer :: k_int, quadrant

      k_real = anint(x*INV_PI_OVER_2)
      k_int = nint(k_real)
      r = (x - k_real*PI_OVER_2_HI) - k_real*PI_OVER_2_LO
      quadrant = modulo(k_int, 4)

      select case (quadrant)
      case (0)
         y = sin_reduced(r)
      case (1)
         y = cos_reduced(r)
      case (2)
         y = -sin_reduced(r)
      case default  ! 3
         y = -cos_reduced(r)
      end select
   end function safe_sin_polynomial

   elemental function safe_cos_polynomial(x) result(y)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! `cos(x)` via the same Cody-Waite reduction; the quadrant
      !! map is shifted by 1 vs sin:
      !!
      !!   k mod 4 = 0:  cos(x) =  cos(r)
      !!   k mod 4 = 1:  cos(x) = -sin(r)
      !!   k mod 4 = 2:  cos(x) = -cos(r)
      !!   k mod 4 = 3:  cos(x) =  sin(r)
      !!
      !! `cos(0) == 1` exact.
      real(wp), intent(in) :: x
      real(wp) :: y
      real(wp) :: k_real, r
      integer :: k_int, quadrant

      k_real = anint(x*INV_PI_OVER_2)
      k_int = nint(k_real)
      r = (x - k_real*PI_OVER_2_HI) - k_real*PI_OVER_2_LO
      quadrant = modulo(k_int, 4)

      select case (quadrant)
      case (0)
         y = cos_reduced(r)
      case (1)
         y = -sin_reduced(r)
      case (2)
         y = -cos_reduced(r)
      case default  ! 3
         y = sin_reduced(r)
      end select
   end function safe_cos_polynomial

   elemental function sin_reduced(r) result(y)
      !! Horner sin Taylor on |r| ≤ π/4.  `sin(r) = r * P(r²)` where
      !! P is the odd-power polynomial in r² — avoids computing r¹,
      !! r³, r⁵... separately.  9 terms (highest power r¹⁷) gives
      !! < 1e-18 truncation error.
      real(wp), intent(in) :: r
      real(wp) :: y
      real(wp) :: u, p

      u = r*r
      p = SIN_C8
      p = p*u + SIN_C7
      p = p*u + SIN_C6
      p = p*u + SIN_C5
      p = p*u + SIN_C4
      p = p*u + SIN_C3
      p = p*u + SIN_C2
      p = p*u + SIN_C1
      p = p*u + SIN_C0
      y = r*p
   end function sin_reduced

   elemental function cos_reduced(r) result(y)
      !! Horner cos Taylor on |r| ≤ π/4.  9 even-power terms.
      real(wp), intent(in) :: r
      real(wp) :: y
      real(wp) :: u, p

      u = r*r
      p = COS_C8
      p = p*u + COS_C7
      p = p*u + COS_C6
      p = p*u + COS_C5
      p = p*u + COS_C4
      p = p*u + COS_C3
      p = p*u + COS_C2
      p = p*u + COS_C1
      p = p*u + COS_C0
      y = p
   end function cos_reduced

   ! =================================================================
   ! safe_sqrt / safe_pow — intrinsic for now
   ! =================================================================

   elemental function safe_sqrt(x) result(y)
      !! IEEE-correctly-rounded by mandate, so the intrinsic is
      !! already bit-identical across compliant compilers.
      real(wp), intent(in) :: x
      real(wp) :: y
      y = sqrt(x)
   end function safe_sqrt

   elemental function safe_pow(base, exponent) result(y)
      !! Non-integer exponent only.  For integer exponents
      !! (`a**2`, `a**3`, ...), write the multiplication out — the
      !! compiler unrolls those and they're already deterministic.
      !! Plain elemental inline around the intrinsic; see `safe_exp`
      !! and the module docstring.
      real(wp), intent(in) :: base
      real(wp), intent(in) :: exponent
      real(wp) :: y
      y = base**exponent
   end function safe_pow

   elemental function safe_pow_polynomial(base, exponent) result(y)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! `base^exponent` via the identity
      !!   base^exponent = exp(exponent * log(base))
      !! Routes through the polynomial log and exp paths so the
      !! result is deterministic across builds.
      !!
      !! Domain:
      !!   * `base > 0` — proper case; returns the polynomial value.
      !!   * `base == 0` — returns 0 if `exponent > 0`, `huge()` if
      !!     `exponent < 0`, `1` if `exponent == 0` (the conventional
      !!     IEEE-754 `pow(0, 0) = 1` rule).
      !!   * `base < 0` — non-integer exponent is undefined for
      !!     reals; returns `-huge()` as a sentinel.  Callers that
      !!     need integer `base**n` for negative `base` should use
      !!     repeated multiplication, not `safe_pow`.
      !!
      !! Accuracy: composes two ~ULP-level polynomial passes plus
      !! one multiplication, so expect ~5-10 ULP relative error.
      real(wp), intent(in) :: base, exponent
      real(wp) :: y
      real(wp) :: log_base

      if (base > 0.0_wp) then
         log_base = safe_log_polynomial(base)
         y = safe_exp_polynomial(exponent*log_base)
      else if (base == 0.0_wp) then
         if (exponent > 0.0_wp) then
            y = 0.0_wp
         else if (exponent < 0.0_wp) then
            y = huge(y)
         else
            y = 1.0_wp
         end if
      else
         y = -huge(y)
      end if
   end function safe_pow_polynomial

end module rki_safe_math
