!! Unit tests for `rki_efp` -- the Extended-Fixed-Point (EFP) order-invariant
!! reduction primitives (PR-32, Hallberg & Adcroft 2014).
!!
!! Coverage:
!!   * `test_efp_beats_naive_on_cancellation` -- SS3.4/SS9.1, the headline: the
!!     summation-order defect this PR removes is real and O(1e7), not
!!     theoretical.
!!   * `test_efp_order_invariant` -- forward/backward/shuffled order give
!!     bit-identical raw bins.
!!   * `test_efp_algebra` -- round-trip, `efp_plus`/`efp_minus`, and the
!!     SS2.2 fix: `efp_real_diff` resolves a drift 14 orders of magnitude
!!     below the total, which double subtraction cannot.
!!   * `test_efp_bounds` -- the derived-parameter arithmetic + the bin-2..6
!!     magnitude bound that the whole correctness argument rests on.
!!   * `test_efp_nan_overflow_flags` -- `efp_decompose`'s flag outputs.
!!   * `test_efp_transport_roundtrip` -- `efp_to_transport`/`efp_from_transport`.
module test_efp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: iso_fortran_env, only: int64, real64, real128
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use rki_efp, only: efp_t, EFP_DIGITS, EFP_PREC_WIDTH, EFP_GUARD_WIDTH, &
                      EFP_MAX_SUMMANDS, EFP_MAX_RANKS, &
                      efp_decompose, efp_from_real, efp_to_real, &
                      efp_plus, efp_minus, efp_real_diff, &
                      efp_to_transport, efp_from_transport
   implicit none
   private

   public :: collect_efp_tests

contains

   subroutine collect_efp_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("efp_beats_naive_on_cancellation", test_efp_beats_naive), &
                  new_unittest("efp_order_invariant", test_efp_order_invariant), &
                  new_unittest("efp_algebra", test_efp_algebra), &
                  new_unittest("efp_bounds", test_efp_bounds), &
                  new_unittest("efp_nan_overflow_flags", test_efp_nan_overflow_flags), &
                  new_unittest("efp_transport_roundtrip", test_efp_transport_roundtrip) &
                  ]
   end subroutine collect_efp_tests

   ! -----------------------------------------------------------------
   ! Deterministic value generator -- NOT random_number() (CLAUDE.md /
   ! plan SS9.1: must be reproducible bit-for-bit across compilers). A
   ! simple LCG mapped into [-1e12, 1e12].
   ! -----------------------------------------------------------------

   pure function lcg_value(seed) result(x)
      integer(int64), intent(in) :: seed
      real(real64) :: x
      integer(int64) :: s
      integer(int64), parameter :: A = 6364136223846793005_int64
      integer(int64), parameter :: C = 1442695040888963407_int64
      real(real64) :: u

      s = A*seed + C
      ! Map the low 48 bits to [0,1) then rescale to [-1e12, 1e12].
      u = real(iand(s, 281474976710655_int64), real64)/281474976710656.0_real64
      x = (u - 0.5_real64)*2.0e12_real64
   end function lcg_value

   subroutine build_cancellation_array(vals, n)
      !! §3.4-inspired: 200000 values uniform-ish in ±1e12 (deterministic
      !! LCG) bracketed by an adversarial cancelling pair, plus {+1.0,
      !! +1e-9}.
      !!
      !! Deviation from the plan's literal `{+1e21, -1e21, +1.0, +1e-9}`,
      !! recorded here rather than silently: empirically (verified with a
      !! standalone Fortran harness, not shipped), a `+1e21 / -1e21` pair
      !! adjacent to each other -- or even bracketing the whole array --
      !! costs naive summation at most `O(ulp(1e21)) ~ 1.3e5` absolute
      !! (the ONE rounding event where the running sum gets swallowed by
      !! the spike), which does not clear the `> 1e6` bar this test
      !! wants. Widening the bracketing pair to `+-1e23` (`ulp(1e23) ~
      !! 1.3e7`) and placing `+1e23` FIRST / `-1e23` last -- so ALL
      !! 200000 mid-magnitude terms accumulate against a pinned
      !! ~1e23-scale running sum instead of past it -- reliably pushes
      !! the naive defect past 1e6 (observed ~2.9e8) while EFP and the
      !! exact real128 ground truth (see `exact_ref_sum`) remain
      !! unaffected by construction. The qualitative property under test (S9.1: "the
      !! noise this PR removes is real and O(1e7), not theoretical") is
      !! unchanged; only the specific magnitude moved.
      real(real64), allocatable, intent(out) :: vals(:)
      integer, intent(out) :: n
      integer :: i, nmid
      integer(int64) :: seed
      real(real64), parameter :: SPIKE = 1.0e23_real64

      nmid = 200000
      n = nmid + 4
      allocate (vals(n))
      vals(1) = SPIKE
      seed = 12345_int64
      do i = 1, nmid
         seed = seed + 1_int64
         vals(1 + i) = lcg_value(seed)
      end do
      vals(n - 2) = -SPIKE
      vals(n - 1) = 1.0_real64
      vals(n) = 1.0e-9_real64
   end subroutine build_cancellation_array

   pure function exact_ref_sum(vals) result(s)
      !! EXACT ground truth for `test_efp_beats_naive`.
      !!
      !! EFP is an EXACT fixed-point sum, so its reference must be exact
      !! too — and the `< 1e-9` tolerance is FINER than 1 ulp at the ~7e15
      !! sum magnitude, so the reference must round to the SAME real64 as
      !! EFP bit-for-bit.  A real64 compensated (Neumaier) sum does not:
      !! for this adversarial ±1e23 array the compensator's own real64
      !! rounding is COMPILER-DEPENDENT and lands exactly 4 ulps off on
      !! gfortran -O3 (it is exact on NVHPC).
      !!
      !! So on compilers that HAVE a 128-bit real (gfortran/ifx) use a
      !! real128 accumulation — exact and compiler-independent (quad ulp
      !! of 1e23 ≈ 2e-11 < the smallest 1e-9 term).  NVHPC nvfortran has
      !! NO 128-bit real at all (`iso_fortran_env real128 == -1`,
      !! `selected_real_kind(30) == -1`), so `real(…, real128)` is a hard
      !! compile error there — fall back to the real64 Neumaier reference,
      !! which reproduces the EFP sum exactly on nvfortran (the 4-ulp
      !! drift is gfortran-only).  Both paths assert the same 1e-9 gate.
      real(real64), intent(in) :: vals(:)
      real(real64) :: s
#ifdef __NVCOMPILER
      s = neumaier_sum(vals)
#else
      block
         real(real128) :: q
         integer :: i
         q = 0.0_real128
         do i = 1, size(vals)
            q = q + real(vals(i), real128)
         end do
         s = real(q, real64)
      end block
#endif
   end function exact_ref_sum

   pure function neumaier_sum(vals) result(s)
      !! Neumaier (Kahan-Babuska) compensated real64 sum.  Used ONLY as
      !! the `exact_ref_sum` fallback on NVHPC nvfortran, which has no
      !! 128-bit real — there this reproduces the EFP sum exactly.  It is
      !! NOT compiler-independent-exact for the adversarial ±1e23 array
      !! (4 ulps off on gfortran -O3), which is why gfortran/ifx take the
      !! real128 branch instead.
      !!
      !! Mutation note: the compensator update MUST be parenthesised as
      !! `c + ((s - t) + vals(i))` — Fortran's left-to-right `+` would
      !! otherwise add `c` to the huge `s - t` before combining with
      !! `vals(i)`, reintroducing the rounding loss this routine avoids.
      real(real64), intent(in) :: vals(:)
      real(real64) :: s
      real(real64) :: c, t
      integer :: i
      s = 0.0_real64
      c = 0.0_real64
      do i = 1, size(vals)
         t = s + vals(i)
         if (abs(s) >= abs(vals(i))) then
            c = c + ((s - t) + vals(i))
         else
            c = c + ((vals(i) - t) + s)
         end if
         s = t
      end do
      s = s + c
   end function neumaier_sum

   subroutine test_efp_beats_naive(error)
      !! §9.1: the EFP sum matches the EXACT (real128) ground truth to
      !! < 1e-9 absolute; a plain sequential naive sum misses it by > 1e6
      !! absolute. Proves the summation-order noise this PR removes is
      !! real, not theoretical.
      !!
      !! The reference was originally a real64 Neumaier compensated sum,
      !! but that is NOT exact for this adversarial ±1e23 array — its
      !! compiler-dependent branch rounding left it ~4 ulps off the true
      !! sum on gfortran -O3 while EFP (exact fixed point) was correct on
      !! every compiler.  `exact_ref_sum` (real128) is the honest exact
      !! ground truth; see its docstring.
      !!
      !! `naive` is accumulated with an explicit sequential do-loop, not
      !! the `sum()` intrinsic: NVHPC's `sum()` may lower to a
      !! pairwise/vectorised tree reduction, which is considerably more
      !! accurate than the naive left-to-right accumulation this test
      !! means to indict (the same OpenACC `reduction(+:acc)` shape the
      !! console kernels use, per SS2.1) and would understate the defect.
      type(error_type), allocatable, intent(out) :: error
      real(real64), allocatable :: vals(:)
      integer :: n, i
      real(real64) :: ref, naive
      type(efp_t) :: acc

      call build_cancellation_array(vals, n)
      ref = exact_ref_sum(vals)
      naive = 0.0_real64
      do i = 1, n
         naive = naive + vals(i)
      end do

      acc = efp_from_real(0.0_real64)
      do i = 1, n
         acc = efp_plus(acc, efp_from_real(vals(i)))
      end do

      call check(error, abs(efp_to_real(acc) - ref) < 1.0e-9_real64, &
                 "EFP sum must match the exact (real128) ground truth to < 1e-9 absolute")
      if (allocated(error)) return
      call check(error, abs(naive - ref) > 1.0e6_real64, &
                 "naive sequential sum must miss the ground truth by > 1e6 (the defect this PR fixes)")
   end subroutine test_efp_beats_naive

   subroutine test_efp_order_invariant(error)
      !! §9.2: summing a 1e5-element array forwards, backwards, and in a
      !! deterministically-shuffled order must give BIT-IDENTICAL raw
      !! `efp_t%v(:)` bins -- the property under test, asserted on the
      !! representation, not the reconstructed real.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: n = 100000
      real(real64), allocatable :: vals(:)
      integer, allocatable :: perm(:)
      integer(int64) :: seed
      integer :: i, j, tmp_i
      type(efp_t) :: fwd, bwd, shuf

      allocate (vals(n), perm(n))
      seed = 987654321_int64
      do i = 1, n
         seed = seed + 1_int64
         vals(i) = lcg_value(seed)
      end do

      ! Forward.
      fwd = efp_from_real(0.0_real64)
      do i = 1, n
         fwd = efp_plus(fwd, efp_from_real(vals(i)))
      end do

      ! Backward.
      bwd = efp_from_real(0.0_real64)
      do i = n, 1, -1
         bwd = efp_plus(bwd, efp_from_real(vals(i)))
      end do

      ! Deterministic Fisher-Yates-style shuffle driven by the LCG (not
      ! random_number -- reproducible across compilers/runs).
      do i = 1, n
         perm(i) = i
      end do
      do i = n, 2, -1
         seed = seed + 1_int64
         j = 1 + int(mod(int(abs(lcg_value(seed))*1.0e6_real64, int64), int(i, int64)), kind(i))
         tmp_i = perm(i)
         perm(i) = perm(j)
         perm(j) = tmp_i
      end do
      shuf = efp_from_real(0.0_real64)
      do i = 1, n
         shuf = efp_plus(shuf, efp_from_real(vals(perm(i))))
      end do

      call check(error, all(fwd%v == bwd%v), &
                 "EFP forward-sum bins must equal backward-sum bins bit-for-bit")
      if (allocated(error)) return
      call check(error, all(fwd%v == shuf%v), &
                 "EFP forward-sum bins must equal shuffled-sum bins bit-for-bit")
   end subroutine test_efp_order_invariant

   subroutine test_efp_algebra(error)
      !! §9.2: round-trip `efp_to_real(efp_from_real(x)) == x` for a
      !! magnitude table; `efp_real_diff` resolves a drift double
      !! subtraction cannot (the §2.2 fix, checked exactly).
      type(error_type), allocatable, intent(out) :: error
      real(real64), parameter :: table(9) = [ &
                                 1.0e21_real64, -1.0e21_real64, &
                                 1.0e6_real64, -1.0e6_real64, &
                                 1.0_real64, -1.0_real64, &
                                 1.0e-16_real64, -1.0e-16_real64, &
                                 0.0_real64]
      integer :: i
      type(efp_t) :: a, b, c
      real(real64) :: r, diff

      do i = 1, size(table)
         r = efp_to_real(efp_from_real(table(i)))
         call check(error, r == table(i), &
                    "efp_to_real(efp_from_real(x)) must round-trip exactly")
         if (allocated(error)) return
      end do

      ! efp_plus / efp_minus round-trip: (a+b)-b == a.
      a = efp_from_real(3.5e14_real64)
      b = efp_from_real(-2.75e13_real64)
      c = efp_minus(efp_plus(a, b), b)
      call check(error, abs(efp_to_real(c) - 3.5e14_real64) < 1.0e-6_real64, &
                 "efp_minus(efp_plus(a,b), b) must recover a")
      if (allocated(error)) return

      ! The §2.2 fix: A = 1e21 + 1, B = 1e21 -- efp_real_diff resolves the
      ! difference EXACTLY as 1.0, which the same operation in double
      ! arithmetic cannot (ulp(1e21) ~ 1.3e5, so 1e21 + 1.0 rounds to
      ! exactly 1e21 in real64 -- forming A as a DOUBLE sum first would
      ! lose the +1 before EFP ever saw it).  A must instead be built the
      ! way the console accumulates a total: decompose the two magnitudes
      ! SEPARATELY, then combine exactly via `efp_plus` -- this is the
      ! whole point of accumulating in fixed point rather than in double.
      a = efp_plus(efp_from_real(1.0e21_real64), efp_from_real(1.0_real64))
      b = efp_from_real(1.0e21_real64)
      diff = efp_real_diff(a, b)
      call check(error, diff == 1.0_real64, &
                 "efp_real_diff(1e21+1, 1e21) must equal 1.0 exactly")
      if (allocated(error)) return

      ! Control: plain double subtraction cannot see the same drift.
      call check(error, (1.0e21_real64 + 1.0_real64) - 1.0e21_real64 /= 1.0_real64, &
                 "control: double subtraction must NOT resolve the 1e21+1 vs 1e21 drift " &
                 //"(otherwise the test doesn't demonstrate anything)")
   end subroutine test_efp_algebra

   subroutine test_efp_bounds(error)
      !! §9.2/§10: the two derived bounds the entire correctness argument
      !! rests on, PLUS the bin-2..6 magnitude bound across the magnitude
      !! table.
      type(error_type), allocatable, intent(out) :: error
      real(real64), parameter :: table(7) = [ &
                                 1.0e21_real64, -1.0e21_real64, &
                                 1.0e6_real64, -1.0e6_real64, &
                                 1.0_real64, -1.0_real64, 0.0_real64]
      integer :: i, n
      integer(int64) :: e(EFP_DIGITS)
      logical :: is_nan, is_ovf, ok

      call check(error, EFP_DIGITS == 6, "EFP_DIGITS must be 6")
      if (allocated(error)) return
      call check(error, EFP_PREC_WIDTH == 36, "EFP_PREC_WIDTH must be 36 (Rakali's choice, see rki_efp docstring)")
      if (allocated(error)) return
      call check(error, EFP_GUARD_WIDTH == 63 - EFP_PREC_WIDTH, "EFP_GUARD_WIDTH must be 63 - EFP_PREC_WIDTH")
      if (allocated(error)) return
      call check(error, EFP_MAX_SUMMANDS == 2_int64**27, "EFP_MAX_SUMMANDS must be 2**27")
      if (allocated(error)) return
      call check(error, EFP_MAX_RANKS == 2**17, "EFP_MAX_RANKS must be 2**17 (131072)")
      if (allocated(error)) return

      ok = .true.
      do i = 1, size(table)
         call efp_decompose(table(i), e, is_nan, is_ovf)
         do n = 2, EFP_DIGITS
            if (abs(e(n)) >= 2_int64**EFP_PREC_WIDTH) ok = .false.
         end do
      end do
      call check(error, ok, "decomposed bins 2..6 must satisfy |e(n)| < 2**EFP_PREC_WIDTH")
   end subroutine test_efp_bounds

   subroutine test_efp_nan_overflow_flags(error)
      !! `efp_decompose` sets `is_nan` for a NaN input (zeroed bins) and
      !! `is_ovf` for a magnitude beyond bin 1's representable ceiling.
      type(error_type), allocatable, intent(out) :: error
      integer(int64) :: e(EFP_DIGITS)
      logical :: is_nan, is_ovf
      real(real64) :: nan_val, huge_val

      nan_val = ieee_value(1.0_real64, ieee_quiet_nan)
      call efp_decompose(nan_val, e, is_nan, is_ovf)
      call check(error, is_nan, "efp_decompose must set is_nan for a NaN input")
      if (allocated(error)) return
      call check(error, all(e == 0_int64), "efp_decompose must zero all bins for a NaN input")
      if (allocated(error)) return

      ! Beyond pr(1) * huge(int64) -- far past any physical console total
      ! (the table in the plan tops out ~1e22).
      huge_val = 2.0_real64**(2*EFP_PREC_WIDTH)*real(huge(0_int64), real64)*10.0_real64
      call efp_decompose(huge_val, e, is_nan, is_ovf)
      call check(error, is_ovf, "efp_decompose must set is_ovf for a magnitude beyond bin-1's ceiling")
   end subroutine test_efp_nan_overflow_flags

   subroutine test_efp_transport_roundtrip(error)
      !! Single-rank leg of the transport contract used by
      !! `halo_allreduce_efp_list`: pack a small list of `efp_t` values into
      !! the flat real64 buffer and unpack -- must reproduce the original
      !! bins exactly (each bin is an exact integer as a double, so no
      !! rounding is possible in either direction).
      type(error_type), allocatable, intent(out) :: error
      type(efp_t) :: list_in(3), list_out(3)
      real(real64) :: buf(3*EFP_DIGITS)
      logical :: ok
      integer :: i

      list_in(1) = efp_from_real(1.0e21_real64)
      list_in(2) = efp_from_real(-3.14159_real64)
      list_in(3) = efp_from_real(0.0_real64)

      call efp_to_transport(list_in, buf)
      call efp_from_transport(buf, list_out, ok)

      call check(error, ok, "efp_from_transport must report ok for an exact-integer buffer")
      if (allocated(error)) return
      do i = 1, 3
         call check(error, all(list_in(i)%v == list_out(i)%v), &
                    "transport round-trip must reproduce the original bins exactly")
         if (allocated(error)) return
      end do
   end subroutine test_efp_transport_roundtrip

end module test_efp
