!! Extended-Fixed-Point (EFP) order-invariant reduction primitives.
module rki_efp
   !! Order-invariant real-to-integer summation (Hallberg & Adcroft, 2014,
   !! *An Order-invariant Real-to-Integer Conversion Sum*, Parallel Computing
   !! 40(5-6), doi:10.1016/j.parco.2014.04.007).
   !!
   !! Floating-point addition is not associative, so `Sum x_i` depends on the
   !! order of summation -- the reduction tree shape, the rank count, the
   !! decomposition.  This module maps each real onto a fixed-point integer
   !! vector (`efp_t`) so accumulation becomes plain `integer(int64)`
   !! addition, which IS associative and exact.  A base exponent `P`
   !! (`EFP_PREC_WIDTH`) splits the mantissa into `EFP_DIGITS = 6` bins of
   !! weight `2**(P*(3-n))` for `n = 1..6`; the decomposition
   !! (`efp_decompose`) is a pure function of the input value, so two ranks
   !! decomposing the same value produce identical bins, and summing bins
   !! (`efp_plus`) is then order-invariant by construction.
   !!
   !! This module is host-only, pure-arithmetic, and depends on nothing but
   !! `rki_constants` (for `wp`) and `iso_fortran_env` -- no MPI, no grid, no
   !! state -- so `src/comm/` can consume it (`halo_allreduce_efp_list`)
   !! with no dependency cycle, and `rki_console_stats` (shared coastal +
   !! ocean) can too.
   !!
   !! **Parameter choice**: `EFP_PREC_WIDTH = 36` (MOM6 uses 46).  Rakali
   !! picks a smaller `P` because the cross-rank transport (see
   !! `halo_allreduce_efp_list` in `src/comm/`) has no `integer(int64)`
   !! MPI allreduce available (`pic_mpi_lib` binds only `dp`/`sp`/`i32`
   !! allreduce overloads) and instead transports the bins as exactly-
   !! represented `real64` values summed by `MPI_SUM` -- exact only while
   !! every partial sum stays within the 53-bit double mantissa.  At
   !! `P = 36` that bounds the rank count to `EFP_MAX_RANKS = 2**17 =
   !! 131072`, far beyond any Rakali run, while `EFP_MAX_SUMMANDS = 2**27
   !! ~= 1.34e8` per local reduction block comfortably covers a single
   !! k-slab of an 11500^2 grid.  See `docs/CAPABILITIES_AND_LIMITATIONS.md`
   !! for the full bound table.
   !!
   !! **Do NOT** widen `EFP_PREC_WIDTH` without re-deriving `EFP_MAX_RANKS`
   !! and `EFP_MAX_SUMMANDS` -- both are `parameter`s computed FROM it, and
   !! `test_efp_bounds` pins the arithmetic.
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   private

   public :: efp_t
   public :: EFP_DIGITS, EFP_PREC_WIDTH, EFP_GUARD_WIDTH
   public :: EFP_MAX_SUMMANDS, EFP_MAX_RANKS
   public :: efp_decompose, efp_carry, efp_regularize
   public :: efp_from_real, efp_to_real
   public :: efp_plus, efp_minus, efp_real_diff
   public :: efp_to_transport, efp_from_transport
   public :: efp_bin1_within_transport_bound

   integer, parameter :: EFP_DIGITS = 6
      !! Number of fixed-point bins per `efp_t` value.
   integer, parameter :: EFP_PREC_WIDTH = 36
      !! Bits per bin (`P`).  See the module docstring for the derivation.
   integer, parameter :: EFP_GUARD_WIDTH = 63 - EFP_PREC_WIDTH
      !! = 27.  `int64` has 63 usable bits (sign-magnitude use here, not
      !! two's-complement range); this is the headroom for accumulating
      !! multiple summands into one bin before it could overflow.
   integer(int64), parameter :: EFP_MAX_SUMMANDS = 2_int64**EFP_GUARD_WIDTH
      !! = 134217728 (~1.34e8).  Local (single-rank, single k-slab) upper
      !! bound on the number of values that may be accumulated into one
      !! bin via `efp_plus` before an `efp_carry` is required to keep
      !! `|bin| < 2**63`.  Enforced fail-loud at the kernel call site
      !! (`rki_ocean_console_stats`), not silently.
   integer, parameter :: EFP_MAX_RANKS = 2**(53 - EFP_PREC_WIDTH)
      !! = 131072.  Upper bound on the number of MPI ranks the double-
      !! precision transport (`halo_allreduce_efp_list`) can combine
      !! exactly: after a local `efp_carry`, bins 2..6 satisfy
      !! `|e(n)| < 2**P`, so a partial sum over `EFP_MAX_RANKS` ranks stays
      !! `<= 2**53`, the largest exactly-representable double integer.

   ! -- fixed-point bin weights pr(n) = 2**(P*(3-n)), n = 1..6, and their
   ! exact reciprocals I_pr(n) = 2**(-P*(3-n)).  Six NAMED SCALARS (not an
   ! array) so the decomposition unrolls with no array indexing -- see the
   ! in-module `!$acc routine seq` duplicate `efp_decompose_impl` in
   ! `rki_ocean_console_stats`, which needs the same scalar shape to inline
   ! into a device reduction loop (OpenACC has no portable array
   ! reduction, and MOM6's `!$omp declare target(pr, I_pr)` device-constant-
   ! array route is not available here).
   real(real64), parameter :: EFP_PR1 = 2.0_real64**(2*EFP_PREC_WIDTH)
   real(real64), parameter :: EFP_PR2 = 2.0_real64**(1*EFP_PREC_WIDTH)
   real(real64), parameter :: EFP_PR3 = 1.0_real64
   real(real64), parameter :: EFP_PR4 = 2.0_real64**(-1*EFP_PREC_WIDTH)
   real(real64), parameter :: EFP_PR5 = 2.0_real64**(-2*EFP_PREC_WIDTH)
   real(real64), parameter :: EFP_PR6 = 2.0_real64**(-3*EFP_PREC_WIDTH)
   real(real64), parameter :: EFP_IPR1 = 1.0_real64/EFP_PR1
   real(real64), parameter :: EFP_IPR2 = 1.0_real64/EFP_PR2
   real(real64), parameter :: EFP_IPR3 = 1.0_real64/EFP_PR3
   real(real64), parameter :: EFP_IPR4 = 1.0_real64/EFP_PR4
   real(real64), parameter :: EFP_IPR5 = 1.0_real64/EFP_PR5
   real(real64), parameter :: EFP_IPR6 = 1.0_real64/EFP_PR6

   integer(int64), parameter :: EFP_PREC_I64 = 2_int64**EFP_PREC_WIDTH
      !! `2**P` as an int64 -- the bin-2..6 magnitude bound used by
      !! `efp_carry` / `efp_regularize`.

   type :: efp_t
      !! An order-invariant fixed-point value: `Sum_n v(n) * pr(n)`.  The
      !! component is PUBLIC (unlike MOM6's private `EFP_type%v`) so the
      !! comm facade (`halo_allreduce_efp_list`) can pack/unpack it without
      !! an accessor procedure -- Rakali splits the arithmetic module
      !! (here) from the collective (in `src/comm/`), so the layering
      !! requires `v` to be reachable from both.  Tests asserting bit-
      !! identity compare `v(:)` directly, never the reconstructed real.
      integer(int64) :: v(EFP_DIGITS) = 0_int64
   end type efp_t

contains

   pure subroutine efp_decompose(r, e, is_nan, is_ovf)
      !! Greedy sign-magnitude fixed-point decomposition of `r` into six
      !! `int64` bins of weight `pr(n)`, `n = 1..6`.  Unrolled (no loop, no
      !! array indexing over `pr`/`I_pr`) so this is a template a
      !! `!$acc routine seq` in-module duplicate can mirror exactly (see
      !! `rki_ocean_console_stats::efp_decompose_impl`, which pins against
      !! this procedure in `test_efp_impl_matches_canonical`).
      !!
      !! `is_nan` is set (and `e = 0`) for a NaN input, using
      !! `ieee_is_nan`-equivalent structural comparison so the module stays
      !! independent of `ieee_arithmetic`'s import list elsewhere.
      !! `is_ovf` is set when `|r| >= pr(1) * huge(1_int64)` -- the largest
      !! magnitude bin 1 can represent -- and `e` is truncated to that
      !! bound rather than silently wrapping.
      real(real64), intent(in) :: r
      integer(int64), intent(out) :: e(EFP_DIGITS)
      logical, intent(out) :: is_nan
      logical, intent(out) :: is_ovf
      real(real64) :: rs, s
      real(real64), parameter :: MAX_E1 = real(huge(0_int64), real64)
         !! `huge(1_int64)` widened to real64 (rounds to the nearest
         !! representable double, ~9.223372036854776e18) -- the largest
         !! magnitude bin 1 can safely hold; used as the overflow ceiling
         !! below.

      is_nan = ieee_is_nan(r)
      is_ovf = .false.
      e = 0_int64
      if (is_nan) return

      s = 1.0_real64
      rs = r
      if (rs < 0.0_real64) then
         s = -1.0_real64
         rs = -rs
      end if

      if (rs*EFP_IPR1 >= MAX_E1) then
         is_ovf = .true.
         e(1) = int(sign(MAX_E1, s), int64)
         return
      end if

      e(1) = int(s*aint(rs*EFP_IPR1), int64)
      rs = rs - real(abs(e(1)), real64)*EFP_PR1

      e(2) = int(s*aint(rs*EFP_IPR2), int64)
      rs = rs - real(abs(e(2)), real64)*EFP_PR2

      e(3) = int(s*aint(rs*EFP_IPR3), int64)
      rs = rs - real(abs(e(3)), real64)*EFP_PR3

      e(4) = int(s*aint(rs*EFP_IPR4), int64)
      rs = rs - real(abs(e(4)), real64)*EFP_PR4

      e(5) = int(s*aint(rs*EFP_IPR5), int64)
      rs = rs - real(abs(e(5)), real64)*EFP_PR5

      e(6) = int(s*aint(rs*EFP_IPR6), int64)
      ! Remainder `rs - |e(6)|*pr(6)` is the deterministic dropped quantum,
      ! `0 <= remainder < pr(6) = 2**-3P` -- not retained (matches MOM6:
      ! only six bins are kept).
   end subroutine efp_decompose

   pure subroutine efp_carry(e)
      !! Renormalise bins 6..2 into `(-2**P, 2**P)`, propagating the excess
      !! into the next-more-significant bin, without changing the
      !! represented value.  Bin 1 is left untouched (unbounded by
      !! construction; see the module docstring on `EFP_MAX_RANKS`).
      !! Mirrors MOM6 `carry_overflow`, which loops
      !! `EFP_DIGITS..2` for the same reason.
      integer(int64), intent(inout) :: e(EFP_DIGITS)
      integer :: n
      integer(int64) :: carry_amt

      do n = EFP_DIGITS, 2, -1
         if (e(n) >= EFP_PREC_I64 .or. e(n) <= -EFP_PREC_I64) then
            carry_amt = e(n)/EFP_PREC_I64
            e(n) = e(n) - carry_amt*EFP_PREC_I64
            e(n - 1) = e(n - 1) + carry_amt
         end if
      end do
   end subroutine efp_carry

   pure subroutine efp_regularize(e)
      !! `efp_carry` plus: force every bin to share the overall sign, so a
      !! single well-conditioned FP accumulation (`efp_to_real`) can form
      !! `Sum pr(n)*e(n)` without alternating-sign cancellation error.
      !! Mirrors MOM6 `regularize_ints`.
      !!
      !! After `efp_carry`, bins 2..6 satisfy `|e(n)| < 2**P`.  The overall
      !! sign is that of the first (most-significant) nonzero bin.  A
      !! single descending pass then borrows/carries any bin of the
      !! opposite sign into `[0, 2**P)` (or `(-2**P, 0]`) from its
      !! next-more-significant neighbour -- exact, since
      !! `e(n-1)*2**P + e(n)` is invariant under `(e(n-1) -+ 1, e(n) +- 2**P)`
      !! -- and cascades: a neighbour driven negative/positive by one
      !! borrow is itself fixed on the next loop iteration.
      integer(int64), intent(inout) :: e(EFP_DIGITS)
      integer :: n
      logical :: positive

      call efp_carry(e)

      positive = .true.
      do n = 1, EFP_DIGITS
         if (e(n) /= 0_int64) then
            positive = e(n) > 0_int64
            exit
         end if
      end do

      if (positive) then
         do n = EFP_DIGITS, 2, -1
            if (e(n) < 0_int64) then
               e(n) = e(n) + EFP_PREC_I64
               e(n - 1) = e(n - 1) - 1_int64
            end if
         end do
      else
         do n = EFP_DIGITS, 2, -1
            if (e(n) > 0_int64) then
               e(n) = e(n) - EFP_PREC_I64
               e(n - 1) = e(n - 1) + 1_int64
            end if
         end do
      end if
   end subroutine efp_regularize

   pure function efp_from_real(r) result(a)
      !! `real64 -> efp_t`.  Wraps `efp_decompose`, discarding the NaN /
      !! overflow flags (callers needing them should call `efp_decompose`
      !! directly).
      real(real64), intent(in) :: r
      type(efp_t) :: a
      logical :: is_nan, is_ovf
      call efp_decompose(r, a%v, is_nan, is_ovf)
   end function efp_from_real

   pure function efp_to_real(a) result(r)
      !! `efp_t -> real64`.  Regularises a LOCAL COPY of `a` (never mutates
      !! the argument) -- `pure` with `intent(in)`, per
      !! `FORTRAN_STYLE.md`'s "default new procedures to pure" (MOM6's
      !! `EFP_to_real` instead mutates its `intent(inout)` argument).
      type(efp_t), intent(in) :: a
      real(real64) :: r
      integer(int64) :: e(EFP_DIGITS)
      e = a%v
      call efp_regularize(e)
      r = real(e(1), real64)*EFP_PR1 + real(e(2), real64)*EFP_PR2 &
          + real(e(3), real64)*EFP_PR3 + real(e(4), real64)*EFP_PR4 &
          + real(e(5), real64)*EFP_PR5 + real(e(6), real64)*EFP_PR6
   end function efp_to_real

   pure function efp_plus(a, b) result(c)
      !! Exact bin-wise integer addition, REGULARISED (not merely
      !! carried).  Order-invariant BIT-FOR-BIT:
      !! `efp_plus(a,b)%v == efp_plus(b,a)%v`, and a running fold
      !! `acc = efp_plus(acc, x_i)` over any permutation of the `x_i`
      !! converges to the SAME raw bins (not merely the same reconstructed
      !! real) -- `test_efp_order_invariant` asserts this on `%v(:)`
      !! directly, per the plan's §9.2.
      !!
      !! `efp_carry` alone is NOT sufficient here: it only bounds each
      !! bin's MAGNITUDE (`|e(n)| < 2**P`), and a magnitude-bounded but
      !! mixed-sign bin vector is a REDUNDANT (non-unique) representation
      !! of a given value -- e.g. `[e2=1, e3=-(2**P-100)]` and
      !! `[e2=0, e3=100]` both represent the value `100`, but differ
      !! bit-for-bit.  Two different accumulation orders can land on
      !! different members of that redundant family even though the
      !! represented value is identical.  `efp_regularize` additionally
      !! forces every bin to share the overall sign, which -- like
      !! ordinary sign-magnitude fixed-radix representations -- IS
      !! unique for a given value, closing that gap.
      type(efp_t), intent(in) :: a, b
      type(efp_t) :: c
      c%v = a%v + b%v
      call efp_regularize(c%v)
   end function efp_plus

   pure function efp_minus(a, b) result(c)
      !! Exact bin-wise integer subtraction, regularised.  See `efp_plus`
      !! for why regularisation (not mere carry) is required for
      !! bit-for-bit order invariance.
      type(efp_t), intent(in) :: a, b
      type(efp_t) :: c
      c%v = a%v - b%v
      call efp_regularize(c%v)
   end function efp_minus

   pure function efp_real_diff(a, b) result(r)
      !! `real64 = efp_to_real(efp_minus(a, b))` -- the difference of two
      !! ~1e21-scale EFP totals resolved to the EFP quantum (`2**-3P`), NOT
      !! to `ulp(1e21)` as a double subtraction would give.  This is the
      !! fix documented in the plan's SS2.2: an implementer who converts
      !! both operands to `real64` FIRST and subtracts loses the whole
      !! benefit of this module.
      type(efp_t), intent(in) :: a, b
      real(real64) :: r
      r = efp_to_real(efp_minus(a, b))
   end function efp_real_diff

   pure subroutine efp_to_transport(list, buf)
      !! Pack a list of `efp_t` values into a flat `real64(6*n)` buffer for
      !! a single collective (`halo_allreduce_efp_list`).  Each bin is
      !! transported as an EXACTLY-representable double (see
      !! `EFP_MAX_RANKS`): `buf((i-1)*EFP_DIGITS + n) = real(list(i)%v(n))`.
      type(efp_t), intent(in) :: list(:)
      real(real64), intent(out) :: buf(:)
      integer :: i, n
      do i = 1, size(list)
         do n = 1, EFP_DIGITS
            buf((i - 1)*EFP_DIGITS + n) = real(list(i)%v(n), real64)
         end do
      end do
   end subroutine efp_to_transport

   pure subroutine efp_from_transport(buf, list, ok)
      !! Inverse of `efp_to_transport`: unpack a flat `real64(6*n)` buffer
      !! (post-collective, still exact integers as doubles) back into
      !! `efp_t` values, converting each bin back to `int64` and carrying.
      !! `ok = .false.` iff any unpacked double is not an exact integer
      !! (would indicate the transport-exactness bound was violated) --
      !! the caller (`halo_allreduce_efp_list`) turns that into a fail-loud
      !! `error stop`, never a silent truncation.
      real(real64), intent(in) :: buf(:)
      type(efp_t), intent(out) :: list(:)
      logical, intent(out) :: ok
      integer :: i, n
      real(real64) :: val

      ok = .true.
      do i = 1, size(list)
         do n = 1, EFP_DIGITS
            val = buf((i - 1)*EFP_DIGITS + n)
            if (val /= aint(val)) ok = .false.
            list(i)%v(n) = int(val, int64)
         end do
         ! Regularise (not merely carry) so the post-combine bins are the
         ! SAME canonical representation a single-rank EFP sum of the
         ! whole field would produce -- see `efp_plus`'s docstring for why
         ! carry alone is not bit-for-bit order-invariant.
         call efp_regularize(list(i)%v)
      end do
   end subroutine efp_from_transport

   pure function efp_bin1_within_transport_bound(a, nranks) result(ok)
      !! Bin 1 is NOT bounded by `efp_carry` (only bins 2..6 are -- see the
      !! module docstring), so the double-precision transport needs its
      !! own guard: after summing `nranks` local values, the partial sum
      !! in bin 1 must stay `<= 2**53 / nranks` for the MPI_SUM-on-doubles
      !! combine to remain exact (the analogue of MOM6's
      !! `prec_error = huge(1_int64) / num_PEs`, but for the int64-as-
      !! double transport rather than int64 transport).  `.false.` ⇒ the
      !! caller must `error stop`, never silently proceed.
      type(efp_t), intent(in) :: a
      integer, intent(in) :: nranks
      logical :: ok
      integer(int64) :: bound
      if (nranks <= 0) then
         ok = .true.
         return
      end if
      bound = (2_int64**53)/int(nranks, int64)
      ok = abs(a%v(1)) <= bound
   end function efp_bin1_within_transport_bound

end module rki_efp
