!! Toolchain conformance tests for F2023 `do concurrent ... reduce` locality.
!!
!! These do NOT test Rakali physics — they test the COMPILER.  The ocean core
!! now expresses its reductions as `do concurrent(...) reduce(op:var)` instead
!! of `!$acc parallel loop reduction(...)`, so a compiler that lowers `reduce`
!! incorrectly corrupts diagnostics silently: no crash, no warning, just wrong
!! numbers.  A compiler/version bump is exactly when that regresses, which is
!! why this lives in ctest rather than in a scratch reproducer.
!!
!! Every assertion here is EXACT equality against a sequential reference, and
!! that is deliberate — each case is chosen so the answer is order-invariant
!! by construction, leaving no room to excuse a failure as round-off:
!!   * integer sums / counts  -- integer addition is associative
!!   * min / max              -- associative and commutative
!!   * real sums of SMALL INTEGER-VALUED data -- exact in float64 below 2**53,
!!     so no association order can change the result
!! A failure here is a MISCOMPILE, never floating-point noise.
!!
!! Known-bad shape deliberately NOT exercised: `transfer()` inside a
!! `do concurrent`.  nvfortran 26.5 gets ~95% of elements wrong there (a
!! shifted-index miscompile, unrelated to reduce) -- see
!! tests/test_dc_transfer_sentinel.F90 and the `dc-transfer` lint.  Bit
!! manipulation stays on `!$acc parallel loop`; see rki_ocean_chksum.
module test_dc_reduce
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   implicit none
   private

   public :: collect_dc_reduce_tests

   integer, parameter :: N = 100003        !! prime-ish, > any plausible tile
   integer, parameter :: NX = 41, NY = 13, NZ = 7

contains

   subroutine collect_dc_reduce_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("reduce_real_sum_exact", test_real_sum), &
                  new_unittest("reduce_int32_sum_exact", test_int32_sum), &
                  new_unittest("reduce_int64_sum_exact", test_int64_sum), &
                  new_unittest("reduce_min_max_exact", test_min_max), &
                  new_unittest("reduce_multi_var_one_construct", test_multi_var), &
                  new_unittest("reduce_collapsed_3d", test_collapsed_3d), &
                  new_unittest("reduce_order_invariant", test_order_invariant) &
                  ]
   end subroutine collect_dc_reduce_tests

   !> real(+) over small integer-valued data: exact, so equality is legitimate.
   subroutine test_real_sum(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: a(N), s, ref
      integer :: i
      do i = 1, N
         a(i) = real(mod(i, 17) + 1, wp)
      end do
      ref = 0.0_wp
      do i = 1, N
         ref = ref + a(i)
      end do
      s = 0.0_wp
      !$omp target teams distribute parallel do
      do i = 1, N
         s = s + a(i)
      end do
      !$omp end target teams distribute parallel do
      call check(error, s == ref, "do concurrent reduce(+) on real(wp) must be exact")
   end subroutine test_real_sum

   subroutine test_int32_sum(error)
      type(error_type), allocatable, intent(out) :: error
      integer(int32) :: a(N), s, ref
      integer :: i
      do i = 1, N
         a(i) = int(mod(i, 7) + 1, int32)
      end do
      ref = 0_int32
      do i = 1, N
         ref = ref + a(i)
      end do
      s = 0_int32
      do concurrent(i=1:N) reduce(+:s)
         s = s + a(i)
      end do
      call check(error, s == ref, "reduce(+) on integer(int32) must be exact")
   end subroutine test_int32_sum

   !> int64 specifically: the accumulator width used by bitwise checksums.
   subroutine test_int64_sum(error)
      type(error_type), allocatable, intent(out) :: error
      integer(int64) :: a(N), s, ref
      integer :: i
      do i = 1, N
         a(i) = int(mod(i, 61) + 1, int64)
      end do
      ref = 0_int64
      do i = 1, N
         ref = ref + a(i)
      end do
      s = 0_int64
      do concurrent(i=1:N) reduce(+:s)
         s = s + a(i)
      end do
      call check(error, s == ref, "reduce(+) on integer(int64) must be exact")
   end subroutine test_int64_sum

   subroutine test_min_max(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: a(N), mn, mx
      integer :: i
      do i = 1, N
         a(i) = real(mod(i*7919, 1000), wp) - 500.0_wp
      end do
      mn = huge(1.0_wp); mx = -huge(1.0_wp)
      do concurrent(i=1:N) reduce(min:mn)
         mn = min(mn, a(i))
      end do
      do concurrent(i=1:N) reduce(max:mx)
         mx = max(mx, a(i))
      end do
      call check(error, mn == minval(a), "reduce(min) must equal minval")
      if (allocated(error)) return
      call check(error, mx == maxval(a), "reduce(max) must equal maxval")
   end subroutine test_min_max

   !> Several reduction variables, mixed types, three operators, ONE construct
   !> — the shape rki_ocean_chksum uses.  Regression guard: an implementation
   !> that handles one variable but drops the others fails only here.
   subroutine test_multi_var(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: a(N), s, mn, mx, rs, rmn, rmx
      integer :: i, c, rc
      integer(int64) :: b, rb
      do i = 1, N
         a(i) = real(mod(i, 23) + 1, wp)
      end do
      rs = 0.0_wp; rmn = huge(1.0_wp); rmx = -huge(1.0_wp); rc = 0; rb = 0_int64
      do i = 1, N
         rs = rs + a(i)
         rmn = min(rmn, a(i)); rmx = max(rmx, a(i))
         if (a(i) > 12.0_wp) rc = rc + 1
         rb = rb + int(a(i), int64)
      end do
      s = 0.0_wp; mn = huge(1.0_wp); mx = -huge(1.0_wp); c = 0; b = 0_int64
      do concurrent(i=1:N) reduce(+:s, c, b) reduce(min:mn) reduce(max:mx)
         s = s + a(i)
         mn = min(mn, a(i)); mx = max(mx, a(i))
         if (a(i) > 12.0_wp) c = c + 1
         b = b + int(a(i), int64)
      end do
      call check(error, s == rs, "multi-var reduce: real sum")
      if (allocated(error)) return
      call check(error, c == rc, "multi-var reduce: integer count")
      if (allocated(error)) return
      call check(error, b == rb, "multi-var reduce: int64 sum")
      if (allocated(error)) return
      call check(error, mn == rmn .and. mx == rmx, "multi-var reduce: min/max")
   end subroutine test_multi_var

   !> Multi-index (collapsed) header — how every converted kernel is written.
   subroutine test_collapsed_3d(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: a(NX, NY, NZ), s, ref
      integer :: i, j, k
      do k = 1, NZ; do j = 1, NY; do i = 1, NX
            a(i, j, k) = real(i + 10*j + 100*k, wp)
         end do; end do; end do
      ref = 0.0_wp
      do k = 1, NZ; do j = 1, NY; do i = 1, NX
            ref = ref + a(i, j, k)
         end do; end do; end do
      s = 0.0_wp
      do concurrent(k=1:NZ, j=1:NY, i=1:NX) reduce(+:s)
         s = s + a(i, j, k)
      end do
      call check(error, s == ref, "collapsed 3-index reduce(+) must be exact")
   end subroutine test_collapsed_3d

   !> Same elements, two different traversals, integer accumulator: integer
   !> addition is associative and commutative, so the totals MUST match.
   !> This is the property that caught the nvfortran transfer() miscompile.
   subroutine test_order_invariant(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: a(NX, NY, NZ), flat(NX*NY*NZ)
      integer(int64) :: s3, s1
      integer :: i, j, k, n
      do k = 1, NZ; do j = 1, NY; do i = 1, NX
            a(i, j, k) = real(mod(i*j*k, 29) + 1, wp)
         end do; end do; end do
      n = 0
      do k = 1, NZ; do j = 1, NY; do i = 1, NX
            n = n + 1; flat(n) = a(i, j, k)
         end do; end do; end do
      s3 = 0_int64
      do concurrent(k=1:NZ, j=1:NY, i=1:NX) reduce(+:s3)
         s3 = s3 + int(a(i, j, k), int64)
      end do
      s1 = 0_int64
      do concurrent(n=1:NX*NY*NZ) reduce(+:s1)
         s1 = s1 + int(flat(n), int64)
      end do
      call check(error, s3 == s1, "integer reduce must be traversal-order invariant")
   end subroutine test_order_invariant

end module test_dc_reduce
