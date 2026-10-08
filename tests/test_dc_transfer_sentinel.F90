!! SENTINEL: `transfer()` inside a `do concurrent` must be correct.
!!
!! This test asserts the RIGHT behaviour.  It is EXPECTED TO FAIL on NVHPC
!! GPU builds (`-stdpar=gpu`) today, and tests/CMakeLists.txt marks it
!! `WILL_FAIL TRUE` for that configuration only.  So:
!!
!!   * NVHPC GPU build, bug present -> binary fails -> WILL_FAIL -> ctest PASSES (quiet)
!!   * NVHPC GPU build, bug FIXED   -> binary passes -> WILL_FAIL -> ctest FAILS
!!         ^ that is the notification.  When this goes red on an NVHPC GPU
!!           build, NVIDIA has fixed it: drop the WILL_FAIL, delete the
!!           `dc-transfer` lint, and rki_ocean_chksum can move to
!!           `do concurrent`.
!!   * every other configuration    -> no WILL_FAIL -> must pass normally.
!!         That includes an NVHPC CPU build, where `do concurrent` lowers to
!!         ordinary serial loops and is correct, so a regression there is
!!         caught the ordinary way.
!!
!! The bug (nvfortran 26.5, -stdpar=gpu -gpu=cc70,mem:separate): a plain
!! elementwise `got(i) = transfer(x(i), 0_int64)` in a `do concurrent` gets
!! ~95% of elements wrong (94645/100000), reading a SHIFTED element --
!! `got 0x4008000000000000` (3.0) where `0x4000000000000000` (2.0) was
!! expected.  No reduction is involved; it is an index/scoping miscompile in
!! the do-concurrent lowering.  The identical expression under
!! `!$acc parallel loop` is CORRECT, which is why rki_ocean_chksum stays on
!! OpenACC.  This test IS the reproducer -- it needs no reduction, no
!! int64 accumulator and no popcnt, only `transfer` in a `do concurrent`.
module test_dc_transfer_sentinel
   use, intrinsic :: iso_fortran_env, only: int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   implicit none
   private

   public :: collect_dc_transfer_sentinel_tests

   integer, parameter :: N = 100000

contains

   subroutine collect_dc_transfer_sentinel_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("transfer_in_do_concurrent_is_correct", test_transfer_dc) &
                  ]
   end subroutine collect_dc_transfer_sentinel_tests

   subroutine test_transfer_dc(error)
      !! Bit-identical reinterpretation, computed sequentially and in a
      !! `do concurrent`.  `transfer` is a pure bit reinterpretation, so the
      !! two MUST agree element-for-element -- there is no ordering,
      !! association or rounding freedom to appeal to.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: x(N)
      integer(int64) :: got(N), want(N)
      integer :: i, bad

      do i = 1, N
         x(i) = real(mod(i, 17) + 1, wp)
      end do
      do i = 1, N
         want(i) = transfer(x(i), 0_int64)
      end do

      got = 0_int64
      !$omp parallel do
      do i = 1, N
         got(i) = transfer(x(i), 0_int64)   ! dc-transfer-ok: this IS the probe
      end do
      !$omp end parallel do

      bad = count(got /= want)
      call check(error, bad == 0, "transfer() in do concurrent must be bit-exact")
   end subroutine test_transfer_dc

end module test_dc_transfer_sentinel
