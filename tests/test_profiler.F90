!! Unit tests for the wall-clock profiler
module test_profiler
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_profiler, only: profiler_init, profiler_end, &
                           profiler_start, profiler_stop, &
                           profiler_enable, profiler_disable, &
                           profiler_reset, profiler_report, &
                           profiler_get_time
   implicit none
   private

   public :: collect_profiler_tests

contains

   subroutine collect_profiler_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("start_stop_accumulates", test_start_stop_accumulates), &
                  new_unittest("multiple_regions", test_multiple_regions), &
                  new_unittest("repeated_calls_accumulate", test_repeated_calls), &
                  new_unittest("disable_skips_timing", test_disable_skips_timing), &
                  new_unittest("reset_clears_times", test_reset_clears_times), &
                  new_unittest("report_runs", test_report_runs), &
                  new_unittest("get_time_unknown_region_zero", test_unknown_region), &
                  new_unittest("nvtx_only_flag", test_nvtx_only_flag) &
                  ]
   end subroutine collect_profiler_tests

   ! Spin until system_clock advances by a few ticks, so the profiler always
   ! records a non-zero interval. `dummy` is `volatile` so NVHPC -fast -O3
   ! cannot dead-code-eliminate the spin (its result is otherwise unread).
   subroutine spin_work(out)
      use, intrinsic :: iso_fortran_env, only: int64
      real(wp), intent(out) :: out
      volatile :: out
      integer(int64) :: t0, t1, rate
      integer :: i
      real(wp) :: s
      s = 0.0_wp
      call system_clock(t0, rate)
      i = 0
      do
         i = i + 1
         s = s + sqrt(real(i, wp))
         if (iand(i, 1023) == 0) then
            call system_clock(t1)
            if (t1 - t0 >= 4 .or. i >= 10000000) exit
         end if
      end do
      out = s
   end subroutine spin_work

   subroutine test_start_stop_accumulates(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: t, dummy
      checks: block

         call profiler_init(.true.)
         call profiler_start("region_a")
         call spin_work(dummy)
         call profiler_stop("region_a")

         t = profiler_get_time("region_a")
         call check(error, t > 0.0_wp, "expected positive accumulated time")
         if (allocated(error)) exit checks

      end block checks
      continue
      call profiler_end()
   end subroutine test_start_stop_accumulates

   subroutine test_multiple_regions(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: ta, tb, dummy
      checks: block

         call profiler_init(.true.)
         call profiler_start("region_a")
         call spin_work(dummy)
         call profiler_stop("region_a")

         call profiler_start("region_b")
         call spin_work(dummy)
         call profiler_stop("region_b")

         ta = profiler_get_time("region_a")
         tb = profiler_get_time("region_b")
         call check(error, ta > 0.0_wp, "region_a should have time")
         if (allocated(error)) exit checks
         call check(error, tb > 0.0_wp, "region_b should have time")
         if (allocated(error)) exit checks

      end block checks
      continue
      call profiler_end()
   end subroutine test_multiple_regions

   subroutine test_repeated_calls(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: t1, t2, dummy
      checks: block

         call profiler_init(.true.)
         call profiler_start("region_x")
         call spin_work(dummy)
         call profiler_stop("region_x")
         t1 = profiler_get_time("region_x")

         call profiler_start("region_x")
         call spin_work(dummy)
         call profiler_stop("region_x")
         t2 = profiler_get_time("region_x")

         call check(error, t2 >= t1, "second call must accumulate (>=)")
         if (allocated(error)) exit checks

         call check(error, t2 > 0.0_wp, "total must be positive")
         if (allocated(error)) exit checks

      end block checks
      continue
      call profiler_end()
   end subroutine test_repeated_calls

   subroutine test_disable_skips_timing(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: t, dummy
      checks: block

         call profiler_init(.true.)
         call profiler_disable()

         call profiler_start("region_off")
         call spin_work(dummy)
         call profiler_stop("region_off")

         t = profiler_get_time("region_off")
         call check(error, t == 0.0_wp, "disabled profiler must not record any time")
         if (allocated(error)) exit checks

         call profiler_enable()
         call profiler_start("region_on")
         call spin_work(dummy)
         call profiler_stop("region_on")
         t = profiler_get_time("region_on")
         call check(error, t > 0.0_wp, "re-enabled profiler must record time")
         if (allocated(error)) exit checks

      end block checks
      continue
      call profiler_end()
   end subroutine test_disable_skips_timing

   subroutine test_reset_clears_times(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: t, dummy
      checks: block

         call profiler_init(.true.)
         call profiler_start("region_r")
         call spin_work(dummy)
         call profiler_stop("region_r")

         call profiler_reset()
         t = profiler_get_time("region_r")
         call check(error, t == 0.0_wp, "reset must zero accumulated time")
         if (allocated(error)) exit checks

      end block checks
      continue
      call profiler_end()
   end subroutine test_reset_clears_times

   subroutine test_report_runs(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dummy

      call profiler_init(.true.)
      call profiler_start("rpt")
      call spin_work(dummy)
      call profiler_stop("rpt")

      ! Should run without crashing both with and without an explicit title
      call profiler_report()
      call profiler_report("test_title")
      call profiler_report("with_root", root_region="rpt")

      call check(error, .true., "report should not crash")

      call profiler_end()
   end subroutine test_report_runs

   subroutine test_unknown_region(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: t

      call profiler_init(.true.)
      t = profiler_get_time("never_started")
      call check(error, t == 0.0_wp, "unknown region must return 0")
      call profiler_end()
   end subroutine test_unknown_region

   subroutine test_nvtx_only_flag(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: t, dummy

      call profiler_init(.true.)
      call profiler_start("nvtx_region", nvtx_only=.true.)
      call spin_work(dummy)
      call profiler_stop("nvtx_region")

      t = profiler_get_time("nvtx_region")
      call check(error, t > 0.0_wp, "nvtx_only region still records timing")

      ! Report should still run (nvtx_only entries are skipped from the table
      ! but the call must not crash).
      call profiler_report("with_nvtx_only")

      call profiler_end()
   end subroutine test_nvtx_only_flag

end module test_profiler
