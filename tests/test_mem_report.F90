!! Unit tests for the memory-budget reporter (rki_mem_report).
!!
!! Cases:
!!   * `mem_host_rss_bytes()` is positive where /proc/self/status exists
!!     (Linux); on platforms without /proc (macOS dev boxes) it returns
!!     the documented -1 sentinel, so the RSS cases skip rather than fail.
!!   * Allocating + touching a ~64 MB array grows RSS by a detectable
!!     amount.  The check is deliberately loose (> 16 MB) because RSS
!!     reflects page-in semantics — only touched pages are resident, and
!!     the allocator/OS may reuse freed pages — so an exact figure would
!!     be flaky.  The array is freed afterwards.
!!   * The device queries must not crash on a CPU build: they return -1
!!     (no OpenACC runtime) OR a positive total (GPU build).
!!   * On a GPU build, after a device mapping: used > 0 and
!!     total >= used (device-wide sanity of `used = total - free`).
!!   * The full report sequence (budget -> actuals -> growth, i.e. the
!!     reconciliation lines) executes without crashing on BOTH builds.
!!   * The GB/MB formatter produces the documented spot values.
module test_mem_report
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: iso_fortran_env, only: int64
   use rki_constants, only: wp
   use rki_mem_report, only: mem_host_rss_bytes, mem_device_total_bytes, &
                             mem_device_free_bytes, mem_device_used_bytes, &
                             mem_format_bytes, mem_log_state_budget, &
                             mem_log_device_actuals, mem_log_device_growth, &
                             arr_bytes, mem_log_computed_budget, mem_log_computed_line
   implicit none
   private

   public :: collect_mem_report_tests

contains

   subroutine collect_mem_report_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("rss_positive_on_linux", test_rss_positive), &
                  new_unittest("rss_grows_with_allocation", test_rss_grows), &
                  new_unittest("device_query_no_crash", test_device_query), &
                  new_unittest("device_used_sane", test_device_used), &
                  new_unittest("reconciliation_logs_no_crash", test_reconciliation_logs), &
                  new_unittest("arr_bytes_exact_counts", test_arr_bytes_exact), &
                  new_unittest("arr_bytes_unallocated_zero", test_arr_bytes_unalloc), &
                  new_unittest("computed_budget_no_crash", test_computed_budget), &
                  new_unittest("format_bytes_spot_values", test_format_bytes) &
                  ]
   end subroutine collect_mem_report_tests

   subroutine test_arr_bytes_exact(error)
      !! arr_bytes must return size * element-width for allocated arrays of
      !! each supported type/rank (the exact counted footprint).
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: r2(:, :), r3(:, :, :)
      integer, allocatable :: i1(:)
      logical, allocatable :: l2(:, :)
      integer(int64) :: wp_w, i_w, l_w

      wp_w = int(storage_size(1.0_wp)/8, int64)
      i_w = int(storage_size(0)/8, int64)
      l_w = int(storage_size(.true.)/8, int64)

      allocate (r2(10, 20))
      allocate (r3(4, 5, 6))
      allocate (i1(7))
      allocate (l2(3, 3))

      call check(error, arr_bytes(r2) == 200_int64*wp_w, "r2 = 200 elements * wp width")
      if (allocated(error)) return
      call check(error, arr_bytes(r3) == 120_int64*wp_w, "r3 = 120 elements * wp width")
      if (allocated(error)) return
      call check(error, arr_bytes(i1) == 7_int64*i_w, "i1 = 7 elements * int width")
      if (allocated(error)) return
      call check(error, arr_bytes(l2) == 9_int64*l_w, "l2 = 9 elements * logical width")
   end subroutine test_arr_bytes_exact

   subroutine test_arr_bytes_unalloc(error)
      !! An unallocated array counts 0 — the property that lets a gated-off
      !! slot contribute nothing to a state's bytes() total.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: r2(:, :), r3(:, :, :)
      integer, allocatable :: i2(:, :)

      call check(error, arr_bytes(r2) == 0_int64, "unallocated real(:,:) counts 0")
      if (allocated(error)) return
      call check(error, arr_bytes(r3) == 0_int64, "unallocated real(:,:,:) counts 0")
      if (allocated(error)) return
      call check(error, arr_bytes(i2) == 0_int64, "unallocated integer(:,:) counts 0")
   end subroutine test_arr_bytes_unalloc

   subroutine test_computed_budget(error)
      !! The upfront counted-budget loggers must run without crashing; the
      !! per-component line is skipped for a 0-byte (gated-off) slot.
      type(error_type), allocatable, intent(out) :: error

      call mem_log_computed_budget("test state", 128_int64*1024_int64*1024_int64)
      call mem_log_computed_line("barotropic", 64_int64*1024_int64*1024_int64)
      call mem_log_computed_line("gated-off closure", 0_int64)  ! skipped, no crash
      call check(error, .true.)
   end subroutine test_computed_budget

   subroutine test_rss_positive(error)
      type(error_type), allocatable, intent(out) :: error
      integer(int64) :: rss

      rss = mem_host_rss_bytes()
      ! Where /proc/self/status exists (Linux) VmRSS parses to a positive
      ! byte count; on platforms without /proc (macOS) the query returns the
      ! documented -1 sentinel and there is nothing to assert — skip.
      if (rss == -1_int64) return
      call check(error, rss > 0_int64, &
                 "VmRSS should parse to a positive byte count where /proc is available")
   end subroutine test_rss_positive

   subroutine test_rss_grows(error)
      type(error_type), allocatable, intent(out) :: error
      integer(int64) :: rss_before, rss_after, growth
      integer(int64), parameter :: N = 8_int64*1024_int64*1024_int64  ! 64 MB of int64.
      integer(int64), parameter :: MIN_GROWTH = 16_int64*1024_int64*1024_int64
      integer(int64), allocatable :: big(:)
      integer(int64) :: i

      rss_before = mem_host_rss_bytes()
      ! RSS unsupported (no /proc, e.g. macOS) — the growth assertion is
      ! Linux-only; skip on the -1 sentinel rather than fail.
      if (rss_before == -1_int64) return
      allocate (big(N))
      ! Touch every page so the kernel actually backs it with RAM.
      do i = 1_int64, N
         big(i) = i
      end do
      ! Defeat dead-store elimination: force a read of the buffer.
      if (big(N) /= N) then
         call check(error, .false., "touched buffer readback mismatch")
         return
      end if
      rss_after = mem_host_rss_bytes()
      growth = rss_after - rss_before
      deallocate (big)

      call check(error, growth > MIN_GROWTH, &
                 "RSS should grow by > 16 MB after touching a 64 MB buffer")
   end subroutine test_rss_grows

   subroutine test_device_query(error)
      type(error_type), allocatable, intent(out) :: error
      integer(int64) :: total, free

      ! On a CPU build these return -1; on a GPU build a positive total.
      ! Must never crash either way.
      total = mem_device_total_bytes()
      free = mem_device_free_bytes()
      call check(error, total == -1_int64 .or. total > 0_int64, &
                 "device total must be -1 (CPU) or positive (GPU)")
      if (allocated(error)) return
      call check(error, free == -1_int64 .or. free > 0_int64, &
                 "device free must be -1 (CPU) or positive (GPU)")
   end subroutine test_device_query

   subroutine test_device_used(error)
      !! On a GPU build, after mapping an array (forces context + a live
      !! allocation): used > 0 and total >= used.  On a CPU build the
      !! directive is inert and the query returns -1 — asserts skipped.
      type(error_type), allocatable, intent(out) :: error
      integer(int64) :: used, total
      real(wp), allocatable :: probe(:)

      allocate (probe(1024*1024))  ! 4 MB
      probe = 0.0_wp
      ! [acc->omp-cpu] dropped device directive (host no-op): enter data

      used = mem_device_used_bytes()
      total = mem_device_total_bytes()

      if (total > 0_int64) then
         ! GPU build: the mapping above forced CUDA context creation, so
         ! device usage is strictly positive and bounded by the total.
         call check(error, used > 0_int64, &
                    "device used must be > 0 on a GPU build after a mapping")
         if (.not. allocated(error)) then
            call check(error, total >= used, "device total must be >= used")
         end if
      else
         ! CPU build: no OpenACC runtime — the sentinel is the contract.
         call check(error, used == -1_int64, &
                    "device used must be -1 when no OpenACC runtime is present")
      end if

      ! [acc->omp-cpu] dropped device directive (host no-op): exit data
      deallocate (probe)
   end subroutine test_device_used

   subroutine test_reconciliation_logs(error)
      !! The full report sequence a run emits — budget (latches
      !! free-before-mapping) -> mapping -> actuals (measured delta +
      !! post-setup baseline) -> growth (reconciliation, loud + quiet) —
      !! must execute without crashing on both CPU and GPU builds.
      !! Logger text isn't captured (test-drive has no stream capture);
      !! this guards the code paths, incl. the device-absent early
      !! returns.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: workspace(:)

      call mem_log_state_budget("test state", 4_int64*1024_int64*1024_int64)

      allocate (workspace(2*1024*1024))  ! 8 MB
      workspace = 1.0_wp
      ! [acc->omp-cpu] dropped device directive (host no-op): enter data

      call mem_log_device_actuals("test state mapped")
      call mem_log_device_growth("test status", quiet=.true.)
      call mem_log_device_growth("test end of run")

      ! [acc->omp-cpu] dropped device directive (host no-op): exit data
      deallocate (workspace)

      ! Reaching here without a crash is the assertion.
      call check(error, .true.)
   end subroutine test_reconciliation_logs

   subroutine test_format_bytes(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, mem_format_bytes(8589934592_int64) == "8.00 GB", &
                 "8 GiB should format as 8.00 GB")
      if (allocated(error)) return
      call check(error, mem_format_bytes(536870912_int64) == "512.00 MB", &
                 "512 MiB should format as 512.00 MB")
      if (allocated(error)) return
      ! Just under 1 GiB stays in MB.
      call check(error, mem_format_bytes(1073741823_int64) == "1024.00 MB", &
                 "just under 1 GiB should format in MB")
   end subroutine test_format_bytes

end module test_mem_report
