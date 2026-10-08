!! Setup-time memory-budget reporting for the device-mapped state.
module rki_mem_report
   !! Free procedures that estimate, log, and verify the device memory
   !! footprint of the simulation state **before** the OpenACC
   !! `enter_data` mapping is attempted — so an over-budget run reports
   !! a readable diagnostic instead of dying opaquely inside the CUDA
   !! runtime (the motivating failure: a 10M-cell ocean config that
   !! exhausted device memory at `enter_data` with no setup-time hint).
   !!
   !! Design decision (no hand-maintained per-slot byte counts — those
   !! rot the moment a slot is added): the ONLY trustworthy figures are
   !! the OpenACC runtime's own free-memory queries.  The device
   !! footprint of the state mapping is **measured**, not estimated, as
   !! the drop in free device memory across `enter_data`
   !! (`mem_log_state_budget` snapshots free-before;
   !! `mem_log_device_actuals` reports the delta).  Host RSS growth
   !! across state init is still logged, but ONLY as an advisory
   !! host-side figure + a conservative upper bound for the pre-flight
   !! OOM warning — it includes host-only allocations (NetCDF/HDF5
   !! library pages, I/O buffers, host mirrors of the device arrays) and
   !! measured 6-31x larger than the true device map (2026-07 audit), so
   !! it must never be read as a device estimate.
   !!
   !! Three structural offsets to keep in mind when reconciling these
   !! figures against `nvidia-smi` (audited 2026-07 on V100):
   !!   1. **CUDA context + CUBIN load**: ~300-500 MB per process,
   !!      created at the first OpenACC API call / directive.  It is
   !!      inside `used = total - free` and inside nvidia-smi's
   !!      per-process figure, but is NOT state — a tiny config still
   !!      shows ~300 MB in use.
   !!   2. **Device-wide accounting + pool granularity**:
   !!      `used = total - free` counts EVERY process on the device
   !!      (shared-device activity perturbs it and the measured mapping
   !!      delta), and the NVHPC pool allocator rounds small maps up by
   !!      a few MB.
   !!   3. **One-off CUDA runtime/kernel reservation at first launch**:
   !!      the first time-loop step grows per-process device usage by
   !!      GBs (measured: ~3.0 GB coastal, ~8.6 GB ocean on V100) —
   !!      constant thereafter, independent of step count and of
   !!      `NV_ACC_MEM_MANAGE`, and NOT proportional to grid size.
   !!      Consistent with the driver's local-memory/stack backing for
   !!      stack-heavy kernels (fixed-size `local()` arrays), which
   !!      scales with the binary's kernels x device SM count.  It is
   !!      real per-process usage (visible in nvidia-smi) but is not a
   !!      leak; the growth WARNING therefore measures against the
   !!      steady post-first-sample baseline, not the post-`enter_data`
   !!      one.
   !!
   !! Post-setup growth (`mem_log_device_growth`) makes drift
   !! self-announcing: kernels that lazy-allocate device workspaces on
   !! first call (the coastal `*_workspace_ensure` pattern — measured at
   !! ~33% of the coastal device footprint) land AFTER `enter_data`, so
   !! the periodic status re-checks and the end-of-run summary reports
   !! the growth, warning above `GROWTH_WARN_FRAC`.
   !!
   !! Device queries use the OpenACC property API
   !! (`acc_get_property` with `acc_property_memory` /
   !! `acc_property_free_memory`).  They are compiled only when the
   !! OpenACC runtime module is available (`RAKALI_HAS_OPENACC_RUNTIME`,
   !! set by `cmake/compiler_flags.cmake` for NVHPC + Cray — covers both
   !! the `-acc=gpu` and `-mp=gpu` NVHPC backends).  On every other
   !! toolchain (gfortran, ifx) the device functions return `-1` and the
   !! loggers degrade to the host-only line.
   use, intrinsic :: iso_fortran_env, only: int64
   use rki_constants, only: wp
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
#ifdef RAKALI_HAS_OPENACC_RUNTIME
   use openacc, only: acc_get_device_num, acc_get_device_type, &
                      acc_get_property, acc_property_memory, &
                      acc_property_free_memory
   use, intrinsic :: iso_c_binding, only: c_size_t
#endif
   implicit none
   private

   public :: mem_host_rss_bytes
   public :: mem_device_total_bytes, mem_device_free_bytes
   public :: mem_device_used_bytes
   public :: mem_log_state_budget, mem_log_device_actuals
   public :: mem_log_device_growth
   public :: mem_format_bytes
   public :: arr_bytes
   public :: mem_set_counted_budget, mem_log_computed_budget, mem_log_computed_line

   interface arr_bytes
      !! Byte footprint of an allocatable array (0 when unallocated), so a
      !! gated-off slot naturally contributes nothing to a counted total.
      !! The COUNTED companion to the runtime's MEASURED free-memory
      !! delta: each `*_state_t` sums its own arrays through this, the
      !! driver reports the total **before** `enter_data` (works on CPU
      !! builds too, where there is no device query), and
      !! `mem_log_device_actuals` reconciles the count against the
      !! measured mapping — a state array added without a matching
      !! `bytes()` term makes the measured map exceed the count and
      !! self-announces the drift (the anti-rot guard for the otherwise
      !! rot-prone per-slot count).
      module procedure arr_bytes_r1, arr_bytes_r2, arr_bytes_r3, arr_bytes_r4
      module procedure arr_bytes_i1, arr_bytes_i2, arr_bytes_i3
      module procedure arr_bytes_l1, arr_bytes_l2, arr_bytes_l3
   end interface arr_bytes

   real(wp), parameter :: COUNT_RECONCILE_FRAC = 0.25_wp
      !! The measured `enter_data` mapping may exceed the COUNTED state
      !! footprint by this fraction (NVHPC pool granularity + any
      !! device-only scratch mapped alongside the state) before the
      !! reconciliation warns.  Overshoot past it means the device mapped
      !! significantly more than every `bytes()` term accounts for — a
      !! state array is almost certainly missing from the count.

   integer(int64), parameter :: COUNT_RECONCILE_FLOOR = 64_int64*1024_int64*1024_int64
      !! Absolute slack (64 MiB) added to the reconciliation tolerance so
      !! tiny configs — where pool rounding dwarfs the relative margin —
      !! never trip the drift warning.

   integer(int64), parameter :: BYTES_PER_KIB = 1024_int64
   integer(int64), parameter :: BYTES_PER_MIB = 1024_int64*1024_int64
   integer(int64), parameter :: BYTES_PER_GIB = 1024_int64*1024_int64*1024_int64

   real(wp), parameter :: GROWTH_WARN_FRAC = 0.15_wp
      !! Device-memory growth above this fraction of the STEADY
      !! post-setup baseline (first post-setup sample — after the
      !! one-off first-launch reservation and lazy workspaces) triggers
      !! a one-shot logger warning — progressive growth (a leak)
      !! becomes self-announcing instead of silently skewing the
      !! setup-time report.

   ! Snapshots latched by the rank-0 logging path (module state, same
   ! lifetime/ownership model as the profiler's timers).  `-1` = not yet
   ! sampled / no device.
   integer(int64) :: free_before_map = -1_int64
      !! Free device bytes just before `enter_data` (latched by
      !! `mem_log_state_budget`); baseline for the measured mapping delta.
   integer(int64) :: used_after_map = -1_int64
      !! Device bytes in use just after `enter_data` (latched by
      !! `mem_log_device_actuals`); baseline for the "grew ~X after
      !! setup" attribution line.
   integer(int64) :: used_steady = -1_int64
      !! Device bytes in use at the FIRST post-setup growth sample
      !! (first status fire or end of run) — after the one-off
      !! first-launch runtime reservation and the lazy kernel
      !! workspaces have landed.  Baseline for the leak WARNING:
      !! growth past this is progressive, not the expected one-offs.
   logical :: growth_warned = .false.
      !! One-shot latch so the growth warning doesn't spam every status
      !! fire once the threshold is crossed.
   integer(int64) :: counted_state_bytes = -1_int64
      !! Counted state-array footprint latched by `mem_log_computed_budget`
      !! (summed from every `bytes()` term BEFORE `enter_data`); reconciled
      !! against the measured mapping in `mem_log_device_actuals`.  `-1` =
      !! not counted this run.

contains

   function mem_host_rss_bytes() result(bytes)
      !! Resident set size of the current process, in bytes, parsed from
      !! `VmRSS` in `/proc/self/status` (reported in kB → scaled to
      !! bytes).  Returns `-1` if the file is unreadable or the field is
      !! absent (non-Linux platforms).  Never aborts.
      integer(int64) :: bytes
      integer :: unit, ios
      character(len=256) :: line
      character(len=256) :: iomsg
      character(len=32) :: kb_str
      integer(int64) :: kb

      bytes = -1_int64
      open (newunit=unit, file="/proc/self/status", status="old", &
            action="read", iostat=ios, iomsg=iomsg)
      if (ios /= 0) return

      do
         read (unit, "(a)", iostat=ios, iomsg=iomsg) line
         if (ios /= 0) exit
         if (line(1:6) == "VmRSS:") then
            ! Layout: "VmRSS:\t   12345 kB" — strip the label + "kB".
            kb_str = adjustl(line(7:))
            ! Drop the trailing " kB" unit by reading the leading integer.
            read (kb_str, *, iostat=ios, iomsg=iomsg) kb
            if (ios == 0) bytes = kb*BYTES_PER_KIB
            exit
         end if
      end do

      close (unit)
   end function mem_host_rss_bytes

   function mem_device_total_bytes() result(bytes)
      !! Total memory of the device that the OpenACC runtime is bound to,
      !! in bytes.  Returns `-1` when no OpenACC runtime is present
      !! (CPU / gfortran / ifx builds).
      integer(int64) :: bytes
#ifdef RAKALI_HAS_OPENACC_RUNTIME
      integer(c_size_t) :: prop
      prop = acc_get_property(acc_get_device_num(acc_get_device_type()), &
                              acc_get_device_type(), acc_property_memory)
      bytes = int(prop, int64)
#else
      bytes = -1_int64
#endif
   end function mem_device_total_bytes

   function mem_device_free_bytes() result(bytes)
      !! Free memory on the bound device, in bytes.  Returns `-1` when no
      !! OpenACC runtime is present.
      integer(int64) :: bytes
#ifdef RAKALI_HAS_OPENACC_RUNTIME
      integer(c_size_t) :: prop
      prop = acc_get_property(acc_get_device_num(acc_get_device_type()), &
                              acc_get_device_type(), acc_property_free_memory)
      bytes = int(prop, int64)
#else
      bytes = -1_int64
#endif
   end function mem_device_free_bytes

   function mem_device_used_bytes() result(bytes)
      !! Device memory currently in use, in bytes (`total - free`).
      !! DEVICE-WIDE: includes this process's CUDA context (~300-500 MB)
      !! and every other process on a shared device — see the module
      !! header's structural offsets.  Returns `-1` when no OpenACC
      !! runtime is present.
      integer(int64) :: bytes
      integer(int64) :: free_b, total_b

      bytes = -1_int64
      free_b = mem_device_free_bytes()
      total_b = mem_device_total_bytes()
      if (free_b < 0_int64 .or. total_b < 0_int64) return
      bytes = total_b - free_b
   end function mem_device_used_bytes

   pure function mem_format_bytes(bytes) result(str)
      !! Human-readable byte count: two decimals, GB at/above 1 GiB,
      !! MB below.  Used for every figure in the budget log so the
      !! units stay consistent.
      integer(int64), intent(in) :: bytes
      character(len=:), allocatable :: str
      real(wp) :: gib, mib

      if (bytes >= BYTES_PER_GIB) then
         gib = real(bytes, wp)/real(BYTES_PER_GIB, wp)
         str = two_dp(gib)//" GB"
      else
         mib = real(bytes, wp)/real(BYTES_PER_MIB, wp)
         str = two_dp(mib)//" MB"
      end if
   end function mem_format_bytes

   pure function two_dp(x) result(str)
      !! Format a real to exactly two decimal places without leading
      !! blanks (`write` with an `f` edit descriptor, then `adjustl` +
      !! trim).
      real(wp), intent(in) :: x
      character(len=:), allocatable :: str
      character(len=32) :: buf
      write (buf, "(f0.2)") x
      ! `f0.2` of a value < 1 emits ".50" on some compilers; prefix a 0.
      if (buf(1:1) == ".") then
         str = "0"//trim(adjustl(buf))
      else
         str = trim(adjustl(buf))
      end if
   end function two_dp

   ! --- Counted footprint: byte size of an allocatable (0 if unallocated) ---
   ! One specific per (type, rank).  `size(a, kind=int64)` is only
   ! referenced inside the `allocated` guard, so an unallocated slot is
   ! never touched.  The element width comes from a same-kind literal
   ! (`storage_size` of the array itself would need it allocated).

   pure function arr_bytes_r1(a) result(b)
      real(wp), allocatable, intent(in) :: a(:)
      integer(int64) :: b
      b = 0_int64
      if (allocated(a)) b = size(a, kind=int64)*int(storage_size(1.0_wp)/8, int64)
   end function arr_bytes_r1

   pure function arr_bytes_r2(a) result(b)
      real(wp), allocatable, intent(in) :: a(:, :)
      integer(int64) :: b
      b = 0_int64
      if (allocated(a)) b = size(a, kind=int64)*int(storage_size(1.0_wp)/8, int64)
   end function arr_bytes_r2

   pure function arr_bytes_r3(a) result(b)
      real(wp), allocatable, intent(in) :: a(:, :, :)
      integer(int64) :: b
      b = 0_int64
      if (allocated(a)) b = size(a, kind=int64)*int(storage_size(1.0_wp)/8, int64)
   end function arr_bytes_r3

   pure function arr_bytes_r4(a) result(b)
      real(wp), allocatable, intent(in) :: a(:, :, :, :)
      integer(int64) :: b
      b = 0_int64
      if (allocated(a)) b = size(a, kind=int64)*int(storage_size(1.0_wp)/8, int64)
   end function arr_bytes_r4

   pure function arr_bytes_i1(a) result(b)
      integer, allocatable, intent(in) :: a(:)
      integer(int64) :: b
      b = 0_int64
      if (allocated(a)) b = size(a, kind=int64)*int(storage_size(0)/8, int64)
   end function arr_bytes_i1

   pure function arr_bytes_i2(a) result(b)
      integer, allocatable, intent(in) :: a(:, :)
      integer(int64) :: b
      b = 0_int64
      if (allocated(a)) b = size(a, kind=int64)*int(storage_size(0)/8, int64)
   end function arr_bytes_i2

   pure function arr_bytes_i3(a) result(b)
      integer, allocatable, intent(in) :: a(:, :, :)
      integer(int64) :: b
      b = 0_int64
      if (allocated(a)) b = size(a, kind=int64)*int(storage_size(0)/8, int64)
   end function arr_bytes_i3

   pure function arr_bytes_l1(a) result(b)
      logical, allocatable, intent(in) :: a(:)
      integer(int64) :: b
      b = 0_int64
      if (allocated(a)) b = size(a, kind=int64)*int(storage_size(.true.)/8, int64)
   end function arr_bytes_l1

   pure function arr_bytes_l2(a) result(b)
      logical, allocatable, intent(in) :: a(:, :)
      integer(int64) :: b
      b = 0_int64
      if (allocated(a)) b = size(a, kind=int64)*int(storage_size(.true.)/8, int64)
   end function arr_bytes_l2

   pure function arr_bytes_l3(a) result(b)
      logical, allocatable, intent(in) :: a(:, :, :)
      integer(int64) :: b
      b = 0_int64
      if (allocated(a)) b = size(a, kind=int64)*int(storage_size(.true.)/8, int64)
   end function arr_bytes_l3

   subroutine mem_set_counted_budget(total_bytes)
      !! Latch the COUNTED state-array footprint for the `enter_data`
      !! reconciliation, WITHOUT logging.  Must run before
      !! `mem_log_device_actuals` (the reconciliation reads this), even
      !! when the human-readable budget is printed later (e.g. folded into
      !! the driver's post-setup banner).  `total_bytes` is the sum of
      !! every state `bytes()` term (0 for gated-off slots).
      integer(int64), intent(in) :: total_bytes
      counted_state_bytes = total_bytes
   end subroutine mem_set_counted_budget

   subroutine mem_log_computed_budget(label, total_bytes)
      !! Log the COUNTED state-array footprint.  Unlike the RSS figure
      !! (host-only I/O buffers + library pages) and the device query
      !! (needs an OpenACC runtime), this is an exact sum of the state's
      !! allocatable arrays — available on every toolchain, so a CPU build
      !! gets a real footprint line too.  Latching for the reconciliation
      !! is `mem_set_counted_budget` (separate, so the print can move to
      !! the banner while the latch stays before `enter_data`).
      character(len=*), intent(in) :: label
      integer(int64), intent(in) :: total_bytes

      ! Two-space indent + four-space breakdown so the block nests cleanly
      ! under the driver's post-setup banner (top-level items at two spaces,
      ! sub-items at four).
      call logger%info("  Memory:   state arrays counted ~ "// &
                       mem_format_bytes(total_bytes)//" ("//trim(label)// &
                       ", exact allocatable footprint)")
   end subroutine mem_log_computed_budget

   subroutine mem_log_computed_line(label, bytes)
      !! Log one indented breakdown line for a component of the counted
      !! footprint (barotropic / layers / closures / …).  Skipped when the
      !! component is empty so gated-off slots don't clutter the report.
      character(len=*), intent(in) :: label
      integer(int64), intent(in) :: bytes

      if (bytes <= 0_int64) return
      call logger%info("    "//trim(label)//": "//mem_format_bytes(bytes))
   end subroutine mem_log_computed_line

   subroutine mem_log_state_budget(label, host_growth_bytes)
      !! Log the pre-mapping memory budget.  Always logs the host-side
      !! allocation growth (ADVISORY — includes host-only I/O buffers +
      !! library pages, NOT a device estimate); when the device queries
      !! succeed it also logs the device number + free/total and latches
      !! the free-memory snapshot that `mem_log_device_actuals` turns
      !! into the measured mapping delta.  Warns when the host growth
      !! (a conservative upper bound on what `enter_data` can map)
      !! exceeds 95% of free device memory — the regime where
      !! `enter_data` is likely to abort.
      character(len=*), intent(in) :: label
      integer(int64), intent(in) :: host_growth_bytes

      integer(int64) :: free_b, total_b
      integer :: devnum

      call logger%info("Memory budget ("//trim(label)//"): host allocations across setup ~ "// &
                       mem_format_bytes(host_growth_bytes)// &
                       " (RSS growth; advisory, not a device estimate)")

      free_b = mem_device_free_bytes()
      total_b = mem_device_total_bytes()
      if (free_b < 0_int64 .or. total_b < 0_int64) return  ! CPU build.

      ! Baseline for the measured mapping delta.  The property query
      ! above already forced CUDA context creation, so the delta across
      ! enter_data excludes the ~300-500 MB context cost.
      free_before_map = free_b

      devnum = device_num()
      call logger%info("  device "//to_string(devnum)//" free before mapping: "// &
                       mem_format_bytes(free_b)//" of "//mem_format_bytes(total_b)// &
                       " (host-growth upper bound: "// &
                       mem_format_bytes(host_growth_bytes)//")")

      if (real(host_growth_bytes, wp) > 0.95_wp*real(free_b, wp)) then
         call logger%warning("device memory likely insufficient — mapping <= "// &
                             mem_format_bytes(host_growth_bytes)//" with only "// &
                             mem_format_bytes(free_b)//" free")
      end if
   end subroutine mem_log_state_budget

   subroutine mem_log_device_actuals(label)
      !! Log actual device memory after `enter_data`: the device-wide
      !! usage (`used = total - free` — includes the CUDA context and
      !! any other process on a shared device) plus the MEASURED state
      !! mapping delta (free-before-mapping minus free-now; accurate to
      !! NVHPC pool granularity, ~a few MB).  Latches the post-mapping
      !! usage as the baseline for `mem_log_device_growth`.  Silent
      !! no-op when device queries are unavailable (CPU build).
      character(len=*), intent(in) :: label
      integer(int64) :: free_b, total_b, used_b

      free_b = mem_device_free_bytes()
      total_b = mem_device_total_bytes()
      if (free_b < 0_int64 .or. total_b < 0_int64) return

      used_b = total_b - free_b
      call logger%info("device memory in use: "//mem_format_bytes(used_b)// &
                       " of "//mem_format_bytes(total_b)//" ("//trim(label)// &
                       "; includes CUDA context ~300-500 MB + other processes)")
      if (free_before_map >= 0_int64) then
         block
            integer(int64) :: mapped_b, tol_b
            mapped_b = max(free_before_map - free_b, 0_int64)
            call logger%info("  state mapped: ~"//mem_format_bytes(mapped_b)// &
                             " (device-measured; shared-device activity can perturb)")

            ! Reconcile the MEASURED mapping against the COUNTED footprint.
            ! Overshoot past the tolerance means the device mapped more than
            ! every bytes() term accounts for — a state array is almost
            ! certainly missing from the count (the drift the counted total
            ! would otherwise hide).  Under-count from lazy/host-only arrays
            ! is expected and stays silent.
            if (counted_state_bytes > 0_int64) then
               tol_b = counted_state_bytes + &
                       int(COUNT_RECONCILE_FRAC*real(counted_state_bytes, wp), int64) + &
                       COUNT_RECONCILE_FLOOR
               if (mapped_b > tol_b) then
                  call logger%warning("device mapped ~"//mem_format_bytes(mapped_b)// &
                                      " but only ~"//mem_format_bytes(counted_state_bytes)// &
                                      " counted — a state array is likely missing from a "// &
                                      "bytes() term (memory-accounting drift)")
               end if
            end if
         end block
      end if

      ! Baselines for post-setup growth reporting.
      used_after_map = used_b
      used_steady = -1_int64
      growth_warned = .false.
   end subroutine mem_log_device_actuals

   subroutine mem_log_device_growth(label, quiet)
      !! Report device-memory growth since the post-`enter_data`
      !! baseline.  Two allocations land AFTER the setup-time report and
      !! would otherwise stay invisible: the lazy first-call kernel
      !! workspaces (`*_workspace_ensure`, ~33% of the coastal device
      !! footprint) and the one-off CUDA runtime/kernel reservation at
      !! first launch (GBs — module-header offset 3).  Call this in the
      !! end-of-run summary (prints usage + the "grew ~X after setup"
      !! attribution, so even a 10-step run shows the true footprint)
      !! and at the periodic status with `quiet=.true.` (no info lines).
      !! The one-shot leak WARNING measures PROGRESSIVE growth — past
      !! `GROWTH_WARN_FRAC` of the steady baseline latched at the first
      !! post-setup sample — so the expected one-offs never trip it.
      !! Silent no-op on CPU builds or when the baseline was never
      !! latched.
      character(len=*), intent(in) :: label
      logical, intent(in), optional :: quiet

      integer(int64) :: used_b, growth_b
      logical :: quiet_mode

      quiet_mode = .false.
      if (present(quiet)) quiet_mode = quiet

      used_b = mem_device_used_bytes()
      if (used_b < 0_int64 .or. used_after_map < 0_int64) return

      growth_b = used_b - used_after_map

      if (.not. quiet_mode) then
         call logger%info("device memory in use: "//mem_format_bytes(used_b)// &
                          " ("//trim(label)//")")
         if (growth_b > 0_int64) then
            call logger%info("  grew ~"//mem_format_bytes(growth_b)// &
                             " after setup (lazy kernel workspaces + one-off "// &
                             "CUDA kernel/runtime reservation)")
         end if
      end if

      if (used_steady < 0_int64) then
         ! First post-setup sample: the one-off reservations have landed;
         ! this is the steady baseline for progressive-growth (leak)
         ! detection.
         used_steady = used_b
      else if (.not. growth_warned .and. &
               real(used_b - used_steady, wp) > &
               GROWTH_WARN_FRAC*real(used_steady, wp)) then
         call logger%warning("device memory grew ~"// &
                             mem_format_bytes(used_b - used_steady)// &
                             " past the steady post-setup baseline ("// &
                             mem_format_bytes(used_steady)// &
                             ") — possible device-memory leak ("//trim(label)//")")
         growth_warned = .true.
      end if
   end subroutine mem_log_device_growth

   function device_num() result(num)
      !! Active OpenACC device number, or `-1` on CPU builds.  Kept
      !! private — only the budget logger needs it.
      integer :: num
#ifdef RAKALI_HAS_OPENACC_RUNTIME
      num = acc_get_device_num(acc_get_device_type())
#else
      num = -1
#endif
   end function device_num

end module rki_mem_report
