!! Named profiling regions with wall-clock timers and optional NVTX ranges.
!! -DUSE_NVTX enables NVTX (else no-op); -DDISABLE_PROFILER disables everything.
module rki_profiler
   use, intrinsic :: iso_fortran_env, only: dp => real64, int64
   use pic_logger, only: logger => global_logger
#ifdef USE_NVTX
   use nvtx, only: nvtxStartRange, nvtxEndRange
#endif
   implicit none
   private

   integer, parameter :: MAX_REGIONS = 256
   integer, parameter :: MAX_NAME_LEN = 64

   type :: tracked_region
      character(len=MAX_NAME_LEN) :: name = ""
      real(dp) :: start_time = 0.0_dp
      real(dp) :: total_time = 0.0_dp
      integer :: call_count = 0
      logical :: active = .false.
      logical :: nvtx_only = .false.
         !! If true, only shows in NVTX timeline, not in text report
   end type tracked_region

   type :: profiler_state
      logical :: initialized = .false.
      logical :: enabled = .true.
      integer :: num_regions = 0
      type(tracked_region) :: regions(MAX_REGIONS)
   end type profiler_state

   type(profiler_state), save :: state

   public :: profiler_init, profiler_end
   public :: profiler_start, profiler_stop
   public :: profiler_enable, profiler_disable
   public :: profiler_report, profiler_reset
   public :: profiler_get_time

contains

   function get_wall_time() result(t)
      real(dp) :: t
      integer(int64) :: count, count_rate
      call system_clock(count, count_rate)
      t = real(count, dp)/real(count_rate, dp)
   end function get_wall_time

   subroutine nvtx_range_push(name)
      character(len=*), intent(in) :: name
#ifdef USE_NVTX
      call nvtxStartRange(name)
#endif
   end subroutine nvtx_range_push

   subroutine nvtx_range_pop()
#ifdef USE_NVTX
      call nvtxEndRange()
#endif
   end subroutine nvtx_range_pop

   subroutine profiler_init(enabled)
      !! Initialise the profiler
      logical, intent(in), optional :: enabled

      state%initialized = .true.
      state%enabled = .true.
      if (present(enabled)) state%enabled = enabled
      state%num_regions = 0
   end subroutine profiler_init

   subroutine profiler_end()
      !! Finalise the profiler
      state%initialized = .false.
      state%num_regions = 0
   end subroutine profiler_end

   subroutine profiler_enable()
      !! Public for the unit-test suite only.
      state%enabled = .true.
   end subroutine profiler_enable

   subroutine profiler_disable()
      !! Public for the unit-test suite only.
      state%enabled = .false.
   end subroutine profiler_disable

   function find_or_create_region(name) result(idx)
      character(len=*), intent(in) :: name
      integer :: idx
      integer :: i

      do i = 1, state%num_regions
         if (trim(state%regions(i)%name) == trim(name)) then
            idx = i
            return
         end if
      end do

      if (state%num_regions < MAX_REGIONS) then
         state%num_regions = state%num_regions + 1
         idx = state%num_regions
         state%regions(idx)%name = trim(name)
         state%regions(idx)%total_time = 0.0_dp
         state%regions(idx)%call_count = 0
         state%regions(idx)%active = .false.
      else
         idx = -1
      end if
   end function find_or_create_region

   subroutine profiler_start(name, nvtx_only)
      !! Start a named profiling region
      !! If nvtx_only is true, the region only appears in NVTX timeline, not in text report
      character(len=*), intent(in) :: name
      logical, intent(in), optional :: nvtx_only
      integer :: idx

#ifdef DISABLE_PROFILER
      return
#endif

      if (.not. state%enabled) return

      idx = find_or_create_region(name)
      if (idx < 0) return

      if (state%regions(idx)%active) return

      if (present(nvtx_only) .and. state%regions(idx)%call_count == 0) then
         state%regions(idx)%nvtx_only = nvtx_only
      end if

      state%regions(idx)%active = .true.
      call nvtx_range_push(name)
      state%regions(idx)%start_time = get_wall_time()
   end subroutine profiler_start

   subroutine profiler_stop(name)
      !! Stop a named profiling region
      character(len=*), intent(in) :: name
      integer :: idx
      real(dp) :: elapsed

#ifdef DISABLE_PROFILER
      return
#endif

      if (.not. state%enabled) return

      do idx = 1, state%num_regions
         if (trim(state%regions(idx)%name) == trim(name)) exit
      end do

      if (idx > state%num_regions) return
      if (.not. state%regions(idx)%active) return

      elapsed = get_wall_time() - state%regions(idx)%start_time
      call nvtx_range_pop()

      state%regions(idx)%total_time = state%regions(idx)%total_time + elapsed
      state%regions(idx)%call_count = state%regions(idx)%call_count + 1
      state%regions(idx)%active = .false.
   end subroutine profiler_stop

   function profiler_get_time(name) result(t)
      !! Accumulated time (s) for a named region. Public for the unit-test suite only.
      character(len=*), intent(in) :: name
      real(dp) :: t
      integer :: idx

      t = 0.0_dp
      do idx = 1, state%num_regions
         if (trim(state%regions(idx)%name) == trim(name)) then
            t = state%regions(idx)%total_time
            return
         end if
      end do
   end function profiler_get_time

   subroutine profiler_reset()
      !! Reset all timing data
      integer :: i

      do i = 1, state%num_regions
         state%regions(i)%total_time = 0.0_dp
         state%regions(i)%call_count = 0
         state%regions(i)%active = .false.
      end do
   end subroutine profiler_reset

   subroutine profiler_report(title, root_region)
      !! Print profiling report
      !! If root_region is specified, percentages are relative to that region
      character(len=*), intent(in), optional :: title
      character(len=*), intent(in), optional :: root_region
      integer :: i, j, root_idx, n_print, tmp_idx
      integer :: sorted_idx(MAX_REGIONS)
      real(dp) :: total_time, pct
      character(len=256) :: line

      if (state%num_regions == 0) return

      root_idx = -1
      if (present(root_region)) then
         do i = 1, state%num_regions
            if (trim(state%regions(i)%name) == trim(root_region)) then
               root_idx = i
               exit
            end if
         end do
      end if

      if (root_idx > 0) then
         total_time = state%regions(root_idx)%total_time
      else
         total_time = 0.0_dp
         do i = 1, state%num_regions
            total_time = total_time + state%regions(i)%total_time
         end do
      end if

      n_print = 0
      do i = 1, state%num_regions
         if (i /= root_idx .and. .not. state%regions(i)%nvtx_only) then
            n_print = n_print + 1
            sorted_idx(n_print) = i
         end if
      end do

      ! Sort by time descending (insertion sort)
      do i = 2, n_print
         tmp_idx = sorted_idx(i)
         j = i - 1
         ! Guard the sorted_idx(j) read with a SEPARATE bounds check: Fortran
         ! does not short-circuit `.and.`, so a combined
         ! `do while (j >= 1 .and. ...sorted_idx(j)...)` reads sorted_idx(0)
         ! when j reaches 0 (OOB — aborts under -fcheck=bounds).
         do
            if (j < 1) exit
            if (state%regions(sorted_idx(j))%total_time >= &
                state%regions(tmp_idx)%total_time) exit
            sorted_idx(j + 1) = sorted_idx(j)
            j = j - 1
         end do
         sorted_idx(j + 1) = tmp_idx
      end do

      call logger%info("")
      call logger%info("============================================================")
      if (present(title)) then
         call logger%info("Profiler Report: "//trim(title))
      else
         call logger%info("Profiler Report")
      end if
      call logger%info("============================================================")
      call logger%info("  Region                          Time (s)    Calls    %    ")
      call logger%info("------------------------------------------------------------")

      do i = 1, n_print
         j = sorted_idx(i)
         if (total_time > 0.0_dp) then
            pct = 100.0_dp*state%regions(j)%total_time/total_time
         else
            pct = 0.0_dp
         end if
         write (line, "(A,A32,F12.6,I8,F8.1)") "  ", &
            state%regions(j)%name, &
            state%regions(j)%total_time, &
            state%regions(j)%call_count, &
            pct
         call logger%info(trim(line))
      end do

      call logger%info("------------------------------------------------------------")
      write (line, "(A,F12.6)") "  Total:                        ", total_time
      call logger%info(trim(line))
      call logger%info("============================================================")

#ifdef USE_NVTX
      call logger%info("  (NVTX enabled - use Nsight Systems for GPU timeline)")
#else
      call logger%info("  (NVTX disabled - compile with -DRAKALI_ENABLE_NVTX to enable)")
#endif
   end subroutine profiler_report

end module rki_profiler
