!! Unit tests for the C API module (rki_api) — handle-based.
!!
!! Covers the one-shot and two-phase lifecycles in a single executable
!! (formerly split across test_rki_api + test_rki_api_twophase due to the
!! singleton GPU-mapping constraint), plus the "two handles coexist"
!! test that the singleton design literally couldn't express.
module test_rki_api
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, &
                                                                             c_f_pointer, c_int, c_null_ptr, &
                                                                             c_ptr
   use rki_constants, only: wp
   implicit none
   private

   public :: collect_api_tests
   public :: collect_api_oneshot_test
   public :: collect_api_twophase_test
   public :: collect_api_two_handles_test
   public :: collect_api_bad_handle_test
   public :: collect_api_string_test

   interface
      function rki_create(filename, flen, c_handle) result(ierr) &
         bind(c, name='rki_create')
         import :: c_char, c_int, c_ptr
         character(kind=c_char), intent(in) :: filename(*)
         integer(c_int), intent(in), value :: flen
         type(c_ptr), intent(out) :: c_handle
         integer(c_int) :: ierr
      end function

      function rki_create_from_namelist(filename, flen, c_handle) result(ierr) &
         bind(c, name='rki_create_from_namelist')
         import :: c_char, c_int, c_ptr
         character(kind=c_char), intent(in) :: filename(*)
         integer(c_int), intent(in), value :: flen
         type(c_ptr), intent(out) :: c_handle
         integer(c_int) :: ierr
      end function

      function rki_create_from_namelist_string(buf, blen, c_handle) result(ierr) &
         bind(c, name='rki_create_from_namelist_string')
         import :: c_char, c_int, c_ptr
         character(kind=c_char), intent(in) :: buf(*)
         integer(c_int), intent(in), value :: blen
         type(c_ptr), intent(out) :: c_handle
         integer(c_int) :: ierr
      end function

      function rki_finalize(c_handle) result(ierr) bind(c, name='rki_finalize')
         import :: c_int, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int) :: ierr
      end function

      subroutine rki_destroy(c_handle) bind(c, name='rki_destroy')
         import :: c_ptr
         type(c_ptr), intent(inout) :: c_handle
      end subroutine

      function rki_step(c_handle, n) result(ierr) bind(c, name='rki_step')
         import :: c_int, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int), intent(in), value :: n
         integer(c_int) :: ierr
      end function

      function rki_advance(c_handle, target_time) result(ierr) &
         bind(c, name='rki_advance')
         import :: c_int, c_double, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         real(c_double), intent(in), value :: target_time
         integer(c_int) :: ierr
      end function

      function rki_get_time(c_handle) result(t) bind(c, name='rki_get_time')
         import :: c_double, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         real(c_double) :: t
      end function

      function rki_get_dt(c_handle) result(dt) bind(c, name='rki_get_dt')
         import :: c_double, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         real(c_double) :: dt
      end function

      function rki_get_step_count(c_handle) result(n) &
         bind(c, name='rki_get_step_count')
         import :: c_int, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int) :: n
      end function

      subroutine rki_get_grid_info(c_handle, nx, ny, nghost, dx, dy) &
         bind(c, name='rki_get_grid_info')
         import :: c_int, c_double, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int), intent(out) :: nx, ny, nghost
         real(c_double), intent(out) :: dx, dy
      end subroutine

      function rki_get_h_ptr(c_handle, ptr, nx_total, ny_total) result(ierr) &
         bind(c, name='rki_get_h_ptr')
         import :: c_ptr, c_int
         type(c_ptr), intent(in), value :: c_handle
         type(c_ptr), intent(out) :: ptr
         integer(c_int), intent(out) :: nx_total, ny_total
         integer(c_int) :: ierr
      end function

      function rki_get_hu_ptr(c_handle, ptr, nx_total, ny_total) result(ierr) &
         bind(c, name='rki_get_hu_ptr')
         import :: c_ptr, c_int
         type(c_ptr), intent(in), value :: c_handle
         type(c_ptr), intent(out) :: ptr
         integer(c_int), intent(out) :: nx_total, ny_total
         integer(c_int) :: ierr
      end function

      function rki_get_hv_ptr(c_handle, ptr, nx_total, ny_total) result(ierr) &
         bind(c, name='rki_get_hv_ptr')
         import :: c_ptr, c_int
         type(c_ptr), intent(in), value :: c_handle
         type(c_ptr), intent(out) :: ptr
         integer(c_int), intent(out) :: nx_total, ny_total
         integer(c_int) :: ierr
      end function

      function rki_get_b_ptr(c_handle, ptr, nx_total, ny_total) result(ierr) &
         bind(c, name='rki_get_b_ptr')
         import :: c_ptr, c_int
         type(c_ptr), intent(in), value :: c_handle
         type(c_ptr), intent(out) :: ptr
         integer(c_int), intent(out) :: nx_total, ny_total
         integer(c_int) :: ierr
      end function

      function rki_set_bathymetry(c_handle, b_data, nx_p, ny_p) result(ierr) &
         bind(c, name='rki_set_bathymetry')
         import :: c_int, c_double, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int), intent(in), value :: nx_p, ny_p
         real(c_double), intent(in) :: b_data(nx_p, ny_p)
         integer(c_int) :: ierr
      end function

      function rki_set_initial_h(c_handle, h_data, nx_p, ny_p) result(ierr) &
         bind(c, name='rki_set_initial_h')
         import :: c_int, c_double, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int), intent(in), value :: nx_p, ny_p
         real(c_double), intent(in) :: h_data(nx_p, ny_p)
         integer(c_int) :: ierr
      end function
   end interface

contains

   ! Bundled collector — used by the fpm test driver (tests/tester.F90)
   ! which builds a single executable that runs every test suite back-
   ! to-back in one process. Under that mode the four scenarios share a
   ! process and the OpenMP-target present table on Aurora ifx leaks
   ! between them (handle teardown doesn't fully release nested-
   ! allocatable maps). CMake/CTest avoids the pollution by using the
   ! four single-test collectors below, one binary per scenario.
   subroutine collect_api_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("oneshot_lifecycle", test_oneshot_lifecycle), &
                  new_unittest("twophase_lifecycle", test_twophase_lifecycle), &
                  new_unittest("two_handles_coexist", test_two_handles_coexist), &
                  new_unittest("bad_handle_rejected", test_bad_handle_rejected) &
                  ]
   end subroutine collect_api_tests

   ! Per-scenario collectors — one ctest binary per scenario so each
   ! starts with a virgin present table. Mirrors the canonical pattern
   ! from main-openmp commit 6b884ac.
   subroutine collect_api_oneshot_test(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("oneshot_lifecycle", test_oneshot_lifecycle)]
   end subroutine collect_api_oneshot_test

   subroutine collect_api_twophase_test(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("twophase_lifecycle", test_twophase_lifecycle)]
   end subroutine collect_api_twophase_test

   subroutine collect_api_two_handles_test(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("two_handles_coexist", test_two_handles_coexist)]
   end subroutine collect_api_two_handles_test

   subroutine collect_api_bad_handle_test(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("bad_handle_rejected", test_bad_handle_rejected)]
   end subroutine collect_api_bad_handle_test

   subroutine collect_api_string_test(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("create_from_string_matches_file", &
                                test_create_from_string_matches_file)]
   end subroutine collect_api_string_test

   subroutine write_default_namelist(filename, nx_p, ny_p)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx_p, ny_p
      integer :: io_unit

      open (newunit=io_unit, file=trim(filename), status="replace", action="write")
      write (io_unit, '(a)') "&grid_nml"
      write (io_unit, '("  nx = ", i0, ", ny = ", i0)') nx_p, ny_p
      write (io_unit, '(a)') "  dx = 1.0, dy = 1.0"
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&time_nml"
      write (io_unit, '(a)') "  t_end = 10.0, cfl = 0.45"
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') '&initial_condition_nml'
      write (io_unit, '(a)') '  init_condition = "dam_break"'
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') '&logging_nml'
      write (io_unit, '(a)') '  log_level = "warning"'
      write (io_unit, '(a)') "/"
      close (io_unit)
   end subroutine write_default_namelist

   subroutine write_quiet_namelist(filename, nx_p, ny_p)
      !! Same as default but without an init_condition (caller will set h).
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx_p, ny_p
      integer :: io_unit

      open (newunit=io_unit, file=trim(filename), status="replace", action="write")
      write (io_unit, '(a)') "&grid_nml"
      write (io_unit, '("  nx = ", i0, ", ny = ", i0)') nx_p, ny_p
      write (io_unit, '(a)') "  dx = 1.0, dy = 1.0"
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&time_nml"
      write (io_unit, '(a)') "  t_end = 10.0, cfl = 0.45"
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') '&logging_nml'
      write (io_unit, '(a)') '  log_level = "warning"'
      write (io_unit, '(a)') "/"
      close (io_unit)
   end subroutine write_quiet_namelist

   subroutine test_oneshot_lifecycle(error)
      !! create_from_namelist -> step -> advance -> read fields -> destroy.
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: rc, nx, ny, nghost, step_n, nx_t, ny_t
      real(c_double) :: dx, dy, t
      type(c_ptr) :: c_h, h_ptr, b_ptr
      real(c_double), pointer :: h(:, :)
      character(len=64) :: f

      f = "test_api_oneshot.nml"
      call write_default_namelist(f, 50, 10)

      rc = rki_create_from_namelist(trim(f), len_trim(f, kind=c_int), c_h)
      call check(error, rc == 0, "create_from_namelist should succeed")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if
      call check(error, c_associated(c_h), "handle is non-null")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      call rki_get_grid_info(c_h, nx, ny, nghost, dx, dy)
      call check(error, nx == 50, "nx should be 50")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if
      call check(error, ny == 10, "ny should be 10")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if
      call check(error, nghost == 3, "nghost should be 3 (default)")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      rc = rki_step(c_h, 5_c_int)
      call check(error, rc == 0, "step(5) should succeed")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      t = rki_get_time(c_h)
      call check(error, t > 0.0_c_double, "time > 0 after stepping")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      step_n = rki_get_step_count(c_h)
      call check(error, step_n == 5, "step_count should be 5")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      call check(error, rki_get_dt(c_h) > 0.0_c_double, "dt > 0")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      rc = rki_advance(c_h, 1.0_c_double)
      call check(error, rc == 0, "advance should succeed")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if
      call check(error, rki_get_time(c_h) >= 1.0_c_double - 1.0e-10_c_double, &
                 "time should reach ~1.0")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      rc = rki_get_h_ptr(c_h, h_ptr, nx_t, ny_t)
      call check(error, rc == 0, "get_h_ptr should succeed")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if
      call check(error, c_associated(h_ptr), "h_ptr not null")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      rc = rki_get_b_ptr(c_h, b_ptr, nx_t, ny_t)
      call check(error, rc == 0, "get_b_ptr should succeed")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if
      call check(error, c_associated(b_ptr), "b_ptr not null")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      ! -- Dam break is physically correct (well-balanced, mass-conserving)
      block
         integer :: i, j, nx_int, ny_int
         real(wp) :: h_max_interior, h_min_interior, total_mass
         real(wp) :: expected_mass
         integer :: count_low, count_high

         call c_f_pointer(h_ptr, h, [nx + 2*nghost, ny + 2*nghost])
         nx_int = nx
         ny_int = ny
         h_max_interior = h(nghost + 1, nghost + 1)
         h_min_interior = h(nghost + 1, nghost + 1)
         total_mass = 0.0_wp
         count_low = 0
         count_high = 0
         do j = nghost + 1, nghost + ny_int
            do i = nghost + 1, nghost + nx_int
               h_max_interior = max(h_max_interior, h(i, j))
               h_min_interior = min(h_min_interior, h(i, j))
               total_mass = total_mass + h(i, j)
               if (h(i, j) < 1.3_wp) count_low = count_low + 1
               if (h(i, j) > 1.7_wp) count_high = count_high + 1
            end do
         end do

         call check(error, abs(h_max_interior - 2.0_wp) < 1.0e-3_wp, &
                    "dam break: upstream h_max must still equal H_LEFT = 2.0")
         if (allocated(error)) then; call rki_destroy(c_h); return; end if
         call check(error, h_min_interior >= 1.0_wp - 1.0e-6_wp, &
                    "dam break: h_min must be >= H_RIGHT = 1.0 (no undershoot)")
         if (allocated(error)) then; call rki_destroy(c_h); return; end if
         expected_mass = (25.0_wp*2.0_wp + 25.0_wp*1.0_wp)*real(ny_int, wp)
         call check(error, abs(total_mass - expected_mass)/expected_mass < 1.0e-6_wp, &
                    "dam break: mass conserved to 1e-6 (closed basin)")
         if (allocated(error)) then; call rki_destroy(c_h); return; end if
         call check(error, count_high > 0, &
                    "dam break: some interior cells still at h > 1.7 (high side)")
         if (allocated(error)) then; call rki_destroy(c_h); return; end if
         call check(error, count_low > 0, &
                    "dam break: some interior cells at h < 1.3 (low side)")
         if (allocated(error)) then; call rki_destroy(c_h); return; end if
      end block

      call rki_destroy(c_h)
      call check(error,.not. c_associated(c_h), "destroy nulls handle")
   end subroutine test_oneshot_lifecycle

   subroutine test_twophase_lifecycle(error)
      !! create -> step-blocked -> wrong-size-bathy rejected -> set fields ->
      !! finalize -> double-finalize-blocked -> step still water -> destroy.
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: rc, nx_t, ny_t
      integer, parameter :: NX = 50, NY = 10
      real(c_double) :: bathy(NX, NY), h_init(NX, NY)
      real(c_double) :: bathy_wrong(60, 10)
      type(c_ptr) :: c_h, h_ptr
      character(len=64) :: f

      f = "test_api_twophase.nml"
      call write_quiet_namelist(f, NX, NY)

      rc = rki_create(trim(f), len_trim(f, kind=c_int), c_h)
      call check(error, rc == 0, "create should succeed")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      rc = rki_step(c_h, 1_c_int)
      call check(error, rc == -1, "step before finalize should return -1")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      bathy_wrong = 0.0_c_double
      rc = rki_set_bathymetry(c_h, bathy_wrong, 60_c_int, 10_c_int)
      call check(error, rc == -2, "wrong-size bathymetry should return -2")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      bathy = 0.0_c_double
      h_init = 5.0_c_double
      rc = rki_set_bathymetry(c_h, bathy, int(NX, c_int), int(NY, c_int))
      call check(error, rc == 0, "set_bathymetry should succeed")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      rc = rki_set_initial_h(c_h, h_init, int(NX, c_int), int(NY, c_int))
      call check(error, rc == 0, "set_initial_h should succeed")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      rc = rki_finalize(c_h)
      call check(error, rc == 0, "finalize should succeed")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      rc = rki_finalize(c_h)
      call check(error, rc == -2, "second finalize should return -2")
      if (allocated(error)) then; call rki_destroy(c_h); return; end if

      ! -- Still water remains exactly still
      block
         real(c_double), pointer :: h(:, :), hu(:, :), hv(:, :)
         type(c_ptr) :: hu_ptr, hv_ptr
         integer(c_int) :: nx_total, ny_total, nghost
         real(c_double) :: dx_c, dy_c
         real(wp) :: err_h, err_hu, err_hv
         integer :: i, j, i0, i1, j0, j1

         rc = rki_step(c_h, 3_c_int)
         call check(error, rc == 0, "step should succeed")
         if (allocated(error)) then; call rki_destroy(c_h); return; end if

         call rki_get_grid_info(c_h, nx_total, ny_total, nghost, dx_c, dy_c)

         rc = rki_get_h_ptr(c_h, h_ptr, nx_t, ny_t)
         call check(error, rc == 0, "get_h_ptr should succeed")
         if (allocated(error)) then; call rki_destroy(c_h); return; end if
         rc = rki_get_hu_ptr(c_h, hu_ptr, nx_t, ny_t)
         call check(error, rc == 0, "get_hu_ptr should succeed")
         if (allocated(error)) then; call rki_destroy(c_h); return; end if
         rc = rki_get_hv_ptr(c_h, hv_ptr, nx_t, ny_t)
         call check(error, rc == 0, "get_hv_ptr should succeed")
         if (allocated(error)) then; call rki_destroy(c_h); return; end if

         call c_f_pointer(h_ptr, h, [nx_t, ny_t])
         call c_f_pointer(hu_ptr, hu, [nx_t, ny_t])
         call c_f_pointer(hv_ptr, hv, [nx_t, ny_t])

         i0 = nghost + 1
         i1 = nx_t - nghost
         j0 = nghost + 1
         j1 = ny_t - nghost
         err_h = 0.0_wp
         err_hu = 0.0_wp
         err_hv = 0.0_wp
         do j = j0, j1
            do i = i0, i1
               err_h = max(err_h, abs(h(i, j) - 5.0_wp))
               err_hu = max(err_hu, abs(hu(i, j)))
               err_hv = max(err_hv, abs(hv(i, j)))
            end do
         end do

         call check(error, err_h < 1.0e-12_wp, &
                    "still water: max|h - 5.0| must be machine precision")
         if (allocated(error)) then; call rki_destroy(c_h); return; end if
         call check(error, err_hu < 1.0e-12_wp, &
                    "still water: max|hu| must be machine precision (no motion)")
         if (allocated(error)) then; call rki_destroy(c_h); return; end if
         call check(error, err_hv < 1.0e-12_wp, &
                    "still water: max|hv| must be machine precision (no motion)")
      end block

      call rki_destroy(c_h)
   end subroutine test_twophase_lifecycle

   subroutine test_two_handles_coexist(error)
      !! Two structured simulations in one process — the headline test
      !! that the singleton design literally couldn't express.
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: rc, n_a, n_b
      type(c_ptr) :: c_a, c_b
      character(len=64) :: f_a, f_b

      f_a = "test_api_handle_a.nml"
      f_b = "test_api_handle_b.nml"
      call write_default_namelist(f_a, 30, 8)
      call write_default_namelist(f_b, 50, 10)

      rc = rki_create_from_namelist(trim(f_a), len_trim(f_a, kind=c_int), c_a)
      call check(error, rc == 0, "create handle A succeeds")
      if (allocated(error)) then; call rki_destroy(c_a); return; end if

      rc = rki_create_from_namelist(trim(f_b), len_trim(f_b, kind=c_int), c_b)
      call check(error, rc == 0, "create handle B succeeds (coexists with A)")
      if (allocated(error)) then; call rki_destroy(c_a); call rki_destroy(c_b); return; end if

      call check(error, c_associated(c_a) .and. c_associated(c_b), &
                 "both handles non-null")
      if (allocated(error)) then; call rki_destroy(c_a); call rki_destroy(c_b); return; end if

      !! Step each independently, with different step counts.
      rc = rki_step(c_a, 3_c_int)
      call check(error, rc == 0, "step A by 3 succeeds")
      if (allocated(error)) then; call rki_destroy(c_a); call rki_destroy(c_b); return; end if
      rc = rki_step(c_b, 7_c_int)
      call check(error, rc == 0, "step B by 7 succeeds")
      if (allocated(error)) then; call rki_destroy(c_a); call rki_destroy(c_b); return; end if

      n_a = rki_get_step_count(c_a)
      n_b = rki_get_step_count(c_b)
      call check(error, n_a == 3, "A advanced exactly 3 steps")
      if (allocated(error)) then; call rki_destroy(c_a); call rki_destroy(c_b); return; end if
      call check(error, n_b == 7, "B advanced exactly 7 steps")
      if (allocated(error)) then; call rki_destroy(c_a); call rki_destroy(c_b); return; end if

      !! Destroy A first, B must still work.
      call rki_destroy(c_a)
      call check(error,.not. c_associated(c_a), "A handle nulled after destroy")
      if (allocated(error)) then; call rki_destroy(c_b); return; end if

      rc = rki_step(c_b, 1_c_int)
      call check(error, rc == 0, "B still steppable after A destroyed")
      if (allocated(error)) then; call rki_destroy(c_b); return; end if

      call rki_destroy(c_b)
   end subroutine test_two_handles_coexist

   subroutine test_bad_handle_rejected(error)
      !! Null and stale handles are rejected with the bad-handle code.
      type(error_type), allocatable, intent(out) :: error
      integer(c_int), parameter :: RKI_ERR_BAD_HANDLE = -100
      integer(c_int) :: rc
      type(c_ptr) :: c_h, c_null

      c_null = c_null_ptr
      rc = rki_step(c_null, 1_c_int)
      call check(error, rc == RKI_ERR_BAD_HANDLE, &
                 "null handle returns RKI_ERR_BAD_HANDLE from step")
      if (allocated(error)) return

      rc = rki_finalize(c_null)
      call check(error, rc == RKI_ERR_BAD_HANDLE, &
                 "null handle returns RKI_ERR_BAD_HANDLE from finalize")
      if (allocated(error)) return

      !! Stale handle: create then destroy, copy survives, must be rejected.
      block
         character(len=64) :: f
         type(c_ptr) :: stale
         f = "test_api_bad_handle.nml"
         call write_default_namelist(f, 20, 8)

         rc = rki_create_from_namelist(trim(f), len_trim(f, kind=c_int), c_h)
         if (rc /= 0) then
            call check(error, .false., "could not create handle for stale test")
            return
         end if
         stale = c_h
         call rki_destroy(c_h)
         rc = rki_step(stale, 1_c_int)
         call check(error, rc == RKI_ERR_BAD_HANDLE, &
                    "stale handle (post-destroy) returns RKI_ERR_BAD_HANDLE")
      end block
   end subroutine test_bad_handle_rejected

   subroutine test_create_from_string_matches_file(error)
      !! The in-memory create path must produce a byte-identical run to the
      !! file path for the same namelist text: build the config as a file
      !! AND as a newline-separated buffer, step both, compare h cell-wise.
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: rc, nx_f, ny_f, nx_s, ny_s
      type(c_ptr) :: c_file, c_str, hp_f, hp_s
      real(c_double), pointer :: hf(:, :), hs(:, :)
      character(len=64) :: f
      character(len=:), allocatable :: buf
      character(len=1), parameter :: NL = achar(10)
      integer :: i, j
      real(c_double) :: maxdiff

      f = "test_api_string.nml"
      call write_default_namelist(f, 40, 12)

      ! Same content as write_default_namelist, as one newline-separated buffer.
      buf = "&grid_nml"//NL// &
            "  nx = 40, ny = 12"//NL// &
            "  dx = 1.0, dy = 1.0"//NL// &
            "/"//NL// &
            "&time_nml"//NL// &
            "  t_end = 10.0, cfl = 0.45"//NL// &
            "/"//NL// &
            "&initial_condition_nml"//NL// &
            '  init_condition = "dam_break"'//NL// &
            "/"//NL// &
            "&logging_nml"//NL// &
            '  log_level = "warning"'//NL// &
            "/"//NL

      rc = rki_create_from_namelist(trim(f), len_trim(f, kind=c_int), c_file)
      call check(error, rc == 0, "file create succeeds")
      if (allocated(error)) then; call rki_destroy(c_file); return; end if

      rc = rki_create_from_namelist_string(buf, len(buf, kind=c_int), c_str)
      call check(error, rc == 0, "string create succeeds")
      if (allocated(error)) then
         call rki_destroy(c_file); call rki_destroy(c_str); return
      end if

      rc = rki_step(c_file, 5_c_int)
      rc = rki_step(c_str, 5_c_int)

      rc = rki_get_h_ptr(c_file, hp_f, nx_f, ny_f)
      rc = rki_get_h_ptr(c_str, hp_s, nx_s, ny_s)
      call check(error, nx_f == nx_s .and. ny_f == ny_s, &
                 "string + file grids match")
      if (allocated(error)) then
         call rki_destroy(c_file); call rki_destroy(c_str); return
      end if

      call c_f_pointer(hp_f, hf, [nx_f, ny_f])
      call c_f_pointer(hp_s, hs, [nx_s, ny_s])
      maxdiff = 0.0_c_double
      do j = 1, ny_f
         do i = 1, nx_f
            maxdiff = max(maxdiff, abs(hf(i, j) - hs(i, j)))
         end do
      end do

      call check(error, maxdiff == 0.0_c_double, &
                 "string-config run is byte-identical to the file-config run")

      call rki_destroy(c_file)
      call rki_destroy(c_str)
   end subroutine test_create_from_string_matches_file

end module test_rki_api
