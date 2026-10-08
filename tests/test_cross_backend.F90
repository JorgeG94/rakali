!! Cross-backend coexistence test (Phase 4 of the FFI handle refactor).
!!
!! Creates a structured handle and an unstructured handle in the same
!! process, finalizes both, takes one step on each, reads back fields,
!! then destroys both.  This is the test the singleton design literally
!! could not express.
module test_cross_backend
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, &
                                                                             c_f_pointer, c_int, c_loc, &
                                                                             c_null_ptr, c_ptr
   use rki_constants, only: wp
   implicit none
   private

   public :: collect_cross_tests

   ! ----- Structured C interface (subset) -----
   interface
      function rki_create_from_namelist(filename, flen, c_handle) result(ierr) &
         bind(c, name='rki_create_from_namelist')
         import :: c_char, c_int, c_ptr
         character(kind=c_char), intent(in) :: filename(*)
         integer(c_int), intent(in), value :: flen
         type(c_ptr), intent(out) :: c_handle
         integer(c_int) :: ierr
      end function

      function rki_step(c_handle, n) result(ierr) bind(c, name='rki_step')
         import :: c_int, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int), intent(in), value :: n
         integer(c_int) :: ierr
      end function

      function rki_get_step_count(c_handle) result(n) &
         bind(c, name='rki_get_step_count')
         import :: c_int, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int) :: n
      end function

      subroutine rki_destroy(c_handle) bind(c, name='rki_destroy')
         import :: c_ptr
         type(c_ptr), intent(inout) :: c_handle
      end subroutine
   end interface

   ! ----- Unstructured C interface (subset) -----
   interface
      function rki_unstr_create(nverts, ntri, vx, vy, tri_v, bc_tags, c_handle) &
         result(ierr) bind(c, name='rki_unstr_create')
         import :: c_double, c_int, c_ptr
         integer(c_int), value :: nverts
         integer(c_int), value :: ntri
         real(c_double), intent(in) :: vx(nverts)
         real(c_double), intent(in) :: vy(nverts)
         integer(c_int), intent(in) :: tri_v(3, ntri)
         type(c_ptr), value :: bc_tags
         type(c_ptr), intent(out) :: c_handle
         integer(c_int) :: ierr
      end function

      function rki_unstr_set_bathymetry(c_handle, c_bed) result(ierr) &
         bind(c, name='rki_unstr_set_bathymetry')
         import :: c_int, c_ptr
         type(c_ptr), value :: c_handle
         type(c_ptr), value :: c_bed
         integer(c_int) :: ierr
      end function

      function rki_unstr_set_initial_stage(c_handle, c_stage) result(ierr) &
         bind(c, name='rki_unstr_set_initial_stage')
         import :: c_int, c_ptr
         type(c_ptr), value :: c_handle
         type(c_ptr), value :: c_stage
         integer(c_int) :: ierr
      end function

      function rki_unstr_set_time_scheme(c_handle, c_name, c_len) result(ierr) &
         bind(c, name='rki_unstr_set_time_scheme')
         import :: c_int, c_ptr
         type(c_ptr), value :: c_handle
         type(c_ptr), value :: c_name
         integer(c_int), value :: c_len
         integer(c_int) :: ierr
      end function

      function rki_unstr_finalize(c_handle) result(ierr) &
         bind(c, name='rki_unstr_finalize')
         import :: c_int, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int) :: ierr
      end function

      function rki_unstr_step(c_handle, dt) result(ierr) &
         bind(c, name='rki_unstr_step')
         import :: c_double, c_int, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         real(c_double), intent(in), value :: dt
         integer(c_int) :: ierr
      end function

      function rki_unstr_get_step_count(c_handle) result(n) &
         bind(c, name='rki_unstr_get_step_count')
         import :: c_int, c_ptr
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int) :: n
      end function

      subroutine rki_unstr_destroy(c_handle) bind(c, name='rki_unstr_destroy')
         import :: c_ptr
         type(c_ptr), intent(inout) :: c_handle
      end subroutine
   end interface

contains

   subroutine collect_cross_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("structured_and_unstructured_coexist", &
                               test_both_step), &
                  new_unittest("unstr_set_time_scheme_via_api", &
                               test_unstr_set_time_scheme) &
                  ]
   end subroutine collect_cross_tests

   subroutine test_both_step(error)
      !! Run a structured and an unstructured simulation in the same
      !! process. Each takes its own steps, neither sees the other.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: c_struct, c_unstr
      integer(c_int) :: rc, struct_steps, unstr_steps
      character(len=64) :: fname
      integer :: io_unit
      ! Minimal 4-triangle square mesh (same as test_state_unstr).
      real(c_double) :: vx(5), vy(5)
      integer(c_int) :: tri_v(3, 4)
      real(c_double), target :: bed(4), stage(4)

      ! ---- Set up the structured handle from a tiny namelist ----
      fname = "test_cross_struct.nml"
      open (newunit=io_unit, file=trim(fname), status="replace", action="write")
      write (io_unit, '(a)') "&grid_nml"
      write (io_unit, '(a)') "  nx = 20, ny = 8"
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

      rc = rki_create_from_namelist(trim(fname), &
                                    len_trim(fname, kind=c_int), c_struct)
      call check(error, rc == 0, "structured create_from_namelist succeeds")
      if (allocated(error)) then
         call rki_destroy(c_struct)
         return
      end if
      call check(error, c_associated(c_struct), "structured handle is non-null")
      if (allocated(error)) then
         call rki_destroy(c_struct)
         return
      end if

      ! ---- Set up the unstructured handle from a 4-triangle mesh ----
      vx = [0.0_c_double, 1.0_c_double, 1.0_c_double, 0.0_c_double, 0.5_c_double]
      vy = [0.0_c_double, 0.0_c_double, 1.0_c_double, 1.0_c_double, 0.5_c_double]
      tri_v(:, 1) = [1_c_int, 2_c_int, 5_c_int]
      tri_v(:, 2) = [2_c_int, 3_c_int, 5_c_int]
      tri_v(:, 3) = [3_c_int, 4_c_int, 5_c_int]
      tri_v(:, 4) = [4_c_int, 1_c_int, 5_c_int]

      rc = rki_unstr_create(5_c_int, 4_c_int, vx, vy, tri_v, c_null_ptr, c_unstr)
      call check(error, rc == 0, "unstructured create succeeds alongside structured")
      if (allocated(error)) then
         call rki_destroy(c_struct)
         call rki_unstr_destroy(c_unstr)
         return
      end if
      call check(error, c_associated(c_unstr), "unstructured handle is non-null")
      if (allocated(error)) then
         call rki_destroy(c_struct)
         call rki_unstr_destroy(c_unstr)
         return
      end if

      bed = 0.0_c_double
      stage = 1.0_c_double
      rc = rki_unstr_set_bathymetry(c_unstr, c_loc(bed))
      call check(error, rc == 0, "set bathymetry on unstr handle")
      if (allocated(error)) then
         call rki_destroy(c_struct)
         call rki_unstr_destroy(c_unstr)
         return
      end if
      rc = rki_unstr_set_initial_stage(c_unstr, c_loc(stage))
      call check(error, rc == 0, "set initial stage on unstr handle")
      if (allocated(error)) then
         call rki_destroy(c_struct)
         call rki_unstr_destroy(c_unstr)
         return
      end if
      rc = rki_unstr_finalize(c_unstr)
      call check(error, rc == 0, "unstr finalize succeeds with structured live")
      if (allocated(error)) then
         call rki_destroy(c_struct)
         call rki_unstr_destroy(c_unstr)
         return
      end if

      ! ---- Step each independently ----
      rc = rki_step(c_struct, 3_c_int)
      call check(error, rc == 0, "structured step(3) succeeds")
      if (allocated(error)) then
         call rki_destroy(c_struct)
         call rki_unstr_destroy(c_unstr)
         return
      end if

      rc = rki_unstr_step(c_unstr, 0.001_c_double)
      call check(error, rc == 0, "unstr step succeeds with structured live")
      if (allocated(error)) then
         call rki_destroy(c_struct)
         call rki_unstr_destroy(c_unstr)
         return
      end if

      ! ---- Each handle has its own step counter ----
      struct_steps = rki_get_step_count(c_struct)
      unstr_steps = rki_unstr_get_step_count(c_unstr)
      call check(error, struct_steps == 3, &
                 "structured step counter independent (== 3)")
      if (allocated(error)) then
         call rki_destroy(c_struct)
         call rki_unstr_destroy(c_unstr)
         return
      end if
      call check(error, unstr_steps == 1, &
                 "unstr step counter independent (== 1)")
      if (allocated(error)) then
         call rki_destroy(c_struct)
         call rki_unstr_destroy(c_unstr)
         return
      end if

      ! ---- Destroy in reverse order; structured must still work ----
      call rki_unstr_destroy(c_unstr)
      call check(error,.not. c_associated(c_unstr), &
                 "unstr handle nulled after destroy")
      if (allocated(error)) then
         call rki_destroy(c_struct)
         return
      end if

      rc = rki_step(c_struct, 1_c_int)
      call check(error, rc == 0, &
                 "structured handle still steppable after unstr destroyed")
      if (allocated(error)) then
         call rki_destroy(c_struct)
         return
      end if

      call rki_destroy(c_struct)
      call check(error,.not. c_associated(c_struct), &
                 "structured handle nulled after destroy")
   end subroutine test_both_step

   subroutine test_unstr_set_time_scheme(error)
      !! Select a non-default integrator on an unstructured handle through
      !! the C API (rki_unstr_set_time_scheme) and confirm the run steps.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: c_unstr
      integer(c_int) :: rc, nsteps
      real(c_double) :: vx(5), vy(5)
      integer(c_int) :: tri_v(3, 4)
      real(c_double), target :: bed(4), stage(4)
      character(kind=c_char), dimension(3), target :: nm

      vx = [0.0_c_double, 1.0_c_double, 1.0_c_double, 0.0_c_double, 0.5_c_double]
      vy = [0.0_c_double, 0.0_c_double, 1.0_c_double, 1.0_c_double, 0.5_c_double]
      tri_v(:, 1) = [1_c_int, 2_c_int, 5_c_int]
      tri_v(:, 2) = [2_c_int, 3_c_int, 5_c_int]
      tri_v(:, 3) = [3_c_int, 4_c_int, 5_c_int]
      tri_v(:, 4) = [4_c_int, 1_c_int, 5_c_int]

      rc = rki_unstr_create(5_c_int, 4_c_int, vx, vy, tri_v, c_null_ptr, c_unstr)
      call check(error, rc == 0, "unstr create for time-scheme test")
      if (allocated(error)) return

      ! Select rk3 via the C API before finalize.
      nm(1) = "r"
      nm(2) = "k"
      nm(3) = "3"
      rc = rki_unstr_set_time_scheme(c_unstr, c_loc(nm), 3_c_int)
      call check(error, rc == 0, "rki_unstr_set_time_scheme(rk3) returns 0")
      if (allocated(error)) then
         call rki_unstr_destroy(c_unstr)
         return
      end if

      bed = 0.0_c_double
      stage = 1.0_c_double
      rc = rki_unstr_set_bathymetry(c_unstr, c_loc(bed))
      if (rc == 0) rc = rki_unstr_set_initial_stage(c_unstr, c_loc(stage))
      if (rc == 0) rc = rki_unstr_finalize(c_unstr)
      call check(error, rc == 0, "unstr setup succeeds with rk3 selected")
      if (allocated(error)) then
         call rki_unstr_destroy(c_unstr)
         return
      end if

      ! The rk3 path runs through the C-created handle.
      rc = rki_unstr_step(c_unstr, 0.001_c_double)
      call check(error, rc == 0, "unstr step under rk3 succeeds")
      if (allocated(error)) then
         call rki_unstr_destroy(c_unstr)
         return
      end if

      nsteps = rki_unstr_get_step_count(c_unstr)
      call check(error, nsteps == 1, "rk3 step advanced the counter")

      call rki_unstr_destroy(c_unstr)
   end subroutine test_unstr_set_time_scheme

end module test_cross_backend
