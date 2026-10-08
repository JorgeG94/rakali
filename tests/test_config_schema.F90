!! Unit tests for rki_config_schema: the central Rakali namelist schema.
module test_config_schema
   use rki_constants, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_config, only: config_t, semi_implicit_grid_conflict
   use rki_config_schema, only: build_rakali_schema
   use rki_nml_schema, only: nml_schema_t
   implicit none
   private

   integer, parameter :: rk = wp
      !! Compare parsed config values (which are `real(wp)`) against
      !! `_rk` literals of the same kind, so the tolerance checks hold in
      !! a single-precision build (`0.3_wp` matches the parsed `0.3`).

   public :: collect_config_schema_tests

contains

   subroutine collect_config_schema_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("happy_path_physics_vcoord", test_happy), &
                  new_unittest("typo_key_suggestion", test_typo), &
                  new_unittest("out_of_range", test_range), &
                  new_unittest("doc_short_non_default", test_doc_short), &
                  new_unittest("no_external_groups", test_external), &
                  new_unittest("semi_implicit_needs_unstructured", test_si_grid) &
                  ]
   end subroutine collect_config_schema_tests

   subroutine test_si_grid(error)
      !! Semi-implicit stepping is unstructured-only: the predicate flags
      !! exactly the (semi_implicit, non-unstructured) combination that
      !! `validate_config` rejects.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call check(error,.not. semi_implicit_grid_conflict(cfg), &
                 "default config (rk2, structured) is not a conflict")
      if (allocated(error)) return

      cfg%time_scheme = "semi_implicit"
      cfg%grid_type = "structured"
      call check(error, semi_implicit_grid_conflict(cfg), &
                 "semi_implicit on a structured grid must be flagged")
      if (allocated(error)) return

      cfg%grid_type = "unstructured"
      call check(error,.not. semi_implicit_grid_conflict(cfg), &
                 "semi_implicit on an unstructured grid is allowed")
      if (allocated(error)) return

      cfg%time_scheme = "rk2"
      cfg%grid_type = "structured"
      call check(error,.not. semi_implicit_grid_conflict(cfg), &
                 "explicit scheme on a structured grid is allowed")
   end subroutine test_si_grid

   subroutine write_file(path, body)
      character(len=*), intent(in) :: path, body
      integer :: u, i, start
      open (newunit=u, file=path, status='replace', action='write')
      start = 1
      do i = 1, len(body)
         if (body(i:i) == achar(10)) then
            write (u, '(A)') body(start:i - 1)
            start = i + 1
         end if
      end do
      if (start <= len(body)) write (u, '(A)') body(start:)
      close (u)
   end subroutine write_file

   subroutine rm_file(path)
      character(len=*), intent(in) :: path
      integer :: u, ios
      open (newunit=u, file=path, status='old', iostat=ios)
      if (ios == 0) close (u, status='delete')
   end subroutine rm_file

   logical function any_contains(arr, n, needle)
      character(len=*), intent(in) :: arr(:)
      integer, intent(in) :: n
      character(len=*), intent(in) :: needle
      integer :: i
      any_contains = .false.
      do i = 1, n
         if (index(arr(i), needle) > 0) then
            any_contains = .true.
            return
         end if
      end do
   end function any_contains

   logical function file_contains(path, needle)
      character(len=*), intent(in) :: path, needle
      integer :: u, ios
      character(len=512) :: buf
      file_contains = .false.
      open (newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (u, '(A)', iostat=ios) buf
         if (ios /= 0) exit
         if (index(buf, needle) > 0) then
            file_contains = .true.
            exit
         end if
      end do
      close (u)
   end function file_contains

   subroutine test_happy(error)
      !! Default config + a file setting physics use_hllc/manning_n + vcoord
      !! remap_method="PPM": fields set, zero errors, enum canonicalized.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_cs_happy.nml"

      call build_rakali_schema(cfg, schema)
      call write_file(fn, &
                      "&physics_nml"//achar(10)// &
                      "  use_hllc = .true."//achar(10)// &
                      "  manning_n = 0.3"//achar(10)// &
                      "/"//achar(10)// &
                      "&vcoord_nml  remap_method = 'PPM_H4' /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, stat == 0, "expected no errors")
      if (allocated(error)) return
      call check(error, cfg%use_hllc, "use_hllc set")
      if (allocated(error)) return
      call check(error, abs(cfg%manning_n - 0.3_rk) < 1.0e-12_rk, "manning_n set")
      if (allocated(error)) return
      ! enum canonical spelling is the lowercase "ppm_h4" from the allowed set.
      call check(error, trim(cfg%remap_method) == "ppm_h4", "enum canonical")
   end subroutine test_happy

   subroutine test_typo(error)
      !! Typo'd physics key -> error with a did-you-mean suggestion.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_cs_typo.nml"

      call build_rakali_schema(cfg, schema)
      call write_file(fn, "&physics_nml  mannin_n = 0.3 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, stat >= 1, "expected an error")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "unknown key"), "unknown-key msg")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "manning_n"), "suggestion present")
   end subroutine test_typo

   subroutine test_range(error)
      !! Out-of-range time integer (cfl_interval < 1) -> error.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_cs_range.nml"

      call build_rakali_schema(cfg, schema)
      call write_file(fn, "&time_nml  cfl_interval = 0 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, stat >= 1, "expected an error")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "below min"), "range violation")
   end subroutine test_range

   subroutine test_doc_short(error)
      !! write_doc_short with ONE non-default knob contains exactly it.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_cs_ds.nml"
      character(len=*), parameter :: sf = "scratch_cs_short.nml"

      call build_rakali_schema(cfg, schema)
      call write_file(fn, "&physics_nml  bottom_drag_r = 0.5 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat == 0, "parse ok")
      if (allocated(error)) return

      call schema%write_doc_short(sf)
      call check(error, file_contains(sf, "bottom_drag_r"), "short has the non-default knob")
      if (allocated(error)) return
      ! dry_threshold is still default -> must be omitted from the short doc.
      call check(error,.not. file_contains(sf, "dry_threshold"), "short omits default knob")
      call rm_file(sf)
   end subroutine test_doc_short

   subroutine test_external(error)
      !! Every group is migrated onto the strict schema — there are no
      !! external (skipped) groups left in this build.  Strict groups
      !! validate known keys, reject unknown keys, and reject bad enum
      !! values.
      type(error_type), allocatable, intent(out) :: error
      type(config_t), target :: cfg
      type(nml_schema_t) :: schema
      integer :: stat, n_ext
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_cs_ext.nml"

      call build_rakali_schema(cfg, schema)

      ! No external groups remain in the coastal-only schema.
      n_ext = 0
      if (allocated(schema%external_names)) n_ext = size(schema%external_names)
      call check(error, n_ext == 0, "no external groups remain")
      if (allocated(error)) return

      ! A clean vcoord_nml validates and applies.
      call write_file(fn, "&vcoord_nml  vcoord_type = 'zstar'  zstar_h_min = 0.001 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat == 0, "clean vcoord_nml validates")
      if (allocated(error)) return

      ! vcoord_nml is strict — an unknown key must error.
      call write_file(fn, "&vcoord_nml  zstar_h_min = 0.001  some_garbage = 3 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat == 1, "vcoord_nml rejects unknown key")
      if (allocated(error)) return

      ! A bad enum value (vcoord_type) must error too.
      call write_file(fn, "&vcoord_nml  vcoord_type = 'nonsense' /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat == 1, "vcoord_nml rejects bad vcoord_type enum")
      if (allocated(error)) return

      ! physics_nml is a strict group — an unknown key must error.
      call write_file(fn, "&physics_nml  manning_n = 0.025  some_garbage = 3 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)
      call check(error, stat == 1, "migrated physics_nml rejects unknown key")
   end subroutine test_external

end module test_config_schema
