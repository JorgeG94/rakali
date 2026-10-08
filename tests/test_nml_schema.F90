!! Unit tests for the rki_nml_schema strict namelist validation engine.
module test_nml_schema
   use rki_constants, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_nml_schema, only: nml_schema_t, nml_group_t, &
                             nml_real, nml_int, nml_logical, nml_string, &
                             nml_enum, nml_real_array, group_check_iface
   implicit none
   private

   integer, parameter :: rk = wp
      !! Mirror the schema's working real kind (`wp`) so the pointer
      !! targets below associate in either a single- or double-precision
      !! build (the schema's real keys are `real(wp)`).

   public :: collect_nml_schema_tests

   ! Module-level config target slots shared across helper builders.
   ! Pointers in the schema associate with these.  Each test resets them.
   real(rk), target :: cfg_ri_crit
   integer, target :: cfg_n_iter
   logical, target :: cfg_use_lt
   character(len=32), target :: cfg_form
   character(len=16), target :: cfg_variant
   real(rk), target :: cfg_layers(4)

   ! Cross-check shared state.
   real(rk), target :: cfg_lo, cfg_hi

contains

   subroutine collect_nml_schema_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("happy_path_all_types", test_happy_path), &
                  new_unittest("unknown_key_suggestion", test_unknown_key), &
                  new_unittest("unknown_group_and_external", test_unknown_group), &
                  new_unittest("range_enum_strlen", test_range_enum_strlen), &
                  new_unittest("scalar_multivalue_and_bad_syntax", test_bad_syntax), &
                  new_unittest("duplicate_key_and_group", test_duplicates), &
                  new_unittest("required_key_semantics", test_required), &
                  new_unittest("defaults_tracking", test_defaults), &
                  new_unittest("doc_writers_smoke", test_doc_writers), &
                  new_unittest("cross_check_propagation", test_cross_check) &
                  ]
   end subroutine collect_nml_schema_tests

   ! ---- helpers --------------------------------------------------------

   subroutine reset_targets()
      !! Restore config slots to their default values (the schema
      !! captures these as defaults at registration).
      cfg_ri_crit = 0.25_rk
      cfg_n_iter = 10
      cfg_use_lt = .false.
      cfg_form = "sadourny"
      cfg_variant = "om4"
      cfg_layers = [1000.0_rk, 1000.0_rk, 0.0_rk, 0.0_rk]
   end subroutine reset_targets

   subroutine build_schema(schema)
      !! Build a representative schema over the module config slots.
      type(nml_schema_t), intent(out) :: schema
      type(nml_group_t) :: g
      real(rk), pointer :: pr, pa(:)
      integer, pointer :: pi
      logical, pointer :: pl
      character(len=:), pointer :: ps

      call reset_targets()

      g%name = "demo_mixing"
      g%doc = "Synthetic mixing knobs (schema-engine fixture)."
      pr => cfg_ri_crit
      call g%add(nml_real("ri_crit", pr, "Critical Richardson number", &
                          units="nondim", min=0.0_rk, max=2.0_rk))
      pi => cfg_n_iter
      call g%add(nml_int("n_iter", pi, "Solver iterations", min=1, max=100))
      pl => cfg_use_lt
      call g%add(nml_logical("use_lt", pl, "Enable Langmuir"))
      ps => str_ptr_form()
      call g%add(nml_string("form", ps, "Coriolis form"))
      pa => cfg_layers
      call g%add(nml_real_array("layers", pa, "Layer thicknesses", units="m"))
      call schema%add_group(g)

      block
         type(nml_group_t) :: g2
         character(len=:), pointer :: ps2
         g2%name = "demo_pbl"
         g2%doc = "Synthetic PBL knobs (schema-engine fixture)."
         ps2 => str_ptr_variant()
         call g2%add(nml_enum("mstar_variant", ps2, "mstar scheme", &
                              allowed=[character(len=8) :: "constant", "om4", "rh18"]))
         call schema%add_group(g2)
      end block
   end subroutine build_schema

   function str_ptr_form() result(p)
      !! Deferred-len pointer associated with the fixed-len form slot.
      character(len=:), pointer :: p
      p => cfg_form
   end function str_ptr_form

   function str_ptr_variant() result(p)
      !! Deferred-len pointer associated with the fixed-len variant slot.
      character(len=:), pointer :: p
      p => cfg_variant
   end function str_ptr_variant

   logical function any_contains(arr, n, needle)
      !! True if any of arr(1:n) contains the substring `needle`.
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

   subroutine write_file(path, body)
      !! Write `body` (newline-separated) to a scratch namelist file.
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
   ! ---- tests ----------------------------------------------------------

   subroutine test_happy_path(error)
      !! All key types parse; enum canonicalizes; case-insensitive
      !! groups+keys; inline comments + multi-line group bodies work.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_happy.nml"

      call build_schema(schema)
      call write_file(fn, &
                      "&DEMO_MIXING_NML  ! a group"//achar(10)// &
                      "  RI_CRIT = 0.5  ! inline comment"//achar(10)// &
                      "  n_iter = 20,"//achar(10)// &
                      "  use_lt = .true."//achar(10)// &
                      "  form = 'leith'"//achar(10)// &
                      "  layers = 500.0, 250.0"//achar(10)// &
                      "/"//achar(10)// &
                      "&demo_pbl_nml mstar_variant = 'RH18' /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, stat == 0, "expected no errors")
      if (allocated(error)) return
      call check(error, abs(cfg_ri_crit - 0.5_rk) < 1.0e-12_rk, "ri_crit")
      if (allocated(error)) return
      call check(error, cfg_n_iter == 20, "n_iter")
      if (allocated(error)) return
      call check(error, cfg_use_lt, "use_lt")
      if (allocated(error)) return
      call check(error, trim(cfg_form) == "leith", "form")
      if (allocated(error)) return
      call check(error, abs(cfg_layers(1) - 500.0_rk) < 1.0e-9_rk, "layer1")
      if (allocated(error)) return
      call check(error, abs(cfg_layers(2) - 250.0_rk) < 1.0e-9_rk, "layer2")
      if (allocated(error)) return
      ! enum canonical spelling from allowed list (lowercase "rh18")
      call check(error, trim(cfg_variant) == "rh18", "enum canonical")
   end subroutine test_happy_path

   subroutine test_unknown_key(error)
      !! Unknown key produces an error with a did-you-mean suggestion.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_unk_key.nml"

      call build_schema(schema)
      call write_file(fn, "&demo_mixing_nml  ri_crt = 0.3 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, stat >= 1, "expected an error")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "unknown key"), "unknown-key msg")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "ri_crit"), "suggestion present")
   end subroutine test_unknown_key

   subroutine test_unknown_group(error)
      !! Unknown group -> error + suggestion; external group body skipped.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_unk_grp.nml"

      call build_schema(schema)
      call schema%add_external_group("time")
      call write_file(fn, &
                      "&time  dt = 1200.0  some_garbage = 3 /"//achar(10)// &
                      "&demo_mixng_nml  ri_crit = 0.3 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, stat >= 1, "expected an error")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "unknown group"), "unknown-group msg")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "demo_mixing"), "group suggestion")
      if (allocated(error)) return
      ! external 'time' body skipped: no error mentioning some_garbage
      call check(error,.not. any_contains(errs, stat, "some_garbage"), "external skipped")
   end subroutine test_unknown_group

   subroutine test_range_enum_strlen(error)
      !! min/max range, enum violation, and string-too-long all error.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_range.nml"
      character(len=40) :: long_str

      long_str = repeat('x', 40)  ! longer than cfg_form (len 32)
      call build_schema(schema)
      call write_file(fn, &
                      "&demo_mixing_nml"//achar(10)// &
                      "  ri_crit = -1.0"//achar(10)// &
                      "  n_iter = 999"//achar(10)// &
                      "  form = '"//trim(long_str)//"'"//achar(10)// &
                      "/"//achar(10)// &
                      "&demo_pbl_nml mstar_variant = 'bogus' /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, any_contains(errs, stat, "below min"), "min violation")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "above max"), "max violation")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "exceeds target length"), "strlen")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "not in allowed set"), "enum")
   end subroutine test_range_enum_strlen

   subroutine test_bad_syntax(error)
      !! Scalar-with-2-values, repeat-count, indexed assignment rejected.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_syntax.nml"

      call build_schema(schema)
      call write_file(fn, &
                      "&demo_mixing_nml"//achar(10)// &
                      "  ri_crit = 0.3 0.4"//achar(10)// &
                      "  n_iter = 3*5"//achar(10)// &
                      "  layers(2) = 7.0"//achar(10)// &
                      "/")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, any_contains(errs, stat, "expects 1 value"), "scalar 2 values")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "repeat-count"), "repeat count")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "indexed assignment"), "indexed")
   end subroutine test_bad_syntax
   subroutine test_duplicates(error)
      !! Duplicate key within a group + duplicate group occurrence error.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_dup.nml"

      call build_schema(schema)
      call write_file(fn, &
                      "&demo_mixing_nml"//achar(10)// &
                      "  ri_crit = 0.3"//achar(10)// &
                      "  ri_crit = 0.4"//achar(10)// &
                      "/"//achar(10)// &
                      "&demo_mixing_nml  n_iter = 5 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, any_contains(errs, stat, "duplicate key"), "dup key")
      if (allocated(error)) return
      call check(error, any_contains(errs, stat, "duplicate occurrence"), "dup group")
   end subroutine test_duplicates

   subroutine test_required(error)
      !! Required key: errors when group present but key absent; OK when
      !! the group is absent entirely (defaults rule).
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      type(nml_group_t) :: g
      integer :: stat
      character(len=:), allocatable :: errs(:)
      real(rk), pointer :: pr
      character(len=*), parameter :: fn1 = "scratch_req1.nml"
      character(len=*), parameter :: fn2 = "scratch_req2.nml"

      call reset_targets()
      g%name = "reqgrp"
      pr => cfg_ri_crit
      call g%add(nml_real("must_set", pr, "Mandatory knob", required=.true.))
      call schema%add_group(g)

      ! Group present, key absent -> error.
      call write_file(fn1, "&reqgrp /")
      call schema%parse(fn1, status=stat, errors=errs)
      call rm_file(fn1)
      call check(error, any_contains(errs, stat, "required key"), "required fires")
      if (allocated(error)) return

      ! Group absent entirely -> OK.
      call write_file(fn2, "! nothing here")
      call schema%parse(fn2, status=stat, errors=errs)
      call rm_file(fn2)
      call check(error, stat == 0, "absent group OK")
   end subroutine test_required

   subroutine test_defaults(error)
      !! Untouched keys keep defaults + is_default true; after setting
      !! one knob, only that key reports non-default.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: fn = "scratch_def.nml"

      call build_schema(schema)
      ! Parse a file that sets only n_iter.
      call write_file(fn, "&demo_mixing_nml  n_iter = 42 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, stat == 0, "no errors")
      if (allocated(error)) return
      call check(error, cfg_n_iter == 42, "n_iter set")
      if (allocated(error)) return
      ! ri_crit untouched -> still default value.
      call check(error, abs(cfg_ri_crit - 0.25_rk) < 1.0e-12_rk, "ri_crit default kept")
      if (allocated(error)) return
      ! is_default queried via the key objects.
      call check(error, schema%groups(1)%keys(1)%key%is_default(), "ri_crit is_default")
      if (allocated(error)) return
      call check(error,.not. schema%groups(1)%keys(2)%key%is_default(), "n_iter non-default")
   end subroutine test_defaults

   subroutine test_doc_writers(error)
      !! write_doc_all / write_doc_short / render_markdown smoke checks.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      integer :: stat
      character(len=:), allocatable :: errs(:)
      character(len=*), parameter :: nf = "scratch_dw.nml"
      character(len=*), parameter :: af = "scratch_all.nml"
      character(len=*), parameter :: sf = "scratch_short.nml"
      character(len=*), parameter :: mf = "scratch_md.md"

      call build_schema(schema)
      call write_file(nf, "&demo_mixing_nml  n_iter = 77 /")
      call schema%parse(nf, status=stat, errors=errs)
      call rm_file(nf)
      call check(error, stat == 0, "parse ok")
      if (allocated(error)) return

      call schema%write_doc_all(af)
      call schema%write_doc_short(sf)
      call schema%render_markdown(mf)

      ! doc_all has all keys; ri_crit (default) marked (default).
      call check(error, file_contains(af, "ri_crit"), "all has ri_crit")
      if (allocated(error)) return
      call check(error, file_contains(af, "(default)"), "all marks default")
      if (allocated(error)) return
      ! doc_short has only the non-default key.
      call check(error, file_contains(sf, "n_iter"), "short has n_iter")
      if (allocated(error)) return
      call check(error,.not. file_contains(sf, "ri_crit"), "short omits default")
      if (allocated(error)) return
      ! markdown table heading.
      call check(error, file_contains(mf, "| Knob |"), "md table header")

      call rm_file(af)
      call rm_file(sf)
      call rm_file(mf)
   end subroutine test_doc_writers

   logical function file_contains(path, needle)
      !! True if any line of `path` contains `needle`.
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

   function demo_cross_check() result(msg)
      !! Cross-check: cfg_lo must not exceed cfg_hi.
      character(len=:), allocatable :: msg
      msg = ''
      if (cfg_lo > cfg_hi) msg = "lo exceeds hi"
   end function demo_cross_check

   subroutine test_cross_check(error)
      !! A registered cross_check error is collected.
      type(error_type), allocatable, intent(out) :: error
      type(nml_schema_t) :: schema
      type(nml_group_t) :: g
      integer :: stat
      character(len=:), allocatable :: errs(:)
      real(rk), pointer :: plo, phi
      character(len=*), parameter :: fn = "scratch_xc.nml"

      cfg_lo = 5.0_rk
      cfg_hi = 1.0_rk
      g%name = "xcgrp"
      plo => cfg_lo
      phi => cfg_hi
      call g%add(nml_real("lo", plo, "low"))
      call g%add(nml_real("hi", phi, "high"))
      g%cross_check => demo_cross_check
      call schema%add_group(g)

      call write_file(fn, "&xcgrp  lo = 5.0  hi = 1.0 /")
      call schema%parse(fn, status=stat, errors=errs)
      call rm_file(fn)

      call check(error, any_contains(errs, stat, "lo exceeds hi"), "cross_check msg")
   end subroutine test_cross_check

end module test_nml_schema
