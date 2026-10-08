!! Strict, schema-driven namelist validation engine.
module rki_nml_schema
   !! A host-only, schema-driven namelist parser and validator. Validates a
   !! namelist file against a registered schema of groups and keys,
   !! collecting ALL errors (with `path:line:` prefixes, did-you-mean
   !! suggestions, range/enum/length checks) and rejecting unsupported
   !! syntax (repeat-counts, indexed assignment, derived-type refs).
   !!
   !! Defaults contract (load-bearing): every key constructor captures its
   !! default from the CURRENT target value at registration; the schema
   !! NEVER takes a default argument. Register AFTER config defaults exist,
   !! BEFORE parsing.
   !!
   !! Required-key semantics: a `required=.true.` key must appear whenever
   !! its GROUP appears; an absent group is fine (defaults rule).
   use rki_constants, only: wp
   use pic_logger, only: global_logger
   implicit none
   private

   integer, parameter :: rk = wp
      !! Working real kind — tracks the solver `wp` so `real(wp)` config
      !! targets associate in either single- or double-precision builds.

   integer, parameter :: MAX_SUGGEST_LEN = 64
      !! Cap on names fed to the Levenshtein suggester.

   public :: nml_key_t
   public :: nml_real_key_t, nml_int_key_t, nml_logical_key_t
   public :: nml_string_key_t, nml_enum_key_t, nml_real_array_key_t
   public :: nml_group_t, nml_schema_t, nml_key_box_t
   public :: nml_real, nml_int, nml_logical, nml_string, nml_enum, nml_real_array
   public :: group_check_iface

   type, abstract :: nml_key_t
      !! Abstract base for a single namelist key descriptor.
      character(len=:), allocatable :: name
         !! Key name (matched case-insensitively).
      character(len=:), allocatable :: doc
         !! Human-readable description (FORD-style synthesis doc).
      character(len=:), allocatable :: units
         !! Optional units string, default ''.
      logical :: required = .false.
         !! If true, the key must appear whenever its group appears.
      logical :: found = .false.
         !! Parse-state: set true when the key is seen in the file.
   contains
      procedure(parse_tokens_iface), deferred :: parse_tokens
      procedure(value_string_iface), deferred :: value_string
      procedure(default_string_iface), deferred :: default_string
      procedure(is_default_iface), deferred :: is_default
   end type nml_key_t

   abstract interface
      subroutine parse_tokens_iface(this, values, err_msg)
         !! Validate `values` and assign through the target pointer.
         !! On success `err_msg` is left unallocated.
         import :: nml_key_t
         implicit none
         class(nml_key_t), intent(inout) :: this
         character(len=*), intent(in) :: values(:)
         character(len=:), allocatable, intent(out) :: err_msg
      end subroutine parse_tokens_iface

      function value_string_iface(this) result(str)
         !! Current target value, formatted for namelist output.
         import :: nml_key_t
         implicit none
         class(nml_key_t), intent(in) :: this
         character(len=:), allocatable :: str
      end function value_string_iface

      function default_string_iface(this) result(str)
         !! Captured default value, formatted for namelist output.
         import :: nml_key_t
         implicit none
         class(nml_key_t), intent(in) :: this
         character(len=:), allocatable :: str
      end function default_string_iface

      function is_default_iface(this) result(yes)
         !! True if the current target value equals the captured default.
         import :: nml_key_t
         implicit none
         class(nml_key_t), intent(in) :: this
         logical :: yes
      end function is_default_iface
   end interface

   abstract interface
      function group_check_iface() result(msg)
         !! Cross-key validation callback.  Empty string = OK; a
         !! non-empty string is appended to the collected errors.
         implicit none
         character(len=:), allocatable :: msg
      end function group_check_iface
   end interface

   type :: nml_key_box_t
      !! Boxes a polymorphic key so groups can hold a growable array.
      class(nml_key_t), allocatable :: key
   end type nml_key_box_t

   type, extends(nml_key_t) :: nml_real_key_t
      !! Real-valued key with optional [min,max] range validation.
      real(rk), pointer :: tgt => null()
      real(rk) :: default = 0.0_rk
      real(rk) :: vmin = 0.0_rk, vmax = 0.0_rk
      logical :: has_min = .false., has_max = .false.
   contains
      procedure :: parse_tokens => real_parse
      procedure :: value_string => real_value_string
      procedure :: default_string => real_default_string
      procedure :: is_default => real_is_default
   end type nml_real_key_t

   type, extends(nml_key_t) :: nml_int_key_t
      !! Integer-valued key with optional [min,max] range validation.
      integer, pointer :: tgt => null()
      integer :: default = 0
      integer :: vmin = 0, vmax = 0
      logical :: has_min = .false., has_max = .false.
   contains
      procedure :: parse_tokens => int_parse
      procedure :: value_string => int_value_string
      procedure :: default_string => int_default_string
      procedure :: is_default => int_is_default
   end type nml_int_key_t

   type, extends(nml_key_t) :: nml_logical_key_t
      !! Logical-valued key.
      logical, pointer :: tgt => null()
      logical :: default = .false.
   contains
      procedure :: parse_tokens => logical_parse
      procedure :: value_string => logical_value_string
      procedure :: default_string => logical_default_string
      procedure :: is_default => logical_is_default
   end type nml_logical_key_t

   type, extends(nml_key_t) :: nml_string_key_t
      !! String key.  Deferred-len char pointer associates with a
      !! fixed-len target and assumes its length; parsing errors if the
      !! parsed value is LONGER than len(tgt) (no silent truncation).
      character(len=:), pointer :: tgt => null()
      character(len=:), allocatable :: default
   contains
      procedure :: parse_tokens => string_parse
      procedure :: value_string => string_value_string
      procedure :: default_string => string_default_string
      procedure :: is_default => string_is_default
   end type nml_string_key_t

   type, extends(nml_string_key_t) :: nml_enum_key_t
      !! Enum key: a string restricted to an allowed list, validated
      !! case-insensitively and assigned its canonical spelling.
      character(len=:), allocatable :: allowed(:)
   contains
      procedure :: parse_tokens => enum_parse
   end type nml_enum_key_t

   type, extends(nml_key_t) :: nml_real_array_key_t
      !! Real-array key: accepts 1..size(tgt) values, fills from
      !! element 1, errors if more than size(tgt) values are given.
      real(rk), pointer :: tgt(:) => null()
      real(rk), allocatable :: default(:)
   contains
      procedure :: parse_tokens => real_array_parse
      procedure :: value_string => real_array_value_string
      procedure :: default_string => real_array_default_string
      procedure :: is_default => real_array_is_default
   end type nml_real_array_key_t

   type :: nml_group_t
      !! A namelist group: a name, doc, growable key list, parse-state,
      !! and an optional cross-key validation callback.
      character(len=:), allocatable :: name
         !! Group name (matched case-insensitively).
      character(len=:), allocatable :: doc
         !! Group description.
      type(nml_key_box_t), allocatable :: keys(:)
         !! Registered keys.
      logical :: found = .false.
         !! Parse-state: true when the group appears in the file.
      procedure(group_check_iface), pointer, nopass :: cross_check => null()
         !! Optional cross-key validator (empty string = OK).
   contains
      procedure :: add => group_add
      procedure :: find_key => group_find_key
   end type nml_group_t

   type :: nml_schema_t
      !! The top-level schema: validated groups plus external (skipped)
      !! group names.
      type(nml_group_t), allocatable :: groups(:)
         !! Registered, validated groups.
      character(len=:), allocatable :: external_names(:)
         !! Names of groups whose bodies are skipped without validation.
   contains
      procedure :: add_group => schema_add_group
      procedure :: add_external_group => schema_add_external_group
      procedure :: parse => schema_parse
      procedure :: parse_lines => schema_parse_lines
      procedure :: write_doc_all => schema_write_doc_all
      procedure :: write_doc_short => schema_write_doc_short
      procedure :: render_markdown => schema_render_markdown
      procedure :: find_group => schema_find_group
      procedure :: is_external => schema_is_external
   end type nml_schema_t
contains

   ! ------------------------------------------------------------------
   ! Constructors.  Each captures `default = current target value`.
   ! ------------------------------------------------------------------

   function nml_real(name, tgt, doc, units, required, min, max) result(key)
      !! Construct a real key, capturing the current target as default.
      character(len=*), intent(in) :: name, doc
      real(rk), pointer, intent(in) :: tgt
      character(len=*), intent(in), optional :: units
      logical, intent(in), optional :: required
      real(rk), intent(in), optional :: min, max
      type(nml_real_key_t) :: key

      key%name = name
      key%doc = doc
      key%units = ""
      if (present(units)) key%units = units
      if (present(required)) key%required = required
      key%tgt => tgt
      key%default = tgt
      if (present(min)) then
         key%vmin = min
         key%has_min = .true.
      end if
      if (present(max)) then
         key%vmax = max
         key%has_max = .true.
      end if
   end function nml_real

   function nml_int(name, tgt, doc, units, required, min, max) result(key)
      !! Construct an integer key, capturing the current target as default.
      character(len=*), intent(in) :: name, doc
      integer, pointer, intent(in) :: tgt
      character(len=*), intent(in), optional :: units
      logical, intent(in), optional :: required
      integer, intent(in), optional :: min, max
      type(nml_int_key_t) :: key

      key%name = name
      key%doc = doc
      key%units = ""
      if (present(units)) key%units = units
      if (present(required)) key%required = required
      key%tgt => tgt
      key%default = tgt
      if (present(min)) then
         key%vmin = min
         key%has_min = .true.
      end if
      if (present(max)) then
         key%vmax = max
         key%has_max = .true.
      end if
   end function nml_int

   function nml_logical(name, tgt, doc, units, required) result(key)
      !! Construct a logical key, capturing the current target as default.
      character(len=*), intent(in) :: name, doc
      logical, pointer, intent(in) :: tgt
      character(len=*), intent(in), optional :: units
      logical, intent(in), optional :: required
      type(nml_logical_key_t) :: key

      key%name = name
      key%doc = doc
      key%units = ""
      if (present(units)) key%units = units
      if (present(required)) key%required = required
      key%tgt => tgt
      key%default = tgt
   end function nml_logical

   function nml_string(name, tgt, doc, units, required) result(key)
      !! Construct a string key, capturing the current target as default.
      character(len=*), intent(in) :: name, doc
      character(len=:), pointer, intent(in) :: tgt
      character(len=*), intent(in), optional :: units
      logical, intent(in), optional :: required
      type(nml_string_key_t) :: key

      key%name = name
      key%doc = doc
      key%units = ""
      if (present(units)) key%units = units
      if (present(required)) key%required = required
      key%tgt => tgt
      key%default = tgt
   end function nml_string

   function nml_enum(name, tgt, doc, allowed, units, required) result(key)
      !! Construct an enum key restricted to `allowed`, capturing the
      !! current target as default.
      character(len=*), intent(in) :: name, doc
      character(len=:), pointer, intent(in) :: tgt
      character(len=*), intent(in) :: allowed(:)
      character(len=*), intent(in), optional :: units
      logical, intent(in), optional :: required
      type(nml_enum_key_t) :: key

      key%name = name
      key%doc = doc
      key%units = ""
      if (present(units)) key%units = units
      if (present(required)) key%required = required
      key%tgt => tgt
      key%default = tgt
      ! Explicit allocate + element copy on purpose: the auto-allocating
      ! `key%allowed = allowed` makes NVHPC's str_copy write one byte past
      ! the allocation (heap corruption once enough enums register).
      block
         integer :: i
         allocate (character(len=len(allowed)) :: key%allowed(size(allowed)))
         do i = 1, size(allowed)
            key%allowed(i) = allowed(i)
         end do
      end block
   end function nml_enum

   function nml_real_array(name, tgt, doc, units, required) result(key)
      !! Construct a real-array key, capturing the current target as
      !! default (whole-array compare for is_default).
      character(len=*), intent(in) :: name, doc
      real(rk), pointer, intent(in) :: tgt(:)
      character(len=*), intent(in), optional :: units
      logical, intent(in), optional :: required
      type(nml_real_array_key_t) :: key

      key%name = name
      key%doc = doc
      key%units = ""
      if (present(units)) key%units = units
      if (present(required)) key%required = required
      key%tgt => tgt
      key%default = tgt
   end function nml_real_array
   ! ------------------------------------------------------------------
   ! Formatting helpers.
   ! ------------------------------------------------------------------

   pure function fmt_real(v) result(str)
      !! Format a real in a stable, namelist-valid exponential form.
      real(rk), intent(in) :: v
      character(len=:), allocatable :: str
      character(len=32) :: buf
      write (buf, "(E18.10)") v
      str = trim(adjustl(buf))
   end function fmt_real

   pure function fmt_int(v) result(str)
      !! Format an integer.
      integer, intent(in) :: v
      character(len=:), allocatable :: str
      character(len=32) :: buf
      write (buf, "(I0)") v
      str = trim(adjustl(buf))
   end function fmt_int

   pure function fmt_logical(v) result(str)
      !! Format a logical as namelist `.true.`/`.false.`.
      logical, intent(in) :: v
      character(len=:), allocatable :: str
      if (v) then
         str = ".true."
      else
         str = ".false."
      end if
   end function fmt_logical

   pure function lower(s) result(out)
      !! Lowercase an ASCII string.
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
      integer :: i, c
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) then
            out(i:i) = achar(c + 32)
         else
            out(i:i) = s(i:i)
         end if
      end do
   end function lower

   ! ------------------------------------------------------------------
   ! Real key methods.
   ! ------------------------------------------------------------------

   subroutine real_parse(this, values, err_msg)
      class(nml_real_key_t), intent(inout) :: this
      character(len=*), intent(in) :: values(:)
      character(len=:), allocatable, intent(out) :: err_msg
      real(rk) :: v
      integer :: ios
      character(len=128) :: imsg

      if (size(values) /= 1) then
         err_msg = "key '"//this%name//"' expects 1 value, got "//fmt_int(size(values))
         return
      end if
      read (values(1), *, iostat=ios, iomsg=imsg) v
      if (ios /= 0) then
         err_msg = "key '"//this%name//"': cannot parse real from '"//trim(values(1))//"'"
         return
      end if
      if (this%has_min .and. v < this%vmin) then
         err_msg = "key '"//this%name//"' = "//fmt_real(v)//" below min "//fmt_real(this%vmin)
         return
      end if
      if (this%has_max .and. v > this%vmax) then
         err_msg = "key '"//this%name//"' = "//fmt_real(v)//" above max "//fmt_real(this%vmax)
         return
      end if
      this%tgt = v
   end subroutine real_parse

   function real_value_string(this) result(str)
      class(nml_real_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = fmt_real(this%tgt)
   end function real_value_string

   function real_default_string(this) result(str)
      class(nml_real_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = fmt_real(this%default)
   end function real_default_string

   function real_is_default(this) result(yes)
      class(nml_real_key_t), intent(in) :: this
      logical :: yes
      yes = (this%tgt == this%default)
   end function real_is_default

   ! ------------------------------------------------------------------
   ! Integer key methods.
   ! ------------------------------------------------------------------

   subroutine int_parse(this, values, err_msg)
      class(nml_int_key_t), intent(inout) :: this
      character(len=*), intent(in) :: values(:)
      character(len=:), allocatable, intent(out) :: err_msg
      integer :: v, ios
      character(len=128) :: imsg

      if (size(values) /= 1) then
         err_msg = "key '"//this%name//"' expects 1 value, got "//fmt_int(size(values))
         return
      end if
      read (values(1), *, iostat=ios, iomsg=imsg) v
      if (ios /= 0) then
         err_msg = "key '"//this%name//"': cannot parse integer from '"//trim(values(1))//"'"
         return
      end if
      if (this%has_min .and. v < this%vmin) then
         err_msg = "key '"//this%name//"' = "//fmt_int(v)//" below min "//fmt_int(this%vmin)
         return
      end if
      if (this%has_max .and. v > this%vmax) then
         err_msg = "key '"//this%name//"' = "//fmt_int(v)//" above max "//fmt_int(this%vmax)
         return
      end if
      this%tgt = v
   end subroutine int_parse

   function int_value_string(this) result(str)
      class(nml_int_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = fmt_int(this%tgt)
   end function int_value_string

   function int_default_string(this) result(str)
      class(nml_int_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = fmt_int(this%default)
   end function int_default_string

   function int_is_default(this) result(yes)
      class(nml_int_key_t), intent(in) :: this
      logical :: yes
      yes = (this%tgt == this%default)
   end function int_is_default
   ! ------------------------------------------------------------------
   ! Logical key methods.
   ! ------------------------------------------------------------------

   subroutine logical_parse(this, values, err_msg)
      class(nml_logical_key_t), intent(inout) :: this
      character(len=*), intent(in) :: values(:)
      character(len=:), allocatable, intent(out) :: err_msg
      character(len=:), allocatable :: tok

      if (size(values) /= 1) then
         err_msg = "key '"//this%name//"' expects 1 value, got "//fmt_int(size(values))
         return
      end if
      tok = lower(trim(adjustl(values(1))))
      select case (tok)
      case (".true.", "t", "true")
         this%tgt = .true.
      case (".false.", "f", "false")
         this%tgt = .false.
      case default
         err_msg = "key '"//this%name//"': cannot parse logical from '"//trim(values(1))//"'"
      end select
   end subroutine logical_parse

   function logical_value_string(this) result(str)
      class(nml_logical_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = fmt_logical(this%tgt)
   end function logical_value_string

   function logical_default_string(this) result(str)
      class(nml_logical_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = fmt_logical(this%default)
   end function logical_default_string

   function logical_is_default(this) result(yes)
      class(nml_logical_key_t), intent(in) :: this
      logical :: yes
      yes = (this%tgt .eqv. this%default)
   end function logical_is_default

   ! ------------------------------------------------------------------
   ! String key methods.
   ! ------------------------------------------------------------------

   subroutine string_parse(this, values, err_msg)
      class(nml_string_key_t), intent(inout) :: this
      character(len=*), intent(in) :: values(:)
      character(len=:), allocatable, intent(out) :: err_msg
      character(len=:), allocatable :: v

      if (size(values) /= 1) then
         err_msg = "key '"//this%name//"' expects 1 value, got "//fmt_int(size(values))
         return
      end if
      v = trim(values(1))
      if (len(v) > len(this%tgt)) then
         err_msg = "key '"//this%name//"': value '"//v//"' ("//fmt_int(len(v))// &
                   " chars) exceeds target length "//fmt_int(len(this%tgt))
         return
      end if
      this%tgt = v
   end subroutine string_parse

   function string_value_string(this) result(str)
      class(nml_string_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = '"'//trim(this%tgt)//'"'
   end function string_value_string

   function string_default_string(this) result(str)
      class(nml_string_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      str = '"'//trim(this%default)//'"'
   end function string_default_string

   function string_is_default(this) result(yes)
      class(nml_string_key_t), intent(in) :: this
      logical :: yes
      yes = (trim(this%tgt) == trim(this%default))
   end function string_is_default

   ! ------------------------------------------------------------------
   ! Enum key methods (extends string).
   ! ------------------------------------------------------------------

   subroutine enum_parse(this, values, err_msg)
      class(nml_enum_key_t), intent(inout) :: this
      character(len=*), intent(in) :: values(:)
      character(len=:), allocatable, intent(out) :: err_msg
      character(len=:), allocatable :: v, list
      integer :: i

      if (size(values) /= 1) then
         err_msg = "key '"//this%name//"' expects 1 value, got "//fmt_int(size(values))
         return
      end if
      v = trim(adjustl(values(1)))
      do i = 1, size(this%allowed)
         if (lower(v) == lower(trim(this%allowed(i)))) then
            if (len(trim(this%allowed(i))) > len(this%tgt)) then
               err_msg = "key '"//this%name//"': canonical value '"//trim(this%allowed(i))// &
                         "' exceeds target length "//fmt_int(len(this%tgt))
               return
            end if
            this%tgt = trim(this%allowed(i))
            return
         end if
      end do
      list = ""
      do i = 1, size(this%allowed)
         if (i > 1) list = list//", "
         list = list//"'"//trim(this%allowed(i))//"'"
      end do
      err_msg = "key '"//this%name//"': '"//v//"' not in allowed set {"//list//"}"
   end subroutine enum_parse

   ! ------------------------------------------------------------------
   ! Real-array key methods.
   ! ------------------------------------------------------------------

   subroutine real_array_parse(this, values, err_msg)
      class(nml_real_array_key_t), intent(inout) :: this
      character(len=*), intent(in) :: values(:)
      character(len=:), allocatable, intent(out) :: err_msg
      real(rk) :: v
      integer :: i, ios
      character(len=128) :: imsg

      if (size(values) > size(this%tgt)) then
         err_msg = "key '"//this%name//"' accepts at most "//fmt_int(size(this%tgt))// &
                   " values, got "//fmt_int(size(values))
         return
      end if
      if (size(values) < 1) then
         err_msg = "key '"//this%name//"' expects at least 1 value"
         return
      end if
      do i = 1, size(values)
         read (values(i), *, iostat=ios, iomsg=imsg) v
         if (ios /= 0) then
            err_msg = "key '"//this%name//"': cannot parse real from '"//trim(values(i))// &
                      "' (element "//fmt_int(i)//")"
            return
         end if
         this%tgt(i) = v
      end do
   end subroutine real_array_parse

   function real_array_value_string(this) result(str)
      class(nml_real_array_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      integer :: i
      str = ""
      do i = 1, size(this%tgt)
         if (i > 1) str = str//", "
         str = str//fmt_real(this%tgt(i))
      end do
   end function real_array_value_string

   function real_array_default_string(this) result(str)
      class(nml_real_array_key_t), intent(in) :: this
      character(len=:), allocatable :: str
      integer :: i
      str = ""
      do i = 1, size(this%default)
         if (i > 1) str = str//", "
         str = str//fmt_real(this%default(i))
      end do
   end function real_array_default_string

   function real_array_is_default(this) result(yes)
      class(nml_real_array_key_t), intent(in) :: this
      logical :: yes
      yes = .false.
      if (size(this%tgt) /= size(this%default)) return
      yes = all(this%tgt == this%default)
   end function real_array_is_default
   ! ------------------------------------------------------------------
   ! Group management.
   ! ------------------------------------------------------------------

   subroutine group_add(this, key)
      !! Add a key to the group (sourced-allocate into a box).  Error
      !! stop on duplicate key name (programming error).
      class(nml_group_t), intent(inout) :: this
      class(nml_key_t), intent(in) :: key
      type(nml_key_box_t), allocatable :: tmp(:)
      integer :: n, i

      if (.not. allocated(this%keys)) allocate (this%keys(0))
      n = size(this%keys)
      do i = 1, n
         if (lower(this%keys(i)%key%name) == lower(key%name)) then
            call global_logger%error("duplicate key '"//key%name//"' in group '"//this%name//"'")
            error stop "rki_nml_schema: duplicate key registration"
         end if
      end do
      allocate (tmp(n + 1))
      do i = 1, n
         call move_alloc(this%keys(i)%key, tmp(i)%key)
      end do
      allocate (tmp(n + 1)%key, source=key)
      call move_alloc(tmp, this%keys)
   end subroutine group_add

   function group_find_key(this, name) result(idx)
      !! Index of `name` in the group's key list, 0 if absent.
      class(nml_group_t), intent(in) :: this
      character(len=*), intent(in) :: name
      integer :: idx, i
      idx = 0
      if (.not. allocated(this%keys)) return
      do i = 1, size(this%keys)
         if (lower(this%keys(i)%key%name) == lower(name)) then
            idx = i
            return
         end if
      end do
   end function group_find_key

   ! ------------------------------------------------------------------
   ! Schema management.
   ! ------------------------------------------------------------------

   subroutine schema_add_group(this, g)
      !! MOVE a locally-built group into the schema — `g` is consumed
      !! (allocatables transferred). Move semantics on purpose: intrinsic
      !! assignment deep-copies the polymorphic key boxes, which NVHPC
      !! miscompiles (heap corruption); move_alloc transfers descriptors
      !! only. Error stop on duplicate name.
      class(nml_schema_t), intent(inout) :: this
      type(nml_group_t), intent(inout) :: g
      type(nml_group_t), allocatable :: tmp(:)
      integer :: n, i

      if (.not. allocated(this%groups)) allocate (this%groups(0))
      n = size(this%groups)
      do i = 1, n
         if (lower(this%groups(i)%name) == lower(g%name)) then
            call global_logger%error("duplicate group '"//g%name//"' in schema")
            error stop "rki_nml_schema: duplicate group registration"
         end if
      end do
      allocate (tmp(n + 1))
      do i = 1, n
         call move_group(this%groups(i), tmp(i))
      end do
      call move_group(g, tmp(n + 1))
      call move_alloc(tmp, this%groups)
   end subroutine schema_add_group

   subroutine move_group(src, dst)
      !! Transfer a group's contents without copying any polymorphic
      !! key box (descriptor moves only — see [[schema_add_group]]).
      type(nml_group_t), intent(inout) :: src
      type(nml_group_t), intent(inout) :: dst
      if (allocated(src%name)) call move_alloc(src%name, dst%name)
      if (allocated(src%doc)) call move_alloc(src%doc, dst%doc)
      if (allocated(src%keys)) call move_alloc(src%keys, dst%keys)
      dst%found = src%found
      dst%cross_check => src%cross_check
      src%cross_check => null()
   end subroutine move_group

   subroutine schema_add_external_group(this, name)
      !! Register a group name the schema knows exists but does NOT
      !! validate (its body is skipped silently).
      class(nml_schema_t), intent(inout) :: this
      character(len=*), intent(in) :: name
      character(len=:), allocatable :: tmp(:)
      integer :: n, newlen, i

      if (.not. allocated(this%external_names)) allocate (character(len=0) :: this%external_names(0))
      n = size(this%external_names)
      newlen = max(len(this%external_names), len(name))
      allocate (character(len=newlen) :: tmp(n + 1))
      do i = 1, n
         tmp(i) = this%external_names(i)
      end do
      tmp(n + 1) = name
      call move_alloc(tmp, this%external_names)
   end subroutine schema_add_external_group

   function schema_find_group(this, name) result(idx)
      !! Index of `name` in the schema's groups, 0 if absent.
      class(nml_schema_t), intent(in) :: this
      character(len=*), intent(in) :: name
      integer :: idx, i
      idx = 0
      if (.not. allocated(this%groups)) return
      do i = 1, size(this%groups)
         if (lower(this%groups(i)%name) == lower(strip_nml(name))) then
            idx = i
            return
         end if
      end do
   end function schema_find_group

   pure function strip_nml(name) result(out)
      !! Drop a trailing `_nml` suffix (namelist files spell groups
      !! `&<name>_nml`; the schema stores the bare `<name>`).
      character(len=*), intent(in) :: name
      character(len=:), allocatable :: out
      integer :: n
      n = len_trim(name)
      if (n > 4) then
         if (lower(name(n - 3:n)) == "_nml") then
            out = name(1:n - 4)
            return
         end if
      end if
      out = trim(name)
   end function strip_nml

   function schema_is_external(this, name) result(yes)
      !! True if `name` is a registered external group.
      class(nml_schema_t), intent(in) :: this
      character(len=*), intent(in) :: name
      logical :: yes
      integer :: i
      yes = .false.
      if (.not. allocated(this%external_names)) return
      do i = 1, size(this%external_names)
         if (lower(trim(this%external_names(i))) == lower(strip_nml(name))) then
            yes = .true.
            return
         end if
      end do
   end function schema_is_external

   ! ------------------------------------------------------------------
   ! Levenshtein + did-you-mean.
   ! ------------------------------------------------------------------

   pure function levenshtein(a, b) result(dist)
      !! Levenshtein edit distance between two (capped) strings.
      character(len=*), intent(in) :: a, b
      integer :: dist
      integer :: la, lb, i, j
      integer, allocatable :: prev(:), cur(:)
      character(len=MAX_SUGGEST_LEN) :: aa, bb

      aa = a
      bb = b
      la = min(len_trim(a), MAX_SUGGEST_LEN)
      lb = min(len_trim(b), MAX_SUGGEST_LEN)
      allocate (prev(0:lb), cur(0:lb))
      do j = 0, lb
         prev(j) = j
      end do
      do i = 1, la
         cur(0) = i
         do j = 1, lb
            if (aa(i:i) == bb(j:j)) then
               cur(j) = prev(j - 1)
            else
               cur(j) = 1 + min(prev(j), cur(j - 1), prev(j - 1))
            end if
         end do
         prev = cur
      end do
      dist = prev(lb)
   end function levenshtein

   function suggest(name, candidates) result(msg)
      !! Build a " (did you mean 'x'?)" fragment if a close candidate
      !! exists (distance <= max(2, len/3)); else empty.
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: candidates(:)
      character(len=:), allocatable :: msg
      integer :: i, d, best_d, best_i, thresh

      msg = ""
      best_d = huge(0)
      best_i = 0
      thresh = max(2, len_trim(name)/3)
      do i = 1, size(candidates)
         d = levenshtein(lower(trim(name)), lower(trim(candidates(i))))
         if (d < best_d) then
            best_d = d
            best_i = i
         end if
      end do
      if (best_i > 0 .and. best_d <= thresh) then
         msg = " (did you mean '"//trim(candidates(best_i))//"'?)"
      end if
   end function suggest
   ! ------------------------------------------------------------------
   ! Documentation writers.
   ! ------------------------------------------------------------------

   subroutine write_key_line(unit, box, short_only, wrote_any)
      !! Write one aligned `key = value  ! [units] doc` line.  When
      !! `short_only`, skip default-valued keys (and leave wrote_any).
      integer, intent(in) :: unit
      type(nml_key_box_t), intent(in) :: box
      logical, intent(in) :: short_only
      logical, intent(inout) :: wrote_any
      character(len=:), allocatable :: lhs, comment, ustr
      logical :: is_def

      is_def = box%key%is_default()
      if (short_only .and. is_def) return
      wrote_any = .true.

      lhs = "  "//box%key%name//" = "//box%key%value_string()
      ustr = ""
      if (len_trim(box%key%units) > 0) ustr = "["//trim(box%key%units)//"] "
      if (is_def) then
         comment = "  ! "//ustr//trim(box%key%doc)//" (default)"
      else
         comment = "  ! "//ustr//trim(box%key%doc)// &
                   " *** non-default (default: "//box%key%default_string()//")"
      end if
      write (unit, "(A)") lhs//comment
   end subroutine write_key_line

   subroutine schema_write_doc_all(this, path)
      !! Write every group and key as valid namelist with aligned docs.
      class(nml_schema_t), intent(in) :: this
      character(len=*), intent(in) :: path
      integer :: u, ig, ik, ios
      character(len=256) :: iomsg
      logical :: dummy

      ! Best-effort documentation dump: an unwritable directory must not abort.
      open (newunit=u, file=path, status="replace", action="write", iostat=ios, iomsg=iomsg)
      if (ios /= 0) return
      do ig = 1, size(this%groups)
         write (u, "(A)") "&"//this%groups(ig)%name//"_nml"
         if (allocated(this%groups(ig)%keys)) then
            do ik = 1, size(this%groups(ig)%keys)
               dummy = .false.
               call write_key_line(u, this%groups(ig)%keys(ik), .false., dummy)
            end do
         end if
         write (u, "(A)") "/"
      end do
      close (u)
   end subroutine schema_write_doc_all

   subroutine schema_write_doc_short(this, path)
      !! Write only non-default keys (MOM6 parameter_doc.short).  Groups
      !! with no non-default keys are omitted entirely.
      class(nml_schema_t), intent(in) :: this
      character(len=*), intent(in) :: path
      integer :: u, ig, ik, ios
      character(len=256) :: iomsg
      logical :: wrote_any

      ! Best-effort documentation dump: an unwritable directory must not abort.
      open (newunit=u, file=path, status="replace", action="write", iostat=ios, iomsg=iomsg)
      if (ios /= 0) return
      do ig = 1, size(this%groups)
         wrote_any = .false.
         if (allocated(this%groups(ig)%keys)) then
            ! First pass: does any non-default key exist in this group?
            do ik = 1, size(this%groups(ig)%keys)
               if (.not. this%groups(ig)%keys(ik)%key%is_default()) wrote_any = .true.
            end do
         end if
         if (.not. wrote_any) cycle
         write (u, "(A)") "&"//this%groups(ig)%name//"_nml"
         do ik = 1, size(this%groups(ig)%keys)
            call write_key_line(u, this%groups(ig)%keys(ik), .true., wrote_any)
         end do
         write (u, "(A)") "/"
      end do
      close (u)
   end subroutine schema_write_doc_short

   subroutine schema_render_markdown(this, path)
      !! Render the schema as Markdown: per group a heading, doc line,
      !! and a knob table.
      class(nml_schema_t), intent(in) :: this
      character(len=*), intent(in) :: path
      integer :: u, ig, ik
      character(len=:), allocatable :: ustr

      open (newunit=u, file=path, status="replace", action="write")
      do ig = 1, size(this%groups)
         ! Blank separator BETWEEN groups only: a trailing blank line after
         ! the last group gets re-trimmed by the end-of-file pre-commit hook
         ! on every regeneration (a recurring commit-abort nuisance).
         if (ig > 1) write (u, "(A)") ""
         write (u, "(A)") "### &"//this%groups(ig)%name//"_nml"
         write (u, "(A)") ""
         if (allocated(this%groups(ig)%doc)) then
            if (len_trim(this%groups(ig)%doc) > 0) then
               write (u, "(A)") trim(this%groups(ig)%doc)
               write (u, "(A)") ""
            end if
         end if
         write (u, "(A)") "| Knob | Default | Units | Description |"
         write (u, "(A)") "|------|---------|-------|-------------|"
         if (allocated(this%groups(ig)%keys)) then
            do ik = 1, size(this%groups(ig)%keys)
               associate (k => this%groups(ig)%keys(ik)%key)
                  ustr = ""
                  if (len_trim(k%units) > 0) ustr = trim(k%units)
                  write (u, "(A)") "| `"//k%name//"` | `"//k%default_string()// &
                     "` | "//ustr//" | "//trim(k%doc)//" |"
               end associate
            end do
         end if
      end do
      close (u)
   end subroutine schema_render_markdown
   ! ------------------------------------------------------------------
   ! Error collection helper.
   ! ------------------------------------------------------------------

   subroutine add_err(errors, n_err, msg)
      !! Append a message to a growable error list.
      character(len=:), allocatable, intent(inout) :: errors(:)
      integer, intent(inout) :: n_err
      character(len=*), intent(in) :: msg
      character(len=:), allocatable :: tmp(:)
      integer :: cap, newlen, i

      if (.not. allocated(errors)) allocate (character(len=256) :: errors(8))
      cap = size(errors)
      if (n_err >= cap .or. len(msg) > len(errors)) then
         newlen = max(len(errors), len(msg))
         allocate (character(len=newlen) :: tmp(max(2*cap, n_err + 1)))
         do i = 1, n_err
            tmp(i) = errors(i)
         end do
         call move_alloc(tmp, errors)
      end if
      n_err = n_err + 1
      errors(n_err) = msg
   end subroutine add_err

   ! ------------------------------------------------------------------
   ! The parser.
   ! ------------------------------------------------------------------

   subroutine schema_parse(this, path, status, errors)
      !! Parse and validate the namelist file at `path`.
      !!
      !! Error-return convention: if BOTH `status` and `errors` are present,
      !! collected messages return in `errors` (unallocated when none) and
      !! `status` = error count (0 = clean), no stop. If either is absent and
      !! any error occurs, each is logged via `global_logger%error` and the
      !! routine `error stop`s.
      class(nml_schema_t), intent(inout) :: this
      character(len=*), intent(in) :: path
      integer, intent(out), optional :: status
      character(len=:), allocatable, intent(out), optional :: errors(:)

      character(len=:), allocatable :: err_buf(:)
      integer :: n_err

      n_err = 0
      call reset_state(this)
      call parse_file(this, path, err_buf, n_err)
      call post_parse_validate(this, path, err_buf, n_err, status, errors)
   end subroutine schema_parse

   subroutine schema_parse_lines(this, lines, n_lines, status, errors)
      !! In-memory sibling of `schema_parse`: validate + apply a namelist
      !! already held as a `lines(:)` character array (no file touch).
      !! Same error-return convention as `schema_parse`.
      class(nml_schema_t), intent(inout) :: this
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      integer, intent(out), optional :: status
      character(len=:), allocatable, intent(out), optional :: errors(:)

      character(len=:), allocatable :: err_buf(:)
      integer :: n_err
      character(len=*), parameter :: LABEL = "<config-string>"

      n_err = 0
      call reset_state(this)
      call parse_from_lines(this, lines, n_lines, LABEL, err_buf, n_err)
      call post_parse_validate(this, LABEL, err_buf, n_err, status, errors)
   end subroutine schema_parse_lines

   subroutine post_parse_validate(this, label, err_buf, n_err, status, errors)
      !! Required-key + cross-check pass shared by the file and buffer
      !! parsers, then the error-return / error-stop finalization.
      type(nml_schema_t), intent(inout) :: this
      character(len=*), intent(in) :: label
      character(len=:), allocatable, intent(inout) :: err_buf(:)
      integer, intent(inout) :: n_err
      integer, intent(out), optional :: status
      character(len=:), allocatable, intent(out), optional :: errors(:)

      integer :: i, ig
      logical :: testable

      testable = present(status) .and. present(errors)

      ! Required-key check (only for groups present in the source).
      do ig = 1, size(this%groups)
         if (.not. this%groups(ig)%found) cycle
         if (.not. allocated(this%groups(ig)%keys)) cycle
         do i = 1, size(this%groups(ig)%keys)
            associate (k => this%groups(ig)%keys(i)%key)
               if (k%required .and. .not. k%found) then
                  call add_err(err_buf, n_err, label//": required key '"//k%name// &
                               "' missing in group '"//this%groups(ig)%name//"'")
               end if
            end associate
         end do
      end do

      ! Cross-checks: run for every registered group regardless of presence.
      do ig = 1, size(this%groups)
         if (associated(this%groups(ig)%cross_check)) then
            block
               character(len=:), allocatable :: cmsg
               cmsg = this%groups(ig)%cross_check()
               if (len_trim(cmsg) > 0) then
                  call add_err(err_buf, n_err, label//": group '"// &
                               this%groups(ig)%name//"': "//trim(cmsg))
               end if
            end block
         end if
      end do

      if (testable) then
         status = n_err
         if (n_err > 0) then
            allocate (character(len=len(err_buf)) :: errors(n_err))
            do i = 1, n_err
               errors(i) = err_buf(i)
            end do
         end if
         return
      end if

      if (n_err > 0) then
         do i = 1, n_err
            call global_logger%error(trim(err_buf(i)))
         end do
         error stop "rki_nml_schema: namelist validation failed"
      end if
   end subroutine post_parse_validate

   subroutine reset_state(this)
      !! Reset parse-state (found flags) before a parse.
      type(nml_schema_t), intent(inout) :: this
      integer :: ig, ik
      do ig = 1, size(this%groups)
         this%groups(ig)%found = .false.
         if (allocated(this%groups(ig)%keys)) then
            do ik = 1, size(this%groups(ig)%keys)
               this%groups(ig)%keys(ik)%key%found = .false.
            end do
         end if
      end do
   end subroutine reset_state
   subroutine read_lines(path, lines, n_lines, ok)
      !! Read the whole file into a deferred-len line array.
      character(len=*), intent(in) :: path
      character(len=:), allocatable, intent(out) :: lines(:)
      integer, intent(out) :: n_lines
      logical, intent(out) :: ok
      integer :: u, ios, cap, maxlen
      character(len=4096) :: buf
      character(len=256) :: imsg
      character(len=:), allocatable :: tmp(:)
      integer :: i

      n_lines = 0
      ok = .false.
      open (newunit=u, file=path, status="old", action="read", iostat=ios, iomsg=imsg)
      if (ios /= 0) return
      ok = .true.
      cap = 64
      maxlen = 256
      allocate (character(len=maxlen) :: lines(cap))
      do
         read (u, "(A)", iostat=ios, iomsg=imsg) buf
         if (ios /= 0) exit
         if (n_lines >= cap .or. len_trim(buf) > maxlen) then
            maxlen = max(maxlen, len_trim(buf))
            allocate (character(len=maxlen) :: tmp(max(2*cap, n_lines + 1)))
            do i = 1, n_lines
               tmp(i) = lines(i)
            end do
            call move_alloc(tmp, lines)
            cap = size(lines)
         end if
         n_lines = n_lines + 1
         lines(n_lines) = buf
      end do
      close (u)
   end subroutine read_lines

   pure function strip_comment(line) result(out)
      !! Remove a trailing `!` comment, respecting single/double quotes.
      character(len=*), intent(in) :: line
      character(len=:), allocatable :: out
      integer :: i
      character :: q
      logical :: in_str

      in_str = .false.
      q = " "
      out = line
      do i = 1, len(line)
         if (in_str) then
            if (line(i:i) == q) in_str = .false.
         else
            if (line(i:i) == '"' .or. line(i:i) == "'") then
               in_str = .true.
               q = line(i:i)
            else if (line(i:i) == "!") then
               out = line(1:i - 1)
               return
            end if
         end if
      end do
   end function strip_comment

   subroutine parse_file(this, path, errors, n_err)
      !! Read the namelist file at `path`, then drive the parse.  The file
      !! is the only on-disk touch; the walk runs on the in-memory lines
      !! (shared with the from-buffer path, parse_from_lines).
      type(nml_schema_t), intent(inout) :: this
      character(len=*), intent(in) :: path
      character(len=:), allocatable, intent(inout) :: errors(:)
      integer, intent(inout) :: n_err

      character(len=:), allocatable :: lines(:)
      integer :: n_lines
      logical :: ok

      call read_lines(path, lines, n_lines, ok)
      if (.not. ok) then
         call add_err(errors, n_err, path//": cannot open file")
         return
      end if
      call parse_from_lines(this, lines, n_lines, path, errors, n_err)
   end subroutine parse_file

   subroutine parse_from_lines(this, lines, n_lines, label, errors, n_err)
      !! Walk an in-memory namelist (already split into `lines`): strip
      !! comments, scan for `&group` headers, dispatch keys.  `label`
      !! names the source in error messages — a file path for parse_file,
      !! a sentinel (e.g. "<config-string>") for buffer parses.
      type(nml_schema_t), intent(inout) :: this
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      character(len=*), intent(in) :: label
      character(len=:), allocatable, intent(inout) :: errors(:)
      integer, intent(inout) :: n_err

      integer :: li, col, gline
      character(len=:), allocatable :: gname, prefix
      logical, allocatable :: group_seen(:)

      if (allocated(this%groups)) then
         allocate (group_seen(size(this%groups)))
         group_seen = .false.
      else
         allocate (group_seen(0))
      end if

      li = 1
      col = 1
      do
         ! Scan for the next '&' that opens a group.
         call next_amp(lines, n_lines, li, col)
         if (li > n_lines) exit
         gline = li
         call read_group_name(lines, n_lines, li, col, gname)
         prefix = make_prefix(label, gline)
         if (len_trim(gname) == 0) then
            call add_err(errors, n_err, prefix//"malformed group header (expected name after '&')")
            ! Skip to next '/' to recover.
            call skip_to_slash(lines, n_lines, li, col)
            cycle
         end if
         call dispatch_group(this, lines, n_lines, li, col, gname, prefix, &
                             group_seen, errors, n_err, label)
      end do
   end subroutine parse_from_lines

   function make_prefix(path, line) result(p)
      !! Build a `path:line: ` error prefix.
      character(len=*), intent(in) :: path
      integer, intent(in) :: line
      character(len=:), allocatable :: p
      p = path//":"//fmt_int(line)//": "
   end function make_prefix

   subroutine next_amp(lines, n_lines, li, col)
      !! Advance (li,col) to the character AFTER the next '&'.  On EOF
      !! sets li = n_lines + 1.  Bare non-blank tokens before '&' are
      !! tolerated as inter-group whitespace.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      integer, intent(inout) :: li, col
      character(len=:), allocatable :: s
      integer :: i

      do while (li <= n_lines)
         s = strip_comment(lines(li))
         do i = col, len(s)
            if (s(i:i) == "&") then
               col = i + 1
               return
            end if
         end do
         li = li + 1
         col = 1
      end do
   end subroutine next_amp

   subroutine read_group_name(lines, n_lines, li, col, gname)
      !! Read an identifier starting at (li,col); advance past it.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      integer, intent(inout) :: li, col
      character(len=:), allocatable, intent(out) :: gname
      character(len=:), allocatable :: s
      integer :: i, start

      gname = ""
      if (li > n_lines) return
      s = strip_comment(lines(li))
      ! Skip leading blanks.
      do while (col <= len(s))
         if (s(col:col) /= " " .and. s(col:col) /= achar(9)) exit
         col = col + 1
      end do
      start = col
      do i = col, len(s)
         if (is_ident_char(s(i:i))) then
            col = i + 1
         else
            exit
         end if
      end do
      if (col > start) gname = s(start:col - 1)
   end subroutine read_group_name

   pure function is_ident_char(c) result(yes)
      !! True for identifier characters (alnum + underscore).
      character, intent(in) :: c
      logical :: yes
      integer :: ic
      ic = iachar(c)
      yes = (ic >= iachar("a") .and. ic <= iachar("z")) .or. &
            (ic >= iachar("A") .and. ic <= iachar("Z")) .or. &
            (ic >= iachar("0") .and. ic <= iachar("9")) .or. &
            (c == "_")
   end function is_ident_char

   subroutine skip_to_slash(lines, n_lines, li, col)
      !! Advance past the next '/' (group terminator); recovery aid.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      integer, intent(inout) :: li, col
      character(len=:), allocatable :: s
      integer :: i
      do while (li <= n_lines)
         s = strip_comment(lines(li))
         do i = col, len(s)
            if (s(i:i) == "/") then
               col = i + 1
               return
            end if
         end do
         li = li + 1
         col = 1
      end do
   end subroutine skip_to_slash
   subroutine dispatch_group(this, lines, n_lines, li, col, gname, prefix, &
                             group_seen, errors, n_err, path)
      !! Resolve a group name to known / external / unknown and parse
      !! (or skip) its body up to the terminating '/'.
      type(nml_schema_t), intent(inout) :: this
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      integer, intent(inout) :: li, col
      character(len=*), intent(in) :: gname, prefix, path
      logical, allocatable, intent(inout) :: group_seen(:)
      character(len=:), allocatable, intent(inout) :: errors(:)
      integer, intent(inout) :: n_err
      integer :: gidx

      gidx = this%find_group(gname)
      if (gidx > 0) then
         if (group_seen(gidx)) then
            call add_err(errors, n_err, prefix//"duplicate occurrence of group '"//gname//"'")
            call skip_to_slash(lines, n_lines, li, col)
            return
         end if
         group_seen(gidx) = .true.
         this%groups(gidx)%found = .true.
         call parse_group_body(this%groups(gidx), lines, n_lines, li, col, &
                               path, errors, n_err)
         return
      end if
      if (this%is_external(gname)) then
         call skip_to_slash(lines, n_lines, li, col)
         return
      end if
      ! Unknown group.
      block
         character(len=:), allocatable :: names(:)
         character(len=:), allocatable :: sug
         names = group_name_list(this)
         sug = ""
         if (size(names) > 0) sug = suggest(strip_nml(gname), names)
         call add_err(errors, n_err, prefix//"unknown group '"//gname//"'"//sug)
      end block
      call skip_to_slash(lines, n_lines, li, col)
   end subroutine dispatch_group

   function group_name_list(this) result(names)
      !! All registered + external group names as a char array.
      type(nml_schema_t), intent(in) :: this
      character(len=:), allocatable :: names(:)
      integer :: ng, ne, i, mlen, k

      ng = 0
      if (allocated(this%groups)) ng = size(this%groups)
      ne = 0
      if (allocated(this%external_names)) ne = size(this%external_names)
      mlen = 1
      do i = 1, ng
         mlen = max(mlen, len(this%groups(i)%name))
      end do
      if (ne > 0) mlen = max(mlen, len(this%external_names))
      allocate (character(len=mlen) :: names(ng + ne))
      k = 0
      do i = 1, ng
         k = k + 1
         names(k) = this%groups(i)%name
      end do
      do i = 1, ne
         k = k + 1
         names(k) = this%external_names(i)
      end do
   end function group_name_list

   function group_key_list(g) result(names)
      !! All key names of a group as a char array.
      type(nml_group_t), intent(in) :: g
      character(len=:), allocatable :: names(:)
      integer :: nk, i, mlen

      nk = 0
      if (allocated(g%keys)) nk = size(g%keys)
      mlen = 1
      do i = 1, nk
         mlen = max(mlen, len(g%keys(i)%key%name))
      end do
      allocate (character(len=mlen) :: names(nk))
      do i = 1, nk
         names(i) = g%keys(i)%key%name
      end do
   end function group_key_list

   subroutine parse_group_body(g, lines, n_lines, li, col, path, errors, n_err)
      !! Tokenize and apply `key = value...` pairs until '/'.
      type(nml_group_t), intent(inout) :: g
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      integer, intent(inout) :: li, col
      character(len=*), intent(in) :: path
      character(len=:), allocatable, intent(inout) :: errors(:)
      integer, intent(inout) :: n_err

      character(len=256), allocatable :: tokens(:)
      integer, allocatable :: tok_line(:)
      integer :: n_tok, ti

      call collect_tokens(lines, n_lines, li, col, path, tokens, tok_line, n_tok, errors, n_err)

      ti = 1
      do while (ti <= n_tok)
         call apply_one_pair(g, tokens, tok_line, n_tok, ti, path, errors, n_err)
      end do
   end subroutine parse_group_body
   subroutine collect_tokens(lines, n_lines, li, col, path, tokens, tok_line, n_tok, errors, n_err)
      !! Tokenize a group body up to and including the terminating '/'.
      !! Token kinds are encoded by their text: "=" , "," , a quoted or
      !! bare value, an identifier, or special markers for unsupported
      !! syntax which are flagged here.  The '/' ends collection.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      integer, intent(inout) :: li, col
      character(len=*), intent(in) :: path
      character(len=256), allocatable, intent(out) :: tokens(:)
      integer, allocatable, intent(out) :: tok_line(:)
      integer, intent(out) :: n_tok
      character(len=:), allocatable, intent(inout) :: errors(:)
      integer, intent(inout) :: n_err

      character(len=:), allocatable :: s
      integer :: i, slen
      logical :: closed
      character :: c, q

      n_tok = 0
      allocate (tokens(32))
      allocate (tok_line(32))
      closed = .false.

      do while (li <= n_lines .and. .not. closed)
         s = strip_comment(lines(li))
         slen = len(s)
         i = col
         do while (i <= slen)
            c = s(i:i)
            if (c == "/") then
               closed = .true.
               col = i + 1
               exit
            end if
            if (c == " " .or. c == achar(9)) then
               i = i + 1
            else if (c == "=") then
               call push_tok(tokens, tok_line, n_tok, "=", li)
               i = i + 1
            else if (c == ",") then
               call push_tok(tokens, tok_line, n_tok, ",", li)
               i = i + 1
            else if (c == '"' .or. c == "'") then
               q = c
               block
                  integer :: j
                  j = i + 1
                  do while (j <= slen)
                     if (s(j:j) == q) exit
                     j = j + 1
                  end do
                  if (j > slen) then
                     call add_err(errors, n_err, make_prefix(path, li)//"unterminated string literal")
                     call push_tok(tokens, tok_line, n_tok, s(i + 1:slen), li)
                     i = slen + 1
                  else
                     call push_tok(tokens, tok_line, n_tok, s(i + 1:j - 1), li)
                     i = j + 1
                  end if
               end block
            else
               ! Bare token: identifier or value, runs to delimiter.
               block
                  integer :: start, j
                  character(len=:), allocatable :: word
                  start = i
                  j = i
                  do while (j <= slen)
                     if (s(j:j) == " " .or. s(j:j) == achar(9) .or. s(j:j) == "=" .or. &
                         s(j:j) == "," .or. s(j:j) == "/") exit
                     j = j + 1
                  end do
                  word = s(start:j - 1)
                  call check_unsupported(word, path, li, errors, n_err)
                  call push_tok(tokens, tok_line, n_tok, word, li)
                  i = j
               end block
            end if
         end do
         if (.not. closed) then
            li = li + 1
            col = 1
         end if
      end do

      if (.not. closed) then
         call add_err(errors, n_err, make_prefix(path, max(li - 1, 1))// &
                      "group not terminated with '/' before end of file")
      end if
   end subroutine collect_tokens

   subroutine check_unsupported(word, path, line, errors, n_err)
      !! Flag unsupported namelist syntax embedded in a bare token.
      character(len=*), intent(in) :: word
      character(len=*), intent(in) :: path
      integer, intent(in) :: line
      character(len=:), allocatable, intent(inout) :: errors(:)
      integer, intent(inout) :: n_err

      if (index(word, "%") > 0) then
         call add_err(errors, n_err, make_prefix(path, line)// &
                      "derived-type reference '"//trim(word)//"' is not supported")
      else if (index(word, "(") > 0) then
         call add_err(errors, n_err, make_prefix(path, line)// &
                      "indexed assignment '"//trim(word)//"' is not supported")
      else if (is_repeat_count(word)) then
         call add_err(errors, n_err, make_prefix(path, line)// &
                      "repeat-count syntax '"//trim(word)//"' is not supported")
      end if
   end subroutine check_unsupported

   pure function is_repeat_count(word) result(yes)
      !! True if `word` looks like `N*value` (integer, then '*').
      character(len=*), intent(in) :: word
      logical :: yes
      integer :: star, i
      character :: c
      yes = .false.
      star = index(word, "*")
      if (star <= 1) return
      do i = 1, star - 1
         c = word(i:i)
         if (c < "0" .or. c > "9") return
      end do
      yes = .true.
   end function is_repeat_count

   subroutine push_tok(tokens, tok_line, n_tok, val, line)
      !! Append a token, growing the buffer as needed.
      character(len=256), allocatable, intent(inout) :: tokens(:)
      integer, allocatable, intent(inout) :: tok_line(:)
      integer, intent(inout) :: n_tok
      character(len=*), intent(in) :: val
      integer, intent(in) :: line
      character(len=256), allocatable :: tmp(:)
      integer, allocatable :: tmpl(:)
      integer :: cap, i

      cap = size(tokens)
      if (n_tok >= cap) then
         allocate (tmp(2*cap), tmpl(2*cap))
         do i = 1, n_tok
            tmp(i) = tokens(i)
            tmpl(i) = tok_line(i)
         end do
         call move_alloc(tmp, tokens)
         call move_alloc(tmpl, tok_line)
      end if
      n_tok = n_tok + 1
      tokens(n_tok) = val
      tok_line(n_tok) = line
   end subroutine push_tok
   subroutine apply_one_pair(g, tokens, tok_line, n_tok, ti, path, errors, n_err)
      !! Consume one `key = value [value...]` starting at token `ti`;
      !! advance `ti` past the consumed tokens.  Validates structure,
      !! unknown keys, duplicate keys, and applies via parse_tokens.
      type(nml_group_t), intent(inout) :: g
      character(len=256), intent(in) :: tokens(:)
      integer, intent(in) :: tok_line(:)
      integer, intent(in) :: n_tok
      integer, intent(inout) :: ti
      character(len=*), intent(in) :: path
      character(len=:), allocatable, intent(inout) :: errors(:)
      integer, intent(inout) :: n_err

      character(len=:), allocatable :: keyname, prefix
      integer :: kidx, kline, vstart, vcount, j
      character(len=256), allocatable :: vals(:)

      ! Skip stray separators.
      do while (ti <= n_tok)
         if (trim(tokens(ti)) == ",") then
            ti = ti + 1
         else
            exit
         end if
      end do
      if (ti > n_tok) return

      keyname = trim(tokens(ti))
      kline = tok_line(ti)
      prefix = make_prefix(path, kline)

      if (keyname == "=") then
         call add_err(errors, n_err, prefix//"stray '=' with no key name")
         ti = ti + 1
         return
      end if

      ! Expect '=' next.
      if (ti + 1 > n_tok) then
         call add_err(errors, n_err, prefix//"key '"//keyname//"' has no '=' assignment")
         ti = ti + 1
         return
      end if
      if (trim(tokens(ti + 1)) /= "=") then
         call add_err(errors, n_err, prefix//"expected '=' after key '"//keyname// &
                      "', got '"//trim(tokens(ti + 1))//"'")
         ti = ti + 1
         return
      end if

      ! Gather value tokens (skip commas) until next `ident =` or end.
      vstart = ti + 2
      vcount = 0
      allocate (vals(0))
      j = vstart
      do while (j <= n_tok)
         if (trim(tokens(j)) == ",") then
            j = j + 1
            cycle
         end if
         ! Lookahead: is this token a key (followed by '=')?
         if (j + 1 <= n_tok) then
            if (trim(tokens(j + 1)) == "=") exit
         end if
         call append_val(vals, vcount, tokens(j))
         j = j + 1
      end do

      ti = j

      ! Resolve and apply.
      kidx = g%find_key(keyname)
      if (kidx == 0) then
         block
            character(len=:), allocatable :: names(:)
            character(len=:), allocatable :: sug
            names = group_key_list(g)
            sug = ""
            if (size(names) > 0) sug = suggest(keyname, names)
            call add_err(errors, n_err, prefix//"unknown key '"//keyname// &
                         "' in group '"//g%name//"'"//sug)
         end block
         return
      end if
      if (g%keys(kidx)%key%found) then
         call add_err(errors, n_err, prefix//"duplicate key '"//keyname// &
                      "' in group '"//g%name//"'")
         return
      end if
      g%keys(kidx)%key%found = .true.
      if (vcount == 0) then
         call add_err(errors, n_err, prefix//"key '"//keyname//"' has no value")
         return
      end if
      block
         character(len=:), allocatable :: emsg
         call g%keys(kidx)%key%parse_tokens(vals(1:vcount), emsg)
         if (allocated(emsg)) then
            call add_err(errors, n_err, prefix//trim(emsg))
         end if
      end block
   end subroutine apply_one_pair

   subroutine append_val(vals, vcount, tok)
      !! Append one value token to the values buffer.
      character(len=256), allocatable, intent(inout) :: vals(:)
      integer, intent(inout) :: vcount
      character(len=*), intent(in) :: tok
      character(len=256), allocatable :: tmp(:)
      integer :: i

      if (vcount >= size(vals)) then
         allocate (tmp(max(4, 2*size(vals))))
         do i = 1, vcount
            tmp(i) = vals(i)
         end do
         call move_alloc(tmp, vals)
      end if
      vcount = vcount + 1
      vals(vcount) = tok
   end subroutine append_val

end module rki_nml_schema
