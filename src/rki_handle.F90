!! Opaque handle types for the C-compatible API.
!!
!! Replaces the module-level singletons in `rki_api` and `rki_api_unstr`
!! with heap-allocated handles, so multiple simulations can coexist in one
!! process (e.g. structured + unstructured at once, multi-grid coupling).
!!
!! Layout: each handle owns the data + state-machine flags that used to
!! live as `save` variables at module scope. C/Python callers see only an
!! opaque `c_ptr`; `handle_check_*` resolves it to a Fortran pointer and
!! verifies the magic-number sentinel before any access.
!!
!! Phase 1 (this module): types + create/destroy/check helpers only.
!! Phase 2/3: the existing API entrypoints stop touching module-static
!! state and instead route through `h => handle_check_*(c_handle)`.
module rki_handle
   use, intrinsic :: iso_c_binding, only: c_associated, c_f_pointer, c_loc, &
                                                                             c_null_ptr, c_ptr
   use rki_constants, only: wp
   use rki_config, only: config_t
   use rki_state, only: state_t
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_state_unstr, only: unstr_state_t
   use rki_output_unstr, only: output_unstr_state_t
   use rki_gauges_unstr, only: gauge_unstr_t
   implicit none
   private

   public :: solver_handle_t, solver_handle_unstr_t
   public :: handle_create, handle_destroy, handle_check
   public :: handle_create_unstr, handle_destroy_unstr, handle_check_unstr
   public :: HANDLE_OK, HANDLE_ERR_NULL, HANDLE_ERR_INVALID

   !! Magic-number sentinels. Distinct values per backend so a structured
   !! handle passed to an unstructured entrypoint (or vice versa) is rejected
   !! with the same code path as a fully bogus pointer. The values are
   !! arbitrary; they just need to be unlikely to occur in uninitialised
   !! memory and distinct from each other.
   integer, parameter :: MAGIC_STRUCTURED = int(z'72616B53')   !! 'rakS'
   integer, parameter :: MAGIC_UNSTRUCTURED = int(z'72616B55')   !! 'rakU'
   integer, parameter :: MAGIC_DESTROYED = int(z'DEADBEEF')

   !! handle_check return codes
   integer, parameter :: HANDLE_OK = 0
   integer, parameter :: HANDLE_ERR_NULL = 1
   integer, parameter :: HANDLE_ERR_INVALID = 2

   type :: solver_handle_t
      !! Structured-solver handle. Mirrors the module-static state in
      !! `rki_api` so the Phase 2 cutover is a mechanical rename of
      !! bare `state`/`cfg`/flags to `h%state`/`h%cfg`/etc.
      integer :: magic = 0
      type(config_t)        :: cfg
      type(state_t), pointer :: state => null()
         !! Pointer (not value) so `c_loc(h%state%barotropic%h)` keeps a stable
         !! address across handle moves. Allocated in handle_create.
      real(wp) :: t_current = 0.0_wp
      logical  :: is_initialised = .false.
      logical  :: is_finalised = .false.
      logical  :: on_device = .false.
      logical  :: host_is_current = .false.
      logical  :: user_set_bathymetry = .false.
      logical  :: user_set_h = .false.
      logical  :: user_set_layers = .false.
   end type solver_handle_t

   type :: solver_handle_unstr_t
      !! Unstructured-solver handle. Mirrors the module-static state in
      !! `rki_api_unstr`.
      integer :: magic = 0
      type(mesh_unstr_t), pointer :: mesh => null()
      type(unstr_state_t), pointer :: state => null()
      type(output_unstr_state_t)         :: ostate
      type(gauge_unstr_t)                :: gauges
      logical :: is_initialised = .false.
      logical :: is_finalised = .false.
      logical :: on_device = .false.
      logical :: output_open = .false.
      logical :: gauges_open = .false.
   end type solver_handle_unstr_t

contains

   ! ================================================================
   ! Structured handle lifecycle
   ! ================================================================

   function handle_create() result(c_handle)
      !! Allocate a new structured handle and return its opaque C pointer.
      !! Caller (typically a `bind(C)` create entrypoint in `rki_api`) is
      !! responsible for eventually passing it to `handle_destroy`.
      type(c_ptr) :: c_handle
      type(solver_handle_t), pointer :: h

      allocate (h)
      allocate (h%state)
      h%magic = MAGIC_STRUCTURED
      c_handle = c_loc(h)
   end function handle_create

   subroutine handle_destroy(c_handle)
      !! Free a structured handle. Idempotent on null/already-destroyed
      !! handles — silently no-ops rather than aborting, so Python
      !! `__del__` can call this without first checking state.
      type(c_ptr), intent(inout) :: c_handle
      type(solver_handle_t), pointer :: h

      if (.not. c_associated(c_handle)) return
      call c_f_pointer(c_handle, h)
      if (h%magic /= MAGIC_STRUCTURED) return    !! garbage or already destroyed
      h%magic = MAGIC_DESTROYED
      if (associated(h%state)) deallocate (h%state)
      deallocate (h)
      c_handle = c_null_ptr
   end subroutine handle_destroy

   function handle_check(c_handle, h) result(status)
      !! Resolve an opaque `c_ptr` to a Fortran pointer and verify it.
      !!
      !! On success: `h` is associated and `status == HANDLE_OK`.
      !! On failure: `h` is null and `status` is one of HANDLE_ERR_*.
      !! Callers in `rki_api` should propagate the status to C as a
      !! non-zero return code (mirroring existing error-code conventions).
      type(c_ptr), intent(in)                  :: c_handle
      type(solver_handle_t), pointer, intent(out) :: h
      integer :: status

      h => null()
      if (.not. c_associated(c_handle)) then
         status = HANDLE_ERR_NULL
         return
      end if
      call c_f_pointer(c_handle, h)
      if (h%magic /= MAGIC_STRUCTURED) then
         h => null()
         status = HANDLE_ERR_INVALID
         return
      end if
      status = HANDLE_OK
   end function handle_check

   ! ================================================================
   ! Unstructured handle lifecycle
   ! ================================================================

   function handle_create_unstr() result(c_handle)
      !! Allocate a new unstructured handle.
      type(c_ptr) :: c_handle
      type(solver_handle_unstr_t), pointer :: h

      allocate (h)
      allocate (h%mesh)
      allocate (h%state)
      h%magic = MAGIC_UNSTRUCTURED
      c_handle = c_loc(h)
   end function handle_create_unstr

   subroutine handle_destroy_unstr(c_handle)
      !! Free an unstructured handle. Idempotent on null/destroyed.
      type(c_ptr), intent(inout) :: c_handle
      type(solver_handle_unstr_t), pointer :: h

      if (.not. c_associated(c_handle)) return
      call c_f_pointer(c_handle, h)
      if (h%magic /= MAGIC_UNSTRUCTURED) return
      h%magic = MAGIC_DESTROYED
      if (associated(h%state)) deallocate (h%state)
      if (associated(h%mesh)) deallocate (h%mesh)
      deallocate (h)
      c_handle = c_null_ptr
   end subroutine handle_destroy_unstr

   function handle_check_unstr(c_handle, h) result(status)
      !! Resolve and verify an unstructured handle.
      type(c_ptr), intent(in)                        :: c_handle
      type(solver_handle_unstr_t), pointer, intent(out) :: h
      integer :: status

      h => null()
      if (.not. c_associated(c_handle)) then
         status = HANDLE_ERR_NULL
         return
      end if
      call c_f_pointer(c_handle, h)
      if (h%magic /= MAGIC_UNSTRUCTURED) then
         h => null()
         status = HANDLE_ERR_INVALID
         return
      end if
      status = HANDLE_OK
   end function handle_check_unstr

end module rki_handle
