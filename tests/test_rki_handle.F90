!! Unit tests for the opaque-handle FFI infrastructure.
!!
!! Phase 1 of the FFI handle refactor:
!! verify create/check/destroy semantics and that the magic-number
!! sentinel rejects null, garbage, and stale handles for both
!! structured and unstructured handle types.
module test_rki_handle
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: iso_c_binding, only: c_associated, c_loc, c_null_ptr, &
                                                                             c_ptr
   use rki_handle, only: solver_handle_t, solver_handle_unstr_t, &
                         handle_create, handle_destroy, handle_check, &
                         handle_create_unstr, handle_destroy_unstr, &
                         handle_check_unstr, &
                         HANDLE_OK, HANDLE_ERR_NULL, HANDLE_ERR_INVALID
   implicit none
   private

   public :: collect_handle_tests

contains

   subroutine collect_handle_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("structured_create_check_destroy", &
                               test_struct_lifecycle), &
                  new_unittest("structured_check_rejects_null", &
                               test_struct_null), &
                  new_unittest("structured_check_rejects_destroyed", &
                               test_struct_destroyed), &
                  new_unittest("structured_check_rejects_unstr", &
                               test_struct_rejects_unstr), &
                  new_unittest("structured_destroy_idempotent", &
                               test_struct_destroy_idempotent), &
                  new_unittest("unstr_create_check_destroy", &
                               test_unstr_lifecycle), &
                  new_unittest("unstr_check_rejects_null", &
                               test_unstr_null), &
                  new_unittest("unstr_check_rejects_destroyed", &
                               test_unstr_destroyed), &
                  new_unittest("unstr_check_rejects_struct", &
                               test_unstr_rejects_struct), &
                  new_unittest("two_handles_coexist", &
                               test_two_handles_coexist) &
                  ]
   end subroutine collect_handle_tests

   ! ================================================================
   ! Structured handle
   ! ================================================================

   subroutine test_struct_lifecycle(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: c_h
      type(solver_handle_t), pointer :: h
      integer :: status

      c_h = handle_create()
      call check(error, c_associated(c_h), "create returns non-null")
      if (allocated(error)) return

      status = handle_check(c_h, h)
      call check(error, status == HANDLE_OK, "check accepts fresh handle")
      if (allocated(error)) return
      call check(error, associated(h), "check sets pointer")
      if (allocated(error)) return
      call check(error, associated(h%state), &
                 "create allocates state pointer component")
      if (allocated(error)) return

      call handle_destroy(c_h)
      call check(error,.not. c_associated(c_h), "destroy nulls c_ptr")
   end subroutine test_struct_lifecycle

   subroutine test_struct_null(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: c_h
      type(solver_handle_t), pointer :: h
      integer :: status

      c_h = c_null_ptr
      status = handle_check(c_h, h)
      call check(error, status == HANDLE_ERR_NULL, &
                 "null c_ptr is rejected with HANDLE_ERR_NULL")
      if (allocated(error)) return
      call check(error,.not. associated(h), "null check leaves pointer null")
   end subroutine test_struct_null

   subroutine test_struct_destroyed(error)
      !! After destroy, the c_ptr is nulled — but a stale copy held by
      !! the caller would still alias the freed memory. handle_check
      !! must reject the stale copy via the magic-number sentinel.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: c_h, stale
      type(solver_handle_t), pointer :: h
      integer :: status

      c_h = handle_create()
      stale = c_h            !! copy before destroy
      call handle_destroy(c_h)

      status = handle_check(stale, h)
      !! Memory was freed and magic was overwritten with MAGIC_DESTROYED
      !! before deallocate. The dereferenced magic is now either the
      !! sentinel or random freed-heap content — either way, not
      !! MAGIC_STRUCTURED, so check returns HANDLE_ERR_INVALID.
      call check(error, status == HANDLE_ERR_INVALID, &
                 "destroyed-handle copy is rejected as invalid")
      if (allocated(error)) return
      call check(error,.not. associated(h), &
                 "rejected check leaves pointer null")
   end subroutine test_struct_destroyed

   subroutine test_struct_rejects_unstr(error)
      !! An unstructured handle passed to the structured check must
      !! be rejected — the magic numbers are distinct.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: c_h_unstr
      type(solver_handle_t), pointer :: h
      integer :: status

      c_h_unstr = handle_create_unstr()
      status = handle_check(c_h_unstr, h)
      call check(error, status == HANDLE_ERR_INVALID, &
                 "unstr handle rejected by structured check")
      call handle_destroy_unstr(c_h_unstr)
   end subroutine test_struct_rejects_unstr

   subroutine test_struct_destroy_idempotent(error)
      !! handle_destroy on a null or already-destroyed c_ptr must not
      !! crash — Python __del__ may call it more than once.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: c_h

      c_h = c_null_ptr
      call handle_destroy(c_h)        !! null: must no-op
      call check(error,.not. c_associated(c_h), "null destroy leaves c_ptr null")
      if (allocated(error)) return

      c_h = handle_create()
      call handle_destroy(c_h)        !! first destroy
      call handle_destroy(c_h)        !! second destroy on now-null c_h
      call check(error,.not. c_associated(c_h), "double destroy is safe")
   end subroutine test_struct_destroy_idempotent

   ! ================================================================
   ! Unstructured handle (mirrors structured tests)
   ! ================================================================

   subroutine test_unstr_lifecycle(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: c_h
      type(solver_handle_unstr_t), pointer :: h
      integer :: status

      c_h = handle_create_unstr()
      call check(error, c_associated(c_h), "create_unstr returns non-null")
      if (allocated(error)) return

      status = handle_check_unstr(c_h, h)
      call check(error, status == HANDLE_OK, "check_unstr accepts fresh handle")
      if (allocated(error)) return
      call check(error, associated(h%state), "state pointer allocated")
      if (allocated(error)) return
      call check(error, associated(h%mesh), "mesh pointer allocated")
      if (allocated(error)) return

      call handle_destroy_unstr(c_h)
      call check(error,.not. c_associated(c_h), "destroy_unstr nulls c_ptr")
   end subroutine test_unstr_lifecycle

   subroutine test_unstr_null(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: c_h
      type(solver_handle_unstr_t), pointer :: h
      integer :: status

      c_h = c_null_ptr
      status = handle_check_unstr(c_h, h)
      call check(error, status == HANDLE_ERR_NULL, &
                 "null c_ptr rejected by unstr check")
   end subroutine test_unstr_null

   subroutine test_unstr_destroyed(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: c_h, stale
      type(solver_handle_unstr_t), pointer :: h
      integer :: status

      c_h = handle_create_unstr()
      stale = c_h
      call handle_destroy_unstr(c_h)

      status = handle_check_unstr(stale, h)
      call check(error, status == HANDLE_ERR_INVALID, &
                 "destroyed unstr handle rejected")
   end subroutine test_unstr_destroyed

   subroutine test_unstr_rejects_struct(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: c_h_struct
      type(solver_handle_unstr_t), pointer :: h
      integer :: status

      c_h_struct = handle_create()
      status = handle_check_unstr(c_h_struct, h)
      call check(error, status == HANDLE_ERR_INVALID, &
                 "structured handle rejected by unstr check")
      call handle_destroy(c_h_struct)
   end subroutine test_unstr_rejects_struct

   subroutine test_two_handles_coexist(error)
      !! The whole point of the refactor: two handles in one process,
      !! distinguishable, independently destroyable.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: c_s, c_u
      type(solver_handle_t), pointer :: hs
      type(solver_handle_unstr_t), pointer :: hu
      integer :: status

      c_s = handle_create()
      c_u = handle_create_unstr()

      call check(error, c_associated(c_s) .and. c_associated(c_u), &
                 "both handles allocated")
      if (allocated(error)) return

      status = handle_check(c_s, hs)
      call check(error, status == HANDLE_OK, "structured side validates")
      if (allocated(error)) return

      status = handle_check_unstr(c_u, hu)
      call check(error, status == HANDLE_OK, "unstructured side validates")
      if (allocated(error)) return

      !! Destroy in reverse order; structured handle must still validate
      !! after the unstructured one is gone.
      call handle_destroy_unstr(c_u)
      status = handle_check(c_s, hs)
      call check(error, status == HANDLE_OK, &
                 "structured handle survives unstr destroy")
      if (allocated(error)) return

      call handle_destroy(c_s)
   end subroutine test_two_handles_coexist

end module test_rki_handle
