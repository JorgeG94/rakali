!! Reusable 3D scratch-buffer type.
module rki_scratch_3d
   !! Single owner for `(n1, n2, n3)` scratch arrays.  Every column-
   !! local kernel slot (continuity, coriolis_adv, pressure_force,
   !! lateral-mix vorticity scratch, KPP solve scratch) declares its
   !! per-step workspace as a `scratch_3d_buffer_t` rather than as a
   !! bare `real(wp), allocatable :: foo(:, :, :)`.
   !!
   !! Why centralise this:
   !!   1. **One allocation policy.**  The strided hybrid loop
   !!      pattern (see `src/core/ocean/README.md`) sizes every
   !!      column-local scratch at `(stride, stride, nz_ml)` with
   !!      `stride = ocean_state%column_stride`.  Putting that policy
   !!      in one type means a new slot follows the convention by
   !!      construction, not by remembering to copy-paste it.
   !!   2. **One GPU-mapping site.**  `enter_data` / `exit_data`
   !!      live on the type, so every slot that uses this type
   !!      inherits the mapping for free — no per-slot directive
   !!      copies to keep in sync.
   !!   3. **One resize site.**  If a kernel discovers it needs a
   !!      different shape mid-run (rare — only happens if `stride`
   !!      becomes dynamic), `ensure_size` reallocates and remaps in
   !!      one place.
   !!   4. **Cheap finite-check.**  Phase 6 debug guards walk every
   !!      `scratch_3d_buffer_t` instance after each slot kernel and
   !!      verify no NaN / out-of-range leaked into scratch.
   !!
   !! Phase 2 status: alloc + GPU enter_data / exit_data live.  Phase
   !! 1+ items (`ensure_size`, `assert_finite`) still pending.
   use rki_constants, only: wp
   use, intrinsic :: iso_fortran_env, only: int64
   use rki_mem_report, only: arr_bytes
   implicit none
   private

   public :: scratch_3d_buffer_t
   ! Non-polymorphic device-attach entry points.  Owning slots call these
   ! directly (with a `type(scratch_3d_buffer_t)` actual) so the OpenMP
   ! map base is the heap object, not a polymorphic stack box (the AMD
   ! libomptarget cross-slot-overlap fix).  The type-bound `enter_data` /
   ! `exit_data` are thin select-type wrappers over these.
   public :: scratch_3d_buffer_enter_data_impl
   public :: scratch_3d_buffer_exit_data_impl

   type :: scratch_3d_buffer_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Prefer this to
         !! `allocated(this%data)` — tracks GPU device attachment too.

      ! ---- Declared shape ----
      integer :: n1 = 0
         !! First-dim size, typically `stride` (or `stride+1` on
         !! face-located buffers).
      integer :: n2 = 0
         !! Second-dim size, typically `stride`.
      integer :: n3 = 0
         !! Third-dim size, typically `nz_ml` (or `nz_ml+1` for
         !! layer-interface buffers).

      ! ---- Payload ----
      real(wp), allocatable :: data(:, :, :)
         !! The scratch storage.  Access in a kernel as
         !! `slot%scratch%data(i, j, k)`.  Note: this is a *single*
         !! derived type (not array-of-derived-types), so the
         !! `!$acc` indirection bug from the array-of-DT pattern
         !! does not apply.

      ! ---- Optional name ----
      character(len=32) :: name = ""
         !! Human-readable tag for diagnostics + debug guard
         !! messages ("kpp_bulk_ri_scratch", "ppm_h_face_left_x", ...).
   contains
      procedure, non_overridable :: init => scratch_3d_buffer_init
      procedure, non_overridable :: destroy => scratch_3d_buffer_destroy
      procedure, non_overridable :: enter_data => scratch_3d_buffer_enter_data
      procedure, non_overridable :: exit_data => scratch_3d_buffer_exit_data
      ! Phase 6+: procedure :: assert_finite — debug-mode invariant
      procedure, non_overridable :: bytes => scratch_3d_buffer_bytes
   end type scratch_3d_buffer_t

contains

   subroutine scratch_3d_buffer_init(this, n1, n2, n3, name)
      !! Allocate the host-side storage at (n1, n2, n3), zero-filled.
      !! GPU attachment is a separate step via `enter_data` —
      !! `init` runs before the device exists in the typical
      !! init-then-enter-data sequence.
      class(scratch_3d_buffer_t), intent(inout) :: this
      integer, intent(in) :: n1, n2, n3
      character(len=*), intent(in), optional :: name
      this%n1 = n1
      this%n2 = n2
      this%n3 = n3
      if (present(name)) this%name = name
      allocate (this%data(n1, n2, n3), source=0.0_wp)
      this%is_init = .true.
   end subroutine scratch_3d_buffer_init

   subroutine scratch_3d_buffer_destroy(this)
      class(scratch_3d_buffer_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%data)) deallocate (this%data)
      this%n1 = 0
      this%n2 = 0
      this%n3 = 0
      this%name = ""
   end subroutine scratch_3d_buffer_destroy

   subroutine scratch_3d_buffer_enter_data(this)
      !! Type-bound wrapper — keeps `buf%enter_data()` call sites working.
      !! Delegates to the non-polymorphic impl so the OpenMP map base is
      !! the heap object (see the impl + the public-decl comment above).
      class(scratch_3d_buffer_t), intent(inout) :: this
      select type (this)
      type is (scratch_3d_buffer_t)
         call scratch_3d_buffer_enter_data_impl(this)
      end select
   end subroutine scratch_3d_buffer_enter_data

   subroutine scratch_3d_buffer_enter_data_impl(this)
      !! Attach the scratch payload to the device.  `type(...)` (not
      !! `class`) dummy on purpose: a by-reference non-polymorphic dummy
      !! aliases the heap object, so `create(this%data)` attaches against
      !! a heap base — no polymorphic stack box for AMD to reject.
      !!
      !! No-op when the payload was never allocated.  Slots may GATE a
      !! buffer's `init` on a runtime knob (see
      !! `ocean_pressure_force_t%scratch_gated`) and still call this
      !! unconditionally from their `enter_data` walk; mapping an
      !! unallocated allocatable is not defined behaviour, so the guard
      !! lives here rather than at each call site.
      type(scratch_3d_buffer_t), intent(inout) :: this
      if (.not. allocated(this%data)) return
      ! [acc->omp-cpu] dropped device directive (host no-op): enter data
   end subroutine scratch_3d_buffer_enter_data_impl

   subroutine scratch_3d_buffer_exit_data(this)
      !! Type-bound wrapper — see scratch_3d_buffer_enter_data.
      class(scratch_3d_buffer_t), intent(inout) :: this
      select type (this)
      type is (scratch_3d_buffer_t)
         call scratch_3d_buffer_exit_data_impl(this)
      end select
   end subroutine scratch_3d_buffer_exit_data

   subroutine scratch_3d_buffer_exit_data_impl(this)
      !! Release the scratch payload from the device.  No `copyout`
      !! — scratch contents are per-step intermediates with no
      !! host-side meaning.  No-op when never allocated (mirrors
      !! `scratch_3d_buffer_enter_data_impl`'s gate).
      type(scratch_3d_buffer_t), intent(inout) :: this
      if (.not. allocated(this%data)) return
      ! [acc->omp-cpu] dropped device directive (host no-op): exit data
   end subroutine scratch_3d_buffer_exit_data_impl

   pure function scratch_3d_buffer_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the 3D scratch buffer slot (0 when
      !! unallocated).
      class(scratch_3d_buffer_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%data)
   end function scratch_3d_buffer_bytes

end module rki_scratch_3d
