!! C-compatible API for the Rakali shallow water solver
!!
!! Provides bind(C) functions for Python/Julia/C callers via a shared library.
!! Multi-instance: each call to rki_create returns an opaque `void*` handle
!! that identifies one simulation. Multiple handles may coexist in the same
!! process (the historical singleton restriction is gone).
!!
!! Lifecycle (one-shot):
!!   void* h;
!!   rki_create_from_namelist(filename, len, &h)  -> allocate + init + GPU
!!   rki_step(h, n) / rki_advance(h, target_time)
!!   rki_get_h_ptr(h, ptr, n), etc.
!!   rki_destroy(&h)                              -> clean up; nulls h
!!
!! Lifecycle (two-phase, for programmatic bathymetry/IC):
!!   rki_create(filename, len, &h)
!!   rki_set_bathymetry(h, data, nx, ny)
!!   rki_set_initial_h(h, data, nx, ny)
!!   rki_finalize(h)
!!   rki_step / rki_advance / ...
!!
!! Error codes:
!!   0     success
!!   -1    not initialised (or other entry-specific precondition failed)
!!   -2    grid mismatch / already finalised / out-of-range argument
!!   -3    multilayer mode not active / restart grid mismatch
!!  -100   invalid handle (null or stale c_ptr)
module rki_api
   use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, &
                                                                             c_loc, c_null_ptr, c_ptr
   use rki_constants, only: wp
   use rki_config, only: read_config, read_config_from_string, &
                         validate_config, MAX_TIDAL_CONSTITUENTS
   use rki_config_schema, only: nml_dirname
   use rki_nml_schema, only: nml_schema_t
   use rki_state, only: load_bathymetry, fill_bathymetry_ghosts
   use rki_boundary_types, only: bc_init_from_config, bc_init_chapman
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data, &
                         solver_update_host, solver_update_forcing, &
                         solver_update_state
   use rki_init, only: init_from_config, init_multilayer_from_barotropic
   use rki_output, only: restart_write, restart_read
   use rki_handle, only: solver_handle_t, handle_create, handle_destroy, &
                         handle_check, HANDLE_OK
   use pic_logger, only: logger => global_logger, info_level
   use pic_strings, only: to_string
   use rki_banner, only: print_banner
   use pic_knowledge, only: get_knowledge
   implicit none
   private

   integer(c_int), parameter :: RKI_ERR_BAD_HANDLE = -100

   !! State machine for GPU data lifecycle (h%on_device, h%host_is_current):
   !!   (F, F) — initial / after destroy
   !!   (T, F) — after finalize/step/advance (GPU canonical)
   !!   (T, T) — after rki_sync_host (both valid)
   !!   (F, T) — after ensure_host/exit_data (host only)
   !!
   !! Transitions:
   !!   rki_finalize      -> (T, F)
   !!   rki_step/advance  -> (T, F)
   !!   rki_sync_host     -> (T, T)
   !!   ensure_host       -> if host_is_current: no-op; else (F, T)
   !!   ensure_device     -> (T, F)
   !!   rki_destroy       -> handle freed

contains

   ! ================================================================
   ! Helpers
   ! ================================================================

   subroutine c_to_f_string(c_str, c_len, f_str)
      !! Convert C string (char* + length) to Fortran character.
      integer(c_int), intent(in), value :: c_len
      character(kind=c_char), intent(in) :: c_str(c_len)
      character(len=:), allocatable, intent(out) :: f_str
      integer :: i

      allocate (character(len=c_len) :: f_str)
      do i = 1, c_len
         f_str(i:i) = c_str(i)
      end do
   end subroutine c_to_f_string

   function resolve(c_handle, h) result(ierr)
      !! Resolve an opaque c_ptr to a Fortran handle pointer.
      !! Returns 0 on success, RKI_ERR_BAD_HANDLE if the c_ptr is null
      !! or fails the magic-number check. Callers propagate the code.
      type(c_ptr), intent(in) :: c_handle
      type(solver_handle_t), pointer, intent(out) :: h
      integer(c_int) :: ierr
      integer :: status

      status = handle_check(c_handle, h)
      if (status == HANDLE_OK) then
         ierr = 0
      else
         ierr = RKI_ERR_BAD_HANDLE
      end if
   end function resolve

   subroutine ensure_host(h)
      !! If data is on GPU and host is stale, copy back to host.
      type(solver_handle_t), pointer, intent(inout) :: h

      if (h%host_is_current) return
      if (h%on_device) then
         call solver_exit_data(h%state)
         h%on_device = .false.
         h%host_is_current = .true.
      end if
   end subroutine ensure_host

   subroutine ensure_device(h)
      type(solver_handle_t), pointer, intent(inout) :: h

      if (.not. h%on_device) then
         call solver_enter_data(h%state)
         h%on_device = .true.
      end if
      h%host_is_current = .false.
   end subroutine ensure_device

   subroutine refresh_host(h)
      !! Lightweight host refresh that keeps GPU data mapped.
      type(solver_handle_t), pointer, intent(inout) :: h

      if (.not. h%is_initialised) return
      if (h%host_is_current) return
      if (h%on_device) then
         call solver_update_host(h%state)
         h%host_is_current = .true.
      end if
   end subroutine refresh_host

   subroutine push_state_to_device(h)
      !! Push solution arrays from host to GPU and mark device-canonical.
      type(solver_handle_t), pointer, intent(inout) :: h

      if (.not. h%on_device) return
      call solver_update_state(h%state)
      h%host_is_current = .false.
   end subroutine push_state_to_device

   subroutine api_init_multilayer(h)
      !! Initialise multilayer arrays from the barotropic state.
      type(solver_handle_t), pointer, intent(inout) :: h

      call init_multilayer_from_barotropic(h%state, h%cfg)
   end subroutine api_init_multilayer
   ! ================================================================
   ! Lifecycle
   ! ================================================================

   function rki_working_precision() result(nbytes) &
      bind(c, name="rki_working_precision")
      !! Working-precision width in bytes (8 = double, 4 = single) for the
      !! solver real kind `wp`.  Process-global (compile-time), handle-less:
      !! lets the Python/ctypes layer pick the numpy dtype matching the
      !! zero-copy field views, which are `c_loc` of `real(wp)` arrays.
      !! Shared by the structured and unstructured Python paths (one
      !! symbol in librakali, since `wp` is identical for both).
      integer(c_int) :: nbytes
      nbytes = int(storage_size(1.0_wp)/8, c_int)
   end function rki_working_precision

   function rki_create(c_filename, c_len, c_handle) result(ierr) &
      bind(c, name="rki_create")
      !! Phase 1: allocate handle, read config, allocate state arrays.
      !!
      !! On success, c_handle is set to the new opaque handle and ierr is 0.
      !! On failure, c_handle is null and ierr is non-zero.
      integer(c_int), intent(in), value :: c_len
      character(kind=c_char), intent(in) :: c_filename(c_len)
      type(c_ptr), intent(out) :: c_handle
      integer(c_int) :: ierr

      character(len=:), allocatable :: filename
      type(solver_handle_t), pointer :: h
      type(nml_schema_t) :: schema
      character(len=:), allocatable :: doc_dir
      integer :: status

      c_handle = handle_create()
      status = handle_check(c_handle, h)
      if (status /= HANDLE_OK) then
         ierr = RKI_ERR_BAD_HANDLE
         return
      end if

      call c_to_f_string(c_filename, c_len, filename)
      ! read_config builds + strict-parses the schema internally (unknown
      ! key / range / enum errors error-stop here); `schema` is returned
      ! for the parameter-doc dumps below.
      call read_config(filename, h%cfg, schema)
      call validate_config(h%cfg)

      ! Emit the MOM6-style parameter-doc dumps next to the input
      ! namelist (the FFI path does not bootstrap the output directory),
      ! using the schema returned by read_config (already built + parsed).
      doc_dir = nml_dirname(filename)
      call schema%write_doc_all(doc_dir//"rakali_parameter_doc.all")
      call schema%write_doc_short(doc_dir//"rakali_parameter_doc.short")

      call h%state%init_from_config(h%cfg)

      h%t_current = 0.0_wp
      h%is_initialised = .true.
      h%is_finalised = .false.
      h%user_set_bathymetry = .false.
      h%user_set_h = .false.
      h%user_set_layers = .false.
      call print_banner()
      ierr = 0
   end function rki_create

   function rki_finalize(c_handle) result(ierr) bind(c, name="rki_finalize")
      !! Phase 2: complete setup after user has set fields.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return

      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (h%is_finalised) then
         ierr = -2
         return
      end if

      associate (s => h%state, cfg => h%cfg)
         if (.not. h%user_set_bathymetry .and. len_trim(cfg%bathymetry_file) > 0) then
            call load_bathymetry(trim(cfg%bathymetry_file), s)
         end if
         call fill_bathymetry_ghosts(s)

         if (.not. h%user_set_h) then
            call init_from_config(s, cfg)
         else if (s%use_multilayer .and. .not. h%user_set_layers) then
            call api_init_multilayer(h)
         end if

         call bc_init_from_config(s%bc, cfg)
         call bc_init_chapman(s%bc, s%barotropic%h, s%barotropic%b, &
                              s%grid%nx_total, s%grid%ny_total, s%grid%nghost)

         call solver_enter_data(s)
      end associate
      h%on_device = .true.
      h%host_is_current = .false.
      h%is_finalised = .true.
      ierr = 0
   end function rki_finalize

   function rki_create_from_namelist(c_filename, c_len, c_handle) result(ierr) &
      bind(c, name="rki_create_from_namelist")
      !! One-shot create + finalize from a namelist file.
      integer(c_int), intent(in), value :: c_len
      character(kind=c_char), intent(in) :: c_filename(c_len)
      type(c_ptr), intent(out) :: c_handle
      integer(c_int) :: ierr

      ierr = rki_create(c_filename, c_len, c_handle)
      if (ierr /= 0) return
      ierr = rki_finalize(c_handle)
   end function rki_create_from_namelist

   function rki_create_from_string(c_text, c_len, c_handle) result(ierr) &
      bind(c, name="rki_create_from_string")
      !! Phase 1, fully in-memory: allocate the handle and read config from
      !! the `c_text` namelist buffer (newline-separated) — no temp file,
      !! no `chdir`, no on-disk parameter-doc dumps.  Otherwise identical
      !! to [[rki_create]] (same strict validation, state
      !! allocation).  Pair with the `set_*` setters + [[rki_finalize]].
      integer(c_int), intent(in), value :: c_len
      character(kind=c_char), intent(in) :: c_text(c_len)
      type(c_ptr), intent(out) :: c_handle
      integer(c_int) :: ierr

      character(len=:), allocatable :: text
      type(solver_handle_t), pointer :: h
      integer :: status

      c_handle = handle_create()
      status = handle_check(c_handle, h)
      if (status /= HANDLE_OK) then
         ierr = RKI_ERR_BAD_HANDLE
         return
      end if

      call c_to_f_string(c_text, c_len, text)
      ! Strict parse + apply straight into h%cfg, in memory (error-stops on
      ! any unknown key / range / enum violation, like the file path).
      call read_config_from_string(text, h%cfg)
      call validate_config(h%cfg)

      call h%state%init_from_config(h%cfg)

      h%t_current = 0.0_wp
      h%is_initialised = .true.
      h%is_finalised = .false.
      h%user_set_bathymetry = .false.
      h%user_set_h = .false.
      h%user_set_layers = .false.
      call print_banner()
      ierr = 0
   end function rki_create_from_string

   function rki_create_from_namelist_string(c_text, c_len, c_handle) &
      result(ierr) bind(c, name="rki_create_from_namelist_string")
      !! One-shot create + finalize from an in-memory namelist buffer.
      integer(c_int), intent(in), value :: c_len
      character(kind=c_char), intent(in) :: c_text(c_len)
      type(c_ptr), intent(out) :: c_handle
      integer(c_int) :: ierr

      ierr = rki_create_from_string(c_text, c_len, c_handle)
      if (ierr /= 0) return
      ierr = rki_finalize(c_handle)
   end function rki_create_from_namelist_string

   subroutine rki_destroy(c_handle) bind(c, name="rki_destroy")
      !! Destroy the simulation and free the handle. Nulls c_handle.
      !! Idempotent on null/already-destroyed handles.
      type(c_ptr), intent(inout) :: c_handle
      type(solver_handle_t), pointer :: h
      integer :: status

      status = handle_check(c_handle, h)
      if (status /= HANDLE_OK) then
         !! Still call destroy to handle the null case gracefully.
         call handle_destroy(c_handle)
         return
      end if

      !! Must do full exit_data before deallocate to unmap GPU memory and
      !! avoid "partially present" errors on a future create.
      if (h%on_device) then
         call solver_exit_data(h%state)
         h%on_device = .false.
      end if
      if (h%is_initialised) call h%state%destroy()

      call get_knowledge()
      call handle_destroy(c_handle)
   end subroutine rki_destroy

   ! ================================================================
   ! Stepping
   ! ================================================================

   function rki_step(c_handle, n_steps) result(ierr) bind(c, name="rki_step")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: n_steps
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: i

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_finalised) then
         ierr = -1
         return
      end if

      call ensure_device(h)
      h%host_is_current = .false.

      do i = 1, n_steps
         call solver_step(h%state, &
                          dt_max=(h%cfg%t_end - h%t_current), &
                          time=h%t_current)
         h%t_current = h%t_current + h%state%dt
      end do

      ierr = 0
   end function rki_step

   function rki_advance(c_handle, target_time) result(ierr) &
      bind(c, name="rki_advance")
      type(c_ptr), intent(in), value :: c_handle
      real(c_double), intent(in), value :: target_time
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_finalised) then
         ierr = -1
         return
      end if

      call ensure_device(h)
      h%host_is_current = .false.

      do while (h%t_current < target_time - 1.0e-10_wp)
         call solver_step(h%state, &
                          dt_max=real(target_time - h%t_current, wp), &
                          time=h%t_current)
         h%t_current = h%t_current + h%state%dt
      end do

      ierr = 0
   end function rki_advance

   ! ================================================================
   ! Queries
   ! ================================================================

   function rki_get_time(c_handle) result(t) bind(c, name="rki_get_time")
      type(c_ptr), intent(in), value :: c_handle
      real(c_double) :: t
      type(solver_handle_t), pointer :: h
      integer :: status

      status = handle_check(c_handle, h)
      if (status == HANDLE_OK) then
         t = h%t_current
      else
         t = 0.0_wp
      end if
   end function rki_get_time

   function rki_get_dt(c_handle) result(dt) bind(c, name="rki_get_dt")
      type(c_ptr), intent(in), value :: c_handle
      real(c_double) :: dt
      type(solver_handle_t), pointer :: h
      integer :: status

      status = handle_check(c_handle, h)
      if (status == HANDLE_OK .and. h%is_initialised) then
         dt = h%state%dt
      else
         dt = 0.0_wp
      end if
   end function rki_get_dt

   subroutine rki_get_grid_info(c_handle, nx, ny, nghost, dx, dy) &
      bind(c, name="rki_get_grid_info")
      !! Return physical grid dimensions (excluding ghost cells).
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(out) :: nx, ny, nghost
      real(c_double), intent(out) :: dx, dy
      type(solver_handle_t), pointer :: h
      integer :: status

      status = handle_check(c_handle, h)
      if (status == HANDLE_OK .and. h%is_initialised) then
         nx = h%state%grid%nx_phys
         ny = h%state%grid%ny_phys
         nghost = h%state%grid%nghost
         dx = h%state%grid%dx
         dy = h%state%grid%dy
      else
         nx = 0
         ny = 0
         nghost = 0
         dx = 0.0_wp
         dy = 0.0_wp
      end if
   end subroutine rki_get_grid_info

   function rki_get_step_count(c_handle) result(n) &
      bind(c, name="rki_get_step_count")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int) :: n
      type(solver_handle_t), pointer :: h
      integer :: status

      status = handle_check(c_handle, h)
      if (status == HANDLE_OK .and. h%is_initialised) then
         n = h%state%step_count
      else
         n = 0
      end if
   end function rki_get_step_count

   function rki_get_nz_ml(c_handle) result(nz) bind(c, name="rki_get_nz_ml")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int) :: nz
      type(solver_handle_t), pointer :: h
      integer :: status

      status = handle_check(c_handle, h)
      if (status == HANDLE_OK .and. h%is_initialised .and. h%state%use_multilayer) then
         nz = h%state%multilayer%nz_ml
      else
         nz = 0
      end if
   end function rki_get_nz_ml

   ! ================================================================
   ! Zero-copy field accessors (return c_ptr to Fortran data)
   ! ================================================================

   function rki_get_h_ptr(c_handle, ptr, nx_total, ny_total) result(ierr) &
      bind(c, name="rki_get_h_ptr")
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx_total, ny_total
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx_total = 0
      ny_total = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      call ensure_host(h)
      tmp => h%state%barotropic%h
      ptr = c_loc(tmp(1, 1))
      nx_total = h%state%grid%nx_total
      ny_total = h%state%grid%ny_total
      ierr = 0
   end function rki_get_h_ptr

   function rki_get_hu_ptr(c_handle, ptr, nx_total, ny_total) result(ierr) &
      bind(c, name="rki_get_hu_ptr")
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx_total, ny_total
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx_total = 0
      ny_total = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      call ensure_host(h)
      tmp => h%state%barotropic%hu
      ptr = c_loc(tmp(1, 1))
      nx_total = h%state%grid%nx_total
      ny_total = h%state%grid%ny_total
      ierr = 0
   end function rki_get_hu_ptr

   function rki_get_hv_ptr(c_handle, ptr, nx_total, ny_total) result(ierr) &
      bind(c, name="rki_get_hv_ptr")
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx_total, ny_total
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx_total = 0
      ny_total = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      call ensure_host(h)
      tmp => h%state%barotropic%hv
      ptr = c_loc(tmp(1, 1))
      nx_total = h%state%grid%nx_total
      ny_total = h%state%grid%ny_total
      ierr = 0
   end function rki_get_hv_ptr

   function rki_get_b_ptr(c_handle, ptr, nx_total, ny_total) result(ierr) &
      bind(c, name="rki_get_b_ptr")
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx_total, ny_total
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx_total = 0
      ny_total = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      call ensure_host(h)
      tmp => h%state%barotropic%b
      ptr = c_loc(tmp(1, 1))
      nx_total = h%state%grid%nx_total
      ny_total = h%state%grid%ny_total
      ierr = 0
   end function rki_get_b_ptr

   ! RK2 save-buffer accessors (h^n / hu^n / hv^n captured at start of stage 1
   ! by fused_update_stage1, consumed by fused_update_stage2 for the average).
   ! Needed for mid-RK checkpointing.

   function rki_get_h0_ptr(c_handle, ptr, nx_total, ny_total) result(ierr) &
      bind(c, name="rki_get_h0_ptr")
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx_total, ny_total
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx_total = 0
      ny_total = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      call ensure_host(h)
      tmp => h%state%barotropic%h0
      ptr = c_loc(tmp(1, 1))
      nx_total = h%state%grid%nx_total
      ny_total = h%state%grid%ny_total
      ierr = 0
   end function rki_get_h0_ptr

   function rki_get_hu0_ptr(c_handle, ptr, nx_total, ny_total) result(ierr) &
      bind(c, name="rki_get_hu0_ptr")
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx_total, ny_total
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx_total = 0
      ny_total = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      call ensure_host(h)
      tmp => h%state%barotropic%hu0
      ptr = c_loc(tmp(1, 1))
      nx_total = h%state%grid%nx_total
      ny_total = h%state%grid%ny_total
      ierr = 0
   end function rki_get_hu0_ptr

   function rki_get_hv0_ptr(c_handle, ptr, nx_total, ny_total) result(ierr) &
      bind(c, name="rki_get_hv0_ptr")
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx_total, ny_total
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx_total = 0
      ny_total = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      call ensure_host(h)
      tmp => h%state%barotropic%hv0
      ptr = c_loc(tmp(1, 1))
      nx_total = h%state%grid%nx_total
      ny_total = h%state%grid%ny_total
      ierr = 0
   end function rki_get_hv0_ptr

   ! 3D multilayer field accessors

   function rki_get_h_layer_ptr(c_handle, ptr, nx_total, ny_total, nz) result(ierr) &
      bind(c, name="rki_get_h_layer_ptr")
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx_total, ny_total, nz
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx_total = 0
      ny_total = 0
      nz = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. .not. h%state%use_multilayer) then
         ierr = -1
         return
      end if
      call ensure_host(h)
      tmp => h%state%multilayer%h_layer
      ptr = c_loc(tmp(1, 1, 1))
      nx_total = h%state%grid%nx_total
      ny_total = h%state%grid%ny_total
      nz = h%state%multilayer%nz_ml
      ierr = 0
   end function rki_get_h_layer_ptr

   function rki_get_hu_layer_ptr(c_handle, ptr, nx_total, ny_total, nz) result(ierr) &
      bind(c, name="rki_get_hu_layer_ptr")
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx_total, ny_total, nz
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx_total = 0
      ny_total = 0
      nz = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. .not. h%state%use_multilayer) then
         ierr = -1
         return
      end if
      call ensure_host(h)
      tmp => h%state%multilayer%hu_layer
      ptr = c_loc(tmp(1, 1, 1))
      nx_total = h%state%grid%nx_total
      ny_total = h%state%grid%ny_total
      nz = h%state%multilayer%nz_ml
      ierr = 0
   end function rki_get_hu_layer_ptr

   function rki_get_hv_layer_ptr(c_handle, ptr, nx_total, ny_total, nz) result(ierr) &
      bind(c, name="rki_get_hv_layer_ptr")
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx_total, ny_total, nz
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx_total = 0
      ny_total = 0
      nz = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. .not. h%state%use_multilayer) then
         ierr = -1
         return
      end if
      call ensure_host(h)
      tmp => h%state%multilayer%hv_layer
      ptr = c_loc(tmp(1, 1, 1))
      nx_total = h%state%grid%nx_total
      ny_total = h%state%grid%ny_total
      nz = h%state%multilayer%nz_ml
      ierr = 0
   end function rki_get_hv_layer_ptr

   function rki_get_salinity_ptr(c_handle, ptr, nx_total, ny_total, nz) result(ierr) &
      bind(c, name="rki_get_salinity_ptr")
      !! Returns pointer to hS_layer (conserved salinity = h*S).
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx_total, ny_total, nz
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx_total = 0
      ny_total = 0
      nz = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. .not. h%state%use_multilayer) then
         ierr = -1
         return
      end if
      call ensure_host(h)
      tmp => h%state%multilayer%tracers(h%state%multilayer%idx_salinity)%hTr
      ptr = c_loc(tmp(1, 1, 1))
      nx_total = h%state%grid%nx_total
      ny_total = h%state%grid%ny_total
      nz = h%state%multilayer%nz_ml
      ierr = 0
   end function rki_get_salinity_ptr

   function rki_get_temperature_ptr(c_handle, ptr, nx_total, ny_total, nz) result(ierr) &
      bind(c, name="rki_get_temperature_ptr")
      !! Returns pointer to hT_layer (conserved temperature = h*T).
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx_total, ny_total, nz
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx_total = 0
      ny_total = 0
      nz = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. .not. h%state%use_multilayer) then
         ierr = -1
         return
      end if
      call ensure_host(h)
      tmp => h%state%multilayer%tracers(h%state%multilayer%idx_temperature)%hTr
      ptr = c_loc(tmp(1, 1, 1))
      nx_total = h%state%grid%nx_total
      ny_total = h%state%grid%ny_total
      nz = h%state%multilayer%nz_ml
      ierr = 0
   end function rki_get_temperature_ptr

   ! ================================================================
   ! Setters (pre-finalize fields)
   ! ================================================================

   function rki_set_bathymetry(c_handle, b_data, nx_p, ny_p) result(ierr) &
      bind(c, name="rki_set_bathymetry")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: b_data(nx_p, ny_p)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: ng

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (nx_p /= h%state%grid%nx_phys .or. ny_p /= h%state%grid%ny_phys) then
         ierr = -2
         return
      end if

      if (h%on_device) call ensure_host(h)
      ng = h%state%grid%nghost
      h%state%barotropic%b(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = b_data(:, :)
      h%user_set_bathymetry = .true.

      if (h%is_finalised) call fill_bathymetry_ghosts(h%state)
      ierr = 0
   end function rki_set_bathymetry

   function rki_set_initial_h(c_handle, h_data, nx_p, ny_p) result(ierr) &
      bind(c, name="rki_set_initial_h")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: h_data(nx_p, ny_p)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: ng

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (nx_p /= h%state%grid%nx_phys .or. ny_p /= h%state%grid%ny_phys) then
         ierr = -2
         return
      end if

      if (h%on_device) call ensure_host(h)
      ng = h%state%grid%nghost
      h%state%barotropic%h(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = h_data(:, :)
      h%user_set_h = .true.
      ierr = 0
   end function rki_set_initial_h

   ! ================================================================
   ! Forcing setters
   ! ================================================================

   function rki_set_wind(c_handle, wx_data, wy_data, nx_p, ny_p) result(ierr) &
      bind(c, name="rki_set_wind")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: wx_data(nx_p, ny_p)
      real(c_double), intent(in) :: wy_data(nx_p, ny_p)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: ng

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (nx_p /= h%state%grid%nx_phys .or. ny_p /= h%state%grid%ny_phys) then
         ierr = -2
         return
      end if

      ng = h%state%grid%nghost
      h%state%forcing%wind_x(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = wx_data(:, :)
      h%state%forcing%wind_y(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = wy_data(:, :)
      h%state%forcing%active = .true.

      if (h%on_device) call solver_update_forcing(h%state)
      ierr = 0
   end function rki_set_wind

   function rki_set_atm_pressure(c_handle, p_data, nx_p, ny_p) result(ierr) &
      bind(c, name="rki_set_atm_pressure")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: p_data(nx_p, ny_p)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: ng

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (nx_p /= h%state%grid%nx_phys .or. ny_p /= h%state%grid%ny_phys) then
         ierr = -2
         return
      end if

      ng = h%state%grid%nghost
      h%state%forcing%atm_pressure(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = p_data(:, :)
      h%state%forcing%active = .true.

      if (h%on_device) call solver_update_forcing(h%state)
      ierr = 0
   end function rki_set_atm_pressure

   function rki_set_precipitation(c_handle, p_data, nx_p, ny_p) result(ierr) &
      bind(c, name="rki_set_precipitation")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: p_data(nx_p, ny_p)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: ng

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (nx_p /= h%state%grid%nx_phys .or. ny_p /= h%state%grid%ny_phys) then
         ierr = -2
         return
      end if

      ng = h%state%grid%nghost
      h%state%forcing%precipitation(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = p_data(:, :)
      h%state%forcing%active = .true.

      if (h%on_device) call solver_update_forcing(h%state)
      ierr = 0
   end function rki_set_precipitation

   function rki_set_heat_flux(c_handle, q_data, nx_p, ny_p) result(ierr) &
      bind(c, name="rki_set_heat_flux")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: q_data(nx_p, ny_p)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: ng

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (nx_p /= h%state%grid%nx_phys .or. ny_p /= h%state%grid%ny_phys) then
         ierr = -2
         return
      end if

      ng = h%state%grid%nghost
      h%state%forcing%heat_flux(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = q_data(:, :)

      if (h%on_device) call solver_update_forcing(h%state)
      ierr = 0
   end function rki_set_heat_flux

   function rki_set_evap_minus_precip(c_handle, ep_data, nx_p, ny_p) result(ierr) &
      bind(c, name="rki_set_evap_minus_precip")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: ep_data(nx_p, ny_p)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: ng

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (nx_p /= h%state%grid%nx_phys .or. ny_p /= h%state%grid%ny_phys) then
         ierr = -2
         return
      end if

      ng = h%state%grid%nghost
      h%state%forcing%evap_minus_precip(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = ep_data(:, :)

      if (h%on_device) call solver_update_forcing(h%state)
      ierr = 0
   end function rki_set_evap_minus_precip

   ! ================================================================
   ! Data assimilation setters (mid-simulation state injection)
   ! ================================================================

   function rki_set_h(c_handle, h_data, nx_p, ny_p) result(ierr) &
      bind(c, name="rki_set_h")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: h_data(nx_p, ny_p)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: ng

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (nx_p /= h%state%grid%nx_phys .or. ny_p /= h%state%grid%ny_phys) then
         ierr = -2
         return
      end if

      call refresh_host(h)
      ng = h%state%grid%nghost
      h%state%barotropic%h(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = h_data(:, :)
      h%user_set_h = .true.

      call push_state_to_device(h)
      ierr = 0
   end function rki_set_h

   function rki_set_hu(c_handle, hu_data, nx_p, ny_p) result(ierr) &
      bind(c, name="rki_set_hu")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: hu_data(nx_p, ny_p)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: ng

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (nx_p /= h%state%grid%nx_phys .or. ny_p /= h%state%grid%ny_phys) then
         ierr = -2
         return
      end if

      call refresh_host(h)
      ng = h%state%grid%nghost
      h%state%barotropic%hu(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = hu_data(:, :)

      call push_state_to_device(h)
      ierr = 0
   end function rki_set_hu

   function rki_set_hv(c_handle, hv_data, nx_p, ny_p) result(ierr) &
      bind(c, name="rki_set_hv")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: hv_data(nx_p, ny_p)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: ng

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (nx_p /= h%state%grid%nx_phys .or. ny_p /= h%state%grid%ny_phys) then
         ierr = -2
         return
      end if

      call refresh_host(h)
      ng = h%state%grid%nghost
      h%state%barotropic%hv(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = hv_data(:, :)

      call push_state_to_device(h)
      ierr = 0
   end function rki_set_hv

   function rki_set_salinity(c_handle, s_data, nx_p, ny_p, nz_p) result(ierr) &
      bind(c, name="rki_set_salinity")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p, nz_p
      real(c_double), intent(in) :: s_data(nx_p, ny_p, nz_p)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: ng, i, j, k

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (.not. h%state%use_multilayer) then
         ierr = -3
         return
      end if
      if (nx_p /= h%state%grid%nx_phys .or. ny_p /= h%state%grid%ny_phys &
          .or. nz_p /= h%state%multilayer%nz_ml) then
         ierr = -2
         return
      end if

      call refresh_host(h)
      ng = h%state%grid%nghost
      do k = 1, nz_p
         do j = 1, ny_p
            do i = 1, nx_p
               h%state%multilayer%tracers(h%state%multilayer%idx_salinity)%hTr(ng + i, ng + j, k) = &
                  h%state%multilayer%h_layer(ng + i, ng + j, k)*s_data(i, j, k)
            end do
         end do
      end do

      call push_state_to_device(h)
      ierr = 0
   end function rki_set_salinity

   function rki_set_temperature(c_handle, t_data, nx_p, ny_p, nz_p) result(ierr) &
      bind(c, name="rki_set_temperature")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p, nz_p
      real(c_double), intent(in) :: t_data(nx_p, ny_p, nz_p)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: ng, i, j, k

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (.not. h%state%use_multilayer) then
         ierr = -3
         return
      end if
      if (nx_p /= h%state%grid%nx_phys .or. ny_p /= h%state%grid%ny_phys &
          .or. nz_p /= h%state%multilayer%nz_ml) then
         ierr = -2
         return
      end if

      call refresh_host(h)
      ng = h%state%grid%nghost
      do k = 1, nz_p
         do j = 1, ny_p
            do i = 1, nx_p
               h%state%multilayer%tracers(h%state%multilayer%idx_temperature)%hTr(ng + i, ng + j, k) = &
                  h%state%multilayer%h_layer(ng + i, ng + j, k)*t_data(i, j, k)
            end do
         end do
      end do

      call push_state_to_device(h)
      ierr = 0
   end function rki_set_temperature

   ! ================================================================
   ! Boundary condition setters (mid-simulation forcing updates)
   ! ================================================================

   function rki_set_inflow_velocity(c_handle, v) result(ierr) &
      bind(c, name="rki_set_inflow_velocity")
      type(c_ptr), intent(in), value :: c_handle
      real(c_double), intent(in), value :: v
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%bc%inflow_velocity = real(v, wp)
      ierr = 0
   end function rki_set_inflow_velocity

   function rki_set_inflow_salinity(c_handle, s) result(ierr) &
      bind(c, name="rki_set_inflow_salinity")
      type(c_ptr), intent(in), value :: c_handle
      real(c_double), intent(in), value :: s
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%bc%inflow_salinity = real(s, wp)
      ierr = 0
   end function rki_set_inflow_salinity

   function rki_set_inflow_temperature(c_handle, t) result(ierr) &
      bind(c, name="rki_set_inflow_temperature")
      type(c_ptr), intent(in), value :: c_handle
      real(c_double), intent(in), value :: t
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%bc%inflow_temperature = real(t, wp)
      ierr = 0
   end function rki_set_inflow_temperature

   function rki_set_discharge_flux(c_handle, q) result(ierr) &
      bind(c, name="rki_set_discharge_flux")
      type(c_ptr), intent(in), value :: c_handle
      real(c_double), intent(in), value :: q
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%bc%discharge_flux = real(q, wp)
      ierr = 0
   end function rki_set_discharge_flux

   function rki_set_clamped_elevation(c_handle, eta) result(ierr) &
      bind(c, name="rki_set_clamped_elevation")
      type(c_ptr), intent(in), value :: c_handle
      real(c_double), intent(in), value :: eta
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%bc%clamped_elevation = real(eta, wp)
      ierr = 0
   end function rki_set_clamped_elevation

   function rki_get_clamped_elevation(c_handle) result(eta) &
      bind(c, name="rki_get_clamped_elevation")
      type(c_ptr), intent(in), value :: c_handle
      real(c_double) :: eta
      type(solver_handle_t), pointer :: h
      integer :: status

      status = handle_check(c_handle, h)
      if (status == HANDLE_OK .and. h%is_initialised) then
         eta = h%state%bc%clamped_elevation
      else
         eta = 0.0_wp
      end if
   end function rki_get_clamped_elevation

   function rki_get_manning_n(c_handle) result(n) &
      bind(c, name="rki_get_manning_n")
      type(c_ptr), intent(in), value :: c_handle
      real(c_double) :: n
      type(solver_handle_t), pointer :: h
      integer :: status

      status = handle_check(c_handle, h)
      if (status == HANDLE_OK .and. h%is_initialised) then
         n = h%state%barotropic%manning_n
      else
         n = 0.0_wp
      end if
   end function rki_get_manning_n

   function rki_get_coriolis_f(c_handle) result(f) &
      bind(c, name="rki_get_coriolis_f")
      type(c_ptr), intent(in), value :: c_handle
      real(c_double) :: f
      type(solver_handle_t), pointer :: h
      integer :: status

      status = handle_check(c_handle, h)
      if (status == HANDLE_OK .and. h%is_initialised) then
         f = h%state%barotropic%coriolis_f
      else
         f = 0.0_wp
      end if
   end function rki_get_coriolis_f

   function rki_set_tidal_constituents(c_handle, amp, phase, omega, n) result(ierr) &
      bind(c, name="rki_set_tidal_constituents")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: n
      real(c_double), intent(in) :: amp(n), phase(n), omega(n)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      integer :: k

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (n < 0 .or. n > MAX_TIDAL_CONSTITUENTS) then
         ierr = -2
         return
      end if

      h%state%bc%tidal_amp = 0.0_wp
      h%state%bc%tidal_phase = 0.0_wp
      h%state%bc%tidal_omega = 0.0_wp

      do k = 1, n
         h%state%bc%tidal_amp(k) = real(amp(k), wp)
         h%state%bc%tidal_phase(k) = real(phase(k), wp)
         h%state%bc%tidal_omega(k) = real(omega(k), wp)
      end do
      h%state%bc%n_tidal_constituents = n
      ierr = 0
   end function rki_set_tidal_constituents

   ! ================================================================
   ! GPU sync
   ! ================================================================

   subroutine rki_sync_host(c_handle) bind(c, name="rki_sync_host")
      !! Lightweight device-to-host sync for callbacks.
      type(c_ptr), intent(in), value :: c_handle
      type(solver_handle_t), pointer :: h
      integer :: status

      status = handle_check(c_handle, h)
      if (status /= HANDLE_OK) return
      if (h%is_initialised .and. h%on_device .and. .not. h%host_is_current) then
         call solver_update_host(h%state)
         h%host_is_current = .true.
      end if
   end subroutine rki_sync_host

   subroutine rki_sync_from_device(c_handle) bind(c, name="rki_sync_from_device")
      !! Full device-to-host sync that removes data from GPU.
      type(c_ptr), intent(in), value :: c_handle
      type(solver_handle_t), pointer :: h
      integer :: status

      status = handle_check(c_handle, h)
      if (status /= HANDLE_OK) return
      if (h%is_initialised .and. h%on_device) then
         call solver_exit_data(h%state)
         h%on_device = .false.
         h%host_is_current = .true.
      end if
   end subroutine rki_sync_from_device

   ! ================================================================
   ! Restart save / load
   ! ================================================================

   function rki_save_restart(c_handle, c_filename, c_len) result(ierr) &
      bind(c, name="rki_save_restart")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: c_len
      character(kind=c_char), intent(in) :: c_filename(c_len)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      character(len=:), allocatable :: filename

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (.not. h%is_finalised) then
         ierr = -2
         return
      end if

      call c_to_f_string(c_filename, c_len, filename)
      call refresh_host(h)
      call restart_write(filename, h%state, h%t_current)
      ierr = 0
   end function rki_save_restart

   function rki_load_restart(c_handle, c_filename, c_len) result(ierr) &
      bind(c, name="rki_load_restart")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: c_len
      character(kind=c_char), intent(in) :: c_filename(c_len)
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h
      character(len=:), allocatable :: filename
      real(wp) :: t_restart
      integer :: read_ierr
      logical :: loaded_ml

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (h%is_finalised) then
         ierr = -2
         return
      end if

      call c_to_f_string(c_filename, c_len, filename)
      call restart_read(filename, h%state, t_restart, &
                        ierr=read_ierr, loaded_multilayer=loaded_ml)
      if (read_ierr /= 0) then
         ierr = -3
         return
      end if

      h%t_current = t_restart
      h%user_set_h = .true.
      h%user_set_bathymetry = .true.
      h%user_set_layers = loaded_ml
      ierr = 0
   end function rki_load_restart

   function rki_set_use_ppm_tracer(c_handle, flag) result(ierr) &
      bind(c, name="rki_set_use_ppm_tracer")
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: flag
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%multilayer%use_ppm_tracer = (flag /= 0)
      ierr = 0
   end function rki_set_use_ppm_tracer

   function rki_set_use_hllc(c_handle, flag) result(ierr) &
      bind(c, name="rki_set_use_hllc")
      !! Toggle the structured barotropic Riemann solver between HLL
      !! (`flag == 0`) and HLLC (`flag /= 0`).  Safe to call mid-run.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: flag
      integer(c_int) :: ierr
      type(solver_handle_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%barotropic%use_hllc = (flag /= 0)
      ierr = 0
   end function rki_set_use_hllc

end module rki_api
