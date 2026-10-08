!! C-compatible API for the unstructured Rakali solver
!!
!! Provides bind(C) entrypoints for Python (and any other ctypes/cffi
!! caller) to drive the unstructured central-upwind solver without going
!! through any file I/O.  The mesh and all initial state come from numpy
!! arrays handed across the ctypes boundary.
!!
!! Multi-instance: each call to `rki_unstr_create` returns an opaque
!! `void*` handle.  Multiple unstructured handles can coexist in one
!! process, and a structured + unstructured handle can coexist too
!! (they are separate types tracked by distinct magic numbers).
!!
!! `rki_unstr_set_log_level` intentionally stays module-level: the logger
!! is process-global, not per simulation.
!!
!! Lifecycle:
!!
!!   void* h;
!!   rki_unstr_create(nverts, ntri, vert_x, vert_y, tri_v, bc_tags, &h)
!!   rki_unstr_set_bathymetry(h, bed)             -- centroid array
!!   rki_unstr_set_initial_stage(h, stage)        -- centroid array
!!   rki_unstr_set_initial_xmom(h, xmom)          -- optional
!!   rki_unstr_set_initial_ymom(h, ymom)          -- optional
!!   rki_unstr_set_bc_type(h, group, bc_type)
!!   rki_unstr_set_manning(h, n)                  -- optional
!!   rki_unstr_set_wind(h, tau_x, tau_y)          -- optional, uniform Pa
!!   rki_unstr_set_wind_field(h, wind_x, wind_y)  -- optional, per-cell Pa
!!   rki_unstr_finalize(h)
!!   rki_unstr_step(h, dt) / rki_unstr_advance(h, t_end)
!!   rki_unstr_get_*_ptr(h, ptr, n)
!!   rki_unstr_destroy(&h)
!!
!! Error codes (mirror rki_api):
!!   0     success
!!   -1    not initialised / output not open / etc.
!!   -2    range / already-open / grid-mismatch
!!   -3    multilayer not active / restart read failure
!!  -100   invalid handle
module rki_api_unstr
   use, intrinsic :: iso_c_binding, only: c_associated, c_double, &
                                                                             c_f_pointer, c_int, &
                                                                             c_loc, c_null_ptr, c_ptr
   use rki_constants, only: wp, VCOORD_ZSTAR_FULL
   use rki_mesh_unstr, only: mesh_unstr_t, mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, state_unstr_init, &
                              state_unstr_destroy, state_unstr_enter_data, &
                              state_unstr_exit_data
   use rki_solver_unstr, only: solver_unstr_step, solver_unstr_advance, &
                               solver_unstr_update_host, &
                               solver_unstr_update_forcing, &
                               solver_unstr_init_multilayer
   use rki_output_unstr, only: output_unstr_init, &
                               output_unstr_write_snapshot, output_unstr_close
   use rki_restart_unstr, only: restart_write_unstr, restart_read_unstr
   use rki_gauges_unstr, only: gauges_unstr_init, &
                               gauges_unstr_write, gauges_unstr_close
   use rki_kernel_timestep_unstr, only: compute_timestep_unstr
   use rki_timestepper, only: parse_time_scheme, TIME_RK2
   use rki_boundary_unstr, only: MAX_BC_GROUPS
   use rki_handle, only: solver_handle_unstr_t, handle_create_unstr, &
                         handle_destroy_unstr, handle_check_unstr, HANDLE_OK
   use pic_logger, only: logger => global_logger, &
                         debug_level, verbose_level, info_level, &
                         performance_level, warning_level, error_level
   use rki_banner, only: print_banner
   use pic_knowledge, only: get_knowledge
   implicit none
   private

   integer(c_int), parameter :: RKI_ERR_BAD_HANDLE = -100

contains

   ! ================================================================
   ! Helpers
   ! ================================================================

   function resolve(c_handle, h) result(ierr)
      !! Resolve an opaque c_ptr to a Fortran handle pointer.
      type(c_ptr), intent(in) :: c_handle
      type(solver_handle_unstr_t), pointer, intent(out) :: h
      integer(c_int) :: ierr
      integer :: status

      status = handle_check_unstr(c_handle, h)
      if (status == HANDLE_OK) then
         ierr = 0
      else
         ierr = RKI_ERR_BAD_HANDLE
      end if
   end function resolve

   subroutine ensure_host(h)
      !! Bring state arrays back to the host if they live on the device.
      !! Once called, on_device flips to false and a subsequent step
      !! will need to re-finalize.
      type(solver_handle_unstr_t), pointer, intent(inout) :: h

      if (h%on_device) then
         call state_unstr_exit_data(h%state)
         call mesh_unstr_exit_data(h%mesh)
         h%on_device = .false.
         h%is_finalised = .false.
      end if
   end subroutine ensure_host

   subroutine refresh_host(h)
      !! Lightweight host refresh that keeps the GPU mapping intact.
      !!
      !! Use this from read-only Python getters. Tearing down via
      !! `ensure_host` would force `prepare()` again on the next step.
      type(solver_handle_unstr_t), pointer, intent(inout) :: h

      if (.not. h%is_initialised) return
      if (h%on_device) then
         call solver_unstr_update_host(h%state)
      end if
   end subroutine refresh_host

   ! ================================================================
   ! Process-global logger configuration (no handle)
   ! ================================================================

   function rki_unstr_set_log_level(level) result(ierr) &
      bind(c, name="rki_unstr_set_log_level")
      !! Process-global logger configuration.  No handle needed.
      integer(c_int), value :: level
      integer(c_int) :: ierr
      call logger%configure(int(level))
      ierr = 0
   end function rki_unstr_set_log_level

   ! ================================================================
   ! Lifecycle
   ! ================================================================

   function rki_unstr_create(nverts, ntri, vert_x, vert_y, tri_v, bc_tags, &
                             c_handle) result(ierr) &
      bind(c, name="rki_unstr_create")
      !! Build mesh + state from raw arrays into a new handle.
      !!
      !! `vert_x`, `vert_y` are length `nverts`. `tri_v` is shape
      !! `(3, ntri)` (column-major), 1-based vertex ids in CCW order.
      !! `bc_tags` is shape `(3, ntri)`: integer boundary group id per
      !! edge, 0 for interior or default-wall edges. Pass a NULL pointer
      !! (e.g. None from Python) to default everything to wall.
      integer(c_int), value :: nverts
      integer(c_int), value :: ntri
      real(c_double), intent(in) :: vert_x(nverts)
      real(c_double), intent(in) :: vert_y(nverts)
      integer(c_int), intent(in) :: tri_v(3, ntri)
      type(c_ptr), value :: bc_tags
      type(c_ptr), intent(out) :: c_handle
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      integer :: status

      c_handle = handle_create_unstr()
      status = handle_check_unstr(c_handle, h)
      if (status /= HANDLE_OK) then
         ierr = RKI_ERR_BAD_HANDLE
         return
      end if

      if (c_associated(bc_tags)) then
         block
            integer, pointer :: tags(:, :)
            call c_f_pointer(bc_tags, tags, [3, ntri])
            call mesh_unstr_init_from_triangles(h%mesh, real(vert_x, wp), &
                                                real(vert_y, wp), tri_v, &
                                                bc_tags=tags)
         end block
      else
         call mesh_unstr_init_from_triangles(h%mesh, real(vert_x, wp), &
                                             real(vert_y, wp), tri_v)
      end if

      call state_unstr_init(h%state, h%mesh)

      h%is_initialised = .true.
      h%is_finalised = .false.
      h%on_device = .false.
      ierr = 0
   end function rki_unstr_create

   function rki_unstr_enable_multilayer(c_handle, nz_ml) result(ierr) &
      bind(c, name="rki_unstr_enable_multilayer")
      type(c_ptr), value :: c_handle
      integer(c_int), value :: nz_ml
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (h%is_finalised) then
         ierr = -1
         return
      end if
      if (nz_ml < 1) then
         ierr = -2
         return
      end if

      call state_unstr_init(h%state, h%mesh, nz_ml=nz_ml)
      ierr = 0
   end function rki_unstr_enable_multilayer

   function rki_unstr_get_nz_ml(c_handle) result(nz) &
      bind(c, name="rki_unstr_get_nz_ml")
      type(c_ptr), value :: c_handle
      integer(c_int) :: nz
      type(solver_handle_unstr_t), pointer :: h
      integer :: status

      status = handle_check_unstr(c_handle, h)
      if (status /= HANDLE_OK .or. .not. h%is_initialised) then
         nz = 0
         return
      end if
      if (.not. h%state%use_multilayer) then
         nz = 0
         return
      end if
      nz = h%state%multilayer%nz_ml
   end function rki_unstr_get_nz_ml

   function rki_unstr_finalize(c_handle) result(ierr) &
      bind(c, name="rki_unstr_finalize")
      type(c_ptr), value :: c_handle
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (h%is_finalised) then
         ierr = 0
         return
      end if

      call mesh_unstr_enter_data(h%mesh)
      call state_unstr_enter_data(h%state)

      if (h%state%use_multilayer) then
         call solver_unstr_init_multilayer(h%state, h%mesh)
      end if

      h%on_device = .true.
      h%is_finalised = .true.
      ierr = 0
   end function rki_unstr_finalize

   subroutine rki_unstr_destroy(c_handle) bind(c, name="rki_unstr_destroy")
      !! Release device data, deallocate mesh + state, free the handle
      !! and null the c_ptr.  Idempotent on null/destroyed handles.
      type(c_ptr), intent(inout) :: c_handle
      type(solver_handle_unstr_t), pointer :: h
      integer :: status

      status = handle_check_unstr(c_handle, h)
      if (status /= HANDLE_OK) then
         call handle_destroy_unstr(c_handle)
         return
      end if

      if (h%gauges_open) then
         call gauges_unstr_close(h%gauges)
         h%gauges_open = .false.
      end if
      if (h%output_open) then
         call output_unstr_close(h%ostate)
         h%output_open = .false.
      end if
      call ensure_host(h)
      if (h%is_initialised) then
         call state_unstr_destroy(h%state)
         call mesh_unstr_destroy(h%mesh)
      end if

      call get_knowledge()
      call handle_destroy_unstr(c_handle)
   end subroutine rki_unstr_destroy

   ! ================================================================
   ! Stepping
   ! ================================================================

   function rki_unstr_step(c_handle, dt) result(ierr) &
      bind(c, name="rki_unstr_step")
      type(c_ptr), value :: c_handle
      real(c_double), value :: dt
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_finalised) then
         ierr = -1
         return
      end if
      call solver_unstr_step(h%state, h%mesh, real(dt, wp))
      ierr = 0
   end function rki_unstr_step

   function rki_unstr_advance(c_handle, t_end) result(ierr) &
      bind(c, name="rki_unstr_advance")
      type(c_ptr), value :: c_handle
      real(c_double), value :: t_end
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_finalised) then
         ierr = -1
         return
      end if
      call solver_unstr_advance(h%state, h%mesh, real(t_end, wp))
      ierr = 0
   end function rki_unstr_advance

   ! ================================================================
   ! Queries
   ! ================================================================

   function rki_unstr_get_time(c_handle) result(t) &
      bind(c, name="rki_unstr_get_time")
      type(c_ptr), value :: c_handle
      real(c_double) :: t
      type(solver_handle_unstr_t), pointer :: h
      integer :: status

      status = handle_check_unstr(c_handle, h)
      if (status == HANDLE_OK .and. h%is_initialised) then
         t = real(h%state%time, c_double)
      else
         t = 0.0_c_double
      end if
   end function rki_unstr_get_time

   function rki_unstr_get_dt(c_handle) result(dt) &
      bind(c, name="rki_unstr_get_dt")
      type(c_ptr), value :: c_handle
      real(c_double) :: dt
      type(solver_handle_unstr_t), pointer :: h
      integer :: status

      status = handle_check_unstr(c_handle, h)
      if (status == HANDLE_OK .and. h%is_initialised) then
         dt = real(h%state%dt, c_double)
      else
         dt = 0.0_c_double
      end if
   end function rki_unstr_get_dt

   function rki_unstr_get_step_count(c_handle) result(n) &
      bind(c, name="rki_unstr_get_step_count")
      type(c_ptr), value :: c_handle
      integer(c_int) :: n
      type(solver_handle_unstr_t), pointer :: h
      integer :: status

      status = handle_check_unstr(c_handle, h)
      if (status == HANDLE_OK .and. h%is_initialised) then
         n = h%state%step_count
      else
         n = 0
      end if
   end function rki_unstr_get_step_count

   function rki_unstr_get_ncells_owned(c_handle) result(n) &
      bind(c, name="rki_unstr_get_ncells_owned")
      type(c_ptr), value :: c_handle
      integer(c_int) :: n
      type(solver_handle_unstr_t), pointer :: h
      integer :: status

      status = handle_check_unstr(c_handle, h)
      if (status == HANDLE_OK .and. h%is_initialised) then
         n = h%state%ncells_owned
      else
         n = 0
      end if
   end function rki_unstr_get_ncells_owned

   function rki_unstr_get_ncells_total(c_handle) result(n) &
      bind(c, name="rki_unstr_get_ncells_total")
      type(c_ptr), value :: c_handle
      integer(c_int) :: n
      type(solver_handle_unstr_t), pointer :: h
      integer :: status

      status = handle_check_unstr(c_handle, h)
      if (status == HANDLE_OK .and. h%is_initialised) then
         n = h%state%ncells_total
      else
         n = 0
      end if
   end function rki_unstr_get_ncells_total

   function rki_unstr_get_mass_error(c_handle) result(e) &
      bind(c, name="rki_unstr_get_mass_error")
      type(c_ptr), value :: c_handle
      real(c_double) :: e
      type(solver_handle_unstr_t), pointer :: h
      integer :: status

      status = handle_check_unstr(c_handle, h)
      if (status == HANDLE_OK .and. h%is_initialised) then
         e = real(h%state%mass_error, c_double)
      else
         e = 0.0_c_double
      end if
   end function rki_unstr_get_mass_error

   ! ================================================================
   ! Setters (host-side, before finalize)
   ! ================================================================
   !
   ! Array setters take the buffer as a `c_ptr value` and bind it to a
   ! Fortran pointer with `c_f_pointer` using the cell count from the
   ! resolved handle. The C ABI is the same as `double*` — but `bed(*)`
   ! on a `bind(C)` dummy is flagged as bad style by Fortran linters
   ! (C071 / "assumed-size") because the compiler cannot bounds-check
   ! it. The pointer-then-c_f_pointer pattern avoids the warning while
   ! preserving the byte layout the caller already expects.

   function rki_unstr_set_bathymetry(c_handle, c_bed) result(ierr) &
      bind(c, name="rki_unstr_set_bathymetry")
      type(c_ptr), value :: c_handle
      type(c_ptr), value :: c_bed
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(c_double), pointer :: bed(:)
      integer :: k

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. h%is_finalised) then
         ierr = -1
         return
      end if
      call c_f_pointer(c_bed, bed, [h%state%ncells_owned])
      do k = 1, h%state%ncells_owned
         h%state%barotropic%bed_cv(k) = real(bed(k), wp)
      end do
      ierr = 0
   end function rki_unstr_set_bathymetry

   function rki_unstr_set_initial_stage(c_handle, c_stage) result(ierr) &
      bind(c, name="rki_unstr_set_initial_stage")
      type(c_ptr), value :: c_handle
      type(c_ptr), value :: c_stage
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(c_double), pointer :: stage_in(:)
      integer :: k

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. h%is_finalised) then
         ierr = -1
         return
      end if
      call c_f_pointer(c_stage, stage_in, [h%state%ncells_owned])
      do k = 1, h%state%ncells_owned
         h%state%barotropic%stage_cv(k) = real(stage_in(k), wp)
      end do
      ierr = 0
   end function rki_unstr_set_initial_stage

   function rki_unstr_set_initial_xmom(c_handle, c_xmom) result(ierr) &
      bind(c, name="rki_unstr_set_initial_xmom")
      type(c_ptr), value :: c_handle
      type(c_ptr), value :: c_xmom
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(c_double), pointer :: xm(:)
      integer :: k

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. h%is_finalised) then
         ierr = -1
         return
      end if
      call c_f_pointer(c_xmom, xm, [h%state%ncells_owned])
      do k = 1, h%state%ncells_owned
         h%state%barotropic%xmom_cv(k) = real(xm(k), wp)
      end do
      ierr = 0
   end function rki_unstr_set_initial_xmom

   function rki_unstr_set_initial_ymom(c_handle, c_ymom) result(ierr) &
      bind(c, name="rki_unstr_set_initial_ymom")
      type(c_ptr), value :: c_handle
      type(c_ptr), value :: c_ymom
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(c_double), pointer :: ym(:)
      integer :: k

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. h%is_finalised) then
         ierr = -1
         return
      end if
      call c_f_pointer(c_ymom, ym, [h%state%ncells_owned])
      do k = 1, h%state%ncells_owned
         h%state%barotropic%ymom_cv(k) = real(ym(k), wp)
      end do
      ierr = 0
   end function rki_unstr_set_initial_ymom

   function rki_unstr_set_bc_type(c_handle, group, bc_type) result(ierr) &
      bind(c, name="rki_unstr_set_bc_type")
      type(c_ptr), value :: c_handle
      integer(c_int), value :: group
      integer(c_int), value :: bc_type
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (group < 1 .or. group > MAX_BC_GROUPS) then
         ierr = -2
         return
      end if
      h%state%bc_type(group) = bc_type
      ierr = 0
   end function rki_unstr_set_bc_type

   function rki_unstr_set_bc_eta_ext(c_handle, group, value) result(ierr) &
      bind(c, name="rki_unstr_set_bc_eta_ext")
      type(c_ptr), value :: c_handle
      integer(c_int), value :: group
      real(c_double), value :: value
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (group < 1 .or. group > MAX_BC_GROUPS) then
         ierr = -2
         return
      end if
      h%state%bc_eta_ext(group) = real(value, wp)
      ierr = 0
   end function rki_unstr_set_bc_eta_ext

   function rki_unstr_set_bc_un_ext(c_handle, group, value) result(ierr) &
      bind(c, name="rki_unstr_set_bc_un_ext")
      type(c_ptr), value :: c_handle
      integer(c_int), value :: group
      real(c_double), value :: value
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (group < 1 .or. group > MAX_BC_GROUPS) then
         ierr = -2
         return
      end if
      h%state%bc_un_ext(group) = real(value, wp)
      ierr = 0
   end function rki_unstr_set_bc_un_ext

   function rki_unstr_set_bc_h_ext(c_handle, group, value) result(ierr) &
      bind(c, name="rki_unstr_set_bc_h_ext")
      type(c_ptr), value :: c_handle
      integer(c_int), value :: group
      real(c_double), value :: value
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (group < 1 .or. group > MAX_BC_GROUPS) then
         ierr = -2
         return
      end if
      h%state%bc_h_ext(group) = real(value, wp)
      ierr = 0
   end function rki_unstr_set_bc_h_ext

   function rki_unstr_set_bc_qn_ext(c_handle, group, value) result(ierr) &
      bind(c, name="rki_unstr_set_bc_qn_ext")
      type(c_ptr), value :: c_handle
      integer(c_int), value :: group
      real(c_double), value :: value
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (group < 1 .or. group > MAX_BC_GROUPS) then
         ierr = -2
         return
      end if
      h%state%bc_Qn_ext(group) = real(value, wp)
      ierr = 0
   end function rki_unstr_set_bc_qn_ext

   function rki_unstr_set_manning(c_handle, n_value) result(ierr) &
      bind(c, name="rki_unstr_set_manning")
      type(c_ptr), value :: c_handle
      real(c_double), value :: n_value
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%barotropic%manning_n = real(n_value, wp)
      ierr = 0
   end function rki_unstr_set_manning

   function rki_unstr_set_time_scheme(c_handle, c_name, c_len) result(ierr) &
      bind(c, name="rki_unstr_set_time_scheme")
      !! Select the explicit time integrator on the unstructured 2D
      !! barotropic path.  `c_name` is one of "euler"/"rk2"/"rk3"/"rk43"
      !! (unrecognised falls back to "rk2").  Returns:
      !!   0  on success,
      !!  -1  if the handle is not initialised,
      !!  -2  if a non-rk2 scheme is requested on a multilayer handle
      !!      (multilayer keeps SSP-RK2; the state field is left unchanged).
      type(c_ptr), value :: c_handle
      type(c_ptr), value :: c_name
      integer(c_int), value :: c_len
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      character(len=:), allocatable :: name
      character, pointer :: cstr(:)
      integer :: i, scheme

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if

      call c_f_pointer(c_name, cstr, [c_len])
      allocate (character(len=c_len) :: name)
      do i = 1, c_len
         name(i:i) = cstr(i)
      end do

      scheme = parse_time_scheme(name)
      if (scheme /= TIME_RK2 .and. h%state%use_multilayer) then
         ierr = -2
         return
      end if
      h%state%time_scheme = scheme
      ierr = 0
   end function rki_unstr_set_time_scheme

   function rki_unstr_set_eos(c_handle, rho_0, beta_s, s_ref, alpha_t, t_ref) &
      result(ierr) bind(c, name="rki_unstr_set_eos")
      type(c_ptr), value :: c_handle
      real(c_double), value :: rho_0, beta_s, s_ref, alpha_t, t_ref
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%multilayer%rho_0 = real(rho_0, wp)
      h%state%multilayer%beta_S = real(beta_s, wp)
      h%state%multilayer%S_ref = real(s_ref, wp)
      h%state%multilayer%alpha_T = real(alpha_t, wp)
      h%state%multilayer%T_ref = real(t_ref, wp)
      ierr = 0
   end function rki_unstr_set_eos

   function rki_unstr_set_bpg_jacobian(c_handle, use_jacobian) result(ierr) &
      bind(c, name="rki_unstr_set_bpg_jacobian")
      type(c_ptr), value :: c_handle
      integer(c_int), value :: use_jacobian
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%multilayer%use_bpg_jacobian = (use_jacobian /= 0)
      ierr = 0
   end function rki_unstr_set_bpg_jacobian

   function rki_unstr_set_vcoord(c_handle, coord_type, remap_method, &
                                 depth_transition, blend_width) result(ierr) &
      bind(c, name="rki_unstr_set_vcoord")
      type(c_ptr), value :: c_handle
      integer(c_int), value :: coord_type, remap_method
      real(c_double), value :: depth_transition, blend_width
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. .not. h%state%use_multilayer) then
         ierr = -1
         return
      end if
      h%state%multilayer%vcoord%coord_type = coord_type
      h%state%multilayer%vcoord%remap_method = remap_method
      h%state%multilayer%vcoord%depth_transition = real(depth_transition, wp)
      h%state%multilayer%vcoord%blend_width = real(blend_width, wp)
      call h%state%multilayer%vcoord%init(h%state%multilayer%nz_ml, coord_type, remap_method)
      if (allocated(h%state%multilayer%vcoord%z_ref) .and. allocated(h%state%multilayer%vcoord_z_ref)) then
         h%state%multilayer%vcoord_z_ref = h%state%multilayer%vcoord%z_ref
      end if
      if (h%state%multilayer%vcoord%needs_remap() .and. .not. allocated(h%state%multilayer%remap_dz_old)) then
         allocate (h%state%multilayer%remap_dz_old(h%state%multilayer%nz_ml, h%state%ncells_total), source=0.0_wp)
      end if
      ! VCOORD_ZSTAR_FULL needs a per-cell reference profile.  This
      ! allocation is lazy because the Python flow typically calls
      ! `set_vcoord` AFTER `enable_multilayer` (which allocated state
      ! with the default sigma coord_type, so vcoord_z_ref_col was
      ! skipped at that time).
      if (coord_type == VCOORD_ZSTAR_FULL .and. &
          .not. allocated(h%state%multilayer%vcoord_z_ref_col)) then
         allocate (h%state%multilayer%vcoord_z_ref_col(0:h%state%multilayer%nz_ml, &
                                                       h%state%ncells_total), &
                   source=0.0_wp)
      end if
      ierr = 0
   end function rki_unstr_set_vcoord

   function rki_unstr_set_zstar_full_params(c_handle, h_surf_target, &
                                            h_min, n_surf, stretching) &
      result(ierr) bind(c, name="rki_unstr_set_zstar_full_params")
      !! Set VCOORD_ZSTAR_FULL parameters: surface-layer target thickness
      !! (m), vanishing-layer floor (m), number of fine surface layers
      !! (0 = auto), and stretching mode (0=uniform, 1=log).
      !!
      !! Must be called after `rki_unstr_set_vcoord(..., 'zstar_full')`
      !! and before `rki_unstr_finalize`.
      type(c_ptr), value :: c_handle
      real(c_double), value :: h_surf_target, h_min
      integer(c_int), value :: n_surf, stretching
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. .not. h%state%use_multilayer) then
         ierr = -1
         return
      end if
      h%state%multilayer%vcoord%zstar_h_surf_target = real(h_surf_target, wp)
      h%state%multilayer%vcoord%zstar_h_min = real(h_min, wp)
      h%state%multilayer%vcoord%zstar_n_surf = n_surf
      h%state%multilayer%vcoord%zstar_stretching = stretching
      ierr = 0
   end function rki_unstr_set_zstar_full_params

   function rki_unstr_set_cfl(c_handle, cfl_value) result(ierr) &
      bind(c, name="rki_unstr_set_cfl")
      type(c_ptr), value :: c_handle
      real(c_double), value :: cfl_value
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%cfl = real(cfl_value, wp)
      ierr = 0
   end function rki_unstr_set_cfl

   function rki_unstr_set_betas(c_handle, beta_w, beta_uh, beta_vh) result(ierr) &
      bind(c, name="rki_unstr_set_betas")
      type(c_ptr), value :: c_handle
      real(c_double), value :: beta_w, beta_uh, beta_vh
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%beta_w = real(beta_w, wp)
      h%state%beta_w_dry = real(beta_w, wp)
      h%state%beta_uh = real(beta_uh, wp)
      h%state%beta_uh_dry = real(beta_uh, wp)
      h%state%beta_vh = real(beta_vh, wp)
      h%state%beta_vh_dry = real(beta_vh, wp)
      ierr = 0
   end function rki_unstr_set_betas

   function rki_unstr_set_wind(c_handle, tau_x, tau_y) result(ierr) &
      bind(c, name="rki_unstr_set_wind")
      type(c_ptr), value :: c_handle
      real(c_double), value :: tau_x, tau_y
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%barotropic%wind_tau_x = real(tau_x, wp)
      h%state%barotropic%wind_tau_y = real(tau_y, wp)
      ierr = 0
   end function rki_unstr_set_wind

   function rki_unstr_set_heat_flux(c_handle, c_q) result(ierr) &
      bind(c, name="rki_unstr_set_heat_flux")
      !! Set net surface heat flux (W/m^2, positive into ocean) on every
      !! cell.  Length of `c_q` must equal `state%ncells_owned`.
      type(c_ptr), value :: c_handle
      type(c_ptr), value :: c_q
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(c_double), pointer :: q(:)
      integer :: k

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. h%is_finalised) then
         ierr = -1
         return
      end if
      call c_f_pointer(c_q, q, [h%state%ncells_owned])
      do k = 1, h%state%ncells_owned
         h%state%forcing%heat_flux(k) = real(q(k), wp)
      end do
      if (h%on_device) call solver_unstr_update_forcing(h%state)
      ierr = 0
   end function rki_unstr_set_heat_flux

   function rki_unstr_set_evap_minus_precip(c_handle, c_ep) result(ierr) &
      bind(c, name="rki_unstr_set_evap_minus_precip")
      !! Set virtual surface freshwater flux (m/s, positive for net
      !! evaporation) on every owned cell.  Same sizing convention as
      !! `rki_unstr_set_heat_flux`.
      type(c_ptr), value :: c_handle
      type(c_ptr), value :: c_ep
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(c_double), pointer :: ep(:)
      integer :: k

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. h%is_finalised) then
         ierr = -1
         return
      end if
      call c_f_pointer(c_ep, ep, [h%state%ncells_owned])
      do k = 1, h%state%ncells_owned
         h%state%forcing%evap_minus_precip(k) = real(ep(k), wp)
      end do
      if (h%on_device) call solver_unstr_update_forcing(h%state)
      ierr = 0
   end function rki_unstr_set_evap_minus_precip

   function rki_unstr_set_wind_field(c_handle, c_wx, c_wy) result(ierr) &
      bind(c, name="rki_unstr_set_wind_field")
      !! Set per-cell surface wind STRESS (Pa), x and y components, on every
      !! owned cell.  Added to the barotropic momentum in the fused update
      !! kernel on top of any spatially-uniform stress from
      !! `rki_unstr_set_wind`.  Same sizing convention as
      !! `rki_unstr_set_heat_flux` (length == state%ncells_owned).
      !! Marks the per-cell wind active so the
      !! kernel applies it; the caller converts wind speed to stress.
      type(c_ptr), value :: c_handle
      type(c_ptr), value :: c_wx, c_wy
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(c_double), pointer :: wx(:), wy(:)
      integer :: k

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised .or. h%is_finalised) then
         ierr = -1
         return
      end if
      call c_f_pointer(c_wx, wx, [h%state%ncells_owned])
      call c_f_pointer(c_wy, wy, [h%state%ncells_owned])
      do k = 1, h%state%ncells_owned
         h%state%forcing%wind_x(k) = real(wx(k), wp)
         h%state%forcing%wind_y(k) = real(wy(k), wp)
      end do
      h%state%forcing%wind_active = .true.
      if (h%on_device) call solver_unstr_update_forcing(h%state)
      ierr = 0
   end function rki_unstr_set_wind_field

   function rki_unstr_set_initial_salinity(c_handle, value) result(ierr) &
      bind(c, name="rki_unstr_set_initial_salinity")
      type(c_ptr), value :: c_handle
      real(c_double), value :: value
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%multilayer%initial_salinity = real(value, wp)
      ierr = 0
   end function rki_unstr_set_initial_salinity

   function rki_unstr_set_initial_temperature(c_handle, value) result(ierr) &
      bind(c, name="rki_unstr_set_initial_temperature")
      type(c_ptr), value :: c_handle
      real(c_double), value :: value
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%multilayer%initial_temperature = real(value, wp)
      ierr = 0
   end function rki_unstr_set_initial_temperature

   function rki_unstr_set_bc_salinity(c_handle, group, value) result(ierr) &
      bind(c, name="rki_unstr_set_bc_salinity")
      type(c_ptr), value :: c_handle
      integer(c_int), value :: group
      real(c_double), value :: value
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (group < 1 .or. group > MAX_BC_GROUPS) then
         ierr = -2
         return
      end if
      h%state%multilayer%bc_S_ext(group) = real(value, wp)
      ierr = 0
   end function rki_unstr_set_bc_salinity

   function rki_unstr_set_bc_temperature(c_handle, group, value) result(ierr) &
      bind(c, name="rki_unstr_set_bc_temperature")
      type(c_ptr), value :: c_handle
      integer(c_int), value :: group
      real(c_double), value :: value
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (group < 1 .or. group > MAX_BC_GROUPS) then
         ierr = -2
         return
      end if
      h%state%multilayer%bc_T_ext(group) = real(value, wp)
      ierr = 0
   end function rki_unstr_set_bc_temperature

   function rki_unstr_set_smagorinsky_cs(c_handle, value) result(ierr) &
      bind(c, name="rki_unstr_set_smagorinsky_cs")
      type(c_ptr), value :: c_handle
      real(c_double), value :: value
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%multilayer%smagorinsky_cs = real(value, wp)
      ierr = 0
   end function rki_unstr_set_smagorinsky_cs

   function rki_unstr_set_hdiff_kappa(c_handle, value) result(ierr) &
      bind(c, name="rki_unstr_set_hdiff_kappa")
      type(c_ptr), value :: c_handle
      real(c_double), value :: value
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%multilayer%hdiff_kappa = real(value, wp)
      ierr = 0
   end function rki_unstr_set_hdiff_kappa

   function rki_unstr_set_use_ppm_tracer(c_handle, flag) result(ierr) &
      bind(c, name="rki_unstr_set_use_ppm_tracer")
      type(c_ptr), value :: c_handle
      integer(c_int), value :: flag
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%multilayer%use_ppm_tracer = (flag /= 0)
      ierr = 0
   end function rki_unstr_set_use_ppm_tracer

   function rki_unstr_set_vmix(c_handle, nu_max, nu_bg, ri_alpha, &
                               kappa_s_bg, kappa_t_bg) &
      result(ierr) bind(c, name="rki_unstr_set_vmix")
      type(c_ptr), value :: c_handle
      real(c_double), value :: nu_max, nu_bg, ri_alpha, kappa_s_bg, kappa_t_bg
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%multilayer%vmix_nu_max = real(nu_max, wp)
      h%state%multilayer%vmix_nu_bg = real(nu_bg, wp)
      h%state%multilayer%vmix_ri_alpha = real(ri_alpha, wp)
      h%state%multilayer%kappa_S_bg = real(kappa_s_bg, wp)
      h%state%multilayer%kappa_T_bg = real(kappa_t_bg, wp)
      ierr = 0
   end function rki_unstr_set_vmix

   function rki_unstr_set_kpp(c_handle, enable, ri_crit, cs_nonlocal, c_vt2) &
      result(ierr) bind(c, name="rki_unstr_set_kpp")
      !! Enable/configure the KPP surface-boundary-layer overlay on top
      !! of the PP81 interior closure.  The per-cell KPP workspace fields
      !! are allocated unconditionally at init, so this only flips the
      !! switch and sets the scalar coefficients (no device sync needed —
      !! the scalars ride over with the parent struct).
      type(c_ptr), value :: c_handle
      integer(c_int), value :: enable
      real(c_double), value :: ri_crit, cs_nonlocal, c_vt2
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      h%state%multilayer%vmix_use_kpp = (enable /= 0)
      h%state%multilayer%kpp_ri_crit = real(ri_crit, wp)
      h%state%multilayer%kpp_cs_nonlocal = real(cs_nonlocal, wp)
      h%state%multilayer%kpp_c_vt2 = real(c_vt2, wp)
      ierr = 0
   end function rki_unstr_set_kpp

   ! ================================================================
   ! Field getters (zero-copy, host pointers)
   ! ================================================================

   function rki_unstr_get_stage_ptr(c_handle, ptr, n) result(ierr) &
      bind(c, name="rki_unstr_get_stage_ptr")
      type(c_ptr), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: n
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(wp), pointer :: tmp(:)

      ptr = c_null_ptr
      n = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      call refresh_host(h)
      tmp => h%state%barotropic%stage_cv
      ptr = c_loc(tmp(1))
      n = h%state%ncells_owned
      ierr = 0
   end function rki_unstr_get_stage_ptr

   function rki_unstr_get_xmom_ptr(c_handle, ptr, n) result(ierr) &
      bind(c, name="rki_unstr_get_xmom_ptr")
      type(c_ptr), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: n
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(wp), pointer :: tmp(:)

      ptr = c_null_ptr
      n = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      call refresh_host(h)
      tmp => h%state%barotropic%xmom_cv
      ptr = c_loc(tmp(1))
      n = h%state%ncells_owned
      ierr = 0
   end function rki_unstr_get_xmom_ptr

   function rki_unstr_get_ymom_ptr(c_handle, ptr, n) result(ierr) &
      bind(c, name="rki_unstr_get_ymom_ptr")
      type(c_ptr), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: n
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(wp), pointer :: tmp(:)

      ptr = c_null_ptr
      n = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      call refresh_host(h)
      tmp => h%state%barotropic%ymom_cv
      ptr = c_loc(tmp(1))
      n = h%state%ncells_owned
      ierr = 0
   end function rki_unstr_get_ymom_ptr

   function rki_unstr_get_bed_ptr(c_handle, ptr, n) result(ierr) &
      bind(c, name="rki_unstr_get_bed_ptr")
      type(c_ptr), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: n
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(wp), pointer :: tmp(:)

      ptr = c_null_ptr
      n = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      tmp => h%state%barotropic%bed_cv
      ptr = c_loc(tmp(1))
      n = h%state%ncells_owned
      ierr = 0
   end function rki_unstr_get_bed_ptr

   function rki_unstr_get_height_ptr(c_handle, ptr, n) result(ierr) &
      bind(c, name="rki_unstr_get_height_ptr")
      type(c_ptr), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: n
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(wp), pointer :: tmp(:)

      ptr = c_null_ptr
      n = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      call refresh_host(h)
      tmp => h%state%barotropic%height_cv
      ptr = c_loc(tmp(1))
      n = h%state%ncells_owned
      ierr = 0
   end function rki_unstr_get_height_ptr

   function rki_unstr_get_h_layer_ptr(c_handle, ptr, nz, ncells) result(ierr) &
      bind(c, name="rki_unstr_get_h_layer_ptr")
      type(c_ptr), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nz, ncells
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nz = 0
      ncells = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (.not. h%state%use_multilayer) then
         ierr = -2
         return
      end if
      call refresh_host(h)
      tmp => h%state%multilayer%h_layer
      ptr = c_loc(tmp(1, 1))
      nz = h%state%multilayer%nz_ml
      ncells = h%state%ncells_owned
      ierr = 0
   end function rki_unstr_get_h_layer_ptr

   function rki_unstr_get_hu_layer_ptr(c_handle, ptr, nz, ncells) result(ierr) &
      bind(c, name="rki_unstr_get_hu_layer_ptr")
      type(c_ptr), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nz, ncells
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nz = 0
      ncells = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (.not. h%state%use_multilayer) then
         ierr = -2
         return
      end if
      call refresh_host(h)
      tmp => h%state%multilayer%hu_layer
      ptr = c_loc(tmp(1, 1))
      nz = h%state%multilayer%nz_ml
      ncells = h%state%ncells_owned
      ierr = 0
   end function rki_unstr_get_hu_layer_ptr

   function rki_unstr_get_hv_layer_ptr(c_handle, ptr, nz, ncells) result(ierr) &
      bind(c, name="rki_unstr_get_hv_layer_ptr")
      type(c_ptr), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nz, ncells
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nz = 0
      ncells = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (.not. h%state%use_multilayer) then
         ierr = -2
         return
      end if
      call refresh_host(h)
      tmp => h%state%multilayer%hv_layer
      ptr = c_loc(tmp(1, 1))
      nz = h%state%multilayer%nz_ml
      ncells = h%state%ncells_owned
      ierr = 0
   end function rki_unstr_get_hv_layer_ptr

   function rki_unstr_get_hS_layer_ptr(c_handle, ptr, nz, ncells) result(ierr) &
      bind(c, name="rki_unstr_get_hS_layer_ptr")
      type(c_ptr), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nz, ncells
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nz = 0
      ncells = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (.not. h%state%use_multilayer) then
         ierr = -2
         return
      end if
      call refresh_host(h)
      tmp => h%state%multilayer%tracers(h%state%multilayer%idx_salinity)%hTr
      ptr = c_loc(tmp(1, 1))
      nz = h%state%multilayer%nz_ml
      ncells = h%state%ncells_owned
      ierr = 0
   end function rki_unstr_get_hS_layer_ptr

   function rki_unstr_get_hT_layer_ptr(c_handle, ptr, nz, ncells) result(ierr) &
      bind(c, name="rki_unstr_get_hT_layer_ptr")
      type(c_ptr), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nz, ncells
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nz = 0
      ncells = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (.not. h%state%use_multilayer) then
         ierr = -2
         return
      end if
      call refresh_host(h)
      tmp => h%state%multilayer%tracers(h%state%multilayer%idx_temperature)%hTr
      ptr = c_loc(tmp(1, 1))
      nz = h%state%multilayer%nz_ml
      ncells = h%state%ncells_owned
      ierr = 0
   end function rki_unstr_get_hT_layer_ptr

   function rki_unstr_get_centroid_x_ptr(c_handle, ptr, n) result(ierr) &
      bind(c, name="rki_unstr_get_centroid_x_ptr")
      type(c_ptr), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: n
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(wp), pointer :: tmp(:)

      ptr = c_null_ptr
      n = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      tmp => h%mesh%centroid_x
      ptr = c_loc(tmp(1))
      n = h%mesh%ncells_owned
      ierr = 0
   end function rki_unstr_get_centroid_x_ptr

   function rki_unstr_get_centroid_y_ptr(c_handle, ptr, n) result(ierr) &
      bind(c, name="rki_unstr_get_centroid_y_ptr")
      type(c_ptr), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: n
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(wp), pointer :: tmp(:)

      ptr = c_null_ptr
      n = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      tmp => h%mesh%centroid_y
      ptr = c_loc(tmp(1))
      n = h%mesh%ncells_owned
      ierr = 0
   end function rki_unstr_get_centroid_y_ptr

   function rki_unstr_get_area_ptr(c_handle, ptr, n) result(ierr) &
      bind(c, name="rki_unstr_get_area_ptr")
      type(c_ptr), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: n
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h
      real(wp), pointer :: tmp(:)

      ptr = c_null_ptr
      n = 0
      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      tmp => h%mesh%area
      ptr = c_loc(tmp(1))
      n = h%mesh%ncells_owned
      ierr = 0
   end function rki_unstr_get_area_ptr

   ! ================================================================
   ! Restart files
   ! ================================================================

   function rki_unstr_save_restart(c_handle, c_filename, c_len) result(ierr) &
      bind(c, name="rki_unstr_save_restart")
      type(c_ptr), value :: c_handle
      type(c_ptr), value :: c_filename
      integer(c_int), value :: c_len
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      character(len=:), allocatable :: fname
      character, pointer :: cstr(:)
      integer :: i

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if

      call c_f_pointer(c_filename, cstr, [c_len])
      allocate (character(len=c_len) :: fname)
      do i = 1, c_len
         fname(i:i) = cstr(i)
      end do

      if (h%on_device) call solver_unstr_update_host(h%state)

      call restart_write_unstr(fname, h%state, h%mesh, h%state%time)
      ierr = 0
   end function rki_unstr_save_restart

   function rki_unstr_load_restart(c_handle, c_filename, c_len) result(ierr) &
      bind(c, name="rki_unstr_load_restart")
      type(c_ptr), value :: c_handle
      type(c_ptr), value :: c_filename
      integer(c_int), value :: c_len
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      character(len=:), allocatable :: fname
      character, pointer :: cstr(:)
      integer :: i, read_ierr
      real(wp) :: t_restart

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

      call c_f_pointer(c_filename, cstr, [c_len])
      allocate (character(len=c_len) :: fname)
      do i = 1, c_len
         fname(i:i) = cstr(i)
      end do

      call restart_read_unstr(fname, h%state, t_restart, read_ierr)
      if (read_ierr /= 0) then
         ierr = -3
         return
      end if

      ierr = 0
   end function rki_unstr_load_restart

   ! ================================================================
   ! UGRID NetCDF output
   ! ================================================================

   function rki_unstr_output_init(c_handle, filename, namelen, deflate_level) &
      result(ierr) bind(c, name="rki_unstr_output_init")
      type(c_ptr), value :: c_handle
      type(c_ptr), value :: filename
      integer(c_int), value :: namelen
      integer(c_int), value :: deflate_level
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      character(len=:), allocatable :: fname
      character, pointer :: cstr(:)
      integer :: i

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (h%output_open) then
         ierr = -2
         return
      end if

      if (h%on_device) then
         call solver_unstr_update_host(h%state)
      end if

      call c_f_pointer(filename, cstr, [namelen])
      allocate (character(len=namelen) :: fname)
      do i = 1, namelen
         fname(i:i) = cstr(i)
      end do
      call output_unstr_init(h%ostate, fname, h%state, h%mesh, deflate_level)

      h%output_open = .true.
      ierr = 0
   end function rki_unstr_output_init

   function rki_unstr_output_write(c_handle) result(ierr) &
      bind(c, name="rki_unstr_output_write")
      type(c_ptr), value :: c_handle
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%output_open) then
         ierr = -1
         return
      end if

      if (h%on_device) then
         call solver_unstr_update_host(h%state)
      end if

      call output_unstr_write_snapshot(h%ostate, h%state, h%state%time)

      ierr = 0
   end function rki_unstr_output_write

   function rki_unstr_output_close(c_handle) result(ierr) &
      bind(c, name="rki_unstr_output_close")
      type(c_ptr), value :: c_handle
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%output_open) then
         ierr = -1
         return
      end if

      call output_unstr_close(h%ostate)

      h%output_open = .false.
      ierr = 0
   end function rki_unstr_output_close

   ! ================================================================
   ! Gauge output
   ! ================================================================

   function rki_unstr_gauges_init(c_handle, c_filename, c_len, n_gauges, gx, gy) &
      result(ierr) bind(c, name="rki_unstr_gauges_init")
      type(c_ptr), value :: c_handle
      type(c_ptr), value :: c_filename
      integer(c_int), value :: c_len
      integer(c_int), value :: n_gauges
      real(c_double), intent(in) :: gx(n_gauges)
      real(c_double), intent(in) :: gy(n_gauges)
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      character(len=:), allocatable :: fname
      character, pointer :: cstr(:)
      integer :: i

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%is_initialised) then
         ierr = -1
         return
      end if
      if (h%gauges_open) then
         ierr = -2
         return
      end if

      call c_f_pointer(c_filename, cstr, [c_len])
      allocate (character(len=c_len) :: fname)
      do i = 1, c_len
         fname(i:i) = cstr(i)
      end do

      call gauges_unstr_init(h%gauges, fname, h%mesh, int(n_gauges), &
                             real(gx, wp), real(gy, wp))
      h%gauges_open = .true.
      ierr = 0
   end function rki_unstr_gauges_init

   function rki_unstr_gauges_write(c_handle) result(ierr) &
      bind(c, name="rki_unstr_gauges_write")
      type(c_ptr), value :: c_handle
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%gauges_open) then
         ierr = -1
         return
      end if

      if (h%on_device) call solver_unstr_update_host(h%state)

      call gauges_unstr_write(h%gauges, h%state, h%state%time)
      ierr = 0
   end function rki_unstr_gauges_write

   function rki_unstr_gauges_close(c_handle) result(ierr) &
      bind(c, name="rki_unstr_gauges_close")
      type(c_ptr), value :: c_handle
      integer(c_int) :: ierr
      type(solver_handle_unstr_t), pointer :: h

      ierr = resolve(c_handle, h)
      if (ierr /= 0) return
      if (.not. h%gauges_open) then
         ierr = -1
         return
      end if

      call gauges_unstr_close(h%gauges)
      h%gauges_open = .false.
      ierr = 0
   end function rki_unstr_gauges_close

end module rki_api_unstr
