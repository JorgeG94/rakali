!! NetCDF solution output and restart file I/O
module rki_output
   !! Writes solution snapshots to NetCDF with CF conventions and
   !! provides restart file read/write for continuing simulations.
   use rki_constants, only: wp, DRY_TOLERANCE
   use rki_config, only: config_t
   use rki_state, only: state_t
   use rki_io_netcdf, only: nc_check, nc_create_file, nc_open_read, &
                            nc_open_write, nc_close, nc_def_dim, nc_def_var_2d, nc_def_var_3d, &
                            nc_def_var_4d, nc_put_att, nc_put_att_real, nc_put_att_global, &
                            nc_enddef, nc_put_var_2d, &
                            nc_put_var_3d_slice, &
                            nc_put_var_4d_slice, nc_get_var_2d, nc_get_dim_len, &
                            nc_get_varid, NC_WP
   use netcdf, only: nf90_unlimited, nf90_def_var, nf90_put_var, &
                     nf90_get_var, nf90_inq_varid, nf90_redef, nf90_get_att, nf90_noerr, &
                     nf90_inquire_dimension, nf90_inq_dimid, nf90_sync, nf90_double, &
                     nf90_fill_double, nf90_fill_float
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: output_state_t
   public :: output_init
   public :: output_write_snapshot
   public :: output_close
   public :: restart_write
   public :: restart_read
   public :: output_rank_filename
   public :: output_write_snapshot_from_arrays
   public :: write_ascii_output

#ifdef RAKALI_DOUBLE_PRECISION
   real(wp), parameter :: FILL_VALUE = nf90_fill_double
      !! NetCDF fill value matching working precision
#else
   real(wp), parameter :: FILL_VALUE = nf90_fill_float
      !! NetCDF fill value matching working precision
#endif

   type :: output_state_t
      !! Tracks the state of an open output NetCDF file
      integer :: ncid = -1
         !! NetCDF file ID (-1 = not open)
      integer :: time_dimid = -1
         !! Time dimension ID
      integer :: time_varid = -1
         !! Time variable ID
      integer :: h_varid = -1
         !! Water depth variable ID
      integer :: hu_varid = -1
         !! x-momentum variable ID
      integer :: hv_varid = -1
         !! y-momentum variable ID
      integer :: time_index = 0
         !! Current time record index

      ! Non-hydrostatic output
      logical :: has_nh = .false.
         !! Whether NH fields are present
      integer :: nz = 0
         !! Number of sigma layers (cached for output)
      integer :: u_layer_varid = -1
         !! Per-layer x-velocity variable ID
      integer :: v_layer_varid = -1
         !! Per-layer y-velocity variable ID
      integer :: w_varid = -1
         !! Vertical velocity variable ID
      integer :: q_nh_varid = -1
         !! NH pressure variable ID
      integer :: vort_x_varid = -1
         !! x-vorticity (dw/dy - dv/dz) variable ID
      integer :: vort_y_varid = -1
         !! y-vorticity (du/dz - dw/dx) variable ID
      real(wp), allocatable :: dsig(:)

      ! Multilayer output
      logical :: has_ml = .false.
         !! Whether multilayer fields are present
      integer :: nz_ml = 0
         !! Number of multilayer sigma layers
      integer :: ml_h_layer_varid = -1
         !! Layer thickness variable ID
      integer :: ml_u_varid = -1
         !! Layer x-velocity variable ID
      integer :: ml_v_varid = -1
         !! Layer y-velocity variable ID
      integer, allocatable :: ml_tracer_varids(:)
         !! Per-tracer NetCDF variable IDs, one entry per
         !! `state%multilayer%tracers(:)` slot.  Set to `-1` for
         !! tracers whose output is gated off; otherwise carries the
         !! varid of the corresponding NetCDF variable.  Replaces the
         !! historical `ml_salinity_varid` / `ml_temperature_varid`
         !! pair so adding a passive tracer requires no I/O changes.
      integer :: ml_rho_varid = -1
         !! Layer density variable ID
      integer :: ml_w_varid = -1
         !! Vertical velocity at interfaces variable ID
      integer :: bed_sediment_varid = -1
         !! Bed sediment mass reservoir variable ID (set only with sediment on)
         !! Cached sigma layer thicknesses for vorticity computation
   end type output_state_t

contains

   subroutine output_init(ostate, filename, state, cfg)
      !! Create and initialise a NetCDF output file with CF conventions
      type(output_state_t), intent(inout) :: ostate
         !! Output file state (populated on return)
      character(len=*), intent(in) :: filename
         !! Path to output file
      type(state_t), intent(in) :: state
         !! Simulation state (for grid metadata)
      type(config_t), intent(in) :: cfg
         !! Runtime configuration (for CF metadata)

      integer :: x_dimid, y_dimid, x_varid, y_varid, b_varid
      integer :: dimids_2d(2), dimids_3d(3)
      integer :: z_dimid, zi_dimid, dsig_varid
      integer :: dimids_4d(4), dimids_4d_w(4)
      integer :: i, j, k, ng
      real(wp), allocatable :: x_coords(:), y_coords(:)
      character(len=256) :: time_units, history_str
      integer, parameter :: DT_SIZE = 8
      integer :: dt_values(DT_SIZE)
      integer :: dlevel
         !! Deflate level (0 = no compression)

      call logger%info("Creating output file: "//trim(filename))

      call nc_create_file(filename, ostate%ncid)

      ! Compression level (0 = disabled)
      dlevel = 0
      if (cfg%compress_output) dlevel = cfg%compress_level

      ! Define dimensions
      call nc_def_dim(ostate%ncid, "x", state%grid%nx_phys, x_dimid)
      call nc_def_dim(ostate%ncid, "y", state%grid%ny_phys, y_dimid)
      call nc_def_dim(ostate%ncid, "time", nf90_unlimited, ostate%time_dimid)

      dimids_2d = [x_dimid, y_dimid]
      dimids_3d = [x_dimid, y_dimid, ostate%time_dimid]

      ! Define coordinate variables
      call nc_check(nf90_def_var(ostate%ncid, "x", NC_WP, [x_dimid], x_varid), &
                    "defining x coordinate")
      call nc_put_att(ostate%ncid, x_varid, "units", "m")
      call nc_put_att(ostate%ncid, x_varid, "long_name", "x-coordinate of cell centre")
      call nc_put_att(ostate%ncid, x_varid, "axis", "X")

      call nc_check(nf90_def_var(ostate%ncid, "y", NC_WP, [y_dimid], y_varid), &
                    "defining y coordinate")
      call nc_put_att(ostate%ncid, y_varid, "units", "m")
      call nc_put_att(ostate%ncid, y_varid, "long_name", "y-coordinate of cell centre")
      call nc_put_att(ostate%ncid, y_varid, "axis", "Y")

      ! Time variable — reference time depends on whether this is an
      ! idealised run (reference_time="") or a real-calendar run.
      call nc_check(nf90_def_var(ostate%ncid, "time", nf90_double, &
                                 [ostate%time_dimid], ostate%time_varid), &
                    "defining time variable")
      if (len_trim(cfg%reference_time) > 0) then
         time_units = "seconds since "//trim(cfg%reference_time)
         call nc_put_att(ostate%ncid, ostate%time_varid, "units", trim(time_units))
         call nc_put_att(ostate%ncid, ostate%time_varid, "calendar", "proleptic_gregorian")
      else
         call nc_put_att(ostate%ncid, ostate%time_varid, "units", &
                         "seconds since simulation start")
         call nc_put_att(ostate%ncid, ostate%time_varid, "calendar", "none")
      end if
      call nc_put_att(ostate%ncid, ostate%time_varid, "long_name", "time")
      call nc_put_att(ostate%ncid, ostate%time_varid, "axis", "T")

      ! Bathymetry (static, 2D)
      call nc_def_var_2d(ostate%ncid, "b", dimids_2d, b_varid, dlevel)
      call nc_put_att(ostate%ncid, b_varid, "units", "m")
      call nc_put_att(ostate%ncid, b_varid, "long_name", "bottom elevation (positive up)")
      call nc_put_att(ostate%ncid, b_varid, "coordinates", "x y")

      ! Solution variables (3D: x, y, time)
      call nc_def_var_3d(ostate%ncid, "h", dimids_3d, ostate%h_varid, dlevel)
      call nc_put_att(ostate%ncid, ostate%h_varid, "units", "m")
      call nc_put_att(ostate%ncid, ostate%h_varid, "long_name", "water depth")
      call nc_put_att_real(ostate%ncid, ostate%h_varid, "_FillValue", FILL_VALUE)
      call nc_put_att(ostate%ncid, ostate%h_varid, "coordinates", "x y")

      call nc_def_var_3d(ostate%ncid, "hu", dimids_3d, ostate%hu_varid, dlevel)
      call nc_put_att(ostate%ncid, ostate%hu_varid, "units", "m2 s-1")
      call nc_put_att(ostate%ncid, ostate%hu_varid, "long_name", "x-momentum (h*u)")
      call nc_put_att_real(ostate%ncid, ostate%hu_varid, "_FillValue", FILL_VALUE)
      call nc_put_att(ostate%ncid, ostate%hu_varid, "coordinates", "x y")

      call nc_def_var_3d(ostate%ncid, "hv", dimids_3d, ostate%hv_varid, dlevel)
      call nc_put_att(ostate%ncid, ostate%hv_varid, "units", "m2 s-1")
      call nc_put_att(ostate%ncid, ostate%hv_varid, "long_name", "y-momentum (h*v)")
      call nc_put_att_real(ostate%ncid, ostate%hv_varid, "_FillValue", FILL_VALUE)
      call nc_put_att(ostate%ncid, ostate%hv_varid, "coordinates", "x y")

      ! Non-hydrostatic variables (4D: x, y, z, time)
      if (state%use_nonhydrostatic) then
         ostate%has_nh = .true.
         ostate%nz = state%nh%nz
         allocate (ostate%dsig(state%nh%nz))
         ostate%dsig = state%nh%dsig

         call nc_def_dim(ostate%ncid, "z", state%nh%nz, z_dimid)
         call nc_def_dim(ostate%ncid, "z_interface", state%nh%nz + 1, zi_dimid)

         ! Sigma layer thickness coordinate
         call nc_check(nf90_def_var(ostate%ncid, "dsig", NC_WP, [z_dimid], dsig_varid), &
                       "defining dsig coordinate")
         call nc_put_att(ostate%ncid, dsig_varid, "units", "1")
         call nc_put_att(ostate%ncid, dsig_varid, "long_name", &
                         "sigma layer thickness (fraction of water column)")

         dimids_4d = [x_dimid, y_dimid, z_dimid, ostate%time_dimid]
         dimids_4d_w = [x_dimid, y_dimid, zi_dimid, ostate%time_dimid]

         call nc_def_var_4d(ostate%ncid, "u_layer", dimids_4d, ostate%u_layer_varid, dlevel)
         call nc_put_att(ostate%ncid, ostate%u_layer_varid, "units", "m s-1")
         call nc_put_att(ostate%ncid, ostate%u_layer_varid, "long_name", &
                         "layer x-velocity")
         call nc_put_att_real(ostate%ncid, ostate%u_layer_varid, "_FillValue", FILL_VALUE)

         call nc_def_var_4d(ostate%ncid, "v_layer", dimids_4d, ostate%v_layer_varid, dlevel)
         call nc_put_att(ostate%ncid, ostate%v_layer_varid, "units", "m s-1")
         call nc_put_att(ostate%ncid, ostate%v_layer_varid, "long_name", &
                         "layer y-velocity")
         call nc_put_att_real(ostate%ncid, ostate%v_layer_varid, "_FillValue", FILL_VALUE)

         call nc_def_var_4d(ostate%ncid, "w", dimids_4d_w, ostate%w_varid, dlevel)
         call nc_put_att(ostate%ncid, ostate%w_varid, "units", "m s-1")
         call nc_put_att(ostate%ncid, ostate%w_varid, "long_name", &
                         "vertical velocity at layer interfaces")
         call nc_put_att_real(ostate%ncid, ostate%w_varid, "_FillValue", FILL_VALUE)

         call nc_def_var_4d(ostate%ncid, "q_nh", dimids_4d, ostate%q_nh_varid, dlevel)
         call nc_put_att(ostate%ncid, ostate%q_nh_varid, "units", "m2 s-2")
         call nc_put_att(ostate%ncid, ostate%q_nh_varid, "long_name", &
                         "non-hydrostatic pressure")
         call nc_put_att_real(ostate%ncid, ostate%q_nh_varid, "_FillValue", FILL_VALUE)

         ! Vorticity diagnostics (3D: x, y, time)
         call nc_def_var_3d(ostate%ncid, "vorticity_x", dimids_3d, ostate%vort_x_varid, dlevel)
         call nc_put_att(ostate%ncid, ostate%vort_x_varid, "units", "s-1")
         call nc_put_att(ostate%ncid, ostate%vort_x_varid, "long_name", &
                         "x-component of vorticity (-dv/dz)")
         call nc_put_att_real(ostate%ncid, ostate%vort_x_varid, "_FillValue", FILL_VALUE)

         call nc_def_var_3d(ostate%ncid, "vorticity_y", dimids_3d, ostate%vort_y_varid, dlevel)
         call nc_put_att(ostate%ncid, ostate%vort_y_varid, "units", "s-1")
         call nc_put_att(ostate%ncid, ostate%vort_y_varid, "long_name", &
                         "y-component of vorticity (du/dz)")
         call nc_put_att_real(ostate%ncid, ostate%vort_y_varid, "_FillValue", FILL_VALUE)
      end if

      ! Multilayer variables (4D: x, y, z, time)
      if (state%use_multilayer) then
         ostate%has_ml = .true.
         ostate%nz_ml = state%multilayer%nz_ml

         call nc_def_dim(ostate%ncid, "z_ml", state%multilayer%nz_ml, z_dimid)
         call nc_def_dim(ostate%ncid, "z_ml_interface", state%multilayer%nz_ml + 1, zi_dimid)

         ! Sigma layer thickness coordinate
         call nc_check(nf90_def_var(ostate%ncid, "dsig_ml", NC_WP, [z_dimid], dsig_varid), &
                       "defining dsig_ml coordinate")
         call nc_put_att(ostate%ncid, dsig_varid, "units", "1")
         call nc_put_att(ostate%ncid, dsig_varid, "long_name", &
                         "sigma layer thickness (fraction of water column)")

         dimids_4d = [x_dimid, y_dimid, z_dimid, ostate%time_dimid]
         dimids_4d_w = [x_dimid, y_dimid, zi_dimid, ostate%time_dimid]

         if (cfg%out_h_layer) then
            call nc_def_var_4d(ostate%ncid, "h_layer", dimids_4d, ostate%ml_h_layer_varid, dlevel)
            call nc_put_att(ostate%ncid, ostate%ml_h_layer_varid, "units", "m")
            call nc_put_att(ostate%ncid, ostate%ml_h_layer_varid, "long_name", "layer thickness")
            call nc_put_att_real(ostate%ncid, ostate%ml_h_layer_varid, "_FillValue", FILL_VALUE)
         end if

         if (cfg%out_velocity) then
            call nc_def_var_4d(ostate%ncid, "u_ml", dimids_4d, ostate%ml_u_varid, dlevel)
            call nc_put_att(ostate%ncid, ostate%ml_u_varid, "units", "m s-1")
            call nc_put_att(ostate%ncid, ostate%ml_u_varid, "long_name", "layer x-velocity")
            call nc_put_att_real(ostate%ncid, ostate%ml_u_varid, "_FillValue", FILL_VALUE)

            call nc_def_var_4d(ostate%ncid, "v_ml", dimids_4d, ostate%ml_v_varid, dlevel)
            call nc_put_att(ostate%ncid, ostate%ml_v_varid, "units", "m s-1")
            call nc_put_att(ostate%ncid, ostate%ml_v_varid, "long_name", "layer y-velocity")
            call nc_put_att_real(ostate%ncid, ostate%ml_v_varid, "_FillValue", FILL_VALUE)
         end if

         ! Tracer registry: define one NetCDF variable per registered
         ! tracer, using metadata from the tracer entry itself.  Adding
         ! a sediment / passive scalar to `state%multilayer%tracers(:)`
         ! produces a correctly-labelled output variable with no
         ! changes to this module.  `cfg%out_salinity` /
         ! `cfg%out_temperature` are still honoured (legacy switches);
         ! future tracers default to write-enabled.
         if (allocated(state%multilayer%tracers)) then
            call define_tracer_outputs(ostate, state, cfg, dimids_4d, dlevel)
         end if

         ! Bed sediment reservoir — a time-varying 2D (x, y, time) field, only
         ! when sediment is active.  Lets deposition / erosion be seen; the
         ! suspended load already auto-outputs via the tracer registry above.
         if (state%multilayer%use_sediment) then
            call nc_def_var_3d(ostate%ncid, "bed_sediment", dimids_3d, &
                               ostate%bed_sediment_varid, dlevel)
            call nc_put_att(ostate%ncid, ostate%bed_sediment_varid, "units", "kg m-2")
            call nc_put_att(ostate%ncid, ostate%bed_sediment_varid, "long_name", &
                            "bed sediment mass per unit area")
            call nc_put_att_real(ostate%ncid, ostate%bed_sediment_varid, "_FillValue", FILL_VALUE)
            call nc_put_att(ostate%ncid, ostate%bed_sediment_varid, "coordinates", "x y")
         end if

         if (cfg%out_density) then
            call nc_def_var_4d(ostate%ncid, "rho", dimids_4d, ostate%ml_rho_varid, dlevel)
            call nc_put_att(ostate%ncid, ostate%ml_rho_varid, "units", "kg m-3")
            call nc_put_att(ostate%ncid, ostate%ml_rho_varid, "long_name", "layer density")
            call nc_put_att_real(ostate%ncid, ostate%ml_rho_varid, "_FillValue", FILL_VALUE)
         end if

         if (cfg%out_w) then
            call nc_def_var_4d(ostate%ncid, "w_ml", dimids_4d_w, ostate%ml_w_varid, dlevel)
            call nc_put_att(ostate%ncid, ostate%ml_w_varid, "units", "m s-1")
            call nc_put_att(ostate%ncid, ostate%ml_w_varid, "long_name", &
                            "diasurface (sigma-relative) transport velocity at layer interfaces")
            call nc_put_att_real(ostate%ncid, ostate%ml_w_varid, "_FillValue", FILL_VALUE)
         end if
      end if

      ! Global attributes
      call nc_put_att_global(ostate%ncid, "Conventions", "CF-1.11")
      call nc_put_att_global(ostate%ncid, "title", "Tidal shallow water solver output")
      call nc_put_att_global(ostate%ncid, "source", "Rakali v0.1.0")

      ! History — CF recommends an audit trail with timestamps
      call date_and_time(values=dt_values)
      write (history_str, '(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,&
            &"Z: Created by Rakali v0.1.0")') &
            dt_values(1), dt_values(2), dt_values(3), &
            dt_values(5), dt_values(6), dt_values(7)
      call nc_put_att_global(ostate%ncid, "history", trim(history_str))

      ! Optional institution
      if (len_trim(cfg%institution) > 0) then
         call nc_put_att_global(ostate%ncid, "institution", trim(cfg%institution))
      end if

      call nc_enddef(ostate%ncid)

      ! Write coordinate values
      ng = state%grid%nghost
      allocate (x_coords(state%grid%nx_phys), y_coords(state%grid%ny_phys))
      do i = 1, state%grid%nx_phys
         x_coords(i) = (real(i, wp) - 0.5_wp)*state%grid%dx
      end do
      do j = 1, state%grid%ny_phys
         y_coords(j) = (real(j, wp) - 0.5_wp)*state%grid%dy
      end do
      call nc_check(nf90_put_var(ostate%ncid, x_varid, x_coords), &
                    "writing x coordinates")
      call nc_check(nf90_put_var(ostate%ncid, y_varid, y_coords), &
                    "writing y coordinates")
      deallocate (x_coords, y_coords)

      ! Write sigma layer thicknesses (static)
      if (ostate%has_nh) then
         call nc_check(nf90_put_var(ostate%ncid, dsig_varid, state%nh%dsig), &
                       "writing dsig coordinate")
      end if
      if (ostate%has_ml) then
         call nc_check(nf90_put_var(ostate%ncid, dsig_varid, state%multilayer%dsig_ml), &
                       "writing dsig_ml coordinate")
      end if

      ! Write bathymetry (static field)
      call nc_put_var_2d(ostate%ncid, b_varid, &
                         state%barotropic%b(ng + 1:ng + state%grid%nx_phys, ng + 1:ng + state%grid%ny_phys))

      ostate%time_index = 0

   end subroutine output_init

   subroutine define_tracer_outputs(ostate, state, cfg, dimids_4d, dlevel)
      !! Define one NetCDF variable per entry in
      !! `state%multilayer%tracers(:)`, using each tracer's identity
      !! strings (`name`, `long_name`, `units`, `standard_name`) so
      !! the registry alone determines what gets written.  Sets
      !! `ostate%ml_tracer_varids(it)` to `-1` when the tracer is
      !! gated off (the legacy `cfg%out_salinity` /
      !! `cfg%out_temperature` switches still apply for those two
      !! built-in tracers; new tracers default to write-enabled).
      type(output_state_t), intent(inout) :: ostate
      type(state_t), intent(in)    :: state
      type(config_t), intent(in)    :: cfg
      integer, intent(in)    :: dimids_4d(4)
      integer, intent(in)    :: dlevel

      integer :: it, varid
      logical :: do_write
      character(len=32) :: varname

      allocate (ostate%ml_tracer_varids(size(state%multilayer%tracers)))
      ostate%ml_tracer_varids = -1

      do it = 1, size(state%multilayer%tracers)
         associate (tr => state%multilayer%tracers(it))
            varname = trim(tr%name)
            do_write = .true.
            ! Honour the legacy out_<name> switches for the built-in
            ! tracers; new tracers fall through to "always write".
            if (varname == "salinity") do_write = cfg%out_salinity
            if (varname == "temperature") do_write = cfg%out_temperature
            if (.not. do_write) cycle
            if (len_trim(varname) == 0) cycle  ! skip un-named tracers

            call nc_def_var_4d(ostate%ncid, trim(varname), dimids_4d, varid, dlevel)
            ostate%ml_tracer_varids(it) = varid

            if (len_trim(tr%units) > 0) then
               call nc_put_att(ostate%ncid, varid, "units", trim(tr%units))
            end if
            if (len_trim(tr%long_name) > 0) then
               call nc_put_att(ostate%ncid, varid, "long_name", trim(tr%long_name))
            end if
            if (len_trim(tr%standard_name) > 0) then
               call nc_put_att(ostate%ncid, varid, "standard_name", &
                               trim(tr%standard_name))
            end if
            call nc_put_att_real(ostate%ncid, varid, "_FillValue", FILL_VALUE)
         end associate
      end do
   end subroutine define_tracer_outputs

   subroutine write_tracer_snapshots(ostate, state, nx, ny, ng, nz_ml, buf3d)
      !! Per-tracer host-side write loop.  Divides the conserved
      !! `tracers(it)%hTr` (= h * Tr) by `h_layer` to recover the
      !! primitive variable, then writes through
      !! `ostate%ml_tracer_varids(it)` for each tracer that has output
      !! enabled.
      type(output_state_t), intent(inout) :: ostate
      type(state_t), intent(in)    :: state
      integer, intent(in)    :: nx, ny, ng, nz_ml
      real(wp), intent(inout) :: buf3d(nx, ny, nz_ml)

      integer :: it, ii, jj, kk
      real(wp) :: h_k

      do it = 1, size(state%multilayer%tracers)
         if (ostate%ml_tracer_varids(it) < 0) cycle

         do kk = 1, nz_ml
            do jj = 1, ny
               do ii = 1, nx
                  h_k = state%multilayer%h_layer(ng + ii, ng + jj, kk)
                  if (h_k > DRY_TOLERANCE) then
                     buf3d(ii, jj, kk) = &
                        state%multilayer%tracers(it)%hTr(ng + ii, ng + jj, kk)/h_k
                  else
                     buf3d(ii, jj, kk) = 0.0_wp
                  end if
               end do
            end do
         end do
         call nc_put_var_4d_slice(ostate%ncid, ostate%ml_tracer_varids(it), &
                                  buf3d, ostate%time_index)
      end do
   end subroutine write_tracer_snapshots

   subroutine output_write_snapshot(ostate, state, time)
      !! Write a solution snapshot at the current time
      type(output_state_t), intent(inout) :: ostate
         !! Output file state
      type(state_t), intent(in) :: state
         !! Current simulation state
      real(wp), intent(in) :: time
         !! Current simulation time (s)

      integer :: ng, nx, ny, i, j, nz
      real(wp) :: dz_inv
      real(wp), allocatable :: vort_x(:, :), vort_y(:, :)

      if (ostate%ncid < 0) return

      ostate%time_index = ostate%time_index + 1
      ng = state%grid%nghost
      nx = state%grid%nx_phys
      ny = state%grid%ny_phys

      ! Write time
      call nc_check(nf90_put_var(ostate%ncid, ostate%time_varid, &
                                 time, start=[ostate%time_index]), &
                    "writing time value")

      ! Write interior solution fields
      call nc_put_var_3d_slice(ostate%ncid, ostate%h_varid, &
                               state%barotropic%h(ng + 1:ng + nx, ng + 1:ng + ny), &
                               ostate%time_index)

      call nc_put_var_3d_slice(ostate%ncid, ostate%hu_varid, &
                               state%barotropic%hu(ng + 1:ng + nx, ng + 1:ng + ny), &
                               ostate%time_index)

      call nc_put_var_3d_slice(ostate%ncid, ostate%hv_varid, &
                               state%barotropic%hv(ng + 1:ng + nx, ng + 1:ng + ny), &
                               ostate%time_index)

      ! Non-hydrostatic fields
      if (ostate%has_nh) then
         nz = ostate%nz

         call nc_put_var_4d_slice(ostate%ncid, ostate%u_layer_varid, &
                                  state%nh%u_layer(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz), &
                                  ostate%time_index)

         call nc_put_var_4d_slice(ostate%ncid, ostate%v_layer_varid, &
                                  state%nh%v_layer(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz), &
                                  ostate%time_index)

         call nc_put_var_4d_slice(ostate%ncid, ostate%w_varid, &
                                  state%nh%w(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz + 1), &
                                  ostate%time_index)

         call nc_put_var_4d_slice(ostate%ncid, ostate%q_nh_varid, &
                                  state%nh%q_nh(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz), &
                                  ostate%time_index)

         ! Compute vorticity from vertical shear (host-side, nz >= 2)
         allocate (vort_x(nx, ny), vort_y(nx, ny))
         vort_x = 0.0_wp
         vort_y = 0.0_wp

         if (nz >= 2) then
            do j = 1, ny
               do i = 1, nx
                  dz_inv = 1.0_wp/(0.5_wp*state%barotropic%h(ng + i, ng + j) &
                                   *(ostate%dsig(1) + ostate%dsig(2)))
                  vort_y(i, j) = (state%nh%u_layer(ng + i, ng + j, 2) &
                                  - state%nh%u_layer(ng + i, ng + j, 1))*dz_inv
                  vort_x(i, j) = -(state%nh%v_layer(ng + i, ng + j, 2) &
                                   - state%nh%v_layer(ng + i, ng + j, 1))*dz_inv
               end do
            end do
         end if

         call nc_put_var_3d_slice(ostate%ncid, ostate%vort_x_varid, &
                                  vort_x, ostate%time_index)
         call nc_put_var_3d_slice(ostate%ncid, ostate%vort_y_varid, &
                                  vort_y, ostate%time_index)

         deallocate (vort_x, vort_y)
      end if

      ! Multilayer fields
      if (ostate%has_ml) then
         block
            integer :: nz_ml, ii, jj, kk, it
            real(wp), allocatable :: buf3d(:, :, :)
            real(wp) :: h_k
            logical :: need_vel, need_tracers

            nz_ml = ostate%nz_ml
            need_vel = (ostate%ml_u_varid >= 0)
            need_tracers = .false.
            if (allocated(ostate%ml_tracer_varids)) then
               need_tracers = any(ostate%ml_tracer_varids >= 0)
            end if

            if (ostate%ml_h_layer_varid >= 0) then
               call nc_put_var_4d_slice(ostate%ncid, ostate%ml_h_layer_varid, &
                                        state%multilayer%h_layer(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz_ml), &
                                        ostate%time_index)
            end if

            ! Bed sediment reservoir (2D slice at this time index).
            if (ostate%bed_sediment_varid >= 0) then
               call nc_put_var_3d_slice(ostate%ncid, ostate%bed_sediment_varid, &
                                        state%multilayer%bed_sediment(ng + 1:ng + nx, ng + 1:ng + ny), &
                                        ostate%time_index)
            end if

            ! Compute and write primitive variables that need h division
            if (need_vel .or. need_tracers) then
               allocate (buf3d(nx, ny, nz_ml))

               if (need_vel) then
                  ! u_ml
                  do kk = 1, nz_ml
                     do jj = 1, ny
                        do ii = 1, nx
                           h_k = state%multilayer%h_layer(ng + ii, ng + jj, kk)
                           if (h_k > DRY_TOLERANCE) then
                              buf3d(ii, jj, kk) = state%multilayer%hu_layer(ng + ii, ng + jj, kk)/h_k
                           else
                              buf3d(ii, jj, kk) = 0.0_wp
                           end if
                        end do
                     end do
                  end do
                  call nc_put_var_4d_slice(ostate%ncid, ostate%ml_u_varid, &
                                           buf3d, ostate%time_index)

                  ! v_ml
                  do kk = 1, nz_ml
                     do jj = 1, ny
                        do ii = 1, nx
                           h_k = state%multilayer%h_layer(ng + ii, ng + jj, kk)
                           if (h_k > DRY_TOLERANCE) then
                              buf3d(ii, jj, kk) = state%multilayer%hv_layer(ng + ii, ng + jj, kk)/h_k
                           else
                              buf3d(ii, jj, kk) = 0.0_wp
                           end if
                        end do
                     end do
                  end do
                  call nc_put_var_4d_slice(ostate%ncid, ostate%ml_v_varid, &
                                           buf3d, ostate%time_index)
               end if

               ! Per-tracer write loop — divides hTr by h_layer to
               ! recover the primitive variable (PSU, degC, kg/m^3,
               ! ...) using the same dry-cell guard as velocity.
               if (need_tracers) then
                  call write_tracer_snapshots(ostate, state, &
                                              nx, ny, ng, nz_ml, buf3d)
               end if

               deallocate (buf3d)
            end if

            if (ostate%ml_rho_varid >= 0) then
               call nc_put_var_4d_slice(ostate%ncid, ostate%ml_rho_varid, &
                                        state%multilayer%rho_layer(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz_ml), &
                                        ostate%time_index)
            end if

            if (ostate%ml_w_varid >= 0) then
               call nc_put_var_4d_slice(ostate%ncid, ostate%ml_w_varid, &
                                        state%multilayer%w_ml(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz_ml + 1), &
                                        ostate%time_index)
            end if
         end block
      end if

      call nc_check(nf90_sync(ostate%ncid), "syncing output file")

      call logger%debug("Output snapshot "//to_string(ostate%time_index)// &
                        " at t = "//to_string(time)//" s")

   end subroutine output_write_snapshot

   subroutine output_write_snapshot_from_arrays(ostate, h, hu, hv, b, time)
      !! Write a solution snapshot from pre-copied host arrays
      !!
      !! Used by the async output path: arrays are already in host memory
      !! (interior cells only, no ghost cells) so we write directly.
      !! NH fields are not supported in this path (hydrostatic only).
      type(output_state_t), intent(inout) :: ostate
         !! Output file state
      real(wp), intent(in) :: h(:, :)
         !! Water depth (nx_phys, ny_phys) — interior only
      real(wp), intent(in) :: hu(:, :)
         !! x-momentum (nx_phys, ny_phys) — interior only
      real(wp), intent(in) :: hv(:, :)
         !! y-momentum (nx_phys, ny_phys) — interior only
      real(wp), intent(in) :: b(:, :)
         !! Bathymetry (nx_phys, ny_phys) — interior only
      real(wp), intent(in) :: time
         !! Current simulation time (s)

      integer :: nx, ny

      if (ostate%ncid < 0) return

      ostate%time_index = ostate%time_index + 1
      nx = size(h, 1)
      ny = size(h, 2)

      ! Write time
      call nc_check(nf90_put_var(ostate%ncid, ostate%time_varid, &
                                 time, start=[ostate%time_index]), &
                    "writing time value")

      ! Write solution fields (already interior-only)
      call nc_put_var_3d_slice(ostate%ncid, ostate%h_varid, h, ostate%time_index)
      call nc_put_var_3d_slice(ostate%ncid, ostate%hu_varid, hu, ostate%time_index)
      call nc_put_var_3d_slice(ostate%ncid, ostate%hv_varid, hv, ostate%time_index)

      call nc_check(nf90_sync(ostate%ncid), "syncing output file")

      call logger%debug("Output snapshot "//to_string(ostate%time_index)// &
                        " at t = "//to_string(time)//" s")

   end subroutine output_write_snapshot_from_arrays

   subroutine output_close(ostate)
      !! Close the output NetCDF file
      type(output_state_t), intent(inout) :: ostate

      if (ostate%ncid >= 0) then
         call nc_close(ostate%ncid)
         ostate%ncid = -1
         if (allocated(ostate%dsig)) deallocate (ostate%dsig)
         if (allocated(ostate%ml_tracer_varids)) deallocate (ostate%ml_tracer_varids)
         call logger%info("Output file closed ("// &
                          to_string(ostate%time_index)//" snapshots written)")
      end if

   end subroutine output_close

   ! -----------------------------------------------------------------------
   ! Output filename generation
   ! -----------------------------------------------------------------------

   function output_rank_filename(base_dir, prefix, rank) result(fname)
      !! Generate a filename: <base_dir>/<prefix>_rank_NNNNNN.nc (the rank
      !! suffix is a leftover of the removed MPI support; callers pass 0)
      character(len=*), intent(in) :: base_dir
         !! Output directory
      character(len=*), intent(in) :: prefix
         !! File prefix (e.g. "output" or "restart")
      integer, intent(in) :: rank
         !! Rank suffix (always 0)
      character(len=512) :: fname

      write (fname, '(A,"/",A,"_rank_",I6.6,".nc")') &
         trim(base_dir), trim(prefix), rank

   end function output_rank_filename

   ! -----------------------------------------------------------------------
   ! Restart file I/O
   ! -----------------------------------------------------------------------

   subroutine restart_write(filename, state, time)
      !! Write a restart file containing full state for later continuation
      !!
      !! Saves the barotropic state (h, hu, hv, b) plus the multilayer
      !! arrays (h_layer, hu_layer, hv_layer, hS_layer) when
      !! state%use_multilayer is true.  The presence of the h_layer
      !! variable in the file is the marker for multilayer mode — old
      !! files written before this change simply lack it and round-trip
      !! as barotropic.
      !!
      !! NOTE: NH fields (u_layer, v_layer, w, q_nh) are NOT saved.
      !! A warm start from restart reverts to hydrostatic initial conditions
      !! for the NH solver state.
      character(len=*), intent(in) :: filename
         !! Path to restart file
      type(state_t), intent(in) :: state
         !! Current simulation state
      real(wp), intent(in) :: time
         !! Current simulation time (s)

      integer :: ncid, x_dimid, y_dimid, z_dimid
      integer :: h_varid, hu_varid, hv_varid, b_varid, time_varid
      integer :: h_layer_varid, hu_layer_varid, hv_layer_varid, hS_layer_varid, hT_layer_varid
      integer :: dimids_2d(2), dimids_3d(3)
      integer :: ng, nx, ny, nz

      call logger%info("Writing restart file: "//trim(filename))

      ng = state%grid%nghost
      nx = state%grid%nx_phys
      ny = state%grid%ny_phys

      call nc_create_file(filename, ncid)

      ! Dimensions
      call nc_def_dim(ncid, "x", nx, x_dimid)
      call nc_def_dim(ncid, "y", ny, y_dimid)
      dimids_2d = [x_dimid, y_dimid]

      ! Variables
      call nc_def_var_2d(ncid, "h", dimids_2d, h_varid)
      call nc_def_var_2d(ncid, "hu", dimids_2d, hu_varid)
      call nc_def_var_2d(ncid, "hv", dimids_2d, hv_varid)
      call nc_def_var_2d(ncid, "b", dimids_2d, b_varid)

      ! Multilayer arrays (only if active)
      if (state%use_multilayer) then
         nz = state%multilayer%nz_ml
         call nc_def_dim(ncid, "z", nz, z_dimid)
         dimids_3d = [x_dimid, y_dimid, z_dimid]
         call nc_def_var_3d(ncid, "h_layer", dimids_3d, h_layer_varid)
         call nc_def_var_3d(ncid, "hu_layer", dimids_3d, hu_layer_varid)
         call nc_def_var_3d(ncid, "hv_layer", dimids_3d, hv_layer_varid)
         call nc_def_var_3d(ncid, "hS_layer", dimids_3d, hS_layer_varid)
         call nc_def_var_3d(ncid, "hT_layer", dimids_3d, hT_layer_varid)
      end if

      ! Scalar time
      call nc_check(nf90_def_var(ncid, "time", nf90_double, time_varid), &
                    "defining restart time")

      ! Metadata attributes
      call nc_put_att(ncid, h_varid, "units", "m")
      call nc_put_att(ncid, hu_varid, "units", "m2 s-1")
      call nc_put_att(ncid, hv_varid, "units", "m2 s-1")
      call nc_put_att(ncid, b_varid, "units", "m")
      call nc_put_att(ncid, time_varid, "units", "seconds since simulation start")
      if (state%use_multilayer) then
         call nc_put_att(ncid, h_layer_varid, "units", "m")
         call nc_put_att(ncid, hu_layer_varid, "units", "m2 s-1")
         call nc_put_att(ncid, hv_layer_varid, "units", "m2 s-1")
         call nc_put_att(ncid, hS_layer_varid, "units", "m PSU")
         call nc_put_att(ncid, hT_layer_varid, "units", "m degC")
      end if
      call nc_put_att_global(ncid, "Conventions", "CF-1.11")
      call nc_put_att_global(ncid, "title", "Tidal restart file")
      call nc_put_att_global(ncid, "source", "Rakali v0.1.0")
      call nc_put_att_global(ncid, "dx", to_string(state%grid%dx))
      call nc_put_att_global(ncid, "dy", to_string(state%grid%dy))
      call nc_put_att_global(ncid, "nghost", to_string(state%grid%nghost))

      call nc_enddef(ncid)

      ! Write data (interior only)
      call nc_put_var_2d(ncid, h_varid, state%barotropic%h(ng + 1:ng + nx, ng + 1:ng + ny))
      call nc_put_var_2d(ncid, hu_varid, state%barotropic%hu(ng + 1:ng + nx, ng + 1:ng + ny))
      call nc_put_var_2d(ncid, hv_varid, state%barotropic%hv(ng + 1:ng + nx, ng + 1:ng + ny))
      call nc_put_var_2d(ncid, b_varid, state%barotropic%b(ng + 1:ng + nx, ng + 1:ng + ny))

      if (state%use_multilayer) then
         call nc_check(nf90_put_var(ncid, h_layer_varid, &
                                    state%multilayer%h_layer(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz)), &
                       "writing restart h_layer")
         call nc_check(nf90_put_var(ncid, hu_layer_varid, &
                                    state%multilayer%hu_layer(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz)), &
                       "writing restart hu_layer")
         call nc_check(nf90_put_var(ncid, hv_layer_varid, &
                                    state%multilayer%hv_layer(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz)), &
                       "writing restart hv_layer")
         call nc_check(nf90_put_var(ncid, hS_layer_varid, &
                               state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz)), &
                       "writing restart hS_layer")
         call nc_check(nf90_put_var(ncid, hT_layer_varid, &
                            state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz)), &
                       "writing restart hT_layer")
      end if

      call nc_check(nf90_put_var(ncid, time_varid, time), "writing restart time")

      call nc_close(ncid)

      call logger%info("Restart file written at t = "//to_string(time)//" s")

   end subroutine restart_write

   subroutine restart_read(filename, state, time, ierr, loaded_multilayer)
      !! Read a restart file and restore the simulation state
      !!
      !! Grid dimensions in the restart file must match the current state.
      !! The time variable is returned so the caller can resume the time loop.
      !!
      !! Multilayer fields (h_layer, hu_layer, hv_layer, hS_layer) are
      !! restored when both the file and the destination state have
      !! multilayer mode active and the z dimension matches.  Files
      !! without these variables (older restart format) round-trip as
      !! barotropic; mismatched cases log a warning and skip the layer
      !! fields, leaving them at their post-allocation values.
      !!
      !! NOTE: NH fields (u_layer, v_layer, w, q_nh) are NOT restored.
      character(len=*), intent(in) :: filename
         !! Path to restart file
      type(state_t), intent(inout) :: state
         !! Simulation state (arrays must be pre-allocated)
      real(wp), intent(out) :: time
         !! Simulation time from the restart file (s)
      integer, intent(out), optional :: ierr
         !! 0 on success, non-zero on failure (currently: 1 = grid mismatch)
      logical, intent(out), optional :: loaded_multilayer
         !! .true. iff multilayer fields were actually populated from
         !! the file.  Callers that auto-initialise layers from the
         !! barotropic state should skip that step when this is true.

      integer :: ncid, varid, file_nx, file_ny, file_nz, status
      integer :: h_layer_varid
      integer :: ng, nx, ny, nz
      logical :: file_has_multilayer

      call logger%info("Reading restart file: "//trim(filename))

      if (present(ierr)) ierr = 0
      if (present(loaded_multilayer)) loaded_multilayer = .false.

      ng = state%grid%nghost
      nx = state%grid%nx_phys
      ny = state%grid%ny_phys

      call nc_open_read(filename, ncid)

      ! Validate dimensions
      call nc_get_dim_len(ncid, "x", file_nx)
      call nc_get_dim_len(ncid, "y", file_ny)

      if (file_nx /= nx .or. file_ny /= ny) then
         call logger%error("Restart grid mismatch: file has "// &
                           to_string(file_nx)//" x "//to_string(file_ny)// &
                           " but simulation expects "// &
                           to_string(nx)//" x "//to_string(ny))
         call nc_close(ncid)
         time = 0.0_wp
         if (present(ierr)) ierr = 1
         return
      end if

      ! Read state variables into interior
      call nc_get_varid(ncid, "h", varid)
      call nc_get_var_2d(ncid, varid, state%barotropic%h(ng + 1:ng + nx, ng + 1:ng + ny))

      call nc_get_varid(ncid, "hu", varid)
      call nc_get_var_2d(ncid, varid, state%barotropic%hu(ng + 1:ng + nx, ng + 1:ng + ny))

      call nc_get_varid(ncid, "hv", varid)
      call nc_get_var_2d(ncid, varid, state%barotropic%hv(ng + 1:ng + nx, ng + 1:ng + ny))

      call nc_get_varid(ncid, "b", varid)
      call nc_get_var_2d(ncid, varid, state%barotropic%b(ng + 1:ng + nx, ng + 1:ng + ny))

      ! Multilayer fields: presence of h_layer is the marker
      status = nf90_inq_varid(ncid, "h_layer", h_layer_varid)
      file_has_multilayer = (status == nf90_noerr)

      if (file_has_multilayer .and. state%use_multilayer) then
         nz = state%multilayer%nz_ml
         call nc_get_dim_len(ncid, "z", file_nz)
         if (file_nz /= nz) then
            call logger%warning("Restart multilayer mismatch: file has "// &
                                to_string(file_nz)//" layers but state expects "// &
                                to_string(nz)//"; skipping layer fields")
         else
            call nc_check(nf90_get_var(ncid, h_layer_varid, &
                                       state%multilayer%h_layer(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz)), &
                          "reading restart h_layer")
            call nc_get_varid(ncid, "hu_layer", varid)
            call nc_check(nf90_get_var(ncid, varid, &
                                       state%multilayer%hu_layer(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz)), &
                          "reading restart hu_layer")
            call nc_get_varid(ncid, "hv_layer", varid)
            call nc_check(nf90_get_var(ncid, varid, &
                                       state%multilayer%hv_layer(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz)), &
                          "reading restart hv_layer")
            call nc_get_varid(ncid, "hS_layer", varid)
            call nc_check(nf90_get_var(ncid, varid, &
                               state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz)), &
                          "reading restart hS_layer")
            ! hT_layer added in the temperature PR; older restart files
            ! lack this variable.  Fall back to initialising hT from
            ! initial_temperature so they remain usable.
            status = nf90_inq_varid(ncid, "hT_layer", varid)
            if (status == nf90_noerr) then
               call nc_check(nf90_get_var(ncid, varid, &
                            state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz)), &
                             "reading restart hT_layer")
            else
               call logger%warning("Restart file lacks hT_layer; initialising "// &
                                   "temperature to initial_temperature = "// &
                                   to_string(state%multilayer%initial_temperature)//" degC")
               state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz) = &
                  state%multilayer%h_layer(ng + 1:ng + nx, ng + 1:ng + ny, 1:nz)* &
                  state%multilayer%initial_temperature
            end if
            if (present(loaded_multilayer)) loaded_multilayer = .true.
         end if
      else if (file_has_multilayer .and. .not. state%use_multilayer) then
         call logger%warning("Restart file contains multilayer fields but "// &
                             "state has multilayer disabled; skipping layer fields")
      else if (.not. file_has_multilayer .and. state%use_multilayer) then
         call logger%warning("Restart file is barotropic but state has "// &
                             "multilayer enabled; layer fields will be re-initialised by "// &
                             "the caller")
      end if

      ! Read time
      call nc_get_varid(ncid, "time", varid)
      call nc_check(nf90_get_var(ncid, varid, time), "reading restart time")

      call nc_close(ncid)

      call logger%info("Restart loaded: t = "//to_string(time)//" s")

   end subroutine restart_read

   subroutine write_ascii_output(s, output_dir)
      !! Write the final water depth field to an ASCII file
      type(state_t), intent(in) :: s
      character(len=*), intent(in) :: output_dir

      integer :: i, j, io_unit, io_stat
      character(len=256) :: io_msg
      character(len=512) :: filename

      filename = trim(output_dir)//"/final_h.dat"

      open (newunit=io_unit, file=trim(filename), status="replace", &
            action="write", iostat=io_stat, iomsg=io_msg)

      if (io_stat /= 0) then
         call logger%warning("Could not open output file "// &
                             trim(filename)//": "//trim(io_msg))
         ! Fall back to writing to current directory
         filename = "final_h.dat"
         open (newunit=io_unit, file=trim(filename), status="replace", &
               action="write", iostat=io_stat, iomsg=io_msg)
         if (io_stat /= 0) then
            call logger%error("Could not write output file: "//trim(io_msg))
            return
         end if
      end if

      ! Write header
      write (io_unit, "(A)") "# x  h"

      ! Write 1D profile along centreline (j = middle)
      j = s%grid%nghost + s%grid%ny_phys/2 + 1
      do i = s%grid%nghost + 1, s%grid%nx_total - s%grid%nghost
         write (io_unit, "(2ES20.12)") &
            (real(i - s%grid%nghost, wp) - 0.5_wp)*s%grid%dx, s%barotropic%h(i, j)
      end do

      close (io_unit)
      call logger%info("Output written to: "//trim(filename))

   end subroutine write_ascii_output

end module rki_output
