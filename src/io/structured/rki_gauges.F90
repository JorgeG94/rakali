!! Gauge (station) point output — high-frequency time series at specific locations
!!
!! Reads gauge locations from a text file, computes bilinear interpolation
!! weights at initialisation, and writes interpolated eta/u/v to a NetCDF
!! time series file at a configurable interval.
module rki_gauges
   use rki_constants, only: wp
   use rki_config, only: config_t
   use rki_state, only: state_t
   use rki_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, &
                            nc_put_att, nc_put_att_real, nc_put_att_global, &
                            nc_check, NC_WP
   use netcdf, only: nf90_def_var, nf90_enddef, nf90_put_var, nf90_sync, &
                     nf90_double, nf90_char, nf90_unlimited, nf90_global, &
                     nf90_noerr
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: gauge_output_t
   public :: gauges_init
   public :: gauges_write
   public :: gauges_destroy

   integer, parameter :: MAX_NAME_LEN = 64
      !! Maximum length of a gauge name string
   integer, parameter :: MAX_GAUGES = 100
      !! Maximum number of gauges that can be read from file
   integer, parameter :: MAX_LINE_LEN = 512
      !! Maximum length of a line in the gauge file

   type :: gauge_output_t
      !! Manages gauge point output to NetCDF time series

      logical :: active = .false.
         !! Whether gauge output is enabled
      integer :: n_gauges = 0
         !! Number of gauge points
      character(len=MAX_NAME_LEN), allocatable :: names(:)
         !! Gauge names
      real(wp), allocatable :: gauge_x(:)
         !! Gauge x-coordinates (metres, in global domain)
      real(wp), allocatable :: gauge_y(:)
         !! Gauge y-coordinates (metres, in global domain)

      ! Bilinear interpolation (cell indices and weights)
      integer, allocatable :: ix0(:)
         !! Lower-left cell i-index for each local gauge
      integer, allocatable :: iy0(:)
         !! Lower-left cell j-index for each local gauge
      real(wp), allocatable :: wx(:)
         !! x-interpolation weight (0-1)
      real(wp), allocatable :: wy(:)
         !! y-interpolation weight (0-1)

      ! Gauge subset (gauges that fall inside the domain)
      integer :: n_local = 0
         !! Number of gauges inside the domain
      integer, allocatable :: local_idx(:)
         !! Maps local gauge index to global gauge index

      ! NetCDF output state
      integer :: ncid = -1
         !! NetCDF file ID
      integer :: time_varid = -1
         !! Time variable ID
      integer :: eta_varid = -1
         !! Water surface elevation variable ID
      integer :: u_varid = -1
         !! Depth-averaged u-velocity variable ID
      integer :: v_varid = -1
         !! Depth-averaged v-velocity variable ID
      integer :: time_index = 0
         !! Current time record index
   end type gauge_output_t

contains

   subroutine gauges_init(gauges, cfg, state, rank)
      !! Initialise gauge output: read locations, compute weights, create file
      type(gauge_output_t), intent(inout) :: gauges
         !! Gauge output state
      type(config_t), intent(in) :: cfg
         !! Simulation configuration
      type(state_t), intent(in) :: state
         !! Solution state (for grid metadata and bathymetry)
      integer, intent(in) :: rank
         !! File-name rank suffix (always 0; the suffix is a leftover of the removed MPI support)

      gauges%active = .false.

      ! If no gauge file specified, nothing to do
      if (len_trim(cfg%gauge_file) == 0) return

      ! Read gauge locations from text file
      call read_gauge_file(cfg%gauge_file, gauges%names, gauges%gauge_x, &
                           gauges%gauge_y, gauges%n_gauges)

      if (gauges%n_gauges == 0) then
         call logger%info("Gauges: no gauge locations found in file")
         return
      end if

      ! Find which gauges fall inside the domain
      call compute_local_gauges(gauges, state)

      if (gauges%n_local == 0) then
         call logger%info("Gauges: no gauge locations inside the domain")
         return
      end if

      ! Compute bilinear interpolation weights for local gauges
      call compute_interp_weights(gauges, state)

      ! Create the NetCDF output file
      call create_gauge_netcdf(gauges, state, cfg, rank)

      gauges%active = .true.
      call logger%info("Gauges: "//trim(to_string(gauges%n_local))// &
                       " of "//trim(to_string(gauges%n_gauges))// &
                       " gauges active on this rank")
   end subroutine gauges_init

   ! ---------------------------------------------------------------------------
   ! Private implementation routines
   ! ---------------------------------------------------------------------------

   subroutine read_gauge_file(filename, names, gauge_x, gauge_y, n_gauges)
      !! Parse a gauge locations text file
      !!
      !! Format: one gauge per line, ``name  x  y``.  Lines starting with
      !! ``#`` or blank lines are skipped.
      character(len=*), intent(in) :: filename
         !! Path to gauge locations file
      character(len=MAX_NAME_LEN), allocatable, intent(out) :: names(:)
         !! Gauge names
      real(wp), allocatable, intent(out) :: gauge_x(:)
         !! Gauge x-coordinates (m)
      real(wp), allocatable, intent(out) :: gauge_y(:)
         !! Gauge y-coordinates (m)
      integer, intent(out) :: n_gauges
         !! Number of gauges read

      character(len=MAX_LINE_LEN) :: line
      character(len=MAX_NAME_LEN) :: tmp_names(MAX_GAUGES)
      real(wp) :: tmp_x(MAX_GAUGES), tmp_y(MAX_GAUGES)
      integer :: iunit, ios, n
      character(len=256) :: ios_msg
      logical :: file_exists

      n_gauges = 0

      inquire (file=trim(filename), exist=file_exists)
      if (.not. file_exists) then
         call logger%error("Gauges: file not found: "//trim(filename))
         return
      end if

      open (newunit=iunit, file=trim(filename), status="old", action="read", &
            iostat=ios, iomsg=ios_msg)
      if (ios /= 0) then
         call logger%error("Gauges: cannot open file "//trim(filename)// &
                           ": "//trim(ios_msg))
         return
      end if

      n = 0
      do
         read (iunit, "(A)", iostat=ios, iomsg=ios_msg) line
         if (ios /= 0) exit

         ! Skip comments and blank lines
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == "#") cycle

         n = n + 1
         if (n > MAX_GAUGES) then
            call logger%warning("Gauges: maximum of "// &
                                trim(to_string(MAX_GAUGES))// &
                                " gauges exceeded, ignoring remainder")
            n = MAX_GAUGES
            exit
         end if

         read (line, *, iostat=ios, iomsg=ios_msg) tmp_names(n), tmp_x(n), tmp_y(n)
         if (ios /= 0) then
            call logger%warning("Gauges: skipping malformed line: "// &
                                trim(line)//" ("//trim(ios_msg)//")")
            n = n - 1
            cycle
         end if
      end do

      close (iunit)

      n_gauges = n
      allocate (names(n), gauge_x(n), gauge_y(n))
      names(1:n) = tmp_names(1:n)
      gauge_x(1:n) = tmp_x(1:n)
      gauge_y(1:n) = tmp_y(1:n)

      call logger%info("Gauges: read "//trim(to_string(n_gauges))// &
                       " gauge locations from "//trim(filename))
   end subroutine read_gauge_file

   subroutine compute_local_gauges(gauges, state)
      !! Determine which gauges fall inside the domain
      type(gauge_output_t), intent(inout) :: gauges
         !! Gauge output state
      type(state_t), intent(in) :: state
         !! Solution state (for grid metadata)

      integer :: g, n, ig, jg
      integer :: i_start, j_start, nx_phys, ny_phys
      integer :: tmp_local(MAX_GAUGES)

      i_start = 1
      j_start = 1
      nx_phys = state%grid%nx_phys
      ny_phys = state%grid%ny_phys

      n = 0
      do g = 1, gauges%n_gauges
         ! Global physical cell containing this gauge (1-based)
         ig = floor(gauges%gauge_x(g)/state%grid%dx) + 1
         jg = floor(gauges%gauge_y(g)/state%grid%dy) + 1

         ! Check if inside the domain
         if (ig >= i_start .and. ig <= i_start + nx_phys - 1 .and. &
             jg >= j_start .and. jg <= j_start + ny_phys - 1) then
            n = n + 1
            tmp_local(n) = g
         end if
      end do

      gauges%n_local = n
      if (n > 0) then
         allocate (gauges%local_idx(n))
         gauges%local_idx(1:n) = tmp_local(1:n)
      end if
   end subroutine compute_local_gauges

   subroutine compute_interp_weights(gauges, state)
      !! Compute bilinear interpolation indices and weights for local gauges
      type(gauge_output_t), intent(inout) :: gauges
         !! Gauge output state
      type(state_t), intent(in) :: state
         !! Solution state

      integer :: g, gi, ig, jg, il, jl, ia, ja
      real(wp) :: x_centre, y_centre, gx, gy, dx, dy
      integer :: nghost, i_start, j_start, nx_phys, ny_phys

      dx = state%grid%dx
      dy = state%grid%dy
      nghost = state%grid%nghost
      i_start = 1
      j_start = 1
      nx_phys = state%grid%nx_phys
      ny_phys = state%grid%ny_phys

      allocate (gauges%ix0(gauges%n_local))
      allocate (gauges%iy0(gauges%n_local))
      allocate (gauges%wx(gauges%n_local))
      allocate (gauges%wy(gauges%n_local))

      do g = 1, gauges%n_local
         gi = gauges%local_idx(g)
         gx = gauges%gauge_x(gi)
         gy = gauges%gauge_y(gi)

         ! Global physical cell containing gauge (1-based)
         ig = floor(gx/dx) + 1
         jg = floor(gy/dy) + 1

         ! Local physical cell (1-based, 1..nx_phys)
         il = ig - i_start + 1
         jl = jg - j_start + 1

         ! Array index (with ghost offset)
         ia = il + nghost
         ja = jl + nghost

         ! Cell centre in global coordinates
         x_centre = (real(ig, wp) - 0.5_wp)*dx
         y_centre = (real(jg, wp) - 0.5_wp)*dy

         ! Determine interpolation stencil: if gauge is left of centre,
         ! interpolate between (ia-1, ia); if right, between (ia, ia+1)
         if (gx < x_centre) then
            ! Gauge is left of cell centre
            if (il > 1) then
               gauges%ix0(g) = ia - 1
               gauges%wx(g) = (gx - (x_centre - dx))/dx
            else
               ! At left edge of the domain, clamp
               gauges%ix0(g) = ia
               gauges%wx(g) = 0.0_wp
            end if
         else
            ! Gauge is right of (or at) cell centre
            if (il < nx_phys) then
               gauges%ix0(g) = ia
               gauges%wx(g) = (gx - x_centre)/dx
            else
               ! At right edge of the domain, clamp
               gauges%ix0(g) = ia
               gauges%wx(g) = 0.0_wp
            end if
         end if

         if (gy < y_centre) then
            if (jl > 1) then
               gauges%iy0(g) = ja - 1
               gauges%wy(g) = (gy - (y_centre - dy))/dy
            else
               gauges%iy0(g) = ja
               gauges%wy(g) = 0.0_wp
            end if
         else
            if (jl < ny_phys) then
               gauges%iy0(g) = ja
               gauges%wy(g) = (gy - y_centre)/dy
            else
               gauges%iy0(g) = ja
               gauges%wy(g) = 0.0_wp
            end if
         end if

         ! Clamp weights to [0, 1] for safety
         gauges%wx(g) = max(0.0_wp, min(1.0_wp, gauges%wx(g)))
         gauges%wy(g) = max(0.0_wp, min(1.0_wp, gauges%wy(g)))
      end do
   end subroutine compute_interp_weights

   subroutine create_gauge_netcdf(gauges, state, cfg, rank)
      !! Create the gauge output NetCDF file with dimensions and static data
      type(gauge_output_t), intent(inout) :: gauges
         !! Gauge output state
      type(state_t), intent(in) :: state
         !! Solution state (for bathymetry interpolation)
      type(config_t), intent(in) :: cfg
         !! Simulation configuration
      integer, intent(in) :: rank
         !! File-name rank suffix (always 0; the suffix is a leftover of the removed MPI support)

      character(len=512) :: filename
      integer :: ncid, station_dimid, time_dimid, strlen_dimid
      integer :: name_varid, sx_varid, sy_varid, sb_varid
      integer :: dimids_1d(1), dimids_2d(2), dimids_name(2)
      integer :: g, gi, i0, j0, status
      real(wp) :: wxg, wyg, val
      real(wp), allocatable :: station_x(:), station_y(:), station_b(:)
      character(len=MAX_NAME_LEN) :: padded_name

      ! Build filename
      if (rank > 0) then
         write (filename, "(A,A,I6.6,A)") trim(cfg%output_dir), &
            "/gauges_rank_", rank, ".nc"
      else
         filename = trim(cfg%output_dir)//"/gauges.nc"
      end if

      call nc_create_file(trim(filename), ncid)
      gauges%ncid = ncid

      ! Define dimensions
      call nc_def_dim(ncid, "station", gauges%n_local, station_dimid)
      call nc_def_dim(ncid, "time", nf90_unlimited, time_dimid)
      call nc_def_dim(ncid, "name_strlen", MAX_NAME_LEN, strlen_dimid)

      ! Define time variable
      dimids_1d = [time_dimid]
      status = nf90_def_var(ncid, "time", nf90_double, dimids_1d, &
                            gauges%time_varid)
      call nc_check(status, "def_var time")
      if (len_trim(cfg%reference_time) > 0) then
         call nc_put_att(ncid, gauges%time_varid, "units", &
                         "seconds since "//trim(cfg%reference_time))
      else
         call nc_put_att(ncid, gauges%time_varid, "units", "seconds")
      end if
      call nc_put_att(ncid, gauges%time_varid, "calendar", "proleptic_gregorian")

      ! Define station_name variable (name_strlen, station)
      dimids_name = [strlen_dimid, station_dimid]
      status = nf90_def_var(ncid, "station_name", nf90_char, dimids_name, &
                            name_varid)
      call nc_check(status, "def_var station_name")

      ! Define station coordinate variables
      dimids_1d = [station_dimid]
      status = nf90_def_var(ncid, "station_x", nf90_double, dimids_1d, &
                            sx_varid)
      call nc_check(status, "def_var station_x")
      call nc_put_att(ncid, sx_varid, "units", "m")
      call nc_put_att(ncid, sx_varid, "long_name", "gauge x-coordinate")

      status = nf90_def_var(ncid, "station_y", nf90_double, dimids_1d, &
                            sy_varid)
      call nc_check(status, "def_var station_y")
      call nc_put_att(ncid, sy_varid, "units", "m")
      call nc_put_att(ncid, sy_varid, "long_name", "gauge y-coordinate")

      status = nf90_def_var(ncid, "station_b", nf90_double, dimids_1d, &
                            sb_varid)
      call nc_check(status, "def_var station_b")
      call nc_put_att(ncid, sb_varid, "units", "m")
      call nc_put_att(ncid, sb_varid, "long_name", &
                      "interpolated bottom elevation")

      ! Define time-varying output variables (station, time)
      dimids_2d = [station_dimid, time_dimid]
      status = nf90_def_var(ncid, "eta", nf90_double, dimids_2d, &
                            gauges%eta_varid)
      call nc_check(status, "def_var eta")
      call nc_put_att(ncid, gauges%eta_varid, "units", "m")
      call nc_put_att(ncid, gauges%eta_varid, "long_name", &
                      "water surface elevation")

      status = nf90_def_var(ncid, "u", nf90_double, dimids_2d, &
                            gauges%u_varid)
      call nc_check(status, "def_var u")
      call nc_put_att(ncid, gauges%u_varid, "units", "m/s")
      call nc_put_att(ncid, gauges%u_varid, "long_name", &
                      "depth-averaged x-velocity")

      status = nf90_def_var(ncid, "v", nf90_double, dimids_2d, &
                            gauges%v_varid)
      call nc_check(status, "def_var v")
      call nc_put_att(ncid, gauges%v_varid, "units", "m/s")
      call nc_put_att(ncid, gauges%v_varid, "long_name", &
                      "depth-averaged y-velocity")

      ! Global attributes
      call nc_put_att_global(ncid, "title", "Rakali gauge time series")
      call nc_put_att_global(ncid, "Conventions", "CF-1.8")
      if (len_trim(cfg%institution) > 0) then
         call nc_put_att_global(ncid, "institution", trim(cfg%institution))
      end if

      ! End define mode
      status = nf90_enddef(ncid)
      call nc_check(status, "enddef gauges")

      ! Write static data: station names, coordinates, interpolated bathymetry
      allocate (station_x(gauges%n_local))
      allocate (station_y(gauges%n_local))
      allocate (station_b(gauges%n_local))

      do g = 1, gauges%n_local
         gi = gauges%local_idx(g)
         station_x(g) = gauges%gauge_x(gi)
         station_y(g) = gauges%gauge_y(gi)

         ! Write station name as fixed-length padded string
         padded_name = gauges%names(gi)
         status = nf90_put_var(ncid, name_varid, trim(padded_name), &
                               start=[1, g], count=[len_trim(padded_name), 1])
         call nc_check(status, "put_var station_name")

         ! Interpolate bathymetry at gauge location
         i0 = gauges%ix0(g)
         j0 = gauges%iy0(g)
         wxg = gauges%wx(g)
         wyg = gauges%wy(g)
         station_b(g) = (1.0_wp - wxg)*(1.0_wp - wyg)*state%barotropic%b(i0, j0) &
                        + wxg*(1.0_wp - wyg)*state%barotropic%b(i0 + 1, j0) &
                        + (1.0_wp - wxg)*wyg*state%barotropic%b(i0, j0 + 1) &
                        + wxg*wyg*state%barotropic%b(i0 + 1, j0 + 1)
      end do

      status = nf90_put_var(ncid, sx_varid, station_x)
      call nc_check(status, "put_var station_x")
      status = nf90_put_var(ncid, sy_varid, station_y)
      call nc_check(status, "put_var station_y")
      status = nf90_put_var(ncid, sb_varid, station_b)
      call nc_check(status, "put_var station_b")

      deallocate (station_x, station_y, station_b)

      status = nf90_sync(ncid)
      call nc_check(status, "sync gauges init")

      call logger%info("Gauges: created "//trim(filename))
   end subroutine create_gauge_netcdf

   subroutine gauges_write(gauges, state, time)
      !! Write a gauge time series record (interpolated eta, u, v)
      type(gauge_output_t), intent(inout) :: gauges
         !! Gauge output state
      type(state_t), intent(in) :: state
         !! Solution state
      real(wp), intent(in) :: time
         !! Current simulation time (s)

      integer :: g, i0, j0, status
      real(wp) :: wxg, wyg, h_interp, hu_interp, hv_interp, b_interp
      real(wp) :: h_safe
      real(wp), allocatable :: eta_buf(:), u_buf(:), v_buf(:)
      real(wp), parameter :: MIN_DEPTH = 1.0e-6_wp
         !! Minimum depth for velocity computation

      if (.not. gauges%active) return

      gauges%time_index = gauges%time_index + 1

      allocate (eta_buf(gauges%n_local))
      allocate (u_buf(gauges%n_local))
      allocate (v_buf(gauges%n_local))

      do g = 1, gauges%n_local
         i0 = gauges%ix0(g)
         j0 = gauges%iy0(g)
         wxg = gauges%wx(g)
         wyg = gauges%wy(g)

         ! Bilinear interpolation of h, hu, hv, b
         h_interp = (1.0_wp - wxg)*(1.0_wp - wyg)*state%barotropic%h(i0, j0) &
                    + wxg*(1.0_wp - wyg)*state%barotropic%h(i0 + 1, j0) &
                    + (1.0_wp - wxg)*wyg*state%barotropic%h(i0, j0 + 1) &
                    + wxg*wyg*state%barotropic%h(i0 + 1, j0 + 1)

         hu_interp = (1.0_wp - wxg)*(1.0_wp - wyg)*state%barotropic%hu(i0, j0) &
                     + wxg*(1.0_wp - wyg)*state%barotropic%hu(i0 + 1, j0) &
                     + (1.0_wp - wxg)*wyg*state%barotropic%hu(i0, j0 + 1) &
                     + wxg*wyg*state%barotropic%hu(i0 + 1, j0 + 1)

         hv_interp = (1.0_wp - wxg)*(1.0_wp - wyg)*state%barotropic%hv(i0, j0) &
                     + wxg*(1.0_wp - wyg)*state%barotropic%hv(i0 + 1, j0) &
                     + (1.0_wp - wxg)*wyg*state%barotropic%hv(i0, j0 + 1) &
                     + wxg*wyg*state%barotropic%hv(i0 + 1, j0 + 1)

         b_interp = (1.0_wp - wxg)*(1.0_wp - wyg)*state%barotropic%b(i0, j0) &
                    + wxg*(1.0_wp - wyg)*state%barotropic%b(i0 + 1, j0) &
                    + (1.0_wp - wxg)*wyg*state%barotropic%b(i0, j0 + 1) &
                    + wxg*wyg*state%barotropic%b(i0 + 1, j0 + 1)

         ! Water surface elevation: eta = h + b
         eta_buf(g) = h_interp + b_interp

         ! Depth-averaged velocities with thin-layer guard
         h_safe = max(h_interp, MIN_DEPTH)
         u_buf(g) = hu_interp/h_safe
         v_buf(g) = hv_interp/h_safe
      end do

      ! Write time
      status = nf90_put_var(gauges%ncid, gauges%time_varid, time, &
                            start=[gauges%time_index])
      call nc_check(status, "put_var gauge time")

      ! Write eta, u, v (station, time) — write all stations at this time index
      status = nf90_put_var(gauges%ncid, gauges%eta_varid, eta_buf, &
                            start=[1, gauges%time_index], &
                            count=[gauges%n_local, 1])
      call nc_check(status, "put_var gauge eta")

      status = nf90_put_var(gauges%ncid, gauges%u_varid, u_buf, &
                            start=[1, gauges%time_index], &
                            count=[gauges%n_local, 1])
      call nc_check(status, "put_var gauge u")

      status = nf90_put_var(gauges%ncid, gauges%v_varid, v_buf, &
                            start=[1, gauges%time_index], &
                            count=[gauges%n_local, 1])
      call nc_check(status, "put_var gauge v")

      deallocate (eta_buf, u_buf, v_buf)

      status = nf90_sync(gauges%ncid)
      call nc_check(status, "sync gauges write")
   end subroutine gauges_write

   subroutine gauges_destroy(gauges)
      !! Close the gauge NetCDF file and deallocate all arrays
      type(gauge_output_t), intent(inout) :: gauges
         !! Gauge output state

      if (gauges%ncid >= 0) then
         call nc_close(gauges%ncid)
         gauges%ncid = -1
      end if

      if (allocated(gauges%names)) deallocate (gauges%names)
      if (allocated(gauges%gauge_x)) deallocate (gauges%gauge_x)
      if (allocated(gauges%gauge_y)) deallocate (gauges%gauge_y)
      if (allocated(gauges%ix0)) deallocate (gauges%ix0)
      if (allocated(gauges%iy0)) deallocate (gauges%iy0)
      if (allocated(gauges%wx)) deallocate (gauges%wx)
      if (allocated(gauges%wy)) deallocate (gauges%wy)
      if (allocated(gauges%local_idx)) deallocate (gauges%local_idx)

      gauges%active = .false.
      gauges%n_gauges = 0
      gauges%n_local = 0
      gauges%time_index = 0
   end subroutine gauges_destroy

end module rki_gauges
