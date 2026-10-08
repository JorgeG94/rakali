!! Atmospheric forcing from NetCDF files (GFS/ERA5 style)
module rki_forcing
   !! Reads time-varying atmospheric forcing from a CF-convention NetCDF
   !! file (dims lon, lat, time) and interpolates onto the model grid
   !! (bilinear in space, linear in time). Fields:
   !!   - 10m wind (u10, v10) → surface stress (Pa)
   !!   - Mean sea-level pressure (msl or sp) → Pa
   !!   - Total precipitation (tp) → m/s
   use rki_constants, only: wp
   use rki_config, only: config_t
   use rki_state, only: state_t
   use rki_io_netcdf, only: nc_check, nc_open_read, nc_close, &
                            nc_get_dim_len, nc_get_varid
   use netcdf, only: nf90_get_var, nf90_noerr, nf90_inq_varid
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: forcing_t
   public :: forcing_init
   public :: forcing_init_from_config
   public :: forcing_update
   public :: forcing_destroy

   real(wp), parameter :: RHO_AIR = 1.225_wp
      !! Air density (kg/m^3) for wind stress conversion

   type :: forcing_t
      !! Container for atmospheric forcing data and interpolation state

      ! NetCDF file handle
      integer :: ncid = -1
         !! NetCDF file ID (kept open for lazy reads)
      logical :: active = .false.
         !! Whether forcing is enabled

      ! Forcing grid coordinates
      integer :: nlon = 0
         !! Number of longitude points in forcing file
      integer :: nlat = 0
         !! Number of latitude points in forcing file
      integer :: ntime = 0
         !! Number of time steps in forcing file
      real(wp), allocatable :: lon(:)
         !! Longitude coordinates (degrees or metres)
      real(wp), allocatable :: lat(:)
         !! Latitude coordinates (degrees or metres)
      real(wp), allocatable :: time(:)
         !! Time coordinates (seconds since reference)

      ! Variable IDs (-1 = not present)
      integer :: varid_u10 = -1
      integer :: varid_v10 = -1
      integer :: varid_msl = -1
      integer :: varid_tp = -1

      ! Two time-bracketing snapshots on the forcing grid
      integer :: idx_lo = 0
         !! Time index of the earlier snapshot (1-based)
      integer :: idx_hi = 0
         !! Time index of the later snapshot (1-based)
      real(wp), allocatable :: u10_lo(:, :), u10_hi(:, :)
         !! 10m wind x-component (m/s)
      real(wp), allocatable :: v10_lo(:, :), v10_hi(:, :)
         !! 10m wind y-component (m/s)
      real(wp), allocatable :: msl_lo(:, :), msl_hi(:, :)
         !! Atmospheric pressure (Pa)
      real(wp), allocatable :: tp_lo(:, :), tp_hi(:, :)
         !! Precipitation rate (m/s)

      ! Precomputed bilinear interpolation indices and weights, per model cell
      integer, allocatable :: ix0(:, :)
         !! Lower-left forcing x-index for each model cell
      integer, allocatable :: iy0(:, :)
         !! Lower-left forcing y-index for each model cell
      real(wp), allocatable :: wx(:, :)
         !! x-direction weight (0-1)
      real(wp), allocatable :: wy(:, :)
         !! y-direction weight (0-1)

      ! Model grid origin (for coordinate mapping)
      real(wp) :: x_origin = 0.0_wp
         !! Model grid x-origin (same units as forcing lon)
      real(wp) :: y_origin = 0.0_wp
         !! Model grid y-origin (same units as forcing lat)

   end type forcing_t

contains

   subroutine forcing_init(forc, filename, state, x_origin, y_origin)
      !! Initialise forcing from a NetCDF file: open it, read grid
      !! coordinates, precompute interpolation weights, load the first two
      !! snapshots. Public only for the unit-test suite.
      type(forcing_t), intent(inout) :: forc
         !! Forcing state to initialise
      character(len=*), intent(in) :: filename
         !! Path to NetCDF forcing file
      type(state_t), intent(in) :: state
         !! Model state (for grid dimensions)
      real(wp), intent(in), optional :: x_origin
         !! Model grid x-origin in forcing coordinates
      real(wp), intent(in), optional :: y_origin
         !! Model grid y-origin in forcing coordinates

      integer :: ncid, status

      if (len_trim(filename) == 0) return

      forc%x_origin = 0.0_wp
      forc%y_origin = 0.0_wp
      if (present(x_origin)) forc%x_origin = x_origin
      if (present(y_origin)) forc%y_origin = y_origin

      call logger%info("Loading atmospheric forcing from: "//trim(filename))

      call nc_open_read(filename, ncid)
      forc%ncid = ncid

      ! Read grid dimensions and coordinates
      call read_forcing_grid(forc, ncid)

      ! Find variable IDs (optional: missing vars get -1)
      call find_forcing_vars(forc, ncid)

      ! Allocate snapshot arrays on forcing grid
      call allocate_snapshots(forc)

      ! Precompute spatial interpolation weights
      call compute_interp_weights(forc, state)

      ! Load first two snapshots
      if (forc%ntime >= 2) then
         forc%idx_lo = 1
         forc%idx_hi = 2
         call load_snapshot(forc, ncid, 1, &
                            forc%u10_lo, forc%v10_lo, forc%msl_lo, forc%tp_lo)
         call load_snapshot(forc, ncid, 2, &
                            forc%u10_hi, forc%v10_hi, forc%msl_hi, forc%tp_hi)
      else if (forc%ntime == 1) then
         forc%idx_lo = 1
         forc%idx_hi = 1
         call load_snapshot(forc, ncid, 1, &
                            forc%u10_lo, forc%v10_lo, forc%msl_lo, forc%tp_lo)
         forc%u10_hi = forc%u10_lo
         forc%v10_hi = forc%v10_lo
         forc%msl_hi = forc%msl_lo
         forc%tp_hi = forc%tp_lo
      end if

      forc%active = .true.

      call logger%info("Forcing grid: "//to_string(forc%nlon)//" x "// &
                       to_string(forc%nlat)//", "// &
                       to_string(forc%ntime)//" time steps")
      call logger%info("Forcing time range: "//to_string(forc%time(1))// &
                       " to "//to_string(forc%time(forc%ntime))//" s")

   end subroutine forcing_init

   subroutine forcing_init_from_config(forc, cfg, state)
      !! Initialise forcing from config parameters
      type(forcing_t), intent(inout) :: forc
      type(config_t), intent(in) :: cfg
      type(state_t), intent(in) :: state

      call forcing_init(forc, trim(cfg%forcing_file), state, &
                        cfg%forcing_x_origin, cfg%forcing_y_origin)
   end subroutine forcing_init_from_config

   subroutine forcing_update(forc, state, sim_time)
      !! Interpolate forcing to `sim_time` and fill state%forcing%wind_x,
      !! wind_y, atm_pressure, precipitation (advances time bracket as needed).
      type(forcing_t), intent(inout) :: forc
         !! Forcing state
      type(state_t), intent(inout) :: state
         !! Model state with forcing arrays to fill
      real(wp), intent(in) :: sim_time
         !! Current simulation time (s)

      real(wp) :: alpha
      integer :: i, j, i0, j0
      real(wp) :: w00, w10, w01, w11
      real(wp) :: u10_val, v10_val, msl_val, tp_val
      real(wp) :: wind_speed, cd, tau_x, tau_y
      integer :: nx, ny, ng
      integer :: i0p, j0p

      if (.not. forc%active) return

      ! Advance time bracket if needed
      call advance_bracket(forc, sim_time)

      ! Compute temporal interpolation weight
      if (forc%idx_lo == forc%idx_hi) then
         alpha = 0.0_wp
      else
         alpha = (sim_time - forc%time(forc%idx_lo))/ &
                 (forc%time(forc%idx_hi) - forc%time(forc%idx_lo))
         alpha = max(0.0_wp, min(1.0_wp, alpha))
      end if

      nx = state%grid%nx_total
      ny = state%grid%ny_total
      ng = state%grid%nghost

      do j = ng + 1, ny - ng
         do i = ng + 1, nx - ng
            i0 = forc%ix0(i, j)
            j0 = forc%iy0(i, j)
            i0p = min(i0 + 1, forc%nlon)
            j0p = min(j0 + 1, forc%nlat)

            ! Bilinear weights
            w00 = (1.0_wp - forc%wx(i, j))*(1.0_wp - forc%wy(i, j))
            w10 = forc%wx(i, j)*(1.0_wp - forc%wy(i, j))
            w01 = (1.0_wp - forc%wx(i, j))*forc%wy(i, j)
            w11 = forc%wx(i, j)*forc%wy(i, j)

            ! Interpolate u10 in space, then time
            u10_val = bilinear_time(forc%u10_lo, forc%u10_hi, &
                                    i0, j0, i0p, j0p, &
                                    w00, w10, w01, w11, alpha)
            v10_val = bilinear_time(forc%v10_lo, forc%v10_hi, &
                                    i0, j0, i0p, j0p, &
                                    w00, w10, w01, w11, alpha)
            msl_val = bilinear_time(forc%msl_lo, forc%msl_hi, &
                                    i0, j0, i0p, j0p, &
                                    w00, w10, w01, w11, alpha)
            tp_val = bilinear_time(forc%tp_lo, forc%tp_hi, &
                                   i0, j0, i0p, j0p, &
                                   w00, w10, w01, w11, alpha)

            ! Convert 10m wind to stress: tau = rho_air * C_d * |W| * W
            wind_speed = sqrt(u10_val*u10_val + v10_val*v10_val)
            cd = drag_coefficient(wind_speed)
            tau_x = RHO_AIR*cd*wind_speed*u10_val
            tau_y = RHO_AIR*cd*wind_speed*v10_val

            state%forcing%wind_x(i, j) = tau_x
            state%forcing%wind_y(i, j) = tau_y
            state%forcing%atm_pressure(i, j) = msl_val
            state%forcing%precipitation(i, j) = tp_val
         end do
      end do

   end subroutine forcing_update

   subroutine forcing_destroy(forc)
      !! Clean up forcing state and close NetCDF file
      type(forcing_t), intent(inout) :: forc

      if (forc%ncid >= 0) then
         call nc_close(forc%ncid)
         forc%ncid = -1
      end if

      if (allocated(forc%lon)) deallocate (forc%lon)
      if (allocated(forc%lat)) deallocate (forc%lat)
      if (allocated(forc%time)) deallocate (forc%time)
      if (allocated(forc%u10_lo)) deallocate (forc%u10_lo)
      if (allocated(forc%u10_hi)) deallocate (forc%u10_hi)
      if (allocated(forc%v10_lo)) deallocate (forc%v10_lo)
      if (allocated(forc%v10_hi)) deallocate (forc%v10_hi)
      if (allocated(forc%msl_lo)) deallocate (forc%msl_lo)
      if (allocated(forc%msl_hi)) deallocate (forc%msl_hi)
      if (allocated(forc%tp_lo)) deallocate (forc%tp_lo)
      if (allocated(forc%tp_hi)) deallocate (forc%tp_hi)
      if (allocated(forc%ix0)) deallocate (forc%ix0)
      if (allocated(forc%iy0)) deallocate (forc%iy0)
      if (allocated(forc%wx)) deallocate (forc%wx)
      if (allocated(forc%wy)) deallocate (forc%wy)

      forc%active = .false.

   end subroutine forcing_destroy

   ! ================================================================
   ! Private helpers
   ! ================================================================

   subroutine read_forcing_grid(forc, ncid)
      !! Read lon, lat, time coordinate arrays from the forcing file
      type(forcing_t), intent(inout) :: forc
      integer, intent(in) :: ncid

      integer :: varid

      ! Read dimensions
      call nc_get_dim_len(ncid, "lon", forc%nlon)
      call nc_get_dim_len(ncid, "lat", forc%nlat)
      call nc_get_dim_len(ncid, "time", forc%ntime)

      ! Read coordinate arrays
      allocate (forc%lon(forc%nlon))
      allocate (forc%lat(forc%nlat))
      allocate (forc%time(forc%ntime))

      call nc_get_varid(ncid, "lon", varid)
      call nc_check(nf90_get_var(ncid, varid, forc%lon), "reading lon")

      call nc_get_varid(ncid, "lat", varid)
      call nc_check(nf90_get_var(ncid, varid, forc%lat), "reading lat")

      call nc_get_varid(ncid, "time", varid)
      call nc_check(nf90_get_var(ncid, varid, forc%time), "reading time")

   end subroutine read_forcing_grid

   subroutine find_forcing_vars(forc, ncid)
      !! Locate forcing variable IDs, tolerating missing variables
      type(forcing_t), intent(inout) :: forc
      integer, intent(in) :: ncid

      integer :: status

      ! Wind components
      status = nf90_inq_varid(ncid, "u10", forc%varid_u10)
      if (status /= nf90_noerr) forc%varid_u10 = -1

      status = nf90_inq_varid(ncid, "v10", forc%varid_v10)
      if (status /= nf90_noerr) forc%varid_v10 = -1

      ! Pressure: try "msl" first, then "sp"
      status = nf90_inq_varid(ncid, "msl", forc%varid_msl)
      if (status /= nf90_noerr) then
         status = nf90_inq_varid(ncid, "sp", forc%varid_msl)
         if (status /= nf90_noerr) forc%varid_msl = -1
      end if

      ! Precipitation
      status = nf90_inq_varid(ncid, "tp", forc%varid_tp)
      if (status /= nf90_noerr) forc%varid_tp = -1

      if (forc%varid_u10 >= 0) then
         call logger%info("  Found wind component: u10")
      end if
      if (forc%varid_v10 >= 0) then
         call logger%info("  Found wind component: v10")
      end if
      if (forc%varid_msl >= 0) then
         call logger%info("  Found pressure variable")
      end if
      if (forc%varid_tp >= 0) then
         call logger%info("  Found precipitation: tp")
      end if

   end subroutine find_forcing_vars

   subroutine allocate_snapshots(forc)
      !! Allocate the two time-bracket snapshot arrays
      type(forcing_t), intent(inout) :: forc

      integer :: nl, na

      nl = forc%nlon
      na = forc%nlat

      allocate (forc%u10_lo(nl, na), source=0.0_wp)
      allocate (forc%u10_hi(nl, na), source=0.0_wp)
      allocate (forc%v10_lo(nl, na), source=0.0_wp)
      allocate (forc%v10_hi(nl, na), source=0.0_wp)
      allocate (forc%msl_lo(nl, na), source=0.0_wp)
      allocate (forc%msl_hi(nl, na), source=0.0_wp)
      allocate (forc%tp_lo(nl, na), source=0.0_wp)
      allocate (forc%tp_hi(nl, na), source=0.0_wp)

   end subroutine allocate_snapshots

   subroutine compute_interp_weights(forc, state)
      !! Precompute bilinear interpolation indices and weights: map each
      !! model cell centre to the forcing grid, store lower-left index and
      !! fractional weights.
      type(forcing_t), intent(inout) :: forc
      type(state_t), intent(in) :: state

      integer :: i, j, nx, ny, ng
      real(wp) :: x_model, y_model, fx, fy

      nx = state%grid%nx_total
      ny = state%grid%ny_total
      ng = state%grid%nghost

      allocate (forc%ix0(nx, ny), source=1)
      allocate (forc%iy0(nx, ny), source=1)
      allocate (forc%wx(nx, ny), source=0.0_wp)
      allocate (forc%wy(nx, ny), source=0.0_wp)

      do j = ng + 1, ny - ng
         do i = ng + 1, nx - ng
            ! Model cell centre in forcing coordinates
            x_model = forc%x_origin + (real(i - ng, wp) - 0.5_wp)*state%grid%dx
            y_model = forc%y_origin + (real(j - ng, wp) - 0.5_wp)*state%grid%dy

            ! Find position in forcing grid
            call find_bracket(forc%lon, forc%nlon, x_model, forc%ix0(i, j), fx)
            call find_bracket(forc%lat, forc%nlat, y_model, forc%iy0(i, j), fy)

            forc%wx(i, j) = fx
            forc%wy(i, j) = fy
         end do
      end do

   end subroutine compute_interp_weights

   pure subroutine find_bracket(coords, n, val, idx, frac)
      !! Find the bracketing interval in a sorted coordinate array:
      !! idx with coords(idx) <= val < coords(idx+1) and
      !! frac = (val - coords(idx)) / (coords(idx+1) - coords(idx)).
      !! Clamps to boundaries if val is outside the grid.
      real(wp), intent(in) :: coords(:)
      integer, intent(in) :: n
      real(wp), intent(in) :: val
      integer, intent(out) :: idx
      real(wp), intent(out) :: frac

      integer :: k

      if (n <= 1) then
         idx = 1
         frac = 0.0_wp
         return
      end if

      ! Clamp below
      if (val <= coords(1)) then
         idx = 1
         frac = 0.0_wp
         return
      end if

      ! Clamp above
      if (val >= coords(n)) then
         idx = n
         frac = 0.0_wp
         return
      end if

      ! Linear search (forcing grids are typically small, <1000 points)
      do k = 1, n - 1
         if (val >= coords(k) .and. val < coords(k + 1)) then
            idx = k
            frac = (val - coords(k))/(coords(k + 1) - coords(k))
            return
         end if
      end do

      ! Fallback (should not reach here)
      idx = n - 1
      frac = 1.0_wp

   end subroutine find_bracket

   subroutine load_snapshot(forc, ncid, tidx, u10, v10, msl, tp)
      !! Read one time snapshot of all forcing fields from NetCDF
      type(forcing_t), intent(in) :: forc
      integer, intent(in) :: ncid
      integer, intent(in) :: tidx
         !! 1-based time index to read
      real(wp), intent(out) :: u10(:, :), v10(:, :), msl(:, :), tp(:, :)

      integer :: start3(3), count3(3)

      start3 = [1, 1, tidx]
      count3 = [forc%nlon, forc%nlat, 1]

      u10 = 0.0_wp
      v10 = 0.0_wp
      msl = 0.0_wp
      tp = 0.0_wp

      if (forc%varid_u10 >= 0) then
         call nc_check(nf90_get_var(ncid, forc%varid_u10, u10, &
                                    start=start3, count=count3), "reading u10")
      end if

      if (forc%varid_v10 >= 0) then
         call nc_check(nf90_get_var(ncid, forc%varid_v10, v10, &
                                    start=start3, count=count3), "reading v10")
      end if

      if (forc%varid_msl >= 0) then
         call nc_check(nf90_get_var(ncid, forc%varid_msl, msl, &
                                    start=start3, count=count3), "reading msl/sp")
      end if

      if (forc%varid_tp >= 0) then
         call nc_check(nf90_get_var(ncid, forc%varid_tp, tp, &
                                    start=start3, count=count3), "reading tp")
      end if

   end subroutine load_snapshot

   subroutine advance_bracket(forc, sim_time)
      !! Advance the time bracket so that time(idx_lo) <= sim_time < time(idx_hi)
      !!
      !! Reads new snapshots from disk as needed.
      type(forcing_t), intent(inout) :: forc
      real(wp), intent(in) :: sim_time

      ! Advance forward if sim_time has passed the current bracket
      do while (forc%idx_hi < forc%ntime .and. &
                sim_time >= forc%time(forc%idx_hi))
         ! Shift hi -> lo
         forc%idx_lo = forc%idx_hi
         forc%u10_lo = forc%u10_hi
         forc%v10_lo = forc%v10_hi
         forc%msl_lo = forc%msl_hi
         forc%tp_lo = forc%tp_hi

         ! Load new hi
         forc%idx_hi = forc%idx_hi + 1
         call load_snapshot(forc, forc%ncid, forc%idx_hi, &
                            forc%u10_hi, forc%v10_hi, &
                            forc%msl_hi, forc%tp_hi)
      end do

   end subroutine advance_bracket

   pure function bilinear_time(fld_lo, fld_hi, i0, j0, i1, j1, &
                               w00, w10, w01, w11, alpha) result(val)
      !! Bilinear spatial interpolation + linear temporal interpolation
      real(wp), intent(in) :: fld_lo(:, :), fld_hi(:, :)
      integer, intent(in) :: i0, j0, i1, j1
      real(wp), intent(in) :: w00, w10, w01, w11
      real(wp), intent(in) :: alpha
      real(wp) :: val

      real(wp) :: val_lo, val_hi

      val_lo = w00*fld_lo(i0, j0) + w10*fld_lo(i1, j0) + &
               w01*fld_lo(i0, j1) + w11*fld_lo(i1, j1)
      val_hi = w00*fld_hi(i0, j0) + w10*fld_hi(i1, j0) + &
               w01*fld_hi(i0, j1) + w11*fld_hi(i1, j1)

      val = (1.0_wp - alpha)*val_lo + alpha*val_hi

   end function bilinear_time

   pure function drag_coefficient(wind_speed) result(cd)
      !! Large & Pond (1981) drag coefficient for 10m wind, capped at
      !! 2.5e-3 for extreme winds.
      real(wp), intent(in) :: wind_speed
         !! 10m wind speed magnitude (m/s)
      real(wp) :: cd

      if (wind_speed < 11.0_wp) then
         cd = 1.2e-3_wp
      else
         cd = (0.49_wp + 0.065_wp*wind_speed)*1.0e-3_wp
         cd = min(cd, 2.5e-3_wp)
      end if

   end function drag_coefficient

end module rki_forcing
