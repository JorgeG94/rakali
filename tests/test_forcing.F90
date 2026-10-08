!! Unit tests for the atmospheric forcing module
module test_forcing
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_forcing, only: forcing_t, forcing_init, forcing_update, forcing_destroy
   use netcdf, only: nf90_create, nf90_def_dim, nf90_def_var, nf90_enddef, &
                     nf90_put_var, nf90_close, nf90_clobber, nf90_netcdf4, nf90_double, &
                     nf90_noerr
   implicit none
   private

   public :: collect_forcing_tests

contains

   subroutine collect_forcing_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("init_destroy_lifecycle", test_init_destroy), &
                  new_unittest("uniform_wind_stress", test_uniform_wind), &
                  new_unittest("uniform_pressure", test_uniform_pressure), &
                  new_unittest("uniform_precipitation", test_uniform_precip), &
                  new_unittest("time_interpolation", test_time_interp), &
                  new_unittest("spatial_interpolation", test_spatial_interp), &
                  new_unittest("no_forcing_file_noop", test_no_file) &
                  ]
   end subroutine collect_forcing_tests

   ! ================================================================
   ! Helper: create a minimal NetCDF forcing file
   ! ================================================================

   subroutine create_forcing_file(filename, nlon, nlat, ntime, &
                                  lon_vals, lat_vals, time_vals, &
                                  u10_vals, v10_vals, msl_vals, tp_vals)
      !! Create a test NetCDF file with forcing data
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nlon, nlat, ntime
      real(wp), intent(in) :: lon_vals(nlon)
      real(wp), intent(in) :: lat_vals(nlat)
      real(wp), intent(in) :: time_vals(ntime)
      real(wp), intent(in) :: u10_vals(nlon, nlat, ntime)
      real(wp), intent(in) :: v10_vals(nlon, nlat, ntime)
      real(wp), intent(in) :: msl_vals(nlon, nlat, ntime)
      real(wp), intent(in) :: tp_vals(nlon, nlat, ntime)

      integer :: ncid, dim_lon, dim_lat, dim_time
      integer :: var_lon, var_lat, var_time
      integer :: var_u10, var_v10, var_msl, var_tp
      integer :: status, dimids(3)

      status = nf90_create(trim(filename), ior(nf90_clobber, nf90_netcdf4), ncid)

      ! Define dimensions
      status = nf90_def_dim(ncid, "lon", nlon, dim_lon)
      status = nf90_def_dim(ncid, "lat", nlat, dim_lat)
      status = nf90_def_dim(ncid, "time", ntime, dim_time)

      ! Define coordinate variables
      status = nf90_def_var(ncid, "lon", nf90_double, [dim_lon], var_lon)
      status = nf90_def_var(ncid, "lat", nf90_double, [dim_lat], var_lat)
      status = nf90_def_var(ncid, "time", nf90_double, [dim_time], var_time)

      ! Define data variables (lon, lat, time)
      dimids = [dim_lon, dim_lat, dim_time]
      status = nf90_def_var(ncid, "u10", nf90_double, dimids, var_u10)
      status = nf90_def_var(ncid, "v10", nf90_double, dimids, var_v10)
      status = nf90_def_var(ncid, "msl", nf90_double, dimids, var_msl)
      status = nf90_def_var(ncid, "tp", nf90_double, dimids, var_tp)

      status = nf90_enddef(ncid)

      ! Write data
      status = nf90_put_var(ncid, var_lon, lon_vals)
      status = nf90_put_var(ncid, var_lat, lat_vals)
      status = nf90_put_var(ncid, var_time, time_vals)
      status = nf90_put_var(ncid, var_u10, u10_vals)
      status = nf90_put_var(ncid, var_v10, v10_vals)
      status = nf90_put_var(ncid, var_msl, msl_vals)
      status = nf90_put_var(ncid, var_tp, tp_vals)

      status = nf90_close(ncid)

   end subroutine create_forcing_file

   ! ================================================================
   ! Tests
   ! ================================================================

   subroutine test_init_destroy(error)
      !! Forcing module init/destroy should not crash with valid file
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      type(forcing_t) :: forc
      character(len=*), parameter :: FNAME = "test_forc_lifecycle.nc"
      integer, parameter :: NL = 3, NA = 3, NT = 2
      real(wp) :: lon(NL), lat(NA), time(NT)
      real(wp) :: u10(NL, NA, NT), v10(NL, NA, NT)
      real(wp) :: msl(NL, NA, NT), tp(NL, NA, NT)

      lon = [0.0_wp, 50.0_wp, 100.0_wp]
      lat = [0.0_wp, 50.0_wp, 100.0_wp]
      time = [0.0_wp, 3600.0_wp]
      u10 = 5.0_wp; v10 = 0.0_wp; msl = 101325.0_wp; tp = 0.0_wp

      call create_forcing_file(FNAME, NL, NA, NT, lon, lat, time, &
                               u10, v10, msl, tp)

      call state%init(10, 10, 2, 10.0_wp, 10.0_wp, 0.45_wp)
      call forcing_init(forc, FNAME, state)

      call check(error, forc%active, "Forcing should be active after init")
      if (allocated(error)) return

      call check(error, forc%nlon == NL, "nlon should match file")
      if (allocated(error)) return

      call check(error, forc%nlat == NA, "nlat should match file")
      if (allocated(error)) return

      call check(error, forc%ntime == NT, "ntime should match file")
      if (allocated(error)) return

      call forcing_destroy(forc)
      call check(error,.not. forc%active, "Forcing should be inactive after destroy")

      call state%destroy()

   end subroutine test_init_destroy

   subroutine test_uniform_wind(error)
      !! Uniform 10 m/s wind should produce correct stress via drag formula
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      type(forcing_t) :: forc
      character(len=*), parameter :: FNAME = "test_forc_wind.nc"
      integer, parameter :: NL = 2, NA = 2, NT = 1
      real(wp) :: lon(NL), lat(NA), time(NT)
      real(wp) :: u10(NL, NA, NT), v10(NL, NA, NT)
      real(wp) :: msl(NL, NA, NT), tp(NL, NA, NT)
      real(wp) :: expected_stress
      real(wp), parameter :: RHO_AIR = 1.225_wp
      real(wp), parameter :: CD_LOW = 1.2e-3_wp
      real(wp), parameter :: WIND_SPEED = 5.0_wp
      integer :: ic, jc

      ! Forcing grid covers the model domain
      lon = [0.0_wp, 200.0_wp]
      lat = [0.0_wp, 200.0_wp]
      time = [0.0_wp]
      u10 = WIND_SPEED
      v10 = 0.0_wp
      msl = 101325.0_wp
      tp = 0.0_wp

      call create_forcing_file(FNAME, NL, NA, NT, lon, lat, time, &
                               u10, v10, msl, tp)

      call state%init(10, 10, 2, 10.0_wp, 10.0_wp, 0.45_wp)
      call forcing_init(forc, FNAME, state)
      call forcing_update(forc, state, 0.0_wp)

      ! tau = rho_air * Cd * |W| * W = 1.225 * 1.2e-3 * 5 * 5
      expected_stress = RHO_AIR*CD_LOW*WIND_SPEED*WIND_SPEED

      ! Check interior cell
      ic = state%grid%nghost + state%grid%nx_phys/2
      jc = state%grid%nghost + state%grid%ny_phys/2

      call check(error, abs(state%forcing%wind_x(ic, jc) - expected_stress) < 1.0e-10_wp, &
                 "wind_x should match bulk drag formula")
      if (allocated(error)) return

      call check(error, abs(state%forcing%wind_y(ic, jc)) < 1.0e-14_wp, &
                 "wind_y should be zero for x-only wind")

      call forcing_destroy(forc)
      call state%destroy()

   end subroutine test_uniform_wind

   subroutine test_uniform_pressure(error)
      !! Uniform pressure should fill atm_pressure array correctly
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      type(forcing_t) :: forc
      character(len=*), parameter :: FNAME = "test_forc_pressure.nc"
      integer, parameter :: NL = 2, NA = 2, NT = 1
      real(wp) :: lon(NL), lat(NA), time(NT)
      real(wp) :: u10(NL, NA, NT), v10(NL, NA, NT)
      real(wp) :: msl(NL, NA, NT), tp(NL, NA, NT)
      real(wp), parameter :: P_ATM = 100800.0_wp
      integer :: ic, jc

      lon = [0.0_wp, 200.0_wp]
      lat = [0.0_wp, 200.0_wp]
      time = [0.0_wp]
      u10 = 0.0_wp; v10 = 0.0_wp; msl = P_ATM; tp = 0.0_wp

      call create_forcing_file(FNAME, NL, NA, NT, lon, lat, time, &
                               u10, v10, msl, tp)

      call state%init(10, 10, 2, 10.0_wp, 10.0_wp, 0.45_wp)
      call forcing_init(forc, FNAME, state)
      call forcing_update(forc, state, 0.0_wp)

      ic = state%grid%nghost + state%grid%nx_phys/2
      jc = state%grid%nghost + state%grid%ny_phys/2

      call check(error, abs(state%forcing%atm_pressure(ic, jc) - P_ATM) < 1.0e-6_wp, &
                 "atm_pressure should match forcing file value")

      call forcing_destroy(forc)
      call state%destroy()

   end subroutine test_uniform_pressure

   subroutine test_uniform_precip(error)
      !! Uniform precipitation should fill precipitation array correctly
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      type(forcing_t) :: forc
      character(len=*), parameter :: FNAME = "test_forc_precip.nc"
      integer, parameter :: NL = 2, NA = 2, NT = 1
      real(wp) :: lon(NL), lat(NA), time(NT)
      real(wp) :: u10(NL, NA, NT), v10(NL, NA, NT)
      real(wp) :: msl(NL, NA, NT), tp(NL, NA, NT)
      real(wp), parameter :: RAIN = 0.001_wp  ! 1 mm/s
      integer :: ic, jc

      lon = [0.0_wp, 200.0_wp]
      lat = [0.0_wp, 200.0_wp]
      time = [0.0_wp]
      u10 = 0.0_wp; v10 = 0.0_wp; msl = 101325.0_wp; tp = RAIN

      call create_forcing_file(FNAME, NL, NA, NT, lon, lat, time, &
                               u10, v10, msl, tp)

      call state%init(10, 10, 2, 10.0_wp, 10.0_wp, 0.45_wp)
      call forcing_init(forc, FNAME, state)
      call forcing_update(forc, state, 0.0_wp)

      ic = state%grid%nghost + state%grid%nx_phys/2
      jc = state%grid%nghost + state%grid%ny_phys/2

      call check(error, abs(state%forcing%precipitation(ic, jc) - RAIN) < 1.0e-14_wp, &
                 "precipitation should match forcing file value")

      call forcing_destroy(forc)
      call state%destroy()

   end subroutine test_uniform_precip

   subroutine test_time_interp(error)
      !! Linear time interpolation between two snapshots
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      type(forcing_t) :: forc
      character(len=*), parameter :: FNAME = "test_forc_time.nc"
      integer, parameter :: NL = 2, NA = 2, NT = 2
      real(wp) :: lon(NL), lat(NA), time(NT)
      real(wp) :: u10(NL, NA, NT), v10(NL, NA, NT)
      real(wp) :: msl(NL, NA, NT), tp(NL, NA, NT)
      real(wp) :: expected_msl
      integer :: ic, jc

      lon = [0.0_wp, 200.0_wp]
      lat = [0.0_wp, 200.0_wp]
      time = [0.0_wp, 3600.0_wp]
      u10 = 0.0_wp; v10 = 0.0_wp; tp = 0.0_wp
      ! Pressure: 101000 at t=0, 102000 at t=3600
      msl(:, :, 1) = 101000.0_wp
      msl(:, :, 2) = 102000.0_wp

      call create_forcing_file(FNAME, NL, NA, NT, lon, lat, time, &
                               u10, v10, msl, tp)

      call state%init(10, 10, 2, 10.0_wp, 10.0_wp, 0.45_wp)
      call forcing_init(forc, FNAME, state)

      ! At t = 1800 s (halfway), expect 101500 Pa
      call forcing_update(forc, state, 1800.0_wp)

      ic = state%grid%nghost + state%grid%nx_phys/2
      jc = state%grid%nghost + state%grid%ny_phys/2
      expected_msl = 101500.0_wp

      call check(error, abs(state%forcing%atm_pressure(ic, jc) - expected_msl) < 1.0e-6_wp, &
                 "Pressure should be linearly interpolated at t=1800s")
      if (allocated(error)) return

      ! At t = 900 s (quarter), expect 101250 Pa
      call forcing_update(forc, state, 900.0_wp)
      expected_msl = 101250.0_wp

      call check(error, abs(state%forcing%atm_pressure(ic, jc) - expected_msl) < 1.0e-6_wp, &
                 "Pressure should be linearly interpolated at t=900s")

      call forcing_destroy(forc)
      call state%destroy()

   end subroutine test_time_interp

   subroutine test_spatial_interp(error)
      !! Bilinear spatial interpolation with a linear pressure gradient
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      type(forcing_t) :: forc
      character(len=*), parameter :: FNAME = "test_forc_spatial.nc"
      integer, parameter :: NL = 3, NA = 3, NT = 1
      real(wp) :: lon(NL), lat(NA), time(NT)
      real(wp) :: u10(NL, NA, NT), v10(NL, NA, NT)
      real(wp) :: msl(NL, NA, NT), tp(NL, NA, NT)
      integer :: il, jl
      integer :: ic, jc
      real(wp) :: x_cell, expected_p

      ! Forcing grid: 3 points at x = 0, 50, 100
      lon = [0.0_wp, 50.0_wp, 100.0_wp]
      lat = [0.0_wp, 50.0_wp, 100.0_wp]
      time = [0.0_wp]
      u10 = 0.0_wp; v10 = 0.0_wp; tp = 0.0_wp

      ! Linear pressure gradient in x: p = 100000 + 10 * x
      do il = 1, NL
         do jl = 1, NA
            msl(il, jl, 1) = 100000.0_wp + 10.0_wp*lon(il)
         end do
      end do

      call create_forcing_file(FNAME, NL, NA, NT, lon, lat, time, &
                               u10, v10, msl, tp)

      ! Model grid: 10x10 cells, dx=dy=10m, covering 0-100m
      call state%init(10, 10, 2, 10.0_wp, 10.0_wp, 0.45_wp)
      call forcing_init(forc, FNAME, state)
      call forcing_update(forc, state, 0.0_wp)

      ! Check a cell at x = 45 m (cell index 5, centre at 45 m with 0-based)
      ic = state%grid%nghost + 5
      jc = state%grid%nghost + 5
      x_cell = (5.0_wp - 0.5_wp)*10.0_wp  ! = 45 m
      expected_p = 100000.0_wp + 10.0_wp*x_cell  ! = 100450

      call check(error, abs(state%forcing%atm_pressure(ic, jc) - expected_p) < 1.0_wp, &
                 "Bilinear interpolation should reproduce linear field")

      call forcing_destroy(forc)
      call state%destroy()

   end subroutine test_spatial_interp

   subroutine test_no_file(error)
      !! Empty forcing file path should leave forcing inactive
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      type(forcing_t) :: forc

      call state%init(10, 10, 2, 10.0_wp, 10.0_wp, 0.45_wp)
      call forcing_init(forc, "", state)

      call check(error,.not. forc%active, &
                 "Forcing should be inactive with empty filename")

      call forcing_destroy(forc)
      call state%destroy()

   end subroutine test_no_file

end module test_forcing
