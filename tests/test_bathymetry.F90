!! Unit tests for bathymetry NetCDF loading
module test_bathymetry
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t, load_bathymetry, fill_bathymetry_ghosts
   use netcdf, only: nf90_create, nf90_close, nf90_def_dim, nf90_def_var, &
                     nf90_enddef, nf90_put_var, nf90_clobber, nf90_netcdf4, &
                     nf90_double, nf90_float
   implicit none
   private

   public :: collect_bathymetry_tests

   integer, parameter :: NX = 6, NY = 4, NG = 2
   real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp

contains

   subroutine collect_bathymetry_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("load_b_round_trip", test_load_b), &
                  new_unittest("ghost_cells_extrapolated", test_ghost_extrapolation), &
                  new_unittest("load_elevation_variable", test_load_elevation), &
                  new_unittest("load_depth_variable", test_load_depth), &
                  new_unittest("fill_bathymetry_ghosts_direct", test_fill_ghosts_direct) &
                  ]
   end subroutine collect_bathymetry_tests

   subroutine delete_file(filename)
      character(len=*), intent(in) :: filename
      integer :: u
      logical :: exists
      inquire (file=trim(filename), exist=exists)
      if (exists) then
         open (newunit=u, file=trim(filename), status='old')
         close (u, status='delete')
      end if
   end subroutine delete_file

   ! Write a NetCDF file with dims (x, y) and a variable named `varname`.
   subroutine write_bathy_file(filename, varname, data)
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      real(wp), intent(in) :: data(:, :)
      integer :: ncid, x_dimid, y_dimid, varid, status
      integer :: nc_type
#ifdef RAKALI_DOUBLE_PRECISION
      nc_type = nf90_double
#else
      nc_type = nf90_float
#endif
      ! Force double for round-trip exactness regardless of build precision —
      ! storage type only needs to match the read precision.
      if (wp == kind(1.0d0)) then
         nc_type = nf90_double
      else
         nc_type = nf90_float
      end if

      status = nf90_create(trim(filename), &
                           ior(nf90_clobber, nf90_netcdf4), ncid)
      status = nf90_def_dim(ncid, "x", size(data, 1), x_dimid)
      status = nf90_def_dim(ncid, "y", size(data, 2), y_dimid)
      status = nf90_def_var(ncid, varname, nc_type, [x_dimid, y_dimid], varid)
      status = nf90_enddef(ncid)
      status = nf90_put_var(ncid, varid, data)
      status = nf90_close(ncid)
   end subroutine write_bathy_file

   subroutine test_load_b(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      character(len=*), parameter :: fname = "test_bathy_b.nc"
      real(wp) :: bathy(NX, NY)
      integer :: i, j
      checks: block

         ! Make a simple ramp: b(i,j) = i + 0.1*j
         do j = 1, NY
            do i = 1, NX
               bathy(i, j) = real(i, wp) + 0.1_wp*real(j, wp)
            end do
         end do
         call write_bathy_file(fname, "b", bathy)

         call s%init(NX, NY, NG, DX, DY, 0.45_wp)
         call load_bathymetry(fname, s)

         ! Interior round-trip
         do j = 1, NY
            do i = 1, NX
               call check(error, abs(s%barotropic%b(NG + i, NG + j) - bathy(i, j)) &
                          < 1.0e-6_wp, "interior bathymetry mismatch")
               if (allocated(error)) exit checks
            end do
         end do

      end block checks
      continue
      call s%destroy()
      call delete_file(fname)
   end subroutine test_load_b

   subroutine test_ghost_extrapolation(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      character(len=*), parameter :: fname = "test_bathy_ghosts.nc"
      real(wp) :: bathy(NX, NY)
      integer :: i, j
      checks: block

         do j = 1, NY
            do i = 1, NX
               bathy(i, j) = 10.0_wp*real(i, wp) + real(j, wp)
            end do
         end do
         call write_bathy_file(fname, "b", bathy)

         call s%init(NX, NY, NG, DX, DY, 0.45_wp)
         call load_bathymetry(fname, s)

         ! West ghost columns equal first interior column
         do j = 1, s%grid%ny_total
            do i = 1, NG
               call check(error, &
                          abs(s%barotropic%b(i, j) - s%barotropic%b(NG + 1, j)) &
                          < 1.0e-12_wp, "west ghost should equal first interior column")
               if (allocated(error)) exit checks
            end do
         end do

         ! East ghost columns equal last interior column
         do j = 1, s%grid%ny_total
            do i = 1, NG
               call check(error, &
                          abs(s%barotropic%b(NG + NX + i, j) - s%barotropic%b(NG + NX, j)) &
                          < 1.0e-12_wp, "east ghost should equal last interior column")
               if (allocated(error)) exit checks
            end do
         end do

         ! South ghost rows equal first interior row
         do j = 1, NG
            do i = 1, s%grid%nx_total
               call check(error, &
                          abs(s%barotropic%b(i, j) - s%barotropic%b(i, NG + 1)) &
                          < 1.0e-12_wp, "south ghost should equal first interior row")
               if (allocated(error)) exit checks
            end do
         end do

         ! North ghost rows equal last interior row
         do j = 1, NG
            do i = 1, s%grid%nx_total
               call check(error, &
                          abs(s%barotropic%b(i, NG + NY + j) - s%barotropic%b(i, NG + NY)) &
                          < 1.0e-12_wp, "north ghost should equal last interior row")
               if (allocated(error)) exit checks
            end do
         end do

      end block checks
      continue
      call s%destroy()
      call delete_file(fname)
   end subroutine test_ghost_extrapolation

   subroutine test_load_elevation(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      character(len=*), parameter :: fname = "test_bathy_elev.nc"
      real(wp) :: bathy(NX, NY)
      integer :: i, j

      bathy = 0.0_wp
      do j = 1, NY
         do i = 1, NX
            bathy(i, j) = -5.0_wp + 0.5_wp*real(i + j, wp)
         end do
      end do
      call write_bathy_file(fname, "elevation", bathy)

      call s%init(NX, NY, NG, DX, DY, 0.45_wp)
      call load_bathymetry(fname, s)

      call check(error, &
                 abs(s%barotropic%b(NG + 3, NG + 2) - bathy(3, 2)) < 1.0e-6_wp, &
                 "loading via 'elevation' name failed")

      call s%destroy()
      call delete_file(fname)
   end subroutine test_load_elevation

   subroutine test_load_depth(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      character(len=*), parameter :: fname = "test_bathy_depth.nc"
      real(wp) :: bathy(NX, NY)
      integer :: i, j

      do j = 1, NY
         do i = 1, NX
            bathy(i, j) = -100.0_wp + real(i*j, wp)
         end do
      end do
      call write_bathy_file(fname, "depth", bathy)

      call s%init(NX, NY, NG, DX, DY, 0.45_wp)
      call load_bathymetry(fname, s)

      call check(error, &
                 abs(s%barotropic%b(NG + 4, NG + 3) - bathy(4, 3)) < 1.0e-6_wp, &
                 "loading via 'depth' name failed")

      call s%destroy()
      call delete_file(fname)
   end subroutine test_load_depth

   subroutine test_fill_ghosts_direct(error)
      !! Exercise fill_bathymetry_ghosts directly without round-tripping through
      !! NetCDF (covers a callable that's also reused by other modules).
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      integer :: i, j
      checks: block

         call s%init(NX, NY, NG, DX, DY, 0.45_wp)

         ! Reset all to a poison value, set interior only
         s%barotropic%b = -999.0_wp
         do j = 1, NY
            do i = 1, NX
               s%barotropic%b(NG + i, NG + j) = real(100*i + j, wp)
            end do
         end do

         call fill_bathymetry_ghosts(s)

         ! West/east edges must now equal nearest interior column
         do j = 1, s%grid%ny_total
            call check(error, abs(s%barotropic%b(1, j) - s%barotropic%b(NG + 1, j)) &
                       < 1.0e-12_wp, "west ghost extrapolation failed")
            if (allocated(error)) exit checks
         end do

         do i = 1, s%grid%nx_total
            call check(error, abs(s%barotropic%b(i, 1) - s%barotropic%b(i, NG + 1)) &
                       < 1.0e-12_wp, "south ghost extrapolation failed")
            if (allocated(error)) exit checks
         end do

      end block checks
      continue
      call s%destroy()
   end subroutine test_fill_ghosts_direct

end module test_bathymetry
