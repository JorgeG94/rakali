!! Unit tests for the gauge output module
module test_gauges
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_config, only: config_t
   use rki_gauges, only: gauge_output_t, gauges_init, gauges_write, gauges_destroy
   use netcdf, only: nf90_open, nf90_close, nf90_nowrite, nf90_noerr, &
                     nf90_inq_varid, nf90_get_var, nf90_inq_dimid, &
                     nf90_inquire_dimension
   implicit none
   private

   public :: collect_gauges_tests

contains

   subroutine collect_gauges_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("no_gauge_file", test_no_gauge_file), &
                  new_unittest("gauge_file_parsing", test_gauge_file_parsing), &
                  new_unittest("bilinear_interpolation", test_bilinear_interpolation), &
                  new_unittest("write_multiple_snapshots", test_write_multiple_snapshots) &
                  ]
   end subroutine collect_gauges_tests

   ! ================================================================
   ! Helper: delete a file if it exists
   ! ================================================================

   subroutine delete_file(filename)
      character(len=*), intent(in) :: filename
      integer :: unit_num
      logical :: exists

      inquire (file=trim(filename), exist=exists)
      if (exists) then
         open (newunit=unit_num, file=trim(filename), status='old')
         close (unit_num, status='delete')
      end if
   end subroutine delete_file

   ! ================================================================
   ! Helper: create a gauge location text file
   ! ================================================================

   subroutine write_gauge_file(filename, names, x_coords, y_coords, n)
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: names(:)
      real(wp), intent(in) :: x_coords(:)
      real(wp), intent(in) :: y_coords(:)
      integer, intent(in) :: n

      integer :: unit_num, k

      open (newunit=unit_num, file=trim(filename), status='replace', action='write')
      write (unit_num, '(A)') '# Gauge locations for testing'
      do k = 1, n
         write (unit_num, '(A, 2(2X, F12.4))') trim(names(k)), x_coords(k), y_coords(k)
      end do
      close (unit_num)
   end subroutine write_gauge_file

   ! ================================================================
   ! Test 1: No gauge file — gauges should be inactive
   ! ================================================================

   subroutine test_no_gauge_file(error)
      type(error_type), allocatable, intent(out) :: error

      type(gauge_output_t) :: gauges
      type(config_t) :: cfg
      type(state_t) :: state

      ! Default config has gauge_file = "" (empty)
      call state%init(20, 10, 2, 1.0_wp, 1.0_wp, 0.45_wp)
      call gauges_init(gauges, cfg, state, 0)

      call check(error,.not. gauges%active, "gauges should be inactive with no file")
      if (allocated(error)) return

      call gauges_destroy(gauges)
      call state%destroy()

   end subroutine test_no_gauge_file

   ! ================================================================
   ! Test 2: Gauge file parsing — verify n_gauges, names, coordinates
   ! ================================================================

   subroutine test_gauge_file_parsing(error)
      type(error_type), allocatable, intent(out) :: error

      type(gauge_output_t) :: gauges
      type(config_t) :: cfg
      type(state_t) :: state
      character(len=64) :: names(2)
      real(wp) :: x_coords(2), y_coords(2)
      character(len=*), parameter :: gauge_filename = "test_gauge_parse.txt"
      character(len=*), parameter :: nc_filename = "gauges.nc"
      checks: block

         ! Create gauge file with 2 stations
         names(1) = "fort_denison"
         names(2) = "point_2"
         x_coords = [12.5_wp, 5.0_wp]
         y_coords = [8.3_wp, 3.0_wp]
         call write_gauge_file(gauge_filename, names, x_coords, y_coords, 2)

         ! Set up config
         cfg%gauge_file = gauge_filename
         cfg%gauge_output_interval = 1.0_wp
         cfg%output_dir = "."

         ! State: 20x10 cells, dx=dy=1.0, so domain is [0, 20] x [0, 10]
         call state%init(20, 10, 2, 1.0_wp, 1.0_wp, 0.45_wp)
         ! Set h to something nonzero so gauges are in wet cells
         state%barotropic%h(:, :) = 1.0_wp

         call gauges_init(gauges, cfg, state, 0)

         ! Verify parsing
         call check(error, gauges%active, "gauges should be active with valid file")
         if (allocated(error)) exit checks

         call check(error, gauges%n_gauges == 2, "expected 2 gauges")
         if (allocated(error)) exit checks

         call check(error, abs(gauges%gauge_x(1) - 12.5_wp) < 1.0e-10_wp, &
                    "gauge 1 x-coordinate mismatch")
         if (allocated(error)) exit checks

         call check(error, abs(gauges%gauge_y(1) - 8.3_wp) < 1.0e-10_wp, &
                    "gauge 1 y-coordinate mismatch")
         if (allocated(error)) exit checks

         call check(error, abs(gauges%gauge_x(2) - 5.0_wp) < 1.0e-10_wp, &
                    "gauge 2 x-coordinate mismatch")
         if (allocated(error)) exit checks

         call check(error, abs(gauges%gauge_y(2) - 3.0_wp) < 1.0e-10_wp, &
                    "gauge 2 y-coordinate mismatch")
         if (allocated(error)) exit checks

      end block checks
      continue
      call gauges_destroy(gauges)
      call state%destroy()
      call delete_file(gauge_filename)
      call delete_file(nc_filename)

   end subroutine test_gauge_file_parsing

   ! ================================================================
   ! Test 3: Bilinear interpolation — known linear field
   ! ================================================================

   subroutine test_bilinear_interpolation(error)
      type(error_type), allocatable, intent(out) :: error

      type(gauge_output_t) :: gauges
      type(config_t) :: cfg
      type(state_t) :: state
      character(len=64) :: names(1)
      real(wp) :: x_coords(1), y_coords(1)
      character(len=*), parameter :: gauge_filename = "test_gauge_interp.txt"
      character(len=*), parameter :: nc_filename = "gauges.nc"
      integer :: i, j, ng
      real(wp) :: expected_eta
      integer :: ncid, varid, status, dimid, ntime
      real(wp) :: eta_read(1)
      checks: block

         ! Grid: 20x10 cells, dx=dy=1.0, nghost=2
         ! Cell centres at x = (i - nghost - 0.5)*dx for physical cell i
         ! Physical cells: i from nghost+1 to nghost+nx_phys
         ! Cell 1 (physical) centre: x = 0.5*dx = 0.5
         ! Cell k (physical) centre: x = (k - 0.5)*dx
         call state%init(20, 10, 2, 1.0_wp, 1.0_wp, 0.45_wp)
         ng = state%grid%nghost

         ! Set h(i,j) = real(i - ng, wp) for a field that varies linearly in x
         ! b = 0 so eta = h + b = h
         state%barotropic%b(:, :) = 0.0_wp
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               state%barotropic%h(i, j) = real(i - ng, wp)
            end do
         end do

         ! Place gauge at x = 5.5 (between physical cell 5 and 6)
         ! Physical cell 5 has array index ng+5 = 7, centre at x = 4.5
         ! Physical cell 6 has array index ng+6 = 8, centre at x = 5.5
         ! h at cell 5: 5.0, h at cell 6: 6.0
         ! At x = 5.5, the bilinear interpolation of h should be:
         !   weight towards cell 6 = (5.5 - 4.5) / 1.0 = 1.0 => h = 6.0
         ! Actually let's use x = 5.0 (midpoint between cells 5 and 6):
         !   h = 0.5 * 5.0 + 0.5 * 6.0 = 5.5
         ! Use y = 3.5 (cell 3 centre at 2.5, cell 4 centre at 3.5 => weight = 1.0)
         ! h varies only in i, so y-interpolation is trivial
         names(1) = "test_gauge"
         x_coords(1) = 5.0_wp
         y_coords(1) = 3.5_wp

         call write_gauge_file(gauge_filename, names, x_coords, y_coords, 1)

         cfg%gauge_file = gauge_filename
         cfg%gauge_output_interval = 0.0_wp
         cfg%output_dir = "."

         call gauges_init(gauges, cfg, state, 0)

         call check(error, gauges%active, "gauges should be active")
         if (allocated(error)) exit checks

         ! Write one snapshot at time = 0
         call gauges_write(gauges, state, 0.0_wp)

         ! Close the NetCDF file so we can read it back
         call gauges_destroy(gauges)

         ! Read back the eta value from the NetCDF file
         status = nf90_open(nc_filename, nf90_nowrite, ncid)
         call check(error, status == nf90_noerr, "failed to open gauges.nc for reading")
         if (allocated(error)) exit checks

         status = nf90_inq_varid(ncid, "eta", varid)
         call check(error, status == nf90_noerr, "eta variable not found in gauges.nc")
         if (allocated(error)) then
            status = nf90_close(ncid)
            exit checks
         end if

         ! Read eta(time=1, station=1)
         status = nf90_get_var(ncid, varid, eta_read, start=[1, 1], count=[1, 1])
         call check(error, status == nf90_noerr, "failed to read eta from gauges.nc")
         if (allocated(error)) then
            status = nf90_close(ncid)
            exit checks
         end if
         status = nf90_close(ncid)

         ! Expected: bilinear interpolation of h at x=5.0
         ! Cell centres: cell 5 at x=4.5 (h=5), cell 6 at x=5.5 (h=6)
         ! wx = (5.0 - 4.5) / 1.0 = 0.5
         ! eta = h + b = (1-0.5)*5.0 + 0.5*6.0 + 0 = 5.5
         expected_eta = 5.5_wp

         ! Bilinear interpolation reproduces a LINEAR field exactly, so
         ! the tolerance here should be roundoff. The previous 0.1 bound
         ! (~2% relative error) would have hidden a broken interpolation
         ! that silently returned the nearest cell centre.
         call check(error, abs(eta_read(1) - expected_eta) < 1.0e-12_wp, &
                    "bilinear interpolation of a linear field must be exact")
         if (allocated(error)) exit checks

      end block checks
      continue
      call state%destroy()
      call delete_file(gauge_filename)
      call delete_file(nc_filename)

   end subroutine test_bilinear_interpolation

   ! ================================================================
   ! Test 4: Write multiple snapshots — verify time dimension grows
   ! ================================================================

   subroutine test_write_multiple_snapshots(error)
      type(error_type), allocatable, intent(out) :: error

      type(gauge_output_t) :: gauges
      type(config_t) :: cfg
      type(state_t) :: state
      character(len=64) :: names(1)
      real(wp) :: x_coords(1), y_coords(1)
      character(len=*), parameter :: gauge_filename = "test_gauge_multi.txt"
      character(len=*), parameter :: nc_filename = "gauges.nc"
      integer :: ncid, varid, dimid, status, ntime
      real(wp) :: time_vals(3), eta_vals(3)
      checks: block

         ! Set up a simple 20x10 grid
         call state%init(20, 10, 2, 1.0_wp, 1.0_wp, 0.45_wp)
         state%barotropic%h(:, :) = 2.0_wp
         state%barotropic%b(:, :) = 0.0_wp
         state%barotropic%hu(:, :) = 0.0_wp
         state%barotropic%hv(:, :) = 0.0_wp

         ! One gauge at the domain centre
         names(1) = "centre"
         x_coords(1) = 10.0_wp
         y_coords(1) = 5.0_wp
         call write_gauge_file(gauge_filename, names, x_coords, y_coords, 1)

         cfg%gauge_file = gauge_filename
         cfg%gauge_output_interval = 0.0_wp
         cfg%output_dir = "."

         call gauges_init(gauges, cfg, state, 0)

         call check(error, gauges%active, "gauges should be active")
         if (allocated(error)) exit checks

         ! Write 3 snapshots at different times with varying h
         state%barotropic%h(:, :) = 2.0_wp
         call gauges_write(gauges, state, 0.0_wp)

         state%barotropic%h(:, :) = 3.0_wp
         call gauges_write(gauges, state, 1.0_wp)

         state%barotropic%h(:, :) = 4.0_wp
         call gauges_write(gauges, state, 2.0_wp)

         ! Close and read back
         call gauges_destroy(gauges)

         status = nf90_open(nc_filename, nf90_nowrite, ncid)
         call check(error, status == nf90_noerr, "failed to open gauges.nc")
         if (allocated(error)) exit checks

         ! Check time dimension length
         status = nf90_inq_dimid(ncid, "time", dimid)
         call check(error, status == nf90_noerr, "time dimension not found")
         if (allocated(error)) then
            status = nf90_close(ncid)
            exit checks
         end if

         status = nf90_inquire_dimension(ncid, dimid, len=ntime)
         call check(error, ntime == 3, "expected 3 time steps in gauge output")
         if (allocated(error)) then
            status = nf90_close(ncid)
            exit checks
         end if

         ! Check time values
         status = nf90_inq_varid(ncid, "time", varid)
         status = nf90_get_var(ncid, varid, time_vals)
         call check(error, abs(time_vals(1) - 0.0_wp) < 1.0e-10_wp, "time(1) should be 0")
         if (allocated(error)) then
            status = nf90_close(ncid)
            exit checks
         end if

         call check(error, abs(time_vals(2) - 1.0_wp) < 1.0e-10_wp, "time(2) should be 1")
         if (allocated(error)) then
            status = nf90_close(ncid)
            exit checks
         end if

         call check(error, abs(time_vals(3) - 2.0_wp) < 1.0e-10_wp, "time(3) should be 2")
         if (allocated(error)) then
            status = nf90_close(ncid)
            exit checks
         end if

         ! Check eta values (h + b = h since b = 0)
         status = nf90_inq_varid(ncid, "eta", varid)
         status = nf90_get_var(ncid, varid, eta_vals, start=[1, 1], count=[1, 3])
         call check(error, abs(eta_vals(1) - 2.0_wp) < 1.0e-6_wp, "eta(t=0) should be 2.0")
         if (allocated(error)) then
            status = nf90_close(ncid)
            exit checks
         end if

         call check(error, abs(eta_vals(2) - 3.0_wp) < 1.0e-6_wp, "eta(t=1) should be 3.0")
         if (allocated(error)) then
            status = nf90_close(ncid)
            exit checks
         end if

         call check(error, abs(eta_vals(3) - 4.0_wp) < 1.0e-6_wp, "eta(t=2) should be 4.0")
         if (allocated(error)) then
            status = nf90_close(ncid)
            exit checks
         end if

         status = nf90_close(ncid)

      end block checks
      continue
      call state%destroy()
      call delete_file(gauge_filename)
      call delete_file(nc_filename)

   end subroutine test_write_multiple_snapshots

end module test_gauges
