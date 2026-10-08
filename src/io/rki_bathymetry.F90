!! Bathymetry loading from NetCDF files
module rki_bathymetry
   !! Reads bathymetry (bottom elevation) from a NetCDF file and populates
   !! the state%barotropic%b array. Expects a 2D variable on a structured grid.
   !! If no file is configured, falls back to flat bottom.
   use rki_constants, only: wp
   use rki_grid, only: hgrid_t
   use rki_io_netcdf, only: nc_check, nc_open_read, nc_close, &
                            nc_get_dim_len, nc_get_varid, nc_get_var_2d
   use netcdf, only: nf90_inquire_variable, nf90_inquire_dimension
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: load_bathymetry_into_array
   public :: fill_bathymetry_ghosts_array

contains

   subroutine load_bathymetry_into_array(filename, b, grid)
      !! Load bathymetry from a NetCDF file into a target 2D array
      !! shaped `(grid%nx_total, grid%ny_total)`.
      !!
      !! The file must contain:
      !!   - Dimensions: x (nx_phys), y (ny_phys)
      !!   - Variable: b(x, y) or elevation(x, y) or depth(x, y)
      !!
      !! Handles NetCDF C/Fortran dimension reversal: files written by
      !! Python/C store b(x,y) in C order, which Fortran reads as b(y,x).
      !! If the direct dimensions don't match but the transposed ones do,
      !! we read transposed and copy correctly.
      !!
      !! Bathymetry lands in the interior; ghost cells are filled by
      !! constant extrapolation from the nearest interior cell via
      !! `fill_bathymetry_ghosts_array`.
      character(len=*), intent(in) :: filename
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid

      integer :: ncid, varid, file_nx, file_ny
      integer :: ng, i, j
      integer :: var_dimids(2)
      integer :: var_ndims, dim1_len, dim2_len
      integer :: status
      character(len=64) :: dim1_name, dim2_name
      logical :: needs_transpose
      real(wp), allocatable :: b_interior(:, :)

      call logger%info("Loading bathymetry from: "//trim(filename))

      call nc_open_read(filename, ncid)

      ! Read grid dimensions (by name — these are the logical sizes)
      call nc_get_dim_len(ncid, "x", file_nx)
      call nc_get_dim_len(ncid, "y", file_ny)

      ! Try variable names in order: b, elevation, depth
      call try_get_bathymetry_var(ncid, varid)

      ! Query the variable's actual dimension order in the file.
      ! Files written by Python/C store b(x,y) in C-order. Fortran's
      ! NetCDF library reverses this, so the variable appears as b(y,x)
      ! in Fortran. We detect this by checking the first dimension name.
      call nc_check(nf90_inquire_variable(ncid, varid, ndims=var_ndims, &
                                          dimids=var_dimids), "querying bathymetry variable")
      call nc_check(nf90_inquire_dimension(ncid, var_dimids(1), &
                                           name=dim1_name, len=dim1_len), "querying dim 1")
      call nc_check(nf90_inquire_dimension(ncid, var_dimids(2), &
                                           name=dim2_name, len=dim2_len), "querying dim 2")

      call logger%info("Bathymetry variable: "//trim(dim1_name)//"="// &
                       to_string(dim1_len)//" x "//trim(dim2_name)//"="// &
                       to_string(dim2_len)//" (Fortran order)")

      ! In Fortran, dim1 is the first array index (varies fastest).
      ! If dim1 is "x" and matches nx, no transpose needed.
      ! If dim1 is "y" (C-order file reversed), we need to transpose.
      needs_transpose = (trim(dim1_name) == "y")

      ! Validate dimensions
      if (needs_transpose) then
         if (dim1_len /= grid%ny_phys .or. dim2_len /= grid%nx_phys) then
            call logger%error("Bathymetry grid mismatch: file has "// &
                              to_string(dim2_len)//" x "//to_string(dim1_len)// &
                              " but simulation expects "// &
                              to_string(grid%nx_phys)//" x "// &
                              to_string(grid%ny_phys))
            call nc_close(ncid)
            error stop "Bathymetry grid mismatch"
         end if
      else
         if (dim1_len /= grid%nx_phys .or. dim2_len /= grid%ny_phys) then
            call logger%error("Bathymetry grid mismatch: file has "// &
                              to_string(dim1_len)//" x "//to_string(dim2_len)// &
                              " but simulation expects "// &
                              to_string(grid%nx_phys)//" x "// &
                              to_string(grid%ny_phys))
            call nc_close(ncid)
            error stop "Bathymetry grid mismatch"
         end if
      end if

      ng = grid%nghost

      ! Read into temporary array matching the Fortran storage order
      allocate (b_interior(dim1_len, dim2_len))
      call nc_get_var_2d(ncid, varid, b_interior)
      call nc_close(ncid)

      if (needs_transpose) then
         ! b_interior is (y, x) in Fortran — transpose to (x, y) for b
         do j = 1, grid%ny_phys
            do i = 1, grid%nx_phys
               b(ng + i, ng + j) = b_interior(j, i)
            end do
         end do
      else
         b(ng + 1:ng + grid%nx_phys, ng + 1:ng + grid%ny_phys) = b_interior
      end if

      deallocate (b_interior)

      ! Fill ghost cells by constant extrapolation from nearest interior cell.
      ! Without this, ghost cells retain b=0 which creates artificial cliffs
      ! against real bathymetry (e.g. b=-50m interior vs b=0 ghost), generating
      ! extreme velocities and tiny CFL timesteps.
      call fill_bathymetry_ghosts_array(b, grid)

      call logger%info("Bathymetry loaded: min = "// &
                       to_string(minval(b(ng + 1:ng + grid%nx_phys, &
                                          ng + 1:ng + grid%ny_phys)))// &
                       " m, max = "// &
                       to_string(maxval(b(ng + 1:ng + grid%nx_phys, &
                                          ng + 1:ng + grid%ny_phys)))// &
                       " m")

   end subroutine load_bathymetry_into_array

   subroutine try_get_bathymetry_var(ncid, varid)
      !! Try to find the bathymetry variable by common names
      use netcdf, only: nf90_noerr, nf90_inq_varid
      integer, intent(in) :: ncid
      integer, intent(out) :: varid

      integer :: status

      ! Try "b" first (our convention)
      status = nf90_inq_varid(ncid, "b", varid)
      if (status == nf90_noerr) return

      ! Try "elevation"
      status = nf90_inq_varid(ncid, "elevation", varid)
      if (status == nf90_noerr) return

      ! Try "depth"
      status = nf90_inq_varid(ncid, "depth", varid)
      if (status == nf90_noerr) return

      ! None found
      call logger%error("No bathymetry variable found (tried: b, elevation, depth)")

   end subroutine try_get_bathymetry_var

   subroutine fill_bathymetry_ghosts_array(b, grid)
      !! Fill ghost cell bathymetry by constant extrapolation from the
      !! nearest interior cell.  This ensures boundary flux computations
      !! see a consistent bottom elevation across the ghost-interior
      !! interface.  Works on a target 2D array directly so both the
      !! coastal and ocean barotropic slots can reuse it.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid

      integer :: ng, i, j, nx, ny

      ng = grid%nghost
      nx = grid%nx_phys
      ny = grid%ny_phys

      ! West and east ghost columns
      do j = 1, grid%ny_total
         do i = 1, ng
            b(i, j) = b(ng + 1, j)                  ! west
            b(ng + nx + i, j) = b(ng + nx, j)       ! east
         end do
      end do

      ! South and north ghost rows (corners already filled above)
      do j = 1, ng
         do i = 1, grid%nx_total
            b(i, j) = b(i, ng + 1)                  ! south
            b(i, ng + ny + j) = b(i, ng + ny)       ! north
         end do
      end do

   end subroutine fill_bathymetry_ghosts_array

end module rki_bathymetry
