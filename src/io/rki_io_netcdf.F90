!! Thin wrapper around netcdf-fortran (nf90_*) with error checking
module rki_io_netcdf
   !! Provides NetCDF read/write primitives and an error-checking helper.
   !! Used by bathymetry loader and output writer.
   use, intrinsic :: iso_fortran_env, only: real32
   use netcdf, only: nf90_noerr, nf90_strerror, nf90_create, nf90_open, &
                     nf90_close, nf90_def_dim, nf90_def_var, nf90_put_att, nf90_enddef, &
                     nf90_put_var, nf90_get_var, nf90_inq_dimid, nf90_inquire_dimension, &
                     nf90_inq_varid, nf90_get_att, nf90_def_var_deflate, &
                     nf90_clobber, nf90_netcdf4, nf90_nowrite, nf90_write, &
                     nf90_unlimited, nf90_double, nf90_float, nf90_int, nf90_global, &
                     nf90_redef, nf90_sync
   use rki_constants, only: wp
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: nc_check
   public :: nc_create_file
   public :: nc_open_read
   public :: nc_open_write
   public :: nc_close
   public :: nc_def_dim
   public :: nc_def_var_2d
   public :: nc_def_var_3d
   public :: nc_put_att
   public :: nc_put_att_real
   public :: nc_put_att_int
   public :: nc_put_att_global
   public :: nc_put_att_global_int
   public :: nc_get_att_int
   public :: nc_enddef
   public :: nc_put_var_2d
   public :: nc_put_var_3d_slice
   public :: nc_get_var_2d
   public :: nc_get_var_3d
   public :: nc_get_var_1d
   public :: nc_get_dim_len
   public :: nc_get_varid
   public :: rki_def_var_1d
   public :: rki_put_var_1d
   public :: rki_put_var_2d_slice
   public :: nc_def_var_4d
   public :: nc_put_var_4d_slice
   public :: nc_get_var_slab_3d
   public :: nc_get_att_text
   public :: nc_put_att_real_r4
   public :: nc_put_var_3d_slice_r4
   public :: nc_put_var_4d_slice_r4
   public :: NC_WP
   public :: NC_R4

#ifdef RAKALI_DOUBLE_PRECISION
   integer, parameter :: NC_WP = nf90_double
      !! NetCDF type matching working precision
#else
   integer, parameter :: NC_WP = nf90_float
      !! NetCDF type matching working precision
#endif

   integer, parameter :: NC_R4 = nf90_float
      !! Explicit 32-bit NetCDF element type.  Used ONLY by the ocean
      !! diagnostic stream when `&ocean_diag_nml output_precision="single"`
      !! asks for a halved-width output file.  Restarts, gauges, coastal
      !! output and console/conservation totals stay on `NC_WP` — a lossy
      !! restart must not be requestable.  The 64-bit counterpart is
      !! `NC_WP` itself — deliberately NOT given an `NC_R8` alias, so
      !! "double" has exactly one spelling and cannot drift away from the
      !! build's working precision.

contains

   subroutine nc_check(status, context)
      !! Check NetCDF return status and log error if it failed
      integer, intent(in) :: status
         !! Return code from nf90_* call
      character(len=*), intent(in) :: context
         !! Description of what was being attempted

      if (status /= nf90_noerr) then
         call logger%error("NetCDF error in "//trim(context)//": "// &
                           trim(nf90_strerror(status)))
         error stop "NetCDF operation failed"
      end if

   end subroutine nc_check

   subroutine nc_create_file(filename, ncid)
      !! Create a new NetCDF-4 file (overwrites if exists)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ncid

      call nc_check(nf90_create(trim(filename), &
                                ior(nf90_clobber, nf90_netcdf4), ncid), &
                    "creating "//trim(filename))

   end subroutine nc_create_file

   subroutine nc_open_read(filename, ncid)
      !! Open an existing NetCDF file for reading
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ncid

      call nc_check(nf90_open(trim(filename), nf90_nowrite, ncid), &
                    "opening "//trim(filename))

   end subroutine nc_open_read

   subroutine nc_open_write(filename, ncid)
      !! Open an existing NetCDF file for writing
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ncid

      call nc_check(nf90_open(trim(filename), nf90_write, ncid), &
                    "opening for write "//trim(filename))

   end subroutine nc_open_write

   subroutine nc_close(ncid)
      !! Close a NetCDF file
      integer, intent(in) :: ncid

      call nc_check(nf90_close(ncid), "closing file")

   end subroutine nc_close

   subroutine nc_def_dim(ncid, name, length, dimid)
      !! Define a dimension (use nf90_unlimited for unlimited)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: length
      integer, intent(out) :: dimid

      call nc_check(nf90_def_dim(ncid, name, length, dimid), &
                    "defining dimension "//trim(name))

   end subroutine nc_def_dim

   subroutine rki_def_var_1d(ncid, name, dimid, varid, deflate_level)
      !! Define a 1D variable with working precision type
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: dimid
      integer, intent(out) :: varid
      integer, intent(in), optional :: deflate_level
         !! Compression level (0=none, 1-9=deflate). Enables shuffle when > 0.

      call nc_check(nf90_def_var(ncid, name, NC_WP, [dimid], varid), &
                    "defining variable "//trim(name))

      if (present(deflate_level)) then
         if (deflate_level > 0) then
            call nc_check(nf90_def_var_deflate(ncid, varid, &
                                               shuffle=1, deflate=1, deflate_level=deflate_level), &
                          "enabling compression for "//trim(name))
         end if
      end if

   end subroutine rki_def_var_1d

   subroutine rki_put_var_1d(ncid, varid, data)
      !! Write a full 1D array
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(in) :: data(:)

      call nc_check(nf90_put_var(ncid, varid, data), "writing 1D variable")

   end subroutine rki_put_var_1d

   subroutine rki_put_var_2d_slice(ncid, varid, data, time_index)
      !! Write a 1D slice into a 2D variable at a given time index
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(in) :: data(:)
      integer, intent(in) :: time_index

      call nc_check(nf90_put_var(ncid, varid, data, &
                                 start=[1, time_index], &
                                 count=[size(data), 1]), &
                    "writing 2D slice")

   end subroutine rki_put_var_2d_slice

   subroutine nc_def_var_2d(ncid, name, dimids, varid, deflate_level)
      !! Define a 2D variable with working precision type
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: dimids(2)
      integer, intent(out) :: varid
      integer, intent(in), optional :: deflate_level
         !! Compression level (0=none, 1-9=deflate). Enables shuffle when > 0.

      call nc_check(nf90_def_var(ncid, name, NC_WP, dimids, varid), &
                    "defining variable "//trim(name))

      if (present(deflate_level)) then
         if (deflate_level > 0) then
            call nc_check(nf90_def_var_deflate(ncid, varid, &
                                               shuffle=1, deflate=1, deflate_level=deflate_level), &
                          "enabling compression for "//trim(name))
         end if
      end if

   end subroutine nc_def_var_2d

   subroutine nc_def_var_3d(ncid, name, dimids, varid, deflate_level, xtype)
      !! Define a 3D variable (x, y, time) with working precision type
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: dimids(3)
      integer, intent(out) :: varid
      integer, intent(in), optional :: deflate_level
         !! Compression level (0=none, 1-9=deflate). Enables shuffle when > 0.
      integer, intent(in), optional :: xtype
         !! NetCDF element type.  Absent => `NC_WP` (working precision), so
         !! every legacy call site is unchanged.  The ocean diagnostic
         !! stream passes `NC_R4` when single-precision output is requested;
         !! the caller is then responsible for handing `nf90_put_var` real32
         !! data (see `nc_put_var_3d_slice_r4`).
      integer :: xt

      xt = NC_WP
      if (present(xtype)) xt = xtype

      call nc_check(nf90_def_var(ncid, name, xt, dimids, varid), &
                    "defining variable "//trim(name))

      if (present(deflate_level)) then
         if (deflate_level > 0) then
            call nc_check(nf90_def_var_deflate(ncid, varid, &
                                               shuffle=1, deflate=1, deflate_level=deflate_level), &
                          "enabling compression for "//trim(name))
         end if
      end if

   end subroutine nc_def_var_3d

   subroutine nc_put_att(ncid, varid, name, value)
      !! Write a character attribute to a variable
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: value

      call nc_check(nf90_put_att(ncid, varid, name, value), &
                    "writing attribute "//trim(name))

   end subroutine nc_put_att

   subroutine nc_put_att_real(ncid, varid, name, value)
      !! Write a real-valued attribute to a variable
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: value

      call nc_check(nf90_put_att(ncid, varid, name, value), &
                    "writing attribute "//trim(name))

   end subroutine nc_put_att_real

   subroutine nc_put_att_real_r4(ncid, varid, name, value)
      !! Write a 32-bit real attribute to a variable.  NetCDF rejects a
      !! `_FillValue` whose type differs from the variable's, so a variable
      !! defined with `NC_R4` must have its `_FillValue` / `missing_value`
      !! written through here, not through `nc_put_att_real` (real64).
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      character(len=*), intent(in) :: name
      real(real32), intent(in) :: value

      call nc_check(nf90_put_att(ncid, varid, name, value), &
                    "writing real32 attribute "//trim(name))

   end subroutine nc_put_att_real_r4

   subroutine nc_put_att_global(ncid, name, value)
      !! Write a global character attribute
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: value

      call nc_check(nf90_put_att(ncid, nf90_global, name, value), &
                    "writing global attribute "//trim(name))

   end subroutine nc_put_att_global

   subroutine nc_put_att_int(ncid, varid, name, value)
      !! Write an integer attribute to a variable
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      character(len=*), intent(in) :: name
      integer, intent(in) :: value

      call nc_check(nf90_put_att(ncid, varid, name, value), &
                    "writing attribute "//trim(name))

   end subroutine nc_put_att_int

   subroutine nc_put_att_global_int(ncid, name, value)
      !! Write a global integer attribute
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: value

      call nc_check(nf90_put_att(ncid, nf90_global, name, value), &
                    "writing global attribute "//trim(name))

   end subroutine nc_put_att_global_int

   subroutine nc_get_att_int(ncid, name, value)
      !! Read a global integer attribute
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: value

      call nc_check(nf90_get_att(ncid, nf90_global, name, value), &
                    "reading global attribute "//trim(name))

   end subroutine nc_get_att_int

   subroutine nc_enddef(ncid)
      !! End define mode, switch to data mode
      integer, intent(in) :: ncid

      call nc_check(nf90_enddef(ncid), "ending define mode")

   end subroutine nc_enddef

   subroutine nc_put_var_2d(ncid, varid, data)
      !! Write a full 2D array
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(in) :: data(:, :)

      call nc_check(nf90_put_var(ncid, varid, data), "writing 2D variable")

   end subroutine nc_put_var_2d

   subroutine nc_put_var_3d_slice(ncid, varid, data, time_index)
      !! Write a 2D slice into a 3D variable at a given time index
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(in) :: data(:, :)
      integer, intent(in) :: time_index

      call nc_check(nf90_put_var(ncid, varid, data, &
                                 start=[1, 1, time_index], &
                                 count=[size(data, 1), size(data, 2), 1]), &
                    "writing 3D slice")

   end subroutine nc_put_var_3d_slice

   subroutine nc_put_var_3d_slice_r4(ncid, varid, data, time_index)
      !! Write a 32-bit 2D slice into a 3D variable at a given time index.
      !! The real32 twin of `nc_put_var_3d_slice`: handing real64 data to
      !! an `NF90_FLOAT` variable is legal (netcdf converts) but writes the
      !! same bytes at twice the host->library traffic, so the diag writer
      !! stages the conversion itself and calls this.
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(real32), intent(in) :: data(:, :)
      integer, intent(in) :: time_index

      call nc_check(nf90_put_var(ncid, varid, data, &
                                 start=[1, 1, time_index], &
                                 count=[size(data, 1), size(data, 2), 1]), &
                    "writing 3D slice (real32)")

   end subroutine nc_put_var_3d_slice_r4

   subroutine nc_put_var_4d_slice_r4(ncid, varid, data, time_index)
      !! Write a 32-bit 3D slice into a 4D variable at a given time index.
      !! The real32 twin of `nc_put_var_4d_slice`.
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(real32), intent(in) :: data(:, :, :)
      integer, intent(in) :: time_index

      call nc_check(nf90_put_var(ncid, varid, data, &
                                 start=[1, 1, 1, time_index], &
                                 count=[size(data, 1), size(data, 2), size(data, 3), 1]), &
                    "writing 4D slice (real32)")

   end subroutine nc_put_var_4d_slice_r4

   subroutine nc_get_var_2d(ncid, varid, data)
      !! Read a full 2D array
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(inout) :: data(:, :)

      call nc_check(nf90_get_var(ncid, varid, data), "reading 2D variable")

   end subroutine nc_get_var_2d

   subroutine nc_get_var_3d(ncid, varid, data)
      !! Read a full 3D array.  Mirrors `nc_get_var_2d`; the target's
      !! shape must match the file variable's Fortran storage order
      !! (the caller handles any C/Fortran dimension reversal, as
      !! `rki_bathymetry` does for 2D).
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(inout) :: data(:, :, :)

      call nc_check(nf90_get_var(ncid, varid, data), "reading 3D variable")

   end subroutine nc_get_var_3d

   subroutine nc_get_var_1d(ncid, varid, data)
      !! Read a full 1D array (the read sibling of `rki_put_var_1d`).
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(inout) :: data(:)

      call nc_check(nf90_get_var(ncid, varid, data), "reading 1D variable")

   end subroutine nc_get_var_1d

   subroutine nc_get_dim_len(ncid, name, length)
      !! Get the length of a named dimension
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: length

      integer :: dimid

      call nc_check(nf90_inq_dimid(ncid, name, dimid), &
                    "finding dimension "//trim(name))
      call nc_check(nf90_inquire_dimension(ncid, dimid, len=length), &
                    "querying dimension "//trim(name))

   end subroutine nc_get_dim_len

   subroutine nc_def_var_4d(ncid, name, dimids, varid, deflate_level, xtype)
      !! Define a 4D variable (x, y, z, time) with working precision type
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: dimids(4)
      integer, intent(out) :: varid
      integer, intent(in), optional :: deflate_level
         !! Compression level (0=none, 1-9=deflate). Enables shuffle when > 0.
      integer, intent(in), optional :: xtype
         !! NetCDF element type.  Absent => `NC_WP` (working precision).
         !! See `nc_def_var_3d` — same contract, 4D shape.
      integer :: xt

      xt = NC_WP
      if (present(xtype)) xt = xtype

      call nc_check(nf90_def_var(ncid, name, xt, dimids, varid), &
                    "defining variable "//trim(name))

      if (present(deflate_level)) then
         if (deflate_level > 0) then
            call nc_check(nf90_def_var_deflate(ncid, varid, &
                                               shuffle=1, deflate=1, deflate_level=deflate_level), &
                          "enabling compression for "//trim(name))
         end if
      end if

   end subroutine nc_def_var_4d

   subroutine nc_put_var_4d_slice(ncid, varid, data, time_index)
      !! Write a 3D slice into a 4D variable at a given time index
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(in) :: data(:, :, :)
      integer, intent(in) :: time_index

      call nc_check(nf90_put_var(ncid, varid, data, &
                                 start=[1, 1, 1, time_index], &
                                 count=[size(data, 1), size(data, 2), size(data, 3), 1]), &
                    "writing 4D slice")

   end subroutine nc_put_var_4d_slice

   subroutine nc_get_var_slab_3d(ncid, varid, start, count, data)
      !! Strided slab read: `start`/`count` are the NetCDF start/count
      !! vectors in FORTRAN dimension order, with length matching the
      !! target variable's rank (3 for a 2D-field-plus-time variable, 4
      !! for a 3D-field-plus-time variable).  `data` is always a plain
      !! 3D Fortran array; for the 2D-plus-time case pass a target with
      !! a trailing singleton extent (`data(:,:,1)`) — the netcdf-fortran
      !! generic interface dispatches on the array RANK of `data`, not
      !! the length of `start`/`count`, so one routine covers both
      !! shapes. This is the first strided (non-full-variable) read in
      !! the tree — mirrors `nc_put_var_4d_slice` for the start/count
      !! shape, on the read side.
      integer, intent(in) :: ncid, varid
      integer, intent(in) :: start(:), count(:)
      real(wp), intent(inout) :: data(:, :, :)

      call nc_check(nf90_get_var(ncid, varid, data, start=start, count=count), &
                    "reading NetCDF slab")

   end subroutine nc_get_var_slab_3d

   subroutine nc_get_att_text(ncid, varid, name, value, ok)
      !! Read a character (text) attribute from a variable — used for the
      !! CF `units` attribute on a time axis.  Non-fail-loud: `ok` is
      !! `.false.` (and `value` blank) when the attribute is absent or the
      !! wrong type, so callers can fall back to a default rather than
      !! aborting on a missing/optional attribute.
      integer, intent(in) :: ncid, varid
      character(len=*), intent(in) :: name
      character(len=:), allocatable, intent(out) :: value
         !! Deferred-length (not assumed-length `intent(out)`, which risks
         !! silent truncation — fortitude C072): filled from a fixed local
         !! buffer sized for a CF `units` string.
      logical, intent(out) :: ok

      integer :: status
      character(len=1024) :: buf

      buf = ""
      status = nf90_get_att(ncid, varid, name, buf)
      ok = (status == nf90_noerr)
      if (ok) then
         value = trim(buf)
      else
         value = ""
      end if

   end subroutine nc_get_att_text

   subroutine nc_get_varid(ncid, name, varid)
      !! Get the variable ID for a named variable
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: varid

      call nc_check(nf90_inq_varid(ncid, name, varid), &
                    "finding variable "//trim(name))

   end subroutine nc_get_varid

end module rki_io_netcdf
