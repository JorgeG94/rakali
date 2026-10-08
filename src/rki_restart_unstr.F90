!! Restart file I/O for the unstructured solver
module rki_restart_unstr
   !! Writes and reads checkpoint files for the unstructured solver.
   !! Stores `(stage_cv, xmom_cv, ymom_cv, bed_cv)` for owned cells,
   !! plus simulation time and step count as global attributes.
   !!
   !! File format: NetCDF-4 with one dimension `face = ncells_owned`.
   !! Cell global ids stored if available.
   use rki_constants, only: wp
   use rki_state_unstr, only: unstr_state_t
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_io_netcdf, only: nc_check, nc_create_file, nc_open_read, &
                            nc_close, nc_def_dim, rki_def_var_1d, &
                            nc_put_att, nc_put_att_real, &
                            nc_put_att_global, nc_put_att_global_int, &
                            nc_enddef, rki_put_var_1d, nc_get_dim_len, &
                            nc_get_varid, NC_WP
   use netcdf, only: nf90_int64, nf90_def_var, nf90_put_var, &
                     nf90_get_var, nf90_global, nf90_put_att, &
                     nf90_get_att, nf90_noerr, nf90_inquire_attribute
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: restart_write_unstr
   public :: restart_read_unstr

contains

   subroutine restart_write_unstr(filename, state, mesh, time)
      !! Write a checkpoint file. Caller must ensure state is on host.
      character(len=*), intent(in) :: filename
      type(unstr_state_t), intent(in) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: time

      integer :: ncid, face_dimid
      integer :: stage_varid, xmom_varid, ymom_varid, bed_varid
      integer :: cgid_varid
      integer :: n

      n = mesh%ncells_owned

      call nc_create_file(trim(filename), ncid)

      ! Global attributes
      call nc_put_att_global(ncid, "Conventions", "CF-1.11")
      call nc_put_att_global(ncid, "title", "Rakali unstructured restart")
      call nc_put_att_global(ncid, "source", "Rakali")
      call nc_put_att_global_int(ncid, "ncells_owned", n)
      call nc_put_att_global_int(ncid, "step_count", state%step_count)
      call nc_check(nf90_put_att(ncid, nf90_global, "time", time), &
                    "writing time attribute")

      ! Dimension
      call nc_def_dim(ncid, "face", n, face_dimid)

      ! Solution fields
      call rki_def_var_1d(ncid, "stage", face_dimid, stage_varid)
      call nc_put_att(ncid, stage_varid, "long_name", "Free-surface elevation")
      call nc_put_att(ncid, stage_varid, "units", "m")

      call rki_def_var_1d(ncid, "xmom", face_dimid, xmom_varid)
      call nc_put_att(ncid, xmom_varid, "long_name", "x-momentum (h*u)")
      call nc_put_att(ncid, xmom_varid, "units", "m2/s")

      call rki_def_var_1d(ncid, "ymom", face_dimid, ymom_varid)
      call nc_put_att(ncid, ymom_varid, "long_name", "y-momentum (h*v)")
      call nc_put_att(ncid, ymom_varid, "units", "m2/s")

      call rki_def_var_1d(ncid, "bed", face_dimid, bed_varid)
      call nc_put_att(ncid, bed_varid, "long_name", "Bed elevation")
      call nc_put_att(ncid, bed_varid, "units", "m")

      if (allocated(mesh%cell_global_id)) then
         call nc_check(nf90_def_var(ncid, "cell_global_id", nf90_int64, &
                                    [face_dimid], cgid_varid), &
                       "defining cell_global_id")
      end if

      call nc_enddef(ncid)

      ! Write data
      call rki_put_var_1d(ncid, stage_varid, state%barotropic%stage_cv(1:n))
      call rki_put_var_1d(ncid, xmom_varid, state%barotropic%xmom_cv(1:n))
      call rki_put_var_1d(ncid, ymom_varid, state%barotropic%ymom_cv(1:n))
      call rki_put_var_1d(ncid, bed_varid, state%barotropic%bed_cv(1:n))

      if (allocated(mesh%cell_global_id)) then
         call nc_check(nf90_put_var(ncid, cgid_varid, &
                                    mesh%cell_global_id(1:n)), &
                       "writing cell_global_id")
      end if

      call nc_close(ncid)

      call logger%info("Restart written: "//trim(filename)// &
                       " (t = "//to_string(time)//" s, "// &
                       to_string(n)//" cells)")

   end subroutine restart_write_unstr

   subroutine restart_read_unstr(filename, state, time, ierr)
      !! Read a checkpoint and populate state arrays.
      !!
      !! Must be called before finalize (arrays on host). Validates
      !! that the file's cell count matches `state%ncells_owned`.
      character(len=*), intent(in) :: filename
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(out) :: time
      integer, intent(out) :: ierr
         !! 0 = success, 1 = cell count mismatch

      integer :: ncid, varid, file_ncells, n, sc, att_stat
      real(wp), allocatable :: tmp(:)

      n = state%ncells_owned

      call nc_open_read(trim(filename), ncid)

      ! Validate cell count
      call nc_get_dim_len(ncid, "face", file_ncells)
      if (file_ncells /= n) then
         call logger%error("Restart cell count mismatch: file has "// &
                           to_string(file_ncells)//" but state expects "// &
                           to_string(n))
         call nc_close(ncid)
         ierr = 1
         return
      end if

      ! Read time from global attribute
      call nc_check(nf90_get_att(ncid, nf90_global, "time", time), &
                    "reading time attribute")

      ! Read step_count if present
      att_stat = nf90_inquire_attribute(ncid, nf90_global, "step_count")
      if (att_stat == nf90_noerr) then
         call nc_check(nf90_get_att(ncid, nf90_global, "step_count", sc), &
                       "reading step_count")
         state%step_count = sc
      end if

      ! Read fields
      allocate (tmp(n))

      call nc_get_varid(ncid, "stage", varid)
      call nc_check(nf90_get_var(ncid, varid, tmp), "reading stage")
      state%barotropic%stage_cv(1:n) = tmp

      call nc_get_varid(ncid, "xmom", varid)
      call nc_check(nf90_get_var(ncid, varid, tmp), "reading xmom")
      state%barotropic%xmom_cv(1:n) = tmp

      call nc_get_varid(ncid, "ymom", varid)
      call nc_check(nf90_get_var(ncid, varid, tmp), "reading ymom")
      state%barotropic%ymom_cv(1:n) = tmp

      call nc_get_varid(ncid, "bed", varid)
      call nc_check(nf90_get_var(ncid, varid, tmp), "reading bed")
      state%barotropic%bed_cv(1:n) = tmp

      deallocate (tmp)
      call nc_close(ncid)

      state%time = time
      ierr = 0

      call logger%info("Restart loaded: "//trim(filename)// &
                       " (t = "//to_string(time)//" s, "// &
                       to_string(n)//" cells)")

   end subroutine restart_read_unstr

end module rki_restart_unstr
