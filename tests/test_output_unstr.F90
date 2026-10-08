!! Unit tests for the unstructured UGRID output writer
module test_output_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t, &
                             mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy
   use rki_state_unstr, only: unstr_state_t, &
                              state_unstr_init, state_unstr_destroy
   use rki_output_unstr, only: output_unstr_state_t, &
                               output_unstr_init, &
                               output_unstr_init_from_arrays, &
                               output_unstr_write_snapshot, &
                               output_unstr_write_snapshot_from_arrays, &
                               output_unstr_close
   use netcdf, only: nf90_open, nf90_close, nf90_nowrite, nf90_noerr, &
                     nf90_inq_varid, nf90_inq_dimid, nf90_inquire_dimension
   implicit none
   private

   public :: collect_output_unstr_tests

contains

   subroutine collect_output_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("barotropic_round_trip", test_barotropic_round_trip), &
                  new_unittest("multilayer_round_trip", test_multilayer_round_trip), &
                  new_unittest("init_from_arrays", test_init_from_arrays), &
                  new_unittest("deflate_levels", test_deflate_levels) &
                  ]
   end subroutine collect_output_unstr_tests

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

   subroutine build_test_mesh(mesh)
      !! Same 4-triangle square mesh used across the unstructured test suite.
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp) :: vx(5), vy(5)
      integer :: tv(3, 4)

      vx = [0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.5_wp]
      vy = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.5_wp]
      tv(:, 1) = [1, 2, 5]
      tv(:, 2) = [2, 3, 5]
      tv(:, 3) = [3, 4, 5]
      tv(:, 4) = [4, 1, 5]

      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_test_mesh

   subroutine test_barotropic_round_trip(error)
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      type(output_unstr_state_t) :: ostate
      character(len=*), parameter :: fname = "test_unstr_baro.nc"
      integer :: ncid, dimid, varid, status, ntime, k
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh)

         do k = 1, mesh%ncells_owned
            state%barotropic%stage_cv(k) = 1.0_wp + 0.1_wp*real(k, wp)
            state%barotropic%xmom_cv(k) = 0.5_wp*real(k, wp)
            state%barotropic%ymom_cv(k) = -0.25_wp*real(k, wp)
            state%barotropic%bed_cv(k) = -10.0_wp
         end do

         call output_unstr_init(ostate, fname, state, mesh)
         call output_unstr_write_snapshot(ostate, state, 0.0_wp)
         call output_unstr_write_snapshot(ostate, state, 1.5_wp)
         call output_unstr_close(ostate)

         status = nf90_open(fname, nf90_nowrite, ncid)
         call check(error, status == nf90_noerr, "open ugrid output")
         if (allocated(error)) exit checks

         status = nf90_inq_dimid(ncid, "time", dimid)
         call check(error, status == nf90_noerr, "find time dim")
         if (allocated(error)) then
            status = nf90_close(ncid); exit checks
         end if

         status = nf90_inquire_dimension(ncid, dimid, len=ntime)
         call check(error, ntime == 2, "expected 2 time steps")
         if (allocated(error)) then
            status = nf90_close(ncid); exit checks
         end if

         status = nf90_inq_varid(ncid, "stage", varid)
         call check(error, status == nf90_noerr, "stage variable should be present")
         if (allocated(error)) then
            status = nf90_close(ncid); exit checks
         end if
         status = nf90_close(ncid)

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
      call delete_file(fname)
   end subroutine test_barotropic_round_trip

   subroutine test_multilayer_round_trip(error)
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      type(output_unstr_state_t) :: ostate
      character(len=*), parameter :: fname = "test_unstr_ml.nc"
      integer, parameter :: NZ = 3
      integer :: ncid, varid, status, k, kk
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)

         do k = 1, mesh%ncells_owned
            state%barotropic%stage_cv(k) = 1.0_wp
            state%barotropic%bed_cv(k) = -10.0_wp
            state%barotropic%xmom_cv(k) = 0.0_wp
            state%barotropic%ymom_cv(k) = 0.0_wp
            do kk = 1, NZ
               state%multilayer%h_layer(kk, k) = 11.0_wp/real(NZ, wp)
               state%multilayer%hu_layer(kk, k) = 0.0_wp
               state%multilayer%hv_layer(kk, k) = 0.0_wp
               state%multilayer%rho_layer(kk, k) = 1025.0_wp
            end do
            do kk = 1, NZ + 1
               state%multilayer%w_interface(kk, k) = 0.0_wp
            end do
         end do

         call output_unstr_init(ostate, fname, state, mesh)
         call output_unstr_write_snapshot(ostate, state, 0.0_wp)
         call output_unstr_close(ostate)

         status = nf90_open(fname, nf90_nowrite, ncid)
         call check(error, status == nf90_noerr, "open ml ugrid output")
         if (allocated(error)) exit checks

         status = nf90_inq_varid(ncid, "h_layer", varid)
         call check(error, status == nf90_noerr, "h_layer variable should be present")
         if (allocated(error)) then
            status = nf90_close(ncid); exit checks
         end if
         status = nf90_close(ncid)

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
      call delete_file(fname)
   end subroutine test_multilayer_round_trip

   subroutine test_init_from_arrays(error)
      type(error_type), allocatable, intent(out) :: error
      type(output_unstr_state_t) :: ostate
      character(len=*), parameter :: fname = "test_unstr_arrays.nc"
      integer, parameter :: NV = 4, NC = 2
      real(wp) :: vx(NV), vy(NV)
      integer :: tv(3, NC)
      real(wp) :: cx(NC), cy(NC), area(NC), bed(NC)
      real(wp) :: stage(NC), xmom(NC), ymom(NC)
      integer :: ncid, dimid, status, ntime
      checks: block

         vx = [0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp]
         vy = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp]
         tv(:, 1) = [1, 2, 3]
         tv(:, 2) = [1, 3, 4]
         cx = [0.66_wp, 0.33_wp]
         cy = [0.33_wp, 0.66_wp]
         area = [0.5_wp, 0.5_wp]
         bed = [-5.0_wp, -5.0_wp]
         stage = [1.0_wp, 1.0_wp]
         xmom = [0.0_wp, 0.0_wp]
         ymom = [0.0_wp, 0.0_wp]

         call output_unstr_init_from_arrays(ostate, fname, NV, NC, vx, vy, tv, &
                                            cx, cy, area, bed)
         call output_unstr_write_snapshot_from_arrays(ostate, stage, xmom, ymom, 0.0_wp)
         call output_unstr_write_snapshot_from_arrays(ostate, stage, xmom, ymom, 1.0_wp)
         call output_unstr_close(ostate)

         status = nf90_open(fname, nf90_nowrite, ncid)
         call check(error, status == nf90_noerr, "open from-arrays output")
         if (allocated(error)) exit checks

         status = nf90_inq_dimid(ncid, "time", dimid)
         status = nf90_inquire_dimension(ncid, dimid, len=ntime)
         call check(error, ntime == 2, "expected 2 time steps from arrays path")
         if (allocated(error)) then
            status = nf90_close(ncid); exit checks
         end if
         status = nf90_close(ncid)

      end block checks
      continue
      call delete_file(fname)
   end subroutine test_init_from_arrays

   subroutine test_deflate_levels(error)
      !! Init with two different deflate levels and confirm both succeed.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      type(output_unstr_state_t) :: ostate0, ostate5
      character(len=*), parameter :: fname0 = "test_unstr_deflate0.nc"
      character(len=*), parameter :: fname5 = "test_unstr_deflate5.nc"

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh)
      state%barotropic%stage_cv = 1.0_wp
      state%barotropic%bed_cv = -10.0_wp

      call output_unstr_init(ostate0, fname0, state, mesh, deflate_level=0)
      call output_unstr_write_snapshot(ostate0, state, 0.0_wp)
      call output_unstr_close(ostate0)

      call output_unstr_init(ostate5, fname5, state, mesh, deflate_level=5)
      call output_unstr_write_snapshot(ostate5, state, 0.0_wp)
      call output_unstr_close(ostate5)

      call check(error, .true., "both deflate levels should succeed")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
      call delete_file(fname0)
      call delete_file(fname5)
   end subroutine test_deflate_levels

end module test_output_unstr
