!! Unit tests for the unstructured gauge output module
module test_gauges_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t, &
                             mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy
   use rki_state_unstr, only: unstr_state_t, &
                              state_unstr_init, state_unstr_destroy
   use rki_gauges_unstr, only: gauge_unstr_t, gauges_unstr_init, &
                               gauges_unstr_write, gauges_unstr_close
   use netcdf, only: nf90_open, nf90_close, nf90_nowrite, nf90_noerr, &
                     nf90_inq_varid, nf90_inq_dimid, nf90_inquire_dimension
   implicit none
   private

   public :: collect_gauges_unstr_tests

contains

   subroutine collect_gauges_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("init_write_close", test_init_write_close), &
                  new_unittest("multiple_gauges", test_multiple_gauges), &
                  new_unittest("gauge_outside_skipped", test_outside_mesh) &
                  ]
   end subroutine collect_gauges_unstr_tests

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

   subroutine test_init_write_close(error)
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      type(gauge_unstr_t) :: gauges
      character(len=*), parameter :: fname = "test_gauges_unstr.nc"
      real(wp) :: gx(1), gy(1)
      integer :: ncid, dimid, varid, status, ntime, k
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh)

         do k = 1, mesh%ncells_owned
            state%barotropic%stage_cv(k) = 0.5_wp
            state%barotropic%bed_cv(k) = -1.0_wp
            state%barotropic%xmom_cv(k) = 0.5_wp
            state%barotropic%ymom_cv(k) = -0.25_wp
         end do

         ! One gauge in the centre of the unit square (the apex vertex)
         gx(1) = 0.5_wp
         gy(1) = 0.4_wp

         call gauges_unstr_init(gauges, fname, mesh, 1, gx, gy)
         call gauges_unstr_write(gauges, state, 0.0_wp)
         call gauges_unstr_write(gauges, state, 0.5_wp)
         call gauges_unstr_close(gauges)

         status = nf90_open(fname, nf90_nowrite, ncid)
         call check(error, status == nf90_noerr, "open gauge ouput")
         if (allocated(error)) exit checks

         status = nf90_inq_dimid(ncid, "time", dimid)
         status = nf90_inquire_dimension(ncid, dimid, len=ntime)
         call check(error, ntime == 2, "expected 2 gauge time records")
         if (allocated(error)) then
            status = nf90_close(ncid); exit checks
         end if

         status = nf90_inq_varid(ncid, "eta", varid)
         call check(error, status == nf90_noerr, "eta variable should exist")
         if (allocated(error)) then
            status = nf90_close(ncid); exit checks
         end if
         status = nf90_close(ncid)

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
      call delete_file(fname)
   end subroutine test_init_write_close

   subroutine test_multiple_gauges(error)
      !! Place gauges in three different triangles of the 4-triangle mesh.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      type(gauge_unstr_t) :: gauges
      character(len=*), parameter :: fname = "test_gauges_multi.nc"
      real(wp) :: gx(3), gy(3)
      integer :: k
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh)

         do k = 1, mesh%ncells_owned
            state%barotropic%stage_cv(k) = real(k, wp)
            state%barotropic%bed_cv(k) = -1.0_wp
            state%barotropic%xmom_cv(k) = 0.0_wp
            state%barotropic%ymom_cv(k) = 0.0_wp
         end do

         ! Triangle 1 (bottom): vertices (0,0)-(1,0)-(0.5,0.5) — sample at (0.5, 0.1)
         ! Triangle 2 (right):  (1,0)-(1,1)-(0.5,0.5)           — sample at (0.9, 0.5)
         ! Triangle 3 (top):    (1,1)-(0,1)-(0.5,0.5)           — sample at (0.5, 0.9)
         gx = [0.5_wp, 0.9_wp, 0.5_wp]
         gy = [0.1_wp, 0.5_wp, 0.9_wp]

         call gauges_unstr_init(gauges, fname, mesh, 3, gx, gy)
         call check(error, gauges%n_local == 3, "expected 3 gauges to land in owned cells")
         if (allocated(error)) exit checks

         call gauges_unstr_write(gauges, state, 0.0_wp)
         call gauges_unstr_close(gauges)

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
      call delete_file(fname)
   end subroutine test_multiple_gauges

   subroutine test_outside_mesh(error)
      !! A gauge well outside the unit square should be skipped silently.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      type(gauge_unstr_t) :: gauges
      character(len=*), parameter :: fname = "test_gauges_outside.nc"
      real(wp) :: gx(2), gy(2)
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh)

         ! First gauge inside, second outside the unit square
         gx = [0.5_wp, 10.0_wp]
         gy = [0.4_wp, 10.0_wp]

         call gauges_unstr_init(gauges, fname, mesh, 2, gx, gy)
         call check(error, gauges%n_local == 1, "out-of-mesh gauge should be skipped")
         if (allocated(error)) exit checks

         call gauges_unstr_close(gauges)

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
      call delete_file(fname)
   end subroutine test_outside_mesh

end module test_gauges_unstr
