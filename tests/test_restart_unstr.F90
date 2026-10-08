!! Unit tests for the unstructured restart writer/reader
module test_restart_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t, &
                             mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy
   use rki_state_unstr, only: unstr_state_t, &
                              state_unstr_init, state_unstr_destroy
   use rki_restart_unstr, only: restart_write_unstr, restart_read_unstr
   implicit none
   private

   public :: collect_restart_unstr_tests

contains

   subroutine collect_restart_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("barotropic_round_trip", test_barotropic_round_trip), &
                  new_unittest("step_count_persists", test_step_count) &
                  ]
   end subroutine collect_restart_unstr_tests

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

   subroutine test_barotropic_round_trip(error)
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: src, dst
      character(len=*), parameter :: fname = "test_restart_unstr.nc"
      real(wp), parameter :: T_IN = 7.5_wp
      real(wp) :: t_out
      integer :: ierr, k
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(src, mesh)
         call state_unstr_init(dst, mesh)

         do k = 1, mesh%ncells_owned
            src%barotropic%stage_cv(k) = 1.0_wp + 0.1_wp*real(k, wp)
            src%barotropic%xmom_cv(k) = 0.5_wp*real(k, wp)
            src%barotropic%ymom_cv(k) = -0.25_wp*real(k, wp)
            src%barotropic%bed_cv(k) = -10.0_wp - 0.05_wp*real(k, wp)
         end do

         call restart_write_unstr(fname, src, mesh, T_IN)
         call restart_read_unstr(fname, dst, t_out, ierr)

         call check(error, ierr == 0, "restart_read_unstr should succeed")
         if (allocated(error)) exit checks

         call check(error, abs(t_out - T_IN) < 1.0e-12_wp, "time round-trips")
         if (allocated(error)) exit checks

         do k = 1, mesh%ncells_owned
            call check(error, abs(src%barotropic%stage_cv(k) - dst%barotropic%stage_cv(k)) < 1.0e-12_wp, &
                       "stage_cv round-trip")
            if (allocated(error)) exit checks
            call check(error, abs(src%barotropic%xmom_cv(k) - dst%barotropic%xmom_cv(k)) < 1.0e-12_wp, &
                       "xmom_cv round-trip")
            if (allocated(error)) exit checks
            call check(error, abs(src%barotropic%ymom_cv(k) - dst%barotropic%ymom_cv(k)) < 1.0e-12_wp, &
                       "ymom_cv round-trip")
            if (allocated(error)) exit checks
            call check(error, abs(src%barotropic%bed_cv(k) - dst%barotropic%bed_cv(k)) < 1.0e-12_wp, &
                       "bed_cv round-trip")
            if (allocated(error)) exit checks
         end do

      end block checks
      continue
      call state_unstr_destroy(src)
      call state_unstr_destroy(dst)
      call mesh_unstr_destroy(mesh)
      call delete_file(fname)
   end subroutine test_barotropic_round_trip

   subroutine test_step_count(error)
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: src, dst
      character(len=*), parameter :: fname = "test_restart_unstr_step.nc"
      real(wp) :: t_out
      integer :: ierr
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(src, mesh)
         call state_unstr_init(dst, mesh)

         src%barotropic%stage_cv = 1.0_wp
         src%barotropic%xmom_cv = 0.0_wp
         src%barotropic%ymom_cv = 0.0_wp
         src%barotropic%bed_cv = -1.0_wp
         src%step_count = 42

         call restart_write_unstr(fname, src, mesh, 1.0_wp)
         call restart_read_unstr(fname, dst, t_out, ierr)

         call check(error, ierr == 0, "restart_read should succeed")
         if (allocated(error)) exit checks

         call check(error, dst%step_count == 42, "step_count must round-trip")
         if (allocated(error)) exit checks

      end block checks
      continue
      call state_unstr_destroy(src)
      call state_unstr_destroy(dst)
      call mesh_unstr_destroy(mesh)
      call delete_file(fname)
   end subroutine test_step_count

end module test_restart_unstr
