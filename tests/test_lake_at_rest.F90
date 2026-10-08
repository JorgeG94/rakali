!! Lake at rest (well-balanced) test
module test_lake_at_rest
   !! Verifies the C-property: a flat water surface over non-flat
   !! bathymetry should remain exactly at rest.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_lake_at_rest_tests

   integer, parameter :: NX = 100
   integer, parameter :: NY = 100
   integer, parameter :: NGHOST = 2
   real(wp), parameter :: DX = 1.0_wp
   real(wp), parameter :: DY = 1.0_wp
   real(wp), parameter :: CFL_NUM = 0.45_wp
   integer, parameter :: N_STEPS_RUN = 100
   real(wp), parameter :: ETA_SURFACE = 1.0_wp
   real(wp), parameter :: BATHY_MAX = 0.5_wp
   real(wp), parameter :: VELOCITY_TOL = 1.0e-10_wp
   real(wp), parameter :: ETA_TOL = 1.0e-10_wp

contains

   subroutine collect_lake_at_rest_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("lake_at_rest_no_spurious_velocity", test_no_velocity), &
                  new_unittest("lake_at_rest_flat_surface", test_flat_surface) &
                  ]
   end subroutine collect_lake_at_rest_tests

   subroutine setup_and_run(state)
      !! Set up parabolic bathymetry with flat surface and run 100 steps
      type(state_t), intent(inout) :: state

      integer :: i, j, n
      real(wp) :: x_center, y_center, x_pos, y_pos, r_sq

      call state%init(NX, NY, NGHOST, DX, DY, CFL_NUM)

      x_center = real(NX, wp)*DX/2.0_wp
      y_center = real(NY, wp)*DY/2.0_wp

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            x_pos = (real(i - NGHOST, wp) - 0.5_wp)*DX
            y_pos = (real(j - NGHOST, wp) - 0.5_wp)*DY

            r_sq = ((x_pos - x_center)/x_center)**2 &
                   + ((y_pos - y_center)/y_center)**2
            state%barotropic%b(i, j) = BATHY_MAX*max(1.0_wp - r_sq, 0.0_wp)
            state%barotropic%h(i, j) = max(ETA_SURFACE - state%barotropic%b(i, j), 0.0_wp)
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
         end do
      end do

      call solver_enter_data(state)
      do n = 1, N_STEPS_RUN
         call solver_step(state)
      end do
      call solver_exit_data(state)
   end subroutine setup_and_run

   subroutine test_no_velocity(error)
      !! No spurious momentum should develop
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      integer :: i, j
      real(wp) :: max_hu, max_hv

      call setup_and_run(state)

      max_hu = 0.0_wp
      max_hv = 0.0_wp
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            max_hu = max(max_hu, abs(state%barotropic%hu(i, j)))
            max_hv = max(max_hv, abs(state%barotropic%hv(i, j)))
         end do
      end do

      call check(error, max_hu < VELOCITY_TOL, &
                 "Spurious x-momentum in lake at rest")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      call check(error, max_hv < VELOCITY_TOL, &
                 "Spurious y-momentum in lake at rest")

      call state%destroy()
   end subroutine test_no_velocity

   subroutine test_flat_surface(error)
      !! Free surface should remain flat
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      integer :: i, j
      real(wp) :: max_eta_err, eta_now

      call setup_and_run(state)

      max_eta_err = 0.0_wp
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            eta_now = state%barotropic%h(i, j) + state%barotropic%b(i, j)
            max_eta_err = max(max_eta_err, abs(eta_now - ETA_SURFACE))
         end do
      end do

      call check(error, max_eta_err < ETA_TOL, &
                 "Surface elevation changed in lake at rest")

      call state%destroy()
   end subroutine test_flat_surface

end module test_lake_at_rest
