!! Unit tests for the solver driver (SSP-RK2 integration)
module test_solver
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_solver_tests

contains

   subroutine collect_solver_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("quiescent_state_unchanged", test_quiescent), &
                  new_unittest("mass_conserved_dam_break", test_mass_conservation), &
                  new_unittest("dt_max_respected", test_dt_max), &
                  new_unittest("no_negative_depth", test_no_negative_depth), &
                  new_unittest("dt_positive_after_step", test_dt_positive), &
                  new_unittest("solver_1d_ny1", test_solver_1d_ny1) &
                  ]
   end subroutine collect_solver_tests

   subroutine test_quiescent(error)
      !! A uniform quiescent state should remain unchanged after stepping
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      integer :: i, j, n
      real(wp) :: max_change

      call s%init(10, 10, 2, 1.0_wp, 1.0_wp, 0.45_wp)

      s%barotropic%h = 1.0_wp
      s%barotropic%hu = 0.0_wp
      s%barotropic%hv = 0.0_wp
      s%barotropic%b = 0.0_wp

      call solver_enter_data(s)

      do n = 1, 5
         call solver_step(s)
      end do

      call solver_exit_data(s)

      max_change = 0.0_wp
      do j = s%grid%nghost + 1, s%grid%ny_total - s%grid%nghost
         do i = s%grid%nghost + 1, s%grid%nx_total - s%grid%nghost
            max_change = max(max_change, abs(s%barotropic%h(i, j) - 1.0_wp))
            max_change = max(max_change, abs(s%barotropic%hu(i, j)))
            max_change = max(max_change, abs(s%barotropic%hv(i, j)))
         end do
      end do

      call check(error, max_change < 1.0e-12_wp, &
                 "Quiescent state should remain unchanged after stepping")

      call s%destroy()
   end subroutine test_quiescent

   subroutine test_mass_conservation(error)
      !! Total mass should be conserved after several steps with wall BCs
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      integer :: i, j, n
      real(wp) :: mass_before, mass_after

      call s%init(20, 10, 2, 1.0_wp, 1.0_wp, 0.45_wp)

      ! Dam break initial condition
      s%barotropic%b = 0.0_wp
      s%barotropic%hu = 0.0_wp
      s%barotropic%hv = 0.0_wp
      do j = 1, s%grid%ny_total
         do i = 1, s%grid%nx_total
            if (i <= s%grid%nghost + s%grid%nx_phys/2) then
               s%barotropic%h(i, j) = 2.0_wp
            else
               s%barotropic%h(i, j) = 1.0_wp
            end if
         end do
      end do

      ! Compute initial mass
      mass_before = 0.0_wp
      do j = s%grid%nghost + 1, s%grid%ny_total - s%grid%nghost
         do i = s%grid%nghost + 1, s%grid%nx_total - s%grid%nghost
            mass_before = mass_before + s%barotropic%h(i, j)
         end do
      end do
      mass_before = mass_before*s%grid%dx*s%grid%dy

      call solver_enter_data(s)
      do n = 1, 10
         call solver_step(s)
      end do
      call solver_exit_data(s)

      mass_after = 0.0_wp
      do j = s%grid%nghost + 1, s%grid%ny_total - s%grid%nghost
         do i = s%grid%nghost + 1, s%grid%nx_total - s%grid%nghost
            mass_after = mass_after + s%barotropic%h(i, j)
         end do
      end do
      mass_after = mass_after*s%grid%dx*s%grid%dy

      call check(error, abs(mass_after - mass_before)/mass_before < 1.0e-10_wp, &
                 "Mass should be conserved to machine precision")

      call s%destroy()
   end subroutine test_mass_conservation

   subroutine test_dt_max(error)
      !! dt_max argument should cap the timestep
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s

      call s%init(10, 10, 2, 100.0_wp, 100.0_wp, 0.45_wp)

      s%barotropic%h = 1.0_wp
      s%barotropic%hu = 0.0_wp
      s%barotropic%hv = 0.0_wp
      s%barotropic%b = 0.0_wp

      call solver_enter_data(s)
      call solver_step(s, dt_max=0.001_wp)
      call solver_exit_data(s)

      call check(error, s%dt <= 0.001_wp + 1.0e-15_wp, &
                 "dt should be capped by dt_max")

      call s%destroy()
   end subroutine test_dt_max

   subroutine test_no_negative_depth(error)
      !! No cells should have negative depth after stepping
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      integer :: i, j, n
      real(wp) :: min_h

      call s%init(20, 10, 2, 1.0_wp, 1.0_wp, 0.4_wp)

      ! Dam break into dry bed
      s%barotropic%b = 0.0_wp
      s%barotropic%hu = 0.0_wp
      s%barotropic%hv = 0.0_wp
      s%barotropic%h = 0.0_wp
      do j = 1, s%grid%ny_total
         do i = 1, s%grid%nghost + s%grid%nx_phys/2
            s%barotropic%h(i, j) = 2.0_wp
         end do
      end do

      call solver_enter_data(s)
      do n = 1, 20
         call solver_step(s)
      end do
      call solver_exit_data(s)

      min_h = 0.0_wp
      do j = s%grid%nghost + 1, s%grid%ny_total - s%grid%nghost
         do i = s%grid%nghost + 1, s%grid%nx_total - s%grid%nghost
            min_h = min(min_h, s%barotropic%h(i, j))
         end do
      end do

      call check(error, min_h >= 0.0_wp, &
                 "No cell should have negative depth after stepping")

      call s%destroy()
   end subroutine test_no_negative_depth

   subroutine test_dt_positive(error)
      !! dt should be positive after a solver step
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s

      call s%init(10, 10, 2, 1.0_wp, 1.0_wp, 0.45_wp)

      s%barotropic%h = 1.0_wp
      s%barotropic%hu = 0.0_wp
      s%barotropic%hv = 0.0_wp
      s%barotropic%b = 0.0_wp

      call solver_enter_data(s)
      call solver_step(s)
      call solver_exit_data(s)

      call check(error, s%dt > 0.0_wp, "dt should be positive after step")
      if (allocated(error)) return
      call check(error, s%dt < 1.0e10_wp, "dt should be finite after step")

      call s%destroy()
   end subroutine test_dt_positive

   subroutine test_solver_1d_ny1(error)
      !! Run a few steps with ny=1 to verify no crash and mass conservation
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      integer :: i, j, n
      real(wp) :: mass_before, mass_after

      call s%init(20, 1, 2, 1.0_wp, 1.0_wp, 0.45_wp)

      ! Dam break IC
      s%barotropic%b = 0.0_wp
      s%barotropic%hu = 0.0_wp
      s%barotropic%hv = 0.0_wp
      do j = 1, s%grid%ny_total
         do i = 1, s%grid%nx_total
            if (i <= s%grid%nghost + s%grid%nx_phys/2) then
               s%barotropic%h(i, j) = 2.0_wp
            else
               s%barotropic%h(i, j) = 1.0_wp
            end if
         end do
      end do

      ! Initial mass
      mass_before = 0.0_wp
      do j = s%grid%nghost + 1, s%grid%ny_total - s%grid%nghost
         do i = s%grid%nghost + 1, s%grid%nx_total - s%grid%nghost
            mass_before = mass_before + s%barotropic%h(i, j)*s%grid%dx*s%grid%dy
         end do
      end do

      call solver_enter_data(s)
      do n = 1, 5
         call solver_step(s)
      end do
      call solver_exit_data(s)

      mass_after = 0.0_wp
      do j = s%grid%nghost + 1, s%grid%ny_total - s%grid%nghost
         do i = s%grid%nghost + 1, s%grid%nx_total - s%grid%nghost
            mass_after = mass_after + s%barotropic%h(i, j)*s%grid%dx*s%grid%dy
         end do
      end do

      call check(error, abs(mass_after - mass_before)/mass_before < 1.0e-10_wp, &
                 "Mass conserved with ny=1 (1D grid)")

      call s%destroy()
   end subroutine test_solver_1d_ny1

end module test_solver
