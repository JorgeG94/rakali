!! Dam break validation test against the exact Riemann solution
module test_dam_break
   !! Runs a 1D dam break on a 2D grid and compares against the
   !! exact Riemann solution (Toro Ch.5). Reports L2 error and
   !! verifies mass conservation.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use ritter_dambreak, only: exact_riemann_dambreak, exact_riemann_sample
   implicit none
   private

   public :: collect_dam_break_tests

   integer, parameter :: NX = 400
   integer, parameter :: NY = 4
   integer, parameter :: NGHOST = 2
   real(wp), parameter :: DX = 0.5_wp
   real(wp), parameter :: DY = 0.5_wp
   real(wp), parameter :: CFL_NUM = 0.45_wp
   real(wp), parameter :: T_END = 1.0_wp
   real(wp), parameter :: H_LEFT = 2.0_wp
   real(wp), parameter :: H_RIGHT = 1.0_wp
   real(wp), parameter :: X_DAM = real(NX, wp)*DX*0.5_wp

contains

   subroutine collect_dam_break_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("dam_break_mass_conserved", test_mass), &
                  new_unittest("dam_break_l2_accuracy", test_l2_accuracy) &
                  ]
   end subroutine collect_dam_break_tests

   subroutine run_dam_break(state, t_final)
      !! Run the dam break simulation to t_final
      type(state_t), intent(inout) :: state
      real(wp), intent(out) :: t_final

      integer :: i, j, mid_x

      call state%init(NX, NY, NGHOST, DX, DY, CFL_NUM)

      mid_x = NGHOST + NX/2
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%b(i, j) = 0.0_wp
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
            if (i <= mid_x) then
               state%barotropic%h(i, j) = H_LEFT
            else
               state%barotropic%h(i, j) = H_RIGHT
            end if
         end do
      end do

      call solver_enter_data(state)
      t_final = 0.0_wp
      do while (t_final < T_END)
         call solver_step(state, dt_max=(T_END - t_final))
         t_final = t_final + state%dt
      end do
      call solver_exit_data(state)

   end subroutine run_dam_break

   function compute_mass(s) result(mass)
      type(state_t), intent(in) :: s
      real(wp) :: mass
      integer :: i, j

      mass = 0.0_wp
      do j = s%grid%nghost + 1, s%grid%ny_total - s%grid%nghost
         do i = s%grid%nghost + 1, s%grid%nx_total - s%grid%nghost
            mass = mass + s%barotropic%h(i, j)*s%grid%dx*s%grid%dy
         end do
      end do
   end function compute_mass

   subroutine test_mass(error)
      !! Mass should be conserved to near machine precision
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: mass_initial, mass_final, mass_error, t_final
      integer :: i, j, mid_x

      call state%init(NX, NY, NGHOST, DX, DY, CFL_NUM)

      mid_x = NGHOST + NX/2
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%b(i, j) = 0.0_wp
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
            if (i <= mid_x) then
               state%barotropic%h(i, j) = H_LEFT
            else
               state%barotropic%h(i, j) = H_RIGHT
            end if
         end do
      end do

      mass_initial = compute_mass(state)

      call solver_enter_data(state)
      t_final = 0.0_wp
      do while (t_final < T_END)
         call solver_step(state, dt_max=(T_END - t_final))
         t_final = t_final + state%dt
      end do
      call solver_exit_data(state)

      mass_final = compute_mass(state)
      mass_error = abs(mass_final - mass_initial)/mass_initial

      call check(error, mass_error < 1.0e-10_wp, &
                 "Dam break mass conservation violated: "// &
                 "relative error should be < 1e-10")

      call state%destroy()
   end subroutine test_mass

   subroutine test_l2_accuracy(error)
      !! L2 error against exact solution should be below 5%
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: h_star, u_star, t_final
      real(wp) :: l2_err, x, h_exact, u_exact, err_val
      integer :: i, j_mid

      call exact_riemann_dambreak(H_LEFT, H_RIGHT, h_star, u_star)
      call run_dam_break(state, t_final)

      j_mid = NGHOST + NY/2 + 1
      l2_err = 0.0_wp
      do i = NGHOST + 1, state%grid%nx_total - NGHOST
         x = (real(i - NGHOST, wp) - 0.5_wp)*DX
         call exact_riemann_sample(H_LEFT, H_RIGHT, h_star, u_star, &
                                   X_DAM, t_final, x, h_exact, u_exact)
         err_val = abs(state%barotropic%h(i, j_mid) - h_exact)
         l2_err = l2_err + err_val*err_val
      end do
      l2_err = sqrt(l2_err/real(NX, wp))

      ! 1D Ritter dam break on 400 cells at t=1s with first-order HLL
      ! reliably lands around 1-2% L2 relative error. The earlier 5%
      ! bound accepted a doubling of numerical dissipation without
      ! complaint — 2.5% still has 25-50% headroom above the noise
      ! floor without being a rubber stamp.
      call check(error, l2_err/((H_LEFT + H_RIGHT)/2.0_wp) < 0.025_wp, &
                 "Dam break L2 error too large (> 2.5%)")

      call state%destroy()
   end subroutine test_l2_accuracy

end module test_dam_break
