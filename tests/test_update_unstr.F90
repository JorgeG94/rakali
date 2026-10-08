!! Unit tests for the explicit update / saxpy / backup-restore kernels
module test_update_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t, &
                             mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, &
                             mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, &
                              state_unstr_init, &
                              state_unstr_destroy, &
                              state_unstr_enter_data, &
                              state_unstr_exit_data
   use rki_kernel_update_unstr, only: backup_conserved_unstr, &
                                      restore_conserved_unstr, &
                                      apply_explicit_update_unstr, &
                                      saxpy_conserved_unstr
   implicit none
   private

   public :: collect_update_tests

   real(wp), parameter :: TOL = 1.0e-12_wp

contains

   subroutine collect_update_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("explicit_update_hand_computed", test_explicit_update_hand), &
                  new_unittest("explicit_update_zero_dt_no_op", test_explicit_update_zero_dt), &
                  new_unittest("backup_then_restore_round_trip", test_backup_restore), &
                  new_unittest("saxpy_identity_a1_b0", test_saxpy_identity), &
                  new_unittest("saxpy_restore_a0_b1", test_saxpy_restore), &
                  new_unittest("saxpy_average_a05_b05", test_saxpy_average), &
                  new_unittest("ssp_rk2_average_idempotent", test_ssp_rk2_average) &
                  ]
   end subroutine collect_update_tests

   ! ----------------------------------------------------------------------
   ! Helpers
   ! ----------------------------------------------------------------------

   subroutine build_two_triangle_square(mesh)
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp) :: vx(4), vy(4)
      integer :: tv(3, 2)

      vx = [0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp]
      vy = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp]
      tv(:, 1) = [1, 2, 3]
      tv(:, 2) = [1, 3, 4]
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_two_triangle_square

   subroutine run_explicit_update(state, mesh, dt)
      !! Wrap apply_explicit_update_unstr with the data dance.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: dt

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call apply_explicit_update_unstr(state, mesh, dt)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_explicit_update

   subroutine run_backup(state, mesh)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call backup_conserved_unstr(state, mesh)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_backup

   subroutine run_restore(state, mesh)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call restore_conserved_unstr(state, mesh)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_restore

   subroutine run_saxpy(state, mesh, a, b)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: a, b

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call saxpy_conserved_unstr(state, mesh, a, b)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_saxpy

   ! ----------------------------------------------------------------------
   ! Tests
   ! ----------------------------------------------------------------------

   subroutine test_explicit_update_hand(error)
      !! Q := Q + dt * Q_eu, hand-computed on a 2-cell mesh.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: dt

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Initial state.
      state%barotropic%stage_cv = [1.0_wp, 2.0_wp]
      state%barotropic%xmom_cv = [0.5_wp, -0.5_wp]
      state%barotropic%ymom_cv = [0.0_wp, 0.3_wp]

      ! Pre-fill explicit-update arrays directly. These are now
      ! map(to:) so the host pre-fill flows to the device.
      state%barotropic%stage_eu = [0.1_wp, -0.2_wp]
      state%barotropic%xmom_eu = [-0.4_wp, 0.6_wp]
      state%barotropic%ymom_eu = [1.0_wp, -1.0_wp]

      dt = 0.5_wp

      call run_explicit_update(state, mesh, dt)

      ! Q_new = Q_old + dt * Q_eu
      call check(error, abs(state%barotropic%stage_cv(1) - (1.0_wp + 0.5_wp*0.1_wp)) < TOL, &
                 "stage_cv(1) = 1 + 0.5*0.1 = 1.05")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%stage_cv(2) - (2.0_wp + 0.5_wp*(-0.2_wp))) < TOL, &
                 "stage_cv(2) = 2 + 0.5*(-0.2) = 1.9")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(1) - (0.5_wp + 0.5_wp*(-0.4_wp))) < TOL, &
                 "xmom_cv(1) = 0.5 - 0.2 = 0.3")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(2) - (-0.5_wp + 0.5_wp*0.6_wp)) < TOL, &
                 "xmom_cv(2) = -0.5 + 0.3 = -0.2")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(1) - (0.0_wp + 0.5_wp*1.0_wp)) < TOL, &
                 "ymom_cv(1) = 0 + 0.5 = 0.5")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(2) - (0.3_wp + 0.5_wp*(-1.0_wp))) < TOL, &
                 "ymom_cv(2) = 0.3 - 0.5 = -0.2")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_explicit_update_hand

   subroutine test_explicit_update_zero_dt(error)
      !! dt = 0 must leave the conserved state untouched, even if eu
      !! is nonzero.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = [3.0_wp, 4.0_wp]
      state%barotropic%xmom_cv = [0.7_wp, 0.8_wp]
      state%barotropic%ymom_cv = [-0.1_wp, -0.2_wp]
      state%barotropic%stage_eu = [99.0_wp, -99.0_wp]
      state%barotropic%xmom_eu = [99.0_wp, -99.0_wp]
      state%barotropic%ymom_eu = [99.0_wp, -99.0_wp]

      call run_explicit_update(state, mesh, 0.0_wp)

      call check(error, abs(state%barotropic%stage_cv(1) - 3.0_wp) < TOL, "stage_cv(1) unchanged")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%stage_cv(2) - 4.0_wp) < TOL, "stage_cv(2) unchanged")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(1) - 0.7_wp) < TOL, "xmom_cv(1) unchanged")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(2) + 0.2_wp) < TOL, "ymom_cv(2) unchanged")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_explicit_update_zero_dt

   subroutine test_backup_restore(error)
      !! backup -> mutate -> restore should return to original values.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = [5.0_wp, 6.0_wp]
      state%barotropic%xmom_cv = [0.1_wp, 0.2_wp]
      state%barotropic%ymom_cv = [0.3_wp, 0.4_wp]

      call run_backup(state, mesh)

      ! Mutate the conserved state on the host.
      state%barotropic%stage_cv = [99.0_wp, 99.0_wp]
      state%barotropic%xmom_cv = [99.0_wp, 99.0_wp]
      state%barotropic%ymom_cv = [99.0_wp, 99.0_wp]

      call run_restore(state, mesh)

      call check(error, abs(state%barotropic%stage_cv(1) - 5.0_wp) < TOL, "stage_cv(1) restored")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%stage_cv(2) - 6.0_wp) < TOL, "stage_cv(2) restored")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(1) - 0.1_wp) < TOL, "xmom_cv(1) restored")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(2) - 0.4_wp) < TOL, "ymom_cv(2) restored")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_backup_restore

   subroutine test_saxpy_identity(error)
      !! saxpy(a=1, b=0): Q := 1*Q + 0*backup = Q. Identity.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = [1.5_wp, 2.5_wp]
      state%barotropic%xmom_cv = [0.7_wp, -0.3_wp]
      state%barotropic%ymom_cv = [0.1_wp, 0.4_wp]
      state%barotropic%stage_backup = [99.0_wp, 99.0_wp]
      state%barotropic%xmom_backup = [99.0_wp, 99.0_wp]
      state%barotropic%ymom_backup = [99.0_wp, 99.0_wp]

      call run_saxpy(state, mesh, 1.0_wp, 0.0_wp)

      call check(error, abs(state%barotropic%stage_cv(1) - 1.5_wp) < TOL, "stage_cv(1) identity")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%stage_cv(2) - 2.5_wp) < TOL, "stage_cv(2) identity")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(1) - 0.7_wp) < TOL, "xmom_cv(1) identity")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_saxpy_identity

   subroutine test_saxpy_restore(error)
      !! saxpy(a=0, b=1): Q := 0*Q + 1*backup = backup. Pure restore.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = [99.0_wp, 99.0_wp]
      state%barotropic%xmom_cv = [99.0_wp, 99.0_wp]
      state%barotropic%ymom_cv = [99.0_wp, 99.0_wp]
      state%barotropic%stage_backup = [1.0_wp, 2.0_wp]
      state%barotropic%xmom_backup = [0.5_wp, -0.5_wp]
      state%barotropic%ymom_backup = [0.1_wp, -0.1_wp]

      call run_saxpy(state, mesh, 0.0_wp, 1.0_wp)

      call check(error, abs(state%barotropic%stage_cv(1) - 1.0_wp) < TOL, "stage_cv(1) = backup")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%stage_cv(2) - 2.0_wp) < TOL, "stage_cv(2) = backup")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(1) - 0.5_wp) < TOL, "xmom_cv(1) = backup")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(2) + 0.1_wp) < TOL, "ymom_cv(2) = backup")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_saxpy_restore

   subroutine test_saxpy_average(error)
      !! saxpy(a=0.5, b=0.5): the SSP-RK2 final-average step.
      !! Q := 0.5*Q + 0.5*backup -- exactly halfway between current
      !! and saved state.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = [3.0_wp, 5.0_wp]
      state%barotropic%stage_backup = [1.0_wp, 7.0_wp]
      state%barotropic%xmom_cv = [0.4_wp, 0.6_wp]
      state%barotropic%xmom_backup = [-0.4_wp, 1.4_wp]
      state%barotropic%ymom_cv = [0.0_wp, 0.0_wp]
      state%barotropic%ymom_backup = [0.0_wp, 0.0_wp]

      call run_saxpy(state, mesh, 0.5_wp, 0.5_wp)

      call check(error, abs(state%barotropic%stage_cv(1) - 2.0_wp) < TOL, &
                 "0.5*3 + 0.5*1 = 2")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%stage_cv(2) - 6.0_wp) < TOL, &
                 "0.5*5 + 0.5*7 = 6")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(1) - 0.0_wp) < TOL, &
                 "0.5*0.4 + 0.5*(-0.4) = 0")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(2) - 1.0_wp) < TOL, &
                 "0.5*0.6 + 0.5*1.4 = 1.0")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_saxpy_average

   subroutine test_ssp_rk2_average(error)
      !! Smoke test of the full SSP-RK2 average pattern: backup, then
      !! pretend the kernel doubled the state (Q -> 2*Q), then average
      !! with backup. Result should be (2*Q + Q)/2 = 1.5*Q.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = [4.0_wp, 6.0_wp]
      state%barotropic%xmom_cv = [0.0_wp, 0.0_wp]
      state%barotropic%ymom_cv = [0.0_wp, 0.0_wp]

      call run_backup(state, mesh)

      ! Pretend two RK substages doubled the stage.
      state%barotropic%stage_cv = [8.0_wp, 12.0_wp]

      call run_saxpy(state, mesh, 0.5_wp, 0.5_wp)

      call check(error, abs(state%barotropic%stage_cv(1) - 6.0_wp) < TOL, &
                 "(8 + 4)/2 = 6")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%stage_cv(2) - 9.0_wp) < TOL, &
                 "(12 + 6)/2 = 9")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_ssp_rk2_average

end module test_update_unstr
