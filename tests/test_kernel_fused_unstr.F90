!! Unit tests for the fused unstructured update kernel.
!!
!! Drives `fused_update_stage1_unstr` and `fused_update_stage2_unstr`
!! against a small triangular mesh with synthetic Euler-update arrays
!! (stage_eu, xmom_eu, ymom_eu) so the friction / wind / protect /
!! RK2-average branches are exercised.
module test_kernel_fused_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, RHO_WATER
   use rki_mesh_unstr, only: mesh_unstr_t, &
                             mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy
   use rki_state_unstr, only: unstr_state_t, &
                              state_unstr_init, state_unstr_destroy
   use rki_kernel_fused_unstr, only: fused_update_stage1_unstr, &
                                     fused_update_stage2_unstr
   implicit none
   private

   public :: collect_kernel_fused_unstr_tests

contains

   subroutine collect_kernel_fused_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("stage1_save_and_euler_update", test_stage1_basic), &
                  new_unittest("stage1_friction_damps_momentum", test_friction), &
                  new_unittest("stage1_wind_adds_momentum", test_wind), &
                  new_unittest("stage1_wind_field_per_cell", test_wind_field_per_cell), &
                  new_unittest("stage1_protect_dries_thin_cells", test_protect_dry), &
                  new_unittest("stage1_negative_h_clamped_to_bed", test_protect_negative), &
                  new_unittest("stage2_rk2_average", test_stage2_rk2) &
                  ]
   end subroutine collect_kernel_fused_unstr_tests

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

   subroutine prime_state(state, stage, xmom, ymom, bed)
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(in) :: stage, xmom, ymom, bed
      state%barotropic%stage_cv = stage
      state%barotropic%xmom_cv = xmom
      state%barotropic%ymom_cv = ymom
      state%barotropic%bed_cv = bed
      state%barotropic%stage_eu = 0.0_wp
      state%barotropic%xmom_eu = 0.0_wp
      state%barotropic%ymom_eu = 0.0_wp
      state%barotropic%manning_n = 0.0_wp
      state%barotropic%wind_tau_x = 0.0_wp
      state%barotropic%wind_tau_y = 0.0_wp
   end subroutine prime_state

   ! ----------------------------------------------------------------------

   subroutine test_stage1_basic(error)
      !! Zero Euler-update with no friction/wind ⇒ state unchanged + backups
      !! contain the input values.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh)
         call prime_state(state, 1.0_wp, 0.5_wp, -0.25_wp, -10.0_wp)

         call fused_update_stage1_unstr(state, mesh, 0.1_wp)

         call check(error, all(abs(state%barotropic%stage_cv(1:mesh%ncells_owned) - 1.0_wp) < 1.0e-12_wp), &
                    "stage_cv unchanged with zero Euler update")
         if (allocated(error)) exit checks
         call check(error, all(abs(state%barotropic%stage_backup(1:mesh%ncells_owned) - 1.0_wp) < 1.0e-12_wp), &
                    "stage_backup snapshot")
         if (allocated(error)) exit checks
         call check(error, all(abs(state%barotropic%xmom_backup(1:mesh%ncells_owned) - 0.5_wp) < 1.0e-12_wp), &
                    "xmom_backup snapshot")
         if (allocated(error)) exit checks

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_stage1_basic

   subroutine test_friction(error)
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh)
      call prime_state(state, 1.0_wp, 1.0_wp, 0.0_wp, -10.0_wp)
      state%barotropic%manning_n = 0.05_wp

      call fused_update_stage1_unstr(state, mesh, 0.5_wp)

      call check(error, state%barotropic%xmom_cv(1) < 1.0_wp .and. &
                 state%barotropic%xmom_cv(1) > 0.0_wp, &
                 "Manning friction must reduce xmom")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_friction

   subroutine test_wind(error)
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh)
      call prime_state(state, 1.0_wp, 0.0_wp, 0.0_wp, -10.0_wp)
      state%barotropic%wind_tau_x = 1.0_wp

      call fused_update_stage1_unstr(state, mesh, 0.1_wp)

      call check(error, state%barotropic%xmom_cv(1) > 0.0_wp, &
                 "wind stress must add xmom")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_wind

   subroutine test_wind_field_per_cell(error)
      !! A per-cell wind-stress field (forcing%wind_x/y + wind_active) drives
      !! spatially-varying momentum: with the uniform wind_tau left at zero,
      !! each owned cell gains exactly dt * wind(cell) / RHO_WATER, so the
      !! increment differs cell-to-cell.  Exercises the gridded-wind path
      !! added for unstructured parity (the uniform `test_wind` above only
      !! covers barotropic%wind_tau_*).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp), parameter :: DT = 0.1_wp
      integer :: k
      real(wp) :: wx, wy, expect_x, expect_y
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh)
         call prime_state(state, 1.0_wp, 0.0_wp, 0.0_wp, -10.0_wp)

         ! Per-cell (spatially varying) wind stress; uniform wind_tau stays 0.
         do k = 1, state%ncells_total
            state%forcing%wind_x(k) = 0.2_wp*real(k, wp)
            state%forcing%wind_y(k) = -0.1_wp*real(k, wp)
         end do
         state%forcing%wind_active = .true.

         call fused_update_stage1_unstr(state, mesh, DT)

         ! xmom starts at 0, no friction, h = 11 > min_h, uniform tau = 0,
         ! so xmom(k) = dt * wind_x(k) / RHO_WATER exactly.
         do k = 1, mesh%ncells_owned
            wx = 0.2_wp*real(k, wp)
            wy = -0.1_wp*real(k, wp)
            expect_x = DT*(wx/RHO_WATER)
            expect_y = DT*(wy/RHO_WATER)
            call check(error, abs(state%barotropic%xmom_cv(k) - expect_x) < 1.0e-12_wp, &
                       "per-cell wind: xmom = dt*wind_x/RHO_WATER")
            if (allocated(error)) exit checks
            call check(error, abs(state%barotropic%ymom_cv(k) - expect_y) < 1.0e-12_wp, &
                       "per-cell wind: ymom = dt*wind_y/RHO_WATER")
            if (allocated(error)) exit checks
         end do

         ! Strictly spatial: distinct cells get distinct momentum.
         call check(error, abs(state%barotropic%xmom_cv(1) - state%barotropic%xmom_cv(2)) > 1.0e-9_wp, &
                    "per-cell wind drives spatially-varying momentum")

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_wind_field_per_cell

   subroutine test_protect_dry(error)
      !! h = stage - bed below DRY_TOLERANCE ⇒ zero momentum.  (Thin-but-
      !! wet cells above DRY_TOLERANCE are desingularized instead — see
      !! test_protect_unstr — so this uses a genuinely dry depth.)
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh)
         ! stage just barely above bed (below DRY_TOLERANCE = 1e-6)
         call prime_state(state, -10.0_wp + 5.0e-7_wp, 1.0_wp, 1.0_wp, -10.0_wp)

         call fused_update_stage1_unstr(state, mesh, 0.1_wp)

         call check(error, state%barotropic%xmom_cv(1) == 0.0_wp, "thin cell xmom zeroed")
         if (allocated(error)) exit checks
         call check(error, state%barotropic%ymom_cv(1) == 0.0_wp, "thin cell ymom zeroed")
         if (allocated(error)) exit checks

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_protect_dry

   subroutine test_protect_negative(error)
      !! Negative h after the Euler update is clamped: stage := bed,
      !! mass_error increases by the deficit.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh)
         call prime_state(state, -9.5_wp, 0.0_wp, 0.0_wp, -10.0_wp)
         ! Negative stage_eu drains the column past bed
         state%barotropic%stage_eu = -1000.0_wp
         state%mass_error = 0.0_wp

         call fused_update_stage1_unstr(state, mesh, 0.1_wp)

         call check(error, all(abs(state%barotropic%stage_cv(1:mesh%ncells_owned) - (-10.0_wp)) < 1.0e-12_wp), &
                    "negative h must reset stage to bed")
         if (allocated(error)) exit checks
         call check(error, state%mass_error > 0.0_wp, "mass_error accumulates the deficit")
         if (allocated(error)) exit checks

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_protect_negative

   subroutine test_stage2_rk2(error)
      !! With zero flux/forcing, stage2 should average current with backup.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh)

         ! After stage1 with non-zero stage_eu the snapshots were stage=2, xmom=0
         ! before update; here we synthesise that scenario directly.
         state%barotropic%bed_cv = -10.0_wp
         state%barotropic%stage_cv = 4.0_wp
         state%barotropic%xmom_cv = 6.0_wp
         state%barotropic%ymom_cv = 0.0_wp
         state%barotropic%stage_backup = 2.0_wp
         state%barotropic%xmom_backup = 0.0_wp
         state%barotropic%ymom_backup = 0.0_wp
         state%barotropic%stage_eu = 0.0_wp
         state%barotropic%xmom_eu = 0.0_wp
         state%barotropic%ymom_eu = 0.0_wp
         state%barotropic%manning_n = 0.0_wp
         state%barotropic%wind_tau_x = 0.0_wp
         state%barotropic%wind_tau_y = 0.0_wp

         call fused_update_stage2_unstr(state, mesh, 0.1_wp)

         call check(error, abs(state%barotropic%stage_cv(1) - 3.0_wp) < 1.0e-12_wp, &
                    "RK2: stage = (4 + 2)/2 = 3")
         if (allocated(error)) exit checks
         call check(error, abs(state%barotropic%xmom_cv(1) - 3.0_wp) < 1.0e-12_wp, &
                    "RK2: xmom = (6 + 0)/2 = 3")
         if (allocated(error)) exit checks

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_stage2_rk2

end module test_kernel_fused_unstr
