!! Unit tests for the unstructured multilayer dynamics kernels
!! (rki_ml_dynamics_unstr).  Drives the public lifecycle pieces directly
!! against a 4-triangle test mesh.
!!
!! Every kernel here uses `!$acc parallel loop present(state)` (or
!! similar), so each test must call `state_unstr_enter_data` after host
!! IC setup and `state_unstr_exit_data` before reading the result from
!! the host arrays.  Mirrors the pattern in `test_multilayer_unstr`.
module test_ml_dynamics_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t, &
                             mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, &
                             mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, &
                              state_unstr_init, state_unstr_destroy, &
                              state_unstr_enter_data, state_unstr_exit_data
   use rki_ml_dynamics_unstr, only: ml_init_layers_unstr, &
                                    ml_save_barotropic_pre_unstr, &
                                    ml_distribute_barotropic_unstr, &
                                    ml_compute_barotropic_momentum_unstr, &
                                    ml_compute_total_depth_unstr, &
                                    ml_save_layer_state_unstr, &
                                    ml_rk2_average_layers_unstr, &
                                    ml_compute_vertical_velocity_unstr, &
                                    ml_vertical_exchange_unstr
   implicit none
   private

   public :: collect_ml_dynamics_unstr_tests

   integer, parameter :: NZ = 3

contains

   subroutine collect_ml_dynamics_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("init_layers_sigma", test_init_layers), &
                  new_unittest("save_barotropic_pre", test_save_baro_pre), &
                  new_unittest("distribute_barotropic_uniform_accel", test_distribute), &
                  new_unittest("compute_barotropic_momentum_sums_layers", test_recompute_baro), &
                  new_unittest("compute_total_depth_sums_layers", test_total_depth), &
                  new_unittest("save_layer_state_snapshots_q0", test_save_layer), &
                  new_unittest("rk2_average_layers_blends_q_q0", test_rk2_avg), &
                  new_unittest("compute_vertical_velocity_zero_div", test_vert_vel_zero), &
                  new_unittest("compute_vertical_velocity_nonzero_div", test_vert_vel_nonzero), &
                  new_unittest("vertical_exchange_runs", test_vert_exchange) &
                  ]
   end subroutine collect_ml_dynamics_unstr_tests

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

   subroutine init_constant_column(state, mesh, H, U, V)
      !! Set up a uniform column on host: height H, depth-averaged velocity (U, V).
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: H, U, V
      integer :: k

      state%barotropic%bed_cv = -10.0_wp
      state%barotropic%stage_cv = state%barotropic%bed_cv + H
      state%barotropic%height_cv = H
      state%barotropic%xmom_cv = U*H
      state%barotropic%ymom_cv = V*H
      do k = 1, mesh%ncells_total
         state%multilayer%h_layer(:, k) = H/real(NZ, wp)
         state%multilayer%hu_layer(:, k) = U*H/real(NZ, wp)
         state%multilayer%hv_layer(:, k) = V*H/real(NZ, wp)
      end do
   end subroutine init_constant_column

   subroutine map_in(mesh, state)
      type(mesh_unstr_t), intent(inout) :: mesh
      type(unstr_state_t), intent(inout) :: state
      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
   end subroutine map_in

   subroutine map_out(mesh, state)
      type(mesh_unstr_t), intent(inout) :: mesh
      type(unstr_state_t), intent(inout) :: state
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine map_out

   ! ----------------------------------------------------------------------
   ! Tests
   ! ----------------------------------------------------------------------

   subroutine test_init_layers(error)
      !! Pure-sigma init should split H into nz equal layers.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)

         state%barotropic%bed_cv = -5.0_wp
         state%barotropic%stage_cv = 5.0_wp     ! H = 10
         state%barotropic%xmom_cv = 3.0_wp     ! HU
         state%barotropic%ymom_cv = 6.0_wp     ! HV

         call map_in(mesh, state)
         call ml_init_layers_unstr(state)
         call map_out(mesh, state)

         do k = 1, NZ
            call check(error, abs(state%multilayer%h_layer(k, 1) - 10.0_wp/real(NZ, wp)) < 1.0e-12_wp, &
                       "uniform sigma layer thickness")
            if (allocated(error)) exit checks
            call check(error, abs(state%multilayer%hu_layer(k, 1) - 1.0_wp) < 1.0e-12_wp, &
                       "uniform sigma hu_k = dsig*HU")
            if (allocated(error)) exit checks
            call check(error, abs(state%multilayer%hv_layer(k, 1) - 2.0_wp) < 1.0e-12_wp, &
                       "uniform sigma hv_k = dsig*HV")
            if (allocated(error)) exit checks
         end do

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_init_layers

   subroutine test_save_baro_pre(error)
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)
         call init_constant_column(state, mesh, 10.0_wp, 0.5_wp, -0.25_wp)

         call map_in(mesh, state)
         call ml_save_barotropic_pre_unstr(state)
         ! xmom_pre / ymom_pre are mapped `create` (workspace), so the kernel
         ! writes them on the device but exit_data won't copy them back.
         ! [acc->omp-cpu] dropped device directive (host no-op): update
         call map_out(mesh, state)

         call check(error, all(abs(state%multilayer%xmom_pre - state%barotropic%xmom_cv) < 1.0e-12_wp), &
                    "xmom_pre snapshot equals xmom_cv")
         if (allocated(error)) exit checks
         call check(error, all(abs(state%multilayer%ymom_pre - state%barotropic%ymom_cv) < 1.0e-12_wp), &
                    "ymom_pre snapshot equals ymom_cv")
         if (allocated(error)) exit checks

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_save_baro_pre

   subroutine test_distribute(error)
      !! After save_pre + bumping HU, the distribute kernel should add the
      !! same du = dHU/H to every layer's velocity.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp), parameter :: H = 10.0_wp, U0 = 0.5_wp, V0 = 0.0_wp
      real(wp), parameter :: DHU = 2.0_wp
      real(wp) :: u_k, expected_u
      integer :: k, c
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)
         call init_constant_column(state, mesh, H, U0, V0)

         ! Pre-fill xmom_pre / ymom_pre with the pre-update barotropic momentum,
         ! and bump xmom_cv to simulate the post-HLL state.  xmom_pre is
         ! `create`-mapped on enter_data so we must push the host values to the
         ! device explicitly with `update device`.
         state%multilayer%xmom_pre = state%barotropic%xmom_cv
         state%multilayer%ymom_pre = state%barotropic%ymom_cv
         state%barotropic%xmom_cv = state%barotropic%xmom_cv + DHU

         call map_in(mesh, state)
         ! [acc->omp-cpu] dropped device directive (host no-op): update
         call ml_distribute_barotropic_unstr(state)
         call map_out(mesh, state)

         ! Each layer's velocity should equal (U0 + DHU/H) since stage was unchanged
         expected_u = U0 + DHU/H
         do c = 1, mesh%ncells_total
            do k = 1, NZ
               u_k = state%multilayer%hu_layer(k, c)/state%multilayer%h_layer(k, c)
               call check(error, abs(u_k - expected_u) < 1.0e-12_wp, &
                          "every layer should pick up the same du = dHU/H")
               if (allocated(error)) exit checks
            end do
         end do

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_distribute

   subroutine test_recompute_baro(error)
      !! HU_baro = sum_k hu_k after recompute.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: c, k
      real(wp) :: expected_HU, expected_HV
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)
         call init_constant_column(state, mesh, 10.0_wp, 0.0_wp, 0.0_wp)

         ! Inject a non-uniform momentum profile so the sum is non-trivial
         do c = 1, mesh%ncells_total
            do k = 1, NZ
               state%multilayer%hu_layer(k, c) = 0.1_wp*real(k, wp)
               state%multilayer%hv_layer(k, c) = -0.05_wp*real(k, wp)
            end do
         end do

         call map_in(mesh, state)
         call ml_compute_barotropic_momentum_unstr(state)
         call map_out(mesh, state)

         expected_HU = 0.0_wp
         expected_HV = 0.0_wp
         do k = 1, NZ
            expected_HU = expected_HU + 0.1_wp*real(k, wp)
            expected_HV = expected_HV - 0.05_wp*real(k, wp)
         end do

         do c = 1, mesh%ncells_owned
            call check(error, abs(state%barotropic%xmom_cv(c) - expected_HU) < 1.0e-12_wp, &
                       "xmom_cv = sum hu_layer")
            if (allocated(error)) exit checks
            call check(error, abs(state%barotropic%ymom_cv(c) - expected_HV) < 1.0e-12_wp, &
                       "ymom_cv = sum hv_layer")
            if (allocated(error)) exit checks
         end do

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_recompute_baro

   subroutine test_total_depth(error)
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: c, k
      real(wp) :: expected_H
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)

         state%barotropic%bed_cv = -7.0_wp
         do c = 1, mesh%ncells_total
            do k = 1, NZ
               state%multilayer%h_layer(k, c) = real(k, wp)   ! 1, 2, 3
            end do
         end do

         call map_in(mesh, state)
         call ml_compute_total_depth_unstr(state)
         call map_out(mesh, state)

         expected_H = 1.0_wp + 2.0_wp + 3.0_wp

         do c = 1, mesh%ncells_owned
            call check(error, abs(state%barotropic%height_cv(c) - expected_H) < 1.0e-12_wp, &
                       "height_cv = sum h_layer")
            if (allocated(error)) exit checks
            call check(error, abs(state%barotropic%stage_cv(c) - (-7.0_wp + expected_H)) < 1.0e-12_wp, &
                       "stage_cv = bed + sum h")
            if (allocated(error)) exit checks
         end do

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_total_depth

   subroutine test_save_layer(error)
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: c, k
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)

         do c = 1, mesh%ncells_total
            do k = 1, NZ
               state%multilayer%h_layer(k, c) = 0.5_wp*real(k, wp)
               state%multilayer%hu_layer(k, c) = 0.1_wp*real(k, wp)
               state%multilayer%hv_layer(k, c) = -0.05_wp*real(k, wp)
            end do
         end do

         call map_in(mesh, state)
         call ml_save_layer_state_unstr(state)
         ! h_layer0 / hu_layer0 / hv_layer0 are workspace (delete-mapped on exit),
         ! so explicitly pull them back before the host-side assertion.
         ! [acc->omp-cpu] dropped device directive (host no-op): update
         call map_out(mesh, state)

         call check(error, all(abs(state%multilayer%h_layer - state%multilayer%h_layer0) < 1.0e-12_wp), &
                    "h_layer0 = h_layer after save")
         if (allocated(error)) exit checks
         call check(error, all(abs(state%multilayer%hu_layer - state%multilayer%hu_layer0) < 1.0e-12_wp), &
                    "hu_layer0 = hu_layer after save")
         if (allocated(error)) exit checks

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_save_layer

   subroutine test_rk2_avg(error)
      !! After RK2 average: q := 0.5 * (q + q0)
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: c, k
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)

         do c = 1, mesh%ncells_total
            do k = 1, NZ
               state%multilayer%h_layer(k, c) = 4.0_wp
               state%multilayer%h_layer0(k, c) = 2.0_wp
               state%multilayer%hu_layer(k, c) = 6.0_wp
               state%multilayer%hu_layer0(k, c) = 0.0_wp
            end do
         end do

         call map_in(mesh, state)
         call ml_rk2_average_layers_unstr(state)
         call map_out(mesh, state)

         call check(error, all(abs(state%multilayer%h_layer - 3.0_wp) < 1.0e-12_wp), &
                    "h_layer averaged to 3")
         if (allocated(error)) exit checks
         call check(error, all(abs(state%multilayer%hu_layer - 3.0_wp) < 1.0e-12_wp), &
                    "hu_layer averaged to 3")
         if (allocated(error)) exit checks

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_rk2_avg

   subroutine test_vert_vel_zero(error)
      !! Zero mass flux ⇒ zero divergence ⇒ w stays zero.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)
      call init_constant_column(state, mesh, 10.0_wp, 0.0_wp, 0.0_wp)
      state%barotropic%mass_flux_e = 0.0_wp

      call map_in(mesh, state)
      ! mass_flux_e is `create`-mapped — push the host zero into the device buffer.
      ! [acc->omp-cpu] dropped device directive (host no-op): update
      call ml_compute_vertical_velocity_unstr(state, mesh)
      call map_out(mesh, state)

      call check(error, all(abs(state%multilayer%w_interface) < 1.0e-12_wp), &
                 "zero mass flux must give zero w at every interface")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_vert_vel_zero

   subroutine test_vert_vel_nonzero(error)
      !! Non-zero edge flux must produce non-zero w in upper interfaces
      !! and a strictly-zero w at the bottom (no penetration).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: c
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)
         call init_constant_column(state, mesh, 10.0_wp, 0.0_wp, 0.0_wp)

         ! Inject a non-trivial divergence pattern on the owned cells
         do c = 1, mesh%ncells_owned
            state%barotropic%mass_flux_e(:, c) = 0.5_wp
         end do

         call map_in(mesh, state)
         ! mass_flux_e is `create`-mapped on enter_data — push the host pattern
         ! to the device buffer so the kernel actually sees a non-zero divergence.
         ! [acc->omp-cpu] dropped device directive (host no-op): update
         call ml_compute_vertical_velocity_unstr(state, mesh)
         call map_out(mesh, state)

         do c = 1, mesh%ncells_owned
            call check(error, abs(state%multilayer%w_interface(1, c)) < 1.0e-12_wp, &
                       "w(1) must be zero (no penetration at bottom)")
            if (allocated(error)) exit checks
            call check(error, abs(state%multilayer%w_interface(NZ + 1, c)) > 1.0e-12_wp, &
                       "w(nz+1) must be non-zero with non-zero flux")
            if (allocated(error)) exit checks
         end do

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_vert_vel_nonzero

   subroutine test_vert_exchange(error)
      !! Just confirm the exchange runs without error and respects the
      !! "use_multilayer = .false. -> no-op" guard.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: hu_before(NZ, 4)
      integer :: c, k
      checks: block

         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)
         call init_constant_column(state, mesh, 10.0_wp, 0.0_wp, 0.0_wp)
         state%multilayer%w_interface = 0.0_wp

         ! Snapshot momentum
         do c = 1, mesh%ncells_total
            do k = 1, NZ
               hu_before(k, c) = state%multilayer%hu_layer(k, c)
            end do
         end do

         call map_in(mesh, state)
         call ml_vertical_exchange_unstr(state, 0.1_wp)
         call map_out(mesh, state)

         ! With w=0, there is no transport, so momentum is unchanged
         do c = 1, mesh%ncells_owned
            do k = 1, NZ
               call check(error, abs(state%multilayer%hu_layer(k, c) - hu_before(k, c)) < 1.0e-12_wp, &
                          "zero w must leave layer momentum unchanged")
               if (allocated(error)) exit checks
            end do
         end do

      end block checks
      continue
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_vert_exchange

end module test_ml_dynamics_unstr
