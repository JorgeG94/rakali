!! Kato-Phillips-style validation for the unstructured k-epsilon closure.
!!
!! Drives the keps vmix kernel sequence in the exact solver-dispatch order
!! (friction velocities -> k/eps update with wall functions -> implicit
!! momentum apply -> implicit tracer apply) on a stratified, sheared column
!! with surface wind, on the 4-triangle test mesh.  This is the only place the
!! keps-ON path runs end-to-end on device; it confirms the closure activates
!! turbulence, mixes the shear, conserves column momentum + tracer content, and
!! stays bounded.  Mirrors the structured `test_kato_phillips` intent and the
!! `test_coastal_kpp_unstr` mesh/lifecycle pattern.
module test_keps_kato_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, KEPS_STAB_CONSTANT
   use rki_mesh_unstr, only: mesh_unstr_t, &
                             mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, &
                             mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, &
                              state_unstr_init, state_unstr_destroy, &
                              state_unstr_enter_data, state_unstr_exit_data
   use rki_ml_keps_unstr, only: ml_keps_friction_velocities_unstr, &
                                ml_keps_update_unstr, &
                                ml_keps_apply_momentum_unstr, &
                                ml_keps_apply_diffusion_tracer_unstr
   implicit none
   private

   public :: collect_keps_kato_unstr_tests

   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: KAPPA_BG = 1.0e-5_wp

contains

   subroutine collect_keps_kato_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("keps_kato_unstr_wind_mixing", test_wind_mixing), &
                  new_unittest("keps_kato_unstr_conserves", test_conservation) &
                  ]
   end subroutine collect_keps_kato_unstr_tests

   subroutine build_test_mesh(mesh)
      type(mesh_unstr_t), intent(out) :: mesh
      integer, parameter :: NV = 5, NTRI = 4
      real(wp) :: vx(NV), vy(NV)
      integer :: tv(3, NTRI)
      vx = [0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.5_wp]
      vy = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.5_wp]
      tv(:, 1) = [1, 2, 5]
      tv(:, 2) = [2, 3, 5]
      tv(:, 3) = [3, 4, 5]
      tv(:, 4) = [4, 1, 5]
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_test_mesh

   subroutine setup_keps_column(state, mesh, nz, h_layer, tau)
      !! keps ON, stably stratified (densest at bed) + vertically sheared
      !! column, with a constant surface wind stress.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer, tau
      integer :: c, k, idx_s, idx_t
      real(wp) :: h_total

      ! keps fields are allocated unconditionally, so the knob can be flipped
      ! after init (mirrors test_coastal_kpp_unstr enabling KPP post-init).
      state%multilayer%vmix_use_keps = .true.
      state%multilayer%keps_stability = KEPS_STAB_CONSTANT
      state%multilayer%rho_0 = RHO0
      state%multilayer%kappa_S_bg = KAPPA_BG
      state%multilayer%kappa_T_bg = KAPPA_BG
      state%barotropic%manning_n = 0.0_wp
      state%barotropic%wind_tau_x = tau
      state%barotropic%wind_tau_y = 0.0_wp

      h_total = real(nz, wp)*h_layer
      state%barotropic%bed_cv = -h_total
      state%barotropic%stage_cv = 0.0_wp
      state%barotropic%height_cv = h_total
      idx_s = state%multilayer%idx_salinity
      idx_t = state%multilayer%idx_temperature
      do c = 1, mesh%ncells_total
         do k = 1, nz
            state%multilayer%h_layer(k, c) = h_layer
            state%multilayer%hu_layer(k, c) = h_layer*0.1_wp*real(k - 1, wp)   ! shear
            state%multilayer%hv_layer(k, c) = 0.0_wp
            state%multilayer%rho_layer(k, c) = 1030.0_wp - 0.5_wp*real(k - 1, wp)   ! dense at bed
            state%multilayer%tracers(idx_s)%hTr(k, c) = h_layer*(35.0_wp - 0.1_wp*real(k - 1, wp))
            state%multilayer%tracers(idx_t)%hTr(k, c) = h_layer*(15.0_wp + 0.2_wp*real(k - 1, wp))
         end do
      end do
   end subroutine setup_keps_column

   subroutine run_keps_unstr(mesh, state, dt, nsteps)
      !! The keps vmix sequence in dispatch order, nsteps times, inside one
      !! device data region (data stays resident across steps, as in the solver).
      type(mesh_unstr_t), intent(inout) :: mesh
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(in) :: dt
      integer, intent(in) :: nsteps
      integer :: step, idx_s, idx_t

      idx_s = state%multilayer%idx_salinity
      idx_t = state%multilayer%idx_temperature
      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      do step = 1, nsteps
         call ml_keps_friction_velocities_unstr(state)
         call ml_keps_update_unstr(state, dt, .true.)
         call ml_keps_apply_momentum_unstr(state, dt)
         call ml_keps_apply_diffusion_tracer_unstr(state, &
                                                   state%multilayer%tracers(idx_s)%hTr, KAPPA_BG, dt)
         call ml_keps_apply_diffusion_tracer_unstr(state, &
                                                   state%multilayer%tracers(idx_t)%hTr, KAPPA_BG, dt)
      end do
      ! keps_nut/tke are delete-mapped (not copyout) -- pull them back for the
      ! host-side assertions before exit_data tears the device copy down.
      !$omp target update from(state%multilayer%keps_nut, state%multilayer%keps_tke)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_keps_unstr

   subroutine test_wind_mixing(error)
      !! Wind on a stratified, sheared column: keps activates turbulence (nu_t
      !! grows well above the uniform seed), mixes the shear toward uniform, and
      !! stays bounded (clamped to [NUT_MIN, NUT_MAX]).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer, parameter :: NZ = 8
      real(wp), parameter :: H_LAYER = 20.0_wp, DT = 10.0_wp, TAU = 0.1_wp
      integer :: k
      real(wp) :: spread0, spread1
      real(wp) :: u_k(NZ)
      checks: block
         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)
         call setup_keps_column(state, mesh, NZ, H_LAYER, TAU)

         do k = 1, NZ
            u_k(k) = state%multilayer%hu_layer(k, 1)/state%multilayer%h_layer(k, 1)
         end do
         spread0 = maxval(u_k) - minval(u_k)

         call run_keps_unstr(mesh, state, DT, 30)

         ! Bounded + positive everywhere (catches NaN/blow-up: NaN fails both).
         call check(error, all(state%multilayer%keps_nut > 0.0_wp) .and. &
                    all(state%multilayer%keps_nut <= 10.0_wp), &
                    "keps nu_t stays positive and bounded (no blow-up)")
         if (allocated(error)) exit checks
         ! Turbulence activated: the surface wall + shear grew nu_t above the seed.
         call check(error, maxval(state%multilayer%keps_nut(:, 1)) > 5.0e-4_wp, &
                    "wind + shear activated keps turbulence (nu_t grew above seed)")
         if (allocated(error)) exit checks
         ! The vertical shear was mixed down.
         do k = 1, NZ
            u_k(k) = state%multilayer%hu_layer(k, 1)/state%multilayer%h_layer(k, 1)
         end do
         spread1 = maxval(u_k) - minval(u_k)
         call check(error, spread1 < spread0, &
                    "keps viscosity mixed the vertical shear")
      end block checks
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_wind_mixing

   subroutine test_conservation(error)
      !! The keps apply kernels (no-flux Thomas) conserve column momentum and
      !! tracer content to round-off through the full sequence.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer, parameter :: NZ = 6
      real(wp), parameter :: H_LAYER = 15.0_wp, DT = 20.0_wp, TAU = 0.05_wp
      integer :: idx_s
      real(wp) :: sum_u0, sum_s0, sum_u1, sum_s1
      checks: block
         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)
         call setup_keps_column(state, mesh, NZ, H_LAYER, TAU)
         idx_s = state%multilayer%idx_salinity

         sum_u0 = sum(state%multilayer%hu_layer(:, 1))
         sum_s0 = sum(state%multilayer%tracers(idx_s)%hTr(:, 1))

         call run_keps_unstr(mesh, state, DT, 20)

         sum_u1 = sum(state%multilayer%hu_layer(:, 1))
         sum_s1 = sum(state%multilayer%tracers(idx_s)%hTr(:, 1))

         call check(error, abs(sum_u1 - sum_u0) < 1.0e-9_wp*(abs(sum_u0) + 1.0_wp), &
                    "column x-momentum conserved through the keps sequence")
         if (allocated(error)) exit checks
         call check(error, abs(sum_s1 - sum_s0) < 1.0e-9_wp*(abs(sum_s0) + 1.0_wp), &
                    "column salinity content conserved through the keps sequence")
      end block checks
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_conservation

end module test_keps_kato_unstr
