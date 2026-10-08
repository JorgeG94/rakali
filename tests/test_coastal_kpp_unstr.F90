!! Unit tests for the unstructured KPP overlay (`rki_ml_kpp_unstr`).
!!
!! Mirrors `test_coastal_kpp` (structured A-grid) on the triangular
!! grid: the KPP kernels fuse the boundary-layer overlay with the
!! explicit two-pass diffusive flux, so the tests observe the effect on
!! hu/hv/hS/hT and the stored boundary-layer depth.  Drives the kernels
!! directly against a 4-triangle test mesh, mapping state on/off device
!! per the usual unstructured enter_data / exit_data lifecycle.
!!
!! Coverage: interior (KPP with no forcing reproduces PP81 bit-for-bit),
!! BL-depth (bounded under stratification, fills an unstratified column),
!! convective activation with no wind, and non-local redistribution.
module test_coastal_kpp_unstr
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
   use rki_ml_kpp_unstr, only: ml_kpp_boundary_layer_unstr, &
                               ml_kpp_viscosity_unstr, &
                               ml_kpp_diffusion_tracers_unstr
   use rki_ml_vertical_viscosity_unstr, only: ml_vertical_viscosity_unstr
   use rki_ml_tracers_unstr, only: ml_vertical_diffusion_tracer_unstr
   implicit none
   private

   public :: collect_coastal_kpp_unstr_tests

   real(wp), parameter :: RHO0 = 1000.0_wp
   real(wp), parameter :: ALPHA_T = 0.17_wp
   real(wp), parameter :: BETA_S = 0.78_wp
   real(wp), parameter :: NU_MAX = 0.01_wp
   real(wp), parameter :: NU_BG = 1.0e-4_wp
   real(wp), parameter :: KAPPA_BG = 1.0e-5_wp
   real(wp), parameter :: RI_ALPHA = 5.0_wp
   real(wp), parameter :: RI_CRIT = 0.3_wp
   real(wp), parameter :: CS_NONLOCAL = 6.3_wp
   real(wp), parameter :: C_VT2 = 1.8_wp

contains

   subroutine collect_coastal_kpp_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("kpp_unstr_no_forcing_equals_pp81", test_no_forcing_equals_pp81), &
                  new_unittest("kpp_unstr_stratified_bl_in_column", test_strat_bl_depth), &
                  new_unittest("kpp_unstr_unstratified_fills_column", test_unstrat_fills_column), &
                  new_unittest("kpp_unstr_convective_activates_no_wind", test_convective_no_wind), &
                  new_unittest("kpp_unstr_nonlocal_redistributes_conserves", test_nonlocal_redistributes) &
                  ]
   end subroutine collect_coastal_kpp_unstr_tests

   ! -----------------------------------------------------------------
   ! Helpers
   ! -----------------------------------------------------------------

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

   subroutine set_physics(state, use_kpp)
      !! Common EOS / vmix / KPP scalar parameters.
      type(unstr_state_t), intent(inout) :: state
      logical, intent(in) :: use_kpp
      state%multilayer%rho_0 = RHO0
      state%multilayer%alpha_T = ALPHA_T
      state%multilayer%beta_S = BETA_S
      state%multilayer%vmix_nu_max = NU_MAX
      state%multilayer%vmix_nu_bg = NU_BG
      state%multilayer%vmix_ri_alpha = RI_ALPHA
      state%multilayer%kappa_S_bg = KAPPA_BG
      state%multilayer%kappa_T_bg = KAPPA_BG
      state%multilayer%kpp_ri_crit = RI_CRIT
      state%multilayer%kpp_cs_nonlocal = CS_NONLOCAL
      state%multilayer%kpp_c_vt2 = C_VT2
      state%multilayer%vmix_use_kpp = use_kpp
   end subroutine set_physics

   subroutine set_depth(state, mesh, nz, h_layer)
      !! Uniform layer thickness; barotropic depth = sum of layers.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer
      integer :: c, k
      real(wp) :: h_total
      h_total = real(nz, wp)*h_layer
      state%barotropic%bed_cv = -h_total
      state%barotropic%stage_cv = 0.0_wp
      state%barotropic%height_cv = h_total
      do c = 1, mesh%ncells_total
         do k = 1, nz
            state%multilayer%h_layer(k, c) = h_layer
         end do
      end do
   end subroutine set_depth

   subroutine run_kpp_unstr(mesh, state, dt)
      type(mesh_unstr_t), intent(inout) :: mesh
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(in) :: dt
      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call ml_kpp_boundary_layer_unstr(state, &
                                       state%multilayer%tracers(state%multilayer%idx_salinity)%hTr)
      call ml_kpp_viscosity_unstr(state, dt)
      call ml_kpp_diffusion_tracers_unstr(state, &
                                          state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                          state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                          state%multilayer%kappa_S_bg, state%multilayer%kappa_T_bg, dt)
      !$acc update self(state%multilayer%kpp_bl_depth, &
      !$acc&            state%multilayer%kpp_u_star, &
      !$acc&            state%multilayer%kpp_buoy_flux)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_kpp_unstr

   subroutine run_pp81_unstr(mesh, state, dt)
      type(mesh_unstr_t), intent(inout) :: mesh
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(in) :: dt
      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call ml_vertical_viscosity_unstr(state, dt)
      call ml_vertical_diffusion_tracer_unstr(state, &
                                              state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                              state%multilayer%kappa_S_bg, dt)
      call ml_vertical_diffusion_tracer_unstr(state, &
                                              state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                              state%multilayer%kappa_T_bg, dt)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_pp81_unstr

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_no_forcing_equals_pp81(error)
      !! KPP with zero wind + zero surface flux ⇒ overlay inactive and
      !! non-local skipped, so the momentum + tracer fluxes reproduce the
      !! PP81 kernels bit-for-bit on a sheared, stratified column.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state_a, state_b
      integer, parameter :: NZ = 5
      real(wp), parameter :: H_LAYER = 10.0_wp, DT = 1.0_wp, TOL = 1.0e-12_wp
      integer :: ncw, idx_s, idx_t
      checks: block
         call build_test_mesh(mesh)
         call state_unstr_init(state_a, mesh, nz_ml=NZ)
         call state_unstr_init(state_b, mesh, nz_ml=NZ)
         call set_physics(state_a, .false.)
         call set_physics(state_b, .true.)
         call setup_strat_sheared(state_a, mesh, NZ, H_LAYER)
         call setup_strat_sheared(state_b, mesh, NZ, H_LAYER)

         call run_pp81_unstr(mesh, state_a, DT)
         call run_kpp_unstr(mesh, state_b, DT)

         ncw = mesh%ncells_owned
         idx_s = state_a%multilayer%idx_salinity
         idx_t = state_a%multilayer%idx_temperature
         call check(error, maxval(abs(state_b%multilayer%hu_layer(:, 1:ncw) - &
                                      state_a%multilayer%hu_layer(:, 1:ncw))) < TOL, &
                    "KPP hu /= PP81 hu (no forcing)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(state_b%multilayer%hv_layer(:, 1:ncw) - &
                                      state_a%multilayer%hv_layer(:, 1:ncw))) < TOL, &
                    "KPP hv /= PP81 hv (no forcing)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(state_b%multilayer%tracers(idx_s)%hTr(:, 1:ncw) - &
                                      state_a%multilayer%tracers(idx_s)%hTr(:, 1:ncw))) < TOL, &
                    "KPP hS /= PP81 hS (no forcing)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(state_b%multilayer%tracers(idx_t)%hTr(:, 1:ncw) - &
                                      state_a%multilayer%tracers(idx_t)%hTr(:, 1:ncw))) < TOL, &
                    "KPP hT /= PP81 hT (no forcing)")
      end block checks
      call state_unstr_destroy(state_a)
      call state_unstr_destroy(state_b)
      call mesh_unstr_destroy(mesh)
   end subroutine test_no_forcing_equals_pp81

   subroutine setup_strat_sheared(state, mesh, nz, h_layer)
      !! Linearly stratified (densest at bed) + vertically sheared column.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer
      integer :: c, k, idx_s, idx_t
      call set_depth(state, mesh, nz, h_layer)
      idx_s = state%multilayer%idx_salinity
      idx_t = state%multilayer%idx_temperature
      do c = 1, mesh%ncells_total
         do k = 1, nz
            state%multilayer%hu_layer(k, c) = h_layer*0.1_wp*real(k - 1, wp)
            state%multilayer%hv_layer(k, c) = 0.0_wp
            state%multilayer%rho_layer(k, c) = 1030.0_wp - 0.5_wp*real(k - 1, wp)
            state%multilayer%tracers(idx_s)%hTr(k, c) = h_layer*(35.0_wp - 0.1_wp*real(k - 1, wp))
            state%multilayer%tracers(idx_t)%hTr(k, c) = h_layer*(15.0_wp + 0.2_wp*real(k - 1, wp))
         end do
      end do
   end subroutine setup_strat_sheared

   subroutine test_strat_bl_depth(error)
      !! Mixed-layer-over-thermocline + wind: BL must extend below the
      !! surface layer and stop short of the bed.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer, parameter :: NZ = 8, NZ_ML = 3
      real(wp), parameter :: H_LAYER = 25.0_wp, DT = 1.0_wp, TAU = 0.1_wp, U_ML = 0.5_wp
      integer :: c, k, idx_s, idx_t
      real(wp) :: h_b
      checks: block
         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)
         call set_physics(state, .true.)
         call set_depth(state, mesh, NZ, H_LAYER)
         idx_s = state%multilayer%idx_salinity
         idx_t = state%multilayer%idx_temperature
         do c = 1, mesh%ncells_total
            do k = 1, NZ
               state%multilayer%hv_layer(k, c) = 0.0_wp
               state%multilayer%tracers(idx_s)%hTr(k, c) = H_LAYER*35.0_wp
               state%multilayer%tracers(idx_t)%hTr(k, c) = H_LAYER*15.0_wp
               if (k > NZ - NZ_ML) then
                  state%multilayer%hu_layer(k, c) = H_LAYER*U_ML
                  state%multilayer%rho_layer(k, c) = 1027.0_wp
               else
                  state%multilayer%hu_layer(k, c) = 0.0_wp
                  state%multilayer%rho_layer(k, c) = 1030.0_wp - 0.5_wp*real(k - 1, wp)
               end if
            end do
         end do
         state%barotropic%wind_tau_x = TAU

         call run_kpp_unstr(mesh, state, DT)

         h_b = state%multilayer%kpp_bl_depth(1)
         call check(error, h_b > H_LAYER, "stratified BL collapsed to a sub-layer")
         if (allocated(error)) exit checks
         call check(error, h_b < real(NZ, wp)*H_LAYER, "stratified BL reached the bed")
      end block checks
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_strat_bl_depth

   subroutine test_unstrat_fills_column(error)
      !! Uniform density + wind: bulk Ri never crosses ⇒ BL = total H.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer, parameter :: NZ = 4
      real(wp), parameter :: H_LAYER = 20.0_wp, DT = 1.0_wp, TAU = 0.05_wp, TOL = 1.0e-9_wp
      integer :: c, k, idx_s, idx_t
      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)
      call set_physics(state, .true.)
      call set_depth(state, mesh, NZ, H_LAYER)
      idx_s = state%multilayer%idx_salinity
      idx_t = state%multilayer%idx_temperature
      do c = 1, mesh%ncells_total
         do k = 1, NZ
            state%multilayer%hu_layer(k, c) = 0.0_wp
            state%multilayer%hv_layer(k, c) = 0.0_wp
            state%multilayer%rho_layer(k, c) = 1030.0_wp
            state%multilayer%tracers(idx_s)%hTr(k, c) = H_LAYER*35.0_wp
            state%multilayer%tracers(idx_t)%hTr(k, c) = H_LAYER*15.0_wp
         end do
      end do
      state%barotropic%wind_tau_x = TAU

      call run_kpp_unstr(mesh, state, DT)

      call check(error, abs(state%multilayer%kpp_bl_depth(1) - real(NZ, wp)*H_LAYER) < TOL, &
                 "unstratified BL did not fill the column")
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_unstrat_fills_column

   subroutine test_convective_no_wind(error)
      !! Surface cooling, NO wind, neutral column: the convective scale
      !! w_* alone (u_* = 0) activates the overlay and mixes the sheared
      !! momentum more than the quiescent PP81 baseline.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state_a, state_b
      integer, parameter :: NZ = 6
      real(wp), parameter :: H_LAYER = 10.0_wp, DT = 30.0_wp, Q_COOL = -300.0_wp
      integer :: c, k, idx_s, idx_t
      real(wp) :: spread_a, spread_b
      real(wp) :: u_a(NZ), u_b(NZ)
      checks: block
         call build_test_mesh(mesh)
         call state_unstr_init(state_a, mesh, nz_ml=NZ)
         call state_unstr_init(state_b, mesh, nz_ml=NZ)
         call set_physics(state_a, .false.)
         call set_physics(state_b, .true.)
         call setup_neutral_sheared(state_a, mesh, NZ, H_LAYER)
         call setup_neutral_sheared(state_b, mesh, NZ, H_LAYER)
         ! PP81 reference (no cooling — PP81 ignores it anyway).
         call run_pp81_unstr(mesh, state_a, DT)
         ! KPP with cooling, no wind.
         do c = 1, mesh%ncells_total
            state_b%forcing%heat_flux(c) = Q_COOL
         end do
         call run_kpp_unstr(mesh, state_b, DT)

         idx_s = state_a%multilayer%idx_salinity
         idx_t = state_a%multilayer%idx_temperature
         do k = 1, NZ
            u_a(k) = state_a%multilayer%hu_layer(k, 1)/state_a%multilayer%h_layer(k, 1)
            u_b(k) = state_b%multilayer%hu_layer(k, 1)/state_b%multilayer%h_layer(k, 1)
         end do
         spread_a = maxval(u_a) - minval(u_a)
         spread_b = maxval(u_b) - minval(u_b)

         call check(error, state_b%multilayer%kpp_buoy_flux(1) < 0.0_wp, &
                    "cooling did not produce B_0 < 0")
         if (allocated(error)) exit checks
         call check(error, spread_b < spread_a, &
                    "convective overlay did not mix momentum without wind")
      end block checks
      call state_unstr_destroy(state_a)
      call state_unstr_destroy(state_b)
      call mesh_unstr_destroy(mesh)
   end subroutine test_convective_no_wind

   subroutine setup_neutral_sheared(state, mesh, nz, h_layer)
      !! Neutral (uniform-density) but vertically sheared column.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer
      integer :: c, k, idx_s, idx_t
      call set_depth(state, mesh, nz, h_layer)
      idx_s = state%multilayer%idx_salinity
      idx_t = state%multilayer%idx_temperature
      do c = 1, mesh%ncells_total
         do k = 1, nz
            state%multilayer%hu_layer(k, c) = h_layer*0.05_wp*real(k - 1, wp)
            state%multilayer%hv_layer(k, c) = 0.0_wp
            state%multilayer%rho_layer(k, c) = 1025.0_wp
            state%multilayer%tracers(idx_s)%hTr(k, c) = h_layer*35.0_wp
            state%multilayer%tracers(idx_t)%hTr(k, c) = h_layer*15.0_wp
         end do
      end do
   end subroutine setup_neutral_sheared

   subroutine test_nonlocal_redistributes(error)
      !! Destabilizing flux (cooling) on a uniform-temperature column:
      !! the diffusive flux is zero (no gradient), so any change in hT is
      !! the KPP non-local term.  It must make the profile non-uniform and
      !! conserve the column-integrated heat.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer, parameter :: NZ = 6
      real(wp), parameter :: H_LAYER = 10.0_wp, DT = 60.0_wp, Q_COOL = -300.0_wp
      integer :: c, k, idx_s, idx_t
      real(wp) :: col0, col1, spread
      real(wp) :: hT0(NZ)
      checks: block
         call build_test_mesh(mesh)
         call state_unstr_init(state, mesh, nz_ml=NZ)
         call set_physics(state, .true.)
         call set_depth(state, mesh, NZ, H_LAYER)
         idx_s = state%multilayer%idx_salinity
         idx_t = state%multilayer%idx_temperature
         do c = 1, mesh%ncells_total
            do k = 1, NZ
               state%multilayer%hu_layer(k, c) = 0.0_wp
               state%multilayer%hv_layer(k, c) = 0.0_wp
               state%multilayer%rho_layer(k, c) = 1025.0_wp
               state%multilayer%tracers(idx_s)%hTr(k, c) = H_LAYER*35.0_wp
               state%multilayer%tracers(idx_t)%hTr(k, c) = H_LAYER*15.0_wp
            end do
            state%forcing%heat_flux(c) = Q_COOL
         end do
         do k = 1, NZ
            hT0(k) = state%multilayer%tracers(idx_t)%hTr(k, 1)
         end do

         call run_kpp_unstr(mesh, state, DT)

         spread = maxval(state%multilayer%tracers(idx_t)%hTr(:, 1)) - &
                  minval(state%multilayer%tracers(idx_t)%hTr(:, 1))
         col0 = sum(hT0)
         col1 = sum(state%multilayer%tracers(idx_t)%hTr(:, 1))

         call check(error, spread > 1.0e-6_wp, &
                    "non-local transport left a uniform column unchanged")
         if (allocated(error)) exit checks
         call check(error, abs(col1 - col0) < 1.0e-9_wp*max(abs(col0), 1.0_wp), &
                    "non-local transport did not conserve column heat")
      end block checks
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_nonlocal_redistributes

end module test_coastal_kpp_unstr
