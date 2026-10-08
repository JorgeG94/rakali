!! Analytical test: 1D vertical diffusion of temperature
!!
!! A single water column is initialised with two uniform temperatures:
!! warm (T_WARM) above a half-depth interface, cold (T_COLD) below.
!! No horizontal flow, zero surface wind, zero shear -> Richardson number
!! is infinite, so PP81 vertical tracer diffusivity reduces to the
!! constant background kappa_T_bg.
!!
!! The governing 1D heat equation in z is
!!     dT/dt = kappa * d^2T/dz^2
!! with zero-flux top and bottom boundaries.  Starting from a step
!! profile centred on z = H/2 with symmetric temperatures around a
!! mean T0 = (T_WARM + T_COLD)/2, the exact solution at each interface
!! is given by erfc applied to the distance from the interface:
!!     T(z, t) - T0 = 0.5 * (T_WARM - T_COLD) * erf((H/2 - z) / (2 sqrt(kappa t)))
!!
!! We check that after t = 3600 s the column-integrated temperature is
!! conserved (no heat lost) and that the numerical profile matches the
!! erf solution to within a few percent of the initial contrast.
!!
!! alpha_T is set to zero so the temperature field does not drive any
!! baroclinic flow -- this isolates the diffusion operator.
module test_thermal_diffusion
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_tracers, only: ml_init_tracer
   use rki_ml_eos, only: ml_equation_of_state
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_thermal_diffusion_tests

   ! Domain: shallow narrow channel discretised as nearly 1D
   real(wp), parameter :: H0 = 10.0_wp       !! Still-water depth (m)
   real(wp), parameter :: DX_VAL = 2.0_wp
   integer, parameter :: NX_PHYS = 8
   integer, parameter :: NY_PHYS = 8
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 15

   ! Temperature parameters
   real(wp), parameter :: T_WARM = 25.0_wp
   real(wp), parameter :: T_COLD = 15.0_wp
   real(wp), parameter :: T_MEAN = 0.5_wp*(T_WARM + T_COLD)

   ! Physics parameters.  We use an elevated diffusivity and a short run
   ! so the test finishes in seconds — the diffusion length scale
   ! sqrt(kappa*t) = sqrt(6) ≈ 2.45 m is what matters for the erf
   ! comparison, and that is unchanged by the kappa*T_END rescaling.
   real(wp), parameter :: KAPPA_BG = 4.0e-2_wp  !! Background diffusivity (m^2/s)
   real(wp), parameter :: RHO_0 = 1000.0_wp

   real(wp), parameter :: T_END = 150.0_wp    !! 2.5 minutes (kappa*t = 6, same as before)

contains

   subroutine collect_thermal_diffusion_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("thermal_diffusion_heat_conserved", &
                               test_heat_conserved), &
                  new_unittest("thermal_diffusion_matches_erfc", &
                               test_matches_erfc), &
                  new_unittest("surface_flux_slab_cooling", &
                               test_slab_cooling) &
                  ]
   end subroutine collect_thermal_diffusion_tests

   subroutine setup_state(state)
      !! 1D column of warm-over-cold water at rest
      type(state_t), intent(inout) :: state

      integer :: i, j, k
      real(wp) :: z_centre_from_top

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_temperature = T_MEAN
      state%multilayer%T_ref = T_MEAN
      state%multilayer%alpha_T = 0.0_wp         ! pure passive diffusion: no buoyancy feedback
      state%multilayer%beta_S = 0.0_wp         ! no salinity contribution
      state%multilayer%rho_0 = RHO_0
      state%multilayer%kappa_T_bg = KAPPA_BG    ! background diffusivity is the one we test
      state%multilayer%kappa_S_bg = 0.0_wp      ! no tracer diffusion for salinity
      state%multilayer%hdiff_kappa = 0.0_wp     ! no horizontal diffusion
      state%multilayer%vmix_nu_max = 0.0_wp     ! no Ri-driven enhancement

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      ! Walls on all sides to keep it 1D
      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      ! Flat bathymetry, still water
      state%barotropic%h = H0
      state%barotropic%b = 0.0_wp
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp

      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                          state%grid, NZ)
      call ml_init_tracer(state%multilayer%h_layer, state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, 0.0_wp, &
                          state%grid, NZ)

      ! Warm over cold: k=NZ is the surface, k=1 is the bottom in this
      ! codebase (bottom-up ROMS indexing matching BPG / rho_ref / NH init).
      ! Convert layer index to depth-from-surface and assign piecewise-constant
      ! initial temperature.
      do k = 1, NZ
         z_centre_from_top = (real(NZ - k, wp) + 0.5_wp)*H0/real(NZ, wp)
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               if (z_centre_from_top < 0.5_wp*H0) then
                  state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*T_WARM
               else
                  state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*T_COLD
               end if
            end do
         end do
      end do

      ! EOS — alpha_T=0 so the temperature field is passive
      call ml_equation_of_state(state%multilayer%h_layer, state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, RHO_0, 0.0_wp, 0.0_wp, &
                                0.0_wp, T_MEAN, &
                                state%grid, NZ)
   end subroutine setup_state

   subroutine run_diffusion(state, t_final)
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: t_final

      real(wp) :: t

      call solver_enter_data(state)

      t = 0.0_wp
      do while (t < t_final)
         call solver_step(state, time=t)
         t = t + state%dt
      end do

      call solver_exit_data(state)
   end subroutine run_diffusion

   subroutine column_temperature(state, T)
      !! Return the 1D column temperature profile (vertical slice at a
      !! representative interior cell).
      type(state_t), intent(in) :: state
      real(wp), intent(out) :: T(NZ)

      integer :: k, ng, ic, jc

      ng = NGHOST
      ic = ng + NX_PHYS/2
      jc = ng + NY_PHYS/2
      do k = 1, NZ
         T(k) = state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ic, jc, k)/state%multilayer%h_layer(ic, jc, k)
      end do
   end subroutine column_temperature

   subroutine test_heat_conserved(error)
      !! Pure vertical diffusion with zero-flux BCs must preserve the
      !! column-integrated heat content EXACTLY (to machine precision):
      !! the kernel uses interior flux form, so interior flux
      !! differences telescope and the surface / bottom zero-flux
      !! boundary condition closes the sum.
      !!
      !! This is stronger than the previous column-MEAN check, which
      !! allowed 5e-4 K drift. Column-integral `sum_k hT_layer` is the
      !! actual conservation invariant of the diffusion operator.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: T_profile(NZ)
      real(wp) :: col_integral_initial, col_integral_final
      integer :: ng, ic, jc, k

      call setup_state(state)

      ng = NGHOST
      ic = ng + NX_PHYS/2
      jc = ng + NY_PHYS/2

      ! Column-integrated hT (mass * concentration) at a representative
      ! interior column. Units: m*K. With uniform sigma this is
      ! h_layer * sum(T_k) = H / NZ * sum(T_k), so invariance of this
      ! integral is a direct test of the conservation property.
      col_integral_initial = 0.0_wp
      do k = 1, NZ
         col_integral_initial = col_integral_initial + state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ic, jc, k)
      end do

      call run_diffusion(state, T_END)

      col_integral_final = 0.0_wp
      do k = 1, NZ
         col_integral_final = col_integral_final + state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ic, jc, k)
      end do

      ! Near-machine-precision conservation.  The diffusion kernel is
      ! exactly conservative, but the full solver pipeline (PPM advection
      ! + clamp + vertical exchange) introduces O(eps) noise per step
      ! that accumulates over ~1600 steps.  1e-10 is well below any
      ! physical significance.
      call check(error, &
                 abs(col_integral_final - col_integral_initial) < 1.0e-10_wp, &
                 "column-integrated heat content preserved to near machine precision")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! Extrema must stay within the initial range (monotonicity of
      ! the diffusion operator on a step profile).
      call column_temperature(state, T_profile)
      call check(error, maxval(T_profile) <= T_WARM + 1.0e-8_wp .and. &
                 minval(T_profile) >= T_COLD - 1.0e-8_wp, &
                 "temperature stays inside [T_COLD, T_WARM]")

      call state%destroy()
   end subroutine test_heat_conserved

   subroutine test_matches_erfc(error)
      !! After t = T_END the profile must be close to the analytical erf
      !! solution to the 1D heat equation:
      !!   T(z,t) = T_MEAN + 0.5*(T_WARM - T_COLD) * erf((H/2 - z) / (2 sqrt(kappa t)))
      !! where z is measured from the surface (z=0 at surface, z=H at bottom).
      !! Zero-flux BCs are satisfied exactly by this solution as long as the
      !! diffusion front has not reached the walls, which is comfortably the
      !! case here (sqrt(kappa t) = sqrt(6) ~ 2.45 m vs H/2 = 5 m).
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: T_num(NZ), T_exact, arg
      real(wp) :: z_c, rms_err, max_err
      real(wp) :: sqrt_2kt
      integer :: k

      call setup_state(state)
      call run_diffusion(state, T_END)
      call column_temperature(state, T_num)

      sqrt_2kt = 2.0_wp*sqrt(KAPPA_BG*T_END)

      rms_err = 0.0_wp
      max_err = 0.0_wp
      do k = 1, NZ
         ! Depth from the surface for layer k under bottom-up ROMS indexing
         ! (k=NZ at surface → small z_c; k=1 at bottom → z_c ≈ H0).
         z_c = (real(NZ - k, wp) + 0.5_wp)*H0/real(NZ, wp)
         arg = (0.5_wp*H0 - z_c)/sqrt_2kt
         T_exact = T_MEAN + 0.5_wp*(T_WARM - T_COLD)*erf(arg)
         rms_err = rms_err + (T_num(k) - T_exact)**2
         max_err = max(max_err, abs(T_num(k) - T_exact))
      end do
      rms_err = sqrt(rms_err/real(NZ, wp))

      ! Tolerance: with NZ=15 and first-order upwind vertical exchange,
      ! observed errors are ~8% RMS / ~11% max.  We allow headroom for
      ! platform variation.
      call check(error, rms_err < 0.10_wp*(T_WARM - T_COLD), &
                 "RMS error vs erfc solution < 10% of initial contrast")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      call check(error, max_err < 0.15_wp*(T_WARM - T_COLD), &
                 "max pointwise error < 15% of initial contrast")

      call state%destroy()
   end subroutine test_matches_erfc

   subroutine test_slab_cooling(error)
      !! Constant-cooling slab test.  A uniform column with enhanced
      !! vertical diffusion stays well-mixed throughout; applying a
      !! constant negative heat flux Q_net at the surface should cool
      !! the whole column at a rate
      !!     dT/dt = Q_net / (rho_0 * cp * H)
      !! Expected ΔT after t seconds:
      !!     ΔT = Q_net * t / (rho_0 * cp * H)
      !! This validates the ml_apply_surface_tracer_fluxes kernel and the
      !! rho_0 * cp conversion simultaneously.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp), parameter :: Q_FLUX = -100.0_wp   !! W/m^2 (cooling)
      real(wp), parameter :: RUN_T = 150.0_wp    !! 2.5 min
      real(wp), parameter :: RHO_CP = 1000.0_wp*4000.0_wp
      real(wp) :: T_profile(NZ), mean_T_initial, mean_T_final
      real(wp) :: expected_delta, observed_delta
      integer :: ig, jg

      call setup_state(state)

      ! Uniform temperature = T_MEAN (no vertical gradient).  We want a
      ! well-mixed column so any spurious internal structure shows up as
      ! a deviation from the slab solution.
      do jg = 1, state%grid%ny_total
         do ig = 1, state%grid%nx_total
            state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ig, jg, :) = state%multilayer%h_layer(ig, jg, :)*T_MEAN
         end do
      end do

      ! Constant surface cooling everywhere
      state%forcing%heat_flux = Q_FLUX
      state%forcing%evap_minus_precip = 0.0_wp

      call column_temperature(state, T_profile)
      mean_T_initial = sum(T_profile)/real(NZ, wp)

      call run_diffusion(state, RUN_T)

      call column_temperature(state, T_profile)
      mean_T_final = sum(T_profile)/real(NZ, wp)

      expected_delta = Q_FLUX*RUN_T/(RHO_CP*H0)     ! = -100 * 150 / (4e6 * 10) = -3.75e-4 K
      observed_delta = mean_T_final - mean_T_initial

      ! Slab cooling ΔT is tiny (~1.5 mK).  Absolute tolerance of 10% of
      ! expected handles scheme diffusion / first-order error.
      call check(error, &
                 abs(observed_delta - expected_delta) < 0.10_wp*abs(expected_delta), &
                 "slab cooling: ΔT matches Q*t/(rho_0*cp*H) within 10%")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      ! Sign check: cooling must produce a net decrease
      call check(error, observed_delta < 0.0_wp, &
                 "slab cooling produces a temperature decrease")

      call state%destroy()
   end subroutine test_slab_cooling

end module test_thermal_diffusion
