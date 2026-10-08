!! Kato-Phillips wind-driven mixing: solver-level integration test for the
!! k-epsilon vertical-mixing closure.
!!
!! A resting, linearly stratified column is forced by a constant surface
!! wind stress (no surface buoyancy flux).  The k-eps surface wall function
!! injects turbulence that erodes the stratification and forms a surface
!! mixed layer.  This exercises the full k-eps vmix path through the solver
!! (friction velocities + wall functions + momentum & tracer mixing) and
!! checks that a stable mixed layer forms.
!!
!! The sqrt(t) entrainment LAW itself, D(t) = c u_* sqrt(t/N0) (Kato &
!! Phillips 1969), is validated in the 1-D Python reference
!! (test_kato_phillips: c constant to <1%).
!! A clean sqrt(t) in the *solver* additionally needs surface-concentrated
!! wind stress (shear production at the ML base); the coastal multilayer
!! currently applies wind to the depth-uniform barotropic mode, so that is
!! tracked as a follow-up.
module test_kato_phillips
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_eos, only: ml_equation_of_state
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_kato_phillips_tests

   integer, parameter :: NX_PHYS = 6, NY_PHYS = 6, NGHOST = 2
   integer, parameter :: NZ = 30
   real(wp), parameter :: H0 = 60.0_wp           !! 60 m water column
   real(wp), parameter :: DX_VAL = 1000.0_wp
   real(wp), parameter :: TAU = 0.5_wp           !! surface wind stress (Pa)
   ! Stable stratification: warm/light surface (k=NZ), cool/dense bed (k=1)
   real(wp), parameter :: T_BOTTOM = 18.0_wp, T_SURFACE = 20.0_wp
   real(wp), parameter :: S_INIT = 35.0_wp
   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: BETA_S = 0.78_wp, ALPHA_T = 0.2_wp
   real(wp), parameter :: T_REF = 19.0_wp, S_REF = 0.0_wp

contains

   subroutine collect_kato_phillips_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("kato_phillips_wind_mixing", test_wind_mixing), &
                  new_unittest("kato_phillips_sqrt_t_deepening", test_sqrt_t_deepening) &
                  ]
   end subroutine collect_kato_phillips_tests

   subroutine setup_state(state, surface_wind)
      !! Resting, linearly stratified column with k-eps on and a constant
      !! surface wind stress.  `surface_wind` selects surface-concentrated
      !! top-layer stress vs the depth-uniform barotropic forcing.
      type(state_t), intent(inout) :: state
      logical, intent(in) :: surface_wind
      integer :: i, j, k
      real(wp) :: T_lin

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%vmix_use_keps = .true.
      state%multilayer%initial_salinity = S_INIT
      state%multilayer%initial_temperature = 0.5_wp*(T_BOTTOM + T_SURFACE)
      state%multilayer%alpha_T = ALPHA_T
      state%multilayer%beta_S = BETA_S
      state%multilayer%rho_0 = RHO_0
      state%multilayer%T_ref = T_REF
      state%multilayer%S_ref = S_REF

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      ! Closed basin, non-rotating: keep the dynamics 1-D so the only
      ! vertical mixing is the k-eps closure.
      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")
      state%barotropic%coriolis_f = 0.0_wp
      state%barotropic%wind_surface_layer = surface_wind

      state%barotropic%h = H0
      state%barotropic%b = 0.0_wp
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp

      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      ! Stable linear T profile (k=1 bottom cold, k=NZ surface warm); uniform
      ! salinity, so the stratification comes purely from T via the EOS.
      do k = 1, NZ
         T_lin = T_BOTTOM + (T_SURFACE - T_BOTTOM)*real(k - 1, wp)/real(NZ - 1, wp)
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = &
                  state%multilayer%h_layer(i, j, k)*T_lin
               state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = &
                  state%multilayer%h_layer(i, j, k)*S_INIT
            end do
         end do
      end do

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, &
                                RHO_0, BETA_S, S_REF, ALPHA_T, T_REF, state%grid, NZ)

      ! Constant surface wind stress (Pa); no surface buoyancy flux.
      state%forcing%active = .true.
      state%forcing%wind_x = TAU
      state%forcing%wind_y = 0.0_wp
      state%forcing%heat_flux = 0.0_wp
      state%forcing%evap_minus_precip = 0.0_wp
   end subroutine setup_state

   subroutine run_to(state, t_final)
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
   end subroutine run_to

   function mixed_layer_depth(state, ic, jc) result(d_ml)
      !! Depth (m below surface) of the base of the wind-mixed layer: the
      !! shallowest layer (scanning down from the surface k=NZ) whose density
      !! exceeds the surface value by more than 20% of the initial
      !! top-to-bottom density contrast.
      type(state_t), intent(in) :: state
      integer, intent(in) :: ic, jc
      real(wp) :: d_ml
      real(wp) :: rho_surf, drho_thresh, dz
      integer :: k, k_base

      rho_surf = state%multilayer%rho_layer(ic, jc, NZ)
      drho_thresh = 0.2_wp*(state%multilayer%rho_layer(ic, jc, 1) - rho_surf)
      dz = H0/real(NZ, wp)
      k_base = 1
      do k = NZ, 1, -1
         if (state%multilayer%rho_layer(ic, jc, k) - rho_surf > drho_thresh) then
            k_base = k
            exit
         end if
      end do
      d_ml = real(NZ - k_base, wp)*dz
   end function mixed_layer_depth

   subroutine test_wind_mixing(error)
      !! Wind forcing must erode the near-surface stratification, form a
      !! surface mixed layer, and keep the column stable.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: T_END = 7200.0_wp     !! 2 h
      integer :: ic, jc
      real(wp) :: d_ml, grad_init, grad_final

      ic = NGHOST + NX_PHYS/2
      jc = NGHOST + NY_PHYS/2

      call setup_state(state, surface_wind=.false.)
      ! near-surface density contrast over the top three layers
      grad_init = state%multilayer%rho_layer(ic, jc, NZ - 3) &
                  - state%multilayer%rho_layer(ic, jc, NZ)

      call run_to(state, T_END)

      d_ml = mixed_layer_depth(state, ic, jc)
      grad_final = state%multilayer%rho_layer(ic, jc, NZ - 3) &
                   - state%multilayer%rho_layer(ic, jc, NZ)

      ! A surface mixed layer formed: several layers deep, not the whole column.
      call check(error, d_ml > 3.0_wp*(H0/real(NZ, wp)) .and. d_ml < 0.8_wp*H0, &
                 "wind forcing forms a surface mixed layer")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! The near-surface stratification is strongly eroded (top is mixed).
      call check(error, grad_final < 0.2_wp*grad_init, &
                 "near-surface stratification eroded by k-eps wind mixing")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! The column remains stably stratified overall (no inversion / blow-up).
      call check(error, state%multilayer%rho_layer(ic, jc, 1) > &
                 state%multilayer%rho_layer(ic, jc, NZ), &
                 "column remains stably stratified (rho_bottom > rho_surface)")

      call state%destroy()
   end subroutine test_wind_mixing

   subroutine test_sqrt_t_deepening(error)
      !! With surface-concentrated wind stress (wind_surface_layer on), the
      !! solver should reproduce the Kato-Phillips sqrt(t) deepening law:
      !! D(4t)/D(t) -> 2.  Two independent runs at a 4x time ratio (robust to
      !! the closure- and ML-diagnostic-dependent prefactor).
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s1, s2
      real(wp), parameter :: T1 = 3600.0_wp, T2 = 14400.0_wp  !! 1 h, 4 h
      integer :: ic, jc
      real(wp) :: d1, d2, ratio

      ic = NGHOST + NX_PHYS/2
      jc = NGHOST + NY_PHYS/2

      call setup_state(s1, surface_wind=.true.)
      call run_to(s1, T1)
      d1 = mixed_layer_depth(s1, ic, jc)

      call setup_state(s2, surface_wind=.true.)
      call run_to(s2, T2)
      d2 = mixed_layer_depth(s2, ic, jc)

      ratio = d2/max(d1, 1.0e-6_wp)

      call check(error, d1 > 2.0_wp*(H0/real(NZ, wp)) .and. d2 < 0.85_wp*H0, &
                 "mixed layer deepened measurably but not bottom-limited")
      if (allocated(error)) then
         call s1%destroy()
         call s2%destroy()
         return
      end if

      call check(error, ratio > 1.6_wp .and. ratio < 2.4_wp, &
                 "Kato-Phillips sqrt(t) scaling with surface wind: D(4t)/D(t) ~ 2")

      call s1%destroy()
      call s2%destroy()
   end subroutine test_sqrt_t_deepening

end module test_kato_phillips
