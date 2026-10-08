!! Cross-validating integration test: surface flux × stratification × BPG.
!!
!! Catches future drift between the surface-flux kernel orientation and the
!! BPG / rho_ref / NH-init bottom-up indexing.  If a future refactor were to
!! flip the surface-flux kernel back to writing at k=1 (treating it as the
!! "top") while the BPG continues to treat k=1 as the bottom, this test
!! fails immediately:
!!   - the warming would land at k=1 (the actual bottom),
!!   - the column would invert (warm at bottom → unstable),
!!   - the rho(top) < rho(bottom) check at the end would fail.
!!
!! The setup is deliberately minimal: closed basin, no horizontal density
!! gradient, no horizontal forcing.  The only thing that can move tracer is
!! (1) the surface flux kernel, (2) very small diffusive coupling between
!! layers.  That makes the orientation cross-check unambiguous.
module test_surface_flux_strat
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, CP_WATER
   use rki_state, only: state_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_eos, only: ml_equation_of_state
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_surface_flux_strat_tests

   ! Domain — small closed basin with non-trivial vertical resolution
   integer, parameter :: NX_PHYS = 6
   integer, parameter :: NY_PHYS = 6
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   real(wp), parameter :: H0 = 4.0_wp        !! 4 m water column
   real(wp), parameter :: DX_VAL = 1.0_wp

   ! Stable stratification: warm/light at the top (k=NZ),
   ! cool/dense at the bottom (k=1)
   real(wp), parameter :: T_BOTTOM = 5.0_wp
   real(wp), parameter :: T_SURFACE = 15.0_wp
   real(wp), parameter :: S_INIT = 35.0_wp

   ! EOS: temperature-active, salinity passive (uniform)
   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: BETA_S = 0.78_wp
   real(wp), parameter :: ALPHA_T = 0.17_wp
   real(wp), parameter :: T_REF = 10.0_wp
   real(wp), parameter :: S_REF = 0.0_wp

contains

   subroutine collect_surface_flux_strat_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("warming_lands_at_top_layer", &
                               test_warming_lands_at_top), &
                  new_unittest("stable_stratification_preserved", &
                               test_stratification_preserved) &
                  ]
   end subroutine collect_surface_flux_strat_tests

   subroutine setup_stratified_state(state, q_surface)
      !! Build a closed-basin state with stable T-stratification and a
      !! prescribed surface heat flux.
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: q_surface

      integer :: i, j, k
      real(wp) :: T_lin

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_salinity = S_INIT
      state%multilayer%initial_temperature = 0.5_wp*(T_BOTTOM + T_SURFACE)
      state%multilayer%alpha_T = ALPHA_T
      state%multilayer%beta_S = BETA_S
      state%multilayer%rho_0 = RHO_0
      state%multilayer%T_ref = T_REF
      state%multilayer%S_ref = S_REF

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      ! Closed basin — no horizontal forcing of any kind
      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      ! Flat bathymetry, still water
      state%barotropic%h = H0
      state%barotropic%b = 0.0_wp
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp

      ! Build uniform sigma layers (h_layer = H0 / NZ everywhere)
      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      ! Stable stratification: linear T profile from T_BOTTOM at k=1 to
      ! T_SURFACE at k=NZ (bottom-up ROMS indexing).  Uniform salinity, so
      ! the vertical gradient comes purely from temperature via EOS.
      do k = 1, NZ
         T_lin = T_BOTTOM + (T_SURFACE - T_BOTTOM) &
                 *real(k - 1, wp)/real(NZ - 1, wp)
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
                                RHO_0, BETA_S, S_REF, ALPHA_T, T_REF, &
                                state%grid, NZ)

      ! Apply uniform positive surface heat flux (warming the ocean).
      state%forcing%heat_flux = q_surface
      state%forcing%evap_minus_precip = 0.0_wp
   end subroutine setup_stratified_state

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

   subroutine test_warming_lands_at_top(error)
      !! Positive surface heat flux must increase T at the SURFACE (k=NZ)
      !! more than at the bottom (k=1).  If the surface-flux kernel ever
      !! drifts to writing at k=1 again, dT_bottom would dominate dT_top
      !! and the warm-on-top gradient would collapse — both assertions
      !! below catch that immediately.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp), parameter :: Q_WARM = 1000.0_wp  !! W/m^2 (strong, for SNR)
      real(wp), parameter :: T_END = 30.0_wp     !! seconds
      integer :: ic, jc
      real(wp) :: T_top_initial, T_bottom_initial
      real(wp) :: T_top_final, T_bottom_final
      real(wp) :: dT_top, dT_bottom
      real(wp) :: gradient_initial, gradient_final

      call setup_stratified_state(state, q_surface=Q_WARM)

      ic = NGHOST + NX_PHYS/2
      jc = NGHOST + NY_PHYS/2

      T_top_initial = state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ic, jc, NZ) &
                      /state%multilayer%h_layer(ic, jc, NZ)
      T_bottom_initial = state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ic, jc, 1) &
                         /state%multilayer%h_layer(ic, jc, 1)
      gradient_initial = T_top_initial - T_bottom_initial

      call run_to(state, T_END)

      T_top_final = state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ic, jc, NZ) &
                    /state%multilayer%h_layer(ic, jc, NZ)
      T_bottom_final = state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(ic, jc, 1) &
                       /state%multilayer%h_layer(ic, jc, 1)
      gradient_final = T_top_final - T_bottom_final

      dT_top = T_top_final - T_top_initial
      dT_bottom = T_bottom_final - T_bottom_initial

      ! Positive heat flux must warm the top
      call check(error, dT_top > 0.0_wp, &
                 "surface heat flux must increase T at the top layer (k=NZ)")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! Top must warm strictly more than the bottom.  Under the bug
      ! (heat flux misdirected to k=1), this inequality flips.
      call check(error, dT_top > dT_bottom, &
                 "top warming must dominate bottom warming")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! The warm-on-top gradient must not collapse.  Mixing can erode it
      ! by a fraction; a misdirected surface flux would crush it (and
      ! flip its sign) within a few steps.  Tolerance: gradient may shrink
      ! by no more than 10% of its initial value over the run.
      call check(error, gradient_final > 0.9_wp*gradient_initial, &
                 "warm-on-top stratification gradient must be preserved")

      call state%destroy()
   end subroutine test_warming_lands_at_top

   subroutine test_stratification_preserved(error)
      !! With a stably-stratified IC and surface warming, the column must
      !! REMAIN stably stratified: rho(top) < rho(bottom) at every probe.
      !! If surface flux misdirects heat to k=1, the bottom warms, becomes
      !! lighter than the surface, the column inverts, and this assertion
      !! fails.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp), parameter :: Q_WARM = 400.0_wp   !! W/m^2 (strong)
      real(wp), parameter :: T_END = 10.0_wp
      integer :: ic, jc
      real(wp) :: rho_top_initial, rho_bottom_initial
      real(wp) :: rho_top_final, rho_bottom_final

      call setup_stratified_state(state, q_surface=Q_WARM)

      ic = NGHOST + NX_PHYS/2
      jc = NGHOST + NY_PHYS/2

      rho_top_initial = state%multilayer%rho_layer(ic, jc, NZ)
      rho_bottom_initial = state%multilayer%rho_layer(ic, jc, 1)

      ! Sanity: IC must already be stably stratified
      call check(error, rho_top_initial < rho_bottom_initial, &
                 "IC sanity: rho(top, k=NZ) < rho(bottom, k=1)")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      call run_to(state, T_END)

      ! Solver internally refreshes rho_layer each step via EOS in the
      ! multilayer pipeline; solver_exit_data syncs it back to the host.
      ! No post-loop EOS call needed (and it would fail under GPU because
      ! the device data has already been released).
      rho_top_final = state%multilayer%rho_layer(ic, jc, NZ)
      rho_bottom_final = state%multilayer%rho_layer(ic, jc, 1)

      call check(error, rho_top_final < rho_bottom_final, &
                 "after surface warming: rho(top) still < rho(bottom)")

      call state%destroy()
   end subroutine test_stratification_preserved

end module test_surface_flux_strat
