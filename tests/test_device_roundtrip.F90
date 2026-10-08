!! Device-data round trip: leaving and re-entering the device must not change a run.
module test_device_roundtrip
   !! `solver_exit_data` / `solver_enter_data` bracket every snapshot and every
   !! restart write in the driver.  Any state that is carried from one step to
   !! the next but is `delete`d on exit (instead of copied back) silently
   !! reverts to its last host value on re-entry, so a run's result would depend
   !! on how often it writes output.
   !!
   !! This was a real defect: the k-epsilon (tke, eps, nu_t, kappa_t), KPP and
   !! Chapman-BC `eta_old_*` arrays were deleted on exit.  A 600 s stratified
   !! k-epsilon run differed by up to 4e-5 relative between "output every
   !! 100 s" and "output only at the end" on a GPU build.  On a CPU build the
   !! exit/enter directives are no-ops, so this test passes trivially there; it
   !! is the GPU build that exercises it.
   !!
   !! The test runs a wind-forced, stratified k-epsilon column on a fixed
   !! timestep twice -- straight through, and with two exit/re-enter round trips
   !! in the middle -- and requires the final state to be bit-identical.  It
   !! also requires k-epsilon to have come back from the device with an evolved
   !! value; with the original `delete` that check is what fails, because the
   !! host never receives k-epsilon at all.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_eos, only: ml_equation_of_state
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_device_roundtrip_tests

   integer, parameter :: NX_PHYS = 6, NY_PHYS = 6, NGHOST = 2
   integer, parameter :: NZ = 12
   real(wp), parameter :: H0 = 30.0_wp
   real(wp), parameter :: DX_VAL = 1000.0_wp
   real(wp), parameter :: DT_FIXED = 5.0_wp
   integer, parameter :: N_STEPS = 60
   real(wp), parameter :: TAU = 0.5_wp
   real(wp), parameter :: T_BOTTOM = 16.0_wp, T_SURFACE = 20.0_wp
   real(wp), parameter :: S_INIT = 35.0_wp
   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: BETA_S = 0.78_wp, ALPHA_T = 0.2_wp
   real(wp), parameter :: T_REF = 18.0_wp, S_REF = 0.0_wp

contains

   subroutine collect_device_roundtrip_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("exit_enter_does_not_change_the_run", test_roundtrip) &
                  ]
   end subroutine collect_device_roundtrip_tests

   subroutine setup_state(state)
      !! Resting, linearly stratified column, k-epsilon on, constant wind.
      type(state_t), intent(inout) :: state
      integer :: i, j, k
      real(wp) :: t_lin

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

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")
      state%barotropic%coriolis_f = 0.0_wp
      state%barotropic%wind_surface_layer = .true.
      state%dt_fixed = DT_FIXED

      state%barotropic%h = H0
      state%barotropic%b = 0.0_wp
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp

      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, &
                                          state%barotropic%hv, state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, state%grid, NZ)

      do k = 1, NZ
         t_lin = T_BOTTOM + (T_SURFACE - T_BOTTOM)*real(k - 1, wp)/real(NZ - 1, wp)
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = &
                  state%multilayer%h_layer(i, j, k)*t_lin
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

      state%forcing%active = .true.
      state%forcing%wind_x = TAU
      state%forcing%wind_y = 0.0_wp
      state%forcing%heat_flux = 0.0_wp
      state%forcing%evap_minus_precip = 0.0_wp
   end subroutine setup_state

   subroutine advance(state, n, t)
      !! `n` fixed steps starting at time `t`; the device data stays mapped
      !! only for the duration of this call (enter ... exit), which is exactly
      !! what the driver does around a snapshot or restart write.
      type(state_t), intent(inout) :: state
      integer, intent(in) :: n
      real(wp), intent(inout) :: t
      integer :: k

      call solver_enter_data(state)
      do k = 1, n
         call solver_step(state, time=t)
         t = t + state%dt
      end do
      call solver_exit_data(state)
   end subroutine advance

   function maxdiff3(a, b) result(d)
      real(wp), intent(in) :: a(:, :, :), b(:, :, :)
      real(wp) :: d
      d = maxval(abs(a - b))
   end function maxdiff3

   subroutine test_roundtrip(error)
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: s_ref, s_cut
      real(wp) :: t_ref, t_cut, d_tke, d_eps, d_hu, d_hT, tke_peak
      integer :: it
      checks: block

         call setup_state(s_ref)
         call setup_state(s_cut)

         ! Straight through
         t_ref = 0.0_wp
         call advance(s_ref, N_STEPS, t_ref)

         ! Same run, leaving and re-entering the device twice in the middle
         t_cut = 0.0_wp
         call advance(s_cut, N_STEPS/4, t_cut)
         call advance(s_cut, N_STEPS/4, t_cut)
         call advance(s_cut, N_STEPS/2, t_cut)

         call check(error, abs(t_ref - t_cut) < 1.0e-12_wp, "runs reached different end times")
         if (allocated(error)) exit checks

         ! k-epsilon must have come back from the device holding its evolved
         ! value.  If `solver_exit_data` deletes it instead of copying it
         ! back, the host keeps the initial value -- both runs then agree with
         ! each other trivially, so this check is what catches that defect
         ! (and it also guards the comparison below against being vacuous).
         tke_peak = maxval(s_ref%multilayer%keps_tke)
         call check(error, tke_peak > 1.0e-6_wp, &
                    "k-epsilon state was not copied back from the device (or never evolved)")
         if (allocated(error)) exit checks

         d_tke = maxdiff3(s_ref%multilayer%keps_tke, s_cut%multilayer%keps_tke)
         d_eps = maxdiff3(s_ref%multilayer%keps_eps, s_cut%multilayer%keps_eps)
         d_hu = maxdiff3(s_ref%multilayer%hu_layer, s_cut%multilayer%hu_layer)
         it = s_ref%multilayer%idx_temperature
         d_hT = maxdiff3(s_ref%multilayer%tracers(it)%hTr, s_cut%multilayer%tracers(it)%hTr)

         call check(error, d_tke == 0.0_wp, &
                    "TKE depends on how often the device data is exited and re-entered")
         if (allocated(error)) exit checks
         call check(error, d_eps == 0.0_wp, &
                    "dissipation depends on how often the device data is exited and re-entered")
         if (allocated(error)) exit checks
         call check(error, d_hu == 0.0_wp, &
                    "layer momentum depends on how often the device data is exited and re-entered")
         if (allocated(error)) exit checks
         call check(error, d_hT == 0.0_wp, &
                    "temperature depends on how often the device data is exited and re-entered")

      end block checks
      continue
      call s_ref%destroy()
      call s_cut%destroy()
   end subroutine test_roundtrip

end module test_device_roundtrip
