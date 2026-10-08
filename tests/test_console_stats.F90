!! Unit tests for the shared console-stats formatter + the coastal
!! kinetic-energy / max-CFL reductions that feed it.
module test_console_stats
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_console_stats, only: console_stats_t, console_stats_report
   use rki_coastal_console_stats, only: coastal_total_ke, coastal_max_cfl
   use rki_solver, only: accumulate_boundary_mass_out, &
                         accumulate_mass_src, &
                         solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_console_stats_tests

contains

   subroutine collect_console_stats_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("coastal_ke_analytical", test_ke), &
                  new_unittest("coastal_maxcfl_analytical", test_cfl), &
                  new_unittest("formatter_latches_reference", test_formatter), &
                  new_unittest("mass_out_accumulation", test_mass_out), &
                  new_unittest("mass_src_precip_sign_and_weight", test_mass_src_precip), &
                  new_unittest("mass_src_evap_mass_sink", test_mass_src_evap), &
                  new_unittest("mass_src_zero_forcing_byteident", test_mass_src_zero) &
                  ]
   end subroutine collect_console_stats_tests

   subroutine test_mass_out(error)
      !! The boundary-outflux accumulator integrates the flux divergence:
      !! for a uniform interior flux_h = C, one RK2 stage (weight 0.5) adds
      !! −0.5·dt·C·(nx·ny·area) to mass_out, and two stages double it.  This
      !! pins the sign + stage weight that make the mass budget close.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      integer, parameter :: nx = 8, ny = 6, ng = 2
      real(wp), parameter :: dx = 10.0_wp, dy = 20.0_wp, C = -0.25_wp, dt = 3.0_wp
      real(wp) :: expect

      call s%init(nx, ny, ng, dx, dy, 0.45_wp)
      s%dt = dt
      s%barotropic%flux_h = C   ! uniform tendency incl. ghosts (interior summed)

      call accumulate_boundary_mass_out(s, 0.5_wp)
      call accumulate_boundary_mass_out(s, 0.5_wp)

      ! Two 0.5 stages ⇒ full −dt·C·Σarea over the interior.
      expect = -dt*C*real(nx*ny, wp)*dx*dy
      call check(error, abs(s%barotropic%mass_out - expect) <= 1.0e-9_wp*abs(expect), &
                 "mass_out must integrate −Σ flux_h·area·dt over the RK stages")
      call check(error, s%barotropic%mass_out_tracked, "accumulation sets the tracked flag")
      call s%destroy()
   end subroutine test_mass_out

   subroutine test_ke(error)
      !! Uniform barotropic state h=H, u=U, v=0 over an nx*ny interior:
      !! KE = Σ 0.5·(hu²)/h·area = nx·ny · 0.5·H·U² · dx·dy.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      integer, parameter :: nx = 8, ny = 6, ng = 2
      real(wp), parameter :: dx = 10.0_wp, dy = 20.0_wp, H = 3.0_wp, U = 0.5_wp
      real(wp) :: ke, expect

      call s%init(nx, ny, ng, dx, dy, 0.45_wp)
      s%barotropic%h = H
      s%barotropic%hu = H*U
      s%barotropic%hv = 0.0_wp

      ! The reduction assumes device residency (present); map the state the
      ! same way production does, so the GPU build reads the on-device arrays.
      call solver_enter_data(s)
      ke = coastal_total_ke(s)
      call solver_exit_data(s)
      expect = real(nx*ny, wp)*0.5_wp*H*U*U*dx*dy
      call check(error, abs(ke - expect) <= 1.0e-9_wp*expect, &
                 "coastal KE must match the analytical Σ 0.5·H·U²·area")
      call s%destroy()
   end subroutine test_ke

   subroutine test_cfl(error)
      !! max advective CFL for u=U, v=0 is |U|·dt/dx (uniform state).
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      integer, parameter :: nx = 8, ny = 6, ng = 2
      real(wp), parameter :: dx = 10.0_wp, dy = 20.0_wp, H = 3.0_wp, U = 0.5_wp
      real(wp), parameter :: dt = 2.0_wp
      real(wp) :: mcfl, expect

      call s%init(nx, ny, ng, dx, dy, 0.45_wp)
      s%barotropic%h = H
      s%barotropic%hu = H*U
      s%barotropic%hv = 0.0_wp

      ! Map the state to device (present-clause contract), as production does.
      call solver_enter_data(s)
      mcfl = coastal_max_cfl(s, dt)
      call solver_exit_data(s)
      expect = U*dt/dx
      call check(error, abs(mcfl - expect) <= 1.0e-12_wp, &
                 "coastal max CFL must equal |U|·dt/dx for a uniform state")
      call s%destroy()
   end subroutine test_cfl

   subroutine test_formatter(error)
      !! The shared formatter must latch the t=0 reference on the first call
      !! (so later drift is relative to it) and run without crashing. Panic
      !! guards left on but fed finite, sub-threshold values.
      type(error_type), allocatable, intent(out) :: error
      type(console_stats_t) :: stats

      call check(error,.not. stats%is_initialised, "starts un-initialised")
      if (allocated(error)) return

      call console_stats_report(stats, t=0.0_wp, step=0, &
                                total_mass=100.0_wp, total_ke=5.0_wp, &
                                mean_S=35.0_wp, mean_T=15.0_wp, max_cfl=0.1_wp, &
                                total_salt=3500.0_wp, total_heat=1500.0_wp, &
                                has_salt=.true., has_temp=.true., &
                                mean_age=0.0_wp, has_age=.false.)

      call check(error, stats%is_initialised, "latches after first report")
      if (allocated(error)) return
      call check(error, abs(stats%mass0 - 100.0_wp) <= 1.0e-12_wp, &
                 "mass0 reference latched from the first call")
      if (allocated(error)) return

      ! Second call must not re-latch (reference stays at the first value).
      call console_stats_report(stats, t=3600.0_wp, step=10, &
                                total_mass=101.0_wp, total_ke=6.0_wp, &
                                mean_S=35.0_wp, mean_T=15.0_wp, max_cfl=0.2_wp, &
                                total_salt=3600.0_wp, total_heat=1500.0_wp, &
                                has_salt=.true., has_temp=.true., &
                                mean_age=0.0_wp, has_age=.false.)
      call check(error, abs(stats%mass0 - 100.0_wp) <= 1.0e-12_wp, &
                 "reference must not re-latch on later calls")
   end subroutine test_formatter

   subroutine test_mass_src_precip(error)
      !! Uniform precipitation field P > 0: `accumulate_mass_src` with
      !! evap_mass_sink=.false. must add `+stage_weight·dt·P·Σarea` per call.
      !! Two 0.5-weight stages sum to `+dt·P·Σarea` — same magnitude as
      !! `mass_out` but OPPOSITE sign direction (it offsets a precip-driven
      !! apparent "leak").
      !!
      !! Mutation check: if the sign were flipped (−P instead of +P), the
      !! accumulator would produce the wrong sign and the budget would fail
      !! this check.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      integer, parameter :: nx = 8, ny = 6, ng = 2
      real(wp), parameter :: dx = 10.0_wp, dy = 20.0_wp
      real(wp), parameter :: P = 1.0e-3_wp   ! 1 mm/s uniform precip
      real(wp), parameter :: dt = 3.0_wp
      real(wp) :: expect

      call s%init(nx, ny, ng, dx, dy, 0.45_wp)
      s%dt = dt
      s%forcing%active = .true.
      s%forcing%precipitation = P      ! uniform field (ghost rows included)
      s%barotropic%h = 1.0_wp

      ! Two SSP-RK2 stages with weight 0.5 each.
      call accumulate_mass_src(s, 0.5_wp, evap_mass_sink=.false.)
      call accumulate_mass_src(s, 0.5_wp, evap_mass_sink=.false.)

      ! Two 0.5 stages → full +dt·P·Σ_interior area.
      expect = dt*P*real(nx*ny, wp)*dx*dy
      call check(error, abs(s%barotropic%mass_src - expect) <= 1.0e-9_wp*abs(expect), &
                 "mass_src must be +dt·P·Σarea for uniform precip (two 0.5-weight stages)")
      if (allocated(error)) return
      call check(error, s%barotropic%mass_src_tracked, &
                 "mass_src_tracked must be set after accumulate_mass_src")
      call s%destroy()
   end subroutine test_mass_src_precip

   subroutine test_mass_src_evap(error)
      !! Uniform evaporation-only field (evap_minus_precip = E > 0, no rain):
      !! `accumulate_mass_src` with `evap_mass_sink=.true.` must accumulate
      !! `−dt·E·Σarea` (negative; evaporation REMOVES mass from the domain).
      !!
      !! This is the budget's complement to `mass_out`: on a forced run with
      !! only evaporation and no open boundaries, `(total − ref)` drifts
      !! negative, and `mass_src` corrects it so the residual stays ~0.
      !!
      !! Mutation check: if the sign were flipped or the evap_mass_sink branch
      !! were ignored, `mass_src` would be positive and the budget residual
      !! would be ~2×drift instead of 0.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      integer, parameter :: nx = 8, ny = 6, ng = 2
      real(wp), parameter :: dx = 10.0_wp, dy = 20.0_wp
      real(wp), parameter :: E = 5.0e-4_wp   ! 0.5 mm/s evap rate (E−P = E)
      real(wp), parameter :: dt = 3.0_wp
      real(wp) :: expect

      call s%init(nx, ny, ng, dx, dy, 0.45_wp)
      s%dt = dt
      s%forcing%active = .true.
      s%forcing%evap_minus_precip = E    ! evap only, no precip → E−P = E > 0
      s%forcing%precipitation = 0.0_wp
      s%barotropic%h = 1.0_wp

      call accumulate_mass_src(s, 0.5_wp, evap_mass_sink=.true.)
      call accumulate_mass_src(s, 0.5_wp, evap_mass_sink=.true.)

      ! Two 0.5 stages → −dt·E·Σ_interior area (evap removes mass).
      expect = -dt*E*real(nx*ny, wp)*dx*dy
      call check(error, abs(s%barotropic%mass_src - expect) <= 1.0e-9_wp*abs(expect), &
                 "mass_src must be −dt·E·Σarea for evap-only forcing (evap_mass_sink=.true.)")
      call s%destroy()
   end subroutine test_mass_src_evap

   subroutine test_mass_src_zero(error)
      !! With forcing off (state%forcing%active=.false.) the accumulator must
      !! be a no-op: mass_src stays 0.  This guards byte-identity on unforced
      !! runs (no precip / no evap / no sponge configured → no change).
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      integer, parameter :: nx = 8, ny = 6, ng = 2
      real(wp), parameter :: dx = 10.0_wp, dy = 20.0_wp, dt = 3.0_wp

      call s%init(nx, ny, ng, dx, dy, 0.45_wp)
      s%dt = dt
      s%forcing%active = .false.    ! Forcing disabled → no precip/evap.
      s%barotropic%h = 1.0_wp

      call accumulate_mass_src(s, 0.5_wp, evap_mass_sink=.false.)
      call accumulate_mass_src(s, 0.5_wp, evap_mass_sink=.false.)

      call check(error, s%barotropic%mass_src == 0.0_wp, &
                 "mass_src must remain 0 when forcing is off (byte-identity gate)")
      if (allocated(error)) return
      call check(error, s%barotropic%mass_src_tracked, &
                 "mass_src_tracked must still be set even when forcing is off")
      call s%destroy()
   end subroutine test_mass_src_zero

end module test_console_stats
