!! Smoke tests for the full state lifecycle.
!!
!! The bug class this catches: integration-level GPU / OpenACC issues
!! that per-kernel unit tests miss because individual kernels never run
!! the full `solver_enter_data → solver_step → diagnostic →
!! solver_exit_data` sequence.  Concrete example: the
!! `state%forcing` partial-presence error introduced by the
!! `surface_forcing_t` extraction would have failed *here* on GPU,
!! despite all 73 prior unit tests passing.
!!
!! The test file is grid-bound but otherwise minimal so a regression
!! gets surfaced in seconds, on every CI run, on host *and* GPU.
module test_driver_smoke
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_init, only: init_dam_break, compute_total_mass
   implicit none
   private

   public :: collect_driver_smoke_tests

contains

   subroutine collect_driver_smoke_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("state_lifecycle_with_diagnostics", &
                               test_state_lifecycle_with_diagnostics), &
                  new_unittest("multiple_steps_no_drift", test_multiple_steps_no_drift) &
                  ]
   end subroutine collect_driver_smoke_tests

   subroutine test_state_lifecycle_with_diagnostics(error)
      !! Reproduces the integration sequence used by `driver_run`:
      !! init → enter_data → step → compute_total_mass → exit_data.
      !!
      !! The `compute_total_mass` call is the trip-wire — it accesses
      !! `s%barotropic%h(i,j)` from inside an `!$acc parallel loop`, which fails
      !! with a partial-presence error if any sub-struct of `state`
      !! has been attached without first attaching `state` itself.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: mass_before, mass_after

      call state%init(8, 8, 2, 1.0_wp, 1.0_wp, 0.45_wp)
      call init_dam_break(state)

      call solver_enter_data(state)
      mass_before = compute_total_mass(state)

      call solver_step(state)
      mass_after = compute_total_mass(state)

      call solver_exit_data(state)

      call check(error, mass_before > 0.0_wp, &
                 "compute_total_mass before step returns positive mass")
      if (allocated(error)) then; call state%destroy(); return; end if

      call check(error, mass_after > 0.0_wp, &
                 "compute_total_mass after step returns positive mass")
      if (allocated(error)) then; call state%destroy(); return; end if

      ! Closed dam-break with wall BCs conserves mass to machine precision
      call check(error, abs(mass_after - mass_before)/mass_before < 1.0e-10_wp, &
                 "mass conserved across one step")

      call state%destroy()
   end subroutine test_state_lifecycle_with_diagnostics

   subroutine test_multiple_steps_no_drift(error)
      !! Same lifecycle, multiple steps, including a mid-loop diagnostic
      !! call (mirrors the status-print path in `driver_run`).
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: mass_initial, mass_mid, mass_final
      integer :: n

      call state%init(16, 16, 2, 1.0_wp, 1.0_wp, 0.45_wp)
      call init_dam_break(state)

      call solver_enter_data(state)
      mass_initial = compute_total_mass(state)

      do n = 1, 5
         call solver_step(state)
      end do
      mass_mid = compute_total_mass(state)   ! mid-loop diagnostic

      do n = 1, 5
         call solver_step(state)
      end do
      mass_final = compute_total_mass(state)

      call solver_exit_data(state)

      call check(error, abs(mass_mid - mass_initial)/mass_initial < 1.0e-10_wp, &
                 "mass conserved over first 5 steps")
      if (allocated(error)) then; call state%destroy(); return; end if

      call check(error, abs(mass_final - mass_initial)/mass_initial < 1.0e-10_wp, &
                 "mass conserved over 10 steps")

      call state%destroy()
   end subroutine test_multiple_steps_no_drift

end module test_driver_smoke
