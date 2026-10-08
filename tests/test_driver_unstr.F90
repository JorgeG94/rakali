!! Smoke test for the namelist-driven unstructured coastal driver path.
!!
!! `grid_type="unstructured"` routes `driver_run` to `driver_run_coastal_unstr`,
!! which generates a native regular triangulation, sets up a multilayer state,
!! and advances a lake-at-rest to `t_end` -- entirely in Fortran, no Python.
!! This exercises that whole path as a unit (mesh-gen -> state_unstr_init ->
!! enter_data -> solver_unstr_init_multilayer -> solver_unstr_advance ->
!! exit_data -> destroy); a control-flow / init-ordering break surfaces here.
!!
!! Inc A: the driver writes no output, so this asserts the path runs to
!! completion (no crash / abort).  Richer assertions (mass conservation,
!! mixed-layer behaviour) land once the driver emits UGRID output (Inc B/C),
!! mirroring `test_driver_coastal`.
module test_driver_unstr
   use rki_constants, only: wp
   use rki_config, only: config_t
   use rki_driver, only: driver_run, configure_log_level
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_driver_unstr_tests

contains

   subroutine collect_driver_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("driver_run_coastal_unstr_e2e", test_unstr_e2e), &
                  new_unittest("driver_run_coastal_unstr_forced", test_unstr_forced), &
                  new_unittest("driver_run_coastal_unstr_stratified", test_unstr_stratified), &
                  new_unittest("driver_run_coastal_unstr_si_2d", test_unstr_si_2d), &
                  new_unittest("driver_run_coastal_unstr_si_ml", test_unstr_si_ml) &
                  ]
   end subroutine collect_driver_unstr_tests

   subroutine test_unstr_si_2d(error)
      !! Namelist-reachability gate for the unstructured SEMI-IMPLICIT mode:
      !! `time_scheme="semi_implicit"` through the real `driver_run` entry
      !! must route into `state_unstr_enable_si` (quality gate + si slot) and
      !! the SI step.  This DISCRIMINATES, not just smokes: without the
      !! driver wiring the run error-stops in `make_timestepper`
      !! (TIME_SEMI_IMPLICIT has no explicit table), so completing to t_end
      !! proves the SI branch engaged.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call configure_log_level("warning")
      cfg%sim_type = "coastal"
      cfg%grid_type = "unstructured"
      cfg%nx = 12
      cfg%ny = 12
      cfg%dx = 100.0_wp
      cfg%dy = 100.0_wp
      cfg%nz_layers = 1
      cfg%h0 = 100.0_wp
      cfg%time_scheme = "semi_implicit"
      cfg%dt_fixed = 20.0_wp
      cfg%t_end = 60.0_wp
      cfg%output_to_file = .false.

      call driver_run(cfg)

      call check(error, .true., &
                 "unstructured SI 2D namelist path runs to t_end without aborting")
   end subroutine test_unstr_si_2d

   subroutine test_unstr_si_ml(error)
      !! The 3D twin of the SI reachability gate: multilayer + semi_implicit
      !! through `driver_run`.  Without the driver's SEMI_IMPLICIT exemption
      !! on the explicit-integrator multilayer guard this error-stops before
      !! the first step.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call configure_log_level("warning")
      cfg%sim_type = "coastal"
      cfg%grid_type = "unstructured"
      cfg%nx = 12
      cfg%ny = 12
      cfg%dx = 100.0_wp
      cfg%dy = 100.0_wp
      cfg%nz_layers = 3
      cfg%h0 = 100.0_wp
      cfg%time_scheme = "semi_implicit"
      cfg%dt_fixed = 20.0_wp
      cfg%t_end = 60.0_wp
      cfg%output_to_file = .false.

      call driver_run(cfg)

      call check(error, .true., &
                 "unstructured SI multilayer namelist path runs to t_end without aborting")
   end subroutine test_unstr_si_ml

   subroutine test_unstr_e2e(error)
      !! A tiny regular-mesh, multilayer lake-at-rest run to t_end through the
      !! real `driver_run` entry point.  Reaching the assertion means the
      !! unstructured driver path completed end-to-end without aborting.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call configure_log_level("warning")
      cfg%sim_type = "coastal"
      cfg%grid_type = "unstructured"
      cfg%nx = 12
      cfg%ny = 12
      cfg%dx = 100.0_wp
      cfg%dy = 100.0_wp
      cfg%nz_layers = 3
      cfg%h0 = 100.0_wp
      cfg%t_end = 5.0_wp
      cfg%output_to_file = .false.

      call driver_run(cfg)

      call check(error, .true., &
                 "unstructured driver path runs to t_end without aborting")
   end subroutine test_unstr_e2e

   subroutine test_unstr_forced(error)
      !! The full forced + mixed path (Inc B): wind stress + Coriolis + Manning
      !! bottom drag + k-epsilon vertical mixing + Smagorinsky lateral mixing,
      !! all plumbed from the namelist onto the unstructured state.  Must run to
      !! t_end without aborting or destabilising.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call configure_log_level("warning")
      cfg%sim_type = "coastal"
      cfg%grid_type = "unstructured"
      cfg%nx = 12
      cfg%ny = 12
      cfg%dx = 100.0_wp
      cfg%dy = 100.0_wp
      cfg%nz_layers = 4
      cfg%h0 = 100.0_wp
      cfg%t_end = 5.0_wp
      cfg%wind_stress_x = 0.1_wp
      cfg%coriolis_f = 1.0e-4_wp
      cfg%manning_n = 0.025_wp
      cfg%vmix_use_keps = .true.
      cfg%smagorinsky_cs = 0.1_wp
      cfg%output_to_file = .false.

      call driver_run(cfg)

      call check(error, .true., &
                 "forced + mixed unstructured run completes without aborting")
   end subroutine test_unstr_forced

   subroutine test_unstr_stratified(error)
      !! Stratified + forced unstr run (Inc B): a linear T(z) profile via
      !! T_init_surface/bottom on top of wind + Coriolis + k-epsilon mixing.
      !! Must run to t_end without aborting or destabilising.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call configure_log_level("warning")
      cfg%sim_type = "coastal"
      cfg%grid_type = "unstructured"
      cfg%nx = 12
      cfg%ny = 12
      cfg%dx = 100.0_wp
      cfg%dy = 100.0_wp
      cfg%nz_layers = 5
      cfg%h0 = 100.0_wp
      cfg%t_end = 5.0_wp
      cfg%T_init_bottom = 12.0_wp
      cfg%T_init_surface = 22.0_wp
      cfg%wind_stress_x = 0.1_wp
      cfg%coriolis_f = 1.0e-4_wp
      cfg%vmix_use_keps = .true.
      cfg%output_to_file = .false.

      call driver_run(cfg)

      call check(error, .true., &
                 "stratified + forced unstructured run completes without aborting")
   end subroutine test_unstr_stratified

end module test_driver_unstr
