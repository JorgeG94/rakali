!! Unit tests for state_init_from_config — covers the multilayer + NH
!! configuration paths in rki_state.F90 that are not exercised by the
!! plain `state%init` constructor.
module test_state_init_from_config
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_config, only: config_t
   implicit none
   private

   public :: collect_state_init_from_config_tests

contains

   subroutine collect_state_init_from_config_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("barotropic_only", test_barotropic_only), &
                  new_unittest("multilayer_path", test_multilayer_path), &
                  new_unittest("multilayer_registers_tracers", test_register_tracers), &
                  new_unittest("nonhydrostatic_path", test_nh_path), &
                  new_unittest("physics_scalars_pulled_from_config", test_physics_scalars), &
                  new_unittest("forcing_scalars_uniform_fill", test_forcing_fill) &
                  ]
   end subroutine collect_state_init_from_config_tests

   subroutine test_barotropic_only(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      checks: block

         cfg%nx = 10; cfg%ny = 8; cfg%nghost = 2
         cfg%dx = 1.0_wp; cfg%dy = 1.0_wp; cfg%cfl = 0.45_wp

         call s%init_from_config(cfg)

         call check(error, s%grid%nx_phys == 10, "nx_phys should match cfg")
         if (allocated(error)) exit checks
         call check(error, s%grid%ny_phys == 8, "ny_phys should match cfg")
         if (allocated(error)) exit checks
         call check(error,.not. s%use_multilayer, "barotropic-only config")
         if (allocated(error)) exit checks
         call check(error,.not. s%use_nonhydrostatic, "barotropic-only config (NH off)")
         if (allocated(error)) exit checks

      end block checks
      continue
      call s%destroy()
   end subroutine test_barotropic_only

   subroutine test_multilayer_path(error)
      !! Drive the use_multilayer branch — populates vcoord + EOS scalars.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      checks: block

         cfg%nx = 6; cfg%ny = 4; cfg%nghost = 2
         cfg%dx = 1.0_wp; cfg%dy = 1.0_wp; cfg%cfl = 0.45_wp
         cfg%use_multilayer = .true.
         cfg%nz_layers = 4
         cfg%vcoord_type = "sigma"
         cfg%remap_method = "ppm"
         cfg%initial_salinity = 34.5_wp
         cfg%S_ref = 0.0_wp
         cfg%beta_S = 0.78_wp
         cfg%initial_temperature = 12.5_wp
         cfg%T_ref = 12.5_wp
         cfg%alpha_T = 0.17_wp
         cfg%rho_0 = 1025.0_wp
         cfg%bpg_method = "jacobian"
         cfg%use_ppm_tracer = .true.

         call s%init_from_config(cfg)

         call check(error, s%use_multilayer, "use_multilayer flag set")
         if (allocated(error)) exit checks
         call check(error, s%multilayer%nz_ml == 4, "nz_ml plumbed through")
         if (allocated(error)) exit checks
         call check(error, abs(s%multilayer%initial_salinity - 34.5_wp) < 1.0e-12_wp, &
                    "initial_salinity plumbed through")
         if (allocated(error)) exit checks
         call check(error, abs(s%multilayer%rho_0 - 1025.0_wp) < 1.0e-12_wp, &
                    "rho_0 plumbed through")
         if (allocated(error)) exit checks
         call check(error, s%multilayer%use_bpg_jacobian, "bpg_method=jacobian")
         if (allocated(error)) exit checks
         call check(error, s%multilayer%use_ppm_tracer, "use_ppm_tracer plumbed through")
         if (allocated(error)) exit checks

      end block checks
      continue
      call s%destroy()
   end subroutine test_multilayer_path

   subroutine test_register_tracers(error)
      !! register_default_tracers should overlay tracer metadata from cfg.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      integer :: idx_S, idx_T
      checks: block

         cfg%nx = 4; cfg%ny = 4; cfg%nghost = 2
         cfg%dx = 1.0_wp; cfg%dy = 1.0_wp; cfg%cfl = 0.45_wp
         cfg%use_multilayer = .true.
         cfg%nz_layers = 2
         cfg%initial_salinity = 32.0_wp
         cfg%initial_temperature = 18.0_wp
         cfg%S_ref = 0.0_wp
         cfg%T_ref = 18.0_wp
         cfg%beta_S = 0.78_wp
         cfg%alpha_T = 0.17_wp

         call s%init_from_config(cfg)

         idx_S = s%multilayer%idx_salinity
         idx_T = s%multilayer%idx_temperature

         call check(error, trim(s%multilayer%tracers(idx_S)%name) == "salinity", &
                    "salinity tracer name registered")
         if (allocated(error)) exit checks
         call check(error, trim(s%multilayer%tracers(idx_T)%name) == "temperature", &
                    "temperature tracer name registered")
         if (allocated(error)) exit checks
         call check(error, abs(s%multilayer%tracers(idx_S)%tr_init - 32.0_wp) < 1.0e-12_wp, &
                    "salinity tr_init plumbed from cfg")
         if (allocated(error)) exit checks
         call check(error, abs(s%multilayer%tracers(idx_T)%tr_init - 18.0_wp) < 1.0e-12_wp, &
                    "temperature tr_init plumbed from cfg")
         if (allocated(error)) exit checks
         call check(error, abs(s%multilayer%tracers(idx_T)%eos_coeff - (-0.17_wp)) < 1.0e-12_wp, &
                    "temperature eos_coeff is -alpha_T")
         if (allocated(error)) exit checks

      end block checks
      continue
      call s%destroy()
   end subroutine test_register_tracers

   subroutine test_nh_path(error)
      !! Drive the use_nonhydrostatic branch — populates NH config scalars.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      checks: block

         cfg%nx = 6; cfg%ny = 4; cfg%nghost = 2
         cfg%dx = 1.0_wp; cfg%dy = 1.0_wp; cfg%cfl = 0.45_wp
         cfg%use_nonhydrostatic = .true.
         cfg%nz_layers = 5
         cfg%nh_tolerance = 1.0e-7_wp
         cfg%nh_max_iterations = 200
         cfg%split_ratio = 4

         call s%init_from_config(cfg)

         call check(error, s%use_nonhydrostatic, "NH flag set")
         if (allocated(error)) exit checks
         call check(error, s%nh%nz == 5, "NH nz plumbed")
         if (allocated(error)) exit checks
         call check(error, abs(s%nh%nh_tolerance - 1.0e-7_wp) < 1.0e-15_wp, "nh_tolerance plumbed")
         if (allocated(error)) exit checks
         call check(error, s%nh%nh_max_iterations == 200, "nh_max_iterations plumbed")
         if (allocated(error)) exit checks
         call check(error, s%nh%split_ratio == 4, "split_ratio plumbed")
         if (allocated(error)) exit checks

      end block checks
      continue
      call s%destroy()
   end subroutine test_nh_path

   subroutine test_physics_scalars(error)
      !! Manning + Coriolis + HLLC flag should round-trip.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      checks: block

         cfg%nx = 4; cfg%ny = 4; cfg%nghost = 2
         cfg%dx = 1.0_wp; cfg%dy = 1.0_wp; cfg%cfl = 0.45_wp
         cfg%manning_n = 0.025_wp
         cfg%coriolis_f = 1.0e-4_wp
         cfg%use_hllc = .true.
         cfg%dt_fixed = 0.1_wp
         cfg%cfl_interval = 5

         call s%init_from_config(cfg)

         call check(error, abs(s%barotropic%manning_n - 0.025_wp) < 1.0e-12_wp, "manning_n plumbed")
         if (allocated(error)) exit checks
         call check(error, abs(s%barotropic%coriolis_f - 1.0e-4_wp) < 1.0e-15_wp, "coriolis_f plumbed")
         if (allocated(error)) exit checks
         call check(error, s%barotropic%use_hllc, "use_hllc plumbed")
         if (allocated(error)) exit checks
         call check(error, abs(s%dt_fixed - 0.1_wp) < 1.0e-12_wp, "dt_fixed plumbed")
         if (allocated(error)) exit checks
         call check(error, s%cfl_interval == 5, "cfl_interval plumbed")
         if (allocated(error)) exit checks

      end block checks
      continue
      call s%destroy()
   end subroutine test_physics_scalars

   subroutine test_forcing_fill(error)
      !! Uniform wind stress should be broadcast to every cell.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      checks: block

         cfg%nx = 4; cfg%ny = 4; cfg%nghost = 2
         cfg%dx = 1.0_wp; cfg%dy = 1.0_wp; cfg%cfl = 0.45_wp
         cfg%wind_stress_x = 0.5_wp
         cfg%wind_stress_y = -0.25_wp

         call s%init_from_config(cfg)

         call check(error, all(abs(s%forcing%wind_x - 0.5_wp) < 1.0e-12_wp), &
                    "wind_x array set to scalar")
         if (allocated(error)) exit checks
         call check(error, all(abs(s%forcing%wind_y - (-0.25_wp)) < 1.0e-12_wp), &
                    "wind_y array set to scalar")
         if (allocated(error)) exit checks

      end block checks
      continue
      call s%destroy()
   end subroutine test_forcing_fill

end module test_state_init_from_config
