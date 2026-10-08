!! Bundled test driver — fpm gcc build entry point.
!!
!! Compiles every `tests/test_*.F90` module into one executable that
!! runs every testdrive suite in sequence.
!! CMake builds the same modules as N standalone executables (one per
!! suite) via a `configure_file` template in `tests/CMakeLists.txt`, so
!! per-test isolation under gdb / sanitizers is preserved on that path.
!!
!! Each `use` line renames `collect_<name>_tests` to a unique local
!! identifier (suffix `_collect`) so the bundled program doesn't trip
!! Fortran's name-collision rules when two test modules export the same
!! `collect_*` symbol (e.g. `collect_dam_break_tests` exists in both
!! `test_dam_break` and `test_dam_break_unstr`).
program rakali_tests
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mem_report, only: mem_report_collect => collect_mem_report_tests
   use test_efp, only: efp_collect => collect_efp_tests
   use test_boundary, only: boundary_collect => collect_boundary_tests
   use test_convergence, only: convergence_collect => collect_convergence_tests
   use test_cross_backend, only: cross_backend_collect => collect_cross_tests
#ifndef RAKALI_NO_NETCDF
#endif
   use test_safe_math, only: safe_math_collect => collect_safe_math_tests
   use test_console_stats, only: console_stats_collect => collect_console_stats_tests
#ifndef RAKALI_NO_NETCDF
#endif
   use test_dam_break_unstr, only: dam_break_unstr_collect => collect_dam_break_tests
   use test_dam_break, only: dam_break_collect => collect_dam_break_tests
   use test_dam_break_dry, only: dam_break_dry_collect => collect_dam_break_dry_tests
   use test_carrier_greenspan, only: carrier_greenspan_collect => collect_carrier_greenspan_tests
   use test_driver_smoke, only: driver_smoke_collect => collect_driver_smoke_tests
   use test_driver_coastal, only: driver_coastal_collect => collect_driver_coastal_tests
   use test_extrapolate_unstr, only: extrapolate_unstr_collect => collect_extrapolate_tests
   use test_flux_unstr, only: flux_unstr_collect => collect_flux_tests
   use test_flux, only: flux_collect => collect_flux_tests
   use test_flux_dump, only: flux_dump_collect => collect_flux_dump_tests
   use test_forcing, only: forcing_collect => collect_forcing_tests
   use test_friction_unstr, only: friction_unstr_collect => collect_friction_tests
   use test_friction, only: friction_collect => collect_friction_tests
   use test_gauges, only: gauges_collect => collect_gauges_tests
   use test_gaussian_hump, only: gaussian_hump_collect => collect_gaussian_hump_tests
   use test_init_diagnostics, only: init_diagnostics_collect => collect_init_diagnostics_tests
   use test_lake_at_rest, only: lake_at_rest_collect => collect_lake_at_rest_tests
   use test_c_property, only: c_property_collect => collect_c_property_tests
   use test_bpg_uniform_density, only: bpg_uniform_density_collect => collect_bpg_uniform_density_tests
   use test_convergence_harness, only: convergence_harness_collect => collect_convergence_harness_tests
   use test_convergence_structured, only: convergence_structured_collect => collect_convergence_structured_tests
   use test_tracer_conservation, only: tracer_conservation_collect => collect_tracer_conservation_tests
   use test_conservation_parity, only: conservation_parity_collect => collect_conservation_parity_tests
   use test_coastal_conservation, only: coastal_conservation_collect => collect_coastal_conservation_tests
   use test_mesh_unstr, only: mesh_unstr_collect => collect_mesh_unstr_tests
   use test_mode_split, only: mode_split_collect => collect_mode_split_tests
   use test_multilayer_unstr, only: multilayer_unstr_collect => collect_multilayer_unstr_tests
   use test_ml_trunc_friction_unstr, only: ml_trunc_friction_unstr_collect => &
                                           collect_ml_trunc_friction_unstr_tests
   use test_multilayer, only: multilayer_collect => collect_multilayer_tests
   use test_nesting, only: nesting_collect => collect_nesting_tests
   use test_nested_unstr, only: nested_unstr_collect => collect_nested_unstr_tests
   use test_sponge_unstr, only: sponge_unstr_collect => collect_sponge_unstr_tests
   use test_nonuniform_sigma, only: nonuniform_sigma_collect => collect_nonuniform_sigma_tests
   use test_pressure, only: pressure_collect => collect_pressure_tests
   use test_protect_unstr, only: protect_unstr_collect => collect_protect_tests
   use test_remap, only: remap_collect => collect_remap_tests
   use test_remap_pqm, only: remap_pqm_collect => collect_remap_pqm_tests
   use test_rki_api, only: rki_api_collect => collect_api_tests
   use test_rki_handle, only: rki_handle_collect => collect_handle_tests
   use test_solver_unstr, only: solver_unstr_collect => collect_solver_tests
   use test_solver, only: solver_collect => collect_solver_tests
   use test_sponge, only: sponge_collect => collect_sponge_tests
   use test_state_unstr, only: state_unstr_collect => collect_state_unstr_tests
   use test_state, only: state_collect => collect_state_tests
   use test_thacker_bowl, only: thacker_bowl_collect => collect_thacker_bowl_tests
   use test_time_schemes, only: time_schemes_collect => collect_time_schemes_tests
   use test_time_schemes_unstr, only: time_schemes_unstr_collect => collect_time_schemes_unstr_tests
   use test_timestep_unstr, only: timestep_unstr_collect => collect_timestep_tests
   use test_timestep, only: timestep_collect => collect_timestep_tests
   use test_update_unstr, only: update_unstr_collect => collect_update_tests
   use test_update, only: update_collect => collect_update_tests
   use test_vcoord, only: vcoord_collect => collect_vcoord_tests
   use test_wetdry, only: wetdry_collect => collect_wetdry_tests
   use test_2d_tracer, only: ts_2d_tracer_collect => collect_2d_tracer_tests
   use test_lock_exchange, only: lock_exchange_collect => collect_lock_exchange_tests
   use test_multilayer_drycell, only: multilayer_drycell_collect => collect_drycell_tests
   use test_river_plume, only: river_plume_collect => collect_river_plume_tests
   use test_coastal_kpp, only: coastal_kpp_collect => collect_coastal_kpp_tests
   use test_coastal_evap, only: coastal_evap_collect => collect_coastal_evap_tests
   use test_solitary_wave, only: solitary_wave_collect => collect_solitary_wave_tests
   use test_thermal_diffusion, only: thermal_diffusion_collect => collect_thermal_diffusion_tests
   use test_profiler, only: profiler_collect => collect_profiler_tests
   use test_init, only: init_collect => collect_init_tests
   use test_ml_eos, only: ml_eos_collect => collect_ml_eos_tests
   use test_coastal_shared_eos, only: coastal_shared_eos_collect => collect_coastal_shared_eos_tests
   use test_bathymetry, only: bathymetry_collect => collect_bathymetry_tests
   use test_output, only: output_collect => collect_output_tests
   use test_output_unstr, only: output_unstr_collect => collect_output_unstr_tests
   use test_gauges_unstr, only: gauges_unstr_collect => collect_gauges_unstr_tests
   use test_restart_unstr, only: restart_unstr_collect => collect_restart_unstr_tests
   use test_boundary_nh, only: boundary_nh_collect => collect_boundary_nh_tests
   use test_state_init_from_config, only: state_init_from_config_collect => collect_state_init_from_config_tests
   use test_vcoord_target, only: vcoord_target_collect => collect_vcoord_target_tests
   use test_vcoord_zstar_full, only: vcoord_zstar_full_collect => collect_vcoord_zstar_full_tests
   use test_zstar_full_shelf_break, only: zstar_full_shelf_break_collect => collect_zstar_full_shelf_break_tests
   use test_zstar_full_shelf_break_unstr, only: zstar_full_shelf_break_unstr_collect => collect_zstar_full_shelf_break_unstr_tests
   use test_zstar_full_tracer_conservation, only: zstar_full_tracer_conservation_collect => collect_zstar_full_tracer_conservation_tests
   use test_ml_tracer_advect_race, only: ml_tracer_advect_race_collect => collect_ml_tracer_advect_race_tests
   use test_sediment, only: sediment_collect => collect_sediment_tests
   use test_ml_dynamics_unstr, only: ml_dynamics_unstr_collect => collect_ml_dynamics_unstr_tests
   use test_coastal_kpp_unstr, only: coastal_kpp_unstr_collect => collect_coastal_kpp_unstr_tests
   use test_kernel_fused, only: kernel_fused_collect => collect_kernel_fused_tests
   use test_kernel_fused_unstr, only: kernel_fused_unstr_collect => collect_kernel_fused_unstr_tests
   implicit none

   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mem_report", mem_report_collect), &
                new_testsuite("efp", efp_collect), &
                new_testsuite("boundary", boundary_collect), &
                new_testsuite("convergence", convergence_collect), &
                new_testsuite("cross_backend", cross_backend_collect), &
#endif
                new_testsuite("safe_math", safe_math_collect), &
                new_testsuite("console_stats", console_stats_collect), &
#ifndef RAKALI_NO_NETCDF
#endif
                new_testsuite("dam_break_unstr", dam_break_unstr_collect), &
                new_testsuite("dam_break", dam_break_collect), &
                new_testsuite("dam_break_dry", dam_break_dry_collect), &
                new_testsuite("carrier_greenspan", carrier_greenspan_collect), &
                new_testsuite("driver_smoke", driver_smoke_collect), &
                new_testsuite("driver_coastal", driver_coastal_collect), &
                new_testsuite("extrapolate_unstr", extrapolate_unstr_collect), &
                new_testsuite("flux_unstr", flux_unstr_collect), &
                new_testsuite("flux", flux_collect), &
                new_testsuite("flux_dump", flux_dump_collect), &
                new_testsuite("forcing", forcing_collect), &
                new_testsuite("friction_unstr", friction_unstr_collect), &
                new_testsuite("friction", friction_collect), &
                new_testsuite("gauges", gauges_collect), &
                new_testsuite("gaussian_hump", gaussian_hump_collect), &
                new_testsuite("init_diagnostics", init_diagnostics_collect), &
                new_testsuite("lake_at_rest", lake_at_rest_collect), &
                new_testsuite("c_property", c_property_collect), &
                new_testsuite("bpg_uniform_density", bpg_uniform_density_collect), &
                new_testsuite("convergence_harness", convergence_harness_collect), &
                new_testsuite("convergence_structured", convergence_structured_collect), &
                new_testsuite("tracer_conservation", tracer_conservation_collect), &
                new_testsuite("conservation_parity", conservation_parity_collect), &
                new_testsuite("coastal_conservation", coastal_conservation_collect), &
                new_testsuite("mesh_unstr", mesh_unstr_collect), &
                new_testsuite("mode_split", mode_split_collect), &
                new_testsuite("multilayer_unstr", multilayer_unstr_collect), &
                new_testsuite("ml_trunc_friction_unstr", ml_trunc_friction_unstr_collect), &
                new_testsuite("multilayer", multilayer_collect), &
                new_testsuite("nesting", nesting_collect), &
                new_testsuite("nested_unstr", nested_unstr_collect), &
                new_testsuite("sponge_unstr", sponge_unstr_collect), &
                new_testsuite("nonuniform_sigma", nonuniform_sigma_collect), &
                new_testsuite("pressure", pressure_collect), &
                new_testsuite("protect_unstr", protect_unstr_collect), &
                new_testsuite("remap", remap_collect), &
                new_testsuite("remap_pqm", remap_pqm_collect), &
                new_testsuite("rki_api", rki_api_collect), &
                new_testsuite("rki_handle", rki_handle_collect), &
                new_testsuite("solver_unstr", solver_unstr_collect), &
                new_testsuite("solver", solver_collect), &
                new_testsuite("sponge", sponge_collect), &
                new_testsuite("state_unstr", state_unstr_collect), &
                new_testsuite("state", state_collect), &
                new_testsuite("thacker_bowl", thacker_bowl_collect), &
                new_testsuite("time_schemes", time_schemes_collect), &
                new_testsuite("time_schemes_unstr", time_schemes_unstr_collect), &
                new_testsuite("timestep_unstr", timestep_unstr_collect), &
                new_testsuite("timestep", timestep_collect), &
                new_testsuite("update_unstr", update_unstr_collect), &
                new_testsuite("update", update_collect), &
                new_testsuite("vcoord", vcoord_collect), &
                new_testsuite("wetdry", wetdry_collect), &
                new_testsuite("2d_tracer", ts_2d_tracer_collect), &
                new_testsuite("lock_exchange", lock_exchange_collect), &
                new_testsuite("multilayer_drycell", multilayer_drycell_collect), &
                new_testsuite("river_plume", river_plume_collect), &
                new_testsuite("coastal_kpp", coastal_kpp_collect), &
                new_testsuite("coastal_evap", coastal_evap_collect), &
                new_testsuite("solitary_wave", solitary_wave_collect), &
                new_testsuite("thermal_diffusion", thermal_diffusion_collect), &
                new_testsuite("profiler", profiler_collect), &
                new_testsuite("init", init_collect), &
                new_testsuite("ml_eos", ml_eos_collect), &
                new_testsuite("coastal_shared_eos", coastal_shared_eos_collect), &
                new_testsuite("bathymetry", bathymetry_collect), &
                new_testsuite("output", output_collect), &
                new_testsuite("output_unstr", output_unstr_collect), &
                new_testsuite("gauges_unstr", gauges_unstr_collect), &
                new_testsuite("restart_unstr", restart_unstr_collect), &
                new_testsuite("boundary_nh", boundary_nh_collect), &
                new_testsuite("state_init_from_config", state_init_from_config_collect), &
                new_testsuite("vcoord_target", vcoord_target_collect), &
                new_testsuite("vcoord_zstar_full", vcoord_zstar_full_collect), &
                new_testsuite("zstar_full_shelf_break", zstar_full_shelf_break_collect), &
                new_testsuite("zstar_full_shelf_break_unstr", zstar_full_shelf_break_unstr_collect), &
                new_testsuite("zstar_full_tracer_conservation", zstar_full_tracer_conservation_collect), &
                new_testsuite("ml_tracer_advect_race", ml_tracer_advect_race_collect), &
                new_testsuite("sediment", sediment_collect), &
                new_testsuite("ml_dynamics_unstr", ml_dynamics_unstr_collect), &
                new_testsuite("coastal_kpp_unstr", coastal_kpp_unstr_collect), &
                new_testsuite("kernel_fused", kernel_fused_collect), &
                new_testsuite("kernel_fused_unstr", kernel_fused_unstr_collect) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop 1
   end if

end program rakali_tests
