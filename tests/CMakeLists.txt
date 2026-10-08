# --------------------------------------------------------------------------
# Test build
#
# Each test source file (`tests/test_*.F90`) is *module-only* — it exports a
# `collect_*_tests` subroutine for the testdrive framework. The matching program
# driver is auto-generated from `tester_template.F90.in` so we preserve the
# per-test executable model (one binary per suite, easy to run under gdb /
# sanitizers / `ctest -R`) without checking in 48 nearly identical driver files.
#
# The same modules are bundled into a single `fpm test` executable via
# `tests/tester.F90`. The flat `tests/` layout is required because fpm's
# `[[test]]` source-dir doesn't recurse subdirs (verified against fpm 0.13).
# Keep the two paths in sync when adding a new test: - drop the module file in
# `tests/` - append it to `RAKALI_TESTS` below (CMake) - append it to
# `tests/tester.F90` (fpm)
# --------------------------------------------------------------------------

# Helper modules used by some unit tests (Ritter dam-break + Thacker bowl
# analytical solutions). Built as a static lib and linked into every test.
add_library(
  rakali_analytical STATIC ritter_dambreak.F90 thacker_bowl.F90
                           convergence_harness.F90 carrier_greenspan.F90)
set_target_properties(rakali_analytical PROPERTIES Fortran_MODULE_DIRECTORY
                                                   ${CMAKE_BINARY_DIR}/modules)
target_include_directories(rakali_analytical
                           PRIVATE "${CMAKE_BINARY_DIR}/modules")
# `ocean_test_metrics` uses `rki_ocean_metrics`, whose own `use` chain reaches
# `rki_io_netcdf` — so the compiler must open `netcdf.mod` while walking the
# .mod tree, even though nothing here calls NetCDF. Every test target picks the
# dir up transitively from ${main_lib}; this target links no core (see below),
# so take the include dirs alone and keep the no-link property intact.
if(RAKALI_ENABLE_NETCDF)
  target_include_directories(
    rakali_analytical
    PRIVATE $<TARGET_PROPERTY:netCDF::netcdff,INTERFACE_INCLUDE_DIRECTORIES>)
endif()
# Same reasoning for pic. `rki_constants` re-exports `pic_types`, so a compiler
# that opens the TRANSITIVE .mod files while walking the chain needs pic's
# module dir on the include path -- nvfortran does, and fails with "Unable to
# open MODULE file pic_types.mod"; gfortran does not, because its .mod files are
# self-contained. Newer pic emits its .mod files into its own
# _deps/pic-build/modules rather than ${CMAKE_BINARY_DIR}/modules. Take the
# include dirs alone, keeping the no-link property the block below depends on.
target_include_directories(
  rakali_analytical
  PRIVATE $<TARGET_PROPERTY:pic::pic,INTERFACE_INCLUDE_DIRECTORIES>)
# Depend on the core for build ordering + the shared .mod dir, but do NOT LINK
# it: a static lib propagates its PRIVATE link deps to consumers, so linking the
# static core here would drag it into EVERY test and defeat the shared-core
# split below. Each test links its own core; analytical's core symbols resolve
# against whichever one the test picks. (Archiving a static lib does no symbol
# resolution, so the missing link is harmless at this stage.)
#
# THIS TARGET IS WHY THE BUILD REQUIRES THE MAKEFILES GENERATOR.  These sources
# `use rki_constants`, so they must not compile until the core's .mod files
# exist -- but `add_dependencies` does NOT express that.  It only chains the
# object-ORDER phonies (cmake_object_order_depends_target_*), which never depend
# on the core's actual outputs.  Makefiles happens to serialise it anyway
# because make recurses directories in order; Ninja's single global graph does
# not, so the analytical objects race ahead of the core and the build dies on
# "Unable to open MODULE file rki_constants.mod".
#
# It is NOT fixable by tightening the dependency: CMake feeds the Ninja Fortran
# dyndep scanner only the module dirs of targets this one LINKS, so
# `add_dependencies` and `$<COMPILE_ONLY:${core_objs}>` are both invisible to it
# (verified: neither produces a .tdi entry, both still race).  Making it
# Ninja-clean means giving this target a real link edge to a core, which is
# exactly what the paragraph above forbids -- a test that links the shared core
# would then drag the static one in too.  Switching generators therefore needs
# that no-link constraint revisited first; don't just flip CMAKE_GENERATOR.
add_dependencies(rakali_analytical ${main_lib})

# Per-test metadata.
#
# Each row encodes one test suite as a single string, fields separated by `|`:
# <test-name>|<module-name>|<collect-fn>|<suite-display>|<regime>
#
# Where: <test-name>     — file stem (also the ctest name and executable name)
# <module-name>   — Fortran module name to `use` (matches <test-name> in every
# current row; column kept so a future test with a different module name can
# plug in trivially) <collect-fn>    — testdrive collector subroutine in that
# module <suite-display> — string passed to testdrive's `new_testsuite()`
# <regime>        — ctest LABEL for the regime the test belongs to: `coastal`
# (A-grid/structured/unstructured/multilayer path), `ocean` (C-grid dyn-core +
# sea-ice), `core` (shared infrastructure), or `ocean+core` (both labels).
# Mandatory — configure fails loud on a missing/invalid regime, and the
# `test-regime-labels` pre-commit hook enforces it too.  Run a single regime
# suite with `ctest -L coastal` / `-L ocean` / `-L core`.
set(RAKALI_TESTS
    "test_mem_report|test_mem_report|collect_mem_report_tests|mem_report|core"
    "test_efp|test_efp|collect_efp_tests|efp|core"
    "test_console_stats|test_console_stats|collect_console_stats_tests|console_stats|coastal"
    "test_boundary|test_boundary|collect_boundary_tests|boundary|coastal"
    "test_nml_schema|test_nml_schema|collect_nml_schema_tests|nml_schema|core"
    "test_config_schema|test_config_schema|collect_config_schema_tests|config_schema|core"
    "test_convergence|test_convergence|collect_convergence_tests|convergence|coastal"
    "test_safe_math|test_safe_math|collect_safe_math_tests|safe_math|core"
    # Split per-subtest (same source) so the heavy 20000-step geostrophic run is
    # its own ctest binary, focusable/skippable in isolation.
    "test_tracer_recon|test_tracer_recon|collect_tracer_recon_tests|tracer_recon|coastal"
    "test_dc_reduce|test_dc_reduce|collect_dc_reduce_tests|dc_reduce|core"
    "test_dc_transfer_sentinel|test_dc_transfer_sentinel|collect_dc_transfer_sentinel_tests|dc_transfer_sentinel|core"
    "test_cross_backend|test_cross_backend|collect_cross_tests|cross_backend|coastal"
    "test_dam_break|test_dam_break|collect_dam_break_tests|dam_break|coastal"
    "test_dam_break_unstr|test_dam_break_unstr|collect_dam_break_tests|dam_break_unstr|coastal"
    "test_dam_break_dry|test_dam_break_dry|collect_dam_break_dry_tests|dam_break_dry|coastal"
    "test_carrier_greenspan|test_carrier_greenspan|collect_carrier_greenspan_tests|carrier_greenspan|coastal"
    "test_driver_smoke|test_driver_smoke|collect_driver_smoke_tests|driver_smoke|coastal"
    "test_driver_coastal|test_driver_coastal|collect_driver_coastal_tests|driver_coastal|coastal"
    "test_driver_unstr|test_driver_unstr|collect_driver_unstr_tests|driver_unstr|coastal"
    "test_extrapolate_unstr|test_extrapolate_unstr|collect_extrapolate_tests|extrapolate_unstr|coastal"
    "test_flux|test_flux|collect_flux_tests|flux|coastal"
    "test_flux_dump|test_flux_dump|collect_flux_dump_tests|flux_dump|coastal"
    "test_flux_unstr|test_flux_unstr|collect_flux_tests|flux_unstr|coastal"
    "test_forcing|test_forcing|collect_forcing_tests|forcing|core"
    "test_friction|test_friction|collect_friction_tests|friction|coastal"
    "test_friction_unstr|test_friction_unstr|collect_friction_tests|friction_unstr|coastal"
    "test_biharmonic|test_biharmonic|collect_biharmonic_tests|biharmonic|coastal"
    "test_keps|test_keps|collect_keps_tests|keps|coastal"
    "test_keps_unstr|test_keps_unstr|collect_keps_unstr_tests|keps_unstr|coastal"
    "test_keps_kato_unstr|test_keps_kato_unstr|collect_keps_kato_unstr_tests|keps_kato_unstr|coastal"
    "test_kato_phillips|test_kato_phillips|collect_kato_phillips_tests|kato_phillips|coastal"
    "test_device_roundtrip|test_device_roundtrip|collect_device_roundtrip_tests|device_roundtrip|coastal"
    "test_gauges|test_gauges|collect_gauges_tests|gauges|core"
    "test_gaussian_hump|test_gaussian_hump|collect_gaussian_hump_tests|gaussian_hump|coastal"
    "test_init_diagnostics|test_init_diagnostics|collect_init_diagnostics_tests|init_diagnostics|core"
    "test_lake_at_rest|test_lake_at_rest|collect_lake_at_rest_tests|lake_at_rest|coastal"
    "test_c_property|test_c_property|collect_c_property_tests|c_property|coastal"
    "test_bpg_uniform_density|test_bpg_uniform_density|collect_bpg_uniform_density_tests|bpg_uniform_density|coastal"
    "test_convergence_harness|test_convergence_harness|collect_convergence_harness_tests|convergence_harness|core"
    "test_convergence_structured|test_convergence_structured|collect_convergence_structured_tests|convergence_structured|coastal"
    "test_tracer_conservation|test_tracer_conservation|collect_tracer_conservation_tests|tracer_conservation|coastal"
    "test_conservation_parity|test_conservation_parity|collect_conservation_parity_tests|conservation_parity|coastal"
    "test_coastal_conservation|test_coastal_conservation|collect_coastal_conservation_tests|coastal_conservation|coastal"
    "test_si_unstr_seiche|test_si_unstr_seiche|collect_si_unstr_seiche_tests|si_unstr_seiche|coastal"
    "test_si_unstr_channel|test_si_unstr_channel|collect_si_unstr_channel_tests|si_unstr_channel|coastal"
    "test_si_unstr_mg|test_si_unstr_mg|collect_si_unstr_mg_tests|si_unstr_mg|coastal"
    "test_si_unstr_ml_baro_limit|test_si_unstr_ml_baro_limit|collect_si_unstr_ml_baro_limit_tests|si_unstr_ml_baro_limit|coastal"
    "test_si_unstr_ml_tidal_tracer|test_si_unstr_ml_tidal_tracer|collect_si_unstr_ml_tidal_tracer_tests|si_unstr_ml_tidal_tracer|coastal"
    "test_si_unstr_wetdry|test_si_unstr_wetdry|collect_si_unstr_wetdry_tests|si_unstr_wetdry|coastal"
    "test_si_unstr_ml_wetdry|test_si_unstr_ml_wetdry|collect_si_unstr_ml_wetdry_tests|si_unstr_ml_wetdry|coastal"
    "test_si_wetdry_bsub|test_si_wetdry_bsub|collect_si_wetdry_bsub_tests|si_wetdry_bsub|coastal"
    "test_si_unstr_estuary_bench|test_si_unstr_estuary_bench|collect_si_unstr_estuary_bench_tests|si_unstr_estuary_bench|coastal"
    "test_si_closures_unstr|test_si_closures_unstr|collect_si_closures_unstr_tests|si_closures_unstr|coastal"
    "test_mesh_unstr|test_mesh_unstr|collect_mesh_unstr_tests|mesh_unstr|coastal"
    "test_mesh_gen_unstr|test_mesh_gen_unstr|collect_mesh_gen_unstr_tests|mesh_gen_unstr|coastal"
    "test_mode_split|test_mode_split|collect_mode_split_tests|mode_split|coastal"
    "test_multilayer|test_multilayer|collect_multilayer_tests|multilayer|coastal"
    "test_multilayer_unstr|test_multilayer_unstr|collect_multilayer_unstr_tests|multilayer_unstr|coastal"
    "test_nesting|test_nesting|collect_nesting_tests|nesting|coastal"
    "test_nonuniform_sigma|test_nonuniform_sigma|collect_nonuniform_sigma_tests|nonuniform_sigma|coastal"
    "test_pressure|test_pressure|collect_pressure_tests|pressure|coastal"
    "test_protect_unstr|test_protect_unstr|collect_protect_tests|protect_unstr|coastal"
    "test_remap|test_remap|collect_remap_tests|remap|coastal"
    "test_remap_pqm|test_remap_pqm|collect_remap_pqm_tests|remap_pqm|coastal"
    # Per-test isolation for the rki_api tests — each create/destroy of a solver
    # handle leaves residual state in the OpenMP-target present table on Aurora
    # ifx that poisons the next test's enter_data. One ctest binary per scenario
    # keeps each in a fresh process. The ${MODULE_NAME}.F90 substitution below
    # lets four binaries share tests/test_rki_api.F90 as their source.
    "test_rki_api_oneshot|test_rki_api|collect_api_oneshot_test|rki_api_oneshot|core"
    "test_rki_api_twophase|test_rki_api|collect_api_twophase_test|rki_api_twophase|core"
    "test_rki_api_two_handles|test_rki_api|collect_api_two_handles_test|rki_api_two_handles|core"
    "test_rki_api_bad_handle|test_rki_api|collect_api_bad_handle_test|rki_api_bad_handle|core"
    "test_rki_api_string|test_rki_api|collect_api_string_test|rki_api_string|core"
    "test_rki_handle|test_rki_handle|collect_handle_tests|rki_handle|core"
    "test_solver|test_solver|collect_solver_tests|solver|coastal"
    "test_solver_unstr|test_solver_unstr|collect_solver_tests|solver_unstr|coastal"
    "test_sponge|test_sponge|collect_sponge_tests|sponge|coastal"
    "test_sponge_unstr|test_sponge_unstr|collect_sponge_unstr_tests|sponge_unstr|coastal"
    "test_nested_unstr|test_nested_unstr|collect_nested_unstr_tests|nested_unstr|coastal"
    "test_state|test_state|collect_state_tests|state|coastal"
    "test_state_unstr|test_state_unstr|collect_state_unstr_tests|state_unstr|coastal"
    "test_thacker_bowl|test_thacker_bowl|collect_thacker_bowl_tests|thacker_bowl|coastal"
    "test_time_schemes|test_time_schemes|collect_time_schemes_tests|time_schemes|coastal"
    "test_time_schemes_unstr|test_time_schemes_unstr|collect_time_schemes_unstr_tests|time_schemes_unstr|coastal"
    "test_timestep|test_timestep|collect_timestep_tests|timestep|coastal"
    "test_timestep_unstr|test_timestep_unstr|collect_timestep_tests|timestep_unstr|coastal"
    "test_update|test_update|collect_update_tests|update|coastal"
    "test_update_unstr|test_update_unstr|collect_update_tests|update_unstr|coastal"
    "test_vcoord|test_vcoord|collect_vcoord_tests|vcoord|core"
    "test_wetdry|test_wetdry|collect_wetdry_tests|wetdry|coastal"
    "test_2d_tracer|test_2d_tracer|collect_2d_tracer_tests|2d_tracer|coastal"
    "test_lock_exchange|test_lock_exchange|collect_lock_exchange_tests|lock_exchange|coastal"
    "test_multilayer_drycell|test_multilayer_drycell|collect_drycell_tests|multilayer_drycell|coastal"
    "test_river_plume|test_river_plume|collect_river_plume_tests|river_plume|coastal"
    "test_solitary_wave|test_solitary_wave|collect_solitary_wave_tests|solitary_wave|coastal"
    "test_thermal_diffusion|test_thermal_diffusion|collect_thermal_diffusion_tests|thermal_diffusion|coastal"
    "test_profiler|test_profiler|collect_profiler_tests|profiler|core"
    "test_init|test_init|collect_init_tests|init|core"
    "test_ml_eos|test_ml_eos|collect_ml_eos_tests|ml_eos|coastal"
    "test_coastal_shared_eos|test_coastal_shared_eos|collect_coastal_shared_eos_tests|coastal_shared_eos|coastal"
    "test_bathymetry|test_bathymetry|collect_bathymetry_tests|bathymetry|core"
    "test_output|test_output|collect_output_tests|output|core"
    "test_output_unstr|test_output_unstr|collect_output_unstr_tests|output_unstr|coastal"
    "test_gauges_unstr|test_gauges_unstr|collect_gauges_unstr_tests|gauges_unstr|coastal"
    "test_restart_unstr|test_restart_unstr|collect_restart_unstr_tests|restart_unstr|coastal"
    "test_boundary_nh|test_boundary_nh|collect_boundary_nh_tests|boundary_nh|coastal"
    "test_state_init_from_config|test_state_init_from_config|collect_state_init_from_config_tests|state_init_from_config|core"
    "test_vcoord_target|test_vcoord_target|collect_vcoord_target_tests|vcoord_target|core"
    "test_vcoord_zstar_full|test_vcoord_zstar_full|collect_vcoord_zstar_full_tests|vcoord_zstar_full|core"
    "test_zstar_full_shelf_break|test_zstar_full_shelf_break|collect_zstar_full_shelf_break_tests|zstar_full_shelf_break|coastal"
    "test_zstar_full_shelf_break_unstr|test_zstar_full_shelf_break_unstr|collect_zstar_full_shelf_break_unstr_tests|zstar_full_shelf_break_unstr|coastal"
    # Per-test isolation — module-scope kernel workspaces (PPM/hdiff scratch)
    # accumulate state across these four tests and the host-side allocate trips
    # "array already allocated" on the second test. One ctest binary per
    # scenario sidesteps it.
    "test_zstar_full_seiche|test_zstar_full_tracer_conservation|collect_zstar_full_seiche_test|zstar_full_seiche|coastal"
    "test_zstar_full_vdiff|test_zstar_full_tracer_conservation|collect_zstar_full_vdiff_test|zstar_full_vdiff|coastal"
    "test_zstar_full_hdiff|test_zstar_full_tracer_conservation|collect_zstar_full_hdiff_test|zstar_full_hdiff|coastal"
    "test_zstar_full_step_basin|test_zstar_full_tracer_conservation|collect_zstar_full_step_basin_test|zstar_full_step_basin|coastal"
    "test_zstar_full_passive_tracer|test_zstar_full_tracer_conservation|collect_zstar_full_passive_tracer_test|zstar_full_passive_tracer|coastal"
    "test_nz1_passive_tracer|test_zstar_full_tracer_conservation|collect_nz1_passive_tracer_test|nz1_passive_tracer|coastal"
    "test_ml_tracer_advect_race|test_ml_tracer_advect_race|collect_ml_tracer_advect_race_tests|ml_tracer_advect_race|coastal"
    "test_sediment|test_sediment|collect_sediment_tests|sediment|coastal"
    "test_ml_dynamics_unstr|test_ml_dynamics_unstr|collect_ml_dynamics_unstr_tests|ml_dynamics_unstr|coastal"
    "test_ml_trunc_friction_unstr|test_ml_trunc_friction_unstr|collect_ml_trunc_friction_unstr_tests|ml_trunc_friction_unstr|coastal"
    "test_coastal_kpp_unstr|test_coastal_kpp_unstr|collect_coastal_kpp_unstr_tests|coastal_kpp_unstr|coastal"
    "test_kernel_fused|test_kernel_fused|collect_kernel_fused_tests|kernel_fused|coastal"
    "test_kernel_fused_unstr|test_kernel_fused_unstr|collect_kernel_fused_unstr_tests|kernel_fused_unstr|coastal"
    "test_surface_flux_strat|test_surface_flux_strat|collect_surface_flux_strat_tests|surface_flux_strat|coastal"
    "test_coastal_kpp|test_coastal_kpp|collect_coastal_kpp_tests|coastal_kpp|coastal"
    "test_coastal_evap|test_coastal_evap|collect_coastal_evap_tests|coastal_evap|coastal"
    "test_inertial_oscillation|test_inertial_oscillation|collect_inertial_oscillation_tests|inertial_oscillation|coastal"
    # Per-test isolation — same workspace-pollution reason as above.
    "test_zstar_ssh_runs|test_zstar_ssh_tracking|collect_zstar_ssh_runs_test|zstar_ssh_runs|coastal"
    "test_zstar_ssh_ratios|test_zstar_ssh_tracking|collect_zstar_ssh_ratios_test|zstar_ssh_ratios|coastal"
)

# Tests that MUST link the STATIC core.
#
# The shared core (.so) is the default test link (RAKALI_TESTS_LINK_SHARED, on
# by default) because it skips the per-exe static-archive pull + RDC device link
# — the slow part of a GPU test build.
#
# A test needs the static core ONLY when its own `do concurrent` (a device
# kernel) CALLS a core `!$acc routine seq`: that call is a device→device
# reference that can only be resolved by device-linking WITHIN one binary, never
# across a .so boundary. A self-contained `do concurrent`, or a host-level call
# to an acc-seq routine, does NOT trigger it — those link fine against the .so.
#
# Empirically that is just test_remap_pqm (its on-device `do concurrent` run
# calls remap_column_pqm). This list is the authority: if a GPU build's device
# link fails with an unresolved device symbol for some test, add it here. On a
# non-GPU (gfortran/ifx) build the distinction is moot — nothing generates
# device code — so a stale entry only costs that test the shared-link speedup.
set(RAKALI_STATIC_CORE_TESTS test_remap_pqm)

# Tests that pull in NetCDF/HDF5 — must run single-threaded because HDF5 is not
# built thread-safe (OMP_NUM_THREADS=1 prevents crashes when the OpenMP runtime
# is active via -mp=multicore or -fopenmp).
set(NETCDF_TESTS
    test_forcing
    test_nesting
    test_solver
    test_dam_break
    test_dam_break_dry
    test_carrier_greenspan
    test_convergence
    test_lake_at_rest
    test_c_property
    test_bpg_uniform_density
    test_convergence_harness
    test_convergence_structured
    test_tracer_conservation
    test_conservation_parity
    test_coastal_conservation
    test_thacker_bowl
    test_gaussian_hump
    test_mode_split
    test_gauges
    test_solitary_wave
    test_lock_exchange
    test_river_plume
    test_multilayer_drycell
    test_2d_tracer
    test_thermal_diffusion
    test_bathymetry
    test_output
    test_output_unstr
    test_gauges_unstr
    test_restart_unstr
    test_surface_flux_strat
    test_coastal_evap
    test_kato_phillips
    test_device_roundtrip
    test_inertial_oscillation
    test_flux_dump
    test_zstar_ssh_tracking
    test_driver_coastal
    test_driver_unstr
    test_rki_api_oneshot
    test_rki_api_twophase
    test_rki_api_two_handles
    test_rki_api_bad_handle
    test_rki_api_string)

# Tests that directly `use` modules from the NetCDF I/O subsystem (rki_forcing,
# rki_nesting, rki_gauges, rki_handle, rki_api).  Skipped when
# -DRAKALI_ENABLE_NETCDF=OFF because those modules aren't built.
set(NETCDF_REQUIRED_TESTS
    test_forcing
    test_nesting
    test_gauges
    test_rki_api
    test_rki_handle
    test_cross_backend
    test_bathymetry
    test_output
    test_output_unstr
    test_gauges_unstr
    test_restart_unstr
    test_flux_dump
    test_driver_coastal
    test_driver_unstr
    test_rki_api_oneshot
    test_rki_api_twophase
    test_rki_api_two_handles
    test_rki_api_bad_handle
    test_rki_api_string)

# ---------------------------------------------------------------------------
# LFortran 0.64 (alpha) known failures.  Two categories, handled differently:
#
# (A) UNBUILDABLE — do not COMPILE under lfortran, so they must be excluded from
# the build (otherwise `cmake --build` itself fails; they can't be "tests that
# fail"): test_config   - LLVM-backend codegen ICE ("x is nullptr") when a unit
# consumes the .mod of a module that (transitively) calls rki_config's
# read_config_from_string. test_profiler - VOLATILE attribute statement,
# rejected by lfortran 0.64 ("Attribute declaration not supported").
#
# (B) RUNTIME FAILURES — BUILD fine but crash when run.  Deliberately LEFT IN
# the suite as live regression markers: they will start passing once lfortran
# fixes the bug, so re-running the suite after an upgrade tells us for free.
# test_nml_schema, test_config_schema - SIGSEGV on allocate(dst, source=src) of
# a type with a character- pointer component (the nml_string/nml_enum key
# `tgt`).
#
# MREs for every bug live in tmp_local_artifacts/lf_nz_mre/; track upstream at
# lfortran/lfortran.  Extend either list as more surface.
# ---------------------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID STREQUAL "LFortran")
  set(RAKALI_LFORTRAN_RUNTIME_FAILURES test_nml_schema test_config_schema)

  # (A) remove the unbuildable ones from the build
  foreach(_row ${RAKALI_TESTS})
    string(REPLACE "|" ";" _fields "${_row}")
    list(GET _fields 0 _tname)
    list(FIND RAKALI_LFORTRAN_UNBUILDABLE "${_tname}" _skip_idx)
    if(NOT _skip_idx EQUAL -1)
      list(REMOVE_ITEM RAKALI_TESTS "${_row}")
    endif()
  endforeach()
  # (B) runtime failures stay in RAKALI_TESTS on purpose (regression markers).

  message(STATUS "")
  message(
    STATUS "LFortran 0.64 known failures (re-check after any lfortran upgrade;")
  message(STATUS "  MREs in tmp_local_artifacts/lf_nz_mre/):")
  message(STATUS "    excluded from build (do not compile):")
  foreach(_t ${RAKALI_LFORTRAN_UNBUILDABLE})
    message(STATUS "      ${_t}")
  endforeach()
  message(STATUS "    built, expected to FAIL at runtime (kept as markers):")
  foreach(_t ${RAKALI_LFORTRAN_RUNTIME_FAILURES})
    message(STATUS "      ${_t}")
  endforeach()
  message(STATUS "")
endif()

# Valid values for the <regime> field of a RAKALI_TESTS row (see below).
# `ocean+core` expands to BOTH labels (cross-regime tests that stay closed-side
# but also exercise shared infrastructure).
set(RAKALI_TEST_REGIMES coastal ocean core ocean+core)

# Build per-test executables from the shared template.
foreach(row ${RAKALI_TESTS})
  string(REPLACE "|" ";" fields "${row}")
  list(LENGTH fields n_fields)
  if(NOT n_fields EQUAL 5)
    message(
      FATAL_ERROR
        "RAKALI_TESTS row '${row}' has ${n_fields} fields; expected 5: "
        "<test-name>|<module-name>|<collect-fn>|<suite-display>|<regime> "
        "with <regime> one of: ${RAKALI_TEST_REGIMES}")
  endif()
  list(GET fields 0 test_name)
  list(GET fields 1 MODULE_NAME)
  list(GET fields 2 COLLECT_FN)
  list(GET fields 3 SUITE_NAME)
  list(GET fields 4 TEST_REGIME)
  if(NOT TEST_REGIME IN_LIST RAKALI_TEST_REGIMES)
    message(
      FATAL_ERROR "RAKALI_TESTS row '${row}': invalid regime '${TEST_REGIME}' "
                  "(must be one of: ${RAKALI_TEST_REGIMES})")
  endif()
  # `+`-separated regimes become multiple ctest LABELS.
  string(REPLACE "+" ";" TEST_REGIME_LABELS "${TEST_REGIME}")

  if(NOT RAKALI_ENABLE_NETCDF AND ${test_name} IN_LIST NETCDF_REQUIRED_TESTS)
    continue()
  endif()

  set(TESTER_NAME "tester_${test_name}")
  set(driver_path "${CMAKE_CURRENT_BINARY_DIR}/generated/${TESTER_NAME}.F90")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tester_template.F90.in
                 ${driver_path} @ONLY)

  # Source file is named after the MODULE, not the binary -- this lets multiple
  # per-test binaries share one source file when present-table isolation
  # requires running each test in its own process (e.g. OpenMP target offload
  # with module-scope allocatables). For existing rows where binary == module
  # name this is a no-op.
  add_executable(${test_name} ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}.F90
                              ${driver_path})
  # Per-target module dir so multiple binaries sharing one source (rki_api
  # split, zstar_full_tracer_conservation split, etc.) don't race to write the
  # same .mod file under parallel make. Each target gets its own scratch dir;
  # the shared modules/ dir is still on the include path so the test source can
  # use core_rakali / analytical.
  set_target_properties(
    ${test_name} PROPERTIES Fortran_MODULE_DIRECTORY
                            "${CMAKE_CURRENT_BINARY_DIR}/modules/${test_name}")
  target_include_directories(
    ${test_name} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/modules/${test_name}"
                         "${CMAKE_BINARY_DIR}/modules")
  # Core selection: the SHARED core is the default fast link; only tests in the
  # RAKALI_STATIC_CORE_TESTS authority list (device code in their own source)
  # fall back to the static core. When the shared core wasn't built
  # (RAKALI_TESTS_LINK_SHARED=OFF) every test links static.
  set(_core_for_test ${main_lib})
  if(TARGET ${main_lib}_shared AND NOT (${test_name} IN_LIST
                                        RAKALI_STATIC_CORE_TESTS))
    set(_core_for_test ${main_lib}_shared)
  endif()
  target_link_libraries(
    ${test_name} PRIVATE ${_core_for_test} rakali_analytical
                         test-drive::test-drive)
  add_test(NAME rakali_${test_name} COMMAND ${test_name})
  set_tests_properties(rakali_${test_name} PROPERTIES LABELS
                                                      "${TEST_REGIME_LABELS}")
  if(${test_name} IN_LIST NETCDF_TESTS)
    set_tests_properties(rakali_${test_name} PROPERTIES ENVIRONMENT
                                                        "OMP_NUM_THREADS=1")
  endif()
  # SENTINEL: nvfortran 26.5 miscompiles `transfer()` inside `do concurrent`
  # when it lowers the loop for the GPU (-stdpar=gpu; ~95% of elements wrong, a
  # shifted-index bug; the sentinel source documents it in full). The test
  # asserts the CORRECT behaviour, so it fails there today; invert the verdict
  # for NVHPC GPU builds ONLY, so the suite is quiet while the bug exists AND
  # goes red the moment an NVHPC release fixes it.  When that happens: drop this
  # block, drop the `dc-transfer` pre-commit lint, and rki_ocean_chksum can move
  # to `do concurrent`.  Every other configuration -- including an NVHPC CPU
  # build, whose serial `do concurrent` is correct -- runs it as an ordinary
  # test, so a regression there is caught the normal way.
  if(${test_name} STREQUAL "test_dc_transfer_sentinel"
     AND CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC"
     AND RAKALI_ENABLE_GPU)
    set_tests_properties(rakali_${test_name} PROPERTIES WILL_FAIL TRUE)
  endif()
endforeach()
