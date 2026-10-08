# Core infrastructure sources
set(CORE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rki_constants.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rki_console_stats.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rki_nml_schema.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rki_array_utils.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/ALE/rki_vcoord.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rki_timestepper.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rki_si_common.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rki_config.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rki_config_schema.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/structured/rki_grid.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/structured/rki_tracer.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/structured/rki_surface_forcing.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/structured/rki_barotropic_state.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/framework/rki_scratch_3d.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/framework/rki_safe_math.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/framework/rki_efp.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/framework/rki_mem_report.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/structured/rki_multilayer_state.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/structured/rki_nonhydrostatic_state.F90
    # Ocean-regime god state + per-kernel state types.  All Phase 0d empty
    # shells; declared here (not under solver/kernel sources) because rki_state
    # composes ocean_state_t.
    ${CMAKE_CURRENT_SOURCE_DIR}/equation_of_state/structured/rki_eos.F90
    # Ocean-side EOS shim over the C-grid multilayer state — keeps rki_eos
    # itself free of ocean-state dependencies (the shared EOS is consumed by the
    # coastal path too); depends on rki_eos + rki_multilayer_cgrid_state (both
    # declared above). Sea-ice PR 2: pure device-callable enthalpy library
    # (T<->E inversions); leaf module. Precedes rki_ice_optics/_mass/_column (PR
    # 3a), which all depend on it. Sea-ice PR 3a: shortwave optics (CSIM4
    # branch) + column mass bookkeeping (bottom-freeze, melt peel, rebalance) —
    # both leaf modules `rki_ice_column` depends on. Sea-ice PR 3a:
    # single-column Winton NkIce-layer thermodynamics (SEB + conduction solve,
    # per-column orchestrator, do-concurrent driver). Exercised only by the new
    # test suite — nothing in the driver calls it yet. Sea-ice slot scaffold
    # (SIS2 port, PR 0) — composed by rki_ocean_state, so it must precede it in
    # the list. PR 3a: now depends on rki_ice_column for
    # ICE_BULK_SALINITY/ICE_NK_MAX (the Winton column prognostics' init values).
    # Sea-ice PR 24: analytic initial-condition path (uniform/latitudes
    # concentration seeding). Host-side setup code, no device kernel; depends on
    # rki_ice_state, rki_ice_column (ICE_RHO_ICE/ICE_RHO_SNOW/
    # ICE_BULK_SALINITY), rki_ice_enthalpy (ice_enth_from_ts), rki_grid,
    # rki_multilayer_cgrid_state, rki_ocean_metrics (all declared above).
    # Sea-ice PR 4a: multi-category ITD restore (adjust_ice_categories); depends
    # on rki_ice_state, rki_ice_column (ICE_RHO_ICE et al, declared above),
    # rki_constants, rki_grid, rki_multilayer_cgrid_state. Sea-ice PR 4b:
    # category ice/snow transport + compress_ice; depends on rki_ice_state,
    # rki_ice_itd (ice_adjust_categories, Phase 5), rki_continuity (the 5
    # promoted PPM helpers, declared earlier in this file at the
    # ocean/kernels/structured/continuity_ppm entry),
    # rki_multilayer_cgrid_state, rki_ocean_metrics. Sea-ice PR 5: C-grid EVP
    # ice dynamics; depends on rki_ice_state (ice_cell_concentration_impl),
    # rki_ice_column (ICE_RHO_ICE), rki_ocean_metrics,
    # rki_multilayer_cgrid_state, rki_ocean_periodic (the face/centre ghost-wrap
    # helpers, declared earlier in this file at the ocean/boundary entry).
    # Sea-ice PR 1: ocean-side frazil accumulator (freezing-point clamp), called
    # from the driver after each ocean outer step. Sea-ice PR 3c: coupleable
    # atmospheric-forcing seam filler (v1 scalar restoring); depends only on
    # rki_ice_state (already declared above). Sea-ice PR 3c: ocean -> ice basal
    # heat flux + the sst/s_surf/tfw sample seam; depends on rki_ice_state +
    # rki_ocean_surface_flux (SEAWATER_CP, declared earlier in this file) +
    # rki_ice_column (ICE_RHO_ICE, for the ice-presence gate — declared above at
    # rki_ice_column.F90). Sea-ice PR 3b: frazil-bank spend (new-ice formation)
    # + the ice->ocean brine-rejection coupler (Q_salt refresh); both called
    # from the driver at thermo cadence, after rki_ice_frazil in the loop.
    # Sea-ice PR 3c: driver-facing column step + ocean-coupling-diag reduction;
    # depends on rki_ice_column (Winton column) + rki_ice_state. Sea-ice PR 26:
    # snowfall's ice-free-share ocean contributor (heat + virtual-salt); depends
    # on rki_ice_state + rki_ice_enthalpy (declared above). Runs AFTER
    # rki_ice_thermo_driver, BEFORE rki_ice_ocean_coupler in the driver's
    # thermo-cadence block (mandated order).
    ${CMAKE_CURRENT_SOURCE_DIR}/core/structured/rki_state.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/structured/rki_boundary_types.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rki_init.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/rki_coastal_console_stats.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/rki_profiler.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/unstructured/rki_mesh_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/unstructured/rki_mesh_gen_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/unstructured/rki_boundary_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/unstructured/rki_tracer_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/unstructured/rki_barotropic_state_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/unstructured/rki_multilayer_state_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/unstructured/rki_surface_forcing_unstr.F90
    # Aggregation-MG preconditioner for the unstructured SI eta-solve; used by
    # rki_si_unstr_state, so it must precede it in the list.
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/barotropic/rki_si_unstr_mg.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/unstructured/rki_si_unstr_state.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/unstructured/rki_state_unstr.F90)

# Atmospheric forcing + nesting both load NetCDF inputs.  Skip them when the I/O
# subsystem is disabled so the kernels/solver/benchmarks still build.
if(RAKALI_ENABLE_NETCDF)
  list(APPEND CORE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/core/rki_forcing.F90
       ${CMAKE_CURRENT_SOURCE_DIR}/core/rki_nesting.F90)
endif()

# Solver drivers. * core/coastal/solver/rki_solver{,_unstr}.F90 — HLLC/KNP
# SSP-RK2. The ocean path (sim_type='ocean') dispatches through
# `driver_run_ocean` (src/driver/rki_driver.F90) directly into the split-RK2
# dyn-core (core/ocean/dynamics/split_rk2/) — there is no separate ocean
# solver-driver module (the former Phase-0 stub `rki_ocean_solver.F90` was
# removed as dead code; PR-8).
set(SOLVER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/solver/rki_solver.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/solver/rki_solver_unstr.F90)

# High-level driver — owns the full compute-rank lifecycle.  Lives above the
# solver, below `app/main.F90`.  Test-callable so unit tests can run a full
# lifecycle and catch integration-level bugs. Pulls in
# output/gauges/forcing/nesting, so it requires NetCDF.
set(DRIVER_SOURCES)
if(RAKALI_ENABLE_NETCDF)
  list(APPEND DRIVER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/driver/rki_driver.F90)
endif()

set(KERNEL_SOURCES
    # ----- Coastal dynamics core (HLL/HLLC + KNP + NH) -----
    # Structured barotropic — pure 2D HLL coastal/flooding kernels
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/structured/barotropic/rki_kernel_flux.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/structured/barotropic/rki_kernel_update.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/structured/barotropic/rki_kernel_timestep.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/structured/barotropic/rki_kernel_friction.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/structured/barotropic/rki_kernel_sponge.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/structured/barotropic/rki_boundary.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/structured/barotropic/rki_kernel_fused.F90
    # Structured multilayer — coastal-coupled (HLLC-flavoured mode-split)
    ${CMAKE_CURRENT_SOURCE_DIR}/core/multilayer/structured/rki_ml_common.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/structured/multilayer/rki_ml_dynamics.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/tracer/structured/rki_ml_tracers.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/tracer/structured/rki_recon_weno.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/equation_of_state/structured/rki_ml_eos.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/pressure_force/structured/rki_ml_pressure.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/structured/rki_ml_vertical_viscosity.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/structured/rki_ml_kpp.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/rki_keps_stability.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/structured/rki_ml_keps.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/structured/rki_ml_friction.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/lateral/structured/rki_ml_horizontal_viscosity.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/multilayer/structured/rki_ml_truncate.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/structured/multilayer/rki_kernel_fused_ml.F90
    # Structured non-hydrostatic option (CG Poisson)
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/structured/nonhydrostatic/rki_kernel_pressure.F90
    # Unstructured barotropic — KNP central-upwind kernels
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/barotropic/rki_kernel_extrapolate_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/barotropic/rki_kernel_flux_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/barotropic/rki_kernel_update_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/barotropic/rki_kernel_protect_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/barotropic/rki_kernel_timestep_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/barotropic/rki_kernel_friction_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/barotropic/rki_kernel_sponge_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/barotropic/rki_kernel_fused_unstr.F90
    # Unstructured semi-implicit free-surface option (TPFA Helmholtz CG
    # eta-solve)
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/barotropic/rki_si_unstr_helmholtz.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/barotropic/rki_si_unstr_wetdry.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/barotropic/rki_si_unstr_step.F90
    # Unstructured multilayer — coastal-coupled
    ${CMAKE_CURRENT_SOURCE_DIR}/core/multilayer/unstructured/rki_ml_common_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/multilayer/unstructured/rki_ml_truncate_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/equation_of_state/unstructured/rki_ml_eos_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/multilayer/rki_ml_dynamics_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/tracer/unstructured/rki_ml_tracers_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/pressure_force/unstructured/rki_ml_pressure_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/unstructured/rki_ml_vertical_viscosity_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/unstructured/rki_ml_kpp_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/unstructured/rki_ml_keps_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/vertical/unstructured/rki_ml_friction_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameterizations/lateral/unstructured/rki_ml_horizontal_viscosity_unstr.F90
    # Unstructured semi-implicit 3D multilayer (column solve + dsig tracers)
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/multilayer/rki_si_unstr_column.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/multilayer/rki_si_unstr_closures.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/core/coastal/kernels/unstructured/multilayer/rki_si_unstr_step_ml.F90
    # Shared remap kernel (used by both coastal and ocean) + ocean- side
    # orchestrator that wires it across the tracer registry — both live under
    # ALE/ since they're vcoord-machinery, not split-RK2- specific.
    ${CMAKE_CURRENT_SOURCE_DIR}/ALE/rki_remap_column.F90)
if(RAKALI_PARALLEL_BACKEND STREQUAL "openmp")
  set(_kernel_directive_label "OpenMP target")
else()
  set(_kernel_directive_label "OpenACC")
endif()
if(RAKALI_ENABLE_GPU)
  message(
    STATUS
      "Kernel backend: GPU (do concurrent + ${_kernel_directive_label} offload)"
  )
else()
  message(
    STATUS
      "Kernel backend: CPU (do concurrent serial, ${_kernel_directive_label} directives ignored)"
  )
endif()

# I/O sources (NetCDF)
set(IO_SOURCES)
if(RAKALI_ENABLE_NETCDF)
  list(
    APPEND
    IO_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/io/rki_io_netcdf.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/io/rki_banner.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/io/rki_bathymetry.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/io/structured/rki_output.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/io/structured/rki_flux_dump.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/io/unstructured/rki_output_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/io/unstructured/rki_restart_unstr.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/io/structured/rki_gauges.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/io/unstructured/rki_gauges_unstr.F90)
endif()

# C-compatible API (handle-based, multi-instance) for Python/Julia/C callers.
# rki_api / rki_api_unstr / rki_handle pull in the NetCDF-backed
# output/gauges/restart modules, so they only build with NetCDF enabled.
set(API_SOURCES)
if(RAKALI_ENABLE_NETCDF)
  list(APPEND API_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/api/rki_handle.F90
       ${CMAKE_CURRENT_SOURCE_DIR}/api/rki_api.F90
       ${CMAKE_CURRENT_SOURCE_DIR}/api/rki_api_unstr.F90)
endif()

target_sources(
  ${core_objs} PRIVATE ${CORE_SOURCES} ${IO_SOURCES} ${SOLVER_SOURCES}
                       ${KERNEL_SOURCES} ${API_SOURCES} ${DRIVER_SOURCES})
