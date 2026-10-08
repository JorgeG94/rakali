### &sim_nml

Simulation-regime selector.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `sim_type` | `"coastal"` |  | Simulation regime: coastal |

### &grid_nml

Structured-grid geometry.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `grid_type` | `"structured"` |  | Grid discretization: structured (Cartesian) or unstructured (triangular) |
| `nx` | `200` |  | Number of physical cells in x |
| `ny` | `1` |  | Number of physical cells in y |
| `dx` | `0.1000000000E+01` | m | Cell size in x |
| `dy` | `0.1000000000E+01` | m | Cell size in y |
| `nghost` | `3` |  | Ghost cells on each side |

### &time_nml

Time-integration controls.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `t_end` | `0.1000000000E+01` |  | Simulation end time (in time_unit) |
| `cfl` | `0.4500000000E+00` |  | CFL number for adaptive timestep |
| `dt_max` | `0.1000000000E+11` | s | Maximum allowable timestep |
| `dt_fixed` | `0.0000000000E+00` | s | Fixed timestep (0 = adaptive CFL) |
| `cfl_interval` | `1` |  | Recompute CFL timestep every N steps |
| `time_scheme` | `"rk2"` |  | Time integrator (coastal 2D): euler/rk2/rk3/rk43 (explicit) or semi_implicit (unstructured grids only) |
| `time_unit` | `"s"` |  | Unit for the long-time fields: s/min/hr/day/year |

### &si_nml

Semi-implicit free-surface CG-solver controls.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `precond` | `"none"` |  | SI Helmholtz preconditioner: none (plain CG), jacobi (diagonal PCG), or mg (aggregation-multigrid V-cycle PCG) |
| `rtol` | `0.1000000000E-09` |  | SI CG relative residual tolerance (||r||^2 < rtol^2*||rhs||^2) |
| `wetdry` | `.false.` |  | Casulli subgrid wetting/drying (V(eta) Newton diagonal); default off = byte-identical all-wet |
| `wetdry_nsub` | `1` |  | Subgrid sub-bottoms per cell for the wet/dry volume function (1 = flat cell bottom) |
| `wetdry_bsub_fill` | `"flat"` |  | Sub-bottom table fill when nsub>=2: 'flat' (default, bit-identical) or 'linear' (analytic subgrid bed ramp) |
| `newton_max` | `10` |  | Casulli wet/dry Newton outer-iteration cap (only consulted when wetdry=.true.); default 10 = bit-identical |
| `newton_tol` | `0.1000000000E-10` |  | Casulli wet/dry Newton relative residual tolerance (||F||^2 < tol^2*||rhs||^2); default 1e-11 is tight (tuned for smooth cases) — real bathymetry typically wants ~1e-7 |

### &metadata_nml

CF-output global attributes.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `reference_time` | `""` |  | Time origin for CF output (empty = idealized run) |
| `institution` | `""` |  | Institution name for the CF global attribute (empty = omitted) |

### &logging_nml

Logger verbosity + status cadence.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `log_level` | `"info"` |  | Log verbosity |
| `status_interval` | `0.0000000000E+00` |  | Status-print cadence (in time_unit; 0 = every 100 steps) |

### &vcoord_nml

Vertical-coordinate + ALE-remap controls.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `vcoord_type` | `"sigma"` |  | Vertical coordinate type |
| `remap_method` | `"ppm"` |  | Vertical remapping method |
| `vcoord_depth_transition` | `0.1000000000E+03` | m | z-sigma hybrid: depth where blending starts |
| `vcoord_blend_width` | `0.5000000000E+02` | m | z-sigma hybrid: blending-zone width |
| `zstar_h_surf_target` | `0.0000000000E+00` | m | z*-full: target surface-layer thickness (0 = auto) |
| `zstar_stretching` | `"log"` |  | z*-full: surface-concentration stretching |
| `zstar_h_min` | `0.1000000000E-03` | m | z*-full: vanishing-layer floor |
| `zstar_n_surf` | `0` |  | z*-full: number of fine near-surface layers (0 = auto) |
| `rho_ref_pressure` | `0.2000000000E+08` | Pa | rho-coord: reference pressure for potential density |
| `rho_target_light` | `0.1020000000E+04` | kg/m^3 | rho-coord: lightest (surface) target density |
| `rho_target_dense` | `0.1030000000E+04` | kg/m^3 | rho-coord: densest (bed) target density |

### &physics_nml

Barotropic physics: bottom drag, wind stress, Coriolis.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `manning_n` | `0.0000000000E+00` |  | Manning roughness coefficient |
| `dry_threshold` | `0.1000000000E-02` | m | Depth below which a cell is dry |
| `wind_stress_x` | `0.0000000000E+00` | Pa | Surface wind stress in x |
| `wind_stress_y` | `0.0000000000E+00` | Pa | Surface wind stress in y |
| `coriolis_f` | `0.0000000000E+00` | 1/s | Coriolis parameter f |
| `use_hllc` | `.false.` |  | Use HLLC contact-wave Riemann solver instead of HLL |
| `use_linear_bottom_drag` | `.false.` |  | Linear (Rayleigh) bottom drag instead of quadratic Manning |
| `bottom_drag_r` | `0.0000000000E+00` | 1/s | Linear bottom-drag rate |
| `wind_surface_layer` | `.false.` |  | Apply wind stress to the top layer only (multilayer) |

### &output_nml

File output, forcing/restart/gauge files.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `output_to_file` | `.false.` |  | Enable file output (NetCDF snapshots) |
| `output_dir` | `"./output"` |  | Directory for output files |
| `output_interval` | `0.0000000000E+00` | s | Time between output snapshots |
| `output_format` | `"netcdf"` |  | Output format |
| `restart_interval` | `0.0000000000E+00` | s | Time between restart writes (0 = none) |
| `async_output` | `.false.` |  | Async ring-buffered output (CUDA Fortran) |
| `async_buffer_slots` | `4` |  | Async ring-buffer slot count |
| `compress_output` | `.false.` |  | Deflate-compress NetCDF output |
| `compress_level` | `1` |  | Deflate level (1=fast, 9=max) |
| `bathymetry_file` | `""` |  | NetCDF bathymetry file (empty = flat) |
| `restart_file` | `""` |  | Restart file for warm start (empty = cold) |
| `forcing_file` | `""` |  | NetCDF atmospheric forcing file (empty = uniform) |
| `forcing_x_origin` | `0.0000000000E+00` |  | Model x-origin in forcing-file coordinates |
| `forcing_y_origin` | `0.0000000000E+00` |  | Model y-origin in forcing-file coordinates |
| `gauge_file` | `""` |  | Gauge-locations text file (empty = none) |
| `gauge_output_interval` | `0.0000000000E+00` | s | Gauge output interval (0 = every step) |
| `out_h_layer` | `.true.` |  | Output layer thickness |
| `out_velocity` | `.true.` |  | Output layer velocity (u_ml, v_ml) |
| `out_salinity` | `.true.` |  | Output layer salinity |
| `out_temperature` | `.true.` |  | Output layer temperature |
| `out_density` | `.false.` |  | Output layer density |
| `out_w` | `.false.` |  | Output vertical velocity at interfaces |

### &boundary_nml

Boundary conditions, tidal forcing, inflow/discharge/sponge/nesting.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `bc_west` | `"wall"` |  | West boundary type |
| `bc_east` | `"wall"` |  | East boundary type |
| `bc_south` | `"wall"` |  | South boundary type |
| `bc_north` | `"wall"` |  | North boundary type |
| `tidal_amplitude` | `0.0000000000E+00` | m | Tidal elevation amplitude (single-constituent) |
| `tidal_period` | `0.4471200000E+05` | s | Tidal period (single-constituent, default M2) |
| `n_tidal_constituents` | `0` |  | Active tidal constituents (0 = legacy single) |
| `tidal_amp` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | m | Constituent amplitudes |
| `tidal_phase` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | rad | Constituent phases |
| `tidal_omega` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | rad/s | Constituent angular frequencies |
| `inflow_velocity` | `0.0000000000E+00` | m/s | Prescribed normal velocity for inflow BC |
| `inflow_salinity` | `-0.1000000000E+01` | PSU | Inflow salinity (<0 = zero-gradient) |
| `inflow_temperature` | `-0.9990000000E+03` | degC | Inflow temperature (<0 = zero-gradient) |
| `discharge_flux` | `0.0000000000E+00` | m^2/s | Prescribed volume flux for discharge BC |
| `clamped_elevation` | `0.0000000000E+00` | m | Prescribed surface elevation for clamped BC |
| `sponge_width` | `0` |  | Sponge layer width in cells |
| `sponge_strength` | `0.0000000000E+00` | 1/s | Sponge relaxation rate |
| `sponge_h_ref` | `0.0000000000E+00` | m | Sponge reference depth |
| `nesting_file` | `""` |  | Parent-domain NetCDF output (empty = no nesting) |
| `nesting_x_origin` | `0.0000000000E+00` | m | Child x-origin in parent coordinates |
| `nesting_y_origin` | `0.0000000000E+00` | m | Child y-origin in parent coordinates |

### &tracer_nml

Salinity + temperature IC/EOS/bounds, sediment transport.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `initial_salinity` | `0.3500000000E+02` | PSU | Initial salinity (uniform IC) |
| `S_ref` | `0.0000000000E+00` | PSU | EOS reference salinity |
| `beta_S` | `0.7800000000E+00` | kg/m^3/PSU | Haline contraction coefficient |
| `S_min` | `0.0000000000E+00` | PSU | Lower physical bound for salinity |
| `S_max` | `0.4000000000E+02` | PSU | Upper physical bound for salinity |
| `kappa_S_bg` | `0.1000000000E-04` | m^2/s | Background vertical salinity diffusivity |
| `S_init_surface` | `0.0000000000E+00` | PSU | Initial surface salinity (stratified IC) |
| `S_init_bottom` | `0.0000000000E+00` | PSU | Initial bed salinity (stratified IC) |
| `initial_temperature` | `0.1500000000E+02` | degC | Initial temperature (uniform IC) |
| `T_ref` | `0.1500000000E+02` | degC | EOS reference temperature |
| `alpha_T` | `0.1700000000E+00` | kg/m^3/degC | Thermal expansion coefficient |
| `T_min` | `-0.2000000000E+01` | degC | Lower physical bound for temperature |
| `T_max` | `0.4000000000E+02` | degC | Upper physical bound for temperature |
| `kappa_T_bg` | `0.1000000000E-04` | m^2/s | Background vertical temperature diffusivity |
| `eos_form` | `"linear"` |  | Coastal EOS: linear (legacy, default), wright or roquet_spv (shared nonlinear) |
| `T_init_surface` | `0.0000000000E+00` | degC | Initial surface temperature (stratified IC) |
| `T_init_bottom` | `0.0000000000E+00` | degC | Initial bed temperature (stratified IC) |
| `use_ppm_tracer` | `.true.` |  | High-order tracer advection (PPM/MUSCL) vs upwind |
| `tracer_recon` | `"ppm"` |  | Tracer face reconstruction (structured coastal): ppm|weno5|weno7|weno9 |
| `use_sediment` | `.false.` |  | Master switch for suspended-sediment transport |
| `sediment_w_s` | `0.0000000000E+00` | m/s | Settling velocity (downward; 0 = none) |
| `sediment_initial_conc` | `0.0000000000E+00` | kg/m^3 | Uniform initial suspended concentration |
| `sediment_bed_init` | `0.0000000000E+00` | kg/m^2 | Uniform initial bed sediment mass |
| `sediment_tau_crit_erosion` | `0.0000000000E+00` | Pa | Critical erosion bed shear stress (0 = off) |
| `sediment_erosion_rate` | `0.0000000000E+00` | kg/m^2/s | Ariathurai erosion-rate coefficient E0 |
| `sediment_porosity` | `0.0000000000E+00` |  | Bed porosity (0-1); erosion scaled by (1-porosity) |

### &initial_condition_nml

Coastal initial-condition selector + parameters.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `init_condition` | `"dam_break"` |  | Initial condition |
| `eta_initial` | `0.0000000000E+00` | m | Initial free-surface elevation (still_water IC) |
| `h0` | `0.1000000000E+01` | m | Background depth (gaussian_hump IC) |
| `hump_amplitude` | `0.5000000000E-01` | m | Perturbation amplitude (gaussian_hump IC) |
| `hump_sigma` | `0.5000000000E+01` | m | Gaussian width (gaussian_hump IC) |
| `hump_x_offset` | `0.0000000000E+00` | m | Hump centre x-offset from domain centre |
| `hump_y_offset` | `0.0000000000E+00` | m | Hump centre y-offset from domain centre |
| `init_u_shear` | `0.0000000000E+00` | m/s | Initial vertical shear (top +val, bottom -val) |

### &nonhydrostatic_nml

Non-hydrostatic / multilayer: CG-Poisson, PP81/KPP/k-eps vmix, Smagorinsky, BPG, mode split.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `use_nonhydrostatic` | `.false.` |  | Enable multi-layer NH pressure correction |
| `nz_layers` | `2` |  | Number of sigma layers |
| `nh_tolerance` | `0.1000000000E-05` |  | CG convergence tolerance for NH pressure |
| `nh_max_iterations` | `100` |  | Maximum CG iterations for NH pressure |
| `nh_preconditioner` | `.false.` |  | Vertical-line preconditioner for CG solver |
| `nh_solve_interval` | `1` |  | Solve NH pressure every N steps |
| `split_ratio` | `1` |  | Mode-split ratio: 0=auto, 1=unsplit, N>1=fixed |
| `split_averaging` | `.true.` |  | ROMS-style averaging over barotropic substeps |
| `use_multilayer` | `.false.` |  | Enable coupled hydrostatic vertical layers |
| `rho_0` | `0.1000000000E+04` | kg/m^3 | Reference density for EOS |
| `vmix_nu_max` | `0.1000000000E-01` | m^2/s | PP81 max turbulent viscosity |
| `vmix_nu_bg` | `0.1000000000E-03` | m^2/s | Background vertical viscosity |
| `vmix_ri_alpha` | `0.5000000000E+01` |  | PP81 Richardson-number scaling coefficient |
| `vmix_use_kpp` | `.false.` |  | Enable KPP surface-boundary-layer overlay on PP81 |
| `vmix_use_keps` | `.false.` |  | Enable two-equation k-epsilon vertical mixing |
| `keps_si_picard` | `1` |  | Picard passes of the k<->eps coupling per semi-implicit step (SI+keps only) |
| `kpp_ri_crit` | `0.3000000000E+00` |  | Critical bulk Richardson number for KPP BL-depth |
| `kpp_cs_nonlocal` | `0.6300000000E+01` |  | KPP non-local (counter-gradient) transport coefficient |
| `kpp_c_vt2` | `0.1800000000E+01` |  | KPP V_t^2 unresolved-turbulence coefficient (0 = off) |
| `hdiff_kappa` | `0.0000000000E+00` | m^2/s | Horizontal tracer diffusion coefficient |
| `smagorinsky_cs` | `0.0000000000E+00` |  | Smagorinsky horizontal-viscosity coefficient |
| `biharmonic_nu4` | `0.0000000000E+00` | m^4/s | Biharmonic lateral viscosity (0 = off) |
| `evap_mass_sink` | `.false.` |  | Treat net evaporation as a real surface mass sink (default off = virtual salt flux) |
| `evap_latent_cooling` | `.false.` |  | Evaporative latent-heat surface cooling (requires evap_mass_sink; default off) |
| `ml_min_depth` | `0.5000000000E+00` | m | Minimum total depth for vertical-layer operations |
| `bpg_method` | `"jacobian"` |  | BPG computation method |
| `keps_stability` | `"constant"` |  | k-epsilon stability functions |
