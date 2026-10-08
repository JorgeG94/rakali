!! Runtime configuration via Fortran namelists
module rki_config
   !! Reads simulation parameters from a namelist input file
#ifdef LFORTRAN_PASSING
   use rki_constants, only: wp, DEG2RAD, &
                            nz_stack_required, nz_stack_is_sufficient
#else
   use rki_constants, only: NZ_STACK_MAX, wp, DEG2RAD, &
                            nz_stack_required, nz_stack_is_sufficient
#endif
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use rki_nml_schema, only: nml_schema_t, nml_group_t, &
                             nml_real, nml_int, nml_logical, nml_string, &
                             nml_enum, nml_real_array, nml_real_array_key_t
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rki_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: config_t
   public :: read_config
   public :: read_config_from_string
   public :: validate_config
   public :: semi_implicit_grid_conflict
   public :: build_rakali_schema
   public :: MAX_TIDAL_CONSTITUENTS

   integer, parameter :: MAX_TIDAL_CONSTITUENTS = 10
      !! Maximum number of tidal constituents

   type :: config_t
      !! Container for all runtime simulation parameters

      ! Simulation regime — selects the dynamical core path:
      !   "coastal" → HLL/HLLC (structured) or KNP central-upwind (unstructured)
      !               + nonhydrostatic option (structured only). Shock-capturing,
      !               wet/dry friendly, default for the existing test suite.
      character(len=16) :: sim_type = "coastal"
         !! Simulation regime: "coastal" is the only supported value

      ! Grid parameters
      character(len=16) :: grid_type = "structured"
         !! Grid discretization: "structured" (Cartesian) or "unstructured"
         !! (triangular).  "unstructured" routes the coastal regime through the
         !! native mesh generator + unstructured solver (Python-free).
      integer :: nx = 200
         !! Number of physical cells in x-direction
      integer :: ny = 1
         !! Number of physical cells in y-direction
      real(wp) :: dx = 1.0_wp
         !! Cell size in x-direction (m)
      real(wp) :: dy = 1.0_wp
         !! Cell size in y-direction (m)
      integer :: nghost = 3
         !! Number of ghost cells on each side.  Default 3 (not 2): a
         !! 5-point PPM reconstruction stencil means a rank-seam face's
         !! ghost-side donor needs two neighbours beyond itself —
         !! full-order reconstruction AT the seam requires three ghost
         !! columns.  At nghost=2 the PPM
         !! local-array-edge fallback fires at the seam and degrades the
         !! face to first order, producing a seam-inconsistent tracer mass
         !! flux: it conserves a uniform tracer (constant field cancels)
         !! but LEAKS a structured one (advecting T with rank-A-out /=
         !! rank-B-in on the shared seam face drifts hT), invisible to the
         !! mass/salt closure gates.  Measured on seamount_bench_full
         !! (np2 x-split): heat closure grew to +5e-11 by day 2 at nghost=2,
         !! stays at round-off (~5e-14, no growth) at nghost=3.  Single-rank
         !! runs are unaffected by the wider halo (same interior); the extra
         !! ghost column only costs a little memory.

      ! Time parameters
      real(wp) :: t_end = 1.0_wp
         !! Simulation end time, interpreted in `time_unit` (default s).
      real(wp) :: cfl = 0.45_wp
         !! CFL number for adaptive timestep
      real(wp) :: dt_max = 1.0e10_wp
         !! Maximum allowable timestep (s) — always in seconds, NOT
         !! affected by `time_unit`.  Timesteps are naturally short
         !! and CFL-bounded; writing them in `"day"` would force
         !! awkward fractional values.
      real(wp) :: dt_fixed = 0.0_wp
         !! Fixed timestep (s), 0 = adaptive CFL.  Same convention
         !! as `dt_max` — always seconds.
      integer :: cfl_interval = 1
         !! Recompute CFL timestep every N steps (1 = every step)
      character(len=16) :: time_scheme = "rk2"
         !! Explicit time integrator (coastal 2D barotropic path):
         !! "euler" (1st), "rk2" (default, 2nd), "rk3" (3rd), "rk43"
         !! (3rd, SSP coeff 2 -> 2x dt).  All SSP for wet/dry positivity.
      logical :: semi_implicit = .false.
         !! Derived (not a namelist key): .true. when time_scheme=="semi_implicit".
         !! Set in `apply_time_unit_cascade` (runs at the end of the config read,
         !! before `validate_config`) so downstream code branches on the logical,
         !! not the string.
      character(len=8) :: time_unit = "s"
         !! Unit applied to the long-time fields: `t_end` and
         !! `status_interval` (logging cadence).  `dt_fixed` and `dt_max` stay in
         !! seconds — see their docstrings.  Read at namelist parse,
         !! then multiplied through so the rest of the code keeps
         !! using seconds internally.  Valid values: `"s"`, `"min"`,
         !! `"hr"`, `"day"`, `"year"`.  Year is 365.25 days (Julian,
         !! matches MOM6 / UDUNITS convention).  Default `"s"`
         !! preserves backward compatibility — every existing nml is
         !! bit-identical.

      ! Semi-implicit free-surface CG solver (`&si_nml`).  Only consulted
      ! when `time_scheme="semi_implicit"`; ignored otherwise.
      character(len=8) :: si_precond = "none"
         !! Preconditioner for the SI Helmholtz eta-solve: `"none"`
         !! (unpreconditioned CG, default), `"jacobi"` (diagonal PCG,
         !! `z = eta_diag*r`), or `"mg"` (aggregation-multigrid V-cycle
         !! PCG).
         !! Fail-loud allowlist.  Default `"none"` is byte-identical to the
         !! pre-knob plain-CG path.
      real(wp) :: si_rtol = 1.0e-10_wp
         !! CG RELATIVE residual tolerance: converge when
         !! `||r||^2 < si_rtol^2 * ||rhs||^2`.  Default `1e-10` reproduces
         !! the historical hardcoded `SI_TOL_REL` bit-for-bit; relaxable
         !! (1e-8/1e-6) to cut iterations — `si%mass_res` is the audit
         !! instrument for the resulting mass consistency.
      logical :: si_wetdry = .false.
         !! Casulli (2009) subgrid wetting/drying for the SI free surface:
         !! replaces the linear `1.0` cell-mass term by the exact wetted-
         !! volume `V(eta)` (per-cell wet-fraction Newton diagonal) so
         !! intertidal flats dry at minute-scale dt.  Default `.false.` is
         !! byte-identical (all-wet configs unchanged).  Structured coastal
         !! SI: 2D barotropic + 3D sigma multilayer.  Unstructured coastal
         !! SI: 2D barotropic + 3D sigma multilayer (the triangular wet/dry
         !! Newton).  See `si_wetdry_bsub_fill` for the subgrid bed option.
      integer :: si_wetdry_nsub = 1
         !! Subgrid sub-bottoms per cell for the wet/dry volume function
         !! (>=1).  `1` (default) = v1 flat cell bottoms (V piecewise-linear,
         !! one breakpoint per cell).  `>=2` = the subgrid variant (a per-cell
         !! depth distribution; see `si_wetdry_bsub_fill`).  Only consulted
         !! when `si_wetdry=.true.`.
      character(len=32) :: si_wetdry_bsub_fill = "flat"
         !! Sub-bottom table fill mode when `si_wetdry_nsub>=2`:
         !!   `"flat"` (default) — all sub-bottoms = the plain bed `b`, so the
         !!     volume reduces to the flat-bottom case ⇒ BIT-IDENTICAL to
         !!     `nsub=1` (and to today's behaviour).
         !!   `"linear"` — the analytic subgrid variant: `n_sub` equal-area
         !!     sub-bottoms evenly spaced across a per-cell relief band
         !!     `[b - half, b + half]` where `half = 0.5*max` neighbour bed
         !!     difference (a linear sub-cell bed ramp; Casulli 2009 subgrid
         !!     bathymetry).  The mean sub-bottom is exactly `b`, so the
         !!     still-water volume is unchanged; only the partial-flooding
         !!     regime near `eta≈b` gains breakpoints ⇒ a coarse tidal flat
         !!     drains/floods SMOOTHLY (wet fraction ramps 1/n_sub per
         !!     breakpoint instead of a binary 0→1 cell jump).  On flat bed
         !!     (`half=0`) ⇒ identical to `"flat"`.  Only consulted when
         !!     `si_wetdry=.true.` and `si_wetdry_nsub>=2`.
      integer :: si_newton_max = 10
         !! Casulli wet/dry Newton outer-iteration cap (`&si_nml newton_max`).
         !! Default `10` matches the historical hardcoded value (bit-identical).
         !! Real intertidal bathymetry with a shuffling active set may need
         !! more.  Only consulted when `si_wetdry=.true.`.
      real(wp) :: si_newton_tol = 1.0e-11_wp
         !! Casulli wet/dry Newton RELATIVE residual tolerance
         !! (`&si_nml newton_tol`): converged when
         !! `||F||^2 < newton_tol^2 * ||rhs||^2` AND the active set is frozen.
         !! Default `1e-11` matches the historical hardcoded value
         !! (bit-identical) — tight, tuned for smooth validation cases; real
         !! bathymetry typically converges to ~1e-7 and wants that here.
         !! Only consulted when `si_wetdry=.true.`.
      ! Physics parameters
      real(wp) :: manning_n = 0.0_wp
         !! Manning roughness coefficient
      logical :: use_linear_bottom_drag = .false.
         !! Use linear (Rayleigh) bottom drag instead of quadratic Manning.
      real(wp) :: bottom_drag_r = 0.0_wp
         !! Linear bottom-drag rate (1/s); used when use_linear_bottom_drag.
      logical :: wind_surface_layer = .false.
         !! Apply surface wind stress to the top layer only (surface-
         !! concentrated, generates vertical shear) instead of the
         !! depth-uniform barotropic mode.  Multilayer (nz>1) only; off by
         !! default (= bit-identical depth-uniform forcing).
      real(wp) :: dry_threshold = 0.001_wp
         !! Depth below which a cell is considered dry (m)
      real(wp) :: wind_stress_x = 0.0_wp
         !! Surface wind stress in x-direction (Pa)
      real(wp) :: wind_stress_y = 0.0_wp
         !! Surface wind stress in y-direction (Pa)
      real(wp) :: coriolis_f = 0.0_wp
         !! Coriolis parameter f (1/s), typically 2*Omega*sin(lat)
      logical :: use_hllc = .false.
         !! Use the HLLC approximate Riemann solver for the structured
         !! barotropic flux instead of HLL.  HLLC adds a contact wave so
         !! transverse momentum (`hv` across an x-face, `hu` across a
         !! y-face) is upwinded across the shear interface; reduces
         !! numerical diffusion for shear-aligned flow at the cost of
         !! a few extra flops per face.  Default `.false.` (HLL).

      ! Output parameters
      logical :: output_to_file = .false.
         !! Enable file output (NetCDF snapshots)
      character(len=256) :: output_dir = "./output"
         !! Directory for output files
      real(wp) :: output_interval = 0.0_wp
         !! Time between output snapshots (s)
      character(len=16) :: output_format = "netcdf"
         !! Output format: "ascii" or "netcdf"
      real(wp) :: restart_interval = 0.0_wp
         !! Time between restart file writes (s), 0 = no restarts
      logical :: async_output = .false.
         !! Use async ring-buffered output (CUDA Fortran only, ignored otherwise)
      integer :: async_buffer_slots = 4
         !! Number of ring buffer slots for async output (default: 4)
      logical :: compress_output = .false.
         !! Enable deflate compression for NetCDF output (NetCDF4/HDF5)
      integer :: compress_level = 1
         !! Deflate compression level (1=fast, 9=max). 1 is usually optimal.

      ! Output variable selection (multilayer)
      logical :: out_h_layer = .true.
         !! Output layer thickness
      logical :: out_velocity = .true.
         !! Output layer velocity (u_ml, v_ml)
      logical :: out_salinity = .true.
         !! Output layer salinity
      logical :: out_temperature = .true.
         !! Output layer temperature
      logical :: out_density = .false.
         !! Output layer density
      logical :: out_w = .false.
         !! Output vertical velocity at interfaces

      ! I/O file parameters
      character(len=256) :: bathymetry_file = ""
         !! Path to NetCDF bathymetry file (empty = flat bottom)
      character(len=256) :: restart_file = ""
         !! Path to restart file for warm start (empty = cold start)
      character(len=256) :: forcing_file = ""
         !! Path to NetCDF atmospheric forcing file (empty = uniform config values)
      real(wp) :: forcing_x_origin = 0.0_wp
         !! Model grid x-origin in forcing file coordinates
      real(wp) :: forcing_y_origin = 0.0_wp
         !! Model grid y-origin in forcing file coordinates
      character(len=256) :: gauge_file = ""
         !! Path to gauge locations text file (empty = no gauge output)
      real(wp) :: gauge_output_interval = 0.0_wp
         !! Gauge output interval in seconds (0 = every timestep)

      ! Nesting parameters
      character(len=256) :: nesting_file = ""
         !! Path to parent domain NetCDF output (empty = no nesting)
      real(wp) :: nesting_x_origin = 0.0_wp
         !! Child grid x-origin in parent coordinates (m)
      real(wp) :: nesting_y_origin = 0.0_wp
         !! Child grid y-origin in parent coordinates (m)

      ! Boundary parameters
      character(len=16) :: bc_west = "wall"
         !! West boundary: "wall", "open", "tidal", "nested"
      character(len=16) :: bc_east = "wall"
         !! East boundary: "wall", "open", "tidal", "nested"
      character(len=16) :: bc_south = "wall"
         !! South boundary: "wall", "open", "tidal", "nested"
      character(len=16) :: bc_north = "wall"
         !! North boundary: "wall", "open", "tidal", "nested"

      ! Tidal forcing parameters
      real(wp) :: tidal_amplitude = 0.0_wp
         !! Tidal elevation amplitude (m) — single-constituent legacy
      real(wp) :: tidal_period = 44712.0_wp
         !! Tidal period (s), default M2 = 12.42 h — single-constituent legacy

      ! Multi-constituent tidal forcing
      integer :: n_tidal_constituents = 0
         !! Number of active tidal constituents (0 = use legacy single)
      real(wp) :: tidal_amp(MAX_TIDAL_CONSTITUENTS) = 0.0_wp
         !! Constituent amplitudes (m)
      real(wp) :: tidal_phase(MAX_TIDAL_CONSTITUENTS) = 0.0_wp
         !! Constituent phases (radians)
      real(wp) :: tidal_omega(MAX_TIDAL_CONSTITUENTS) = 0.0_wp
         !! Constituent angular frequencies (rad/s)

      ! Inflow/discharge/clamped BC parameters
      real(wp) :: inflow_velocity = 0.0_wp
         !! Prescribed normal velocity for inflow BC (m/s)
      real(wp) :: inflow_salinity = -1.0_wp
         !! Prescribed salinity for inflow BC (PSU); <0 = zero-gradient
      real(wp) :: inflow_temperature = -999.0_wp
         !! Prescribed temperature for inflow BC (degC); <0 = zero-gradient
      real(wp) :: discharge_flux = 0.0_wp
         !! Prescribed volume flux for discharge BC (m^2/s)
      real(wp) :: clamped_elevation = 0.0_wp
         !! Prescribed surface elevation for clamped BC (m)

      ! Sponge layer parameters
      integer :: sponge_width = 0
         !! Sponge layer width in cells
      real(wp) :: sponge_strength = 0.0_wp
         !! Sponge relaxation rate (1/s)
      real(wp) :: sponge_h_ref = 0.0_wp
         !! Sponge reference depth (m)

      ! Metadata parameters
      character(len=64) :: reference_time = ""
         !! Time origin for CF output, e.g. "2026-01-01 00:00:00" (empty = idealized run)
      character(len=256) :: institution = ""
         !! Institution name for CF global attribute (empty = omitted)

      ! Non-hydrostatic pressure correction
      logical :: use_nonhydrostatic = .false.
         !! Enable multi-layer NH pressure correction
      integer :: nz_layers = 2
         !! Number of sigma layers for NH solver

      ! Coupled hydrostatic multi-layer
      logical :: use_multilayer = .false.
         !! Enable coupled hydrostatic vertical layers (alternative to NH)
      real(wp) :: rho_0 = 1000.0_wp
         !! Reference density for EOS (kg/m^3)
      real(wp) :: vmix_nu_max = 0.01_wp
         !! Pacanowski-Philander max turbulent viscosity (m^2/s)
      real(wp) :: vmix_nu_bg = 1.0e-4_wp
         !! Background vertical viscosity (m^2/s)
      real(wp) :: vmix_ri_alpha = 5.0_wp
         !! Richardson number scaling coefficient (shared between momentum and tracers)
      logical :: vmix_use_kpp = .false.
         !! Enable the KPP (Large et al. 1994) surface-boundary-layer
         !! overlay on top of the PP81 interior closure.  Default off:
         !! the coastal multilayer path runs pure PP81 (bit-identical
         !! to pre-KPP).  When on, the solver swaps the three PP81
         !! vmix calls for the KPP kernels (BL-depth + shape-function
         !! enhanced diffusivity + non-local transport).
      logical :: vmix_use_keps = .false.
         !! Enable the two-equation k-epsilon vertical-mixing closure
         !! (rki_ml_keps).  Default off = PP81.  When on it overrides
         !! both PP81 and the KPP overlay: the solver advances k/eps and
         !! applies the k-epsilon eddy viscosity to momentum (tracer
         !! mixing is wired in a later increment).
      character(len=16) :: keps_stability = "constant"
         !! k-epsilon stability functions: "constant" (default; c_mu=0.09,
         !! Pr_t=1), "galperin" (quasi-equilibrium), or "canuto" (Canuto-A).
      integer :: keps_si_picard = 1
         !! Picard passes of the k<->eps coupling per semi-implicit step
         !! (`time_scheme="semi_implicit"` + `vmix_use_keps` only; the
         !! explicit path ignores it).  1 (default) = single pass, the
         !! shipped explicit cadence — steady-state-exact per the Q1
         !! prototype (local_archive/prototypes/si_closures).  2 (the certified
         !! setting; up to 3 accepted) re-runs the k/eps solve with the
         !! refreshed nu_t production coefficient for transient-fidelity
         !! studies.  Fixed pass count (GPU column-parallel, no per-column
         !! divergence).
      real(wp) :: kpp_ri_crit = 0.3_wp
         !! Critical bulk Richardson number for the KPP BL-depth sweep.
      real(wp) :: kpp_cs_nonlocal = 6.3_wp
         !! Non-local (counter-gradient) transport coefficient C_s
         !! (LMD94 eq 20, limit value).
      real(wp) :: kpp_c_vt2 = 1.8_wp
         !! Unresolved-turbulence coefficient for the V_t^2 term in the
         !! bulk-Ri denominator (LMD94 eq 23).  Set 0 to disable V_t^2.
      real(wp) :: hdiff_kappa = 0.0_wp
         !! Horizontal diffusion coefficient for tracers (m^2/s).
         !! Smooths sharp density fronts to reduce BPG overshoot.
         !! Typical values: 1-10 m^2/s depending on grid resolution.

      ! Tracer parameters (salinity + temperature)
      real(wp) :: initial_salinity = 35.0_wp
         !! Initial salinity for all layers (PSU).  Used as a uniform IC
         !! unless `S_init_surface` and `S_init_bottom` are BOTH set
         !! non-zero, in which case the seed builds a linear S(z) profile
         !! from `S_init_bottom` at k=1 (bed) to `S_init_surface` at
         !! k=nz (surface).  NOTE the stable polarity is the inverse of
         !! temperature: salty/dense water belongs at the bed, so a stable
         !! haline column has `S_init_bottom > S_init_surface` (fresh
         !! surface, e.g. a river-plume column).
      real(wp) :: S_init_surface = 0.0_wp
         !! Initial salinity at the surface layer (k=nz_ml, PSU).
         !! Activates linear S(z) stratification when set non-zero
         !! together with `S_init_bottom`.
      real(wp) :: S_init_bottom = 0.0_wp
         !! Initial salinity at the bed (k=1, PSU).
      real(wp) :: S_ref = 0.0_wp
         !! Reference salinity in EOS (PSU); density anomaly = beta_S*(S - S_ref)
      real(wp) :: beta_S = 0.78_wp
         !! Haline contraction coefficient (kg/m^3 per PSU)
      real(wp) :: S_min = 0.0_wp
         !! Lower physical bound for salinity (PSU)
      real(wp) :: S_max = 40.0_wp
         !! Upper physical bound for salinity (PSU)
      real(wp) :: kappa_S_bg = 1.0e-5_wp
         !! Background vertical salinity diffusivity (m^2/s)
      real(wp) :: initial_temperature = 15.0_wp
         !! Initial potential temperature for all layers (degC).  Used
         !! as a uniform IC unless `T_init_surface` and `T_init_bottom`
         !! are BOTH set non-zero, in which case the seed builds a
         !! linear T(z) profile from `T_init_bottom` at k=1 (bed) to
         !! `T_init_surface` at k=nz (surface).  Stratified IC is the
         !! prerequisite for baroclinic-instability-driven eddies.
      real(wp) :: T_init_surface = 0.0_wp
         !! Initial temperature at the surface layer (k=nz_ml, degC).
         !! Activates linear T(z) stratification when set non-zero
         !! together with `T_init_bottom`.
      real(wp) :: T_init_bottom = 0.0_wp
         !! Initial temperature at the bed (k=1, degC).
      real(wp) :: T_ref = 15.0_wp
         !! Reference temperature in EOS (degC); density anomaly = -alpha_T*(T - T_ref)
      real(wp) :: alpha_T = 0.17_wp
         !! Thermal expansion coefficient (kg/m^3 per degC); ~0.17 near 15 degC.
         !! Set to 0 to decouple temperature from dynamics (pure passive tracer).
      real(wp) :: T_min = -2.0_wp
         !! Lower physical bound for temperature (degC); seawater freezing
      real(wp) :: T_max = 40.0_wp
         !! Upper physical bound for temperature (degC)
      real(wp) :: kappa_T_bg = 1.0e-5_wp
         !! Background vertical temperature diffusivity (m^2/s)
      character(len=16) :: eos_form = "linear"
         !! COASTAL equation-of-state selector (both backends):
         !! "linear" (default) keeps the legacy coastal linear kernel
         !! `rho = rho_0 + beta_S*(S-S_ref) - alpha_T*(T-T_ref)` —
         !! bit-identical by construction for every existing config.
         !! "wright" (Wright 1997 rational) / "roquet_spv" (Roquet et
         !! al. 2015 TEOS-10 polynomial) route density through the
         !! SHARED nonlinear EOS (`rki_eos`), evaluated at p=0
         !! (surface/potential density); alpha_T/beta_S/S_ref/T_ref
         !! are then ignored by the density computation.  CAVEAT: the
         !! coastal KPP surface buoyancy-flux estimate `B_0 =
         !! (g/rho_0)(alpha_T F_T - beta_S F_S)` still reads the
         !! namelist alpha_T/beta_S as constant expansion coefficients
         !! — set them to locally representative values when pairing
         !! eos_form with vmix_use_kpp (interior Ri mixing reads
         !! rho_layer directly and is fully consistent).  The ocean
         !! regime ignores this knob (`&ocean_pgf_nml eos=` instead).
      logical :: use_ppm_tracer = .true.
         !! Use high-order tracer advection: PPM (structured) or MUSCL
         !! (unstructured). When false, falls back to first-order upwind.
      character(len=8) :: tracer_recon = "ppm"
         !! Tracer horizontal-advection face reconstruction (structured coastal):
         !! "ppm" (default, CW84 PPM), "weno5" (Jiang-Shu 1996 + Borges 2008,
         !! nghost>=3), "weno7" (Balsara-Shu 2000, nghost>=4), "weno9"
         !! (Balsara-Shu 2000, nghost>=5). "pqm" is deferred.
         !! Bit-identical to ppm when tracer_recon="ppm".
      real(wp) :: smagorinsky_cs = 0.0_wp
         !! Smagorinsky horizontal viscosity coefficient for layer momentum.
      real(wp) :: biharmonic_nu4 = 0.0_wp
         !! Biharmonic (del^4) lateral viscosity coefficient (m^4/s) for layer
         !! momentum; 0 = off.  Capped at the explicit-step stability limit.
         !! nu_h = (Cs * dx)^2 * |S| where |S| is the strain rate.
         !! Provides grid-aware dissipation at sharp velocity gradients.
      logical :: evap_mass_sink = .false.
         !! Treat net evaporation (evap_minus_precip > 0) as a real surface
         !! water-mass sink at the top layer (k=nz) instead of a virtual salt
         !! flux.  Removes water: h_top -= E*dt; salinity mass hS is left
         !! unchanged so concentration rises; heat leaves with the water
         !! (hT_top -= E*dt*T_top) so temperature concentration is preserved.
         !! Default .false. = bit-identical (legacy virtual-salt-flux path).
      logical :: evap_latent_cooling = .false.
         !! Latent-heat cooling of the surface from evaporation (v2).  When
         !! true (and evap_mass_sink is also true), the water that evaporates
         !! withdraws its latent heat of vaporization from the top layer:
         !! an ADDITIONAL surface heat sink Q_lat = rho_water*E_eff*L_v [W/m²]
         !! applied to k=nz, where E_eff is the effective (clamped) evaporation
         !! rate actually removed.  This is distinct from — and does NOT
         !! double-count — the sensible heat that already leaves with the
         !! removed water mass (that is handled by the proportional hT scaling
         !! in evap_mass_sink).  Requires evap_mass_sink=.true.; no-op
         !! otherwise (fail-loud at configure).  Default .false. = bit-identical.

      ! Sediment transport (single class, gated; default off = bit-identical).
      logical :: use_sediment = .false.
         !! Master switch for suspended-sediment transport.
      real(wp) :: sediment_w_s = 0.0_wp
         !! Settling (sinking) velocity (m/s, downward; 0 = no settling).
      real(wp) :: sediment_initial_conc = 0.0_wp
         !! Uniform initial suspended-sediment concentration (kg/m^3).
      real(wp) :: sediment_bed_init = 0.0_wp
         !! Uniform initial bed sediment mass (kg/m^2).
      real(wp) :: sediment_tau_crit_erosion = 0.0_wp
         !! Critical bed shear stress for erosion (Pa). 0 = erosion off.
      real(wp) :: sediment_erosion_rate = 0.0_wp
         !! Ariathurai erosion-rate coefficient E0 (kg/m^2/s). 0 = erosion off.
      real(wp) :: sediment_porosity = 0.0_wp
         !! Bed porosity (0-1); erosion scaled by (1 - porosity).
         !! 0 = off, typical 0.1-0.2.
      character(len=16) :: bpg_method = "jacobian"
         !! BPG computation method: "jacobian" (Shchepetkin-McWilliams,
         !! DEFAULT) or "sigma" (cheaper, but see below).
         !!
         !! `"sigma"` differences the integrated density anomaly at constant
         !! SIGMA INDEX.  Over sloping bathymetry layer k sits at a different
         !! physical depth in each column, so a spatially UNIFORM density
         !! whose value differs from `rho_0` yields a non-zero gradient --
         !! a fictitious force proportional to (uniform anomaly) x (bed
         !! slope).  Measured on a closed box over a sinusoidal bed, 400
         !! steps: relative mass change 1.7e-6 with "sigma" vs 1.3e-14 with
         !! "jacobian".  Retained only for reproducing legacy runs; the
         !! identity it violates is pinned by `test_bpg_uniform_density`.
      real(wp) :: ml_min_depth = 0.5_wp
         !! Minimum total water depth for vertical layer operations (m).
         !! Columns shallower than this skip vertical exchange, mixing,
         !! salinity advection, and BPG.  Increase if instabilities appear
         !! near wetting/drying fronts.
      real(wp) :: nh_tolerance = 1.0e-6_wp
         !! CG convergence tolerance for NH pressure
      integer :: nh_max_iterations = 100
         !! Maximum CG iterations for NH pressure
      logical :: nh_preconditioner = .false.
         !! Enable vertical-line preconditioner for CG solver
      integer :: nh_solve_interval = 1
         !! Solve NH pressure every N steps (1 = every step, N > 1 = reuse q_nh)

      ! Initial condition selection
      character(len=32) :: init_condition = "dam_break"
         !! Initial condition: "dam_break", "gaussian_hump", "castro_dam_break", or "still_water"
      real(wp) :: eta_initial = 0.0_wp
         !! Initial free-surface elevation for still_water IC (m, relative to bathymetry datum)
      real(wp) :: h0 = 1.0_wp
         !! Background depth for gaussian_hump IC (m)
      real(wp) :: hump_amplitude = 0.05_wp
         !! Perturbation amplitude for gaussian_hump IC (m)
      real(wp) :: hump_sigma = 5.0_wp
         !! Gaussian width for gaussian_hump IC (m)
      real(wp) :: hump_x_offset = 0.0_wp
         !! Hump centre x-offset from domain centre (m)
      real(wp) :: hump_y_offset = 0.0_wp
         !! Hump centre y-offset from domain centre (m)
      real(wp) :: init_u_shear = 0.0_wp
         !! Initial vertical shear: top layer gets +val, bottom gets -val (m/s)

      ! Mode splitting parameters
      integer :: split_ratio = 1
         !! Barotropic-baroclinic split ratio: 0=auto, 1=unsplit, N>1=fixed
      logical :: split_averaging = .true.
         !! ROMS-style time averaging over barotropic substeps

      ! Vertical coordinate parameters
      character(len=16) :: vcoord_type = "sigma"
         !! Vertical coordinate type: "sigma", "zsigma", "zstar", etc.
      character(len=16) :: remap_method = "ppm"
         !! Vertical remapping method: "pcm", "plm", "ppm", "ppm_h4", "pqm".  Default PPM
         !! (piecewise parabolic, Colella & Woodward 1984) — the state-of-the-art
         !! conservative vertical-remap stencil for the coastal ALE coords. "pqm" is
         !! piecewise quartic (White & Adcroft 2008); it falls back to PPM for nz < 5.
      real(wp) :: vcoord_depth_transition = 100.0_wp
         !! z-sigma hybrid: depth (m) where blending starts
      real(wp) :: vcoord_blend_width = 50.0_wp
         !! z-sigma hybrid: blending zone width (m)
      ! Full MOM6 z* parameters (VCOORD_ZSTAR_FULL)
      real(wp) :: zstar_h_surf_target = 0.0_wp
         !! Target physical thickness of the surface layer (m).
         !! 0 => auto: fall back to the lite (uniform sigma) profile
         !! per column (effectively zstar-lite per-column).  Set > 0 to
         !! anchor the surface layer at a fixed thickness regardless of H.
      character(len=16) :: zstar_stretching = "log"
         !! Surface-concentration stretching: "log" | "uniform"
      real(wp) :: zstar_h_min = 1.0e-4_wp
         !! Vanishing-layer floor (m).  Layers that would land below the
         !! local bed get clipped to this thickness rather than going to zero.
      integer :: zstar_n_surf = 0
         !! Number of "fine" near-surface layers using stretching.
         !! 0 => auto (use max(1, nz/3))
      ! Isopycnal coordinate (VCOORD_RHO) parameters
      real(wp) :: rho_ref_pressure = 2.0e7_wp
         !! Reference pressure (Pa, default 2e7 = 2000 dbar) for the
         !! potential density that defines the VCOORD_RHO coordinate.
      real(wp) :: rho_target_light = 1020.0_wp
         !! Lightest target interface potential density (kg/m³) — the
         !! surface interface for VCOORD_RHO.  With `rho_target_dense`
         !! builds the uniform light→dense target-density linspace
         !! (MOM6 `ALE_COORDINATE_CONFIG=UNIFORM` analogue for the RHO
         !! mode).  A stretched (RFNC1-style) profile is a P3 follow-up.
      real(wp) :: rho_target_dense = 1030.0_wp
         !! Densest target interface potential density (kg/m³) — the bed
         !! interface for VCOORD_RHO.  See `rho_target_light`.

      ! Logging parameters
      character(len=16) :: log_level = "info"
         !! Log verbosity: "debug", "verbose", "info", "performance", "warning", "error"
      real(wp) :: status_interval = 0.0_wp
         !! Print status every N seconds of simulation time (0 = default every 100 steps)

   end type config_t

contains

   subroutine read_config(filename, cfg, schema)
      !! Read simulation configuration from a namelist file.
      !!
      !! Groups still read natively (`read(nml=)`) are seeded from the
      !! `cfg` defaults, read, and copied back as before.  Groups that
      !! have migrated to the strict schema (`rki_nml_schema`) are
      !! validated and applied here, AFTER the native reads, by building
      !! a local schema (defaults captured from `cfg`'s pristine field
      !! initialisers) and strict-parsing the file — so `error stop` on
      !! any unknown key / range / enum violation fires for every caller
      !! (driver, FFI, unit tests), not just the production entry points.
      !!
      !! When `schema` is present the built schema is returned for the
      !! caller's parameter-doc dumps.  Its key pointers alias into
      !! `cfg`; both entry points keep their `cfg` alive (and `target`),
      !! so returning it is safe.
      character(len=*), intent(in) :: filename
         !! Path to the namelist input file
      type(config_t), target, intent(out) :: cfg
         !! Populated configuration
      type(nml_schema_t), intent(out), optional :: schema
         !! Built + parsed schema (aliases into `cfg`); for doc dumps.

      ! Dispatch so the schema is built DIRECTLY into the caller's
      ! variable when requested — never copied.  Intrinsic assignment of
      ! nml_schema_t deep-copies polymorphic allocatable key boxes, which
      ! NVHPC miscompiles (heap corruption in the FFI tests); building in
      ! place sidesteps it and is cheaper anyway.
      if (present(schema)) then
         call read_config_impl(filename, cfg, schema)
      else
         block
            type(nml_schema_t) :: sch
            call read_config_impl(filename, cfg, sch)
         end block
      end if

   end subroutine read_config

   subroutine read_config_impl(filename, cfg, sch)
      !! Body of [[read_config]]: build the strict schema (capturing
      !! defaults from the pristine `cfg` field initialisers), parse the
      !! file (validate + apply straight into `cfg`), run the post-parse
      !! time-unit cascade.
      character(len=*), intent(in) :: filename
         !! Path to the namelist input file
      type(config_t), target, intent(out) :: cfg
         !! Populated configuration
      type(nml_schema_t), intent(out) :: sch
         !! Built + parsed schema (aliases into `cfg`)

      ! `cfg` arrives at its pristine field-initialiser defaults.  Every
      ! namelist group is now on the strict schema: build it (capturing
      ! those defaults), then parse — which validates AND applies every
      ! override straight into `cfg` via the schema's key pointers.  No
      ! native readers, no flat shadow locals, no copy-back remain.

      ! Build the strict schema so its key constructors capture defaults
      ! from the pristine `cfg` field initialisers.
      call build_rakali_schema(cfg, sch)

      ! A missing namelist file is not an error — the pristine `cfg`
      ! defaults stand (production main.F90 fails loudly on a missing
      ! file earlier; tests + the FFI rely on this defaults-on-absent
      ! behaviour).  Guard the strict parse so it never error-stops on
      ! "cannot open file".
      block
         logical :: file_exists
         inquire (file=trim(filename), exist=file_exists)
         if (.not. file_exists) then
            call logger%debug("Namelist file not found, keeping defaults: "//trim(filename))
            return
         end if
      end block

      ! Strict schema parse: validate + APPLY every group straight into
      ! `cfg`.  Strict (no status/errors args) → on any unknown key /
      ! range / enum violation it logs the full report and error-stops,
      ! for EVERY caller (driver, FFI, unit tests).

      call sch%parse(trim(filename))

      call apply_time_unit_cascade(cfg)

   end subroutine read_config_impl

   subroutine read_config_from_string(text, cfg, schema)
      !! In-memory sibling of [[read_config]]: build + strict-parse +
      !! apply a namelist held entirely in the `text` buffer (newline-
      !! separated), with NO filesystem touch — no temp file, no `chdir`,
      !! no cleanup.  Same validation + error-stop behaviour as the
      !! file path (the schema parser runs from the in-memory line array).
      character(len=*), intent(in) :: text
         !! Whole namelist as one string (records separated by '\n').
      type(config_t), target, intent(out) :: cfg
      type(nml_schema_t), intent(out), optional :: schema

      if (present(schema)) then
         call read_config_from_string_impl(text, cfg, schema)
      else
         block
            type(nml_schema_t) :: sch
            call read_config_from_string_impl(text, cfg, sch)
         end block
      end if
   end subroutine read_config_from_string

   subroutine read_config_from_string_impl(text, cfg, sch)
      !! Body of [[read_config_from_string]] — mirrors read_config_impl
      !! with the line array replacing the file unit.
      character(len=*), intent(in) :: text
      type(config_t), target, intent(out) :: cfg
      type(nml_schema_t), intent(out) :: sch

      character(len=:), allocatable :: lines(:)
      integer :: n_lines

      call build_rakali_schema(cfg, sch)

      ! Empty buffer → keep pristine defaults (mirrors the missing-file
      ! branch of read_config_impl).
      if (len_trim(text) == 0) return

      call split_to_lines(text, lines, n_lines)

      ! Strict schema parse from the same lines: validate + APPLY every
      ! migrated group straight into `cfg`.
      call sch%parse_lines(lines, n_lines)

      call apply_time_unit_cascade(cfg)
   end subroutine read_config_from_string_impl

   subroutine split_to_lines(text, lines, n_lines)
      !! Split a newline-separated buffer into a fixed-len character array
      !! (one record per line, trailing CR stripped) for internal-file
      !! reads + the schema line walker.
      character(len=*), intent(in) :: text
      character(len=:), allocatable, intent(out) :: lines(:)
      integer, intent(out) :: n_lines

      integer :: i, start, maxlen, cur, nl
      character(len=1) :: c

      ! First pass: count records + longest line.
      nl = 1
      maxlen = 0
      cur = 0
      do i = 1, len(text)
         c = text(i:i)
         if (c == achar(10)) then
            nl = nl + 1
            maxlen = max(maxlen, cur)
            cur = 0
         else if (c /= achar(13)) then
            cur = cur + 1
         end if
      end do
      maxlen = max(maxlen, cur, 1)

      allocate (character(len=maxlen) :: lines(nl))
      lines = ""

      ! Second pass: copy each record.
      n_lines = 0
      start = 1
      cur = 0
      do i = 1, len(text)
         c = text(i:i)
         if (c == achar(10)) then
            n_lines = n_lines + 1
            if (cur > 0) lines(n_lines) = text(start:start + cur - 1)
            start = i + 1
            cur = 0
         else if (c /= achar(13)) then
            ! CR is skipped (a trailing CR before LF is excluded because
            ! `cur` counts only content chars and `start` marks the first).
            if (cur == 0) start = i
            cur = cur + 1
         end if
      end do
      ! Trailing record (no final newline).
      n_lines = n_lines + 1
      if (cur > 0) lines(n_lines) = text(start:start + cur - 1)
   end subroutine split_to_lines

   subroutine apply_time_unit_cascade(cfg)
      !! Post-parse time-unit fixup shared by the file + string config
      !! paths.  `t_end` and `status_interval` are given
      !! in `time_unit` from `&time_nml`; multiply through to SI seconds.
      !! `dt_fixed`/`dt_max` stay in seconds; default `time_unit="s"`
      !! gives a 1× factor (bit-identical to before).
      type(config_t), intent(inout) :: cfg

      real(wp) :: time_factor

      select case (trim(cfg%time_unit))
      case ("s", "sec", "second", "seconds")
         time_factor = 1.0_wp
      case ("min", "minute", "minutes")
         time_factor = 60.0_wp
      case ("hr", "hour", "hours", "h")
         time_factor = 3600.0_wp
      case ("day", "days", "d")
         time_factor = 86400.0_wp
      case ("year", "years", "yr")
         ! Julian year (365.25 d) — matches MOM6 / UDUNITS.
         time_factor = 365.25_wp*86400.0_wp
      case default
         call logger%error("Invalid time_unit = '"//trim(cfg%time_unit)// &
                           "': must be one of s/min/hr/day/year")
         time_factor = 1.0_wp  ! fallthrough so validation reports it
      end select
      cfg%t_end = cfg%t_end*time_factor
      cfg%status_interval = cfg%status_interval*time_factor

      ! Derive the logical flag so downstream code never compares strings.
      cfg%semi_implicit = (trim(adjustl(cfg%time_scheme)) == "semi_implicit")
   end subroutine apply_time_unit_cascade

   subroutine validate_config(cfg)
      !! Validate configuration parameters after reading
      !!
      !! Logs errors for invalid values that would crash the solver,
      !! and warnings for suspicious but non-fatal settings.
      use rki_ml_common_unstr, only: BPG_NZ_MAX, bpg_nz_required, &
                                     bpg_nz_is_sufficient
      type(config_t), intent(in) :: cfg

      logical :: has_error

      has_error = .false.

      ! Simulation regime — this build carries the coastal path only;
      ! any other regime tag must fail loud, never silently fall through.
      if (trim(cfg%sim_type) /= "coastal") then
         call logger%error("Invalid sim_type = '"//trim(cfg%sim_type)// &
                           "': must be 'coastal' (the only regime in this build)")
         has_error = .true.
      end if

      ! Grid parameters
      if (cfg%nx < 1) then
         call logger%error("Invalid nx = "//to_string(cfg%nx)//": must be >= 1")
         has_error = .true.
      end if
      if (cfg%ny < 1) then
         call logger%error("Invalid ny = "//to_string(cfg%ny)//": must be >= 1")
         has_error = .true.
      end if
      if (cfg%nghost < 1) then
         call logger%error("Invalid nghost = "//to_string(cfg%nghost)//": must be >= 1")
         has_error = .true.
      end if
      if (cfg%dx <= 0.0_wp) then
         call logger%error("Invalid dx = "//to_string(cfg%dx)//": must be > 0")
         has_error = .true.
      end if
      if (cfg%dy <= 0.0_wp) then
         call logger%error("Invalid dy = "//to_string(cfg%dy)//": must be > 0")
         has_error = .true.
      end if

      ! Time parameters
      if (cfg%cfl <= 0.0_wp .or. cfg%cfl > 1.0_wp) then
         call logger%error("Invalid cfl = "//to_string(cfg%cfl)//": must be in (0, 1]")
         has_error = .true.
      end if
      if (cfg%t_end <= 0.0_wp) then
         call logger%error("Invalid t_end = "//to_string(cfg%t_end)//": must be > 0")
         has_error = .true.
      end if

      ! Multilayer parameters
      if (cfg%use_multilayer) then
         if (cfg%nz_layers < 1) then
            call logger%error("Invalid nz_layers = "//to_string(cfg%nz_layers)// &
                              ": must be >= 1 for multilayer")
            has_error = .true.
         end if
         if (cfg%use_nonhydrostatic) then
            call logger%error("use_multilayer and use_nonhydrostatic are mutually exclusive")
            has_error = .true.
         end if
         ! Unstructured multilayer: the per-column stack workspaces are
         ! sized by BPG_NZ_MAX (rki_ml_common_unstr, default 20), NOT by
         ! NZ_STACK_MAX.  The binding consumer (the k-eps interface
         ! tridiagonal) indexes plain BPG_NZ_MAX arrays at nz + 1, so
         ! the requirement is nz + 1 <= BPG_NZ_MAX.  Exceeding it is not
         ! a crash: the unstructured BPG kernels silently return (no
         ! baroclinic pressure gradient — plausible-looking wrong
         ! physics) and the unguarded column kernels silently overrun
         ! thread-local storage.  Refuse here, fail-loud.
         if (trim(cfg%grid_type) == "unstructured" .and. &
             .not. bpg_nz_is_sufficient(cfg%nz_layers)) then
            call logger%error("Unstructured multilayer nz_layers = "// &
                              to_string(cfg%nz_layers)//" needs BPG_NZ_MAX >= "// &
                              to_string(bpg_nz_required(cfg%nz_layers))// &
                              " but this build has BPG_NZ_MAX = "// &
                              to_string(BPG_NZ_MAX)//" (the k-eps interface "// &
                              "tridiagonal indexes per-column stack arrays at "// &
                              "nz + 1). Reduce nz_layers to <= "// &
                              to_string(BPG_NZ_MAX - 1)//" or rebuild with "// &
                              "-DRAKALI_BPG_NZ_MAX="// &
                              to_string(bpg_nz_required(cfg%nz_layers))//" or larger.")
            has_error = .true.
         end if
      end if

      ! Per-column stack-workspace bound.  FAIL LOUD.
      !
      ! Every layered kernel (BPG, ALE remap, the diag-remap fills)
      ! carries fixed-size `NZ_STACK_MAX` thread-local column arrays.
      ! Overrunning them corrupts thread-local storage SILENTLY — wrong
      ! answers, no crash — so this must refuse the run, not warn.
      !
      ! The requirement is `nz + 1` — see `nz_stack_required`.
      if (cfg%use_multilayer) then
         if (.not. nz_stack_is_sufficient(cfg%nz_layers)) then
            call logger%error("nz_layers = "//to_string(cfg%nz_layers)// &
                              " needs NZ_STACK_MAX >= "// &
                              to_string(nz_stack_required(cfg%nz_layers))// &
                              " but this binary was compiled with NZ_STACK_MAX = "// &
                              to_string(NZ_STACK_MAX)// &
                              ". Per-column stack kernels (BPG, ALE remap) would "// &
                              "overrun thread-local storage and "// &
                              "silently produce wrong answers. Rebuild with "// &
                              "-DRAKALI_NZ_STACK_MAX="// &
                              to_string(nz_stack_required(cfg%nz_layers))//" (or larger).")
            has_error = .true.
         end if
      end if

      ! tracer_recon knob guards.  The "is this rung supported on this
      ! config" decision runs for EVERY config (fail-loud): a weno* rung is
      ! silently ignored by anything but the coastal multilayer solver, so
      ! it must abort on non-multilayer configs (2D barotropic / NH)
      ! rather than parse-and-drop.  The nghost /
      ! use_ppm_tracer quality guards apply only to the valid weno path.
      block
         use rki_recon_weno, only: parse_tracer_recon, &
                                   tracer_recon_required_nghost, &
                                   tracer_recon_support_status, &
                                   TRACER_RECON_PPM, &
                                   TRACER_RECON_OK, &
                                   TRACER_RECON_REJECT_PQM, &
                                   TRACER_RECON_REJECT_UNKNOWN, &
                                   TRACER_RECON_REJECT_OCEAN, &
                                   TRACER_RECON_REJECT_NO_ML
         integer :: recon_code_local, support_status
         recon_code_local = parse_tracer_recon(trim(cfg%tracer_recon))
         ! Second argument is the "is this a non-coastal regime" flag of the
         ! shared predicate — always .false. here (coastal is the only regime).
         support_status = tracer_recon_support_status(recon_code_local, &
                                                      .false., &
                                                      cfg%use_multilayer)
         select case (support_status)
         case (TRACER_RECON_REJECT_PQM)
            call logger%error("tracer_recon='pqm' is deferred: dominated by "// &
                              "weno5-Z under practical limiters; use weno5/7/9.")
            has_error = .true.
         case (TRACER_RECON_REJECT_UNKNOWN)
            call logger%error("tracer_recon='"//trim(cfg%tracer_recon)// &
                              "' not recognised; valid: ppm, weno5, weno7, weno9")
            has_error = .true.
         case (TRACER_RECON_REJECT_OCEAN)
            ! Unreachable here (the regime flag above is hard .false.);
            ! kept for select-case totality over the shared predicate's codes.
            call logger%error("tracer_recon='"//trim(cfg%tracer_recon)// &
                              "' is only supported by the coastal multilayer solver")
            has_error = .true.
         case (TRACER_RECON_REJECT_NO_ML)
            call logger%error("tracer_recon='"//trim(cfg%tracer_recon)// &
                              "' requires multilayer tracers (use_multilayer=.true.)")
            has_error = .true.
         case (TRACER_RECON_OK)
            ! Valid path.  For a weno rung on coastal multilayer, apply the
            ! multilayer-specific quality guards.
            if (recon_code_local /= TRACER_RECON_PPM) then
               if (cfg%nghost < tracer_recon_required_nghost(recon_code_local)) then
                  call logger%error("tracer_recon='"//trim(cfg%tracer_recon)// &
                                    "' requires nghost >= "// &
                                    to_string(tracer_recon_required_nghost(recon_code_local))// &
                                    " but nghost = "//to_string(cfg%nghost))
                  has_error = .true.
               end if
               if (.not. cfg%use_ppm_tracer) then
                  call logger%error("tracer_recon='"//trim(cfg%tracer_recon)// &
                                    "' requires use_ppm_tracer=.true.")
                  has_error = .true.
               end if
            end if
         case default
            ! Unreachable: tracer_recon_support_status returns only the
            ! codes above.  Fail loud if a new code is ever added unhandled.
            call logger%error("internal error: unhandled tracer_recon support status "// &
                              to_string(support_status))
            has_error = .true.
         end select
      end block

      ! NH parameters
      if (cfg%use_nonhydrostatic) then
         if (cfg%nz_layers < 1) then
            call logger%error("Invalid nz_layers = "//to_string(cfg%nz_layers)// &
                              ": must be >= 1")
            has_error = .true.
         end if
         if (cfg%nh_tolerance <= 0.0_wp) then
            call logger%error("Invalid nh_tolerance = "//to_string(cfg%nh_tolerance)// &
                              ": must be > 0")
            has_error = .true.
         end if
         if (cfg%nh_max_iterations < 1) then
            call logger%error("Invalid nh_max_iterations = "// &
                              to_string(cfg%nh_max_iterations)//": must be >= 1")
            has_error = .true.
         end if
      end if

      ! Semi-implicit free-surface mode: 2D barotropic (Phase 0/1) or 3D
      ! multilayer (Phase 2).  Unstructured grids only; non-hydrostatic is
      ! excluded.
      if (cfg%semi_implicit) then
         if (semi_implicit_grid_conflict(cfg)) then
            call logger%error("time_scheme='semi_implicit' is only supported "// &
                              "on unstructured grids (grid_type='unstructured').")
            has_error = .true.
         end if
         if (cfg%use_nonhydrostatic) then
            call logger%error("time_scheme='semi_implicit' is incompatible with "// &
                              "use_nonhydrostatic=.true. (the implicit free surface "// &
                              "assumes hydrostatic; NH-on-semi-implicit is deferred).")
            has_error = .true.
         end if
         if (cfg%dt_fixed <= 0.0_wp) then
            call logger%error("time_scheme='semi_implicit' requires dt_fixed > 0 "// &
                              "(the semi-implicit scheme bypasses the gravity-wave CFL and takes "// &
                              "the fixed step directly).")
            has_error = .true.
         end if
         ! Every non-wall edge must be a Dirichlet-eta type the SI seed
         ! supports (tidal / clamped).  Any other non-wall type (open/
         ! Flather, chapman, sponge, nested, inflow, discharge) fills its
         ! ghosts with momentum the Helmholtz operator treats as a no-flux
         ! wall (Hf=0) — the RHS divergence then sees flux through a "wall"
         ! with no balancing PGF: a silent mass leak.  Abort instead.
         if (.not. si_bc_supported(cfg%bc_west)) then
            call logger%error("time_scheme='semi_implicit': bc_west='"// &
                              trim(cfg%bc_west)//"' is unsupported (only wall / tidal / "// &
                              "clamped are wired into the implicit eta-solve).")
            has_error = .true.
         end if
         if (.not. si_bc_supported(cfg%bc_east)) then
            call logger%error("time_scheme='semi_implicit': bc_east='"// &
                              trim(cfg%bc_east)//"' is unsupported (only wall / tidal / "// &
                              "clamped are wired into the implicit eta-solve).")
            has_error = .true.
         end if
         if (.not. si_bc_supported(cfg%bc_south)) then
            call logger%error("time_scheme='semi_implicit': bc_south='"// &
                              trim(cfg%bc_south)//"' is unsupported (only wall / tidal / "// &
                              "clamped are wired into the implicit eta-solve).")
            has_error = .true.
         end if
         if (.not. si_bc_supported(cfg%bc_north)) then
            call logger%error("time_scheme='semi_implicit': bc_north='"// &
                              trim(cfg%bc_north)//"' is unsupported (only wall / tidal / "// &
                              "clamped are wired into the implicit eta-solve).")
            has_error = .true.
         end if
         ! Phase-2 (3D multilayer) scope guards — fail loud rather than
         ! silently drop a configured feature the SI step does not run.
         ! The unstructured 3D SI path (P3.2) is sigma-only /
         ! wall-tidal-clamped / PP81-or-closure / no-sediment.
         if (cfg%use_multilayer) then
            if (cfg%nz_layers < 2) then
               call logger%error("time_scheme='semi_implicit' with use_multilayer "// &
                                 "requires nz_layers >= 2 (the nz=1 multilayer fast path "// &
                                 "is not wired to the SI column solve; use the pure 2D "// &
                                 "barotropic SI mode instead).")
               has_error = .true.
            end if
            ! Boundary-type support (wall / tidal / clamped only) is already
            ! enforced for EVERY SI config — 2D and ML alike — by the
            ! all-edges si_bc_supported guard above; no ML-specific re-check
            ! is needed (the seed/operator machinery is shared).
            ! k-epsilon and KPP are wired into the SI implicit column solve
            ! (Q1: closure fills si%nu_iface; the implicit tridiagonal + the
            ! implicit tracer vdiff consume it) — no guard needed for them.
            if (cfg%use_sediment) then
               call logger%error("time_scheme='semi_implicit' with use_multilayer "// &
                                 "does not run the sediment settling/erosion kernels "// &
                                 "(Phase 2); use_sediment would be silently dropped.")
               has_error = .true.
            end if
            if (trim(cfg%vcoord_type) /= "sigma") then
               call logger%error("time_scheme='semi_implicit' with use_multilayer "// &
                                 "targets vcoord_type='sigma' in Phase 2 (the remap "// &
                                 "vcoords need the pre-remap Lagrangian-thickness "// &
                                 "consistency pass, deferred).")
               has_error = .true.
            end if
         end if

         ! nsub >= 2 is live via `wetdry_bsub_fill` ("flat" default =
         ! bit-identical to nsub=1; "linear" = the analytic Casulli 2009
         ! subgrid ramp).  The LINEAR fill is incompatible with the
         ! multilayer wet/dry step on BOTH backends: the eta-solve conserves
         ! the subgrid v_cell(eta) but the ML thickness update consumes
         ! max(eta-b,0) on the plain bed, so mass leaks at the wet/dry
         ! front.  The state-init si_bsub_ml_conflict guards error-stop too;
         ! this is the configure-time twin.
         if (cfg%si_wetdry .and. cfg%si_wetdry_nsub > 1 .and. &
             trim(adjustl(cfg%si_wetdry_bsub_fill)) == "linear" .and. &
             cfg%use_multilayer) then
            call logger%error("time_scheme='semi_implicit': wetdry_bsub_fill='linear' "// &
                              "(nsub>=2) is not supported with use_multilayer — use "// &
                              "'flat' with multilayer, or 'linear' on the 2D barotropic "// &
                              "path.")
            has_error = .true.
         end if

      end if

      ! `keps_si_picard > 1` is consumed only under semi_implicit + ML +
      ! vmix_use_keps; anywhere else it is silently inert => warn from
      ! top-level scope so explicit-scheme / 2D configs get feedback too
      ! (warning not error: default 1, harmless).
      if (cfg%keps_si_picard > 1 .and. .not. (cfg%semi_implicit .and. &
                                              cfg%use_multilayer .and. cfg%vmix_use_keps)) then
         call logger%warning("keps_si_picard > 1 has no effect without "// &
                             "time_scheme='semi_implicit' + use_multilayer + "// &
                             "vmix_use_keps = .true.")
      end if

      ! si_wetdry only takes effect when the semi-implicit slot allocates;
      ! with any other time_scheme the knob would be silently dropped.
      if (cfg%si_wetdry .and. .not. cfg%semi_implicit) then
         call logger%error("si_wetdry=.true. requires time_scheme='semi_implicit' "// &
                           "— the knob would be silently dropped (the semi-implicit "// &
                           "state slot only allocates under time_scheme='semi_implicit').")
         has_error = .true.
      end if

      ! Split ratio
      if (cfg%split_ratio < 0) then
         call logger%error("Invalid split_ratio = "//to_string(cfg%split_ratio)// &
                           ": must be >= 0")
         has_error = .true.
      end if

      ! Sponge width
      if (cfg%sponge_width < 0) then
         call logger%error("Invalid sponge_width = "//to_string(cfg%sponge_width)// &
                           ": must be >= 0")
         has_error = .true.
      end if

      ! Warn about unknown BC strings (they silently default to WALL)
      call warn_unknown_bc(cfg%bc_west, "bc_west")
      call warn_unknown_bc(cfg%bc_east, "bc_east")
      call warn_unknown_bc(cfg%bc_south, "bc_south")
      call warn_unknown_bc(cfg%bc_north, "bc_north")

      ! Warnings for suspicious settings.
      ! split_ratio defaults to 1 (= "unsplit, no subcycling"); only 0
      ! ("auto") or values > 1 actually engage mode splitting, and they
      ! only matter when nonhydrostatic mode is on.
      if (cfg%split_ratio /= 1 .and. .not. cfg%use_nonhydrostatic) then
         call logger%warning("split_ratio = "//to_string(cfg%split_ratio)// &
                             " has no effect without use_nonhydrostatic = .true.")
      end if

      if (cfg%output_to_file .and. cfg%output_interval <= 0.0_wp) then
         call logger%warning("output_to_file = .true. but output_interval <= 0 — no snapshots will be written")
      end if
      if (cfg%async_output) then
         call logger%warning("async_output = .true. is not supported in this build")
      end if

      ! ---- Evaporative latent-heat cooling (coastal multilayer) ----
      ! Latent cooling only makes physical sense when the evaporated water is
      ! actually removed (evap_mass_sink): the sink withdraws the latent heat
      ! of the mass being evaporated.  With evap_mass_sink off, no water leaves
      ! and the cooling would be an arbitrary, unbudgeted heat loss.  Fail loud
      ! rather than silently no-op.
      if (cfg%evap_latent_cooling .and. .not. cfg%evap_mass_sink) then
         call logger%error("evap_latent_cooling=.true. requires evap_mass_sink=.true. "// &
                           "(latent cooling is tied to the evaporated water mass being removed)")
         has_error = .true.
      end if

      if (has_error) then
         error stop "Configuration validation failed — see errors above"
      end if

   end subroutine validate_config

   pure logical function semi_implicit_grid_conflict(cfg) result(conflict)
      !! `.true.` when the config asks for `time_scheme="semi_implicit"` on
      !! anything but an unstructured grid.  Semi-implicit stepping exists
      !! only on the triangular backend; `validate_config` and the state
      !! initialisers reject the combination.
      type(config_t), intent(in) :: cfg
      conflict = (trim(adjustl(cfg%time_scheme)) == "semi_implicit") .and. &
                 (trim(adjustl(cfg%grid_type)) /= "unstructured")
   end function semi_implicit_grid_conflict

   subroutine warn_unknown_bc(bc_str, param_name)
      !! Warn if a BC string does not match any known type
      character(len=*), intent(in) :: bc_str
      character(len=*), intent(in) :: param_name

      select case (trim(bc_str))
      case ("wall", "open", "tidal", "nested", "inflow", &
            "discharge", "clamped", "sponge", "chapman")
         continue
      case default
         call logger%warning("Unknown boundary type '"//trim(bc_str)// &
                             "' for "//param_name//", defaulting to wall")
      end select

   end subroutine warn_unknown_bc

   pure logical function si_bc_supported(bc_str) result(ok)
      !! .true. iff the edge type is wired into the semi-implicit
      !! eta-solve (2D and 3D-ML paths alike): wall (mirror ghost =>
      !! exact zero face flux) or the two Dirichlet-eta seeds (tidal /
      !! clamped).  Everything else (open/Flather, chapman, sponge,
      !! nested, inflow, discharge) fills ghost momentum the Helmholtz
      !! operator walls off => a silent mass leak; validate_config
      !! aborts on those.
      character(len=*), intent(in) :: bc_str
      select case (trim(bc_str))
      case ("wall", "tidal", "clamped")
         ok = .true.
      case default
         ok = .false.
      end select
   end function si_bc_supported

   ! ==================================================================
   ! Strict-schema construction (rki_nml_schema).  Lives here (not in
   ! rki_config_schema) because read_config builds + parses the schema
   ! itself — a module-level `use rki_config_schema` from rki_config
   ! would be a circular dependency.  rki_config_schema re-exports
   ! build_rakali_schema so the production entry points keep their
   ! original import.
   ! ==================================================================

   subroutine build_rakali_schema(cfg, schema)
      !! Register the validated groups + the still-external groups onto
      !! `schema`, capturing defaults from `cfg`.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(out) :: schema

      call register_sim(cfg, schema)
      call register_grid(cfg, schema)
      call register_time(cfg, schema)
      call register_si(cfg, schema)
      call register_metadata(cfg, schema)
      call register_logging(cfg, schema)
      call register_vcoord(cfg, schema)
      call register_physics(cfg, schema)
      call register_output(cfg, schema)
      call register_boundary(cfg, schema)
      call register_tracer(cfg, schema)
      call register_initial_condition(cfg, schema)
      call register_nonhydrostatic(cfg, schema)
   end subroutine build_rakali_schema

   subroutine register_sim(cfg, schema)
      !! `&sim_nml`: the simulation regime selector.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      character(len=:), pointer :: ps

      g%name = "sim"
      g%doc = "Simulation-regime selector."
      ps => cfg%sim_type
      call g%add(nml_enum("sim_type", ps, &
                          "Simulation regime: coastal", &
                          allowed=[character(len=8) :: "coastal"]))
      call schema%add_group(g)
   end subroutine register_sim

   subroutine register_grid(cfg, schema)
      !! `&grid_nml`: structured-grid geometry.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      integer, pointer :: pi
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "grid"
      g%doc = "Structured-grid geometry."
      ps => cfg%grid_type
      call g%add(nml_enum("grid_type", ps, &
                          "Grid discretization: structured (Cartesian) or unstructured (triangular)", &
                          allowed=[character(len=12) :: "structured", "unstructured"]))
      pi => cfg%nx
      call g%add(nml_int("nx", pi, "Number of physical cells in x", min=1))
      pi => cfg%ny
      call g%add(nml_int("ny", pi, "Number of physical cells in y", min=1))
      pr => cfg%dx
      call g%add(nml_real("dx", pr, "Cell size in x", units="m"))
      pr => cfg%dy
      call g%add(nml_real("dy", pr, "Cell size in y", units="m"))
      pi => cfg%nghost
      call g%add(nml_int("nghost", pi, "Ghost cells on each side", min=1))
      call schema%add_group(g)
   end subroutine register_grid

   subroutine register_time(cfg, schema)
      !! `&time_nml`: time-integration controls.  `t_end` is interpreted
      !! in `time_unit` and converted to seconds by the post-parse
      !! cascade in read_config; `dt_fixed`/`dt_max` stay in seconds.
      !! `time_unit` is a plain string (the cascade validates its set).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      integer, pointer :: pi
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "time"
      g%doc = "Time-integration controls."
      pr => cfg%t_end
      call g%add(nml_real("t_end", pr, "Simulation end time (in time_unit)"))
      pr => cfg%cfl
      call g%add(nml_real("cfl", pr, "CFL number for adaptive timestep"))
      pr => cfg%dt_max
      call g%add(nml_real("dt_max", pr, "Maximum allowable timestep", units="s"))
      pr => cfg%dt_fixed
      call g%add(nml_real("dt_fixed", pr, "Fixed timestep (0 = adaptive CFL)", units="s"))
      pi => cfg%cfl_interval
      call g%add(nml_int("cfl_interval", pi, "Recompute CFL timestep every N steps", min=1))
      ps => cfg%time_scheme
      call g%add(nml_enum("time_scheme", ps, &
                          "Time integrator (coastal 2D): euler/rk2/rk3/rk43 (explicit) or "// &
                          "semi_implicit (unstructured grids only)", &
                          allowed=[character(len=13) :: "euler", "rk2", "rk3", "rk43", "semi_implicit"]))
      ps => cfg%time_unit
      call g%add(nml_string("time_unit", ps, &
                            "Unit for the long-time fields: s/min/hr/day/year"))
      call schema%add_group(g)
   end subroutine register_time

   subroutine register_si(cfg, schema)
      !! `&si_nml`: semi-implicit free-surface CG-solver controls.  Only
      !! consulted when `time_scheme="semi_implicit"`; both knobs default
      !! to the pre-knob plain-CG behaviour (byte-identical).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      character(len=:), pointer :: ps
      logical, pointer :: pl
      integer, pointer :: pi

      g%name = "si"
      g%doc = "Semi-implicit free-surface CG-solver controls."
      ps => cfg%si_precond
      call g%add(nml_enum("precond", ps, &
                          "SI Helmholtz preconditioner: none (plain CG), jacobi (diagonal PCG), "// &
                          "or mg (aggregation-multigrid V-cycle PCG)", &
                          allowed=[character(len=6) :: "none", "jacobi", "mg"]))
      pr => cfg%si_rtol
      call g%add(nml_real("rtol", pr, &
                          "SI CG relative residual tolerance (||r||^2 < rtol^2*||rhs||^2)", &
                          min=1.0e-30_wp))
      pl => cfg%si_wetdry
      call g%add(nml_logical("wetdry", pl, &
                             "Casulli subgrid wetting/drying (V(eta) Newton diagonal); "// &
                             "default off = byte-identical all-wet"))
      pi => cfg%si_wetdry_nsub
      call g%add(nml_int("wetdry_nsub", pi, &
                         "Subgrid sub-bottoms per cell for the wet/dry volume function "// &
                         "(1 = flat cell bottom)", min=1))
      ps => cfg%si_wetdry_bsub_fill
      call g%add(nml_enum("wetdry_bsub_fill", ps, &
                          "Sub-bottom table fill when nsub>=2: 'flat' (default, "// &
                          "bit-identical) or 'linear' (analytic subgrid bed ramp)", &
                          allowed=[character(len=6) :: "flat", "linear"]))
      pi => cfg%si_newton_max
      call g%add(nml_int("newton_max", pi, &
                         "Casulli wet/dry Newton outer-iteration cap (only consulted "// &
                         "when wetdry=.true.); default 10 = bit-identical", min=1))
      pr => cfg%si_newton_tol
      call g%add(nml_real("newton_tol", pr, &
                          "Casulli wet/dry Newton relative residual tolerance "// &
                          "(||F||^2 < tol^2*||rhs||^2); default 1e-11 is tight (tuned for "// &
                          "smooth cases) — real bathymetry typically wants ~1e-7", &
                          min=1.0e-30_wp))
      call schema%add_group(g)
   end subroutine register_si

   subroutine register_metadata(cfg, schema)
      !! `&metadata_nml`: CF-output global attributes.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      character(len=:), pointer :: ps

      g%name = "metadata"
      g%doc = "CF-output global attributes."
      ps => cfg%reference_time
      call g%add(nml_string("reference_time", ps, &
                            "Time origin for CF output (empty = idealized run)"))
      ps => cfg%institution
      call g%add(nml_string("institution", ps, &
                            "Institution name for the CF global attribute (empty = omitted)"))
      call schema%add_group(g)
   end subroutine register_metadata

   subroutine register_logging(cfg, schema)
      !! `&logging_nml`: logger verbosity + status cadence.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "logging"
      g%doc = "Logger verbosity + status cadence."
      ps => cfg%log_level
      call g%add(nml_enum("log_level", ps, "Log verbosity", &
                          allowed=[character(len=11) :: "debug", "verbose", "info", &
                                   "performance", "warning", "error"]))
      pr => cfg%status_interval
      call g%add(nml_real("status_interval", pr, &
                          "Status-print cadence (in time_unit; 0 = every 100 steps)"))
      call schema%add_group(g)
   end subroutine register_logging

   subroutine register_vcoord(cfg, schema)
      !! `&vcoord_nml`: vertical-coordinate + ALE-remap controls.
      !! `vcoord_type` / `remap_method` / `zstar_stretching` are enums
      !! with the canonical sets accepted by the rki_vcoord parsers.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      integer, pointer :: pi
      real(wp), pointer :: pr
      character(len=:), pointer :: ps
      logical, pointer :: pl

      g%name = "vcoord"
      g%doc = "Vertical-coordinate + ALE-remap controls."
      ps => cfg%vcoord_type
      call g%add(nml_enum("vcoord_type", ps, "Vertical coordinate type", &
                          allowed=[character(len=12) :: "sigma", "zsigma", &
                                   "zstar", "zstar_full", "zstar_sigma", "z", &
                                   "eulerian_z", "isopycnal", "lagrangian", "gprime", &
                                   "z_fixed", "rho", "hycom"]))
      ps => cfg%remap_method
      call g%add(nml_enum("remap_method", ps, "Vertical remapping method", &
                          allowed=[character(len=6) :: "pcm", "plm", "ppm", "ppm_h4", "pqm"]))
      pr => cfg%vcoord_depth_transition
      call g%add(nml_real("vcoord_depth_transition", pr, &
                          "z-sigma hybrid: depth where blending starts", units="m"))
      pr => cfg%vcoord_blend_width
      call g%add(nml_real("vcoord_blend_width", pr, &
                          "z-sigma hybrid: blending-zone width", units="m"))
      pr => cfg%zstar_h_surf_target
      call g%add(nml_real("zstar_h_surf_target", pr, &
                          "z*-full: target surface-layer thickness (0 = auto)", units="m"))
      ps => cfg%zstar_stretching
      call g%add(nml_enum("zstar_stretching", ps, &
                          "z*-full: surface-concentration stretching", &
                          allowed=[character(len=7) :: "log", "uniform"]))
      pr => cfg%zstar_h_min
      call g%add(nml_real("zstar_h_min", pr, "z*-full: vanishing-layer floor", units="m"))
      pi => cfg%zstar_n_surf
      call g%add(nml_int("zstar_n_surf", pi, &
                         "z*-full: number of fine near-surface layers (0 = auto)", min=0))
      pr => cfg%rho_ref_pressure
      call g%add(nml_real("rho_ref_pressure", pr, &
                          "rho-coord: reference pressure for potential density", units="Pa"))
      pr => cfg%rho_target_light
      call g%add(nml_real("rho_target_light", pr, &
                          "rho-coord: lightest (surface) target density", units="kg/m^3"))
      pr => cfg%rho_target_dense
      call g%add(nml_real("rho_target_dense", pr, &
                          "rho-coord: densest (bed) target density", units="kg/m^3"))
      call schema%add_group(g)
   end subroutine register_vcoord

   subroutine register_physics(cfg, schema)
      !! `&physics_nml`: barotropic physics + bottom drag + wind stress.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      logical, pointer :: pl

      g%name = "physics"
      g%doc = "Barotropic physics: bottom drag, wind stress, Coriolis."
      pr => cfg%manning_n
      call g%add(nml_real("manning_n", pr, "Manning roughness coefficient"))
      pr => cfg%dry_threshold
      call g%add(nml_real("dry_threshold", pr, "Depth below which a cell is dry", units="m"))
      pr => cfg%wind_stress_x
      call g%add(nml_real("wind_stress_x", pr, "Surface wind stress in x", units="Pa"))
      pr => cfg%wind_stress_y
      call g%add(nml_real("wind_stress_y", pr, "Surface wind stress in y", units="Pa"))
      pr => cfg%coriolis_f
      call g%add(nml_real("coriolis_f", pr, "Coriolis parameter f", units="1/s"))
      pl => cfg%use_hllc
      call g%add(nml_logical("use_hllc", pl, "Use HLLC contact-wave Riemann solver instead of HLL"))
      pl => cfg%use_linear_bottom_drag
      call g%add(nml_logical("use_linear_bottom_drag", pl, &
                             "Linear (Rayleigh) bottom drag instead of quadratic Manning"))
      pr => cfg%bottom_drag_r
      call g%add(nml_real("bottom_drag_r", pr, "Linear bottom-drag rate", units="1/s"))
      pl => cfg%wind_surface_layer
      call g%add(nml_logical("wind_surface_layer", pl, &
                             "Apply wind stress to the top layer only (multilayer)"))
      call schema%add_group(g)
   end subroutine register_physics

   subroutine register_output(cfg, schema)
      !! `&output_nml`: file output + I/O + forcing/restart/gauge files.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      integer, pointer :: pi
      logical, pointer :: pl
      character(len=:), pointer :: ps

      g%name = "output"
      g%doc = "File output, forcing/restart/gauge files."
      pl => cfg%output_to_file
      call g%add(nml_logical("output_to_file", pl, "Enable file output (NetCDF snapshots)"))
      ps => cfg%output_dir
      call g%add(nml_string("output_dir", ps, "Directory for output files"))
      pr => cfg%output_interval
      call g%add(nml_real("output_interval", pr, "Time between output snapshots", units="s"))
      ps => cfg%output_format
      call g%add(nml_enum("output_format", ps, "Output format", &
                          allowed=[character(len=6) :: "ascii", "netcdf"]))
      pr => cfg%restart_interval
      call g%add(nml_real("restart_interval", pr, "Time between restart writes (0 = none)", units="s"))
      pl => cfg%async_output
      call g%add(nml_logical("async_output", pl, "Async ring-buffered output (CUDA Fortran)"))
      pi => cfg%async_buffer_slots
      call g%add(nml_int("async_buffer_slots", pi, "Async ring-buffer slot count", min=1))
      pl => cfg%compress_output
      call g%add(nml_logical("compress_output", pl, "Deflate-compress NetCDF output"))
      pi => cfg%compress_level
      call g%add(nml_int("compress_level", pi, "Deflate level (1=fast, 9=max)", min=1, max=9))
      ps => cfg%bathymetry_file
      call g%add(nml_string("bathymetry_file", ps, "NetCDF bathymetry file (empty = flat)"))
      ps => cfg%restart_file
      call g%add(nml_string("restart_file", ps, "Restart file for warm start (empty = cold)"))
      ps => cfg%forcing_file
      call g%add(nml_string("forcing_file", ps, "NetCDF atmospheric forcing file (empty = uniform)"))
      pr => cfg%forcing_x_origin
      call g%add(nml_real("forcing_x_origin", pr, "Model x-origin in forcing-file coordinates"))
      pr => cfg%forcing_y_origin
      call g%add(nml_real("forcing_y_origin", pr, "Model y-origin in forcing-file coordinates"))
      ps => cfg%gauge_file
      call g%add(nml_string("gauge_file", ps, "Gauge-locations text file (empty = none)"))
      pr => cfg%gauge_output_interval
      call g%add(nml_real("gauge_output_interval", pr, "Gauge output interval (0 = every step)", units="s"))
      pl => cfg%out_h_layer
      call g%add(nml_logical("out_h_layer", pl, "Output layer thickness"))
      pl => cfg%out_velocity
      call g%add(nml_logical("out_velocity", pl, "Output layer velocity (u_ml, v_ml)"))
      pl => cfg%out_salinity
      call g%add(nml_logical("out_salinity", pl, "Output layer salinity"))
      pl => cfg%out_temperature
      call g%add(nml_logical("out_temperature", pl, "Output layer temperature"))
      pl => cfg%out_density
      call g%add(nml_logical("out_density", pl, "Output layer density"))
      pl => cfg%out_w
      call g%add(nml_logical("out_w", pl, "Output vertical velocity at interfaces"))
      call schema%add_group(g)
   end subroutine register_output

   subroutine register_boundary(cfg, schema)
      !! `&boundary_nml`: per-side BC types, tidal forcing, inflow/
      !! discharge/clamped/sponge/nesting parameters.  The four `bc_*`
      !! keys are enums over the bc_type_from_string accepted set.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      real(wp), pointer :: pa(:)
      integer, pointer :: pi
      character(len=:), pointer :: ps
      character(len=9), parameter :: bc_allowed(*) = &
                                     [character(len=9) :: "wall", "open", "tidal", "nested", "inflow", &
                                                           "discharge", "clamped", "sponge", "chapman"]

      g%name = "boundary"
      g%doc = "Boundary conditions, tidal forcing, inflow/discharge/sponge/nesting."
      ps => cfg%bc_west
      call g%add(nml_enum("bc_west", ps, "West boundary type", allowed=bc_allowed))
      ps => cfg%bc_east
      call g%add(nml_enum("bc_east", ps, "East boundary type", allowed=bc_allowed))
      ps => cfg%bc_south
      call g%add(nml_enum("bc_south", ps, "South boundary type", allowed=bc_allowed))
      ps => cfg%bc_north
      call g%add(nml_enum("bc_north", ps, "North boundary type", allowed=bc_allowed))
      pr => cfg%tidal_amplitude
      call g%add(nml_real("tidal_amplitude", pr, "Tidal elevation amplitude (single-constituent)", units="m"))
      pr => cfg%tidal_period
      call g%add(nml_real("tidal_period", pr, "Tidal period (single-constituent, default M2)", units="s"))
      pi => cfg%n_tidal_constituents
      call g%add(nml_int("n_tidal_constituents", pi, &
                         "Active tidal constituents (0 = legacy single)", min=0))
      pa => cfg%tidal_amp
      call g%add(nml_real_array("tidal_amp", pa, "Constituent amplitudes", units="m"))
      pa => cfg%tidal_phase
      call g%add(nml_real_array("tidal_phase", pa, "Constituent phases", units="rad"))
      pa => cfg%tidal_omega
      call g%add(nml_real_array("tidal_omega", pa, "Constituent angular frequencies", units="rad/s"))
      pr => cfg%inflow_velocity
      call g%add(nml_real("inflow_velocity", pr, "Prescribed normal velocity for inflow BC", units="m/s"))
      pr => cfg%inflow_salinity
      call g%add(nml_real("inflow_salinity", pr, "Inflow salinity (<0 = zero-gradient)", units="PSU"))
      pr => cfg%inflow_temperature
      call g%add(nml_real("inflow_temperature", pr, "Inflow temperature (<0 = zero-gradient)", units="degC"))
      pr => cfg%discharge_flux
      call g%add(nml_real("discharge_flux", pr, "Prescribed volume flux for discharge BC", units="m^2/s"))
      pr => cfg%clamped_elevation
      call g%add(nml_real("clamped_elevation", pr, "Prescribed surface elevation for clamped BC", units="m"))
      pi => cfg%sponge_width
      call g%add(nml_int("sponge_width", pi, "Sponge layer width in cells", min=0))
      pr => cfg%sponge_strength
      call g%add(nml_real("sponge_strength", pr, "Sponge relaxation rate", units="1/s"))
      pr => cfg%sponge_h_ref
      call g%add(nml_real("sponge_h_ref", pr, "Sponge reference depth", units="m"))
      ps => cfg%nesting_file
      call g%add(nml_string("nesting_file", ps, "Parent-domain NetCDF output (empty = no nesting)"))
      pr => cfg%nesting_x_origin
      call g%add(nml_real("nesting_x_origin", pr, "Child x-origin in parent coordinates", units="m"))
      pr => cfg%nesting_y_origin
      call g%add(nml_real("nesting_y_origin", pr, "Child y-origin in parent coordinates", units="m"))
      call schema%add_group(g)
   end subroutine register_boundary

   subroutine register_tracer(cfg, schema)
      !! `&tracer_nml`: salinity + temperature IC/EOS/bounds + sediment.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      logical, pointer :: pl
      character(len=:), pointer :: ps

      g%name = "tracer"
      g%doc = "Salinity + temperature IC/EOS/bounds, sediment transport."
      pr => cfg%initial_salinity
      call g%add(nml_real("initial_salinity", pr, "Initial salinity (uniform IC)", units="PSU"))
      pr => cfg%S_ref
      call g%add(nml_real("S_ref", pr, "EOS reference salinity", units="PSU"))
      pr => cfg%beta_S
      call g%add(nml_real("beta_S", pr, "Haline contraction coefficient", units="kg/m^3/PSU"))
      pr => cfg%S_min
      call g%add(nml_real("S_min", pr, "Lower physical bound for salinity", units="PSU"))
      pr => cfg%S_max
      call g%add(nml_real("S_max", pr, "Upper physical bound for salinity", units="PSU"))
      pr => cfg%kappa_S_bg
      call g%add(nml_real("kappa_S_bg", pr, "Background vertical salinity diffusivity", units="m^2/s"))
      pr => cfg%S_init_surface
      call g%add(nml_real("S_init_surface", pr, "Initial surface salinity (stratified IC)", units="PSU"))
      pr => cfg%S_init_bottom
      call g%add(nml_real("S_init_bottom", pr, "Initial bed salinity (stratified IC)", units="PSU"))
      pr => cfg%initial_temperature
      call g%add(nml_real("initial_temperature", pr, "Initial temperature (uniform IC)", units="degC"))
      pr => cfg%T_ref
      call g%add(nml_real("T_ref", pr, "EOS reference temperature", units="degC"))
      pr => cfg%alpha_T
      call g%add(nml_real("alpha_T", pr, "Thermal expansion coefficient", units="kg/m^3/degC"))
      pr => cfg%T_min
      call g%add(nml_real("T_min", pr, "Lower physical bound for temperature", units="degC"))
      pr => cfg%T_max
      call g%add(nml_real("T_max", pr, "Upper physical bound for temperature", units="degC"))
      pr => cfg%kappa_T_bg
      call g%add(nml_real("kappa_T_bg", pr, "Background vertical temperature diffusivity", units="m^2/s"))
      ps => cfg%eos_form
      call g%add(nml_enum("eos_form", ps, &
                          "Coastal EOS: linear (legacy, default), wright or roquet_spv (shared nonlinear)", &
                          allowed=[character(len=10) :: "linear", "wright", "roquet_spv"]))
      pr => cfg%T_init_surface
      call g%add(nml_real("T_init_surface", pr, "Initial surface temperature (stratified IC)", units="degC"))
      pr => cfg%T_init_bottom
      call g%add(nml_real("T_init_bottom", pr, "Initial bed temperature (stratified IC)", units="degC"))
      pl => cfg%use_ppm_tracer
      call g%add(nml_logical("use_ppm_tracer", pl, "High-order tracer advection (PPM/MUSCL) vs upwind"))
      ps => cfg%tracer_recon
      call g%add(nml_enum("tracer_recon", ps, &
                          "Tracer face reconstruction (structured coastal): ppm|weno5|weno7|weno9", &
                          allowed=[character(len=8) :: "ppm", "pqm", "weno5", "weno7", "weno9"]))
      pl => cfg%use_sediment
      call g%add(nml_logical("use_sediment", pl, "Master switch for suspended-sediment transport"))
      pr => cfg%sediment_w_s
      call g%add(nml_real("sediment_w_s", pr, "Settling velocity (downward; 0 = none)", units="m/s"))
      pr => cfg%sediment_initial_conc
      call g%add(nml_real("sediment_initial_conc", pr, "Uniform initial suspended concentration", units="kg/m^3"))
      pr => cfg%sediment_bed_init
      call g%add(nml_real("sediment_bed_init", pr, "Uniform initial bed sediment mass", units="kg/m^2"))
      pr => cfg%sediment_tau_crit_erosion
      call g%add(nml_real("sediment_tau_crit_erosion", pr, "Critical erosion bed shear stress (0 = off)", units="Pa"))
      pr => cfg%sediment_erosion_rate
      call g%add(nml_real("sediment_erosion_rate", pr, "Ariathurai erosion-rate coefficient E0", units="kg/m^2/s"))
      pr => cfg%sediment_porosity
      call g%add(nml_real("sediment_porosity", pr, "Bed porosity (0-1); erosion scaled by (1-porosity)"))
      call schema%add_group(g)
   end subroutine register_tracer

   subroutine register_initial_condition(cfg, schema)
      !! `&initial_condition_nml`: coastal IC selector + parameters.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "initial_condition"
      g%doc = "Coastal initial-condition selector + parameters."
      ps => cfg%init_condition
      call g%add(nml_enum("init_condition", ps, "Initial condition", &
                          allowed=[character(len=16) :: "dam_break", "gaussian_hump", &
                                   "castro_dam_break", "still_water"]))
      pr => cfg%eta_initial
      call g%add(nml_real("eta_initial", pr, "Initial free-surface elevation (still_water IC)", units="m"))
      pr => cfg%h0
      call g%add(nml_real("h0", pr, "Background depth (gaussian_hump IC)", units="m"))
      pr => cfg%hump_amplitude
      call g%add(nml_real("hump_amplitude", pr, "Perturbation amplitude (gaussian_hump IC)", units="m"))
      pr => cfg%hump_sigma
      call g%add(nml_real("hump_sigma", pr, "Gaussian width (gaussian_hump IC)", units="m"))
      pr => cfg%hump_x_offset
      call g%add(nml_real("hump_x_offset", pr, "Hump centre x-offset from domain centre", units="m"))
      pr => cfg%hump_y_offset
      call g%add(nml_real("hump_y_offset", pr, "Hump centre y-offset from domain centre", units="m"))
      pr => cfg%init_u_shear
      call g%add(nml_real("init_u_shear", pr, "Initial vertical shear (top +val, bottom -val)", units="m/s"))
      call schema%add_group(g)
   end subroutine register_initial_condition

   subroutine register_nonhydrostatic(cfg, schema)
      !! `&nonhydrostatic_nml`: NH/multilayer switches, CG-Poisson
      !! controls, PP81 vmix, KPP knobs, k-eps, Smagorinsky, BPG, mode
      !! split.  `keps_stability` and `bpg_method` are enums.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      integer, pointer :: pi
      logical, pointer :: pl
      character(len=:), pointer :: ps

      g%name = "nonhydrostatic"
      g%doc = "Non-hydrostatic / multilayer: CG-Poisson, PP81/KPP/k-eps vmix, Smagorinsky, BPG, mode split."
      pl => cfg%use_nonhydrostatic
      call g%add(nml_logical("use_nonhydrostatic", pl, "Enable multi-layer NH pressure correction"))
      pi => cfg%nz_layers
      call g%add(nml_int("nz_layers", pi, "Number of sigma layers", min=1))
      pr => cfg%nh_tolerance
      call g%add(nml_real("nh_tolerance", pr, "CG convergence tolerance for NH pressure"))
      pi => cfg%nh_max_iterations
      call g%add(nml_int("nh_max_iterations", pi, "Maximum CG iterations for NH pressure", min=1))
      pl => cfg%nh_preconditioner
      call g%add(nml_logical("nh_preconditioner", pl, "Vertical-line preconditioner for CG solver"))
      pi => cfg%nh_solve_interval
      call g%add(nml_int("nh_solve_interval", pi, "Solve NH pressure every N steps", min=1))
      pi => cfg%split_ratio
      call g%add(nml_int("split_ratio", pi, "Mode-split ratio: 0=auto, 1=unsplit, N>1=fixed", min=0))
      pl => cfg%split_averaging
      call g%add(nml_logical("split_averaging", pl, "ROMS-style averaging over barotropic substeps"))
      pl => cfg%use_multilayer
      call g%add(nml_logical("use_multilayer", pl, "Enable coupled hydrostatic vertical layers"))
      pr => cfg%rho_0
      call g%add(nml_real("rho_0", pr, "Reference density for EOS", units="kg/m^3"))
      pr => cfg%vmix_nu_max
      call g%add(nml_real("vmix_nu_max", pr, "PP81 max turbulent viscosity", units="m^2/s"))
      pr => cfg%vmix_nu_bg
      call g%add(nml_real("vmix_nu_bg", pr, "Background vertical viscosity", units="m^2/s"))
      pr => cfg%vmix_ri_alpha
      call g%add(nml_real("vmix_ri_alpha", pr, "PP81 Richardson-number scaling coefficient"))
      pl => cfg%vmix_use_kpp
      call g%add(nml_logical("vmix_use_kpp", pl, "Enable KPP surface-boundary-layer overlay on PP81"))
      pl => cfg%vmix_use_keps
      call g%add(nml_logical("vmix_use_keps", pl, "Enable two-equation k-epsilon vertical mixing"))
      pi => cfg%keps_si_picard
      call g%add(nml_int("keps_si_picard", pi, &
                         "Picard passes of the k<->eps coupling per semi-implicit step (SI+keps only)", &
                         min=1, max=3))
      pr => cfg%kpp_ri_crit
      call g%add(nml_real("kpp_ri_crit", pr, "Critical bulk Richardson number for KPP BL-depth"))
      pr => cfg%kpp_cs_nonlocal
      call g%add(nml_real("kpp_cs_nonlocal", pr, "KPP non-local (counter-gradient) transport coefficient"))
      pr => cfg%kpp_c_vt2
      call g%add(nml_real("kpp_c_vt2", pr, "KPP V_t^2 unresolved-turbulence coefficient (0 = off)"))
      pr => cfg%hdiff_kappa
      call g%add(nml_real("hdiff_kappa", pr, &
                          "Horizontal tracer diffusion coefficient", &
                          units="m^2/s"))
      pr => cfg%smagorinsky_cs
      call g%add(nml_real("smagorinsky_cs", pr, "Smagorinsky horizontal-viscosity coefficient"))
      pr => cfg%biharmonic_nu4
      call g%add(nml_real("biharmonic_nu4", pr, "Biharmonic lateral viscosity (0 = off)", units="m^4/s"))
      pl => cfg%evap_mass_sink
      call g%add(nml_logical("evap_mass_sink", pl, &
                             "Treat net evaporation as a real surface mass sink (default off = virtual salt flux)"))
      pl => cfg%evap_latent_cooling
      call g%add(nml_logical("evap_latent_cooling", pl, &
                             "Evaporative latent-heat surface cooling (requires evap_mass_sink; default off)"))
      pr => cfg%ml_min_depth
      call g%add(nml_real("ml_min_depth", pr, "Minimum total depth for vertical-layer operations", units="m"))
      ps => cfg%bpg_method
      call g%add(nml_enum("bpg_method", ps, "BPG computation method", &
                          allowed=[character(len=8) :: "sigma", "jacobian"]))
      ps => cfg%keps_stability
      call g%add(nml_enum("keps_stability", ps, "k-epsilon stability functions", &
                          allowed=[character(len=9) :: "constant", "galperin", "canuto"]))
      call schema%add_group(g)
   end subroutine register_nonhydrostatic

end module rki_config
