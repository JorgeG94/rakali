!! Solution state type for the shallow water solver
module rki_state
   !! Manages the conserved variables, bathymetry, and flux arrays
   !!
   !! All arrays use SoA layout with dimensions (nx_total, ny_total) where
   !! nx_total = nx_physical + 2*nghost. The first index (i / x-direction)
   !! is the fast-varying index for coalesced GPU memory access.
   use, intrinsic :: iso_fortran_env, only: int64
   use rki_constants, only: wp, KEPS_STAB_CONSTANT, KEPS_STAB_GALPERIN, KEPS_STAB_CANUTO
   use rki_config, only: config_t
   use rki_mem_report, only: arr_bytes
   use rki_boundary_types, only: bc_params_t, BC_WALL
   use rki_vcoord, only: parse_vcoord_type, parse_remap_method, parse_stretching_mode
   use rki_timestepper, only: parse_time_scheme, TIME_RK2, TIME_SEMI_IMPLICIT
   use rki_grid, only: hgrid_t
   use rki_surface_forcing, only: surface_forcing_t
   use rki_barotropic_state, only: barotropic_state_t
   use rki_multilayer_state, only: multilayer_state_t
   use rki_nonhydrostatic_state, only: nonhydrostatic_state_t
   use rki_tracer, only: tracer_t
#ifndef RAKALI_NO_NETCDF
   use rki_bathymetry, only: load_bathymetry_into_array, &
                             fill_bathymetry_ghosts_array
#endif
   implicit none
   private

   public :: state_t
   public :: register_default_tracers
   public :: load_bathymetry, fill_bathymetry_ghosts

   type :: state_t
      !! Container for all solution arrays and grid metadata

      ! Barotropic conserved variables, bathymetry, flux scratch, RK2
      ! saves, and barotropic physics scalars (Manning's n, Coriolis f).
      type(barotropic_state_t) :: barotropic

      ! Grid metadata
      type(hgrid_t) :: grid
         !! Cell counts, ghost width, and cell sizes

      ! Time-stepping
      real(wp) :: dt = 0.0_wp
         !! Current timestep (s)
      real(wp) :: cfl = 0.45_wp
         !! CFL number
      real(wp) :: dt_fixed = 0.0_wp
         !! Fixed timestep (s), 0 = adaptive CFL
      integer :: cfl_interval = 1
         !! Recompute CFL timestep every N steps (1 = every step)
      integer :: time_scheme = TIME_RK2
         !! Explicit time integrator (`TIME_*` from rki_timestepper).
         !! Default RK2 reproduces the historical SSP-RK2 path bit-for-bit.
      integer :: step_count = 0
         !! Internal step counter for CFL interval

      ! Surface forcing fields (wind, pressure, precipitation, heat,
      ! virtual freshwater).  Filled by the forcing reader or API setters;
      ! consumed by barotropic friction + multilayer surface-flux kernels.
      type(surface_forcing_t) :: forcing

      ! Boundary condition parameters
      type(bc_params_t) :: bc
         !! Boundary condition configuration

      ! Coupled hydrostatic multi-layer
      logical :: use_multilayer = .false.
         !! Enable coupled hydrostatic vertical layers
      type(multilayer_state_t) :: multilayer
         !! Per-layer arrays, vcoord, tracer registry, EOS, vmix, BPG.

      ! Non-hydrostatic pressure correction (multi-layer) + mode-split
      logical :: use_nonhydrostatic = .false.
         !! Enable NH pressure correction
      type(nonhydrostatic_state_t) :: nh
         !! Per-layer NH velocities/pressure + RK2 saves, CG workspaces,
         !! barotropic pre-HLL momentum saves, mode-splitting arrays.

   contains
      procedure, non_overridable :: init => state_init
      procedure, non_overridable :: init_from_config => state_init_from_config
      procedure, non_overridable :: destroy => state_destroy
      procedure, non_overridable :: bytes => state_bytes
   end type state_t

contains

   subroutine state_init(this, nx_phys, ny_phys, nghost, dx, dy, cfl)
      !! Allocate all arrays and set grid metadata
      class(state_t), intent(inout) :: this
      integer, intent(in) :: nx_phys
         !! Number of physical cells in x
      integer, intent(in) :: ny_phys
         !! Number of physical cells in y
      integer, intent(in) :: nghost
         !! Ghost cell width
      real(wp), intent(in) :: dx
         !! Cell size in x (m)
      real(wp), intent(in) :: dy
         !! Cell size in y (m)
      real(wp), intent(in) :: cfl
         !! CFL number

      integer :: nx, ny

      this%cfl = cfl

      call this%grid%init(nx_phys, ny_phys, nghost, dx, dy)
      nx = this%grid%nx_total
      ny = this%grid%ny_total

      ! Barotropic state arrays (h, hu, hv + RK saves, bathymetry,
      ! flux + mass-flux scratch)
      call this%barotropic%init(this%grid)

      ! Surface forcing arrays (initialised to zero = no forcing)
      call this%forcing%init(this%grid)

      ! Coupled hydrostatic multi-layer arrays + tracer registry.
      ! Caller must set nz_ml + vcoord params first; multilayer%init reads them.
      if (this%use_multilayer) then
         call this%multilayer%init(this%grid)
      end if

      ! Non-hydrostatic arrays + mode-splitting buffers.  Caller must set
      ! nz, CG/NH config, and split_ratio first; nh%init reads them.
      if (this%use_nonhydrostatic) then
         call this%nh%init(this%grid)
      end if

   end subroutine state_init

   subroutine state_init_from_config(this, cfg)
      !! Allocate arrays from config
      class(state_t), intent(inout) :: this
      type(config_t), intent(in) :: cfg
         !! Configuration parameters

      ! Set multilayer parameters from config (read by multilayer%init below)
      this%use_multilayer = cfg%use_multilayer
      if (cfg%use_multilayer) then
         this%multilayer%nz_ml = cfg%nz_layers
         ! Vcoord parameters so multilayer%init picks them up
         this%multilayer%vcoord%coord_type = parse_vcoord_type(cfg%vcoord_type)
         this%multilayer%vcoord%remap_method = parse_remap_method(cfg%remap_method)
         this%multilayer%vcoord%depth_transition = cfg%vcoord_depth_transition
         this%multilayer%vcoord%blend_width = cfg%vcoord_blend_width
         this%multilayer%vcoord%zstar_h_surf_target = cfg%zstar_h_surf_target
         this%multilayer%vcoord%zstar_h_min = cfg%zstar_h_min
         this%multilayer%vcoord%zstar_stretching = parse_stretching_mode(cfg%zstar_stretching)
         this%multilayer%vcoord%zstar_n_surf = cfg%zstar_n_surf
         this%multilayer%initial_salinity = cfg%initial_salinity
         this%multilayer%S_init_surface = cfg%S_init_surface
         this%multilayer%S_init_bottom = cfg%S_init_bottom
         this%multilayer%S_ref = cfg%S_ref
         this%multilayer%beta_S = cfg%beta_S
         this%multilayer%S_min = cfg%S_min
         this%multilayer%S_max = cfg%S_max
         this%multilayer%kappa_S_bg = cfg%kappa_S_bg
         this%multilayer%initial_temperature = cfg%initial_temperature
         this%multilayer%T_init_surface = cfg%T_init_surface
         this%multilayer%T_init_bottom = cfg%T_init_bottom
         this%multilayer%T_ref = cfg%T_ref
         this%multilayer%alpha_T = cfg%alpha_T
         this%multilayer%T_min = cfg%T_min
         this%multilayer%T_max = cfg%T_max
         this%multilayer%kappa_T_bg = cfg%kappa_T_bg
         this%multilayer%rho_0 = cfg%rho_0
         ! Coastal EOS selector (&tracer_nml eos_form): default "linear"
         ! keeps the legacy kernel (eos_use_shared = .false.); "wright" /
         ! "roquet_spv" opt in to the shared nonlinear EOS.  The handle
         ! is a flat POD read on the host at dispatch — no device map.
         block
            use rki_ml_eos, only: ml_parse_eos_form
            call ml_parse_eos_form(cfg%eos_form, &
                                   this%multilayer%eos_use_shared, &
                                   this%multilayer%eos%variant)
         end block
         this%multilayer%eos%rho0 = cfg%rho_0
         this%multilayer%eos%S_ref = cfg%S_ref
         this%multilayer%eos%T_ref = cfg%T_ref
         this%multilayer%eos%alpha_T = cfg%alpha_T
         this%multilayer%eos%beta_S = cfg%beta_S
         this%multilayer%eos%is_init = .true.
         this%multilayer%vmix_nu_max = cfg%vmix_nu_max
         this%multilayer%vmix_nu_bg = cfg%vmix_nu_bg
         this%multilayer%vmix_ri_alpha = cfg%vmix_ri_alpha
         this%multilayer%vmix_use_kpp = cfg%vmix_use_kpp
         this%multilayer%vmix_use_keps = cfg%vmix_use_keps
         this%multilayer%kpp_ri_crit = cfg%kpp_ri_crit
         this%multilayer%kpp_cs_nonlocal = cfg%kpp_cs_nonlocal
         this%multilayer%kpp_c_vt2 = cfg%kpp_c_vt2
         this%multilayer%hdiff_kappa = cfg%hdiff_kappa
         this%multilayer%smagorinsky_cs = cfg%smagorinsky_cs
         this%multilayer%biharmonic_nu4 = cfg%biharmonic_nu4
         this%multilayer%evap_mass_sink = cfg%evap_mass_sink
         this%multilayer%evap_latent_cooling = cfg%evap_latent_cooling
         this%multilayer%use_bpg_jacobian = (trim(cfg%bpg_method) == "jacobian")
         this%multilayer%keps_stability = KEPS_STAB_CONSTANT
         if (trim(cfg%keps_stability) == "galperin") then
            this%multilayer%keps_stability = KEPS_STAB_GALPERIN
         else if (trim(cfg%keps_stability) == "canuto") then
            this%multilayer%keps_stability = KEPS_STAB_CANUTO
         end if
         this%multilayer%use_ppm_tracer = cfg%use_ppm_tracer
         block
            use rki_recon_weno, only: parse_tracer_recon
            this%multilayer%tracer_recon_code = parse_tracer_recon(trim(cfg%tracer_recon))
         end block
         ! Sediment (set before this%init so bed_sediment allocates).
         this%multilayer%use_sediment = cfg%use_sediment
         this%multilayer%sediment_w_s = cfg%sediment_w_s
         this%multilayer%sediment_initial_conc = cfg%sediment_initial_conc
         this%multilayer%sediment_bed_init = cfg%sediment_bed_init
         this%multilayer%sediment_tau_crit_erosion = cfg%sediment_tau_crit_erosion
         this%multilayer%sediment_erosion_rate = cfg%sediment_erosion_rate
         this%multilayer%sediment_porosity = cfg%sediment_porosity
      end if

      ! Parse the time-scheme tag up front.  Semi-implicit is unstructured-only;
      ! a direct-API caller who sets only the string must fail loud here.
      this%time_scheme = parse_time_scheme(cfg%time_scheme)
      if (this%time_scheme == TIME_SEMI_IMPLICIT) then
         error stop "time_scheme='semi_implicit' is only supported on "// &
            "unstructured grids"
      end if

      ! Set NH and mode-splitting parameters from config (read by
      ! nh%init below)
      this%use_nonhydrostatic = cfg%use_nonhydrostatic
      this%nh%nz = cfg%nz_layers
      this%nh%nh_tolerance = cfg%nh_tolerance
      this%nh%nh_max_iterations = cfg%nh_max_iterations
      this%nh%nh_preconditioner = cfg%nh_preconditioner
      this%nh%nh_solve_interval = cfg%nh_solve_interval
      this%nh%split_ratio = cfg%split_ratio
      this%nh%split_averaging = cfg%split_averaging

      ! Delegate to main init (allocates arrays)
      call this%init(cfg%nx, cfg%ny, cfg%nghost, cfg%dx, cfg%dy, cfg%cfl)

      ! Set physics parameters from config
      this%barotropic%manning_n = cfg%manning_n
      this%barotropic%use_linear_bottom_drag = cfg%use_linear_bottom_drag
      this%barotropic%bottom_drag_r = cfg%bottom_drag_r
      this%barotropic%wind_surface_layer = cfg%wind_surface_layer
      this%forcing%wind_stress_x = cfg%wind_stress_x
      this%forcing%wind_stress_y = cfg%wind_stress_y
      this%barotropic%coriolis_f = cfg%coriolis_f
      this%barotropic%use_hllc = cfg%use_hllc
      this%dt_fixed = cfg%dt_fixed
      this%cfl_interval = cfg%cfl_interval

      ! The selectable explicit integrator (euler/rk3/rk43) is wired into
      ! the 2D barotropic fast path only.  Multilayer + non-hydrostatic
      ! keep their hard-wired SSP-RK2 (mode-split / pressure stages sit
      ! between the Euler update and the combine, so they can't share the
      ! generic N-stage loop yet).  Fail loud rather than silently ignore.
      ! The selectable explicit integrators (euler/rk3/rk43) are wired into
      ! the 2D barotropic fast path only.  NH is excluded on all non-rk2
      ! schemes.
      if (this%time_scheme /= TIME_RK2 .and. this%use_nonhydrostatic) then
         error stop "time_scheme: euler/rk3/rk43 are incompatible "// &
            "with use_nonhydrostatic=.true.; non-hydrostatic runs use rk2."
      end if
      if (this%time_scheme /= TIME_RK2 .and. this%use_multilayer) then
         error stop "time_scheme: euler/rk3/rk43 require the 2D barotropic path "// &
            "(use_multilayer=.false.); multilayer runs use rk2."
      end if

      ! Fill spatially-varying forcing arrays from uniform config values
      this%forcing%wind_x(:, :) = cfg%wind_stress_x
      this%forcing%wind_y(:, :) = cfg%wind_stress_y

      ! Register the prognostic salinity + temperature tracers in the
      ! generic registry, accessed via tracers(idx_salinity/idx_temperature)%hTr.
      if (cfg%use_multilayer) then
         call register_default_tracers( &
            this%multilayer%tracers(this%multilayer%idx_salinity), &
            this%multilayer%tracers(this%multilayer%idx_temperature), &
            cfg)
         ! Register the suspended-sediment tracer (single class) when enabled.
         ! Horizontal transport / diffusion / clamp ride the registry; settling
         ! + bed exchange are separate kernels.  eos_coeff = 0 (no density
         ! feedback); concentration in kg/m^3, clamped >= 0 (no upper bound).
         if (cfg%use_sediment) then
            call this%multilayer%register_passive_tracer(this%grid, "sediment", &
                                                         cfg%sediment_initial_conc, &
                                                         0.0_wp, huge(1.0_wp), &
                                                         cfg%kappa_S_bg, cfg%hdiff_kappa)
            this%multilayer%idx_sediment = size(this%multilayer%tracers)
            associate (tr_sed => this%multilayer%tracers(this%multilayer%idx_sediment))
               tr_sed%long_name = "Suspended sediment concentration"
               tr_sed%units = "kg m-3"
               tr_sed%standard_name = "mass_concentration_of_suspended_matter_in_sea_water"
            end associate
         end if
      end if

   end subroutine state_init_from_config

   subroutine register_default_tracers(tr_S, tr_T, cfg)
      !! Overlay scalar config onto the already-allocated salinity /
      !! temperature tracer slots (by reference so coastal + ocean paths
      !! share the overlay).  Does not touch hTr / hTr0 — those are
      !! populated by the path's IC step once h_layer is known.
      type(tracer_t), intent(inout) :: tr_S
      type(tracer_t), intent(inout) :: tr_T
      type(config_t), intent(in)    :: cfg

      ! Salinity
      tr_S%name = "salinity"
      tr_S%long_name = "Sea water salinity"
      tr_S%units = "PSU"
      tr_S%standard_name = "sea_water_salinity"
      tr_S%tr_init = cfg%initial_salinity
      tr_S%tr_inflow = cfg%inflow_salinity
      tr_S%tr_min = cfg%S_min
      tr_S%tr_max = cfg%S_max
      tr_S%kappa_bg = cfg%kappa_S_bg
      tr_S%hdiff_kappa = cfg%hdiff_kappa
      tr_S%eos_coeff = cfg%beta_S
      tr_S%eos_ref = cfg%S_ref

      ! Temperature
      tr_T%name = "temperature"
      tr_T%long_name = "Sea water potential temperature"
      tr_T%units = "degC"
      tr_T%standard_name = "sea_water_potential_temperature"
      tr_T%tr_init = cfg%initial_temperature
      tr_T%tr_inflow = cfg%inflow_temperature
      tr_T%tr_min = cfg%T_min
      tr_T%tr_max = cfg%T_max
      tr_T%kappa_bg = cfg%kappa_T_bg
      tr_T%hdiff_kappa = cfg%hdiff_kappa
      tr_T%eos_coeff = -cfg%alpha_T
      tr_T%eos_ref = cfg%T_ref

   end subroutine register_default_tracers

   subroutine state_destroy(this)
      !! Deallocate all arrays
      class(state_t), intent(inout) :: this

      call this%barotropic%destroy()
      call this%forcing%destroy()
      call this%multilayer%destroy()
      call this%nh%destroy()

      ! Chapman BC arrays
      if (allocated(this%bc%eta_old_west)) deallocate (this%bc%eta_old_west)
      if (allocated(this%bc%eta_old_east)) deallocate (this%bc%eta_old_east)
      if (allocated(this%bc%eta_old_south)) deallocate (this%bc%eta_old_south)
      if (allocated(this%bc%eta_old_north)) deallocate (this%bc%eta_old_north)

   end subroutine state_destroy

   pure function state_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the whole coastal state, summed
      !! from each sub-slot's own `bytes()` (each 0 when its feature is
      !! gated off, so the total tracks conditional allocation directly).
      !! Reported before `enter_data` and reconciled against the measured
      !! device mapping — see `rki_mem_report`.
      class(state_t), intent(in) :: this
      integer(int64) :: nbytes

      nbytes = this%barotropic%bytes() + this%forcing%bytes()
      if (this%use_multilayer) nbytes = nbytes + this%multilayer%bytes()
      if (this%use_nonhydrostatic) nbytes = nbytes + this%nh%bytes()
      ! Chapman persistent-state arrays live directly on the bc slot.
      nbytes = nbytes + arr_bytes(this%bc%eta_old_west) &
               + arr_bytes(this%bc%eta_old_east) &
               + arr_bytes(this%bc%eta_old_south) &
               + arr_bytes(this%bc%eta_old_north)
   end function state_bytes

   subroutine load_bathymetry(filename, state)
      !! Coastal-state entry point — forwards to the array-based loader
      !! in `rki_bathymetry`.  Lives here (not there) to avoid a USE cycle.
      character(len=*), intent(in) :: filename
      type(state_t), intent(inout) :: state
#ifndef RAKALI_NO_NETCDF
      call load_bathymetry_into_array(filename, state%barotropic%b, state%grid)
#else
      ! File-loaded bathymetry needs the NetCDF reader, compiled out when
      ! RAKALI_ENABLE_NETCDF=OFF.  Fail fast; embed `filename` (also avoids
      ! an unused-arg warning).
      if (state%grid%nghost < 0) return     ! refs `state` to silence warn
      error stop "load_bathymetry: requires RAKALI_ENABLE_NETCDF=ON at "// &
         "build time (the NetCDF-backed rki_bathymetry reader is needed "// &
         "to load '"//trim(filename)//"')."
#endif
   end subroutine load_bathymetry

   subroutine fill_bathymetry_ghosts(state)
      !! Wrapper mirror of `load_bathymetry` — populates ghost rows on
      !! `state%barotropic%b` after a custom bathy load (e.g. via the
      !! Python FFI).
      type(state_t), intent(inout) :: state
#ifndef RAKALI_NO_NETCDF
      call fill_bathymetry_ghosts_array(state%barotropic%b, state%grid)
#else
      ! Same gate — the ghost-fill helper also lives in rki_bathymetry,
      ! which is compiled out when RAKALI_ENABLE_NETCDF=OFF.
      if (state%grid%nghost < 0) return     ! refs `state` to silence warn
      error stop "fill_bathymetry_ghosts: requires RAKALI_ENABLE_NETCDF=ON "// &
         "at build time."
#endif
   end subroutine fill_bathymetry_ghosts

end module rki_state
