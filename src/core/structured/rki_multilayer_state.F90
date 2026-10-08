!! Multilayer (3D hydrostatic) solution state for the structured solver.
module rki_multilayer_state
   !! Per-layer prognostic arrays, vertical-coordinate definition,
   !! diagnostic vertical velocity, EOS-derived density + reference
   !! profile, tracer registry, and the scalar physics params the
   !! multilayer kernels read (vmix, Smagorinsky, BPG/PPM flags, EOS).
   !! Lives at `state%multilayer`; allocated by `state_init` only when
   !! `use_multilayer = .true.` (else only scalar defaults).
   use, intrinsic :: iso_fortran_env, only: int64
   use rki_constants, only: wp, VCOORD_ZSTAR_FULL, KEPS_STAB_CONSTANT
   use rki_grid, only: hgrid_t
   use rki_vcoord, only: vcoord_t
   use rki_tracer, only: tracer_t
   use rki_mem_report, only: arr_bytes
   use rki_eos, only: eos_t
   use rki_barotropic_state, only: barotropic_state_t
   use rki_ml_tracers, only: ml_seed_tracer
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic, &
                              ml_build_vcoord_zstar_full
   use pic_logger, only: global_logger
   implicit none
   private

   public :: multilayer_state_t

   type :: multilayer_state_t
      !! Per-layer prognostic arrays + multilayer scalar config.

      ! Layer count
      integer :: nz_ml = 1
         !! Number of multilayer sigma layers

      ! Vertical coordinate definition
      type(vcoord_t) :: vcoord
         !! Coordinate type (sigma / zsigma), remap method, parameters
      real(wp), allocatable :: vcoord_z_ref(:)
         !! Reference z-level interfaces (0:nz_ml), GPU-mapped copy of
         !! vcoord%z_ref.  Top-level allocatable so OpenACC can map it
         !! without derived-type attach issues.
      real(wp), allocatable :: vcoord_z_ref_col(:, :, :)
         !! Per-column reference z-level interfaces for VCOORD_ZSTAR_FULL,
         !! shape (0:nz_ml, nx_total, ny_total).  Built once at solver
         !! init from local bathymetry; static thereafter.  Only
         !! allocated when coord_type == VCOORD_ZSTAR_FULL.
      real(wp), allocatable :: dsig_ml(:)
         !! Sigma layer fractions (nz_ml), sum = 1

      ! Layer prognostic arrays (nx_total, ny_total, nz_ml)
      real(wp), allocatable :: h_layer(:, :, :)
         !! Layer thickness
      real(wp), allocatable :: hu_layer(:, :, :)
         !! Layer x-momentum h_k*u_k
      real(wp), allocatable :: hv_layer(:, :, :)
         !! Layer y-momentum h_k*v_k

      ! RK2 saves
      real(wp), allocatable :: h_layer0(:, :, :)
      real(wp), allocatable :: hu_layer0(:, :, :)
      real(wp), allocatable :: hv_layer0(:, :, :)

      ! Diagnostic vertical velocity at interfaces
      real(wp), allocatable :: w_ml(:, :, :)
         !! (nx_total, ny_total, nz_ml+1)

      ! Workspace: barotropic momentum before HLL update
      real(wp), allocatable :: hu_ml_pre(:, :)
      real(wp), allocatable :: hv_ml_pre(:, :)

      ! Generic tracer registry. Salinity and temperature live as the
      ! first two entries; sediment / passive tracers will append here.
      ! Special-shape kernels (EOS, surface heat / virtual-salt flux)
      ! locate S and T via the named indices below.
      type(tracer_t), allocatable :: tracers(:)
         !! Registered prognostic tracers (allocated when use_multilayer)
      integer :: idx_salinity = 0
         !! Index into tracers(:) for salinity. 0 = not registered.
      integer :: idx_temperature = 0
         !! Index into tracers(:) for temperature. 0 = not registered.
      logical :: registry_locked = .false.
         !! Set by `solver_tracers_enter_data`, cleared by
         !! `solver_tracers_exit_data`.  While locked,
         !! `register_passive_tracer` REFUSES: the device map snapshots
         !! `tracers(:)` element-by-element, so growing the registry
         !! reallocates the array (and deep-copies every element's
         !! `hTr`/`hTr0`), stranding every existing tracer's device
         !! mapping on freed host addresses — an illegal address or a
         !! silent stale-host read on the `mem:separate` GPU build.

      ! EOS-derived density
      real(wp), allocatable :: rho_layer(:, :, :)
         !! (nx_total, ny_total, nz_ml)

      ! Reference-density profile. Horizontally-uniform `rho_ref(z)` on
      ! a uniform z-grid `0..rho_ref_dz*nz_ref`, subtracted in both BPG
      ! kernels so the integrand is the anomaly `rho - rho_ref(z)`
      ! (10-100x smaller than `rho - rho_0` in stratified water).  Built
      ! once at solver init (see `ml_compute_reference_density`).
      integer :: nz_ref = 0
         !! Highest valid index in `rho_ref(0:nz_ref)`. 0 = unset.
      real(wp) :: rho_ref_dz = 0.0_wp
         !! Reference-grid spacing (m). 0 = unset.
      real(wp), allocatable :: rho_ref(:)
         !! Reference density profile sampled at `iz * rho_ref_dz`
         !! for `iz = 0..nz_ref`. Index 0 = surface.

      ! ----- Scalar physics parameters -----
      ! Salinity defaults / EOS / clamp / mixing.  Some are duplicated
      ! on the registered tracers (`tracers(idx_salinity)%eos_ref`,
      ! `%tr_min`, `%kappa_bg`, ...).
      real(wp) :: initial_salinity = 35.0_wp
         !! Initial salinity for all layers (PSU)
      real(wp) :: S_init_surface = 0.0_wp
         !! Linear-S(z) stratification, surface (k=nz) value.  When both
         !! S_init_surface and S_init_bottom are non-zero (and nz>1) the
         !! salinity seed becomes a linear profile from S_init_bottom at
         !! the bed (k=1) to S_init_surface at the surface; otherwise the
         !! uniform initial_salinity is used.  Stable polarity is the
         !! inverse of temperature: S_init_bottom > S_init_surface (salty
         !! bed, fresh surface).
      real(wp) :: S_init_bottom = 0.0_wp
         !! Linear-S(z) stratification, bed (k=1) value.  See S_init_surface.
      real(wp) :: S_ref = 0.0_wp
         !! Reference salinity in EOS (PSU)
      real(wp) :: beta_S = 0.78_wp
         !! Haline contraction coefficient (kg/m^3 per PSU)
      real(wp) :: S_min = 0.0_wp
         !! Lower physical salinity bound (PSU)
      real(wp) :: S_max = 40.0_wp
         !! Upper physical salinity bound (PSU)
      real(wp) :: kappa_S_bg = 1.0e-5_wp
         !! Background vertical salinity diffusivity (m^2/s)

      ! Temperature equivalents
      real(wp) :: initial_temperature = 15.0_wp
         !! Initial temperature for all layers (degC)
      real(wp) :: T_init_surface = 0.0_wp
         !! Linear-T(z) stratification, surface (k=nz) value.  When both
         !! T_init_surface and T_init_bottom are non-zero (and nz>1) the
         !! temperature seed becomes a linear profile from T_init_bottom at
         !! the bed (k=1) to T_init_surface at the surface; otherwise the
         !! uniform initial_temperature is used.  Mirrors the ocean path.
      real(wp) :: T_init_bottom = 0.0_wp
         !! Linear-T(z) stratification, bed (k=1) value.  See T_init_surface.
      real(wp) :: T_ref = 15.0_wp
         !! Reference temperature in EOS (degC)
      real(wp) :: alpha_T = 0.17_wp
         !! Thermal expansion coefficient (kg/m^3 per degC).
         !! Set to 0 to make temperature a passive tracer.
      real(wp) :: T_min = -2.0_wp
         !! Lower physical temperature bound (degC)
      real(wp) :: T_max = 40.0_wp
         !! Upper physical temperature bound (degC)
      real(wp) :: kappa_T_bg = 1.0e-5_wp
         !! Background vertical temperature diffusivity (m^2/s)

      ! EOS reference density (also used as a Boussinesq scaling)
      real(wp) :: rho_0 = 1000.0_wp
         !! Reference density for EOS (kg/m^3)

      ! Shared (ocean) nonlinear EOS — opt-in via `&tracer_nml eos_form`
      logical :: eos_use_shared = .false.
         !! When .true. the EOS call sites route through the shared
         !! `rki_eos` dispatch (`ml_equation_of_state_shared`)
         !! instead of the legacy coastal linear kernel.  Default off ⇒
         !! every existing coastal config is bit-identical by
         !! construction.  Set from `&tracer_nml eos_form` via
         !! `ml_parse_eos_form` (state_init / driver).
      type(eos_t) :: eos
         !! Shared EOS handle (variant tag + scalar coefficients) for
         !! the opt-in path.  Flat POD (no allocatable components):
         !! read on the host at dispatch time and passed by value into
         !! the device kernels, so it needs NO enter_data/exit_data
         !! mapping.  `eos%rho0` mirrors `rho_0` (vanishing-layer
         !! fallback density); `eos%p_ref` stays 0 ⇒ surface/potential
         !! density.

      ! Pacanowski-Philander vertical mixing coefficients
      real(wp) :: vmix_nu_max = 0.01_wp
         !! Pacanowski-Philander max turbulent viscosity (m^2/s)
      real(wp) :: vmix_nu_bg = 1.0e-4_wp
         !! Background vertical viscosity (m^2/s)
      real(wp) :: vmix_ri_alpha = 5.0_wp
         !! Richardson number scaling coefficient (shared S+T in PP81)

      ! KPP (Large et al. 1994) surface-boundary-layer overlay.  When
      ! `vmix_use_kpp` is on, the solver swaps the three PP81 vmix calls
      ! for the KPP kernels (`rki_ml_kpp`).  The 2D workspace fields
      ! below are allocated only when `vmix_use_kpp`.
      logical :: vmix_use_kpp = .false.
         !! Master switch for the KPP overlay (default off = pure PP81).
      logical :: vmix_use_keps = .false.
         !! Master switch for the two-equation k-epsilon closure (default
         !! off = PP81/KPP).  Interface-shaped prognostic fields below
         !! are allocated only when this is set.
      real(wp), allocatable :: keps_tke(:, :, :)
         !! Turbulent kinetic energy at interfaces (m^2/s^2), shape (nx,ny,nz+1).
      real(wp), allocatable :: keps_eps(:, :, :)
         !! TKE dissipation rate at interfaces (m^2/s^3).
      real(wp), allocatable :: keps_nut(:, :, :)
         !! k-eps turbulent viscosity at interfaces (m^2/s); read as the old
         !! value by the next diffusion step, so it persists across steps.
      real(wp), allocatable :: keps_kt(:, :, :)
         !! k-eps turbulent (scalar) diffusivity at interfaces (m^2/s).
      real(wp), allocatable :: keps_ustar_b(:, :)
         !! Bed friction velocity for the k-eps wall functions (m/s), (nx,ny).
         !! Recomputed each step from the bottom-drag law; kept as a field
         !! (not scratch) so it maps to the device like the KPP u_star.
      real(wp), allocatable :: keps_ustar_s(:, :)
         !! Surface friction velocity for the k-eps wall functions (m/s).
      integer :: keps_stability = KEPS_STAB_CONSTANT
         !! Stability-function scheme for nu_t / kappa_t (KEPS_STAB_*): constant
         !! c_mu (default) or Galperin stratification-aware quasi-equilibrium.
      real(wp) :: kpp_ri_crit = 0.3_wp
         !! Critical bulk Richardson number for the BL-depth sweep.
      real(wp) :: kpp_cs_nonlocal = 6.3_wp
         !! Non-local transport coefficient C_s (LMD94 eq 20).
      real(wp) :: kpp_c_vt2 = 1.8_wp
         !! Unresolved-turbulence coefficient for V_t^2 (LMD94 eq 23).
      real(wp), allocatable :: kpp_bl_depth(:, :)
         !! Diagnosed boundary-layer depth (m, positive down).  Persists
         !! across steps — the previous value seeds the V_t^2 velocity
         !! scale (self-bootstraps from zero on the first call).
      real(wp), allocatable :: kpp_u_star(:, :)
         !! Surface friction velocity sqrt(|tau|/rho_0) (m/s), per column.
      real(wp), allocatable :: kpp_buoy_flux(:, :)
         !! Surface buoyancy flux B_0 (m^2/s^3); <0 = destabilizing.

      ! Horizontal mixing
      real(wp) :: hdiff_kappa = 0.0_wp
         !! Horizontal diffusion coefficient for tracers (m^2/s)
      real(wp) :: smagorinsky_cs = 0.0_wp
         !! Smagorinsky horizontal viscosity coefficient
      real(wp) :: biharmonic_nu4 = 0.0_wp
         !! Biharmonic (del^4) lateral viscosity coefficient (m^4/s); 0 = off

      ! Scheme flags
      logical :: use_bpg_jacobian = .true.
         !! Use the pressure-Jacobian BPG (Shchepetkin & McWilliams 2003).
         !! DEFAULT ON.  The sigma alternative differences the integrated
         !! density anomaly at constant SIGMA INDEX, which over sloping
         !! bathymetry produces a fictitious force proportional to
         !! (uniform density anomaly) x (bed slope) -- so a column of
         !! spatially UNIFORM density whose value differs from `rho_0` is
         !! spuriously accelerated, and the resulting flow leaks mass
         !! (measured: 1.7e-6 relative over 400 steps, vs 1.3e-14 with the
         !! Jacobian).  Pinned by `test_bpg_uniform_density`.
      logical :: use_ppm_tracer = .true.
         !! Use PPM for tracer advection (false = first-order upwind)
      logical :: evap_mass_sink = .false.
         !! When true, net evaporation is a real top-layer water-mass sink
         !! (h_top -= E*dt) with salt held (concentration rises) and heat
         !! carried away with the water.  Default false = virtual salt flux.
      logical :: evap_latent_cooling = .false.
         !! When true (needs evap_mass_sink), evaporation also withdraws its
         !! latent heat of vaporization from the surface layer (k=nz) — an
         !! additional cooling on top of the sensible heat that leaves with
         !! the removed water.  Default false = bit-identical.
      integer :: tracer_recon_code = 0
         !! Parsed TRACER_RECON_* code (0 = PPM, 1 = WENO5, 2 = WENO7,
         !! 3 = WENO9).  Set at state init from config%tracer_recon.

      ! Sediment transport (single class, default off = bit-identical).
      ! Suspended sediment is a registered tracer (idx_sediment) so its
      ! advection / diffusion / clamp come from the registry; settling +
      ! bed exchange are separate kernels.  Bed is a 2D mass reservoir.
      logical :: use_sediment = .false.
         !! Master switch for sediment transport (default off).
      real(wp) :: sediment_w_s = 0.0_wp
         !! Settling (sinking) velocity (m/s, downward; 0 = no settling).
      real(wp) :: sediment_initial_conc = 0.0_wp
         !! Uniform initial suspended-sediment concentration (kg/m^3).
      real(wp) :: sediment_bed_init = 0.0_wp
         !! Uniform initial bed sediment mass (kg/m^2).
      real(wp) :: sediment_tau_crit_erosion = 0.0_wp
         !! Critical bed shear stress for erosion (Pa; 0 = erosion off).
      real(wp) :: sediment_erosion_rate = 0.0_wp
         !! Ariathurai erosion-rate coefficient E0 (kg/m^2/s; 0 = off).
      real(wp) :: sediment_porosity = 0.0_wp
         !! Bed porosity (0-1); erosion scaled by (1 - porosity).
      integer :: idx_sediment = 0
         !! Index into tracers(:) for suspended sediment. 0 = not registered.
      real(wp), allocatable :: bed_sediment(:, :)
         !! Bed sediment mass reservoir (kg/m^2), per column.  Allocated when
         !! use_sediment.  Settling deposits into it; total (suspended + bed)
         !! sediment is conserved.

   contains
      procedure, non_overridable :: init => multilayer_state_init
      procedure, non_overridable :: destroy => multilayer_state_destroy
      procedure, non_overridable :: register_passive_tracer => multilayer_register_passive_tracer
      procedure, non_overridable :: init_layers => multilayer_init_layers
      procedure, non_overridable :: build_vcoord_reference => multilayer_build_vcoord_reference
      procedure, non_overridable :: seed_tracers => multilayer_seed_tracers
      procedure, non_overridable :: bytes => multilayer_state_bytes
   end type multilayer_state_t

contains

   subroutine multilayer_state_init(this, grid)
      !! Allocate all multilayer arrays at the grid size and configured
      !! layer count, init the vcoord object, copy `vcoord%z_ref` into
      !! the GPU-mappable `vcoord_z_ref`, allocate the rho-ref profile,
      !! and register S + T.  Caller must set `nz_ml`,
      !! `vcoord%coord_type`, `vcoord%remap_method`, and any other
      !! vcoord parameters first.
      class(multilayer_state_t), intent(inout) :: this
      type(hgrid_t), intent(in)    :: grid

      integer :: nx, ny, nz_ml

      nx = grid%nx_total
      ny = grid%ny_total
      nz_ml = this%nz_ml

      call this%vcoord%init(nz_ml, this%vcoord%coord_type, &
                            this%vcoord%remap_method)
      allocate (this%dsig_ml(nz_ml))
      this%dsig_ml = this%vcoord%dsig_target

      ! GPU-mappable copy of z_ref (avoids derived-type attach issues)
      allocate (this%vcoord_z_ref(0:nz_ml), source=0.0_wp)
      if (allocated(this%vcoord%z_ref)) then
         this%vcoord_z_ref = this%vcoord%z_ref
      end if

      ! Full MOM6 z* per-column z_ref.  Allocated here; the per-cell
      ! contents are filled later from the bathymetry by
      ! `ml_build_vcoord_zstar_full` (called from rki_init).
      if (this%vcoord%coord_type == VCOORD_ZSTAR_FULL) then
         allocate (this%vcoord_z_ref_col(0:nz_ml, nx, ny), source=0.0_wp)
      end if

      allocate (this%h_layer(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%hu_layer(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%hv_layer(nx, ny, nz_ml), source=0.0_wp)

      allocate (this%h_layer0(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%hu_layer0(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%hv_layer0(nx, ny, nz_ml), source=0.0_wp)

      allocate (this%w_ml(nx, ny, nz_ml + 1), source=0.0_wp)

      allocate (this%hu_ml_pre(nx, ny), source=0.0_wp)
      allocate (this%hv_ml_pre(nx, ny), source=0.0_wp)

      allocate (this%rho_layer(nx, ny, nz_ml), source=0.0_wp)

      ! KPP boundary-layer workspace — only when the overlay is enabled
      ! (caller sets `vmix_use_kpp` before this routine).  Left
      ! unallocated on the default PP81 path.
      if (this%vmix_use_kpp) then
         allocate (this%kpp_bl_depth(nx, ny), source=0.0_wp)
         allocate (this%kpp_u_star(nx, ny), source=0.0_wp)
         allocate (this%kpp_buoy_flux(nx, ny), source=0.0_wp)
      end if

      ! k-epsilon interface-shaped prognostics — only when the closure is on.
      ! Seeded with small positive floors so nu_t = c_mu k^2/eps is finite on
      ! the first step (the closure self-bootstraps from these).
      if (this%vmix_use_keps) then
         allocate (this%keps_tke(nx, ny, nz_ml + 1), source=1.0e-6_wp)
         allocate (this%keps_eps(nx, ny, nz_ml + 1), source=1.0e-8_wp)
         allocate (this%keps_nut(nx, ny, nz_ml + 1), source=1.0e-4_wp)
         allocate (this%keps_kt(nx, ny, nz_ml + 1), source=1.0e-4_wp)
         allocate (this%keps_ustar_b(nx, ny), source=0.0_wp)
         allocate (this%keps_ustar_s(nx, ny), source=0.0_wp)
      end if

      ! Bed sediment reservoir — only when sediment transport is enabled
      ! (caller sets use_sediment + sediment_bed_init before this routine).
      if (this%use_sediment) then
         allocate (this%bed_sediment(nx, ny), source=this%sediment_bed_init)
      end if

      ! Reference-density profile: 100 samples, sub-metre resolution.
      this%nz_ref = 100
      allocate (this%rho_ref(0:this%nz_ref), source=0.0_wp)
      this%rho_ref_dz = 0.0_wp

      ! Default tracer registry: salinity at index 1, temperature at
      ! index 2.  Scalar config + CF metadata are overlaid by
      ! `register_default_tracers` (rki_state) on the namelist path.
      allocate (this%tracers(2))
      this%idx_salinity = 1
      this%idx_temperature = 2
      call this%tracers(this%idx_salinity)%init(grid, nz_ml)
      call this%tracers(this%idx_temperature)%init(grid, nz_ml)
      this%tracers(this%idx_salinity)%name = "salinity"
      this%tracers(this%idx_salinity)%long_name = "Sea water salinity"
      this%tracers(this%idx_salinity)%units = "PSU"
      this%tracers(this%idx_salinity)%standard_name = "sea_water_salinity"
      this%tracers(this%idx_temperature)%name = "temperature"
      this%tracers(this%idx_temperature)%long_name = "Sea water potential temperature"
      this%tracers(this%idx_temperature)%units = "degC"
      this%tracers(this%idx_temperature)%standard_name = "sea_water_potential_temperature"
   end subroutine multilayer_state_init

   subroutine multilayer_register_passive_tracer(this, grid, name, tr_init, &
                                                 tr_min, tr_max, kappa_bg, &
                                                 hdiff_kappa)
      !! Append a passive tracer (eos_coeff = 0) to the registry,
      !! growing `tracers(:)` past the default S/T pair.  Must be called
      !! AFTER `init` and BEFORE the solver's `enter_data` (so the new
      !! `hTr` gets device-mapped) — enforced via `registry_locked`:
      !! once the solver has mapped the registry the call logs an error
      !! and returns without touching it.  Caller populates `hTr` once
      !! layer thicknesses are known.  S/T keep their indices.
      class(multilayer_state_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: tr_init, tr_min, tr_max, kappa_bg, hdiff_kappa

      type(tracer_t), allocatable :: tmp(:)
      integer :: n_old, it

      if (this%registry_locked) then
         call global_logger%error("multilayer_state: register_passive_tracer("// &
                                  trim(name)//") refused — registry locked by "// &
                                  "enter_data (growing the registry now would strand "// &
                                  "every existing tracer's device mapping on the "// &
                                  "mem:separate GPU build)")
         return
      end if

      n_old = size(this%tracers)
      call move_alloc(this%tracers, tmp)
      allocate (this%tracers(n_old + 1))
      ! Deep-copy the existing tracers (intrinsic assignment reallocates
      ! the allocatable hTr / hTr0 components and copies their data).
      do it = 1, n_old
         this%tracers(it) = tmp(it)
      end do

      it = n_old + 1
      call this%tracers(it)%init(grid, this%nz_ml)
      this%tracers(it)%name = trim(name)
      this%tracers(it)%long_name = trim(name)
      this%tracers(it)%units = "1"
      this%tracers(it)%standard_name = ""
      this%tracers(it)%tr_init = tr_init
      this%tracers(it)%tr_inflow = tr_init
      this%tracers(it)%tr_min = tr_min
      this%tracers(it)%tr_max = tr_max
      this%tracers(it)%kappa_bg = kappa_bg
      this%tracers(it)%hdiff_kappa = hdiff_kappa
      this%tracers(it)%eos_coeff = 0.0_wp
      this%tracers(it)%eos_ref = 0.0_wp
   end subroutine multilayer_register_passive_tracer

   subroutine multilayer_build_vcoord_reference(this, barotropic, grid)
      !! For VCOORD_ZSTAR_FULL, build the per-column reference z-level
      !! pattern `vcoord_z_ref_col` from the bathymetry.  No-op for the
      !! other vertical coordinates.  Needed both on a cold start (called
      !! by `init_layers`) and after a restart read, which restores the
      !! prognostic layer fields but not this derived array.
      class(multilayer_state_t), intent(inout) :: this
      type(barotropic_state_t), intent(in) :: barotropic
      type(hgrid_t), intent(in) :: grid

      if (this%vcoord%coord_type == VCOORD_ZSTAR_FULL .and. &
          allocated(this%vcoord_z_ref_col)) then
         call ml_build_vcoord_zstar_full(barotropic%b, this%vcoord_z_ref_col, &
                                         grid, this%nz_ml, &
                                         this%vcoord%zstar_h_surf_target, &
                                         this%vcoord%zstar_n_surf, &
                                         this%vcoord%zstar_stretching)
      end if
   end subroutine multilayer_build_vcoord_reference

   subroutine multilayer_init_layers(this, barotropic, grid)
      !! Lay out the initial per-layer distribution from the barotropic
      !! state: for VCOORD_ZSTAR_FULL first build the per-column
      !! reference z-level pattern from the bathymetry, then split the
      !! barotropic (h, hu, hv) across layers via the sigma weights
      !! `dsig_ml`.  Call once after `init` and before `seed_tracers` /
      !! EOS.
      class(multilayer_state_t), intent(inout) :: this
      type(barotropic_state_t), intent(in) :: barotropic
      type(hgrid_t), intent(in) :: grid

      call this%build_vcoord_reference(barotropic, grid)
      call ml_init_layers_from_barotropic(barotropic%h, barotropic%hu, barotropic%hv, &
                                          this%dsig_ml, this%h_layer, this%hu_layer, &
                                          this%hv_layer, grid, this%nz_ml)

   end subroutine multilayer_init_layers

   subroutine multilayer_seed_tracers(this, grid)
      !! Seed salinity + temperature from the IC scalars on the state,
      !! each dispatching uniform-vs-stratified via `ml_seed_tracer`.
      !! Call once after the layer thicknesses are laid out and before
      !! EOS.
      class(multilayer_state_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid

      call ml_seed_tracer(this%h_layer, &
                          this%tracers(this%idx_salinity)%hTr, &
                          this%initial_salinity, &
                          this%S_init_bottom, this%S_init_surface, &
                          grid, this%nz_ml)
      call ml_seed_tracer(this%h_layer, &
                          this%tracers(this%idx_temperature)%hTr, &
                          this%initial_temperature, &
                          this%T_init_bottom, this%T_init_surface, &
                          grid, this%nz_ml)

   end subroutine multilayer_seed_tracers

   subroutine multilayer_state_destroy(this)
      class(multilayer_state_t), intent(inout) :: this

      integer :: it

      call this%vcoord%cleanup()
      if (allocated(this%vcoord_z_ref)) deallocate (this%vcoord_z_ref)
      if (allocated(this%vcoord_z_ref_col)) deallocate (this%vcoord_z_ref_col)
      if (allocated(this%dsig_ml)) deallocate (this%dsig_ml)
      if (allocated(this%h_layer)) deallocate (this%h_layer)
      if (allocated(this%hu_layer)) deallocate (this%hu_layer)
      if (allocated(this%hv_layer)) deallocate (this%hv_layer)
      if (allocated(this%h_layer0)) deallocate (this%h_layer0)
      if (allocated(this%hu_layer0)) deallocate (this%hu_layer0)
      if (allocated(this%hv_layer0)) deallocate (this%hv_layer0)
      if (allocated(this%w_ml)) deallocate (this%w_ml)
      if (allocated(this%hu_ml_pre)) deallocate (this%hu_ml_pre)
      if (allocated(this%hv_ml_pre)) deallocate (this%hv_ml_pre)
      if (allocated(this%rho_layer)) deallocate (this%rho_layer)
      if (allocated(this%keps_tke)) deallocate (this%keps_tke)
      if (allocated(this%keps_eps)) deallocate (this%keps_eps)
      if (allocated(this%keps_nut)) deallocate (this%keps_nut)
      if (allocated(this%keps_kt)) deallocate (this%keps_kt)
      if (allocated(this%keps_ustar_b)) deallocate (this%keps_ustar_b)
      if (allocated(this%keps_ustar_s)) deallocate (this%keps_ustar_s)
      if (allocated(this%kpp_bl_depth)) deallocate (this%kpp_bl_depth)
      if (allocated(this%kpp_u_star)) deallocate (this%kpp_u_star)
      if (allocated(this%kpp_buoy_flux)) deallocate (this%kpp_buoy_flux)
      if (allocated(this%bed_sediment)) deallocate (this%bed_sediment)
      if (allocated(this%rho_ref)) deallocate (this%rho_ref)
      this%nz_ref = 0
      this%rho_ref_dz = 0.0_wp

      if (allocated(this%tracers)) then
         do it = 1, size(this%tracers)
            call this%tracers(it)%destroy()
         end do
         deallocate (this%tracers)
      end if
      this%idx_salinity = 0
      this%idx_temperature = 0
      this%registry_locked = .false.
   end subroutine multilayer_state_destroy

   pure function multilayer_state_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the multilayer slot: the layer
      !! prognostics + RK2 saves + vcoord/rho tables + the tracer registry
      !! (each tracer sums its own arrays) + the gated KPP / k-epsilon /
      !! sediment workspaces (0 when their closure is off).
      class(multilayer_state_t), intent(in) :: this
      integer(int64) :: nbytes
      integer :: it

      nbytes = arr_bytes(this%vcoord_z_ref) + arr_bytes(this%vcoord_z_ref_col) &
               + arr_bytes(this%dsig_ml) &
               + arr_bytes(this%h_layer) + arr_bytes(this%hu_layer) &
               + arr_bytes(this%hv_layer) &
               + arr_bytes(this%h_layer0) + arr_bytes(this%hu_layer0) &
               + arr_bytes(this%hv_layer0) &
               + arr_bytes(this%w_ml) &
               + arr_bytes(this%hu_ml_pre) + arr_bytes(this%hv_ml_pre) &
               + arr_bytes(this%rho_layer) + arr_bytes(this%rho_ref) &
               + arr_bytes(this%keps_tke) + arr_bytes(this%keps_eps) &
               + arr_bytes(this%keps_nut) + arr_bytes(this%keps_kt) &
               + arr_bytes(this%keps_ustar_b) + arr_bytes(this%keps_ustar_s) &
               + arr_bytes(this%kpp_bl_depth) + arr_bytes(this%kpp_u_star) &
               + arr_bytes(this%kpp_buoy_flux) + arr_bytes(this%bed_sediment)

      if (allocated(this%tracers)) then
         do it = 1, size(this%tracers)
            nbytes = nbytes + this%tracers(it)%bytes()
         end do
      end if
   end function multilayer_state_bytes

end module rki_multilayer_state
