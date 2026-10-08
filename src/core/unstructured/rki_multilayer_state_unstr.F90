!! Multilayer (3D hydrostatic) solution state for the unstructured solver.
module rki_multilayer_state_unstr
   !! Per-layer prognostic arrays (`(nz_ml, ncells_total)` axis order —
   !! layer index fastest, matching column-stencil kernels), vcoord
   !! definition + GPU mirror, diagnostic vertical velocity, EOS density
   !! + reference profile, tracer registry, per-step Jacobi/gradient
   !! workspaces, and the scalar physics parameters the kernels read.
   !!
   !! Lives at `state%multilayer`. Allocated by `state_unstr_init` when
   !! `nz_ml` is passed. Mirrors structured `multilayer_state_t` with the
   !! cell-centric axis order swapped.
   use rki_constants, only: wp, VCOORD_ZSTAR_FULL, KEPS_STAB_CONSTANT
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_vcoord, only: vcoord_t
   use rki_tracer_unstr, only: unstr_tracer_t
   use rki_eos, only: eos_t
   use rki_boundary_unstr, only: MAX_BC_GROUPS
   use pic_logger, only: global_logger
   implicit none
   private

   public :: unstr_multilayer_state_t

   type :: unstr_multilayer_state_t
      !! Per-layer prognostic arrays + tracer registry + multilayer scalar config.

      ! Layer count
      integer :: nz_ml = 1
         !! Number of multilayer sigma layers.

      ! Vertical coordinate definition + GPU mirror
      type(vcoord_t) :: vcoord
         !! Coordinate type (sigma / zsigma), remap method, parameters.
      real(wp), allocatable :: dsig_ml(:)
         !! Sigma layer fractions (nz_ml), sum = 1.
      real(wp), allocatable :: vcoord_z_ref(:)
         !! Reference z-level interfaces (0:nz_ml), GPU-mapped copy of
         !! `vcoord%z_ref`.  Top-level allocatable so OpenACC can map
         !! it without nested derived-type pointer issues.
      real(wp), allocatable :: vcoord_z_ref_col(:, :)
         !! Per-cell reference z-level interfaces for VCOORD_ZSTAR_FULL,
         !! shape (0:nz_ml, ncells_total).  Built once at solver init
         !! from local bathymetry; static thereafter.  Only allocated
         !! when coord_type == VCOORD_ZSTAR_FULL.

      ! Layer prognostic arrays (nz_ml, ncells_total)
      real(wp), allocatable :: h_layer(:, :)
      real(wp), allocatable :: hu_layer(:, :)
      real(wp), allocatable :: hv_layer(:, :)

      ! RK2 saves (dynamics only — tracer RK2 saves live on tracers(it)%hTr0)
      real(wp), allocatable :: h_layer0(:, :)
      real(wp), allocatable :: hu_layer0(:, :)
      real(wp), allocatable :: hv_layer0(:, :)

      ! EOS-derived density
      real(wp), allocatable :: rho_layer(:, :)

      ! Diagnostic vertical velocity at layer interfaces
      real(wp), allocatable :: w_interface(:, :)
         !! `(nz_ml + 1, ncells_total)`. Index 1 = bed, nz_ml+1 = surface.

      ! Generic tracer registry. Salinity / temperature live as the
      ! first two entries by default; sediment / passive tracers append.
      type(unstr_tracer_t), allocatable :: tracers(:)
      integer :: idx_salinity = 0
         !! Index into tracers(:) for salinity. 0 = not registered.
      integer :: idx_temperature = 0
         !! Index into tracers(:) for temperature. 0 = not registered.
      logical :: registry_locked = .false.
         !! Set by `tracers_enter_data` (rki_state_unstr), cleared by
         !! `tracers_exit_data`.  While locked,
         !! `register_passive_tracer` REFUSES: the device map snapshots
         !! `tracers(:)` element-by-element, so growing the registry
         !! reallocates the array (and deep-copies every element's
         !! `hTr`/`hTr0`), stranding every existing tracer's device
         !! mapping on freed host addresses — an illegal address or a
         !! silent stale-host read on the `mem:separate` GPU build.

      ! Reference-density profile. Subtracted in both BPG kernels so the
      ! integrand is `rho - rho_ref(z)` rather than `rho - rho_0`
      ! -- 10-100x smaller in stratified water.
      integer :: nz_ref = 0
      real(wp) :: rho_ref_dz = 0.0_wp
      real(wp), allocatable :: rho_ref(:)
         !! Reference density profile sampled at `iz * rho_ref_dz` for
         !! `iz = 0..nz_ref`. Index 0 = surface.

      ! Per-step kernel workspaces (allocated alongside the prognostic
      ! arrays; mapped at enter_data so kernels never allocate per step).
      real(wp), allocatable :: tracer_advect_scratch(:, :)
         !! Jacobi snapshot for `ml_advect_tracer_unstr` / MUSCL.
      real(wp), allocatable :: ml_lgrad_x(:, :)
         !! Limited Green-Gauss gradient (x) for MUSCL advection.
      real(wp), allocatable :: ml_lgrad_y(:, :)
      real(wp), allocatable :: ml_hdiff_hTr_pre(:, :)
         !! Jacobi snapshot for the horizontal-diffusion kernel.
      real(wp), allocatable :: ml_horvisc_hu_pre(:, :)
         !! Jacobi snapshot for the Smagorinsky horizontal-viscosity kernel.
      real(wp), allocatable :: ml_horvisc_hv_pre(:, :)

      ! Barotropic momentum snapshot (saved before HLL update so
      ! `ml_distribute_barotropic_unstr` can spread the increment).
      real(wp), allocatable :: xmom_pre(:)
      real(wp), allocatable :: ymom_pre(:)

      ! Saved layer thicknesses before vertical-coordinate update.
      ! Only allocated when vcoord%needs_remap().
      real(wp), allocatable :: remap_dz_old(:, :)

      ! ----- Scalar physics parameters -----
      real(wp) :: initial_salinity = 35.0_wp
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
      real(wp) :: beta_S = 0.78_wp
      real(wp) :: S_min = 0.0_wp
      real(wp) :: S_max = 40.0_wp
      real(wp) :: kappa_S_bg = 1.0e-5_wp

      real(wp) :: initial_temperature = 15.0_wp
      real(wp) :: T_init_surface = 0.0_wp
         !! Linear-T(z) stratification, surface (k=nz) value.  When both
         !! T_init_surface and T_init_bottom are non-zero (and nz>1) the
         !! temperature seed becomes a linear profile from T_init_bottom at
         !! the bed (k=1) to T_init_surface at the surface; otherwise the
         !! uniform initial_temperature is used.  Mirrors the ocean path.
      real(wp) :: T_init_bottom = 0.0_wp
         !! Linear-T(z) stratification, bed (k=1) value.  See T_init_surface.
      real(wp) :: T_ref = 15.0_wp
      real(wp) :: alpha_T = 0.17_wp
      real(wp) :: T_min = -2.0_wp
      real(wp) :: T_max = 40.0_wp
      real(wp) :: kappa_T_bg = 1.0e-5_wp

      real(wp) :: rho_0 = 1000.0_wp
         !! Reference density for the linear EOS / Boussinesq scaling.

      ! Shared (ocean) nonlinear EOS — opt-in via `&tracer_nml eos_form`
      logical :: eos_use_shared = .false.
         !! When .true. `ml_equation_of_state_unstr` routes through the
         !! shared `rki_eos` point evaluation instead of the
         !! legacy linear kernel.  Default off ⇒ bit-identical.  Set
         !! from `&tracer_nml eos_form` via `ml_parse_eos_form`
         !! (driver).
      type(eos_t) :: eos
         !! Shared EOS handle (variant tag + scalar coefficients) for
         !! the opt-in path.  Flat POD — passed by value into the
         !! device kernel, no enter_data/exit_data mapping needed.
         !! `eos%rho0` mirrors `rho_0`; `eos%p_ref` stays 0 ⇒
         !! surface/potential density.

      ! Pacanowski-Philander vertical mixing
      real(wp) :: vmix_nu_max = 0.01_wp
      real(wp) :: vmix_nu_bg = 1.0e-4_wp
      real(wp) :: vmix_ri_alpha = 5.0_wp

      ! KPP (Large et al. 1994) surface-boundary-layer overlay. When
      ! `vmix_use_kpp` is on, the solver swaps PP81 vmix for the KPP
      ! kernels (`rki_ml_kpp_unstr`): bulk-Ri BL-depth sweep, shape-
      ! function-enhanced diffusivity over PP81, non-local counter-
      ! gradient tracer transport. Per-cell fields allocated
      ! unconditionally (API can flip vmix_use_kpp after init; gating the
      ! allocation would leave them unmapped).
      logical :: vmix_use_kpp = .false.
         !! Master switch for the KPP overlay (default off = pure PP81).
      real(wp) :: kpp_ri_crit = 0.3_wp
         !! Critical bulk Richardson number for the BL-depth sweep.
      real(wp) :: kpp_cs_nonlocal = 6.3_wp
         !! Non-local transport coefficient C_s (LMD94 eq 20).
      real(wp) :: kpp_c_vt2 = 1.8_wp
         !! Unresolved-turbulence coefficient for V_t^2 (LMD94 eq 23).
      real(wp), allocatable :: kpp_bl_depth(:)
         !! Diagnosed boundary-layer depth (m, positive down), per cell.
         !! Persists across steps — the previous value seeds the V_t^2
         !! velocity scale (self-bootstraps from zero on the first call).
      real(wp), allocatable :: kpp_u_star(:)
         !! Surface friction velocity sqrt(|tau|/rho_0) (m/s), per cell.
      real(wp), allocatable :: kpp_buoy_flux(:)
         !! Surface buoyancy flux B_0 (m^2/s^3); <0 = destabilizing.

      ! k-epsilon (two-equation GLS) vertical mixing. When `vmix_use_keps`
      ! is on, the solver swaps PP81/KPP for `rki_ml_keps_unstr`. Allocated
      ! unconditionally (mirrors KPP). Default off = bit-identical.
      logical :: vmix_use_keps = .false.
         !! Master switch for the k-epsilon closure (default off).
      integer :: keps_stability = KEPS_STAB_CONSTANT
         !! Stability-function scheme: constant / galperin / canuto.
      integer :: keps_si_picard = 1
         !! Picard iteration count for the SI k-eps solve (1 = single pass;
         !! 2+ enables save/restore of the pre-step tke/eps and re-runs
         !! ml_keps_update_unstr from the same initial state with the
         !! refreshed nu_t production coefficient).  Mirrors the structured
         !! `keps_si_picard` knob.
      real(wp), allocatable :: keps_tke0(:, :)
         !! Pre-step TKE save slot (nz+1, ncells_total); allocated only when
         !! vmix_use_keps .and. keps_si_picard > 1.
      real(wp), allocatable :: keps_eps0(:, :)
         !! Pre-step dissipation-rate save slot (nz+1, ncells_total); co-
         !! allocated with keps_tke0.
      real(wp), allocatable :: keps_tke(:, :)
         !! Turbulent kinetic energy at interfaces (nz+1, ncell); state.
      real(wp), allocatable :: keps_eps(:, :)
         !! Dissipation rate at interfaces (nz+1, ncell); state.
      real(wp), allocatable :: keps_nut(:, :)
         !! Eddy viscosity at interfaces (nz+1, ncell); carried across steps.
      real(wp), allocatable :: keps_kt(:, :)
         !! Scalar (tracer) diffusivity at interfaces (nz+1, ncell); carried.
      real(wp), allocatable :: keps_ustar_s(:)
         !! Surface friction velocity sqrt(|tau_wind|/rho_0) (m/s), per cell.
      real(wp), allocatable :: keps_ustar_b(:)
         !! Bed friction velocity (m/s), per cell.

      ! Horizontal mixing
      real(wp) :: hdiff_kappa = 0.0_wp
      real(wp) :: smagorinsky_cs = 0.0_wp

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
         !! Use MUSCL for tracer advection (false = first-order upwind).
      logical :: evap_mass_sink = .false.
         !! When true, net evaporation is a real water-mass sink removed
         !! COLUMN-PROPORTIONALLY (salt held, heat carried away with the
         !! water) rather than a virtual salt flux.  Mirrors the structured
         !! `evap_mass_sink`.  Default false = bit-identical (virtual salt).
      logical :: evap_latent_cooling = .false.
         !! When true (needs evap_mass_sink), evaporation additionally
         !! withdraws its latent heat of vaporization from the surface layer
         !! (k=nz).  Mirrors the structured `evap_latent_cooling`.  Default
         !! false = bit-identical.

      ! ----- Per-BC-group tracer inflow values -----
      real(wp) :: bc_S_ext(MAX_BC_GROUPS) = 35.0_wp
      real(wp) :: bc_T_ext(MAX_BC_GROUPS) = 15.0_wp
   contains
      procedure, non_overridable :: init => unstr_multilayer_init
      procedure, non_overridable :: destroy => unstr_multilayer_destroy
      procedure, non_overridable :: enter_data => unstr_multilayer_enter_data
      procedure, non_overridable :: exit_data => unstr_multilayer_exit_data
      procedure, non_overridable :: register_passive_tracer => unstr_multilayer_register_passive_tracer
   end type unstr_multilayer_state_t

contains

   subroutine unstr_multilayer_init(this, mesh, nz_ml)
      !! Allocate every multilayer array at the cell counts in `mesh`
      !! and the configured layer count.  Initialises `vcoord` from
      !! its (already-set) `coord_type` / `remap_method`, allocates
      !! the rho-ref profile, and registers the salinity + temperature
      !! tracers in the tracer registry.
      !!
      !! Caller must have set `this%nz_ml` (or pass it via `nz_ml`),
      !! `this%vcoord%coord_type` / `remap_method`, and any other
      !! vcoord parameters before calling.
      class(unstr_multilayer_state_t), intent(inout) :: this
      type(mesh_unstr_t), intent(in)    :: mesh
      integer, intent(in)    :: nz_ml

      integer :: n, nz

      n = mesh%ncells_total
      nz = nz_ml
      this%nz_ml = nz

      call this%vcoord%init(nz, this%vcoord%coord_type, &
                            this%vcoord%remap_method)
      allocate (this%dsig_ml(nz))
      this%dsig_ml = this%vcoord%dsig_target

      ! GPU-mappable copy of z_ref
      allocate (this%vcoord_z_ref(0:nz), source=0.0_wp)
      if (allocated(this%vcoord%z_ref)) then
         this%vcoord_z_ref = this%vcoord%z_ref
      end if

      ! Full z* per-cell z_ref. Allocated here; per-cell contents filled
      ! from bathymetry by `ml_build_vcoord_zstar_full_unstr` later.
      if (this%vcoord%coord_type == VCOORD_ZSTAR_FULL) then
         allocate (this%vcoord_z_ref_col(0:nz, n), source=0.0_wp)
      end if

      allocate (this%h_layer(nz, n), source=0.0_wp)
      allocate (this%hu_layer(nz, n), source=0.0_wp)
      allocate (this%hv_layer(nz, n), source=0.0_wp)
      allocate (this%rho_layer(nz, n), source=0.0_wp)
      allocate (this%w_interface(nz + 1, n), source=0.0_wp)

      allocate (this%h_layer0(nz, n), source=0.0_wp)
      allocate (this%hu_layer0(nz, n), source=0.0_wp)
      allocate (this%hv_layer0(nz, n), source=0.0_wp)

      ! KPP per-cell boundary-layer workspace (always allocated; tiny).
      allocate (this%kpp_bl_depth(n), source=0.0_wp)
      allocate (this%kpp_u_star(n), source=0.0_wp)
      allocate (this%kpp_buoy_flux(n), source=0.0_wp)

      ! k-epsilon interface state — allocated unconditionally (mirrors KPP).
      ! Seeds mirror the structured path.
      allocate (this%keps_tke(nz + 1, n), source=1.0e-6_wp)
      allocate (this%keps_eps(nz + 1, n), source=1.0e-8_wp)
      allocate (this%keps_nut(nz + 1, n), source=1.0e-4_wp)
      allocate (this%keps_kt(nz + 1, n), source=1.0e-4_wp)
      allocate (this%keps_ustar_s(n), source=0.0_wp)
      allocate (this%keps_ustar_b(n), source=0.0_wp)
      ! Picard save slots: only needed when keps_si_picard > 1.  The caller
      ! (rki_driver / test) sets keps_si_picard BEFORE calling init, so the
      ! allocation can happen here.
      if (this%vmix_use_keps .and. this%keps_si_picard > 1) then
         allocate (this%keps_tke0(nz + 1, n), source=0.0_wp)
         allocate (this%keps_eps0(nz + 1, n), source=0.0_wp)
      end if

      allocate (this%tracer_advect_scratch(nz, n), source=0.0_wp)
      allocate (this%ml_lgrad_x(nz, n), source=0.0_wp)
      allocate (this%ml_lgrad_y(nz, n), source=0.0_wp)
      allocate (this%ml_hdiff_hTr_pre(nz, n), source=0.0_wp)
      allocate (this%ml_horvisc_hu_pre(nz, n), source=0.0_wp)
      allocate (this%ml_horvisc_hv_pre(nz, n), source=0.0_wp)

      allocate (this%xmom_pre(n), source=0.0_wp)
      allocate (this%ymom_pre(n), source=0.0_wp)

      if (this%vcoord%needs_remap()) then
         allocate (this%remap_dz_old(nz, n), source=0.0_wp)
      end if

      ! Reference-density profile: 100 piecewise-linear samples.
      this%nz_ref = 100
      allocate (this%rho_ref(0:this%nz_ref), source=0.0_wp)
      this%rho_ref_dz = 0.0_wp

      ! Default tracer registry: salinity at index 1, temperature at index 2.
      allocate (this%tracers(2))
      this%idx_salinity = 1
      this%idx_temperature = 2
      call this%tracers(this%idx_salinity)%init(mesh, nz)
      call this%tracers(this%idx_temperature)%init(mesh, nz)
      this%tracers(this%idx_salinity)%name = "salinity"
      this%tracers(this%idx_salinity)%long_name = "Sea water salinity"
      this%tracers(this%idx_salinity)%units = "PSU"
      this%tracers(this%idx_salinity)%standard_name = "sea_water_salinity"
      this%tracers(this%idx_temperature)%name = "temperature"
      this%tracers(this%idx_temperature)%long_name = "Sea water potential temperature"
      this%tracers(this%idx_temperature)%units = "degC"
      this%tracers(this%idx_temperature)%standard_name = "sea_water_potential_temperature"
   end subroutine unstr_multilayer_init

   subroutine unstr_multilayer_register_passive_tracer(this, mesh, name, tr_init, &
                                                       tr_min, tr_max, kappa_bg, &
                                                       hdiff_kappa)
      !! Append a passive tracer (eos_coeff = 0) to the registry. Must be
      !! called AFTER `init` and BEFORE `state_unstr_enter_data` (so the
      !! new tracer's `hTr` gets device-mapped) — enforced via
      !! `registry_locked`: once the registry is device-mapped the call
      !! logs an error and returns without touching it. S/T keep their
      !! indices.
      class(unstr_multilayer_state_t), intent(inout) :: this
      type(mesh_unstr_t), intent(in) :: mesh
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: tr_init, tr_min, tr_max, kappa_bg, hdiff_kappa

      type(unstr_tracer_t), allocatable :: tmp(:)
      integer :: n_old, it

      if (this%registry_locked) then
         call global_logger%error("unstr_multilayer_state: register_passive_tracer("// &
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
      call this%tracers(it)%init(mesh, this%nz_ml)
      this%tracers(it)%name = trim(name)
      this%tracers(it)%long_name = trim(name)
      this%tracers(it)%units = "1"
      this%tracers(it)%standard_name = ""
      this%tracers(it)%tr_init = tr_init
      this%tracers(it)%tr_inflow = tr_init
      this%tracers(it)%bc_ext = tr_init
      this%tracers(it)%tr_min = tr_min
      this%tracers(it)%tr_max = tr_max
      this%tracers(it)%kappa_bg = kappa_bg
      this%tracers(it)%hdiff_kappa = hdiff_kappa
      this%tracers(it)%eos_coeff = 0.0_wp
      this%tracers(it)%eos_ref = 0.0_wp
   end subroutine unstr_multilayer_register_passive_tracer

   subroutine unstr_multilayer_destroy(this)
      class(unstr_multilayer_state_t), intent(inout) :: this

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

      if (allocated(this%rho_layer)) deallocate (this%rho_layer)
      if (allocated(this%w_interface)) deallocate (this%w_interface)
      if (allocated(this%kpp_bl_depth)) deallocate (this%kpp_bl_depth)
      if (allocated(this%kpp_u_star)) deallocate (this%kpp_u_star)
      if (allocated(this%kpp_buoy_flux)) deallocate (this%kpp_buoy_flux)
      if (allocated(this%keps_tke)) deallocate (this%keps_tke)
      if (allocated(this%keps_eps)) deallocate (this%keps_eps)
      if (allocated(this%keps_nut)) deallocate (this%keps_nut)
      if (allocated(this%keps_kt)) deallocate (this%keps_kt)
      if (allocated(this%keps_ustar_s)) deallocate (this%keps_ustar_s)
      if (allocated(this%keps_ustar_b)) deallocate (this%keps_ustar_b)
      if (allocated(this%keps_tke0)) deallocate (this%keps_tke0)
      if (allocated(this%keps_eps0)) deallocate (this%keps_eps0)
      if (allocated(this%rho_ref)) deallocate (this%rho_ref)

      if (allocated(this%tracer_advect_scratch)) deallocate (this%tracer_advect_scratch)
      if (allocated(this%ml_lgrad_x)) deallocate (this%ml_lgrad_x)
      if (allocated(this%ml_lgrad_y)) deallocate (this%ml_lgrad_y)
      if (allocated(this%ml_hdiff_hTr_pre)) deallocate (this%ml_hdiff_hTr_pre)
      if (allocated(this%ml_horvisc_hu_pre)) deallocate (this%ml_horvisc_hu_pre)
      if (allocated(this%ml_horvisc_hv_pre)) deallocate (this%ml_horvisc_hv_pre)

      if (allocated(this%xmom_pre)) deallocate (this%xmom_pre)
      if (allocated(this%ymom_pre)) deallocate (this%ymom_pre)
      if (allocated(this%remap_dz_old)) deallocate (this%remap_dz_old)

      this%nz_ml = 1
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
   end subroutine unstr_multilayer_destroy

   subroutine unstr_multilayer_enter_data(this)
      !! Map the multilayer allocatables onto the device. Parent
      !! `state%multilayer` rides over with the parent `copyin(state)`;
      !! only the inner allocatables need explicit attaches. The tracer
      !! registry's two-step attach is handled by `tracers_enter_data`
      !! (in `rki_state_unstr`).
      class(unstr_multilayer_state_t), intent(inout) :: this

      !$omp target enter data map(to: this%dsig_ml, &
      !$omp&                  this%h_layer, this%hu_layer, this%hv_layer, &
      !$omp&                  this%rho_layer, this%w_interface, &
      !$omp&                  this%h_layer0, this%hu_layer0, this%hv_layer0, &
      !$omp&                  this%rho_ref, &
      !$omp&                  this%kpp_bl_depth, this%kpp_u_star, this%kpp_buoy_flux)
      !$omp target enter data map(alloc: this%xmom_pre, this%ymom_pre, &
      !$omp&                  this%tracer_advect_scratch, &
      !$omp&                  this%ml_lgrad_x, this%ml_lgrad_y, &
      !$omp&                  this%ml_hdiff_hTr_pre, &
      !$omp&                  this%ml_horvisc_hu_pre, this%ml_horvisc_hv_pre)
      !$omp target enter data map(to: this%keps_tke, this%keps_eps, this%keps_nut, &
      !$omp&                  this%keps_kt, this%keps_ustar_s, this%keps_ustar_b)
      if (allocated(this%keps_tke0)) then
         !$omp target enter data map(to: this%keps_tke0, this%keps_eps0)
      end if
      if (allocated(this%vcoord_z_ref)) then
         !$omp target enter data map(to: this%vcoord_z_ref)
      end if
      if (allocated(this%vcoord_z_ref_col)) then
         !$omp target enter data map(to: this%vcoord_z_ref_col)
      end if
      if (allocated(this%remap_dz_old)) then
         !$omp target enter data map(alloc: this%remap_dz_old)
      end if
   end subroutine unstr_multilayer_enter_data

   subroutine unstr_multilayer_exit_data(this)
      class(unstr_multilayer_state_t), intent(inout) :: this

      !$omp target exit data map(from: this%h_layer, this%hu_layer, this%hv_layer, &
      !$omp&                  this%rho_layer, this%w_interface)
      !$omp target exit data map(delete: this%dsig_ml, &
      !$omp&                 this%h_layer0, this%hu_layer0, this%hv_layer0, &
      !$omp&                 this%xmom_pre, this%ymom_pre, &
      !$omp&                 this%tracer_advect_scratch, &
      !$omp&                 this%ml_lgrad_x, this%ml_lgrad_y, &
      !$omp&                 this%ml_hdiff_hTr_pre, &
      !$omp&                 this%ml_horvisc_hu_pre, this%ml_horvisc_hv_pre, &
      !$omp&                 this%rho_ref, &
      !$omp&                 this%kpp_bl_depth, this%kpp_u_star, this%kpp_buoy_flux)
      !$omp target exit data map(delete: this%keps_tke, this%keps_eps, this%keps_nut, &
      !$omp&                 this%keps_kt, this%keps_ustar_s, this%keps_ustar_b)
      if (allocated(this%keps_tke0)) then
         !$omp target exit data map(delete: this%keps_tke0, this%keps_eps0)
      end if
      if (allocated(this%vcoord_z_ref_col)) then
         !$omp target exit data map(delete: this%vcoord_z_ref_col)
      end if
      if (allocated(this%vcoord_z_ref)) then
         !$omp target exit data map(delete: this%vcoord_z_ref)
      end if
      if (allocated(this%remap_dz_old)) then
         !$omp target exit data map(delete: this%remap_dz_old)
      end if
   end subroutine unstr_multilayer_exit_data

end module rki_multilayer_state_unstr
