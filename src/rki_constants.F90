!! Physical constants and working precision for the tidal solver
module rki_constants
   !! Working precision kind parameter and physical constants.
   use pic_types, only: dp, sp
   implicit none
   private

   public :: wp
   public :: GRAVITY, RHO_WATER, OMEGA_EARTH
   public :: CP_WATER
   public :: LATENT_HEAT_VAPORIZATION
   public :: DRY_TOLERANCE
   public :: NH_MIN_DEPTH
   public :: THIN_LAYER_THRESHOLD
   public :: DESING_EPS
   public :: H_VANISHED, H_DIV_EPS
   public :: LAND_DEPTH_THRESHOLD
   public :: FROUDE_BACKSTOP
   public :: PI, TWO_PI, DEG2RAD, RAD2DEG
   public :: NZ_STACK_MAX
   public :: BPG_NZ_MAX
   public :: nz_stack_required, nz_stack_is_sufficient
   public :: VCOORD_LAGRANGIAN, VCOORD_EULERIAN_Z, &
             VCOORD_SIGMA, VCOORD_ZSIGMA, VCOORD_ZSTAR, &
             VCOORD_ZSTAR_FULL, VCOORD_ZSTAR_SIGMA, VCOORD_Z_FIXED, &
             VCOORD_RHO, VCOORD_HYCOM
   public :: REMAP_PCM, REMAP_PLM, REMAP_PPM, REMAP_PPM_H4, REMAP_PQM
   public :: KEPS_STAB_CONSTANT, KEPS_STAB_GALPERIN, KEPS_STAB_CANUTO

#ifdef RAKALI_DOUBLE_PRECISION
   integer, parameter :: wp = dp
      !! Working precision kind parameter (double)
#else
   integer, parameter :: wp = sp
      !! Working precision kind parameter (single)
#endif

   real(wp), parameter :: GRAVITY = 9.80665_wp
      !! Gravitational acceleration (m/s^2)
   real(wp), parameter :: RHO_WATER = 1035.0_wp
      !! Seawater density (kg/m^3)
   real(wp), parameter :: OMEGA_EARTH = 7.2921e-5_wp
      !! Earth angular velocity (rad/s)
   real(wp), parameter :: CP_WATER = 4000.0_wp
      !! Specific heat capacity of seawater (J/(kg·K)).  Converts surface
      !! heat flux Q_net [W/m²] to kinematic heat flux Q_net/(rho_0·cp) [m·K/s].
   real(wp), parameter :: LATENT_HEAT_VAPORIZATION = 2.5e6_wp
      !! Latent heat of vaporization of water (J/kg).  A representative
      !! near-surface value (2.501e6 at 0 °C, weakly T-dependent).  Used by
      !! the evaporative latent-heat surface cooling: a net evaporation of
      !! E [m/s] withdraws a heat flux ρ_water·E·L_v [W/m²] from the top of
      !! the water column as the evaporated mass changes phase.
   real(wp), parameter :: DRY_TOLERANCE = 1.0e-6_wp
      !! Minimum depth to consider a cell wet (m)
   real(wp), parameter :: NH_MIN_DEPTH = 1.0e-1_wp
      !! Minimum depth for non-hydrostatic pressure operations (m).
      !! Cells shallower than this are treated as hydrostatic to avoid
      !! 1/h singularities in the NH Poisson operator.
   real(wp), parameter :: THIN_LAYER_THRESHOLD = 1.0e-3_wp
      !! Thin-LAYER depth scale (m) for the MULTILAYER per-layer
      !! stabilisers (skip/zero thresholds in friction, horizontal
      !! viscosity, tracer diffusion, layer distribute/truncate).  The 2D
      !! barotropic update/protect kernels no longer hard-zero momentum at
      !! this depth — they desingularize at `DESING_EPS` instead (the old
      !! hard zero froze wet/dry fronts; see `DESING_EPS`).

   ! ---- Thin-layer constants of record (ocean-path semantic split, D4 taxonomy) ----
   real(wp), parameter :: H_VANISHED = 1.5e-4_wp
      !! Dynamic-vanish threshold (m).  A layer thinner than this is
      !! dynamically VANISHED — skip / merge into a neighbour, do NOT clamp.
      !! When dividing by `h_layer` for a PHYSICAL result, gate on this.
      !! Same value/role as coastal `VANISHING_LAYER_TOL`(_UNSTR).
   real(wp), parameter :: DESING_EPS = H_VANISHED
      !! Kurganov & Petrova (2007) desingularization depth `eps` (m) for
      !! the 2D barotropic update/protect kernels AND the CFL timestep
      !! estimator: the velocity
      !!     u = sqrt(2) h hu / sqrt(h^4 + max(h^4, eps^4))
      !! equals hu/h exactly for h >= eps and rolls off quadratically
      !! below, and momentum is reconstructed as hu := h*u.  Momentum in
      !! cells thinner than eps is smoothly suppressed, NOT zeroed — the
      !! pre-2026-08 hard zero below THIN_LAYER_THRESHOLD (together with
      !! `FROUDE_CAP = 10`, see FROUDE_BACKSTOP) froze thin water at
      !! wet/dry fronts: Ritter dry-bed tip at 1.49 c0 instead of 2 c0,
      !! beach films that never drained.  One value shared by the update
      !! kernels and the dt estimator, so both agree about what "thin"
      !! means.  (Earlier the estimator alone used
      !! `eps = max(1e-3, 0.01*min(dx,dy))`; the grid term is gone — on
      !! coarse grids it made the estimator blind to real velocities in
      !! cells shallower than 1% of the cell width, a latent CFL bug.)
      !!
      !! Value: tied to `H_VANISHED` — the codebase's established
      !! "dynamically vanished" depth below which a column has no
      !! meaningful velocity.  Measured dose-response on the Ritter
      !! dry-bed tip (dx = 0.25 m, exact celerity 2 c0, no-clamp ceiling
      !! 1.888 c0 at this resolution; design sweep at FROUDE_BACKSTOP=100):
      !!     eps = 1e-3    -> 1.617 c0
      !!     eps = 1.5e-4  -> ~1.73 c0   (this value; 1.736 in the final
      !!                                  backstop=50 configuration)
      !!     eps = 3e-5    -> 1.776 c0
      !! Diminishing returns below H_VANISHED, and the front position at
      !! a RESOLVABLE depth (1 cm) is best at 1.5e-4 (lag 0.36 m of
      !! 62.6 m travelled, 0.6%).
   real(wp), parameter :: H_DIV_EPS = 1.0e-20_wp
      !! Pure division-safety epsilon (m), far below any physical thickness:
      !! adding it (`1/(h+H_DIV_EPS)`) only prevents a 1/0, never changes a
      !! well-posed result.  Use where a divisor is ALREADY guaranteed
      !! positive (MOM6 `H_subroundoff` role).
   real(wp), parameter :: LAND_DEPTH_THRESHOLD = 2.0_wp
      !! Cell bathymetry depth (m, positive-down) below which an ocean cell
      !! is treated as LAND for surface-forcing masking.  Consumed by
      !! `ocean_state_seed_from_cfg` to populate `multilayer%wet_mask`.
   real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
      !! Pi
   real(wp), parameter :: TWO_PI = 8.0_wp*atan(1.0_wp)
      !! 2*Pi
   real(wp), parameter :: DEG2RAD = PI/180.0_wp
      !! Degrees → radians conversion factor.
   real(wp), parameter :: RAD2DEG = 180.0_wp/PI
      !! Radians → degrees conversion factor.
   real(wp), parameter :: FROUDE_BACKSTOP = 50.0_wp
      !! Last-resort Froude bound at wet/dry fronts: momentum is rescaled
      !! so |u| <= FROUDE_BACKSTOP*sqrt(g*h).  This is a blowup ARRESTER,
      !! not physics shaping — it replaced `FROUDE_CAP = 10` (2026-08),
      !! which bound exactly where the physics lives: the exact Ritter
      !! dry-bed tongue has Fr = 2(sqrt(h_l/h) - 1), unbounded at the tip,
      !! and the leading wet cell measured at exactly Fr = 10 (tip
      !! celerity 1.49 c0 instead of 2 c0).  A Fr-shaped cap that binds a
      !! draining film also freezes it structurally: saturated cells obey
      !! hu = Fr*sqrt(g)*h^1.5, so a uniform film has zero flux divergence
      !! and dh/dt = 0 (measured: ~40 consecutive beach cells all at
      !! u = -10*sqrt(g h), shoreline pinned).
      !!
      !! Why it exists at all: without it, a thin draining film on a
      !! sloping bed develops a genuine velocity runaway ABOVE the
      !! desingularization depth, where KP07 is an identity — measured
      !! u -> 4e3 m/s at h ~ 3.5e-3 m and a ~2000x CFL dt collapse on the
      !! Carrier-Greenspan beach fixture.
      !!
      !! Why 50 (all measured at DESING_EPS = 1.5e-4, dx = 0.25 m; the
      !! 200 row is from the design sweep, where it matched 100 exactly):
      !!     backstop =  50 -> Ritter tip 1.736 c0, beach max|u| ~10.6 m/s
      !!     backstop = 100 -> Ritter tip 1.744 c0, beach max|u| ~17.8 m/s
      !!     backstop = 200 -> Ritter tip = the 100 value, max|u| ~31.5 m/s
      !! The dry-front tongue is insensitive to the value (the KP07
      !! desingularization governs it), but the Carrier-Greenspan runup
      !! DISCRIMINATES: at backstop 100 the freed runup jet overshoots and
      !! strands a cm-scale wedge on the beach through rundown (excursion
      !! at a 1 cm detector 11.8 m vs exact 20 m); at 50 the excursion is
      !! 18.9 m (5.6% of exact) with the rundown edge at x = 83.9 m vs the
      !! exact 80 m.  50 is also the lowest measured value that leaves the
      !! Ritter tongue unchanged; going much lower re-approaches the
      !! measured cap-10 film-freeze pathology (Fr-shaped caps bind at
      !! drainage velocities once Fr*sqrt(g h) ~ a few m/s at cm depths).
      !! NaN safety: the rescale triggers only on a TRUE comparison, so a
      !! NaN speed passes through un-"corrected" and fails loudly instead
      !! of being laundered into a plausible extreme value (see CLAUDE.md
      !! on NaN-blind clamps).

#ifndef RKI_NZ_STACK_MAX
#define RKI_NZ_STACK_MAX 128
#endif
   integer, parameter :: NZ_STACK_MAX = RKI_NZ_STACK_MAX
      !! Maximum vertical layers for per-column stack workspace inside
      !! `do concurrent`.  Stack arrays in the BPG, remap, kappa-shear,
      !! Redi and diag-remap kernels are dimensioned by this constant so
      !! the compiler can emit fixed-size thread-local storage on GPU.
      !!
      !! **The requirement is `NZ_STACK_MAX >= nz + 1`** — see
      !! `nz_stack_required` below for the derivation and the survey of
      !! every consumer.  Enforced fail-loud by `validate_config`.
      !!
      !! History: 100 -> 256 (2026-06-19), on the belief that the Redi
      !! neutral-diffusion sweep needed `2*nz+2` *of this constant*.  It
      !! does not: Redi's `nsurf` locals are declared `2*NZ_STACK_MAX+2`
      !! (`rki_ocean_redi.F90` `PoLc`/`PoRc`/`KoLc`/`KoRc`/`hEc`), so they
      !! scale with the constant and only need `NZ_STACK_MAX >= nz`.
      !! Lowered 256 -> 128 (2026-08-05) after that measurement: 256 was
      !! ~2x larger than any consumer needed and cost ~4 GB of CUDA
      !! per-thread local memory, which the driver reserves at each
      !! kernel's first launch and never returns.  128 covers the largest
      !! shipped case (nz=90, `acc_channel_kitchensink_xl.nml`) and the
      !! largest test (nz=100) with headroom.
      !!
      !! Heap-allocated workspace arrays (PPM snapshots, remap snapshots)
      !! have no such limit — they are allocated at runtime.

   integer, parameter :: BPG_NZ_MAX = NZ_STACK_MAX
      !! Backward-compatible alias for NZ_STACK_MAX.

   ! ---- Vertical coordinate type constants ----

   integer, parameter :: VCOORD_LAGRANGIAN = -1
      !! Pure Lagrangian / isopycnal — `h_layer` evolves freely and is
      !! never remapped.  `compute_target_h` returns the current `h`
      !! unchanged and the ALE remap step is a no-op.  Distinct from
      !! `VCOORD_EULERIAN_Z` (which targets the static `H · dsig`
      !! decomposition): LAGRANGIAN tracks material surfaces, so for
      !! adiabatic flow it reproduces MOM6's `COORD_CONFIG="gprime"`
      !! NK=2 isopycnal layering.  Parsed from "lagrangian" / "isopycnal".
   integer, parameter :: VCOORD_EULERIAN_Z = 0
      !! Ignore SSH when forming `target_h`; layer thicknesses follow
      !! the static `H · dsig(k)` decomposition.  Used by the ocean
      !! path as the "do nothing, leave the IC layers alone" default
      !! before the user picks a remap-aware coord_type.  Was originally
      !! declared in `rki_ocean_vcoord` as an orphan; moved here so the
      !! VCOORD enum lives in one place.
   integer, parameter :: VCOORD_SIGMA = 1
      !! Pure terrain-following sigma coordinate (default)
   integer, parameter :: VCOORD_ZSIGMA = 2
      !! z-sigma hybrid: z-levels in deep/steep, sigma in shallow
   integer, parameter :: VCOORD_ZSTAR = 4
      !! z-star "lite" (MOM6-like quasi-horizontal): single global
      !! `z_ref` pattern stretched proportionally per column by H/z_ref(nz)
   integer, parameter :: VCOORD_ZSTAR_FULL = 5
      !! z-star "full" MOM6: per-column `z_ref` anchored to local
      !! bathymetry with vanishing layers below the bed.  Surface
      !! layers stay at a fixed physical thickness independent of H
   integer, parameter :: VCOORD_ZSTAR_SIGMA = 6
      !! z*/sigma hybrid: smoothstep blend from pure sigma in shallow
      !! water (H ≤ depth_transition) to z*-lite in deep water
      !! (H ≥ depth_transition + blend_width).  Coastal wet/dry
      !! friendliness in shallow zones, SSH-tracking in deep zones
   integer, parameter :: VCOORD_Z_FIXED = 7
      !! Fixed-z interfaces with vanishing layers in shallow water.
      !! ALE remap pulls `h_layer(k) → h_target(k)` where the targets
      !! come from absolute-depth interfaces `z_target = (0, h_nominal,
      !! 2·h_nominal, ..., h_ref)` with `h_nominal = h_ref/nz_ml`.
      !! In deep cells the interfaces stay locked at those depths; in
      !! shallow cells the bed-side layers vanish to a minimum
      !! thickness `h_min` while the surface layer absorbs the residual.
      !! Mirrors MOM6's `COORD_CONFIG = "gprime"` layer structure and
      !! the existing `seed_h_layer_z_fixed_impl` IC algorithm — used
      !! as a per-step ALE target rather than an init-only seed so the
      !! interface stays anchored as the simulation evolves.
   integer, parameter :: VCOORD_RHO = 8
      !! Isopycnal coordinate (P2).  Layer interfaces are placed on
      !! prescribed potential-density surfaces `rho_target(0:nz)`
      !! (monotone-increasing, referenced to `rho_ref_pressure`).  The
      !! ALE regrid inverts a PPM reconstruction of the column density
      !! profile for the depth where ρ equals each interior target, then
      !! remaps T/S/tracers/velocities onto the new grid.  Bottom-up
      !! convention: `rho_target(0)` (lightest) anchors the surface
      !! interface (k=nz); `rho_target(nz)` (densest) the bed (k=1).
      !! Validation-grade alone (weakly-stratified columns collapse);
      !! HYCOM hybrid is the production follow-on (P3).  See
      !! `rki_ocean_vcoord :: ocean_vcoord_compute_target_h_rho_impl`.
   integer, parameter :: VCOORD_HYCOM = 9
      !! Hybrid z*/isopycnal coordinate (P3, Bleck 2002 / MOM6
      !! `coord_hycom`).  Runs the exact `VCOORD_RHO` density-space
      !! inversion (same `ocean_vcoord_compute_target_h_rho` kernel,
      !! `hybrid=.true.`) but adds two HYCOM deltas around it: a
      !! bottom-up density monotonize before the inversion, and a z*
      !! nominal-floor sweep after it that pushes too-shallow interfaces
      !! down to a fixed-resolution surface band (`dsig·(H+η)`).  Net
      !! effect: a fixed-res near-surface z* band (no surface collapse)
      !! with an isopycnal interior — the production GVC coordinate.
      !! Reuses `rho_target` / `rho_ref_pressure`; no new knobs.

   ! ---- Vertical remapping method constants ----

   integer, parameter :: REMAP_PCM = 1
      !! Piecewise constant (0th order, donor cell)
   integer, parameter :: REMAP_PLM = 2
      !! Piecewise linear with minmod limiter (1st order)
   integer, parameter :: REMAP_PPM = 3
      !! Piecewise parabolic (2nd order, Colella & Woodward)
   integer, parameter :: REMAP_PPM_H4 = 4
      !! Piecewise parabolic with non-uniform 4th-order (H4) edge values —
      !! thickness-weighted edge estimate (White & Adcroft 2008) feeding the
      !! same Colella & Woodward parabola + monotonicity limiter as REMAP_PPM
   integer, parameter :: REMAP_PQM = 5
      !! Piecewise quartic (4th-order, PQM_IH4IH3 — White & Adcroft 2008):
      !! implicit-h4 edge values + implicit-h3 edge slopes + a quartic per-cell
      !! reconstruction with the White & Adcroft monotonicity limiter and a
      !! conservative quartic overlap integral.  Requires nz >= 5; falls back to
      !! REMAP_PPM for thinner columns.  Boundary cells are PCM.

   ! ---- k-epsilon stability-function scheme constants ----

   integer, parameter :: KEPS_STAB_CONSTANT = 0
      !! Constant c_mu (=0.09) / Pr_t (=1) -- standard high-Re k-epsilon
      !! (default; bit-identical to the original closure).
   integer, parameter :: KEPS_STAB_GALPERIN = 1
      !! Galperin et al. (1988) quasi-equilibrium stratification-aware
      !! stability functions (Mellor-Yamada level-2.5 algebraic form).
   integer, parameter :: KEPS_STAB_CANUTO = 2
      !! Canuto et al. (2001) Model A second-moment-closure stability functions
      !! (rational alpha_M / alpha_N form per Umlauf & Burchard 2003).

contains

   pure function nz_stack_required(nz) result(req)
      !! Smallest `NZ_STACK_MAX` that safely covers a run of `nz` layers.
      !!
      !! `nz + 1`.  Derived from a survey of every `NZ_STACK_MAX`-dimensioned
      !! array in the tree; the binding consumers are the ones that index a
      !! plainly-`NZ_STACK_MAX`-declared array to `nz+1` rather than the
      !! `NZ_STACK_MAX+1` an interface array would use:
      !!
      !!   * `rki_ml_keps` `S2`/`N2`/`dz_if`/`a`/`b`/`c`/`rhs`/`cp` — written
      !!     at `ni = nz + 1`.
      !!   * `rki_si_unstr_closures` `kapS`/`kapT` —
      !!     written at `nz + 1`.
      !!   * `rki_ml_tracers` `ml_settling_tracer` — the augmented settling
      !!     column is `nz + 1` deep, and it calls `remap_column_ppm(nz+1, …)`,
      !!     so the remap kernel's own `z_old(0:NZ_STACK_MAX)` / `q_L`/`q_R`/`q6`
      !!     must also hold `nz + 1`.
      !!   * `rki_ocean_diag_fills` `z_iface` — needs `n_bin + 2`, which the
      !!     `n_bin <= NZ_STACK_MAX - 1` output-level guard bounds by `nz + 1`.
      !!
      !! NOT `2*nz + 2`: the Redi neutral-surface locals that motivated the
      !! 256 bump are declared `2*NZ_STACK_MAX + 2` and so scale with the
      !! constant — they need only `NZ_STACK_MAX >= nz`.
      integer, intent(in) :: nz
      integer :: req

      req = nz + 1
   end function nz_stack_required

   pure function nz_stack_is_sufficient(nz) result(ok)
      !! `.true.` when the COMPILED `NZ_STACK_MAX` covers `nz` layers.
      !! Callers that get `.false.` must refuse the run — the overflow is
      !! a silent thread-local-storage overrun, not a crash.
      integer, intent(in) :: nz
      logical :: ok

      ok = nz_stack_required(nz) <= NZ_STACK_MAX
   end function nz_stack_is_sufficient

end module rki_constants
