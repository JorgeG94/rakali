# Rakali — Capabilities and Limitations

A snapshot of what the solver can and can't do today. Sourced from `ROADMAP_COASTAL.md`, `CLAUDE.md`, and the current state of the codebase. Every "can't" item below is something the codebase doesn't do yet — not a hard "won't"; consult the roadmap for project-scope tags and validation criteria for the planned items.

For the full feature checklist (planned, validated, in flight, deferred), see `ROADMAP_COASTAL.md` (coastal + cross-cutting), or `v0.1.0_checklist.md` for near-term. For namelist / config keys + physics, see `REFERENCE.md`. For acronyms (HLL, KNP, CWC, BPG, PP81, MUSCL, ALE, …) and concept definitions (vanishing layer, conservative soft clamp, mode split, …), see `codebase/CONCEPTS.md`.

Rakali is `sim_type='coastal'` — production-mature. HLL/HLLC (structured) + KNP
central-upwind (unstructured) on an A-grid; 2D barotropic, 3D hydrostatic
multilayer (both backends), structured non-hydrostatic. Used in tidal / coastal /
estuary / harbour / dam-break / inundation work.

---

## Time-step constraint — read this first

> **The default path is CFL-limited by surface gravity waves.** The free surface is treated **explicitly** (HLL/KNP Riemann on the barotropic step), so the global timestep is bounded by `Δt < dx / sqrt(g·H)`. For coastal applications (`H ~ 10–100 m`, `dx ~ 10–1000 m`), this typically means `Δt` of **seconds**, not minutes-to-hours. The opt-in `time_scheme="semi_implicit"` mode (next subsection) removes this constraint within its documented scope.

Several individual operators *are* unconditionally stable on their sub-problem (Manning friction is semi-implicit; Coriolis is the exact closed-form rotation), but **none of those relax the global `Δt`** on the explicit path. To unlock large `Δt` in the SCHISM / MOM-split sense you need a semi-implicit free-surface formulation that treats the gravity-wave term implicitly — Rakali now ships this as an **opt-in mode**: `&time_nml time_scheme="semi_implicit"` (**unstructured grids only**: 2D barotropic + 3D sigma multilayer — TPFA on triangle edges with circumcentre η + area-weighted CG, P3.1/P3.2 — see "Time integration and stepping" below and `REFERENCE.md` for the honest limits). The default explicit path is untouched.

If your use case needs `Δt` set by accuracy rather than gravity-wave CFL (basin-scale tides at coarse resolution, multi-year hindcasts) **and** needs physics beyond the SI mode's current scope (non-Dirichlet open boundaries, remap vcoords under SI, or momentum advection inside the implicit step), Rakali in its current form is still the wrong tool — that widening is the SI post-P3 roadmap.

### What this constraint *buys us*

The flip side: by being explicit on the gravity-wave term, Rakali resolves several things that semi-implicit free-surface models (SCHISM, MOM-split-fast-mode) have to filter, smear, or parameterise:

- **Surface and external gravity waves** are propagated faithfully, with no implicit damping or phase error from a Crank-Nicholson wave-equation solve. Tsunami fronts, storm surge propagation, seiches, dam-break / outburst-flood-like transients all carry their high-frequency content through the solver intact.
- **Shocks and hydraulic jumps** are resolved by the well-balanced HLL/KNP Godunov machinery — sharp wet/dry fronts, tidal bores, supercritical-to-subcritical transitions. Semi-implicit ELM-based schemes diffuse these.
- **Lake-at-rest is preserved exactly** under the well-balanced Audusse hydrostatic reconstruction, so quiescent variable-bathymetry initial states don't generate spurious flows. (Holds at the discretisation level, not just to truncation order.)
- **Inertial and internal-wave time scales** are temporally well-resolved without subcycling tricks, because `Δt` is already small enough to resolve them many-points-per-period.
- **Phase-resolving wave dynamics** in the non-hydrostatic path are physically possible (latent capability per ROADMAP_COASTAL.md §12A) — a semi-implicit free surface fundamentally cannot host them.
- **High-frequency forcing coupling** (rapid wind shifts, atmospheric pressure pulses, sub-tidal modulations) is reflected in the dynamics step-by-step rather than time-averaged into the implicit step.

So the choice is a **deliberate tradeoff**, not just a limitation: pay more wall-clock per simulated day, get faithful wave / shock physics in return. Which side wins depends on the application — tidal/coastal/inundation/tsunami-class problems benefit from the explicit treatment; multi-year basin hindcasts where the dynamics are quasi-steady at the implicit `Δt` benefit from SCHISM-style large `Δt`.

---

## Solver backends (coastal path)

**Structured (Cartesian):** 2D shallow water, 3D hydrostatic multilayer (sigma + z-sigma hybrid), and non-hydrostatic (CG Poisson) all available. HLL or HLLC Riemann (selectable via `use_hllc`) + minmod, SSP-RK2, mode-split barotropic-baroclinic.

**Unstructured (triangular):** 2D shallow water and 3D hydrostatic multilayer. Kurganov-Noelle-Petrova central-upwind flux, Green-Gauss gradients, MUSCL tracer advection. **No non-hydrostatic mode on unstructured.**

Both backends share the same explicit time integrator, vertical coordinate types, EOS, and tracer registry. The semi-implicit time scheme is unstructured-only.

---

## Physics (coastal path)

### Shipped
- **Hydrodynamics:** 2D barotropic on both grids; 3D mode-split multilayer on both grids; NH structured-only.
- **Vertical mixing:** Pacanowski-Philander 1981 Richardson-dependent mixing for momentum + tracers (both grids).
- **Horizontal mixing:** Smagorinsky viscosity, Laplacian tracer diffusion (both grids).
- **Bottom friction:** Manning, semi-implicit, on the bottom layer using layer velocity.
- **Coriolis:** closed-form rotation applied **per layer** on both grids — the rotation matrix is the *exact* solution of the inertial ODE so the kernel is unconditionally stable on its own sub-problem (the global `Δt` is still gravity-wave-CFL-limited; see the Time-step section above). NH mode-split path is still depth-averaged-only and pairs with the ALE refactor on the roadmap.
- **Baroclinic pressure gradient:** Jacobian (Shchepetkin-McWilliams) and sigma variants on both grids.
- **Tracers:** salinity + temperature in a registry (`tracers(:)`) on both grids; first-order upwind and second-order (PPM structured / MUSCL unstructured) advection.
- **EOS:** linear two-tracer `rho = rho_0 + beta_S·(S − S_ref) − alpha_T·(T − T_ref)` (default), plus opt-in SHARED nonlinear EOS (`&tracer_nml eos_form = "wright" | "roquet_spv"`, both grids) routing density through the shared `rki_eos` dispatch, evaluated at p = 0 (surface/potential density).
- **Surface tracer fluxes:** `heat_flux` (W/m²) and `evap_minus_precip` (m/s) applied as a top-layer source on both grids; programmatic setters via the C/Python API.
- **Evaporation as a real mass sink** (structured + unstructured multilayer, `evap_mass_sink`, default off = virtual salt flux): net evaporation removes water column-proportionally with salt held (concentration rises, salt mass exact) and sensible heat/momentum carried away with the water (T, u, v invariant); clamps at `ml_min_depth` on drying flats. **Latent-heat cooling** (`evap_latent_cooling`, default off, requires `evap_mass_sink`): the evaporated mass additionally withdraws its latent heat `ρ_water·E·L_v` [W/m²] from the surface layer (distinct from the sensible term — no double-count), cooling the column along `T_mean(t) = T_INIT − L·ln(H0/h(t))`. The two channels are normalized to slightly different departed masses by design — sensible to `ρ_0·dh_eff`, latent to `ρ_water·dh_eff` (matching the standard bulk `Q_lat = ρ_water·E·L_v`); the `ρ_water/ρ_0 ≈ 3.5%` gap is each term using the density natural to its own formulation, not an inconsistency. **Energy-closure caveat:** the whole latent sink lands in the top layer, so it cools `k=nz` by `latent_gain·dh_eff/h_top`; in a thin-surface-layer + strong-evaporation + weak-mixing regime `T(nz)` can hit the `T_min` clamp floor, where the non-conservative tracer clamp silently re-adds heat and the exact energy budget no longer holds (spreading the sink over a surface band is the deferred remedy). Surface-intensified (rather than column-uniform) evaporative salinity structure is a follow-up: it needs the corrected physical vertical velocity `w` (branch `fix/ml-vertical-velocity-w`).
- **Wet/dry:** hydrostatic reconstruction (Audusse); thin-cell momentum via Kurganov & Petrova (2007) desingularized velocity below `DESING_EPS = 1.5e-4 m` (smooth suppression, not zeroing — hard zero only below `DRY_TOLERANCE`) plus a `FROUDE_BACKSTOP = 50` blowup arrester. Known residual: beaches retain a film of up to `slope·dx` per cell (Audusse staircase retention), so mm-scale water strands on fine intertidal slopes; dry-front tips advance at ~1.74 c₀ (structured) / ~1.59 c₀ (unstructured KNP) vs the exact 2 c₀.
- **Velocity truncation:** CFL > 0.8 velocity clipping (both backends).

### Not yet shipped — see ROADMAP_COASTAL.md §1, §3, §4

- **GLS, k-ω / MY2.5** turbulence closures. Shipped: PP81 (default), KPP
  (`vmix_use_kpp`), and two-equation k-ε (`vmix_use_keps`,
  constant/Galperin/Canuto stability functions) — all on BOTH backends
  under the explicit time schemes, and under `semi_implicit` on the
  unstructured backend.
- **Implicit vertical diffusion** on the EXPLICIT path's PP81 route
  (explicit clamped flux caps `nu_max` there). The k-ε apply is already
  implicit (Thomas), and the `semi_implicit` mode solves ALL vertical
  viscosity/diffusion implicitly (PP81 / k-ε / KPP).
- ~~**Nonlinear EOS**~~ — SHIPPED: `&tracer_nml eos_form = "wright" |
  "roquet_spv"` routes density (both backends) through the shared
  `rki_eos` at p = 0. Still open: in-situ (pressure-dependent) evaluation
  and a TEOS-10 GSW branch.
- **Convective adjustment.** No equivalent yet — still only PP81's bounded
  `Ri<0` clip (~1.01e-2 m²/s, not a real convective diffusivity).
- **Double-diffusive mixing, shortwave penetration profile.** The solver
  supports `K_T != K_S` (`vdiff_apply_tracers` takes independent
  `kt_source` / `ks_source`) but no closure produces a nonzero split yet.
- **Biogeochemistry, vegetation drag.** No source terms beyond surface
  fluxes and single-class suspended sediment (structured) today.

### Unstructured-specific gaps vs structured (ROADMAP_COASTAL.md §4)
- CFL-based velocity truncation
- Per-layer bottom friction
- Fused multilayer kernel
- Non-hydrostatic mode

---

## Vertical coordinates

- **Sigma** (terrain-following) — primary, both backends.
- **Z-sigma hybrid** (`VCOORD_ZSIGMA`) — implemented and routed through the conservative remap.
- **Z-star (`VCOORD_ZSTAR`)** — z*-lite implementation (simplified): a global `z_ref` reference layer pattern stretched proportionally per column with the actual water depth `H = h_bed + η`, so `sum(dz) = H` exactly and the relative layer spacing is preserved as η changes. Eliminates the per-step layer-redistribution that distorts the BPG / EOS / vmix at large SSH amplitudes (storm surge, large tides). Goes through the conservative remap path (`ml_distribute_with_remap`).
- **Z-star/sigma hybrid (`VCOORD_ZSTAR_SIGMA`)** — smoothstep blend from pure sigma in shallow water (`H ≤ depth_transition`) to z*-lite in deep water (`H ≥ depth_transition + blend_width`). Coastal sigma's wet/dry friendliness in shoaling zones, z*'s SSH-tracking in deep zones. Both contributing branches independently conserve mass, so the blend conserves by construction. Reuses the same `depth_transition` / `blend_width` knobs as `VCOORD_ZSIGMA`.
- **Z-star full (`VCOORD_ZSTAR_FULL`)** — per-column `z_ref(0:nz, ...)` built from local bathymetry at init via `zstar_full_build_column` (uniform or log surface stretching). Surface layer anchored at `zstar_h_surf_target` regardless of total H — z*-lite cannot do this. Bed-side layers vanish (collapse to `zstar_h_min`) when `H < h_bed_ref`; the surface trim in `vcoord_target_dz_column_zstar_full` keeps `sum(dz) = H` exactly. Per-cell `vcoord_z_ref_col` is GPU-attached. Tracer transport under ZSTAR_FULL uses the **vanishing-layer-gated centred CWC** advection + the **conservative two-pass column-conservative tracer limiter** (both backends) so the per-layer flux split, the remap, and the bounds-enforcement all stay mass-conservative across deep-shallow transitions where vanishing bed layers would otherwise be over-allocated by a centred face flux. Opt-in via `vcoord_type = "zstar_full"`; ZSTAR-lite remains the default for benign-bathymetry coastal cases.
- **S-coordinate (ROMS)** — placeholder only.
- **Isopycnal, partial cells** — not implemented.

### Known limitations on ZSTAR_FULL with intertidal-flat geometry
On bathymetry with frequent wet/dry transitions (e.g. tidal estuary intertidal flats — Lakes Entrance, Sydney Harbour mudflats), ZSTAR_FULL still leaks ~1-2 % salt over a 4 h tidal cycle. The leak is now **wet/dry-driven** (cells crossing `h < DRY_TOLERANCE` zero out their `hTr`; on rewet the soft clamp clips up to `S_min`, creating spurious mass), not advection / remap / clamp-driven inside the column (those are conservative). The same effect hits sigma less severely because its layers are uniform per column. The four ranked options for closing this gap are h-reconstruction at wet/dry edge → tracer salvage at w→d transition → tracer reservoir → sub-cell partial wetting. For coastal-estuary / intertidal-flat cases the recommended vcoord remains sigma or zstar-lite; ZSTAR_FULL is best deployed on shelf-break / steep-bathymetry / global-scale runs where its surface-anchoring property pays for itself.

---

## Boundary conditions

### Shipped (structured)
WALL, OPEN (Flather), TIDAL (multi-constituent), NESTED (one-way from coarser parent), INFLOW, DISCHARGE, CLAMPED, SPONGE, CHAPMAN.

### Shipped (unstructured)
WALL, OPEN, TIDAL, INFLOW, DISCHARGE, CLAMPED, CHAPMAN, OUTFLOW, SPONGE, NESTED
(the last two as of 2026-07).

- `BC_SPONGE`: separate per-cell relaxation kernel (`apply_sponge_unstr`),
  called after each RK2-stage barotropic update. The sponge ghost falls through
  to the wall state in the flux kernel; the relaxation is the separate pass, per
  structured semantics. **v1 limitation**: the sponge band is one cell deep
  (cells with a direct `BC_SPONGE` boundary edge). Graded multi-cell width
  profiles are deferred — a triangular mesh has no regular cell-distance index
  the way the structured `(i,j)` grid does.
- `BC_NESTED`: **smoke-level** — routed through the Flather/`BC_TIDAL` algebra
  with a per-group prescribed external elevation. Live parent-grid coupling (a
  NetCDF reader + spatio-temporal interpolation onto the triangular boundary,
  as the structured `rki_nesting` does) is still deferred.

### Not in either backend
- Open-boundary T,S relaxation/nudging
- 2D tidal harmonics from a file (per-edge constants only)
- Tidal potential / SAL correction
- Wave radiation stress
- Hydraulic structures (weirs, culverts, gates)

---

## Forcing (coastal path)

### Shipped
- Wind stress (NetCDF, spatial + temporal)
- Atmospheric pressure (NetCDF)
- Precipitation (NetCDF)
- Heat flux + virtual freshwater flux (programmatic setters via C/Python API; both grids)
- River discharge at inflow BCs (constant per BC group)

### Not yet shipped
- **NetCDF-driven heat flux / E−P.** Requires the existing `forcing_file` pipeline extension (ROADMAP_COASTAL.md §10).
- **Bulk aerodynamic flux module** (Zeng / COARE). User must pre-compute `τ`, `Q_net`, `E−P` from atmospheric variables today — there's no in-model conversion from `(U10, T_air, q_air, SST, SLP)` to fluxes.
- **Time-varying river `Q(t)` from CSV/NetCDF.**

---

## Coupling and nesting

- **One-way structured nesting** (child reads parent NetCDF) — works.
- **Two-way nesting** — not implemented.
- **Structured ↔ unstructured nesting** — not implemented.
- **Telescoping (>2 levels)** — not implemented.
- **Local sub-stepping** — not implemented; smallest triangle dictates global `dt`.
- **External wave-model coupling** (SWAN/WW3) — not implemented.

---

## Numerics

### Shipped
- HLL and HLLC Riemann (structured, selectable via `use_hllc`), KNP central-upwind (unstructured)
- Minmod limiter (structured), Green-Gauss gradients (unstructured)
- SSP-RK2, CFL-adaptive timestep, mode-split subcycling
- PPM (structured) and MUSCL (unstructured) tracer advection — second-order, monotone

### Not yet shipped
- Roe with entropy fix
- WENO, MUSCL-Hancock
- Limiters beyond minmod (van Leer, van Albada, superbee)
- SSP-RK3 / RK4
- IMEX time stepping (needed once implicit vertical diffusion lands)
- Adaptive timestep on tracer CFL
- MPDATA / FCT for tracers

---

## I/O

### Shipped
- NetCDF output (CF-1.8)
- Optional deflate compression
- Async ring-buffered I/O
- Gauge stations (bilinear, both grids)
- Restart / warm start
- UGRID NetCDF for unstructured

### Not yet shipped
- Diagnostic outputs: vorticity, Okubo-Weiss, Ertel PV, energy budgets
- Lagrangian drifters
- Moving gauges
- Parallel collective NetCDF
- In-situ visualization

---

## Time integration and stepping

**Global `Δt` is set by surface-gravity-wave CFL on the explicit HLL barotropic step** — see the "Time-step constraint" call-out at the top of this document. The per-operator semi-implicit components below don't relax that constraint; the opt-in `time_scheme="semi_implicit"` mode (next subsection) does, within its Phase 0/1 scope.

### What's explicit (and CFL-binding)
- Barotropic free-surface dynamics — HLL Riemann, structured; KNP central-upwind, unstructured. **This is what sets `Δt`.**
- Tracer advection (PPM / MUSCL / upwind).
- Baroclinic pressure gradient (sigma + Jacobian variants).
- SSP-RK2 outer integration on both backends.
- Mode-split barotropic-baroclinic subcycling for the NH path (`split_ratio` in `nonhydrostatic_nml`).

### What's already semi-implicit / unconditionally stable on its sub-problem
- **Manning bottom friction** — semi-implicit `hu_new = hu_old / (1 + Δt · C_f · |u| / h)` (barotropic + per-layer).
- **Coriolis** — closed-form rotation matrix per layer; the *exact* analytical solution of the Coriolis ODE, so it carries no `f·Δt < 1` stability constraint.

These are unconditionally stable as kernels, but the surrounding solver still advances at the gravity-wave-CFL `Δt`, so they buy you accuracy / robustness, not larger time steps.

### Semi-implicit free surface (`time_scheme="semi_implicit"`, opt-in, unstructured grids)
Casulli-style implicit gravity wave: θ-implicit surface PGF (θ=0.6 fixed)
+ matrix-free warm-started CG on the SPD Helmholtz η-system
— **no gravity-wave CFL**; requires `dt_fixed > 0`,
validated at barotropic Courant 10–20+. Internally staggered: the
predictor / back-substitution transports live on triangle edges ⇒ **exact
mass conservation to round-off**; the A-grid cell momenta downstream code
sees are a step-end reconstruction. A per-step mass-consistency
diagnostic (`si%mass_res`, the interior CG residual sum — the only path
mass error can enter the η update) is exposed alongside `cg_iters`.
Unstructured grids only: `time_scheme="semi_implicit"` with a structured
grid is rejected at configure time.

TPFA edge fluxes with η at the triangle **circumcentre** (Casulli &
Walters 2000; the centroid breaks TPFA consistency), per-cell
double-visit edge transports, and an **area-weighted CG** (the FV
operator is self-adjoint only in the cell-area inner product). Guarded by
a **fail-loud mesh-quality gate** (max circumcentre-line deviation > 5°
or non-positive circumcentre distance aborts — clean Delaunay meshes
measure < 1e-2°).

CG knobs (`&si_nml`, all default to plain CG ⇒ byte-identical):
- `precond = "none" | "jacobi" | "mg"` — `jacobi` is a diagonal
  preconditioner (variable cell areas put an O(100-1000)× spread on the
  diagonal even at uniform depth). `mg` selects AGGREGATION
  (agglomeration) multigrid (`rki_si_unstr_mg`; v1: 2D barotropic only —
  mg + wetdry / multilayer aborts at enable): 1526 -> 16 CG iterations on
  the 480-triangle funnel.
- `rtol` (default `1e-10`) — the CG relative residual tolerance;
  relaxing to `1e-8`/`1e-6` cuts iterations with `mass_res` as the audit
  instrument for the resulting mass consistency.
Convergence is gated on the **unpreconditioned** relative residual in
all modes, so iteration counts stay directly comparable.

**2D scope (P3.1)**: wall + tidal/clamped Dirichlet-η boundary edges;
linear predictor (no momentum self-advection).

**3D multilayer scope (P3.2, shipped)**: `use_multilayer` + nz≥2, sigma.
The per-cell column-coupled step in `(nz, ncells)` order: implicit
backward-Euler vertical viscosity (PP81 by default supplies ν) +
**Manning-only** implicit bottom drag (no linear-drag knob on the
unstructured barotropic state), the effective face depth
`Ĥ = −Σ h_k R_k` + implicitly-mixed predictor transport `Ğ = Σ h_k G_k`
feeding the unchanged TPFA eta-solve, **anchor + increment** per-layer
edge recovery (the uniform-density inviscid nz=4 run reproduces the 2D
SI to 1.4e-10 m pointwise), dsig-proportional registry tracer advection
on the exact blended continuity edge flux (uniform-S tidal pumping
invariant to CG tolerance — verified linear in `rtol`), and surface
wind stress on the top layer. **Tidal/clamped Dirichlet-η open
boundaries** are supported on the 3D path too (zero-gradient baroclinic
velocity, zero-gradient outflow / clamped inflow tracers).

**SI vertical-mixing closures (Q1, shipped)**: the column solve
consumes ANY interface ν through the `nu_iface` slot
(`R = solve(−1)` is ν-independent), so `vmix_use_keps` and
`vmix_use_kpp` work under `semi_implicit` too. **k-ε at SI
cadence**: the explicit `ml_keps_update_unstr` is already a
per-column backward-Euler/Patankar integrator (implicit Thomas
diffusion, implicit dissipation sinks, positivity-clamped sources), so
it runs UNCHANGED once per SI step (dt = 120–480 s certified by the Q1
prototype — steady states dt-independent to round-off, k/ε ≥ floors
under shock-off decay; `local_archive/prototypes/si_closures/`); `keps_nut`
fills `nu_iface` (the implicit momentum apply IS the column
tridiagonal), `keps_kt` drives the implicit tracer diffusion. Optional
`keps_si_picard = 2` (default 1) re-runs the k/ε solve with the
refreshed ν_t production coefficient for transient-fidelity studies —
1 pass matches the shipped explicit cadence and is steady-state-exact
per the prototype. **KPP at SI cadence**: `ml_kpp_boundary_layer_unstr`
reused verbatim (diagnostic BL depth / u\*/ B₀), the LMD94
shape-function overlay MAX'd onto PP81 into `nu_iface` with **no
explicit-flux CFL clamp** (the implicit tridiagonal replaces the
coastal explicit clamped apply, which would strangle the overlay at SI
step sizes), S/T tracer diffusion implicit (Thomas) with the overlay
diffusivity + the non-local counter-gradient divergence (conservative;
extra passive tracers fall back to PP81 as on the explicit dispatch).
Closures off ⇒ the untouched PP81 branch (byte-identity by
construction).

**Phase-2 limits
(fail-loud at configure)**: `VCOORD_SIGMA` only, boundaries are
wall/tidal/clamped only (Flather/Chapman/sponge/nested/inflow/discharge
rejected on SI-ML), `nz ≥ 2`, no sediment, baroclinic PGF
explicit (internal
gravity wave stays CFL-limited — the standard split posture). See
`REFERENCE.md` and `docs/semi_implicit_freesurface_plan.md`.

**Wetting/drying (`&si_nml wetdry`, opt-in, shipped, 2D barotropic + 3D
sigma multilayer)**: Casulli-style subgrid wet/dry for the SI free
surface (Casulli 2009, *Int. J. Numer. Meth. Fluids* 60:391-408). The
Helmholtz cell-mass term - the constant `1.0` of the fully-wet
linearisation - is replaced by the exact wetted-volume function `V(eta)`
(per unit cell area), whose derivative is the per-cell wet fraction
`a(eta)` in `[0,1]`; this `a(eta)` becomes the SI operator's variable
diagonal. Because `V(eta)` is convex for open-air wet/dry, a **plain
Newton** outer loop is monotone with finite convergence - no line
search, no nested splitting (measured ~1.2-3.0 outer iterations per step;
hard cap 10, exceeding it emits a fail-loud warning). Face depths use the
frozen max-eta/sill rule (built from `eta^n`, held over the step), with
GATED-MEAN edge conveyance (all-wet reduces exactly to the linear edge
depth — measured all-wet agreement 1.8e-15 m); edges that go dry are
zeroed so a dead edge contributes no transport. The existing CG/Jacobi
machinery does the inner solve unchanged with the `1.0` swapped for
`a(eta)` (`rki_si_unstr_wetdry` + `state_unstr_enable_si(...,
wetdry=.true.)`). The wet/dry state is **stateless across restart**:
`b_sub`/`b_min`/`Hf`/the wet-mask are all reconstructed each step from
`h + b`, so there is nothing new to persist in a checkpoint. Under
wet/dry `si%cg_iters` is the SUM of the inner-CG iterations over the
Newton outer iterates on a step (not a single-solve count), and
`si%mass_res` reports the nonlinear residual sum (the true mass defect at
Newton convergence). **v1 envelope/limits**: the subgrid sub-bottom table
is live (`&si_nml wetdry_nsub` sub-bottoms per cell + `wetdry_bsub_fill`:
`"flat"` default = bit-identical to `nsub=1`, `"linear"` = the analytic
Casulli 2009 per-cell relief ramp — mean sub-bottom exactly `b`,
still-water volume unchanged, a coarse tidal flat drains/floods smoothly
at `1/n_sub` wet-fraction steps). `"linear"` × multilayer is **fenced
fail-loud** at configure and state-init (the η-solve conserves the
subgrid `v_cell(η)` but the ML thickness update consumes the plain bed —
mass would leak at the front; `flat`+ML and `linear`+2D stay allowed);
first-order-upwind tracer margin on the 3D path
(positivity-limited donor cell, not the windowed CW-PPM drain) — the thin
re-wetting front carries a documented ~2.3 PSU tracer overshoot in v1 (a
re-wet cell's `hTr` lags its thickness for one advection pass), confined
to the moving front; the well-wet interior stays uniform to ~0.05 PSU.
Default off (`wetdry=.false.`) leaves the SI operator's constant mass
term in place, so it is byte-identical. The quality gate is hardened for
this phase: `min_edge_delta` tracks the SIGNED circumcentre-line
projection on the outward normal, so a locally non-Delaunay flipped pair
is rejected.
Validated: Thacker-on-triangles (period 1.5%, shoreline excursion 230 m
vs exact 200 mesh-quantised, minV = 0 exactly, mass 1.3e-13 rel, Newton
≤ 4 outer), tidal-flat cycling with the water budget closed to 9e-13 rel
through dry/rewet cycles. ACCEPTANCE BENCH (funnel estuary, nz=4 S/T): SI
amplitude within 0.07–6.2% / phase within 3–4.5% of M2 vs the explicit
KNP; median-of-3 wall-clock speedup — **gfortran-CPU 19.5× (plain) /
48.7× (Jacobi, rtol=1e-8); V100-GPU 1.50× / 3.14× at 480 tris, shrinking
to 1.37×/1.73× at 7680 tris** (the no-MG triangle-CG tier is
launch-latency- and conditioning-bound on GPU; the shipped agglomeration
MG (`precond="mg"`) is the plateau-breaker for that tier — 1526 -> 16
iterations on the 480-triangle funnel).

### Not yet shipped
- **Semi-implicit free surface beyond the P3 campaign scope** — momentum advection in the SI predictor (needs a triangle departure-point walk; the earlier structured-only semi-Lagrangian option was removed along with the structured SI path), SI on structured grids (removed), remap vcoords (z*/zsigma) under SI, non-Dirichlet open boundaries (Flather/Chapman/sponge) on the 3D path.
- **Implicit vertical diffusion** (Thomas-algorithm tridiagonal per column) — would unblock large `nu_max` / coarse vertical grids under KPP / GLS without `Δt` shrinking. Listed in ROADMAP_COASTAL.md §7.
- **IMEX time stepping** — needed once implicit vertical diffusion lands.
- **Local time stepping** on heterogeneous unstructured meshes (sliver triangles currently dictate the global `Δt`).
- **Adaptive `Δt` based on tracer CFL** — currently set by gravity-wave + barotropic CFL only.

---

## API surface

### Shipped
- C-compatible shared library (`librakali.so`) via `iso_c_binding`, opaque-handle FFI.
- Python wrapper (`rakali`) — structured `Rakali` and unstructured `UnstructuredSimulation`.
- Mid-simulation setters: `h`, `hu`, `hv`, `S`, `T`, bathymetry, wind, heat flux, E−P.
- Zero-copy field access via `c_loc` → numpy views.
- Mesh helpers: Hilbert renumberer, coastal-mesh builder.
- YAML config loader (`rakali.yaml_loader`) for Python examples; Fortran-side namelist remains the canonical config format.

### Not yet shipped
- C API setters for EOS coefficients (`alpha_T`, `beta_S`, `rho_0`, refs)
- Schema-validating namelist reader (currently tolerant of unknown keys)
- Formal `SimulationProtocol` for ML surrogate swap-in

---

## Performance and scaling

### Shipped
- Single source `do concurrent` + OpenACC, runs on NVHPC GPU/multicore, gfortran, ifx (OpenMP-target build also exercised).
- Kernel fusion where profitable (`rki_kernel_fused_ml.F90`).
- Hilbert-curve cell renumbering for cache locality (mesh-build time).

### Headline numbers (single V100 unless noted)

| Configuration | Grid | Throughput |
|---|---|---|
| 2D barotropic | 1652×1396 | ~707 Mcells/s |
| 3D multilayer (5 sigma) | 1652×1396 | ~130 Mcells/s |
| 3D multilayer benchmark | 2000×2000 | ~340 Mcells/s |

### Not yet shipped
- Mixed-precision tracers (fp32/bf16)
- Multi-GPU / multi-node execution (MPI support was removed)
- GPU CI runner (GPU tests run locally today)

---

## Validation

Analytical tests shipped: Ritter dam break, Thacker bowl, Gaussian hump, lake-at-rest, solitary wave, thermal diffusion, lock exchange, river plume, multilayer dry-cell. 281+ unit tests under `testdrive`. CPU↔GPU bit-comparison via the gcc + nvhpc CI matrix.

Canonical benchmark configs live under `validation_examples/`:

- `validation_examples/coastal/` — `castro_dam_break/` (Castro et al.
  2011), `river_plume/`, `multilayer_benchmark/`, `native_bench/`. Each has a README with expected behaviour, analytical
  scales (when applicable), and diagnostic interpretation of failure modes.

Not shipped: GRS2018 storm-surge benchmarks, formalised Sydney Harbour observational comparison, golden-output regression harness, performance regression CI.

---

## Data assimilation

- **Shipped:** mid-simulation state setters (bathymetry, h, hu, hv, S, T, wind, heat flux, E−P) via Python — sufficient for nudging-by-hand and bench experiments.
- **Not shipped:** in-model nudging operator, ensemble support, observation operators beyond gauges, 4D-Var adjoint.

---

## Machine learning

Nothing built today. The user-loop API (`set_bc_*` → `advance(t)` → read state) is already the autoregressive-rollout pattern an ML surrogate would need. Scaffolding (`SimulationProtocol`, `mesh_features`, `make_ml_surrogate`, training-data loader) is sketched in ROADMAP_COASTAL.md §17 but deferred until a concrete model exists.

---

## Conservation contract (console Mass / Salt / Heat `Error`)

The periodic `[stats]` console block prints a `Mass : <total>  Error <residual>` line for each conserved quantity.  The `Error` meaning depends on the run configuration.

### What the `Error` measures

For a conserved quantity Q the residual is

```
Error = (Q_total - Q_ref) + out_Q - src_Q
```

where `Q_ref` is the value latched on the first status report, `out_Q` is the cumulative net outflux through open boundaries (positive = left the domain), and `src_Q` is the cumulative surface source (positive = added to the domain).  On a conservative, closed-domain, unforced run all three terms are zero and `Error` is the true numerical-leak residual **subject to the summation-order floor described below**.  With open boundaries or surface forcing the budget correction keeps `Error` a round-off residual rather than a spurious flag for the intended physics.

When a budget term is not instrumented (see limitations below) the console falls back to raw drift `(Q_total - Q_ref) / Q_ref`, byte-identical to the pre-budget behaviour.

### The summation-order floor

Every `Q_total` above is a plain floating-point `!$acc parallel loop reduction(+:acc)` device reduction.  It is not associative or order-deterministic: a changed GPU reduction-tree shape or thread count can change the last few bits of `Q_total`.  Two consequences:

- **`Total mass` / `Total salt` / `Total heat` are NOT bit-reproducible across thread counts or GPU reduction shapes** — the printed totals drift at round-off (~1e-14 to 1e-15 relative), even for bit-identical physics.
- **The `Error` floor is set by summation noise, not physics**, at that same ~1e-14 to 1e-15 relative band.

(An Extended-Fixed-Point reproducing-sum primitive, `src/framework/rki_efp.F90`, exists in the framework and is exercised by `test_efp`, but the console reductions in `rki_ml_common` / `compute_layer_sums` stay FP.)

### Instrumented paths

| Path | Mass | Salt | Heat |
|---|---|---|---|
| 2D barotropic (structured + unstructured) | closed residual | — | — |
| Multilayer `nz > 1` (structured + unstructured) | closed residual | raw drift | raw drift |
| Multilayer `nz = 1` fused fast path | raw drift | raw drift | raw drift |
| Semi-implicit (`time_scheme="semi_implicit"`, unstructured) | raw drift | raw drift | raw drift |

### Known limitations (fall back to raw drift)

- **Evaporation clamp edge case** (coastal multilayer, `evap_mass_sink`): in a thin-surface-layer + strong-evaporation + weak-mixing regime the non-conservative tracer clamp silently re-adds heat, causing a small non-zero residual even on an otherwise-closed budget — this is a known physical-model limitation, not a solver bug.
- **Coastal intertidal salt leak** (coastal multilayer, `ZSTAR_FULL` + intertidal): wet/dry cycling leaks 1–2 % salt per cycle; the Mass budget closes but salt/heat budgets are not instrumented on this path.

---

## What Rakali is not

- **Not a large-`Δt` semi-implicit model by default.** The DEFAULT path treats the free surface explicitly via HLL/KNP, so `Δt` is bounded by gravity-wave CFL `Δt < dx / sqrt(g·H)` — typically seconds at coastal scales. The opt-in `time_scheme="semi_implicit"` mode (structured 2D + sigma-multilayer 3D, plus the unstructured port) removes the gravity-wave CFL in the SCHISM/Casulli sense within its documented scope — see "Semi-implicit free surface" above; multi-rank remains explicit-only.
- **Not a wave model.** The non-hydrostatic solver is latent capability — a Stelling-Zijlema NH path that *could* propagate surface gravity waves, but it's not validated as a wave mode and there are no wave-maker / absorbing BCs tuned for surface gravity waves yet (ROADMAP_COASTAL.md §12A).
- **Not a BGC model.** No biology source terms, no FABM coupling.
- **Not a full morphodynamic model.** Single-class suspended sediment is shipped (structured), but bathymetry is fixed — no bedload / Exner / bed-level feedback yet.
- **No submesoscale / eddy closures.** Out of regime for a coastal solver.

---

## Compiler / platform support

- **NVHPC** (GPU + multicore stdpar) — primary target, exercised in CI.
- **gfortran** — exercised in CI for portability.
- **ifx** — buildable; some OpenMP-target map-clause edge cases observed in vendor tooling, not load-bearing for the canonical NVHPC build.
- **AMD / Intel GPUs** — unsupported today (no HIP / SYCL backend; would be a portability research project per ROADMAP_COASTAL.md §18).

---

## How to read this alongside the roadmap

`ROADMAP_COASTAL.md` carries the full work checklist with project-scope tags (`[E]/[M]/[H]/[V]`) and `Validate:` lines per planned item. This document summarises the *current* state; if there's a discrepancy, the codebase is the authority and the roadmap is the next-best source.
