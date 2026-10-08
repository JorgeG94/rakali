# Rakali — Physics & Numerics Reference

Canonical reference for the physics, numerical schemes, boundary
conditions, validation tests, and namelist configuration. **Codebase
structure, build system, GPU programming model, kernel API, and
Python/Julia interface live in [`CLAUDE.md`](../CLAUDE.md) +
[`docs/codebase/`](codebase/INDEX.md) — this file is the
science-and-config view.**

## 1. Overview

Portable GPU-native coastal solver written in Fortran. Two
grid backends share kernels:

- **Structured** Cartesian — 2D barotropic, 3D hydrostatic
  multilayer, non-hydrostatic. HLL / HLLC Riemann + hydrostatic
  reconstruction + minmod + SSP-RK2.
- **Unstructured** triangular — 2D barotropic + 3D hydrostatic
  multilayer. Kurganov-Noelle-Petrova central-upwind + Green-Gauss
  gradients + SSP-RK2. Cell-centric flux (each triangle accumulates
  into its own arrays only; interior edges computed twice but no
  scatter, no atomics).

## 2. Governing Equations

2D shallow water (conservative form):

```
∂h/∂t + ∂(hu)/∂x + ∂(hv)/∂y = 0

∂(hu)/∂t + ∂(hu² + ½gh²)/∂x + ∂(huv)/∂y =
    -gh B_x + τ_sx/ρ - c_f u √(u² + v²)

∂(hv)/∂t + ∂(huv)/∂x + ∂(hv² + ½gh²)/∂y =
    -gh B_y + τ_sy/ρ - c_f v √(u² + v²)
```

- `h` water depth
- `u, v` depth-averaged velocity components
- `g` gravitational acceleration
- `B` bathymetry (bottom elevation, positive up)
- `τ_s` surface wind stress
- `c_f` bottom friction coefficient
- `ρ` water density

3D hydrostatic multilayer extends this with per-layer momentum + a
prognostic salinity/temperature tracer registry; see §5.

## 3. Numerical Schemes

### Structured

- **Spatial**: finite volume on Cartesian grid
- **Riemann**: HLL approximate solver with hydrostatic reconstruction
  (or HLLC if `use_hllc=.true.` in `&physics_nml`)
- **Reconstruction**: 2nd-order minmod slope limiter
- **Time**: selectable explicit SSP Runge-Kutta — `&time_nml time_scheme`
  = `euler` / `rk2` (default) / `rk3` / `rk43` (see "Time integrator
  menu" below); non-hydrostatic keeps SSP-RK2.  The implicit-gravity-wave
  `semi_implicit` mode exists on the unstructured backend only (see
  "Semi-implicit free surface" below).
- **Wetting/drying**: hydrostatic reconstruction with Kurganov & Petrova
  (2007) desingularized thin-cell velocity (`DESING_EPS`) and a
  high-Froude blowup backstop (`FROUDE_BACKSTOP`)

### Unstructured

- **Spatial**: finite volume on triangular mesh
- **Flux**: Kurganov-Noelle-Petrova central-upwind (no Riemann solver
  needed)
- **Reconstruction**: 2nd-order via per-cell Green-Gauss gradient
- **Time**: selectable explicit SSP Runge-Kutta — same `&time_nml
  time_scheme` menu as the structured path (`euler`/`rk2`/`rk3`/`rk43`,
  default `rk2`; see "Time integrator menu" below). 2D barotropic path
  only; multilayer keeps SSP-RK2.
- **Wetting/drying**: negative-depth protection kernel
- **Mesh**: UGRID-convention NetCDF (node coords, cell-to-node,
  edge-to-cell, boundary edge flags)
- **GPU pattern**: cell-centric double-visit. Each triangle visits
  its three edges and writes only its own arrays. Interior edges
  computed twice. No scatter, no atomics

### Design rationale

Explicit time-stepping maps to GPU parallelism without global solves
or implicit systems. CFL constraint is acceptable for coastal
applications where short gravity waves require small timesteps anyway.
`do concurrent` is portable across GPU vendors with a CPU fallback
without separate kernel implementations.

### Time integrator menu (2D barotropic, structured + unstructured)

`&time_nml time_scheme` selects the explicit integrator on both the
structured (HLL/HLLC) and unstructured (KNP) 2D barotropic paths
(plus the `semi_implicit` mode — next section; unstructured only, 2D
barotropic + 3D sigma multilayer). All
four explicit schemes are strong-stability-preserving (SSP), so each stage is a convex
blend of forward-Euler updates and the wetting/drying positivity of
forward Euler carries through. The scheme is parametrized by a Shu-Osher
coefficient table on a concrete `timestepper_t` (data, not type
hierarchy) and every scheme reuses the same forward-Euler stage
kernel — adding a scheme is a table entry, not a code path. The same
`timestepper_t` drives both backends; each contributes its own
fused stage primitive + convex-combine kernel.

| `time_scheme` | stages | order | SSP coeff `C` | cost / sim-time |
|---|---|---|---|---|
| `euler` | 1 | 1 | 1 | 1× |
| `rk2` *(default)* | 2 | 2 | 1 | 2× |
| `rk3` | 3 | 3 | 1 | 3× |
| `rk43` | 4 | 3 | **2** | **2×** (4 stages ÷ 2× stable dt) |

`C` is the factor by which the stable timestep exceeds the
forward-Euler CFL limit: `rk43` is stable at 2× the dt of the others,
so its CFL-derived `dt` is scaled by 2 and four stages cover the same
simulated time per step as two `rk2` stages. The temporal order is
capped at 3 deliberately — the barotropic flux is 2nd-order
(HLL/HLLC + minmod), so RK4 buys nothing until a higher-order
barotropic reconstruction lands; `rk3` already shrinks the
temporal-error constant. `rk2` is bit-identical to the historical
hard-wired SSP-RK2 (it dispatches to the same fused fast path). The
menu is wired into the 2D barotropic path only; multilayer and
non-hydrostatic runs keep SSP-RK2 (fail-loud if a non-`rk2` scheme is
requested with `use_multilayer`/`use_nonhydrostatic`). Verify the
convergence orders + work-precision with
`benchmark_timestepping_single_gpu`.

### Semi-implicit free surface (`time_scheme="semi_implicit"`)

Unstructured (triangular) grids only: `time_scheme="semi_implicit"` with
`grid_type="structured"` is rejected at configure time.

Casulli-style semi-implicit gravity wave: the free-surface PGF is
treated θ-implicitly (θ=0.6 fixed) and `η^{n+1}` solved by a matrix-free
CG on the SPD Helmholtz system `[1 − gθ²Δt²∇·(Hf∇)]η = RHS`
(warm-started, relative tol `&si_nml rtol`, default 1e-10) — no
gravity-wave CFL, so `dt` is minutes instead of seconds (requires
`dt_fixed > 0`; barotropic Courant 10–20+ validated). Internally the
face transports live on edges and face-flux continuity gives exact mass
conservation to round-off; the cell momenta the rest of the code sees
are reconstructed from the edge transports at step end.

**Geometry**: TPFA edge fluxes (Casulli & Walters 2000) with η at the
triangle **circumcentre** (NOT the KNP centroid: a Delaunay edge is
orthogonal to the circumcentre line by construction, which is what
makes the two-point flux consistent — prototype-measured O(h²)
convergence vs a plateau on centroids), per-cell double-visit edge
transports `(3, ncell)`, and an **area-weighted** CG (the cell-centred
FV operator is self-adjoint only in `⟨x,y⟩ = Σ area_c x_c y_c`; plain
dots stall CG on a variable-area mesh). A **fail-loud mesh-quality
gate** at enable-time rejects meshes with max circumcentre-line
deviation > 5° or non-positive circumcentre distance (obtuse/degenerate
triangles — a clean Delaunay mesh measures < 1e-2°; an axis-aligned
right-triangle split fails BY DESIGN, its circumcentres sit on the
shared hypotenuse).

**Preconditioning (`&si_nml precond`)**: `"jacobi"` does real work even
at uniform depth (variable cell areas put O(100-1000)× spread on the
diagonal); `"mg"` is AGGREGATION (agglomeration) multigrid
(`rki_si_unstr_mg`: deterministic greedy pairwise
strongest-transmissibility aggregation, Galerkin coarse operators,
V(2,2) Jacobi ping-pong — the plateau-breaker; measured 1526 -> 16 CG
iterations on the 480-triangle funnel). 2D barotropic only in v1
(mg + wetdry / multilayer aborts at enable).

**2D barotropic scope (P3.1)**: wall + tidal/clamped Dirichlet-η
boundary edges via the per-group `bc_type`/`bc_eta_ext` tables, linear
predictor (no momentum self-advection). Kernels:
`src/core/coastal/kernels/unstructured/barotropic/rki_si_unstr_*.F90`,
state slot `unstr_state_t%si` (`rki_si_unstr_state`), mesh additions
(circumcentres, `edge_delta`, deviation angles) on `mesh_unstr_t`.

**3D multilayer scope (P3.2)**: `use_multilayer` + `semi_implicit` runs
the column-coupled step — per-cell column solve in `(nz, ncells)` order
(`rki_si_unstr_column`): per-column backward-Euler implicit vertical
viscosity (closure ν — PP81 by default, k-ε / KPP opt-in, see the
closure paragraph below) via Thomas factorisation and three back-solves
(the G/R linear response whose `R_k ≡ −1` row-sum identity keeps the
Helmholtz side ν-independent), implicit **Manning-only** bottom drag on
the bed row (the unstructured barotropic state carries no linear-drag
knob), the effective face depth `Ĥ = −Σ h_k R_k` + implicitly-mixed
transport `Ğ = Σ h_k G_k` feeding the UNCHANGED TPFA eta-solve, and an
**anchor + increment** per-layer edge recovery: the persisted per-layer
edge transports are the anchor and only the implicit-mixing +
explicit-tendency increment `h_k·G_k − hu⁰_k` routes through the cell
average — so the uniform-density inviscid limit reduces EXACTLY to the
edge-persistent 2D step (measured 1.4e-10 m pointwise over 5 seiche
periods). Tracers: registry advection with dsig-proportional first-order
upwind, driven by the EXACT blended edge mass flux the eta update used
(uniform tracers invariant to CG tolerance — verified linear in `rtol`;
monotone only while the FLOW Courant `|u|Δt/Δx ≤ 1`; the gravity-wave
Courant may be 20 but fast advective flow at large dt is outside the
envelope). **Surface wind stress** applies on the top layer (uniform
scalar τ). EOS runs BEFORE the BPG each step (the unstructured house
order — `rho_layer` is zero-initialised, so an EOS-lagged first-step BPG
would see a −ρ₀ anomaly kick). **Tidal/clamped open boundaries** work on
the 3D path (the same Dirichlet-η seed + open-face `Ĥ`; per-layer
zero-gradient baroclinic velocity; tracers zero-gradient on outflow /
clamped inflow); the OBC seed is STATELESS (ghost η re-derived from
`bc%time` + tidal parameters every step). Fail-loud limits: `VCOORD_SIGMA` only,
boundaries wall/tidal/clamped only, `nz ≥ 2`, PP81/k-ε/KPP closures, no
sediment, baroclinic PGF explicit, no predictor momentum-advection.

**Wetting/drying (`&si_nml wetdry`, opt-in)**: Casulli (2009) subgrid
wet/dry replaces the Helmholtz constant `1.0` cell-mass term with the
convex wetted-volume `V(η)` (derivative = wet fraction `a(η) ∈ [0,1]`,
the operator's variable diagonal), solved by a monotone plain-Newton
outer loop (~1.2–3 iterates/step, cap `newton_max`=10) wrapping the
unchanged CG/Jacobi inner solve, with GATED-MEAN frozen edge conveyance
(all-wet ≡ the linear edge depth) and structural dead-edge masking
(per-edge transports zeroed on Hf=0 — verified by the tidal-flat
water-budget closure, 9e-13 rel); `rki_si_unstr_wetdry`, enabled via
`state_unstr_enable_si(..., wetdry=.true., n_sub=...)`. Works on the 2D
barotropic and 3D sigma-multilayer paths, wall + tidal/clamped edges. The
quality gate is hardened for this phase: `min_edge_delta` tracks the
SIGNED circumcentre-line projection on the outward normal (a locally
non-Delaunay flipped pair has a negative projection; the deviation angle
alone cannot detect it). The wet/dry state is **stateless across
restart** (`b_sub`/`b_min`/`Hf`/wet-mask rebuilt each step from
`h + b`). Under wet/dry `si%cg_iters` is the SUM of inner-CG iterations
over the Newton iterates, and `si%mass_res` is the nonlinear residual
sum (true mass defect), not a single CG residual. **v1 limits**: the
subgrid sub-bottom table is live via `&si_nml wetdry_bsub_fill` (`"flat"`
default = bit-identical to `nsub=1`; `"linear"` = the analytic Casulli
2009 per-cell relief ramp, mean sub-bottom exactly `b` so still-water
volume is unchanged) — but `"linear"` is 2D-barotropic only: with
multilayer it is fenced fail-loud at configure AND state-init (the
η-solve conserves the subgrid `v_cell(η)` while the ML thickness update
consumes the plain bed — mass would leak at the front); the
thin re-wetting front carries a documented ~2.3 PSU first-order-upwind
tracer overshoot in v1 (well-wet interior stays uniform to ~0.05 PSU).
Default off leaves the constant mass term in place, byte-identical.

**SI vertical-mixing closures (Q1)**: `vmix_use_keps` and
`vmix_use_kpp` work under `semi_implicit`. The column solve consumes any
interface ν through `nu_iface` (`R = solve(−1)` is ν-independent), so a
closure is a pure fill: **k-ε** runs the unchanged `ml_keps_update_unstr`
once per SI step (it is already backward-Euler/Patankar — implicit
Thomas diffusion, implicit `ε/k`-linearized sinks, positivity-clamped
sources; dt = 120–480 s certified by the prototype in
`local_archive/prototypes/si_closures/`), then `keps_nut → nu_iface`
(momentum mixing happens IN the column tridiagonal;
`ml_keps_apply_momentum` is not called) and `keps_kt` drives the
implicit per-tracer diffusion. `&nonhydrostatic_nml keps_si_picard`
(default 1) optionally re-runs the k/ε solve with the refreshed ν_t
production coefficient (pass 1 is steady-state-exact; pass 2 for
transient-fidelity studies). **KPP** reuses `ml_kpp_boundary_layer_unstr`
verbatim, MAXes the LMD94 overlay `h_b·w_s·G(σ)` onto PP81 into
`nu_iface` with no explicit-flux CFL clamp, and solves S/T diffusion
implicitly with the overlay diffusivity + the conservative non-local
counter-gradient divergence (`rki_si_unstr_closures`); extra passive
tracers fall back to PP81. Closures off ⇒ the untouched PP81 branch
(byte-identity by construction).  Tests: `test_si_closures_unstr`.

Tests: `test_si_unstr_seiche` (triangulated seiche: period + mass +
Jacobi stability + quality-gate accept/reject), `test_si_unstr_channel`
(co-oscillating tidal channel, clamped Dirichlet edge),
`test_si_unstr_mg` (aggregation-MG vs plain CG), `test_si_unstr_ml_baro_limit`
(nz=4 ≡ 2D to 1.4e-10 m), `test_si_unstr_ml_tidal_tracer` (uniform-S
pumping invariance, rtol-scaling verified), `test_si_unstr_wetdry`
(Thacker on triangles + tidal-flat budget + all-wet agreement),
`test_si_unstr_ml_wetdry`, `test_si_wetdry_bsub` (sub-bottom fill),
`test_si_unstr_estuary_bench` (funnel-estuary SI-vs-explicit
head-to-head: amp/phase gates + logged wall-clock).

## 4. Boundary Conditions

Integer codes match across backends so the same value means the same
BC type on either side.

### Structured (9 types)

`rki_boundary` (kernel) + `rki_boundary_types` (codes):

| Type | Code | Description |
|------|------|-------------|
| `BC_WALL` | 1 | Reflective wall with eta mirroring |
| `BC_OPEN` | 2 | Flather radiation |
| `BC_TIDAL` | 3 | Prescribed elevation with Flather momentum |
| `BC_NESTED` | 4 | One-way nesting from coarser parent |
| `BC_INFLOW` | 5 | Prescribed normal velocity, extrapolated depth |
| `BC_DISCHARGE` | 6 | Prescribed volume flux |
| `BC_CLAMPED` | 7 | Dirichlet elevation, zero-gradient velocity |
| `BC_SPONGE` | 8 | Quadratic damping toward reference (separate kernel) |
| `BC_CHAPMAN` | 9 | Orlanski radiation with Flather momentum |

`BC_SPONGE` ghost cells fall through to wall; sponge relaxation is a
separate kernel called after friction. `BC_CHAPMAN` carries persistent
`eta_old_*` arrays in `bc_params_t`, allocated unconditionally by
`bc_init_chapman` and GPU-mapped in `solver_enter_data`.

### Unstructured (7 types)

Dispatched inside `compute_boundary_state` / `rki_kernel_flux_unstr`.
Boundary edges (where `mesh%neighbours(i, k) == -1`) carry a per-edge
group id that maps to a per-group BC type.

| Type | Code | Description | Per-group parameter |
|------|------|-------------|---------------------|
| `BC_WALL` | 1 | Reflective mirror | – |
| `BC_OPEN` | 2 | Transmissive zero-gradient | – |
| `BC_TIDAL` | 3 | Prescribed eta + Flather | `bc_eta_ext(g)` |
| `BC_INFLOW` | 5 | Prescribed outward-normal velocity | `bc_un_ext(g)` |
| `BC_DISCHARGE` | 6 | Prescribed specific discharge | `bc_Qn_ext(g)` |
| `BC_CLAMPED` | 7 | Dirichlet eta, transmissive momentum | `bc_eta_ext(g)` |
| `BC_CHAPMAN` | 9 | Orlanski radiation | per-edge `chapman_eta_old(i, k)` |
| `BC_OUTFLOW` | 10 | Prescribed depth, transmissive momentum | `bc_h_ext(g)` |

`BC_CHAPMAN` carries per-cell-edge state (not per-group) — each
radiating edge needs its own previous-step eta. Per-edge writes are
race-free under cell-centric flux pattern.

`BC_OUTFLOW = 10` is unstructured-only (the code is past the structured
range to keep codes collision-free if ever ported).

### Multilayer

The unstructured tracer advection kernel honours per-group
`bc_S_ext(g)` / `bc_T_ext(g)` at inflow edges so freshwater /
clamped-tracer river inflows work end-to-end.

### Coverage gap

| Type | Code | Notes for porting |
|------|------|-------------------|
| `BC_NESTED` | 4 | Needs unstructured analogue of `rki_nesting` |
| `BC_SPONGE` | 8 | Cell-centric implementation should be straightforward |

## 5. Hydrostatic Multilayer Solver

Enabled via `use_multilayer = .true.` in `&nonhydrostatic_nml`. Sigma
or hybrid (z*-lite, z*-full, z-sigma) coordinates via `vcoord_type`.

### Architecture: mode-split barotropic-then-distribute

Single barotropic HLL step (exactly well-balanced, 2D), then
distribute acceleration uniformly to all layers preserving vertical
shear. Avoids per-layer HLL and its associated sigma-coordinate
pressure-gradient errors.

Each SSP-RK2 stage:

1. Barotropic flux + update (2D, well-balanced)
2. Distribute barotropic acceleration to layers
3. Tracer advection (horizontal) + horizontal diffusion + vertical exchange
4. Equation of state → density
5. Vertical viscosity / diffusion (momentum + tracers)
6. Baroclinic pressure gradient (reduced-gravity cap)
7. Per-layer bottom friction (Manning on bottom layer)
8. Smagorinsky horizontal viscosity
9. Truncation safety net
10. Recompute barotropic momentum from layer sum
11. Per-layer Coriolis rotation

Stage 2 also RK2-averages and clamps tracers.

### Diagnostic vertical velocity `w_ml` (diasurface, structured)

`state%multilayer%w_ml(:,:,k)` is the **diasurface** (sigma-relative)
volume flux per unit area across the moving sigma interface `k` — the
quantity the vertical-exchange kernels advect momentum and tracers
with, **not** the physical Cartesian vertical velocity `w_phys`.

In sigma coordinates the interfaces move with the free surface, so
`∂h_k/∂t = dsig_k·∂H/∂t`.  From per-layer continuity, integrated
bottom-up with `w_ml(:,:,1) = 0` (no bed penetration):

```
w_ml(k+1) = w_ml(k) − div(hu_k) + dsig_k·div_col
div_col   = Σ_m div(hu_m) = −∂H/∂t
```

The `dsig_k·div_col` term is the moving-coordinate correction.  Two
consequences (both verified by `test_w_diasurface_{barotropic,baroclinic}`):

- Pure **barotropic** sigma flow (`hu_k = dsig_k·hu_bt`) ⇒ every
  interface flux is 0 to round-off — no spurious diapycnal transport.
  The former `w_ml(k+1)=w_ml(k)−div(hu_k)` form omitted the correction,
  reported `w_phys`, and produced O(1e-3) spurious interior `w_ml` that
  drove a slow salt leak through the vertical exchange.
- The free surface `w_ml(:,:,nz+1) ≡ 0` **by construction** (the
  moving-top kinematic BC).

**Baroclinic** (sheared) flow still yields sensible nonzero interior
`w_ml`.  The `dsig_k·∂H/∂t` interface-motion model is exact for pure
sigma layering (`VCOORD_SIGMA`, the coastal default); for hybrid z/z*
coordinates it is an approximation.  The unstructured twin
(`ml_compute_vertical_velocity_unstr`) tracks only the barotropic mass
flux and so is a distinct (barotropic-only) formulation — see its
docstring.

### Vertical layer convention (load-bearing)

Both backends are **bottom-up (ROMS-style)**: `k=1` is the bed, `k=nz`
is the surface. Surface forcings land at `(:, :, nz)` / `(nz, cell)`;
bed forcings at `(:, :, 1)` / `(1, cell)`. Three regression gates:
`test_lock_exchange` (dense current at `k=1`), `test_river_plume`
(buoyant plume at `k=nz`), `test_surface_flux_strat` (positive Q warms
`k=nz` and preserves stratification).

### Vertical coordinates

Dispatched via `vcoord_type` (namelist) → `parse_vcoord_type` →
`VCOORD_*` enum:

| Tag | Coord | Status | Notes |
|---|---|---|---|
| `sigma` | Pure sigma | Production | No remap; layers rescaled by `dsig(k) * H`. Default for coastal. |
| `zsigma` | Smoothstep blend sigma→z | Production | Conservative remap. |
| `zstar` | z*-lite | Production | Single global `z_ref` stretched per column. SSH-tracking. |
| `zstar_sigma` | Sigma + z*-lite hybrid | Production | Sigma shallow, z*-lite deep. Conserves by construction. |
| `zstar_full` | Full per-column z* | Production w/ caveat | Per-column `z_ref(0:nz)` from local bathymetry; surface anchored at `zstar_h_surf_target` regardless of H. Bed-side layers can vanish (`zstar_h_min`). Operators that divide by `h_layer` gate on `VANISHING_LAYER_TOL = 1.5e-4 m`. Pairs with gated CWC + conservative clamp. **Caveat**: intertidal estuaries still leak 1-2% salt/cycle from wet/dry destruction. |

### Tracer transport (CWC)

Both backends use a **vanishing-layer-gated centred CWC** per-face
mass-flux split with per-face symmetric redistribution. Guarantees
per-layer conservation under deep-shallow transitions; reduces to
centred CWC for benign bathymetry.

Two prognostic tracers are registered by default:

- **Salinity** `S` (PSU), conserved as `hS = h_layer × S`
- **Temperature** `T` (°C), conserved as `hT = h_layer × T`

Both use the same generic tracer kernels — each is called twice per
solver step, once for S and once for T, with tracer-specific scalars
(inflow value, clamp bounds, background diffusivity) passed as
arguments. Adding a sediment / age / BGC tracer is a registry append;
no new kernel code.

### Equation of state

Linear in both tracers:

```
ρ_k = ρ_0 + β_S × (S_k - S_ref) - α_T × (T_k - T_ref)
```

`β_S` (default 0.78 kg/m³/PSU) is the haline contraction coefficient;
`α_T` (default 0.17 kg/m³/°C) is the thermal expansion coefficient.
Setting `α_T = 0` reduces temperature to a passive tracer with no
dynamical feedback.

### Baroclinic pressure gradient

Two methods via `bpg_method` in `&nonhydrostatic_nml`:

- **Sigma BPG** (`bpg_method = "sigma"`): top-down pressure
  accumulation along sigma surfaces with central differences. Fast
  but produces pressure-gradient error proportional to
  `β_S × S × dH/dx × nz` at steep bathymetry — fundamental
  sigma-coordinate limitation. Unstable with 5 layers on real
  bathymetry after ~20 minutes.
- **Jacobian BPG** (`bpg_method = "jacobian"`): Shchepetkin &
  McWilliams (2003) pressure Jacobian. Interpolates neighbour-column
  pressure to the centre column's physical depth, computing gradients
  along constant-z surfaces. With uniform density, BPG is exactly zero
  regardless of bathymetry. **Recommended for real-world applications.**

Both include a per-step acceleration cap using reduced gravity:
`sqrt(g' × H) / dt` where `g' = g × |ρ_top − ρ_bottom| / ρ_0`. Scales
the cap with actual baroclinic forcing rather than the barotropic
wave speed.

### Per-layer bottom friction

Manning applied to bottom layer (k=1) using the bottom layer velocity
rather than depth-averaged. Semi-implicit:

```
hu_new = hu_old / (1 + dt × C_d × |u_b| / h_k)
```

with `C_d = g × n² / H^{1/3}`. Damps baroclinic shear flows that
depth-averaged friction cannot see.

### Vertical mixing (Pacanowski-Philander 1981)

Richardson-number-dependent:

- Momentum: `ν_t = ν_max / (1 + α_Ri × Ri)² + ν_bg`
- Salinity: `κ_t = ν_max / (1 + α_Ri × Ri)³ + κ_S_bg`
- Temperature: `κ_t = ν_max / (1 + α_Ri × Ri)³ + κ_T_bg`

`α_Ri` shared between momentum and tracers (standard PP81).
Background diffusivities independent. `Ri = g × dρ / (ρ_0 × (du² + dv²) / dz²)`.

### Smagorinsky horizontal viscosity

Grid-aware eddy viscosity for layer momentum:

```
ν_h = (C_s × Δ)² × |S|
|S|  = √( (∂u/∂x − ∂v/∂y)² + (∂u/∂y + ∂v/∂x)² )
Δ    = √(dx × dy)
```

Stability-limited to `0.5 / (dt × (1/dx² + 1/dy²))`. `smagorinsky_cs`
default 0; typical 0.1–0.2.

### Surface tracer fluxes

Prescribed surface forcings applied to top layer (k=nz) only:

- **Heat flux** `heat_flux(x, y)` [W/m²], positive into the water.
  `hT(nz) += dt × Q_net / (ρ_0 × c_p)` with `c_p = 4000 J/(kg·K)`.
- **Virtual freshwater flux** `evap_minus_precip(x, y)` [m/s],
  positive for evaporation. `hS(nz) += dt × (E − P) × S_top` —
  ROMS/MOM virtual-salt-flux convention; concentrates salt without
  removing volume.

Neither forcing changes `h_layer`. Real volume forcing is handled by
`precipitation` (in `&physics_nml`) through `apply_friction` on the
barotropic `h`.

**Evaporation as a real mass sink** (`&nonhydrostatic_nml
evap_mass_sink`, default `.false.` = bit-identical legacy behaviour;
structured multilayer only): when on, net evaporation
(`evap_minus_precip > 0`) is removed as REAL water instead of the
virtual salt flux.  Per step, `dh = E × dt` is taken
column-proportionally — every layer's `h`, `hu`, `hv`, `hT` scales by
`1 − dh/h_col` while `hS` is held.  Consequences (exact for any
profile): salt mass is conserved and the concentration `S = hS/h`
rises (evaporation removes water, not salt); the evaporated water
carries its own heat and momentum away, so `T`, `u`, `v` are
invariant — no spurious warming or acceleration.  The proportional
form preserves the current layer layout (sigma fractions untouched)
so the vcoord distribute/remap machinery sees a consistent column.
Wet/dry: the sink clamps to the water available above `ml_min_depth`;
the unsatisfied demand is dropped honestly (no compensating source).
Precipitation stays on the barotropic `precipitation` path in this
mode.  Ported to **both** backends (structured + unstructured): the
unstructured kernel removes water column-proportionally per cell and
the following `ml_compute_barotropic_momentum_unstr` propagates the
drawdown into `height_cv`/`stage_cv`.

**Latent-heat cooling** (`evap_latent_cooling`, default `.false.` =
bit-identical; requires `evap_mass_sink=.true.`, fail-loud at configure
otherwise): the evaporated mass additionally withdraws its latent heat
of vaporization from the surface layer (`k=nz`) — an ADDITIVE cooling
`d(hT_nz) = −(ρ_water·L_v / (ρ_0·c_p))·dh_eff` [°C·m], with
`L_v = 2.5×10⁶ J/kg` (`LATENT_HEAT_VAPORIZATION`) and `dh_eff` the
effective (clamped) evaporated depth.  This is DISTINCT from — and does
NOT double-count — the sensible heat that leaves with the removed water
mass (that is the proportional `hT` scaling above: phase-change energy
vs. enthalpy of the departing liquid).  A dried-out column withdraws no
latent heat.  The column-mean temperature then follows the continuum
drawdown-cooling law `T_mean(t) = T_INIT − L·ln(H0/h(t))`.

Density-normalization convention: the sensible and latent channels are
normalized to slightly different departed masses by design.  The
sensible term removes `ρ_0·dh_eff` (the model's reference-density water
leaving with its own heat), while the latent term is normalized to
`ρ_water·dh_eff` so it matches the standard bulk latent-heat flux
`Q_lat = ρ_water·E·L_v`.  The two densities differ by
`ρ_water/ρ_0 ≈ 3.5%`; each channel uses the density natural to its own
formulation — this is the intended convention, not an inconsistency.

Energy-closure caveat: the full latent withdrawal lands in the surface
layer alone, so it cools `k=nz` by `latent_gain·dh_eff/h_top`.  For a
THIN surface layer (small z* surface target or strong stratification)
under strong evaporation with weak or absent vertical mixing, `T(nz)`
can be driven below `T_min`, where the NON-CONSERVATIVE tracer clamp
(`rki_solver`) fires and silently re-adds heat.  The exact latent energy
budget therefore holds only while `T(nz)` stays above `T_min`; in that
thin-surface-layer + strong-evaporation + weak-mixing regime the closure
is no longer exact.  Spreading the latent sink over a surface band
(rather than the top layer alone) is the deferred remedy.

Remaining v1 limitations (documented follow-ups): the evaporative
salinity rise is column-uniform, not surface-intensified — a
surface-localised sink needs the corrected physical `w` for w-coupled
tracer transport (branch `fix/ml-vertical-velocity-w`), deferred until
that lands; KPP's boundary-layer detector and non-local transport still
estimate the surface salt buoyancy flux via the virtual-salt proxy
(`evap_minus_precip × S_top`) regardless of the knob — a physically
reasonable proxy for the concentration-rise buoyancy loss, not a salt
double-count, but the two paths are not identical.
Tests: `test_coastal_evap` (structured: analytic drawdown
`h(t) = h0 − E·t`, exact salt-mass constancy with the analytic
`S(t) = h0·S0/h(t)` concentration rise, temperature invariance, a
thin-flat graceful-drying case, exact latent-sink bookkeeping +
energy-closure, and the continuum cooling law through the full solver)
and `test_multilayer_unstr` (unstructured mass-sink drawdown + latent
bookkeeping twins).

### Coriolis

Exact rotation matrix formulation (energy-conserving) applied
per-layer. Latitude → `f` set in `&physics_nml`.

### Wetting/drying stability (multilayer)

Three mechanisms:

1. **`ML_MIN_DEPTH` guard**: columns shallower than `ml_min_depth`
   (default 0.5 m) skip vertical exchange, mixing, tracer advection,
   and BPG.
2. **Thin-layer zeroing**: layer momentum zeroed when
   `h_layer < THIN_LAYER`.
3. **Vertical exchange limiter**: dry-neighbour guard on vertical
   velocity + momentum flux limited to 50% of donor layer.

### CFL-based velocity truncation

Any layer velocity producing CFL > 0.8 is clipped to `0.9 × 0.8 × dx/dt`.
Runs after all multilayer physics, before recomputing barotropic momentum.

### 2D depth-averaged tracer mode

When `nz_layers = 1`, the multilayer pipeline degenerates: vertical
mixing, vertical exchange, and BPG kernels early-return on `nz < 2`,
leaving only barotropic update plus horizontal tracer advection.
Recipe for passive tracer transport (sediment, dye, salt) on top of
pure 2D shallow water without paying for 3D physics.

### Tracer face-reconstruction ladder (structured coastal)

`&tracer_nml tracer_recon` selects the FACE-VALUE reconstruction inside
the existing flux-form horizontal tracer advection (the CWC mass-flux
split, per-layer conservation shell, and downstream clamp are unchanged
— only the upwind face value swaps):

| `tracer_recon` | scheme | nghost | full-period L1 order | notes |
|---|---|---|---|---|
| `ppm` (default) | Colella-Woodward 1984 parabola | 2 | 2 | bit-identical to the pre-existing path |
| `weno5` | Jiang-Shu 1996 + Borges 2008 Z-weights, 3 quadratics | 3 | 3 | the recommended rung |
| `weno7` | Balsara-Shu 2000, 4 cubics, Z-weights | 4 | 4 | best front L1 |
| `weno9` | Balsara-Shu 2000, 5 quartics (simplified 3-term β), Z-weights | 5 | ~6 | documented-expensive |
| `pqm` | — | — | — | **deferred, fail-loud** (see below) |

All rungs use the **CFL-swept-average** form (the face value integrates
the upwind reconstruction over the swept region, Godunov-style) — this
is single-stage stable, unlike the classical point-value WENO which is
weakly unstable in this integrator and needs a multi-stage RK. Z-weights
(Borges et al. 2008) are used throughout: better than Jiang-Shu weights
at coarse resolution, never worse. Full swept-average derivation +
coefficient tables: `local_archive/prototypes/tracer_weno/coefficients.md`.

**Boundary order-reduction ladder.** As a face approaches a wall / open
boundary and the wide stencil no longer fits, the rung degrades per face
`weno9 → weno7 → weno5 → plm → donor` (a stencil may reach at most one
ghost row, since deeper ghosts are zero-gradient copies that carry no
information). Guaranteed in-bounds and conservative; tested in
`test_tracer_recon` (`rung_ladder_table`, `wall_adjacent_weno9_smoke`).

**Per-scheme `nghost` is fail-loud at configure** (weno5→3, weno7→4,
weno9→5): if `&grid_nml nghost` is smaller, `validate_config` aborts
naming both numbers. A weno* `tracer_recon` is also rejected fail-loud on
any config that will not honour it — a run without
multilayer tracers (2D barotropic / NH), with `use_ppm_tracer=.false.`
(contradiction). These
support checks run for every config (not gated on multilayer), so a weno*
rung can never be silently parsed-and-dropped. `pqm`
parses but aborts with a defer message on all configs: it was measured as dominated by
`weno5` under any practical monotonicity limiter, and the faithful
White-Adcroft 2008 inflection-point limiter is a standalone follow-up
(evidence in `local_archive/prototypes/tracer_weno/PROTOTYPE_REPORT.md`).

**GPU cost (V100, sm_70, the single WENO9-inclusive kernel).** 180
registers, **no spills** (16 KB stack frame, 0 spill stores/loads) ⇒
~17% theoretical occupancy — register-bound but not spilling (contrast
the 255-register kappa-shear kernel). Cost rises with the rung's stencil
width; `weno5` is the accuracy-per-cost sweet spot. Default `ppm`
compiles the same kernel but never dispatches into it (bit-identical).

**Known-behaviour note — PPM overshoot.** The default PPM path is NOT
strictly monotone at sharp fronts: at a square-wave corner it overshoots
by ~1.8e-2 (≈3% of the jump), because the Colella-Woodward eq. 1.10
limiter without van-Leer-limited edge interpolation does not bound face
values by neighbour means. This is pre-existing, intended behaviour —
the downstream conservative tracer clamp is the safety net (see the
`test_2d_tracer` comments). The WENO rungs overshoot two orders LESS
(weno5 ≤1e-3 in the same test) and rely on the same clamp; no new
limiter was added.

## 6. Non-Hydrostatic Extension

Multi-layer sigma + conjugate gradient (CG) Poisson solver. Structured
only.

- **Pressure solver**: CG with warm start, configurable tolerance and
  max iterations.
- **Depth guard**: `NH_MIN_DEPTH` (separate from `DRY_TOLERANCE`) to
  avoid 1/h singularities in the Poisson operator.
- **Mode splitting**: barotropic-baroclinic split via `split_ratio`
  (1=unsplit, 0=auto, N=fixed). N fast barotropic substeps (SSP-RK2
  each) per slow baroclinic step. NH pressure gradient injected as
  constant forcing during barotropic subcycling. ROMS-style time
  averaging.
- **4D output**: layer fields (`u_layer`, `v_layer`, `w`, `q_nh`) and
  vorticity diagnostics.

## 7. Validation Test Suite

Tests use the `test-drive` framework. PASS/FAIL against analytical
solutions or regression baselines.

| Test | What it tests | Notes |
|---|---|---|
| `test_dambreak` | Ritter shock + rarefaction, basic conservation | 1D on 2D grid. L2 ~O(dx) near shock, O(dx²) smooth. |
| `test_thacker` | Wetting/drying, moving shoreline, mass conservation | Parabolic bowl, periodic sloshing. Shoreline error < 2 cells. |
| `test_lake_at_rest` | Well-balanced property | Variable bathy, flat surface. Velocities stay at machine zero. |
| `test_convergence` | Spatial convergence rate | 2nd-order in smooth regions. |
| `test_hump_radiation` | Wave dispersion, symmetry, open BC | Circular symmetry preserved at open boundary. |
| `test_solitary_nh` | NH dispersion | NH wave propagates faster than hydrostatic. |
| `test_lock_exchange` | Baroclinic dynamics, EOS, gravity current | Dense current at `k=1`, fresh at `k=nz`. |
| `test_river_plume` | Buoyancy-driven plume + vertical mixing | Surface-trapped buoyant plume at `k=nz`. |
| `test_surface_flux_strat` | Surface heat flux warms top, preserves stratification | Positive Q warms `k=nz`. |
| `test_dry_cell_multilayer` | W/D stability with sigma layers on complex bathy | No blowup at fronts, mass conserved. |
| `test_zstar_full_*` (suite) | z*-full intertidal, conservation, layer protection | See `tests/test_zstar_full_*.F90`. |

Plus ~80 unit tests covering state, config, BCs, flux,
timestep, update, friction, pressure, sponge, w/d, forcing, nesting,
solver integration, mode splitting, multilayer kernels, gauges,
unstructured analogues.

## 8. Runtime Configuration

All runtime parameters via Fortran namelist. Each namelist block is
optional — missing blocks fall back to defaults.

### `&sim_nml`

| Parameter | Default | Description |
|---|---|---|
| `sim_type` | `"coastal"` | `"coastal"` (HLL/HLLC + KNP) |

### `&grid_nml`

| Parameter | Default | Description |
|---|---|---|
| `nx`, `ny` | 200, 1 | Physical cells in x/y |
| `dx`, `dy` | 1.0, 1.0 | Cell size (m) |
| `nghost` | 3 | Ghost cells per side |

### `&time_nml`

| Parameter | Default | Description |
|---|---|---|
| `t_end` | 1.0 | Simulation end time (s) |
| `cfl` | 0.45 | CFL number |
| `dt_max` | 1e10 | Maximum timestep (s) |
| `dt_fixed` | 0.0 | Fixed timestep; 0 = adaptive CFL |
| `cfl_interval` | 1 | Recompute CFL every N steps |

### `&physics_nml`

| Parameter | Default | Description |
|---|---|---|
| `manning_n` | 0.0 | Manning's roughness coefficient |
| `dry_threshold` | 0.001 | Dry depth threshold (m) |
| `wind_stress_x` | 0.0 | Wind stress x-component (Pa) |
| `wind_stress_y` | 0.0 | Wind stress y-component (Pa) |
| `coriolis_f` | 0.0 | Coriolis parameter f (1/s) |
| `use_hllc` | `.false.` | Use HLLC instead of HLL for structured barotropic. **Structured-only**. |

### `&output_nml`

| Parameter | Default | Description |
|---|---|---|
| `output_to_file` | `.false.` | Enable file output |
| `output_dir` | `"./output"` | Output directory |
| `output_interval` | 0.0 | Time between snapshots (s) |
| `output_format` | `"netcdf"` | `"netcdf"` or `"ascii"` |
| `restart_interval` | 0.0 | Restart file interval (s) |
| `compress_output` | `.false.` | Enable deflate compression |
| `compress_level` | 1 | Deflate level (1–9) |
| `bathymetry_file` | `""` | NetCDF bathymetry file |
| `restart_file` | `""` | Warm start from restart file |
| `forcing_file` | `""` | Atmospheric forcing NetCDF file |
| `forcing_x_origin`, `forcing_y_origin` | 0.0 | Forcing grid offset (m) |
| `gauge_file` | `""` | Path to gauge locations text file |
| `gauge_output_interval` | 0.0 | Gauge output interval (s; 0 = every step) |

### `&boundary_nml`

| Parameter | Default | Description |
|---|---|---|
| `bc_west / east / south / north` | `"wall"` | BC type per edge |
| `n_tidal_constituents` | 0 | Number of tidal constituents |
| `tidal_amp(10)` | 0.0 | Constituent amplitudes (m) |
| `tidal_phase(10)` | 0.0 | Constituent phases (rad) |
| `tidal_omega(10)` | 0.0 | Constituent frequencies (rad/s) |
| `inflow_velocity` | 0.0 | Normal velocity for inflow BC (m/s) |
| `inflow_salinity` | -1.0 | Inflow salinity (PSU); <0 → zero-gradient |
| `inflow_temperature` | -999.0 | Inflow temperature (°C); <0 → zero-gradient |
| `discharge_flux` | 0.0 | Volume flux for discharge BC (m²/s) |
| `clamped_elevation` | 0.0 | Dirichlet elevation for clamped BC (m) |
| `sponge_width` | 0 | Sponge layer width (cells) |
| `sponge_strength` | 0.0 | Sponge relaxation rate (1/s) |
| `sponge_h_ref` | 0.0 | Sponge reference depth (m) |
| `nesting_file` | `""` | Parent domain output file |
| `nesting_x_origin`, `nesting_y_origin` | 0.0 | Child origin in parent coords (m) |

### `&nonhydrostatic_nml`

| Parameter | Default | Description |
|---|---|---|
| `use_nonhydrostatic` | `.false.` | Enable NH pressure correction |
| `use_multilayer` | `.false.` | Enable hydrostatic multilayer |
| `nz_layers` | 2 | Number of sigma layers |
| `bpg_method` | `"jacobian"` | `"sigma"` or `"jacobian"` |
| `ml_min_depth` | 0.5 | Minimum depth for vertical operations (m) |
| `nh_tolerance` | 1e-6 | CG solver tolerance (NH only) |
| `nh_max_iterations` | 100 | CG max iterations (NH only) |
| `nh_preconditioner` | `.false.` | Vertical-line preconditioner for CG |
| `nh_solve_interval` | 1 | Solve NH pressure every N steps |
| `split_ratio` | 1 | 0=auto, 1=unsplit, N>1=fixed |
| `split_averaging` | `.true.` | ROMS-style time averaging |
| `rho_0` | 1000.0 | Reference density for EOS (kg/m³) |
| `vmix_nu_max` | 0.01 | Max turbulent viscosity (m²/s) |
| `vmix_nu_bg` | 1.0e-4 | Background vertical viscosity (m²/s) |
| `vmix_ri_alpha` | 5.0 | Ri scaling coefficient (PP81) |
| `hdiff_kappa` | 0.0 | Horizontal tracer diffusion (m²/s) |
| `smagorinsky_cs` | 0.0 | Smagorinsky coefficient (typical 0.1–0.2) |
| `out_h_layer` | `.true.` | Output layer thickness |
| `out_velocity` | `.true.` | Output `u_ml`, `v_ml` |
| `out_salinity` | `.true.` | Output salinity |
| `out_temperature` | `.true.` | Output temperature |
| `out_density` | `.false.` | Output density |
| `out_w` | `.false.` | Output vertical velocity at interfaces |

### `&tracer_nml`

Tracer initial conditions, EOS coefficients, clamp bounds, background
vertical diffusivities. Running without this block is fine — all
fields have sensible defaults.

| Parameter | Default | Description |
|---|---|---|
| `initial_salinity` | 35.0 | Initial salinity (PSU) |
| `S_ref` | 0.0 | Reference salinity in EOS (PSU) |
| `beta_S` | 0.78 | Haline contraction (kg/m³ per PSU) |
| `S_min` | 0.0 | Lower clamp (PSU) |
| `S_max` | 40.0 | Upper clamp (PSU) |
| `kappa_S_bg` | 1.0e-5 | Background vertical S diffusivity (m²/s) |
| `initial_temperature` | 15.0 | Initial potential temperature (°C) |
| `T_ref` | 15.0 | Reference temperature in EOS (°C) |
| `alpha_T` | 0.17 | Thermal expansion (kg/m³ per °C); 0 = passive |
| `T_min` | -2.0 | Lower clamp (°C); seawater freezing |
| `T_max` | 40.0 | Upper clamp (°C) |
| `kappa_T_bg` | 1.0e-5 | Background vertical T diffusivity (m²/s) |
| `eos_form` | `"linear"` | Coastal EOS selector, both backends: `linear` (legacy linear kernel — bit-identical default), `wright` (Wright 1997 rational) or `roquet_spv` (Roquet et al. 2015 TEOS-10 polynomial) routed through the shared `rki_eos`, evaluated at p = 0 (surface/potential density). Non-linear forms ignore `alpha_T`/`beta_S`/`S_ref`/`T_ref` in the density computation; the KPP surface buoyancy-flux estimate `B_0` still reads `alpha_T`/`beta_S` as constant expansion coefficients — set them to locally representative values when combining with `vmix_use_kpp` |

### `&vcoord_nml`

| Parameter | Default | Description |
|---|---|---|
| `vcoord_type` | `"sigma"` | `sigma`, `zsigma`, `zstar`, `zstar_sigma`, `zstar_full`, `z_fixed`, `rho`, `hycom`, `z`, `eulerian_z`, `isopycnal`, `lagrangian`, `gprime` |
| `remap_method` | `"ppm"` | `pcm`, `plm`, `ppm`, `ppm_h4`, `pqm` (conservative remap reconstruction; `pqm` = piecewise-quartic, 4th–5th order, W&A 2008, falls back to PPM for nz<5) |
| `vcoord_depth_transition` | 100.0 | z↔sigma blend depth (m) |
| `vcoord_blend_width` | 50.0 | z↔sigma blend width (m) |
| `zstar_h_surf_target` | 0.0 (auto) | Surface-anchored layer thickness for `zstar_full` (m) |
| `zstar_h_min` | 1.0e-4 | Minimum layer thickness in `zstar_full` (m); below this, the layer vanishes |
| `zstar_stretching` | `"log"` | `"log"` or `"uniform"` |
| `zstar_n_surf` | 0 (auto) | Number of surface-anchored layers (`zstar_full`) |
| `rho_ref_pressure` | 2.0e7 | Reference pressure (Pa) for the potential density that defines `rho` |
| `rho_target_light` | 1020.0 | Lightest (surface) target interface density (kg/m³) for the `rho` light→dense linspace fallback |
| `rho_target_dense` | 1030.0 | Densest (bed) target interface density (kg/m³) for the `rho` linspace fallback |

### `&initial_condition_nml`

| Parameter | Default | Description |
|---|---|---|
| `init_condition` | `"dam_break"` | `dam_break`, `gaussian_hump`, `castro_dam_break`, `still_water` |
| `eta_initial` | 0.0 | Initial elevation for `still_water` (m) |
| `h0` | 1.0 | Background depth for `gaussian_hump` (m) |
| `hump_amplitude` | 0.05 | Perturbation amplitude (m) |
| `hump_sigma` | 5.0 | Gaussian width (m) |
| `hump_x_offset`, `hump_y_offset` | 0.0 | Centre offset (m) |
| `init_u_shear` | 0.0 | Vertical shear (m/s) |

### `&metadata_nml`

| Parameter | Default | Description |
|---|---|---|
| `reference_time` | `""` | CF time origin (e.g. `"2026-01-01 00:00:00"`) |
| `institution` | `""` | Global NetCDF attribute |

### `&logging_nml`

| Parameter | Default | Description |
|---|---|---|
| `log_level` | `"info"` | `debug`, `verbose`, `info`, `performance`, `warning`, `error` |
| `status_interval` | 0.0 | Status print interval (sim seconds) |

### Example

```fortran
&sim_nml
  sim_type = "coastal"
/

&grid_nml
  nx = 1000
  ny = 500
  dx = 100.0
  dy = 100.0
  nghost = 2
/

&time_nml
  t_end = 86400.0
  cfl = 0.45
  dt_max = 10.0
/

&physics_nml
  manning_n = 0.025
  dry_threshold = 0.001
/

&output_nml
  output_to_file = .true.
  output_dir = "./output"
  output_interval = 3600.0
  compress_output = .true.
  compress_level = 1
  bathymetry_file = "bathymetry.nc"
  gauge_file = "gauges.txt"
  gauge_output_interval = 60.0
/

&boundary_nml
  bc_west  = "wall"
  bc_east  = "tidal"
  bc_south = "wall"
  bc_north = "wall"
  n_tidal_constituents = 5
  tidal_amp   = 0.50, 0.12, 0.16, 0.09, 0.11
  tidal_phase = -0.91, -0.79, -0.57, -1.10, -1.14
  tidal_omega = 1.405189e-4, 1.454441e-4, 7.292117e-5, 6.759774e-5, 1.378797e-4
/

&logging_nml
  log_level = "info"
  status_interval = 3600.0
/
```

## 9. Key Design Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Language | Fortran | Numerical performance, scientific community familiarity |
| GPU strategy | `do concurrent` + OpenACC | Single source for GPU, multicore, serial; no vendor lock-in |
| Reductions | `!$acc parallel loop reduction(...)` | Inert comment on non-OpenACC compilers; correct sequential fallback |
| GPU data | `!$acc enter/exit data`, `!$acc update` | Portable across OpenACC-capable compilers; inert elsewhere |
| Compilers | NVHPC, gfortran, ifx, Cray ftn | Mixed-hardware HPC portability |
| Spatial scheme | Finite volume, HLL/HLLC + KNP | Conservative, robust at shocks, GPU-friendly |
| Time scheme | Explicit SSP-RK2 | No global solves; embarrassingly parallel per stage |
| Data layout | SoA, `(nx, ny[, nz])` | Coalesced GPU memory access on first index |
| Grid types | Structured + unstructured side-by-side | Structured for open water (peak throughput); unstructured for coastlines (geometric flexibility) |
| Unstructured GPU pattern | Cell-centric double-visit | Avoids atomics; safe for `do concurrent` |
| Parallelism | Single process; GPU or multicore via `do concurrent` | No MPI |
| I/O strategy | One NetCDF file per run | |
| Python/Julia API | C bindings via `iso_c_binding` + opaque handle | Universal ABI; zero-copy array access; multiple handles per process |
| Build system | CMake | Multi-configuration support, industry standard |
| Configuration | Fortran namelist | Zero dependencies, native to Fortran |
| Dependencies | `pic`, `test-drive`, NetCDF-Fortran | Reuse tested types, logging, testing |
