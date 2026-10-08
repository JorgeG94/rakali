# Executable plan — Wave coupling for the Rakali coastal solver

Status: proposal / plan of record candidate. Scope: this repo is the coastal path
and nothing else (`sim_type='coastal'`; structured Cartesian A-grid + unstructured
triangular; explicit Godunov + Casulli semi-implicit; no ocean dyn-core, no ice).
Verified state (re-confirmed 2026-08-07): rakali has **zero** wave coupling — a
repo-wide grep of `src/` for radiation-stress / rad_stress / Stokes / vortex-force /
WWM / WW3 / SWAN returns **no hits at all**, and
`docs/CAPABILITIES_AND_LIMITATIONS.md:133` lists "Wave radiation stress" under
"Not in either". A clean blank slate, not a partial implementation to extend.

**SI caveat (2026-08-07):** the semi-implicit path is single-rank, sigma-only in
multilayer, and restricted to wall/tidal/clamped boundaries — which rules out the
open-boundary configs most wave-forcing validation needs. See
`schism_parity_report.md` §5 before committing to the SI cells of the four-path
table in §9.

Provenance rule obeyed throughout: SCHISM (`/home/jorge/nci/cdx/schism`,
Apache-2.0) was read only as a *behavioural oracle* for which terms a mature
implementation applies where; **no SCHISM code is reproduced**, and every
algorithm below is cited to primary literature for clean-room re-derivation.
File paths into the SCHISM tree appear only as navigation pointers.

---

## 1. Executive summary + recommended architecture

**Do not port a spectral wave model.** A third-generation action-balance solver
(WWMIII is ~80 kLOC) is a global *implicit spectral-space* solve — a poor
`do concurrent` target and a multi-year build. The value rakali must add is the
**per-cell / per-column wave-forcing source terms**, which are an excellent GPU
fit and map cleanly onto the existing wind-stress / bottom-friction / Coriolis
kernel seams. Split the problem: (a) *couple* to an external wave model for the
wave field; (b) *write* the forcing terms inside rakali.

### Recommended coupling mechanism: model-agnostic NetCDF field ingest (one-way, offline first)

Build a **wave-field state slot + a NetCDF reader** that consumes a
pre-regridded wave field on the rakali grid — significant wave height `Hs`, peak
/ mean period `Tp`/`Tm`, mean direction `θ`, wavenumber `k` (or it is derived
from `Tm` + local depth via the dispersion relation), near-bottom orbital
velocity amplitude `u_orb` + representative period `T_bot`, and wave energy
dissipation `ε_br` (breaking) + `ε_bot` (bottom). This reader reuses the
established forcing-ingest pattern (the ocean `&ocean_dataovr_nml` PR-14 reader;
coastal `rki_forcing.F90` `forcing_update`). One-way, file-driven coupling is
the keystone: it is **immediately useful with offline SWAN/WW3/WWMIII output**,
it is **licence-clean** (rakali never links the wave model — it only reads a
NetCDF the wave model already writes), and it is **model-agnostic**. Online
two-way in-memory coupling is a later, optional phase (§3, Phase 6).

**External-model assessment** (rakali reads their output; never links them):

| Model | Licence | Grid | Fit / recommendation |
|---|---|---|---|
| **WaveWatch III** (NOAA) | Permissive (NOAA WW3 licence; use/redistribute OK) | structured + unstructured-triangle (PDLIB) | **Recommended primary for shelf/regional.** Matches both rakali backends; active GPU efforts; standard NetCDF output. |
| **SWAN** (TU Delft) | GPL-2 | structured + unstructured | **Recommended for surf-zone / nearshore.** GPL is fine for *file* exchange (arm's-length, no linking); do NOT link or vendor. |
| **WWMIII** (in SCHISM tree) | Apache-2.0 | unstructured-native | Usable + it is our behavioural oracle, but entangled in the SCHISM build. Best kept as reference; couple via file if used. |

The ingest schema is the same for all three, so the producer is a deployment
choice, not a code fork.

### Recommended forcing formulation: vortex force (for 3D), radiation stress (for 2D)

**Recommend the vortex-force framework (McWilliams, Restrepo & Lane 2004;
Bennis, Ardhuin & Dumas 2011) as the formulation of record for the 3D
multilayer paths, with the depth-integrated radiation-stress form
(Longuet-Higgins & Stewart 1964) shipped first as the validated 2D-barotropic
stepping stone.** Justification for this architecture specifically:

- The radiation-stress tensor `S_ij` is a **depth-integrated** momentum flux.
  In a 2D barotropic model it is exact and its divergence is a clean single body
  force — cheap, and it validates against the classic closed-form wave
  setup/setdown and longshore-current oracles. So it is the right *first* thing.
- But rakali's 3D paths are **sigma / z\*-lite multilayer**. Distributing a
  depth-integrated `∂S_ij/∂x_j` over sigma layers is ad hoc (there is no
  first-principles vertical shape for it). The **vortex force is
  vertically well-posed**: the Stokes drift `u^{St}(z)` is a physical 3D
  velocity, and the vortex force `−u^{St}×(∇×u)`, the Coriolis-Stokes force
  `−f×u^{St}`, and the Bernoulli-head gradient are all *local* per-layer terms.
  This is exactly why SCHISM offers `RADFLAG=VOR` (Bennis 2011) for 3D and keeps
  `RADFLAG=LON` (Longuet-Higgins) as the depth-integrated option
  (nav: `src/WWMIII/wwm_coupl_selfe.F90`).
- The two are consistent, not rival: in the absence of dissipation the
  vortex-force + Bernoulli terms depth-integrate to the radiation-stress
  divergence, giving a built-in cross-check (Phase 3 test).

So: 2D barotropic gets radiation stress (Phase 1); 3D gets Stokes drift (Phase 2)
+ vortex force / Coriolis-Stokes / Bernoulli (Phase 3). Both consume the same
wave field.

### The one architectural invariant across all four solver paths

Wave momentum forcing is an **explicit, known body force** every step. It
belongs in the **explicit RHS / velocity-predictor accumulator**, exactly where
wind stress already lives — **never** in the semi-implicit Helmholtz operator
(the operator is the gravity-wave term only, and must stay SPD and
forcing-free). The sole exception is **wave-enhanced bottom stress**, which — on
the 3D semi-implicit paths — belongs *inside* the implicit vertical tridiagonal
drag coefficient, again exactly where current-only bottom drag already lives.
This is proven true in-source for all four paths (see §4).

---

## 2. Dependency graph

```
Phase 0  Wave-field slot + NetCDF ingest  (the coupling mechanism; inert)
   │
   ├────────────┬───────────────┬──────────────────┐
   ▼            ▼               ▼                  ▼
Phase 1     Phase 2         Phase 4            Phase 5
Rad-stress  Stokes drift    Wave breaking      Wave-current
(2D baro,   (3D diagnostic  momentum + TKE     bottom stress
 both grids, velocity,       injection         (GM79 / Soulsby97,
 expl+SI)    both grids)     (both grids)       all 4 paths)
   ┊            │
   ┊(consistency▼
   ┊ check)  Phase 3
   └┈┈┈┈┈┈┈┈▶ Vortex force + Coriolis-Stokes + Bernoulli
              (3D, both grids, expl+SI)
                     │
                     ▼
              Phase 6 (optional, deferred)
              Online two-way in-memory coupling
```

Phase 0 gates everything. After it, Phases 1, 2, 4, 5 are independent and can be
built in parallel; Phase 3 needs Phase 2 (Stokes) and cross-checks Phase 1.
Phase 6 is deferred.

---

## 3. Per-phase detail

Conventions used below (verified against source):
- **Wind hook = the template.** Structured surface-concentrated wind:
  `ml_apply_wind_to_surface_layer` (`src/core/coastal/kernels/structured/
  multilayer/rki_ml_dynamics.F90:225`, a `pure` per-cell DC kernel adding
  `dt*wind/RHO_WATER` to `hu_layer(:,:,nz)`). Depth-uniform barotropic wind:
  `apply_friction` (`.../barotropic/rki_kernel_friction.F90`). Unstructured wind
  is depth-uniform inside the fused kernel
  (`.../unstructured/barotropic/rki_kernel_fused_unstr.F90:95-97`).
- **Device residency.** Wave arrays follow the `wind_x` plumbing exactly:
  allocate in the slot `init`; map with `!$acc enter data copyin(...)` in the
  solver `enter_data` (structured template `rki_solver.F90:168`; unstructured
  `rki_surface_forcing_unstr.F90:102`); refresh mid-run with `!$acc update
  device(...)` (`rki_solver.F90:493`, `rki_solver_unstr.F90:678`); `delete` in
  `solver_exit_data` (`rki_solver.F90:340`). The build is `-gpu=mem:separate`
  (no managed memory) — **every** array a kernel touches must be device-present;
  all directives are inert no-ops on host/multicore builds so add them
  unconditionally, and verify on the actual GPU build.
- **Race rule.** Unstructured kernels are cell-centric double-visit — a kernel
  writes only its own cell, never a neighbour. Any gradient/curl (radiation
  stress divergence, vortex force) must be computed by gathering from neighbours
  into the owning cell, reusing the Green-Gauss gradient machinery already used
  by the unstructured BPG — no scatter, no atomics.

---

### Phase 0 — Wave-field state slot + NetCDF ingest (THE coupling mechanism)

**Delivers.** A `wave_field_t` derived type composed onto `state` (mirroring
`surface_forcing_t`, `src/core/structured/rki_surface_forcing.F90:19`), a
matching `unstr_wave_field_t`, and a NetCDF reader that fills them from a
pre-regridded wave file. Device-mapped. **Physics: none yet — fully inert.**
Fields (per-cell): `hs`, `tm`, `dir` (or `dir_x/dir_y` unit components to avoid
angle wraparound on the GPU), `wavenumber` (read or derived from `tm`+depth via
the linear dispersion relation), `u_orb`, `t_bot`, `eps_break`, `eps_bot`.

**Physics / citation.** Linear (Airy) wave theory dispersion `σ² = g k tanh(kH)`
for deriving `k` — **Airy (1841)** / standard (any wave-mechanics text, e.g.
Dean & Dalrymple 1991). No forcing applied.

**Hooks (all four paths).** None into the momentum equations yet. New slot wired
into: the state `init`/`destroy`/`bytes`; the solver `enter_data`/`exit_data`;
a new `wave_update` refresh routine (mirror `solver_update_forcing`
`rki_solver.F90` @493 and `solver_unstr_update_forcing` `rki_solver_unstr.F90`
:666). Reader modelled on `rki_forcing.F90 forcing_update` (@193) and the ocean
`&ocean_dataovr_nml` PR-14 reader.

**GPU strategy.** Pure host-side allocate + fill on read (setup path, no DC);
`enter data copyin` after fill; `update device` on each new wave-file time slice.
No persistent kernel workspace. Arrays are `(nx,ny)` structured / `(ncells)`
unstructured, `wp`.

**Namelist knob.** `&wave_nml enable = .false.` (master, default OFF ⇒
bit-identical) + `wave_file`, `wave_time_mode`. Per-term knobs land in later
phases, all default OFF.

**Test.** Read a two-time-slice synthetic wave NetCDF; assert fields land on
device (run-twice-and-compare pattern for mapped arrays, per the workspace-test
rule); assert a full solver step with `enable=.true.` but all forcing sub-knobs
off is **bit-identical** to `enable=.false.` Derived-`k` unit test against the
dispersion relation at three depths (deep/intermediate/shallow), closed form.

**Size.** Medium (~500 LOC: two types + reader + wiring). **Risk.** Low —
pure plumbing, no numerics in the momentum path. The regridding onto the rakali
grid is assumed offline (see §5, open question).

---

### Phase 1 — Depth-integrated radiation stress (2D barotropic; both grids; explicit + SI)

**Delivers.** Wave setup/setdown and radiation-stress-driven currents on the
**2D barotropic** path. The first *validated* wave forcing.

**Physics / citation.** **Longuet-Higgins & Stewart (1964)** radiation stress
tensor from `Hs`, `k`, `dir`, `H`:
`S_xx = E[n(cos²θ+1) − ½]`, `S_yy = E[n(sin²θ+1) − ½]`,
`S_xy = E n cosθ sinθ`, with `E = (1/16)ρ g Hs²`, `n = ½(1 + 2kH/sinh2kH)`.
Depth-integrated body force on the barotropic momentum:
`F_x = −(1/ρ)(∂S_xx/∂x + ∂S_xy/∂y)`, likewise `F_y`. Longshore-current closure
follows **Longuet-Higgins (1970)**.

**Hooks.**
- **Structured explicit:** new `pure` kernel `ml_apply_radiation_stress` computes
  `S_ij` per cell then a centred-difference divergence (like BPG), adding
  `dt*F/ (ρ·h)`-style increment to the barotropic momentum. Insert as a new call
  co-located with `apply_friction` in `solver_multilayer_rk_stage`
  (`rki_solver.F90` @1935) — i.e. a depth-uniform barotropic body force. At
  `nz=1` it folds into the same stage.
- **Unstructured explicit:** mirror kernel using Green-Gauss gradients (gather
  from neighbours, own-cell write). Insert alongside the fused barotropic update
  (`rki_solver_unstr.F90` @237/@363) or as a standalone barotropic body-force
  call after it.
- **Structured SI (2D):** add `dt*F` to the **face predictor** `hu_star/hv_star`
  right where wind is added (`rki_si_step.F90:262-278`), gated on `Hf>0`. RHS
  assembler `si_compute_eta_rhs` and Helmholtz operator are **untouched**.
- **Unstructured SI (2D):** project `F` onto the edge normal and add to `u_star`
  in the edge predictor (`rki_si_unstr_step.F90:206-210`), gated on
  `Hf_edge>0`. Operator/CG untouched.

**GPU strategy.** Two-pass: (1) per-cell `S_ij` into a small device-resident
scratch slot on the wave state (3 arrays); (2) divergence read from neighbours.
`do concurrent(j,i)` structured, `do concurrent(cell)` unstructured. `S_ij`
scratch mapped in `enter_data` (recompute-on-device, no host prefill). No
managed memory ⇒ scratch must be device-present.

**Namelist knob.** `&wave_nml rad_stress = .false.`.

**Test.** **Wave setup/setdown across a plane beach** — analytical
Longuet-Higgins & Stewart (1964): setdown reaching a minimum at the breakpoint,
then setup with slope `dη/dx = −(1/(ρg(H+η))) dS_xx/dx`; assert the modelled
`η(x)` against the closed form (few-% L2). **Longshore current on a plane beach**
— Longuet-Higgins (1970) closed-form `V(x)` profile under obliquely incident
breaking waves. Both are canonical surf-zone oracles.

**Size.** Medium-large (~700 LOC across 4 paths + 2 tests). **Risk.** Medium —
the divergence stencil at wet/dry and open-boundary edges needs the same care as
BPG (zero the force at dry-neighbour faces); breakpoint location depends on the
wave model's `Hs` already carrying breaking, so garbage-in if the wave field is
unbroken.

---

### Phase 2 — Stokes drift (3D multilayer; both grids)

**Delivers.** A vertically-resolved Stokes drift velocity `u^{St}(z)` diagnosed
per layer from the wave field. A *diagnostic* velocity (no momentum change yet)
consumed by Phase 3 and by tracer transport; it also feeds the Coriolis-Stokes
term in Phase 3.

**Physics / citation.** Monochromatic Stokes profile **Stokes (1847)**:
`u^{St}(z) = (σ a² k) cosh(2k(z+H)) / (2 sinh²(kH)) · (dir)`, with `a=Hs/2`; for
a broad spectrum the **Breivik, Janssen & Bidlot (2014)** monotonic
approximation from `Hs`, `Tm` and surface Stokes `u^{St}_0`. Depth-integrated
Stokes transport `∫u^{St}dz = E/(ρ c)`.

**Hooks.** New `pure` per-column kernel `ml_stokes_drift` writing
`wave%u_stokes(:,:,k)`, `v_stokes` (new device-resident 3D slot on the wave
state). Computed once per step after the vcoord layer geometry is known (after
`ml_distribute_*`), before Phase 3 consumes it. Optionally added to the tracer
**advecting** velocity so scalars feel Lagrangian (bulk) transport — an opt-in
sub-knob; default keeps tracer transport Eulerian (bit-identity). Same insertion
point in both explicit drivers (`rki_solver.F90`, `rki_solver_unstr.F90`) and
both SI-ML drivers (`rki_si_step_ml.F90` Step 2, `rki_si_unstr_step_ml.F90`
Step 2).

**GPU strategy.** `do concurrent(j,i)` (structured) / `do concurrent(cell)`
(unstructured) with an inner `k` loop over layers (contiguous index last:
`(k,j,i)` intent honoured via array shape). New 3D slots `u_stokes/v_stokes`
mapped in `enter_data`. Pure column-local ⇒ no race, no workspace.

**Namelist knob.** `&wave_nml stokes = .false.` (+ `stokes_in_advection`,
default off).

**Test.** **Stokes drift profile** — deep-water exponential decay `u^{St}(z) =
u^{St}_0 exp(2kz)`; assert the per-layer profile against the closed form.
**Depth-integrated Stokes transport** equals `E/(ρc)` to few-%. Both analytical.

**Size.** Medium (~450 LOC). **Risk.** Low-medium — vertical mapping under
sigma/z\* must sample `z` at layer mid-depths consistently with the
bottom-up (k=1 bed) convention; a flipped index silently gives a bed-intensified
"Stokes" profile.

---

### Phase 3 — Vortex force + Coriolis-Stokes + Bernoulli head (3D; both grids; explicit + SI)

**Delivers.** The physically complete 3D wave-current momentum coupling. The
formulation of record for 3D.

**Physics / citation.** **McWilliams, Restrepo & Lane (2004)** / **Bennis,
Ardhuin & Dumas (2011)** vortex-force form: per-layer momentum tendency
`−u^{St}×(∇×u)` (vortex force) `− f×u^{St}` (Coriolis-Stokes) `− ∇J` (Bernoulli
head / wave-induced pressure `J`). Uses `u^{St}` from Phase 2, the relative
vorticity `∇×u` (reuse existing curl/vorticity machinery), and `f` from
`coriolis_f`. SCHISM's `wwave_force(2,nvrt,nsa)` applied at sides across levels
is the behavioural analogue (nav: `src/Hydro/schism_step.F90:875`,
`src/WWMIII/wwm_coupl_selfe.F90:98` RADFLAG=VOR).

**Hooks.** New `pure` per-layer kernel `ml_apply_vortex_force`, slotted **next to
per-layer Coriolis** (which is where a rotational body force belongs):
- **Structured explicit:** call at `rki_solver.F90` @1912 region (right after
  `ml_apply_coriolis`, before recompute-barotropic @1929) — adds a per-layer
  `hu_layer/hv_layer` increment.
- **Unstructured explicit:** call after `ml_apply_coriolis_unstr`
  (`rki_solver_unstr.F90:352` / :492), before recompute (@353/@493).
- **Structured SI-3D:** add as an explicit **layer tendency** in the Step-2
  predictor block (`rki_si_step_ml.F90:191-272`, alongside BPG/Coriolis/wind at
  @267), **before** `si_column_solve` @360. The column solve then integrates it
  into `Gbar`; operator/CG untouched.
- **Unstructured SI-3D:** same, in the Step-2 predictor
  (`rki_si_unstr_step_ml.F90:150-168`), before `si_unstr_column_solve` @239;
  carried automatically by the `Gbar − X0` increment
  (`rki_si_unstr_column.F90:366` docstring: "drag/wind/BPG enter through the
  increment").

**GPU strategy.** Per-layer column-local plus a relative-vorticity read from
neighbours (gather, own-cell write; reuse the vorticity stencil). `do concurrent`
over cells with inner-k. `u_stokes/v_stokes` and the Bernoulli scratch mapped
device-resident. Watch the NaN-launder gotcha on any clamp of the vortex-force
magnitude — guard with `ieee_is_finite`.

**Namelist knob.** `&wave_nml vortex_force = .false.` (+ sub-toggles
`coriolis_stokes`, `bernoulli` for staged validation).

**Test.** **Stokes-Coriolis equilibrium** — with no waves-breaking and steady
`u^{St}`, the Coriolis-Stokes force drives an anti-Stokes Eulerian return flow
`u^{E} = −u^{St}` in steady state (closed form; Hasselmann 1970 / Polton et al.
2005). **Consistency check:** with dissipation off, assert the depth-integrated
vortex-force + Bernoulli terms match the Phase-1 radiation-stress divergence
(cross-formulation oracle).

**Size.** Large (~800 LOC across 4 paths + tests). **Risk.** Medium-high — the
Bernoulli head and the vertical vortex-force component need the vertical Stokes
velocity `w^{St}`; getting the vorticity stencil right on the A-grid vs the
central-upwind triangles is the fiddly part. Stokes drift at open boundaries
needs handling (SCHISM subtracts depth-averaged Stokes at OBC, Bennis 2011) —
flagged in §5.

---

### Phase 4 — Wave-breaking momentum + turbulence injection (both grids)

**Delivers.** The non-conservative surface breaking force into momentum, and a
wave-breaking TKE surface flux into the turbulence closures.

**Physics / citation.** Breaking-induced surface momentum flux from the wave
dissipation `ε_br` (surface-roller / breaking stress ~ `ε_br·dir/c`) —
**Deigaard (1993)** / roller framework; **Battjes & Janssen (1978)** for the
breaking dissipation origin. Wave-enhanced near-surface turbulence via a TKE
surface-flux boundary condition proportional to `ε_br` — **Craig & Banner
(1994)**.

**Hooks.**
- **Momentum (all four paths):** surface-concentrated body force at `k=nz` —
  mirror `ml_apply_wind_to_surface_layer` exactly; same insertion points as the
  wind hook (explicit `rki_solver.F90` @1921 / unstr fused; SI Step-2 predictor).
- **TKE injection:** modify the surface boundary condition of the **k-ε**
  ψ/TKE-flux and the **KPP** surface forcing. Structured k-ε friction velocities
  are computed at `rki_solver.F90:1754` (`ml_keps_friction_velocities`); KPP
  boundary layer at `ml_kpp_boundary_layer` (reads `u_*`, `B_0`). Add a
  `ε_br`-derived surface TKE flux term there (both backends;
  `rki_ml_keps{,_unstr}`, `rki_ml_kpp{,_unstr}`). PP81 has no TKE equation ⇒ no
  injection there (documented limitation).

**GPU strategy.** Per-cell surface-layer add (`do concurrent(j,i)` /
`(cell)`), plus a scalar surface-BC modification inside the existing per-column
closure kernels (no new loop shape). No new workspace.

**Namelist knob.** `&wave_nml breaking = .false.` (+ `breaking_tke`, default off).

**Test.** **Craig & Banner (1994)** wave-enhanced surface layer — near-surface
TKE dissipation elevated above the law-of-the-wall value, with the characteristic
`z`-dependence; assert the modelled near-surface `k`/dissipation profile against
the CB94 analytical surface layer. Secondary: breaking deepens the surface mixed
layer vs a no-wave control.

**Size.** Medium (~500 LOC). **Risk.** Medium — the roller/dissipation split is
wave-model-dependent (some models output a roller field, some only bulk `ε_br`);
the TKE surface-flux constant needs calibration and only lands on k-ε/KPP, not
PP81 (the default) — so the default closure sees only the momentum half.

---

### Phase 5 — Wave-enhanced bottom stress / wave-current BBL (all four paths)

**Delivers.** Combined wave-current bottom stress — the enhanced mean bed stress
and apparent roughness under waves, which under-predicting biases sediment
resuspension and bottom drag on wave-exposed coasts.

**Physics / citation.** **Grant & Madsen (1979)** wave-current bottom boundary
layer (nonlinear interaction of the wave orbital velocity `u_orb` with the mean
current giving an enhanced mean bed stress and an apparent roughness), or the
algebraically simpler **Soulsby (1997)** parameterization of the mean and maximum
combined stress `τ_m`, `τ_max` from the current stress `τ_c`, wave stress `τ_w`,
and the wave-current angle. SCHISM uses a modified Grant-Madsen (nav:
`src/Hydro/misc_subs.F90:4523`, apparent-roughness form).

**Hooks.** The current-only Manning drag is
`cf = g·n²/H^{1/3}`, `speed_b = sqrt(u_b²+v_b²)`, `denom = 1 + dt·cf·speed_b/h_k`
— the wave term enters by replacing `speed_b` with a combined wave-current
friction velocity (Soulsby) and/or scaling `cf` by the wave-current enhancement:
- **Structured explicit:** `ml_apply_bottom_friction`
  (`rki_ml_friction.F90:84-97`) — modify the `speed_b`/`cf` at :86/:93; add a
  `u_orb`/`t_bot` per-cell argument plumbed like `wind_x`. Also the barotropic
  twin in `apply_friction` (`rki_kernel_friction.F90:85-91`).
- **Unstructured explicit:** `ml_apply_bottom_friction_unstr`
  (`rki_ml_friction_unstr.F90:70-78`) — same; wave arrays reached via `state`.
  Also the fused-kernel twin (`rki_kernel_fused_unstr.F90:89-91`).
- **Structured SI-3D:** the drag is **implicit** — modify `lambda_drag` in the
  column solve (`rki_si_column.F90:264-279`), which folds into the bed-row
  diagonal `bk_col(1) = 1 + alph + dt*lambda_drag` (@294). Stays SPD.
- **Unstructured SI-3D:** `lambda_drag` in `si_unstr_column_solve`
  (`rki_si_unstr_column.F90:231-242`), bed diagonal @258.
- (2D SI paths have no bottom-drag term today; adding a wave-current *linear*
  bed drag to the 2D SI predictor is optional and out of this phase's core.)

**GPU strategy.** Purely local per-cell/per-column modification of an existing
`pure` kernel — the ideal DC case, no new loop. New per-cell `u_orb`, `t_bot`
arrays device-resident (wind_x pattern). The combined-stress algebra must be
`ieee_is_finite`-guarded (no NaN launder under `-fast`).

**Namelist knob.** `&wave_nml wave_bbl = .false.` (+ `bbl_form = "soulsby" |
"grant_madsen"`).

**Test.** **Soulsby (1997)** combined wave-current stress — assert the modelled
mean/max bottom stress against the parameterization's own tabulated closed-form
`τ_m/τ_c`, `τ_max/τ_c` for a given wave/current speed ratio and angle
(the reference gives the exact algebraic answer). Secondary: a steady current's
bed stress increases monotonically with imposed `u_orb`.

**Size.** Medium (~450 LOC across 4 paths). **Risk.** Low-medium — mostly local;
the Grant-Madsen apparent-roughness iteration (if that form is chosen over
Soulsby's non-iterative algebra) needs a fixed-iteration cap sized from an
upfront reduction, not an unbounded loop.

---

### Phase 6 — Online two-way in-memory coupling (DEFERRED / optional)

**Delivers.** Live field exchange (rakali → wave model: `η`, currents, water
level; wave model → rakali: the Phase-0 field) without file round-trips, for
feedback-sensitive surf-zone cases. **Deferred** — the file-based one-way path
(Phases 0-5) is sufficient for storm-surge total-water-level and offline
surf-zone work. Mechanism: exchange through the opaque-handle C API (the ocean
path is API-gated, but coastal is API-exposed) or an ESMF/OASIS-style coupler as
a CPU co-process. Not scheduled here; flagged for a separate plan.

---

## 4. The four-path forcing table

How each forcing term enters each of structured/unstructured × explicit/SI. "RHS
predictor" = the explicit velocity accumulator (never the Helmholtz operator).

| Term | Structured **explicit** | Unstructured **explicit** | Structured **semi-implicit** | Unstructured **semi-implicit** |
|---|---|---|---|---|
| **Radiation stress** (2D baro, Phase 1) | New barotropic body-force kernel by `apply_friction`, `rki_solver.F90:1935`; `F=−(1/ρ)∇·S` centred-diff | Green-Gauss `∇·S`, own-cell write; by fused update `rki_solver_unstr.F90:237/363` | `hu_star/hv_star` **face predictor**, add `dt·F` at `rki_si_step.F90:262-278`, gate `Hf>0`; **operator untouched** | Project `F` onto edge normal → `u_star` edge predictor `rki_si_unstr_step.F90:206-210`, gate `Hf_edge>0`; **operator untouched** |
| **Stokes drift** (3D diag, Phase 2) | Per-column kernel → `u_stokes` slot after `ml_distribute_*`; feeds VF + (opt) advection | Same, own-cell | Diagnosed in Step-2 predictor block `rki_si_step_ml.F90:191-272` | Diagnosed in Step-2 predictor `rki_si_unstr_step_ml.F90:150-168` |
| **Vortex force + Coriolis-Stokes + Bernoulli** (3D, Phase 3) | Per-layer `hu_layer/hv_layer` increment next to `ml_apply_coriolis`, `rki_solver.F90:1912` | Per-layer increment after `ml_apply_coriolis_unstr`, `rki_solver_unstr.F90:352/492` | Explicit **layer tendency** in Step-2 predictor (`rki_si_step_ml.F90:191-272`, by wind @267), **before** `si_column_solve` @360 → rides into `Gbar` | Explicit layer tendency Step-2 (`rki_si_unstr_step_ml.F90:150-168`), before column solve @239 → carried by `Gbar−X0` increment |
| **Wave-breaking momentum** (Phase 4) | Surface-layer add at `k=nz`, mirror `ml_apply_wind_to_surface_layer` @1921 | Surface-layer add (new; unstr wind is depth-uniform today) | Surface-layer tendency in Step-2 predictor | Surface-layer tendency in Step-2 predictor |
| **Breaking TKE** (Phase 4) | k-ε surface BC @`rki_solver.F90:1754`; KPP `ml_kpp_boundary_layer` | k-ε/KPP unstr surface BC | Closures run inside the SI ML column path — same surface-BC edit | Same |
| **Wave-current bottom stress** (Phase 5) | Modify `speed_b`/`cf` in `ml_apply_bottom_friction` `rki_ml_friction.F90:84-97` | Modify `ml_apply_bottom_friction_unstr` `rki_ml_friction_unstr.F90:70-78` | **Implicit** `lambda_drag` → bed diagonal `rki_si_column.F90:264-294` | **Implicit** `lambda_drag` → bed diagonal `rki_si_unstr_column.F90:231-258` |

**Where the treatment must differ (the load-bearing distinction):**
- **Explicit vs SI momentum forcing.** In the explicit Godunov update, a body
  force is a direct increment on `hu/hv[_layer]` after the flux update. In the
  SI update it must go into the **explicit velocity predictor** (`hu_star` /
  `u_star` for 2D; the Step-2 per-layer tendency that feeds `Gbar` for 3D) —
  i.e. into the RHS "G" vector — and **must not** enter the Helmholtz operator,
  which stays the pure gravity-wave SPD operator (`si_apply_helmholtz_operator`
  / `si_unstr_apply_helmholtz` read only `Hf` and `η`). Radiation stress and
  vortex force are explicit known terms ⇒ RHS, both regimes.
- **Bottom stress is the exception.** On the explicit paths wave-current bottom
  stress is a semi-implicit divide in the friction kernel (RHS-side). On the SI
  **3D** paths it belongs *inside* the implicit vertical tridiagonal (the
  `lambda_drag` bed-diagonal coefficient) — modifying the operator, not the RHS —
  exactly as current-only drag already does. This is the one term whose SI
  treatment is operator-side.
- **Structured vs unstructured gradients.** Divergence/curl terms use
  centred differences on the structured A-grid but Green-Gauss gather on the
  triangles (own-cell write only, no neighbour scatter — the race rule).
- **Unstructured explicit has no surface-layer wind kernel** (wind is
  depth-uniform in the fused kernel), so Phase 4's surface breaking momentum
  needs a *new* surface-layer kernel there, whereas structured reuses the wind
  template.

---

## 5. Risks and open questions needing a human decision

1. **3D radiation stress vs vortex force — confirm the split.** This plan ships
   radiation stress for 2D (Phase 1) and vortex force for 3D (Phases 2-3),
   because the depth-integrated `S_ij` has no principled vertical distribution on
   sigma layers. If a reviewer wants a single formulation everywhere, that is a
   design change. **Decision: endorse VF-for-3D / rad-stress-for-2D, or mandate
   one framework.**
2. **External wave model + regridding ownership.** Which producer (WW3 for
   shelf/regional, SWAN for surf-zone, WWMIII as reference)? And who does the
   **grid interpolation** from the wave grid to the rakali grid — offline
   pre-regrid (assumed by Phase 0) vs an online regridder? The offline
   assumption is a real constraint on users. **Decision needed.**
3. **Licence boundary (SWAN is GPL-2).** File-based exchange keeps SWAN at
   arm's length (rakali never links it) — this must be signed off. WWMIII is
   Apache-2.0; WW3 permissive. **Confirm the coupling stays file/exchange-only,
   never linked/vendored.**
4. **Stokes drift at open boundaries.** Bennis et al. (2011) / SCHISM subtract
   the depth-averaged Stokes drift at open boundaries to avoid spurious
   recirculation. Phase 3 must decide the OBC treatment for `u^{St}` — this
   interacts with the Flather/Chapman/Orlanski BCs. **Flag for Phase-3 design.**
5. **Breaking dissipation representation.** Roller vs bulk `ε_br` is
   wave-model-dependent, and the TKE injection only reaches k-ε/KPP — the
   **default PP81 closure sees no wave-breaking turbulence**. Acceptable? Or
   should wave coupling flip the default closure recommendation to k-ε?
6. **Two-way coupling cadence & stability** (Phase 6, deferred) — sub-cycling
   ratio and feedback stability are unresolved; deferred deliberately.
7. **Validation-data availability.** The analytical oracles (setup/setdown,
   longshore current, Stokes profile, Soulsby stress) are closed-form and need
   no external data — good. But an end-to-end *coupled* validation (e.g. a
   Duck94-class field case) needs real wave output and is out of the analytical
   test scope.

---

## 6. Explicitly out of scope

- **Porting any spectral wave model** (WWMIII / WW3 / SWAN internals — the
  action-balance spectral solver). Coupling only.
- **Online two-way tight coupling** in v1 (Phase 6 deferred).
- **Phase-resolving wave dynamics** — that is the non-hydrostatic path
  (`ROADMAP_COASTAL.md §12A`), a different capability; this plan is phase-averaged
  wave *forcing*.
- **Wave-driven sediment transport / morphodynamics** — a separate plan; this
  plan delivers the wave-current bottom stress (Phase 5) that a morphodynamics
  effort consumes, but no Exner/bedload here.
- **Infragravity waves, wave setup as a separate elliptic Poisson solve** — setup
  emerges from the momentum body force directly (Phases 1/3), not from a
  dedicated elliptic solver.
- **Ocean-path (`sim_type='ocean'`) wave coupling** — out of scope; the ocean
  EPBL Stokes parameterization there is unrelated.
- **Wave forcing on the barotropic-only 2D SI bottom drag** beyond the optional
  note in Phase 5 (2D SI has no bed-drag term today).
