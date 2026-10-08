# Executable plan — Sediment transport & morphodynamics (coastal path)

Scope: this repo is the coastal path and nothing else (`sim_type='coastal'`; no ocean
dyn-core, no ice). Every capability must land on **all four paths**: {structured,
unstructured} × {explicit Godunov, Casulli semi-implicit}, at 2D barotropic and 3D
sigma/ALE multilayer where the physics is 3D.

**Blocker on the SI half (added 2026-08-07):** `use_sediment` under
`time_scheme='semi_implicit'` with multilayer is a **hard config error**
(`src/core/rki_config.F90:1123`), and the SI multilayer path is sigma-only
(`:1129`) and single-rank (`:1053`). Morphodynamics therefore has *no* semi-implicit
path today, and a moving bed on z\*/ALE has no SI counterpart at all. Either lift
those guards as a prerequisite or declare the two SI cells deferred in every
four-path table below — see `README.md` and `schism_parity_report.md` §5.

Provenance rule (house law): **cite the paper, never the codebase**. SCHISM
(`/home/jorge/nci/cdx/schism`, Apache-2.0) is read only as a *behavioural oracle* for what a
mature implementation sequences and includes; no SCHISM source is transcribed. Every kernel
here is re-derived from the primary literature cited per phase. File paths are navigation aids.

---

## 1. Executive summary

Rakali today has a **single-class suspended-sediment stub**, structured multilayer only:

- `ml_settling_tracer` (`src/tracer/structured/rki_ml_tracers.F90:963`) — conservative PPM
  settling of one suspended tracer into a 2D `bed_sediment(nx,ny)` mass reservoir.
- `ml_erode_bed_sediment` (`:1054`) — Ariathurai & Arulanandan (1978) critical-shear erosion,
  Manning-consistent bed stress, positivity-capped, conservative against the reservoir.
- Wired in `rki_solver.F90:1477` under `use_sediment`, applied **once per step** on the
  RK2-averaged state.

The blocking fact: **`bed_sediment` is a mass reservoir disjoint from the bathymetry
`state%barotropic%b`. The seabed elevation never changes.** There is no Exner equation, no
bedload, no bed layering, no morphodynamic feedback. So the existing piece is a *starting
fragment*, not a foundation for morphology.

**Verdict on the stub: EXTEND, do not replace.** The settling remap (a proven,
GPU-resident, exactly-conservative Colella-Woodward column remap) and the Ariathurai erosion
law are correct and worth keeping. What must change is *ownership and generality*: the reservoir
becomes a per-class **bed-mass inventory** that the Exner update consumes, and the two kernels
move into a dedicated sediment module, generalised over a grain-class registry. Replacing them
would throw away conservation properties already tested (`tests/test_sediment.F90`).

**The keystone item is Phase 3 — the Exner bed evolution.** It is the *only* phase that lets
`b` move, and it is architecturally the hard one: a moving seabed changes the total water depth
`H` that BOTH the explicit Godunov flux (`eta = h + b`, well-balanced reconstruction) and the
implicit η-solve (`si_build_face_depths(h,...)`) depend on, and it re-defines the sigma/ALE
vertical coordinate (which is anchored to the bed) and the wet/dry masks. Every phase after P0–P2
is downstream of it.

**Rough total size:** ~9 PR-sized phases. P3 (keystone) is Large; P1, P2, P4, P5, P8 are
Medium–Large; P0, P6, P7 are Small–Medium. Order-of-magnitude ~8–12 kLOC across all four paths
plus tests — comparable to SCHISM's `Sediment/`+`Sed2d/` combined, but GPU-native and narrower
(no bedform/roughness predictor menu, no wave-asymmetry bedload in v1).

Design invariants (non-negotiable, every phase):
- **Every new knob defaults OFF ⇒ bit-identity** for all existing namelists and tests.
- **GPU:** `do concurrent` per-cell/per-column, persistent workspaces on the state type
  (mapped once via `*_enter_data`), no per-step allocation, `-gpu=mem:separate` (NO managed
  memory) — every array a kernel touches must be device-resident.
- **Conservation:** total sediment mass = suspended (water column) + bed inventory is conserved
  to round-off; morphology moves mass between the bed inventory and `b`, tracked explicitly.
- **A closed-form or standard-benchmark oracle per phase.**

---

## 2. Dependency graph

```
P0  Sediment module + multi-class registry + per-class bed inventory   [foundation]
 │
 ├─► P1  Multi-class suspended transport (settling + erosion/deposition)
 │        │
 │        ▼
 ├─► P2  Bedload transport (MPM / van Rijn 2007 / Soulsby face flux)
 │        │
 │        ▼
 └──────► P3  ★ EXNER BED EVOLUTION (b moves; reconcile h/eta/vcoord/masks) ★  KEYSTONE
                   │        needs P1 (deposition/erosion) + P2 (bedload divergence)
                   │
        ┌──────────┼───────────────┬────────────────┐
        ▼          ▼               ▼                ▼
    P4 Bed       P6 Avalanching   P7 Morfac        (P3 alone already gives
    stratigraphy  (repose-angle    (Roelvink 2006   morphology for the
        │          relaxation)      accel + ramp)   simplest single-layer bed)
        ▼
    P5 Cohesive + consolidation (Partheniades; time-dependent tau_ce)
        │  (consolidation needs the P4 sub-bottom layers)
        ▼
    P8 Wave-current bed stress (Grant-Madsen / Soulsby)
        needs the WAVE PLAN (plan_wave_coupling.md — DOES NOT EXIST YET; flag)
        + upgrades the bed stress feeding P1/P2/P5
```

Critical path: **P0 → P1 → P2 → P3**. Everything else hangs off P3 and can be scheduled in any
order (P5 after P4). P8 is gated on an external deliverable (the wave plan) and should be treated
as provisional until that plan lands.

---

## 3. Per-phase detail

Naming: new module `src/tracer/structured/rki_ml_sediment.F90` +
`src/tracer/unstructured/rki_ml_sediment_unstr.F90` (the settling/erosion kernels migrate out of
`rki_ml_tracers.F90`). Bed-evolution kernels live under a new
`src/core/coastal/kernels/{structured,unstructured}/morphodynamics/`. Namelist group
`&sediment_nml` (one group, per-class arrays), master switch `use_sediment` (already exists).

### P0 — Sediment module + multi-class registry + per-class bed inventory

**Delivers.** A refactor-only PR that (a) moves `ml_settling_tracer` / `ml_erode_bed_sediment`
verbatim into a new `rki_ml_sediment` module, (b) generalises the single `bed_sediment(nx,ny)`
into a per-class bed inventory `bed_sed(nx,ny,nclass)` (structured) / `bed_sed(nclass,ncells)`
(unstructured), and (c) registers `nclass` sediment classes as passive tracers through the
existing `multilayer%register_passive_tracer` seam, each with a settling velocity, density,
diameter, critical-shear and erodibility. With `nclass=1` and default params the run is
**bit-identical** to today (the reservoir is just a rank-3 array with one class).

**Physics + citation.** None new — pure plumbing. Grain-class properties per Soulsby (1997)
(dimensionless grain size `D_* = D[(s-1)g/ν²]^{1/3}`, Shields threshold, settling velocity).

**Hooks (four paths).** Structured/unstructured: the registry already threads every transport
kernel (advection, ALE remap, vdiff, clamp, RK2) — a sediment class rides those for free (this is
the whole point of the existing passive-tracer seam). The settling+erosion call site in
`rki_solver.F90:1477` and its unstructured twin loop over classes. Explicit vs SI: sediment
source terms are applied on the RK2-averaged / SI-updated state *after* the hydrodynamic
sub-step, identically for both time schemes — P0 adds the unstructured settling/erosion kernels
that today are structured-only (the report's "nothing on the unstructured backend" gap), using
the cell-centric double-visit pattern.

**GPU.** `do concurrent(cell) local(...)` with an inner `do c = 1, nclass` class loop (nclass
small, ≤ ~5 — keep it the innermost non-parallel loop). Bed inventory arrays allocated on the
multilayer state, mapped in the existing `*_enter_data`. No new workspaces.

**Knob.** `&sediment_nml n_class`, per-class arrays `d50(:)`, `rho_s(:)`, `w_s(:)`,
`tau_ce(:)`, `e_rate(:)`, `porosity(:)`. Default `n_class=1`, values = current scalars ⇒
bit-identical.

**Test.** `test_sediment` extended: (i) `n_class=1` reproduces the current stub to round-off
(regression gate); (ii) two identical classes carry exactly half the mass each and sum to the
single-class answer (linearity); (iii) unstructured settling conserves suspended+bed to 1e-12
(new-path conservation).

**Size:** Small–Medium. **Risk:** Low (refactor under existing tests + a new-backend port).

### P1 — Multi-class suspended transport (erosion / deposition exchange)

**Delivers.** Per-class erosion and deposition that actually differ by class: class-specific
settling `w_s(c)` (already carried by P0), class-specific pickup. Replaces the single Ariathurai
call with a per-class erosion/deposition exchange, and adds the **van Rijn (1984) reference
concentration** pickup for non-cohesive classes as an alternative to the cohesive
Ariathurai–Krone law (dispatch by class type `cohesive(c)`).

**Physics + citation.**
- Non-cohesive suspended pickup / reference concentration: **van Rijn (1984)** (Part II,
  suspended load) — `c_a = 0.015 · d50/a · T^{1.5}/D_*^{0.3}` at reference height `a`, with
  transport-stage `T = (τ_b - τ_ce)/τ_ce`.
- Settling velocity: **Soulsby (1997)** formula over the full `D_*` range (Stokes → turbulent).
- Cohesive erosion: retain **Ariathurai & Arulanandan (1978)** `E = M(τ_b/τ_ce − 1)_+`;
  deposition **Krone (1962)** `D = w_s·c·(1 − τ_b/τ_cd)_+`.
- Equilibrium suspended profile (test oracle): **Rouse (1937)** `c(z)/c_a =
  [(H−z)/z · a/(H−a)]^{Ro}`, `Ro = w_s/(κ u_*)`.

**Hooks (four paths).** Advection/settling/mixing of each suspended class already flows through
the registry on both backends and under both time schemes (transport is time-scheme-agnostic —
it runs on the post-substep state). The new part is the erosion/deposition **source term**,
computed per cell from the bed shear `τ_b`. For explicit it uses the resolved bottom-layer
velocity as today; for SI it uses the SI-updated bottom-layer velocity — same kernel, called at
the same post-substep point, so no path divergence. Bed shear `τ_b` comes from
`ml_apply_bottom_friction`'s Manning-consistent form in v1 (upgraded in P8).

**GPU.** Per-cell `do concurrent`, inner class loop. Van Rijn / Soulsby closures are pure
`!$acc routine seq` helpers — but per the NVHPC cross-module-no-inline rule, duplicate the hot
ones as `*_impl` bodies inside the kernel module if a profile shows >5% (do not preemptively
duplicate). No new persistent workspace (works on registry + bed inventory, both already mapped).

**Knob.** `&sediment_nml cohesive(:)` (per-class), `tau_cd(:)` (deposition threshold),
`ref_height_a`, `susp_scheme = "vanrijn84" | "ariathurai"`. Default reproduces P0.

**Test.** `test_sediment_rouse` — a 1D vertical column in steady log-law flow, zero net
settling+turbulent-diffusion flux at equilibrium must reproduce the analytical **Rouse profile**
to a few %. Plus a per-class mass-conservation gate (suspended+bed) to 1e-12.

**Size:** Medium. **Risk:** Medium (closure correctness; Rouse needs the vdiff turbulent
diffusivity to be consistent — validate against the existing k-ε/KPP column).

### P2 — Bedload transport

**Delivers.** A bedload mass flux `q_b` (kg/m/s) per non-cohesive class, its divergence feeding
the bed inventory. This is a **new conservative face flux** on both backends; it does NOT touch
the hydrodynamic Riemann solver — it is a separate advective flux of bed material, evaluated at
faces from the cell-centred `τ_b`, and differenced conservatively (structured: face-difference
like the barotropic `flux_h`; unstructured: cell-centric double-visit accumulation, the mandatory
no-neighbour-write pattern).

**Physics + citation.** Menu (dispatch by `bedload_formula`):
- **Meyer-Peter & Müller (1948)** — `Φ = 8(θ − θ_c)^{1.5}_+`, `q_b = Φ√((s−1)g d50³)`. The
  canonical non-cohesive bedload law and the primary test oracle.
- **van Rijn (2007)** — `q_b = 0.5 ρ_s d50 D_*^{-0.3} u_* (τ_b/τ_ce − 1)` (bed-load component of
  the TRANSPOR formulation).
- **Soulsby (1997)** threshold-of-motion Shields curve `θ_c(D_*)` for all classes.
- Bed-slope correction on magnitude and direction: **Bagnold / Koch & Flokstra** gravitational
  deflection (down-slope enhancement), needed so migrating bedforms don't spuriously stall — but
  keep it a gated add-on (`bedload_slope_corr`, default off) to preserve the clean MPM oracle.

**Hooks (four paths).** Bedload flux is a function of bed shear `τ_b` and the local flow
direction. On explicit and SI alike it is evaluated on the post-substep state and its divergence
added to the bed inventory at the same cadence as the suspended source terms — time-scheme
agnostic. Structured face flux reuses the `(i±½,j)` / `(i,j±½)` face-difference machinery;
unstructured reuses the edge-normal accumulation in the KNP flux pattern (compute per edge, add
to own cell only). The divergence lands in the **bed inventory**, not yet in `b` — P3 is what
moves `b`. So P2 is testable in isolation as a bed-inventory transport (a migrating hump of bed
mass) before Exner is wired.

**GPU.** Structured: per-face `do concurrent` producing `qb_x`, `qb_y` (persistent workspaces on
the multilayer state), then a per-cell divergence loop. Unstructured: per-cell double-visit
accumulation into a `qb_div(nclass,ncells)` workspace. All device-resident, mapped in
`*_enter_data`. The Shields/MPM closures are pure inline helpers.

**Knob.** `&sediment_nml enable_bedload`, `bedload_formula = "mpm" | "vanrijn07"`,
`bedload_slope_corr`. Default off ⇒ bit-identical.

**Test.** `test_bedload_mpm` — impose a uniform steady shear above threshold, assert the
computed `q_b` matches the closed-form **Meyer-Peter-Müller** rate to round-off (pure closure
test), and a 1D constant-flux case gives zero bed-inventory divergence in the interior (conserv-
ation). A second test: bedload of a bed-mass hump under uniform flow migrates at the analytical
kinematic celerity `c = dq_b/dz_b` (before Exner; inventory-space).

**Size:** Medium–Large. **Risk:** Medium (face-flux conservation on triangles; slope correction
can destabilise — keep gated).

### P3 — ★ EXNER BED EVOLUTION (the keystone) ★

**Delivers.** The seabed `b` finally moves. The bed inventory changes from P1 (deposition −
erosion) and P2 (−∇·q_b) are converted to a bed-elevation change and applied to
`state%barotropic%b` (structured) / `bed_cv` (unstructured), and the hydrodynamic state is
**reconciled** so the free surface stays physical and all bed-derived quantities are re-derived.

**Physics + citation.**
- **Exner (1925)** sediment continuity: `(1 − p) ∂z_b/∂t = −∇·q_b + (D − E)`, summed over
  classes, with bed porosity `p`. `q_b` from P2; `(D−E)` the net suspended exchange from P1
  (the deposition/erosion already moved mass into/out of the inventory — Exner converts the
  *inventory change* to a *thickness change* via `Δz_b = Δm_bed / (ρ_s (1−p))`, so bookkeeping
  is done in mass and converted once, guaranteeing the water-column↔bed mass balance).
- Optional WENO reconstruction of the bedload divergence for sharp bedforms: **Guérin et al.
  (2016)** — gated add-on, default the plain conservative difference.

**The hard part — sequencing and re-derivation (spelled out per path in §4).** The morphological
update is **operator-split after the hydrodynamic step** (an "asynchronous" morphodynamic
coupling; the standard and stable choice — the bed evolves on the slow morphological timescale,
much slower than the flow, so a split introduces negligible error and avoids a coupled
sediment-hydro Jacobian). Sequence each step:

1. Hydrodynamic sub-step advances `h`/`hu`/`hv` (+ vcoord ALE remap) as today, `b` frozen.
2. Suspended exchange (P1) + bedload divergence (P2) accumulate the per-class **bed-inventory
   change** `Δm_bed`.
3. Exner converts `Δm_bed → Δz_b` and updates `b_new = b_old + Σ_c Δz_b(c)`.
4. **Reconcile the water column.** `eta = h + b` (structured) / `stage = height + bed`
   (unstructured) is the physically continuous quantity; sediment fills the column from below, so
   hold the free surface fixed and set `h_new = eta_old − b_new` (⇒ `h` *decreases* where the bed
   aggrades). This is the mass-consistent choice: it neither creates nor destroys water. Guard
   `h_new ≥ 0` (a cell can dry by infilling — feed the wet/dry mask).
5. **Re-derive everything bed-dependent** (the crux — see §4 table): the vertical coordinate
   (sigma layer thicknesses `dsig·h`, or the ALE `z_ref` for ZSTAR_FULL which is built *from
   bathymetry*), the per-layer `h_layer` (must be re-partitioned to the new `H=h_new` and the
   suspended-tracer masses re-hung on the new thicknesses conservatively — reuse the ALE remap),
   the dry-cell mask (`h_new` vs `DRY_TOLERANCE`), and — for SI — the frozen face depths `Hf`
   that the Helmholtz operator caches.

**Morphological cadence.** Bed update runs every `morph_interval` hydro steps (default 1; >1 is a
cheap sub-cadence, distinct from morfac in P7). Cadence-gated exactly like `is_thermo_step()` for
ALE.

**Hooks (four paths).** See the dedicated §4 table — this is where the four paths genuinely
diverge and each must be stated. In brief: explicit paths re-derive the well-balanced
reconstruction inputs (nothing to re-solve, `b` just enters the next flux); SI paths must in
addition **rebuild the face-depth array** `si_build_face_depths(h_new,...)` (and, if the
multigrid coarse-operator caches depth-dependent coefficients, mark them stale) before the next
η-solve, or the Helmholtz operator is inconsistent with the moved bed.

**GPU.** Per-cell `do concurrent` for the Exner update and the `h = eta − b` reconciliation
(trivially parallel). The vcoord/layer re-partition reuses the existing ALE remap kernels (already
GPU-resident). The SI face-depth rebuild reuses `si_build_face_depths` (already a device kernel).
`b`/`bed_cv` are already device-resident state; the only new persistent arrays are the
per-class inventory deltas (small). No host round-trip — the bed update stays on device;
`b` copies to host only at output cadence.

**Knob.** `&sediment_nml enable_morph`, `morph_interval`. Default off ⇒ `b` static ⇒
bit-identical. (Note: `enable_morph=.false.` with `use_sediment=.true.` is exactly today's
behaviour — passive suspended + inert reservoir.)

**Test.** `test_exner_hump` — the **migrating-trench / migrating-hump** benchmark: a bed
perturbation under steady bedload-only flow migrates downstream at the analytical celerity
`c = 1/(1−p) · dq_b/dz_b` while conserving total bed mass; assert front position and total
(water + suspended + bed) mass to 1e-10. A second, cross-cutting test runs the **same** Exner
case through all four paths (structured/unstructured × explicit/SI) and asserts they agree within
scheme tolerance — this is the gate that proves the four-path reconciliation is consistent.
Standard harder oracle for a later validation example (not a unit test): a **scour-hole /
aggradation** case or the MorphoDyn/NCK migrating-trench experiment.

**Size:** Large (the keystone; touches every backend's depth plumbing). **Risk:** High — the
`h = eta − b` reconciliation interacts with wet/dry and the vcoord in ways that are easy to get
subtly non-conservative; the four-path agreement test is the safety net. This phase is where the
plan most needs a human review of the sequencing decision (see §6).

### P4 — Bed stratigraphy (active layer + sub-bottom layers)

**Delivers.** The bed inventory becomes vertically resolved: an **active (exchange) layer** that
trades with the water column, over a stack of **sub-bottom layers** each carrying per-class mass
fractions. Erosion peels from the active layer (and pulls up from below as it thins); deposition
adds to the active layer (and pushes down / spawns new layers as it thickens). This gives grain
sorting / armouring and a stratigraphic record.

**Physics + citation.** Active-layer bookkeeping: **Hirano (1971)** active-layer concept;
mixing-layer thickness from **van Rijn (2007)** (≈ half the bedform height, or a d50-scaled
floor). Layer merge/split bookkeeping is an algorithm, not a closure — mass-conservative by
construction. Note the unstructured SI state already carries a `bed_sub` sorted sub-bottom array
(`rki_si_unstr_state.F90:154`) — a hook to reuse, not a coincidence.

**Hooks (four paths).** Stratigraphy is a per-cell column of bed layers — identical on all four
paths (it does not touch the hydrodynamic solve at all; it only re-partitions *where in the bed*
the P1/P2/P3 mass exchange comes from). The active-layer top-fraction feeds P1/P2's per-class
`tau_ce` / availability. So the four paths share one kernel; the only per-path concern is that the
sub-bottom arrays are allocated on each backend's multilayer state and mapped in `*_enter_data`.

**GPU.** Per-cell `do concurrent` with an inner bed-layer loop bounded by a compile-time
`N_BED_MAX` (fixed-size stack arrays — NVHPC `local()` needs fixed size, cf. the
`NZ_STACK_MAX` rule). Merge/split writes only within the owning cell's column ⇒ no race. Bed-layer
arrays persistent on the state, mapped once.

**Knob.** `&sediment_nml enable_stratigraphy`, `n_bed_layer`, `active_layer_thickness` (or
`"vanrijn07"` predictor). Default: single active layer = P3 behaviour ⇒ bit-identical.

**Test.** `test_bed_stratigraphy` — a graded (two-class) bed under an erosion event coarsens the
active layer (winnowing / armouring): assert the surface fraction shifts toward the coarse class
and total per-class mass over the whole column is conserved to 1e-12. Plus a deposition case that
spawns and later re-erodes a sub-layer, returning to the initial stratigraphy (reversibility).

**Size:** Medium–Large. **Risk:** Medium (bookkeeping-heavy; conservation is the whole game).

### P5 — Cohesive sediment + consolidation

**Delivers.** Depth- and time-dependent cohesive-bed behaviour: a consolidation model that makes
the mud critical erosion stress `τ_ce` increase with depth in the bed and with residence time
(freshly-deposited mud is easily re-eroded; consolidated mud resists). Builds on the P4 sub-bottom
layers (each layer carries a consolidation state / dry density).

**Physics + citation.** Cohesive erosion **Ariathurai & Arulanandan (1978)** (already have) /
**Partheniades (1965)** excess-shear erosion; deposition **Krone (1962)**. Consolidation:
a **Gibson (1967)** finite-strain consolidation model is the rigorous choice but heavy; the
pragmatic v1 is a **profile/empirical consolidation** — layer-wise `τ_ce(z_bed, t)` and dry
density evolving toward equilibrium on a consolidation timescale (the standard engineering
approach; e.g. Winterwerp & van Kesteren 2004 for the framework). Flocculation-dependent settling
is **explicitly deferred** (SCHISM itself omits it) — keep `w_s` class-constant.

**Hooks (four paths).** Same as P4 — a per-cell bed-column process, four-path identical; it only
modifies the `τ_ce` that P1's cohesive erosion reads. No hydrodynamic-solver contact, so explicit
and SI are identical.

**GPU.** Per-cell `do concurrent`, inner bed-layer loop (reuses P4's fixed-size stack). The
consolidation ODE per layer is a cheap explicit or implicit-per-layer update — no global coupling.

**Knob.** `&sediment_nml enable_consolidation`, `tau_ce_consolidated`, `consolidation_rate`,
`bed_dry_density(:)`. Default off ⇒ constant `τ_ce` = P1 behaviour ⇒ bit-identical.

**Test.** `test_cohesive_consolidation` — deposit mud, hold quiescent for a consolidation time,
then ramp shear: assert the erosion threshold has risen (less mass eroded than an
un-consolidated bed) and mass is conserved. Standard benchmark for a validation example: an
**erosion-deposition tidal cycle** with an estuarine turbidity maximum signature.

**Size:** Medium. **Risk:** Medium (consolidation model choice is a modelling decision — see §6;
v1 empirical is defensible and cheap).

### P6 — Avalanching (bed-slope relaxation)

**Delivers.** Where the bed slope exceeds the **angle of repose**, relax it by moving bed material
downslope until the slope is at (or just below) repose — prevents unphysical vertical bed faces
that Exner alone can produce at steep fronts (scour edges, dune lee sides).

**Physics + citation.** Angle-of-repose avalanching — a slope-limited diffusive relaxation of the
bed; standard in morphodynamic models (e.g. the avalanching in **Roelvink et al. (2009)**
XBeach; underlying granular repose **Bagnold**). Distinct wet vs dry repose angles
(`phi_wet < phi_dry`).

**Hooks (four paths).** A per-cell bed-slope relaxation, iterated a bounded number of sweeps —
identical on all four paths (operates on `b`/`bed_cv` only, after the Exner update, before the
depth reconciliation is finalised — or as its own cadence-gated pass that itself triggers a
re-reconciliation of `h`). Structured uses `(i±1,j),(i,j±1)` neighbour slopes; unstructured uses
edge-neighbour slopes with the double-visit pattern (compute the down-slope flux per edge, apply
to own cell only). Must re-run the P3 `h = eta − b` reconciliation for any cell whose `b` it
moved.

**GPU.** Bounded iteration count (size the max sweeps from a single up-front grid reduction of the
worst slope, cf. the "size max_iter from one reduction" rule — do not fixed-budget over-iterate).
Each sweep is a `do concurrent`; the read-neighbour/write-own race is avoided by computing
down-slope transfers from a pre-sweep snapshot (the standard DC race fix in this codebase). No
neighbour writes on either backend.

**Knob.** `&sediment_nml enable_avalanching`, `repose_angle_dry`, `repose_angle_wet`,
`max_avalanche_sweeps`. Default off ⇒ bit-identical.

**Test.** `test_avalanching` — initialise a bed at super-repose slope, run to rest, assert the
final maximum slope ≤ repose angle everywhere and total bed mass is conserved to 1e-12
(avalanching only *moves* material). A cone-collapse case (a single tall bed spike relaxing to a
repose-angle cone) is a good visual/quantitative oracle.

**Size:** Small–Medium. **Risk:** Low–Medium (iteration convergence; conservation is easy if the
snapshot-flux pattern is used).

### P7 — Morphological acceleration factor (morfac)

**Delivers.** A **morfac** multiplier that scales the per-step bed change so that a short
hydrodynamic simulation represents a long morphological evolution (the flow is quasi-periodic —
e.g. a tidal cycle — and the slow bed can be advanced `morfac×` faster). Plus a spin-up ramp so
morfac engages gradually.

**Physics + citation.** **Roelvink (2006)** "online" morphological acceleration (multiply the bed
change per hydrodynamic time step by morfac; valid when the morphological timescale ≫ the
hydrodynamic timescale and the forcing is tidally periodic). Ramp per the same reference.

**Hooks (four paths).** A scalar multiplier on the Exner `Δz_b` (and on the bed-inventory
exchange, to keep the water-column↔bed mass balance exact under acceleration — the suspended
exchange must be scaled consistently or mass leaks). Path-agnostic: it multiplies the P3 update
identically on all four paths. The one subtlety, stated for the record: under morfac the bed
change per step is `morfac×` larger, so the `h = eta − b` reconciliation and wet/dry can move
faster — the `h_new ≥ 0` guard and the vcoord re-derive in P3 must hold at the accelerated rate
(they do, since they are exact per-step operations, but morfac is where they get stress-tested).

**GPU.** Trivial — a scalar in the existing Exner `do concurrent`. No new kernel.

**Knob.** `&sediment_nml morfac`, `morfac_ramp`. Default `morfac=1.0` ⇒ bit-identical.

**Test.** `test_morfac` — a periodic-forcing case run with `morfac=N` for time `T/N` must
reproduce (to the acceleration approximation's tolerance) the `morfac=1` case run for time `T`:
assert bed-change agreement and exact total-mass conservation at every step (mass conservation is
NOT approximate under morfac — only the time-representativeness is).

**Size:** Small. **Risk:** Low (physics trivial; the risk is entirely the mass-consistency of
scaling both bed and suspended exchange — the test guards it).

### P8 — Wave-current bed stress

**Delivers.** Upgrades the bed shear `τ_b` that feeds P1 (pickup), P2 (bedload) and P5 (cohesive
erosion) from current-only Manning to a **combined wave-current bottom stress**, which is what
actually mobilises sediment in the nearshore (wave orbital velocities dominate resuspension).

**Physics + citation.** **Grant & Madsen (1979)** wave-current boundary-layer interaction (the
nonlinear enhancement of the mean stress by wave orbital motion); **Soulsby (1997)** parametric
`τ_m`/`τ_max` fit (the cheap, GPU-friendly algebraic form — preferred over iterating the full GM
BBL each cell). Needs wave orbital velocity amplitude `u_orb`, period `T_w`, and direction from a
wave field.

**DEPENDENCY — FLAGGED.** This phase consumes a wave field (`u_orb`, `T_w`, wave direction).
The wave plan now exists — **`docs/plans/plan_wave_coupling.md`** — and P8 depends on its
**Phase 5**. P8 therefore has an **unmet cross-plan dependency**: it cannot land in full until
wave coupling (or at minimum a prescribed/offline wave-field ingest path) exists.
Two escape hatches so P8 is not fully blocked:
1. Land the **Soulsby stress algebra** against a *prescribed* wave field (a namelist/analytic
   `u_orb`, `T_w`) — immediately testable, immediately useful with offline wave output, and it is
   the exact seam the future wave plan will fill. This is the recommended P8 v1.
2. Full two-way coupling (wave field from an online wave model) is deferred to the wave plan.

**Hooks (four paths).** The combined stress is a per-cell closure feeding the existing sediment
source terms — path-agnostic, four-path identical (it changes `τ_b`, nothing structural). The
only new input is the wave-field arrays on the forcing state, mapped in `*_enter_data`.

**GPU.** Per-cell `do concurrent`; the Soulsby parametric form is a pure algebraic `*_impl`
helper. No iteration (that is the reason to prefer Soulsby's fit over iterating GM). Wave-field
arrays device-resident.

**Knob.** `&sediment_nml wave_current_stress`, `bbl_form = "soulsby" | "grant_madsen"`. Default
off ⇒ current-only Manning stress = P1 behaviour ⇒ bit-identical.

**Test.** `test_wave_current_stress` — impose a known `(u_current, u_orb, T_w, angle)` and assert
the combined `τ_max` matches the closed-form **Soulsby (1997)** parametric value; check the two
limits (`u_orb=0` recovers the pure-current Manning stress; `u_current=0` recovers the pure-wave
stress).

**Size:** Medium. **Risk:** Medium — but the *scheduling* risk (blocked on the wave plan) is
higher than the implementation risk. Keep P8 v1 to the prescribed-wave-field seam.

---

## 4. The four-path bed-update table

How the Exner update (P3) is sequenced relative to the hydrodynamic step, and what is re-derived
when `b` moves, **for each of the four paths**. Common to all: the morphological update is
operator-split *after* the hydrodynamic sub-step (asynchronous coupling); mass bookkeeping is done
in the bed inventory (mass) and converted to `Δz_b` once; the free surface `eta`/`stage` is held
fixed and the water depth is back-derived (`h = eta − b`), so water mass is untouched.

| | **Structured (Cartesian A-grid)** | **Unstructured (triangles)** |
|---|---|---|
| **Explicit (HLL/HLLC, KNP)** | **Sequence:** RK2 hydro step (`b` frozen, well-balanced `eta=h+b` reconstruction) → P1 exchange + P2 bedload divergence into per-class inventory → Exner `b += Σ Δz_b` → set `h = eta − b`, guard `h≥0`. **Re-derive:** (1) nothing to *re-solve* — `b` simply enters the next step's `eta=h+b` minmod reconstruction (`rki_kernel_flux.F90:148`); (2) per-layer `h_layer` re-partition to new `H` via sigma `dsig·h` or ALE remap; suspended-tracer masses re-hung conservatively by the existing ALE remap; (3) ZSTAR_FULL `z_ref(0:nz)` is rebuilt *from bathymetry* → must be recomputed from `b_new`; (4) dry mask `h` vs `DRY_TOLERANCE`. | **Sequence:** same operator-split. RK2 KNP hydro step (`b`=`bed_cv` frozen, Audusse stage/`bed_ev` reconstruction) → P1/P2 into inventory (cell-centric double-visit, no neighbour writes) → Exner `bed_cv += Σ Δz_b` → `height_cv = max(stage − bed_cv, 0)`. **Re-derive:** (1) `bed_ev` edge reconstruction (`bed_ev = stage_ev − height_ev`) picks up new `bed_cv` next step; (2) per-layer thickness re-partition via unstructured ALE remap (`ml_remap_tracer_unstr`); (3) ZSTAR_FULL `z_ref` (`ml_build_vcoord_zstar_full_unstr`) rebuilt from `bed_cv`; (4) dry mask `height_cv` vs min-depth. |
| **Semi-implicit (Casulli)** | **Sequence:** SI η-solve (Casulli CG-Helmholtz) advances `eta`/`h` (`b` frozen; `si_build_face_depths(h,...)` gives the Helmholtz face depths `Hf`) → P1/P2 into inventory → Exner `b += Σ Δz_b` → `h = eta − b`, guard `h≥0`. **Re-derive: everything the explicit path does, PLUS the SI operator inputs** — the frozen face depths `Hf_x/Hf_y` (`rki_si_helmholtz.F90:61`) are a function of `h`, so after `h` changes they **must be rebuilt by `si_build_face_depths(h_new,...)`** before the next η-solve, or the Helmholtz operator is inconsistent with the moved bed; if the multigrid (`rki_si_mg`) caches depth-dependent coarse-grid coefficients, mark them stale. Wet/dry in SI uses the frozen-face-depth branch (`rki_si_helmholtz.F90:316`) — the morph-driven dry/rewet must go through the same freezing logic. | **Sequence:** same as structured-SI but on triangles: `rki_si_unstr_step{,_ml}` η-solve (`bed_cv` frozen, TPFA Helmholtz face depths) → P1/P2 into inventory → Exner `bed_cv += Σ Δz_b` → `height_cv = max(stage − bed_cv,0)`. **Re-derive:** explicit-unstr set PLUS the unstructured Helmholtz face depths (`rki_si_unstr_helmholtz`) rebuilt from `height_cv`; the SI-unstr wet/dry (`rki_si_unstr_wetdry.F90`) mask refreshed; MG coarse coefficients (`rki_si_unstr_mg`) marked stale if depth-cached. The `bed_sub` sorted sub-bottom array already on `rki_si_unstr_state` (:154) is the natural home for P4 stratigraphy on this path. |

**The single load-bearing difference:** explicit paths need only make `b_new` available to the
next flux reconstruction (the well-balanced scheme does the rest); **SI paths must additionally
rebuild the cached Helmholtz face-depth arrays** (and invalidate any depth-dependent MG operator
cache) because the implicit η-solve froze a depth that the bed just changed. Getting this wrong is
silent: the solve still converges, but to the wrong (stale-depth) operator. The four-path
agreement test in P3 is what catches it.

---

## 5. Conservation & positivity strategy

**Total sediment mass = Σ_columns [ Σ_k hTr_sed(k) (suspended) ] + Σ_classes bed_inventory + bed
material represented as `b`.** The invariant is enforced by doing *all* exchange in **mass units**
through a single reservoir and converting to elevation exactly once:

- **Suspended ↔ bed inventory** (P1): erosion moves `dm` from inventory to the bottom suspended
  layer; deposition/settling moves the column deficit into the inventory (the existing
  `bed(i,j) += mass_old − mass_new` pattern, `rki_ml_tracers.F90:1048`). Both are exact transfers
  — no source/sink. Retain the positivity cap `dm = min(dm, bed_available)`.
- **Bedload divergence** (P2): a conservative face flux — the mass leaving a cell across a face
  equals the mass entering the neighbour (structured face-difference is antisymmetric by
  construction; unstructured interior edges are computed twice but equal, the codebase's
  standard double-visit conservation). Net inventory change = −∇·q_b, globally telescoping to
  boundary fluxes only.
- **Inventory → bed elevation** (P3): `Δz_b = Δm_bed / (ρ_s (1−p))`, a pure unit conversion —
  the mass that left the inventory becomes bed volume, nothing created/destroyed.
- **Water mass** under a moving bed (P3): holding `eta` fixed and setting `h = eta − b`
  conserves water exactly (water surface unchanged; only the bed floor moved). This is *why* the
  reconciliation holds `eta` rather than `h` — the alternative (hold `h`, let `eta` rise) would
  inject spurious potential energy / not represent infilling.
- **Stratigraphy** (P4) and **avalanching** (P6) only *move* bed material between layers / between
  cells — merge/split and down-slope transfer are exact mass redistributions, asserted to 1e-12.
- **Morfac** (P7): scales bed AND suspended-exchange consistently so the water-column↔bed balance
  stays exact under acceleration (the test asserts per-step mass conservation is *exact*, not
  approximate).
- **Global diagnostic.** Add a coastal sediment budget line (total suspended + total inventory +
  Σ ρ_s(1−p)·b·area) to the console stats, mirroring `compute_total_mass` — the report flagged the
  *absence* of a continuous tracer/sediment budget diagnostic as a credibility gap; close it here.

**Positivity / boundedness under wetting-drying** (the known ~1–2 %/tidal-cycle intertidal tracer
leak, gaps #3): suspended sediment is a registry tracer, so it inherits the same dry→`hTr=0`
destruction / rewet→clamp behaviour as salt. Sediment-specific mitigations:
- **Sediment concentration has a natural floor of 0, not a physical min like salinity** — so the
  clamp's *injection* branch (which invents salt to reach `S_min`) must be **disabled for sediment
  classes** (`tr_min = 0`), turning the rewet behaviour from "inject mass" into "start clean". A
  dried cell's suspended load should **settle to the bed inventory *before* the cell dries**, not
  vanish — add a "drain-on-drying" step (when a cell transitions wet→dry, deposit its remaining
  suspended mass into the bed inventory rather than zeroing it). This makes sediment *more*
  conservative across wet/dry than salt, because the bed is a legitimate sink that salt lacks.
- The **Exner-driven drying** (a cell infills until `h < DRY_TOLERANCE`) routes the residual
  suspended mass to the inventory by the same drain-on-drying path, so morphological drying is
  mass-exact even where hydrodynamic drying is not.
- Bed inventory and `b` are never subject to the water clamp, so the *morphological* mass is
  exactly conserved regardless of the wet/dry tracer leak — the leak is bounded to the suspended
  fraction in intertidal cells and further reduced by drain-on-drying. Document the residual
  honestly (as the report does for salt).

---

## 6. Risks & open questions needing a human decision

1. **P3 reconciliation choice (hold `eta` vs hold `h`).** The plan holds the free surface fixed
   and back-derives `h = eta − b` (mass-consistent infilling). This is the physically standard
   choice, but it means aggradation *shoals* the water column and can dry cells — confirm this is
   the intended morphodynamic behaviour (vs. e.g. a fixed-`h` convention used in some idealised
   flume setups). **Decision needed before P3.**
2. **Asynchronous (operator-split) vs synchronous morphodynamic coupling.** The plan splits bed
   evolution after the hydro step (standard, stable, cheap). A fully-coupled sediment-hydro
   Jacobian would be more accurate at large morfac but is a different (and much larger) numerical
   project. Confirm split is acceptable (it is the SCHISM/Delft3D/XBeach norm).
3. **SI Helmholtz-operator staleness.** Whether the SI multigrid caches depth-dependent
   coarse-grid operators that must be invalidated on every bed change — needs a code check during
   P3; if it recomputes face depths every solve anyway, the re-derive is a no-op and cheaper. Flag
   for the P3 implementer.
4. **Consolidation model fidelity (P5).** v1 is empirical (layer-wise `τ_ce(z,t)` relaxation).
   A rigorous Gibson (1967) finite-strain model is much heavier. Confirm empirical is acceptable
   for the target applications (it is standard for engineering estuary work).
5. **P8 blocked on the non-existent wave plan.** `plan_wave_coupling.md` is absent. Decide whether
   to (a) land the prescribed-wave-field seam now (recommended) or (b) defer P8 entirely until the
   wave plan exists. Wave-current stress is the single largest driver of nearshore sediment, so
   the seam is high value even before online waves.
6. **Grain-class count vs performance.** The inner class loop is serial within each cell's `do
   concurrent`. For `nclass` beyond ~5 (rare in coastal practice) the register pressure and stack
   arrays (`local()` fixed-size) may hurt — set a compile-time `N_SED_CLASS_MAX` and document the
   ceiling. Decide the ceiling.
7. **Bed-slope bedload correction (P2) stability.** The gravitational down-slope term can
   destabilise on steep beds; kept gated (`bedload_slope_corr` default off). Confirm it may ship
   off-by-default with avalanching (P6) as the primary steep-slope safeguard.
8. **Interaction with the untallied structured negative-depth clamp** (report §4): the structured
   `h<0→0` clamp adds water untracked. Under morphology-driven drying this could double-count with
   drain-on-drying. The P3 implementer must route morphological drying through the sediment drain
   *before* the hydro clamp sees the cell. Flag.

---

## 7. Out of scope

- **Sea ice** — out of scope for this repo by decision.
- **Spectral wave model** — P8 consumes a wave field but does not build one; wave generation is
  the wave plan's job (and the community norm is to *couple*, not reimplement).
- **Flocculation-dependent settling** (variable `w_s` from floc dynamics) — deferred; `w_s` stays
  class-constant (SCHISM omits it too).
- **Bedform / roughness predictors** (ripple/dune height → hydraulic roughness feedback) — the v1
  bed stress uses Manning (P1–P7) or wave-current (P8); dynamic roughness from bedforms is a later
  refinement, not in this plan.
- **Wave-asymmetry / acceleration-skewness bedload** (Elfrink, Hoefel-Elgar, Dubarbier) — deferred
  to a wave-coupled follow-on.
- **Sediment-density buoyancy feedback** on the EOS (turbidity-current-like effects) — the EOS
  does not see suspended sediment in v1; flag as a possible follow-on if turbidity currents become
  a target.
- **WENO Exner scheme as default** — Guérin et al. (2016) WENO reconstruction is offered as a
  gated add-on in P3, not the default (the plain conservative difference is the default oracle).
- **Water-quality / biogeochemistry** — a separate capability (schism_parity_report gap #3), not
  morphodynamics.
```
