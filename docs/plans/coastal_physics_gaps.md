# Rakali coastal path — physics-credibility assessment

**Question answered:** for the applications this model targets — coastal, estuarine, tidal, storm
surge, river plume, sediment/morphodynamics — what physics is *missing or thin enough* that a
domain scientist would reject it for real work, ranked by scientific consequence.

**Scope / method.** Read-only. This repo is the coastal path and nothing else — `sim_type`
accepts only `'coastal'` (`src/core/rki_config.F90:847`); there is no ocean dyn-core and no ice.
Every claim below names the file/procedure checked in-source; the docs had drifted (the "SI is
structured-2D-barotropic-only" claim is stale — see §6) so **source is authority here**. Companion
report `schism_parity_report.md` is a feature inventory; this is a credibility read of the physics
that *is* present.

---

## 1. Executive summary — could a domain scientist use this today, and for what?

**Yes, for a well-defined slice; no, for several flagship coastal problems as-is.**

Rakali's coastal core is a genuinely strong, well-balanced, GPU-native shallow-water / hydrostatic
multilayer engine: Audusse hydrostatic reconstruction (verified `rki_kernel_flux.F90`), flux-form
finite-volume mass conservation to round-off with *consistent* barotropic/tracer mass fluxes
(`mass_flux_x/y` shared between the h-update and the CWC tracer transport), a column-conservative
tracer clamp, SM03-Jacobian baroclinic PGF with reference-density-anomaly subtraction, PP81 + KPP +
k-ε closures on **both** backends, and an opt-in Casulli semi-implicit free-surface mode that (unlike
the stale doc) runs on **both** backends at **both** 2D barotropic and 3D sigma multilayer.

**Usable today, with care:**
- **Barotropic tidal hydrodynamics and storm surge (wave-free)** — the strongest fit. Wind stress
  (Large & Pond), atmospheric-pressure/inverse-barometer forcing, multi-constituent tidal BCs,
  robust wet/dry, and exact lake-at-rest are all present.
- **Subtidal / partially-mixed estuarine circulation and salt transport on sigma or z\*-lite** —
  credible, provided the domain does not depend on large intertidal flats (see gap #3).
- **River plumes (qualitative)** — buoyancy, KPP, and 1%-conservative salt transport are there.
- **Suspended-sediment transport as a passive tracer** — conservative settling + Ariathurai erosion.

**Not credible as-is (single blocking item each):**
- **Sediment *morphodynamics*** — the seabed never moves (gap #2). No Exner, no bedload.
- **Surf-zone / wave-affected nearshore, and the wave contribution to storm-surge total water
  level** — no wave coupling of any kind (gap #1).
- **Tidal-energy *resource extraction*** — no turbine/energy-sink physics (gap #10).
- **Intertidal-mudflat salt-intrusion studies** — ~1–2 % salt/tidal-cycle leak (gap #3).

The core numerics are sound; the deficits are coupled Earth-system physics (waves, morphodynamics)
and a handful of *thin* spots in existing physics (intertidal tracer conservation, the baroclinic
PGF's wet/dry behaviour, air-sea flux breadth). GPU-first explicit design is a deliberate choice,
not a defect (§6).

---

## 2. Ranked gap table (scientific consequence first)

`M` = missing, `T` = thin/present-but-weak. Size is rough implementation effort.

| # | Gap | Consequence for target apps | Evidence (checked) | M/T | Size | GPU fit | Primary citation |
|---|-----|------------------------------|--------------------|-----|------|---------|------------------|
| 1 | **No wave coupling** — radiation stress, Stokes drift, vortex force, wave setup, wave-enhanced bottom stress | Blocks surf-zone, nearshore, and the wave part of storm-surge total water level (wave setup 10–30 % of surge); blocks wave-driven sediment | repo-wide grep of `src/` for radiation-stress/stokes/vortex/wwm/ww3/swan/wbbl/grant-madsen → **zero hits** | M | Large (couple, don't port the spectral model) | Coupling terms: excellent DC fit | Longuet-Higgins & Stewart (1964); Bennis et al. (2011); Grant & Madsen (1979) |
| 2 | **No morphodynamics** — bed elevation never evolves; single-class suspended only, no bedload | Sediment *morphodynamics* impossible; only passive suspended transport | `ml_erode_bed_sediment`/`ml_settling_tracer` (`rki_ml_tracers.F90:963,1054`) move mass between a suspended tracer and a 2D `bed` reservoir **distinct from bathymetry `b`**; `b` never updated | M | Large | Exner update per-cell = good DC fit; bedload flux reuses HLL machinery | Exner (1925); Meyer-Peter & Müller (1948); van Rijn (2007) |
| 3 | **Salt/tracer non-conservation under intertidal wet/dry cycling** (~1–2 %/tidal cycle) | Blocks intertidal-mudflat salt-intrusion; error accumulates over spring–neap | dry→`hTr=0` (`ml_distribute_tracer:887`), rewet→clamp injects to `S_min` (`ml_clamp_tracer:1219-1248`); documented `CAPABILITIES...:200-201` | T | Medium | N/A (algorithmic) | — (numerical; edge-reconstruction / tracer-salvage fixes known) |
| 4 | **Baroclinic PGF: hard reduced-gravity cap + dry-neighbour skip** | Suppresses legitimate baroclinic accel over steep shelf-break; **zeros baroclinic forcing entirely adjacent to any dry column** → degraded salt-intrusion front dynamics at estuary mouths / intertidal edges | `ml_baroclinic_pressure_gradient_jacobian` caps `|dp/dx|` to `sqrt(g'H)/dt` (`rki_ml_pressure.F90:501-506`); `all_neighbours_wet` → `cycle` (`:333-334, 520-539`) | T | Medium | fine (already DC) | Shchepetkin & McWilliams (2003); Mellor et al. (1994); Haney (1991) |
| 5 | **Air-sea flux breadth** — no bulk (COARE/Zeng) heat/moisture flux; no NetCDF-driven heat/E–P; no time-varying river Q(t); wind drag capped 2.5e-3 (no high-wind rolloff) | Multi-day realistic thermal-estuary and tropical-cyclone surge need these; user must precompute all fluxes offline | `drag_coefficient` = Large & Pond, `min(cd,2.5e-3)` (`rki_forcing.F90:562-576`); heat_flux/E–P programmatic-only (`CAPABILITIES...:244-250`); no `river Q` time-series knob (grep) | T | Medium | per-cell, excellent DC fit | Fairall et al. (2003); Large & Pond (1981); Powell et al. (2003) |
| 6 | **Turbulence: default PP81 is explicit + Ri≥0-clamped (no convective adjustment); no GLS/k-ω/MY2.5** | Under-mixes convective (nocturnal-cooling / dense-water) events unless user opts into k-ε/KPP; no generic-length-scale family | PP81 clamps `Ri<0→0` (`rki_ml_vertical_viscosity.F90:98`, `..._tracer:1685`); flux-capped `nu·dt<0.5h²` (`:104`); GLS "not yet built" (`CLOSURE_MATRIX.md:46`) | T | Medium | k-ε/KPP already DC on both backends | Pacanowski & Philander (1981); Large et al. (1994); Umlauf & Burchard (2003) |
| 7 | **Bottom stress is current-only Manning; no wave-current BBL** | Under-predicts bottom stress / sediment resuspension in wave-exposed coast (ties to #1) | `ml_apply_bottom_friction` = quadratic Manning on lowest active layer (`rki_ml_friction.F90:25-100`); no wave orbital velocity input | M/T | Medium | per-cell DC fit | Grant & Madsen (1979); Soulsby (1997) |
| 8 | **EOS evaluated at p=0** (potential density); linear default | Low consequence for coastal depths; matters only for deep-shelf thermobaric/cabbeling | `ml_equation_of_state_shared` uses `eos%p_ref` default 0 Pa (`rki_ml_eos.F90:56-80`); linear kernel `:82-134` | T | Small | excellent | Wright (1997); Jackett & McDougall (1995); IOC/TEOS-10 (2010) |
| 9 | **f-plane only** — Coriolis is a single scalar `coriolis_f`; no beta, no spatial f | Fine regionally; mis-rotates on large shelf domains spanning latitude | `state%barotropic%coriolis_f` scalar (`rki_barotropic_state.F90:70`); `ml_apply_coriolis(...,f)` scalar (`rki_ml_friction.F90:102`) | T | Small | trivial | standard geophysical FD |
| 10 | **Tidal forcing thin: per-edge constant harmonics; no 2D tidal-harmonic file, no tidal potential/SAL; no energy-extraction sink** | Limits large-domain tidal accuracy; blocks tidal-*energy* (turbine array) work | per-edge `tc_amp/phase/omega` (`rki_boundary.F90:68-70`); no tidal-potential/SAL (grep); no turbine/actuator/momentum-sink (grep → zero) | M/T | Small–Med | per-cell adds are good DC fit | Cartwright & Tayler (1971); Ray (1998); actuator-disc hydraulics |
| 11 | **Non-hydrostatic path is homogeneous-density, structured-only** | No stratified-NH phenomena (internal solitary waves, NH plume fronts) | `nonhydrostatic_state_t` carries only `u_layer`/`q_nh`/`u_layer0`, no density/S/T (`rki_nonhydrostatic_state.F90:30-40`) | M | Large | — | standard NH shallow-water |

---

## 3. Per-application readiness

### Storm surge — **usable today for wave-free regional surge; the strongest fit.**
Present and verified: Audusse well-balanced wet/dry (`rki_kernel_flux.F90`), KP07 thin-cell
velocity desingularization + Froude backstop (`rki_kernel_update.F90`; replaced the 2026-08
hard-zero/`FROUDE_CAP=10` pair that froze wet/dry fronts), multi-constituent tidal BCs
(up to 10, `rki_boundary.F90`), wind stress (Large & Pond) and atmospheric-pressure forcing, and
round-off mass conservation. **Single blocking item:** no **wave setup / radiation stress** (gap #1)
— in energetic events wave setup is a first-order fraction of peak total water level, so surge
without it is biased low nearshore. Secondary: wind-drag cap 2.5e-3 has no high-wind saturation
rolloff (Powell 2003) → over-drag in category-4/5 winds; f-plane (#9) for basin-spanning shelves.
Verdict: **credible for tide+surge hydrodynamics; not for total-water-level where waves matter.**

### Estuarine salt intrusion — **credible for subtidal/partially-mixed on sigma or z\*-lite; blocked for intertidal mudflats.**
The machinery is real: baroclinic SM03-Jacobian PGF with rho-anomaly (`rki_ml_pressure.F90`),
CWC-conservative tracer advection (PPM), KPP/k-ε available, 1%-level salt conservation in the
plume test. **Two thin spots gate the hard cases:** (a) gap #3 — ~1–2 % salt leak per tidal cycle
over intertidal flats from the dry→zero/rewet→clamp-inject cycle (`CAPABILITIES...:200-201`;
the recommended vcoord for these is sigma/z\*-lite, not ZSTAR_FULL); (b) gap #4 — the baroclinic
PGF is **zeroed adjacent to any dry column** (`all_neighbours_wet→cycle`) and hard-capped to the
local internal-wave speed, both of which blunt the salt-wedge front dynamics precisely where the
intrusion is sharpest. Verdict: **credible for deep-channel/subtidal estuaries; a salt-budget study
over a drying-flat estuary would show drift the author already documents.**

### River plume — **credible qualitatively; diffusive on triangles.**
Buoyant transport, KPP surface mixing, and conservative salt (plume test asserts <1 %) are present
and the surface=`k=nz` convention is regression-gated. Limits: on the **unstructured** backend the
"PPM" tracer knob actually selects **MUSCL piecewise-linear (1st-order)** (`CLOSURE_MATRIX.md`
honesty note) → a numerically diffuse plume front; and no wave-current mixing at the plume face
(#1). Verdict: **fine for plume structure/spreading studies; not for sharp-front fidelity on
triangular meshes.**

### Tidal energy — **currents yes; extraction no.**
Barotropic tidal currents are well-resolved (explicit, wave-faithful). **Blocking item:** there is
**no energy-extraction physics** — grep for turbine/actuator/energy-sink/momentum-sink returns
nothing. Also tidal forcing is per-edge constant harmonics with no spatially-varying 2D tidal
field or tidal potential/SAL (#10), limiting large-domain resource accuracy. Verdict: **usable for
tidal-current *characterisation*; not for turbine-array resource/interaction assessment.**

### Sediment / morphodynamics — **passive suspended transport only; NOT a morphodynamic model.**
Verified: `ml_settling_tracer` (conservative PPM settling into a `bed` reservoir) +
`ml_erode_bed_sediment` (Ariathurai-Arulanandan 1978, Manning-consistent bed stress, positivity-
capped). **Blocking item:** the seabed elevation `b` **never changes** — `bed` is a mass reservoir
disjoint from bathymetry, so there is no Exner update and no morphodynamic feedback (gap #2).
Additionally: single grain class, **no bedload**, no bed stratigraphy, no cohesive flocculation,
no wave-current bottom stress (#7), and nothing on the unstructured backend. Verdict: **can carry a
suspended-sediment tracer; cannot do bed evolution / morphodynamics.**

---

## 4. Conservation & boundedness audit (mechanism named)

| Property | Status | Mechanism (file:proc) | Caveat |
|---|---|---|---|
| **Water mass (barotropic)** | Conserved to round-off | Flux-form FV; symmetric per-face HLL fluxes, interior faces computed twice but equal (`rki_kernel_flux.F90:358-376`); global start/end tally `compute_total_mass` (`rki_driver.F90:712`) | Negative-depth clamp `h<0→0` (`rki_kernel_update.F90:45`) **adds** water on structured and is **not tallied**; unstructured tallies it in `state%mass_error` (`protect_unstr`) |
| **Water mass (multilayer)** | Conserved to round-off | Barotropic h distributed by `dsig` (`ml_distribute_barotropic`); total change = `0.5·dt·(flux_h¹+flux_h²)` (`rki_solver.F90:1602`); boundary/precip/sponge sources tracked | same clamp caveat |
| **Momentum** | Not a conservation target; bounded | KP07 desingularization + Froude backstop + `ml_truncate_velocity` (0.9·cfl·dx/dt) are **momentum rescales** (mass untouched) | intentional stabilisers |
| **Tracer (salt/heat) mass** | Conserved in-column; **not** across wet/dry | CWC flux-form advection shares barotropic `mass_flux` (uniform-tracer invariant, test 1e-12); `ml_clamp_tracer` two-pass column-conservative (`rki_ml_tracers.F90:1130`) — redistributes when column-mean in-bounds, **injects/destroys only when column-mean out of range** | dry→`hTr=0` destroys; rewet clamp injects → gap #3 |
| **Tracer boundedness** | Enforced | PPM Colella-Woodward monotonicity limiter (`:485-527`); `ml_clamp_tracer` hard bounds; RK2 `hTr≥0` | hard bound's residual is the non-conservative site |
| **Salt/heat global monitoring** | **Absent on coastal** | only `total_mass` printed; no coastal salt/heat console budget (grep) | user cannot easily see tracer leak |
| **Sediment (suspended+bed)** | Conserved to round-off | `bed += (mass_old−mass_new)`; erosion capped at `min(dm,bed)` (`rki_ml_tracers.F90:1048,1122`) | `bed` ≠ bathymetry (no morpho) |

Net: **mass and per-column tracer conservation are genuinely strong and well-tested** (see §5); the
two real holes are (a) the untallied structured negative-depth clamp and (b) the wet/dry tracer
leak. There is **no continuous salt/heat conservation diagnostic** on the coastal path — a
credibility gap for a user who needs to trust a multi-week salt budget.

---

## 5. Validation coverage & notable absent benchmarks

**Genuinely analytical (exact reference), verified in `tests/`:**
- **Ritter dam break** structured (`test_dam_break.F90`, exact Riemann via `ritter_dambreak.F90`,
  L2 < 2.5 %, mass < 1e-10) and unstructured (`test_dam_break_unstr.F90`, L1 < 3.5 % **plus the only
  mesh-refinement convergence-rate test in the coastal set**, rate ≥ 0.85).
- **Thacker parabolic bowl** (`test_thacker_bowl.F90`, moving shoreline; eta L2 < 8 %, mass < 1e-4,
  no-negative-depth) — but the exact solution is a *simplified* planar form and 8 % is loose.
- **WENO tracer reconstruction** (`test_tracer_recon.F90`, 39 tests) — polynomial-exactness to
  5e-12, conservation 1e-12, mirror symmetry 1e-14: the most rigorous coastal-tracer suite.
- **k-ε** (`test_keps.F90`) — homogeneous-decay power law (2 %) and wall log-law (1e-6): strongest
  closure test.
- **Semi-implicit** seiche/channel/wind-setup analytical checks (`test_semi_implicit_*`).
- Exact-arithmetic kernel suites: `test_friction.F90` (semi-implicit Manning/Coriolis KE-conserving),
  `test_multilayer.F90` (47 tests), `test_pressure.F90` (NH PCG).

**Conservation/boundedness well covered:** mass (barotropic/multilayer/split/NH), salt to 1e-10…1e-12
(`test_zstar_full_tracer_conservation.F90`, `test_2d_tracer.F90`, `test_multilayer_unstr.F90`),
sediment (`test_sediment.F90`), tracer race regression (`test_ml_tracer_advect_race.F90`).

**Weak-tolerance "analytical" (really property tests):** **lock exchange** vs Benjamin front speed
at **±50 %** (`test_lock_exchange.F90`); **Kato-Phillips** entrainment as a **1.6–2.4 ratio band**,
with the clean √t law validated only in an *off-repo Python* reference, not the Fortran solver
(`test_kato_phillips.F90`). These are the two flagship stratified-physics benchmarks and both pass
at order-of-magnitude bands.

**Standard coastal benchmarks with NO coverage found (verified absent):**
- **Carrier-Greenspan** analytical runup (periodic beach oscillation) — Thacker is the *only*
  moving-shoreline analytical case.
- **Dry-bed dam break** (Ritter dry front / Stoker) — both dam-break tests are **wet-bed**; the
  dry leading edge is only smoke-tested (no-NaN), never against analytic front position/celerity.
- ~~**Well-balanced C-property over varying bathymetry** as a standalone structured assertion.~~
  **CORRECTED 2026-08-07 — this claim was wrong.** `tests/test_lake_at_rest.F90` already asserts
  the C-property standalone on the structured backend. It was *weak*, not absent: a smooth,
  fully-submerged parabolic bump, 100 steps, 1e-10 absolute, structured 2D only. The real gap was
  strength and coverage, now closed by `tests/test_c_property.F90` (irregular bed, emergent
  island, 1000 steps, round-off tolerances, both backends + 3D multilayer).
- **Steady transcritical flow over a bump** (classic SWE benchmark with shock).
- **Analytical tidal-flat wetting/drying** (periodic-forcing runup) — only `no_nan_tidal` smoke.
- **Internal-seiche / baroclinic-wave period** vs analytic frequency (seiches run only as
  conservation drivers, never checked against a period).
- **Structured HLL convergence order** — no refinement test at all (only the unstructured KNP one).
- Storm-surge observational/standard cases (e.g. GRS2018), formal Sydney Harbour comparison
  (`CAPABILITIES...:638`).

Demonstration configs (not analytical validation) exist under `validation_examples/coastal/`:
`castro_dam_break`, `river_plume`, `tidal_estuary`, `seamount_sediment`, `multilayer_benchmark`,
`semi_implicit`.

**Bottom line on validation:** per-kernel exact-arithmetic correctness and conservation are
*strong*; **analytical-solution accuracy at realistic tolerance and convergence-order verification
are thin** — a reviewer would want Carrier-Greenspan runup, a dry-bed dam-break front, a
well-balanced C-property assertion, and a tightened lock-exchange/Kato-Phillips before trusting the
stratified and wet/dry physics quantitatively.

---

## 6. "Not gaps" — deliberate design divergence

- **Explicit Godunov (gravity-wave-CFL) vs semi-implicit large-Δt.** A stated tradeoff
  (`CAPABILITIES...:48-70`): explicit HLL/KNP buys faithful tsunami/surge/bore/dam-break/shock
  physics and exact lake-at-rest that semi-implicit ELM schemes smear. Not a defect.
- **Semi-implicit mode is broader than the stale doc says.** Verified in-source: `time_scheme=
  "semi_implicit"` (Casulli CG-Helmholtz) exists on **both** backends at **both** 2D barotropic and
  3D sigma multilayer — structured (`kernels/structured/semi_implicit/rki_si_step{,_ml,_wetdry,
  _closures}.F90`) and unstructured (`kernels/unstructured/.../rki_si_unstr_step{,_ml}.F90`,
  `rki_si_unstr_wetdry.F90`, `rki_si_unstr_closures.F90`). Any report repeating "structured 2D
  barotropic only" is wrong.
- **Triangle-only unstructured** (KNP central-upwind) — a scheme choice, not a hole for most meshes.
- **Modern vcoord stack** (sigma / z-sigma / z\*-lite / z\*-full / z\*-sigma via ALE remap, PCM..PQM)
  is a *different* lineage from SCHISM's LSC2/SZ — comparable generality, not a gap (isopycnal /
  true partial cells honestly absent).
- **No sea ice** — out of scope for this repo by decision. SCHISM has two ice models on its
  unstructured grid; Rakali will not. The one deliberate hole in the parity target.
- **PP81 explicit flux cap** — the `nu·dt<0.5h²` cap looks limiting but is rarely binding because
  the explicit Δt is already seconds (e.g. h=2 m, Δt=2 s ⇒ cap ≈ 1 m²/s), so the real turbulence
  limitation is PP81's Ri≥0 *physics* (#6), not the cap. Design-appropriate given the small Δt.

---

## 7. "Unverified" — what I did NOT check (silence ≠ absence)

- **Did not run anything.** All conservation/leak claims are from source structure, docstrings, and
  the shipped tests' *assertions* — not from executing the solver or the tests on a GPU build.
- **Unstructured backend** was read only where the structured/unstructured split was material
  (tracer clamp, protect_unstr mass tally, MUSCL order); I did not trace the unstructured
  solver pipeline (`rki_solver_unstr.F90`) end-to-end.
- **Semi-implicit path internals** — I confirmed the files exist and the doc's scope, but did not
  audit the Casulli wet/dry Newton conservation or the SI tracer overshoot claims in code.
- **NH path** — confirmed homogeneous-density from the state type; did not read the CG-Poisson
  solver or its convergence/robustness.
- **Boundary conditions** beyond tidal/wall — Chapman/Orlanski/sponge/nested robustness not audited;
  the nested BC's "smoke-level" status is from docs, not re-verified in code.
- **Actual magnitude of the baroclinic-PGF cap's effect** — I verified the cap and dry-skip exist
  but did not quantify how often they bind on a realistic shelf/estuary run.
- **Wind→stress high-wind behaviour** — read the drag-coefficient formula; did not check whether any
  gustiness/Charnock refinement exists elsewhere.
- **Forcing reader** (`rki_forcing.F90`) — confirmed it ingests u10/v10/MSLP/precip and converts
  wind→stress; did not verify heat-flux/E–P are truly absent from *every* file path (docs say
  programmatic-only).
- **Tracer registry passive-tracer path** beyond S/T — assumed from CLAUDE.md that arbitrary passive
  tracers ride the same conservative transport; not independently exercised.
