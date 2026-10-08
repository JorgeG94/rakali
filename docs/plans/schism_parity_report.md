# SCHISM vs Rakali (coastal path) — feature-parity gap analysis

**Author:** read-only exploration agent · **Date:** 2026-08-06
**Revised:** 2026-08-07 — re-verified line-by-line against both trees. Ocean/ice references
struck (this repo is coastal-only), semi-implicit maturity limits added, module sizes corrected.

**Scope.** What SCHISM can do that Rakali cannot. This repo *is* the coastal path and nothing
else: `sim_type` accepts only `'coastal'` (`src/core/rki_config.F90:847`), and there is no
`src/core/ocean/` or `src/core/ice/`. Any claim below of the form "Rakali has X on the ocean
path" was an artefact of the pre-port-over repo and has been removed.

**Target.** SCHISM **physics** parity on **both** time schemes — explicit Godunov *and* Casulli
semi-implicit. Numerical-method divergence (§4) stays deliberate; the physics gaps (§2) are the
roadmap. Sea ice is out of scope **by decision**, not by accident.

## Provenance / licensing boundary (read first)

SCHISM (`/home/jorge/nci/cdx/schism`) is **Apache-2.0** (verified: `LICENSE` is the
Apache 2.0 text; per-file headers e.g. `Ice/ice_fct.F90` carry
"Copyright 2014 College of William and Mary … Licensed under the Apache License, Version 2.0").
Apache-2.0 permits reading, but this report obeys the house rule **cite the paper, not the
codebase**. **No SCHISM source was copied.** Algorithms are described in prose; file paths and
identifiers appear only as navigation pointers; any future Rakali implementation must be
re-derived from the primary literature named below, not transcribed from SCHISM.
Note some sub-modules carry non-Apache lineage inside the SCHISM tree — GOTM is GPL, the
FABM host is dual Apache/GPL, EcoSim/ICM carry ROMS/FERI notices — a further reason to work
only from published equations.

Every "Rakali lacks X" claim below names the SCHISM module implementing X and states where I
looked in Rakali. Rakali coastal + shared physics were read under `src/core/{coastal,structured,
unstructured,multilayer}/`, `src/{ALE,equation_of_state,pressure_force,tracer,parameterizations}/`,
plus `docs/CAPABILITIES_AND_LIMITATIONS.md` and `docs/CLOSURE_MATRIX.md` (both authoritative and
code-drift-checked).

---

## 1. Executive summary (the honest headline)

SCHISM is a 20-year community estuarine/coastal Earth-system model: ~300 kLOC of physics beyond
the bare hydrodynamic core, spanning full 3D morphodynamic sediment transport, a coupled
third-generation spectral wave model, four interchangeable water-quality/biogeochemistry engines
behind a generic 12-slot tracer framework, two sea-ice packages, GOTM turbulence, hydraulic
structures, vegetation/marsh drag, and a rich open-boundary/nudging/tidal-potential apparatus.
Rakali's coastal path is a **fast, GPU-native, wave/shock-faithful hydrodynamic core** with S/T
baroclinic transport and a modern vertical-coordinate/ALE stack — but it has **essentially none of
SCHISM's coupled Earth-system modules**. The largest, most valuable gaps are **sediment/
morphodynamics**, **wave coupling**, and **water-quality/biogeochemistry**; each is a multi-month
capability, not a patch. The good news for Rakali's architecture: the *physics* in all three is
dominated by **per-cell / per-column source terms and reaction kinetics that are embarrassingly
parallel** and map cleanly onto `do concurrent` — the cost is breadth (many formulas, many state
variables, validation), not GPU-hostility. The one genuinely awkward-to-port capability is the
spectral wave model itself (its own global implicit spectral solver); the right move there is
*coupling* to an external wave model, not reimplementation. Numerically SCHISM is semi-implicit
Eulerian-Lagrangian on mixed tri/quad grids; Rakali is explicit Godunov (with an opt-in Casulli
semi-implicit mode) — this is a deliberate divergence, and several SCHISM "features" (large `Δt`,
implicit transport) are the flip side of physics Rakali resolves that SCHISM must filter.

---

## 2. Prioritised gap table (highest value first)

| # | Gap | SCHISM module (nav only) | Size / complexity | GPU-fit (`do concurrent`) | Primary literature (citation of record) |
|---|-----|--------------------------|-------------------|---------------------------|------------------------------------------|
| 1 | **3D + 2D sediment transport & morphodynamics** (multi-class, bedload+suspended, cohesive+non-cohesive, bed stratigraphy, Exner bed evolution, avalanching, roughness predictors, wave-current stress, morph accel factor) | `src/Sediment/` (3D, ~7.1 kLOC), `src/Sed2d/` (2D total-load, ~5.1 kLOC) | Very large (~12 kLOC). Rakali has only a single-class suspended stub | Good — per-cell/column source terms + per-cell Exner update are DC-friendly; bedload flux & avalanching iterations need care | Suspended cohesive: **Ariathurai & Arulanandan (1978)**; **Krone (1962)**. Bedload: **Meyer-Peter & Müller (1948)**, **van Rijn (2007)**, **Soulsby & Damgaard (2005)**, **Wu & Lin (2014)**. Total-load: **Engelund & Hansen (1967)**, **Ackers & White (1973)**, **Soulsby (1997)**, **Camenen & Larson (2011)**. Morphology: **Exner (1925)**; WENO Exner **Guérin et al. (2016)**. Settling/threshold: **Soulsby (1997)** |
| 2 | **Spectral wave model + wave-current coupling** (action-balance spectrum; radiation-stress AND vortex-force forcing; Stokes drift; wave setup; wave-enhanced bottom stress; wave breaking) | `src/WWMIII/` (WWM-III, ~80 kLOC); coupler `wwm_coupl_selfe.F90`; external WW3 field-coupling in `Hydro/misc_subs.F90` (ESMF/OASIS) | Huge for the wave model; small–moderate for the coupling terms | Spectral solver: poor DC fit (global implicit spectral advection) → **couple, don't port**. Coupling terms (rad stress, vortex force, Stokes): good DC fit | Model/source terms: **WAM — Hasselmann et al. (1988) / Komen et al. (1994)**; ST4 **Ardhuin et al. (2010)**; ST6 **Zieger et al. (2015)**; DIA **Hasselmann et al. (1985)**; breaking **Battjes & Janssen (1978)**. Coupling: radiation stress **Longuet-Higgins & Stewart (1964)**; vortex force **Bennis et al. (2011)**, **Uchiyama et al. (2010)** |
| 3 | **Water quality / biogeochemistry ecosystem** behind a generic multi-model tracer framework (12 slots) | `src/ICM/` (~3.7 kLOC), `src/COSINE/` (~2.1 kLOC), `src/EcoSim/` (~6.5 kLOC), `src/Fabm/` (coupler, ~2.0 kLOC), `src/FIB/` (~0.5 kLOC); framework in `Core/schism_glbl.F90` + `Hydro/schism_init.F90` | Large aggregate (~15 kLOC); but Rakali's tracer registry is a real head start | Excellent — reaction kinetics are per-cell ODE integrations, ideal for DC. Transport already exists in Rakali | ICM: **Cerco & Cole (1993)**; sediment diagenesis **DiToro (2001)**. CoSiNE: **Chai et al. (2002)**. EcoSim (bio-optical): **Bissett et al. (1999)**, light **Gregg & Carder (1990)**. FABM host: **Bruggeman & Bolding (2014)**. FIB: **Canteras et al. (1995)**, **Servais et al. (2007)** |
| 4 | **GLS turbulence family (k-ω, Mellor-Yamada 2.5) via GOTM** | GOTM coupling (`itur=4`, `do_turbulence`), vendored `GOTM3.2.5/` **and** `GOTM5.2/`; built-in `itur=3` in `Hydro/misc_subs.F90` | Moderate — incremental on Rakali's existing k-ε | Good — per-column implicit closures; Rakali already ships k-ε (implicit Thomas) both backends | **Umlauf & Burchard (2003)** (GLS); k-ω **Wilcox (1988)**; MY2.5 **Mellor & Yamada (1982)**, **Galperin et al. (1988)**; stability functions **Canuto et al. (2001)**, **Kantha & Clayson (1994)** |
| 5 | **Hydraulic structures** (weirs, orifices, culverts/pipes, radial gates, pumps-as-transfers, compound structures with gate operation) | `Core/hydraulic_structures.F90` (driven by `ihydraulics`/`nhtblocks` in `schism_step.F90`) | Small–moderate (~1 file) | Mixed — per-structure flux laws are trivial; the cross-domain flux transfer between non-adjacent node pairs is a scatter that needs care/atomics on GPU | Standard hydraulics: weir/orifice **Bos (1989)**; radial gate **USACE HEC-RAS** hydraulic reference; compound culvert-over-weir standard energy/momentum forms |
| 6 | **Vegetation & marsh drag** (bendable canopy momentum sink, wave dissipation by vegetation, dynamic marsh growth) | `schism_glbl.F90` (`veg_*`, `isav`, `age_marsh`), applied in `Hydro/misc_subs.F90` + `schism_step.F90`; wave-veg in `wwm_vegdiss.F90` | Small (drag) to moderate (dynamic marsh) | Excellent — a per-cell/per-layer momentum sink, trivially DC | Rigid canopy **Nepf (1999)**; flexible/bent canopy **Luhar & Nepf (2011)**; wave dissipation by vegetation **Mendez & Losada (2004)** |
| 7 | **Richer open boundaries: earth tidal potential + loading tide/SAL; file-driven 2D tidal harmonics; per-tracer OBC nudging** | `Hydro/schism_init.F90` (`ntip`/`tip_dp` tidal potential, `iloadtide` SAL, `iettype`/`ifltype`/`itrtype`, `inu_tr`, `*_nudge.gr3`) | Moderate | Good — potential forcing and nudging are per-cell adds; DC-friendly | Equilibrium tidal potential & Love numbers **Cartwright & Tayler (1971)**; SAL/loading **Ray (1998)**; radiation OBC **Flather (1976)** |
| 8 | **In-situ (pressure-dependent) nonlinear EOS** (UNESCO/TEOS-10; optional sediment-density effect on buoyancy) | `Hydro/misc_subs.F90` `eqstate` (`ieos_type`, `ieos_pres`, `ddensed`) | Small | Excellent — Rakali already has Wright/Roquet at p=0; add the pressure argument | **UNESCO / Fofonoff & Millard (1983)**; **TEOS-10 (IOC et al. 2010)**, **Roquet et al. (2015)**; Wright **Wright (1997)** |
| 9 | **Lagrangian particle tracking / drifters / oil spill** | `Utility/Particle_Tracking/ptrack{3,4}.f90` (offline); `src/OilSpill-Jung/` (spill setup on top of ptrack) | Small–moderate | Moderate — gather/scatter along trajectories; often done offline anyway | Random-walk particle dispersion **Visser (1997)**; oil weathering **NOAA ADIOS / Lehr et al. (2002)** |
| — | **Two independent sea-ice models on the unstructured coastal grid** — **OUT OF SCOPE by decision** | `src/Ice/` (FESOM EVP/mEVP + FCT, ~2.0 kLOC), `src/Multi_ice/` (CICE/Icepack column physics, ~55 kLOC) | Large | n/a | EVP **Hunke & Dukowicz (1997)**; mEVP **Bouillon et al. (2013)**; Icepack/CICE **Hunke et al.** — *see "Not gaps"* |

---

## 3. Per-gap detail (top items)

### Gap 1 — Sediment transport & morphodynamics  (HIGHEST VALUE)

**What SCHISM has.** Two independent, selectable capabilities:
- **`src/Sediment/` (SED3D, ROMS/CSTMS lineage, ~7.1 kLOC).** Multi-class (sand + mud),
  multi-layer **bed stratigraphy** with active-layer bookkeeping, suspended load (PPM
  semi-Lagrangian settling + erosion/deposition exchange with the bed), and **four selectable
  bedload formulas** (Meyer-Peter-Müller, van Rijn 2007, Soulsby-Damgaard, Wu-Lin). **Bed
  morphology** via the Exner equation solved with a JCG or WENO scheme; **avalanching**;
  **bedform/roughness predictors** (current ripples, sand waves, wave ripples); **wave-current
  bottom stress** (Grant-Madsen / Soulsby / Madsen-Warner); wave-asymmetry/acceleration-skewness
  bedload (Elfrink, Hoefel-Elgar, Dubarbier); a **morphological acceleration factor** and ramp.
  Cohesive erosion uses the Ariathurai-Krone law; **no flocculation**.
- **`src/Sed2d/` (SED2D, ~5.1 kLOC).** Depth-averaged total-load morphodynamics: a menu of
  empirical formulas (Engelund-Hansen, Ackers-White, Soulsby-Van Rijn, van Rijn 2007,
  Camenen-Larson 2011, TRANSPOR2004), Exner bed update (FE + WENO), avalanching, bedform
  roughness. No bed layering, no cohesive sediment.

**What Rakali has** (so this is a *narrowing* gap, not a total blank): a **single-class suspended
sediment** capability on the **structured multilayer path only**. In
`src/tracer/structured/rki_ml_tracers.F90`: `ml_settling_tracer` (constant settling velocity,
deposits into a `bed_sediment` mass reservoir) + `ml_erode_bed_sediment` (**Ariathurai-Arulanandan
1978** critical-shear erosion, Manning-consistent bed stress), wired into
`src/core/coastal/solver/rki_solver.F90` (`use_sediment` knob, once-per-step on the RK2-averaged
state, conserving suspended+bed mass). **What is missing vs SCHISM:** multiple grain classes;
bedload of any kind; **bed morphology feedback into the actual bathymetry** (Rakali's `bed_sediment`
is a passive mass reservoir — the hydro bed elevation never changes, so there is **no
morphodynamics / Exner update**); bed stratigraphy; cohesive flocculation; bedform/roughness
predictors; wave-current stress (Rakali has no waves); and any unstructured-backend sediment
(confirmed absent — `grep` of `src/tracer/unstructured/` returns nothing).

**GPU fit.** Good. Erosion/deposition/settling are per-column; the Exner bed update is a per-cell
hyperbolic update; grain-class loops are inner and vectorizable. The already-shipped
`ml_erode_bed_sediment` is a working DC template. The awkward bits are bedload flux divergence
(needs a conservative face flux like the existing HLL/KNP machinery — reusable) and avalanching
(iterative relaxation — cadence-bounded is fine). Morphology feedback requires letting the
bathymetry `b` evolve, which touches the well-balanced reconstruction — a real but bounded change.

### Gap 2 — Wave coupling  (HIGH VALUE)

**What SCHISM has.** `src/WWMIII/` is a full **third-generation spectral wave model** (action
balance over frequency × direction × unstructured space, ~80 kLOC) with WAM/ST4/ST6 source-term
packages, quadruplet + triad nonlinear interactions, depth-induced breaking, bottom friction,
refraction/diffraction, surface roller. It couples to the hydrodynamics (`wwm_coupl_selfe.F90`) via
either **radiation stress** (Longuet-Higgins-Stewart) or a **3D vortex-force** formulation (Bennis
et al. 2011): Stokes drift (3D), conservative vortex force + Bernoulli head, breaking/roller
non-conservative force, wave streaming bottom term, wave setup (an elliptic Poisson solve), and
wave-enhanced bottom stress. Separately, `Hydro/misc_subs.F90` field-couples to external
**WaveWatch III** through ESMF/OASIS (`USE_WW3`). `src/WW3/` itself is an empty stub.

**What Rakali has.** Nothing. `docs/CAPABILITIES_AND_LIMITATIONS.md` states plainly: wave radiation
stress is "Not in either" coastal backend, and external wave-model coupling (SWAN/WW3) is "not
implemented". A **repo-wide** grep of `src/` for radiation-stress / rad_stress / Stokes /
vortex-force / WWM / WW3 / SWAN returns **zero hits** — this is a clean blank slate, not a partial
implementation to extend. Non-hydrostatic Rakali can *phase-resolve* waves, but that is a different
capability from a spectral wave climate.

**GPU fit — split the problem.** The **spectral wave model** is a poor `do concurrent` target (its
value is a global implicit spectral-space solve with its own PETSc/PDLIB machinery) — porting it is
not the right move; the community path is to *couple* an external wave model (WW3 has GPU efforts;
SWAN; or WWM as a CPU co-process) via field exchange. The **coupling terms** Rakali would need to
consume a wave field — radiation-stress divergence, vortex force, Stokes drift, wave-enhanced bottom
drag, wave setup — are all **per-cell/per-layer** and are an excellent DC fit. Recommended framing:
build the wave-forcing *sink terms* + a field-ingest path first (immediately useful with offline
wave output), treat online two-way coupling as a later phase.

### Gap 3 — Water quality / biogeochemistry  (HIGH VALUE, broad)

**What SCHISM has.** A **generic multi-model tracer framework** (`Core/schism_glbl.F90`,
`natrm=12` tracer models sharing one transport core keyed by `irange_tr`), into which plug: **ICM**
(VIMS Integrated Compartment Model water quality — 3 phytoplankton groups, full N/P/Si cycling,
dissolved oxygen, carbonate/pH, a two-layer benthic sediment-diagenesis flux model, SAV, marsh,
clams; ~4.3 kLOC), **CoSiNE** (Chai carbon-silicon-nitrogen ecosystem, 13 tracers), **EcoSim**
(spectrally-resolved bio-optical ecosystem, up to 62 tracers), the **FABM** coupler (host for any
external FABM biogeochemical model), and **FIB** (fecal-indicator bacteria decay). Also generic user
tracers, water age, and salinity-variance dissipation. Oil spill is offline particle tracking.

**What Rakali has.** A genuine **passive-tracer registry** (`multilayer%register_passive_tracer`)
on both backends — every transport kernel (advection, ALE remap, vdiff, clamp, RK2) already iterates
it, so *transport* of arbitrary tracers "just works" (verified: `rki_solver_unstr.F90` carries ~20
`do it = 1, size(state%multilayer%tracers)` loops, including both advection call sites, CWC and
MUSCL). Shipped tracers: **salinity and temperature only** — two slots vs SCHISM's twelve tracer
*models*. **What is missing:** any **reaction kinetics / source-sink biogeochemistry** — no
nutrient cycling, no phytoplankton/DO/carbonate, no sediment diagenesis, no bacteria decay. Grep of
`src/tracer` + `src/core` for biogeochem/phytoplankton/water-quality returns nothing (one incidental
mention in an ocean README).

**GPU fit — excellent.** BGC reaction kinetics are independent per-cell (per-column) stiff-ODE
integrations — the canonical embarrassingly-parallel workload, ideal for `do concurrent`. Rakali's
registry means the transport half is already solved. The cost is **breadth and validation**: many
state variables, many rate expressions, and the numerical care that stiff kinetics + positivity
demand. A FABM-style coupling (compile-time model selection feeding a per-cell kinetics call) would
let Rakali reuse an entire community model library rather than re-deriving each — the highest-leverage
route, modulo FABM's GPL licensing which must be handled at the coupling boundary.

### Gap 4 — GLS turbulence family via GOTM

SCHISM offers the full generic-length-scale suite through GOTM (`itur=4`, `do_turbulence`): k-ε,
k-ω, Mellor-Yamada 2.5, with Canuto/Kantha-Clayson stability functions, plus a built-in `itur=3`
closure. Rakali coastal ships **PP81, KPP, and a clean-room k-ε** (`docs/CLOSURE_MATRIX.md`,
`src/parameterizations/vertical/`), on both backends and under the semi-implicit mode — but the
matrix note is explicit: "The GLS-generic family (k-ω / MY2.5) is **not yet built**." Because
Rakali's k-ε is already a per-column implicit-Thomas integrator with pluggable stability functions,
adding k-ω / MY2.5 is an **incremental, GPU-friendly** extension (the GLS ψ-equation is one
parameter family over the same two-equation machinery) — not a GOTM port. Re-derive from
Umlauf & Burchard (2003); do **not** vendor GOTM (GPL).

### Gap 5 — Hydraulic structures

`Core/hydraulic_structures.F90` supports prescribed flow transfers, rectangular orifices, weirs
(incl. submerged transition), circular pipes/culverts, radial gates (HEC-RAS-simplified), and
compound weir-over-culvert structures, with time-varying gate operation — connecting non-adjacent
node pairs by imposed face fluxes. Rakali has **none** (`docs/CAPABILITIES_AND_LIMITATIONS.md`:
"Hydraulic structures (weirs, culverts, gates)" listed under "Not in either"). This is
small in code but important for engineering/estuary applications (locks, barrages, culverts). GPU
caveat: the structure flux itself is a trivial algebraic law, but transferring flux between two
distant cells is a scatter — on GPU it needs an atomic or a gather-side formulation to avoid a race,
consistent with Rakali's "no neighbour writes" rule.

### Gap 6 — Vegetation & marsh drag

SCHISM carries a bendable-canopy vegetation drag (`veg_*`, `isav`; a per-height momentum sink with
bending), wave dissipation by vegetation (`wwm_vegdiss.F90`), and dynamic marsh growth/mortality
(`age_marsh`, `create_marsh_*`). Rakali has no vegetation term (grep of coastal + vertical params:
absent). The static drag is the easy win: a per-cell/per-layer quadratic momentum sink
`-½ C_d a |u| u`, a trivial `do concurrent` kernel; re-derive from Luhar & Nepf (2011). Dynamic
marsh morphology is a larger, more specialised effort.

### Gap 7 — Open boundaries, tidal potential, nudging

Beyond the nine OBC types Rakali already has (WALL/OPEN-Flather/TIDAL/NESTED/INFLOW/DISCHARGE/
CLAMPED/SPONGE/CHAPMAN — `rki_boundary_types.F90:18-34`),
SCHISM adds: **earth tidal potential** body forcing (`ntip` constituents, `tip_dp`) and **loading
tide / self-attraction-and-loading** (`iloadtide`); **file-driven 2D tidal harmonics** on open
boundaries (Rakali has per-edge constants only, per the capabilities doc); and **per-tracer OBC
nudging** with incoming/outgoing factors + interior relaxation (`inu_tr`, `*_nudge.gr3`). Rakali's
path explicitly lists "Tidal potential / SAL correction", "2D tidal harmonics from a file", and
open-boundary T,S nudging as absent — confirmed by grep (`tidal_potential`, `love_number`,
`self_attraction`: zero hits). Tidal boundary forcing does support multiple constituents per edge
(`MAX_TIDAL_CONSTITUENTS`), but as compiled-in constants, not a harmonics file. All of these are
per-cell additive forcings — good DC fit. Citations: Cartwright & Tayler (1971) for the potential;
Ray (1998) for SAL.

### Gap 8 — In-situ nonlinear EOS

SCHISM's `eqstate` offers full UNESCO nonlinear density with an optional **pressure correction**
(secant bulk modulus) and an optional **sediment-density buoyancy effect**. Rakali coastal ships a
linear EOS plus opt-in **Wright/Roquet nonlinear** density — but evaluated at **p = 0** (surface/
potential density only); `docs/CAPABILITIES_AND_LIMITATIONS.md` lists in-situ pressure-dependent
evaluation and a TEOS-10 branch as still open. Small, excellent GPU fit — the shared `rki_eos`
already computes Wright; adding the pressure argument and a TEOS-10/GSW path is incremental.
Re-derive from Fofonoff & Millard (1983) / IOC et al. (2010) / Roquet et al. (2015).

---

## 4. "Not gaps" — deliberate design divergence, not missing features

- **Semi-implicit ELM large-`Δt` vs explicit Godunov.** SCHISM is semi-implicit finite-element/FV
  with Eulerian-Lagrangian advection (Zhang et al. 2016), so it takes large `Δt` and damps/smears
  gravity waves, shocks, and sharp fronts. Rakali is explicit HLL/HLLC + KNP (gravity-wave-CFL
  bound) **by design**, buying faithful tsunami/surge/dam-break/hydraulic-jump physics and exact
  lake-at-rest. This is a stated tradeoff, not a deficiency — and Rakali now *also* offers an opt-in
  Casulli semi-implicit free-surface mode (2D + 3D, both backends) for when large `Δt` is wanted.
  A "SCHISM does implicit transport / big timesteps" item is therefore **not** a gap.
  **But** that mode is early-phase and heavily restricted — see §5, which is a *prerequisite*
  for the parity target, not a gap in it.
- **Mixed triangle/quad unstructured grids.** SCHISM supports mixed tri/quad elements (`i34`).
  Rakali's unstructured backend is triangle-only (KNP central-upwind). Real but minor; a design
  choice tied to the central-upwind scheme, not a capability hole for most coastal meshes.
- **LSC2 / SZ vertical coordinates.** SCHISM uses LSC2 (localized sigma, shaved cells) + SZ hybrid.
  Rakali fields a *different* modern stack (sigma / z-sigma / z*-lite / z*-full / z*-sigma via ALE
  remap, PCM..PQM). Different lineage, comparable generality — a divergence, not a gap. (Isopycnal /
  true partial cells remain unimplemented in Rakali, noted honestly in its docs.)
- **Sea ice on the coastal grid.** SCHISM has two ice models on its unstructured grid
  (`src/Ice/`, `src/Multi_ice/`). Rakali has none, and **will not** — ice is out of scope for
  this repo by decision. This is the one place where the parity target is deliberately
  incomplete; call it out rather than let it read as an oversight.
- **CPU/MPI-only maturity.** SCHISM's breadth reflects a serial/MPI CPU heritage. Rakali trades
  breadth for a single GPU-native source; the missing modules are a *roadmap*, and the point of the
  divergence is that when they land they run on the GPU.

---

## 5. Semi-implicit maturity — the gating prerequisite

The stated target is physics parity on **both** time schemes. That makes the Casulli
semi-implicit mode's current maturity a first-order planning constraint rather than a footnote:
it is a Phase-0/1/2 capability, and every restriction below is a hard `logger%error` at config
validation, not a silent degradation. Source of record: `src/core/rki_config.F90`.

| Restriction | Where | Consequence for parity work |
|---|---|---|
| **Single-rank only** — `px*py` must be 1; the η-solve CG iterates interior-only with no per-iteration halo exchange | `rki_config.F90:1053` | **The biggest one.** No SI result can be run at production scale or validated against a multi-rank explicit run. Any SI physics added now is untestable beyond one GPU |
| Boundary types limited to **wall / tidal / clamped** on all four edges | `:1065`, `:1071`, `:1077`, `:1083` | Open/Flather, Chapman, sponge, inflow, discharge, nested are unavailable. Rules out most realistic estuary/shelf configs — and Gap 7 (tidal potential, nudging) cannot be exercised on the SI path at all |
| Multilayer requires **`vcoord_type='sigma'`** | `:1129` | The entire z\*/ALE stack is explicit-only. Gap 1 morphodynamics on a moving bed and any z\*-dependent physics have no SI counterpart |
| Multilayer requires **`nz_layers >= 2`** | `:1109` | The nz=1 fused fast path is not wired to the SI column solve |
| **Incompatible with non-hydrostatic** | `:1037` | Expected and probably permanent (implicit free surface assumes hydrostatic) |
| **`use_sediment` is a hard error** under SI+multilayer | `:1123` | Gap 1 has *no* SI path today. Morphodynamics parity requires lifting this first |
| Requires **`dt_fixed > 0`** | `:1043` | By design — the scheme bypasses the gravity-wave CFL |
| Unstructured SI has **no semi-Lagrangian predictor** (`si_advect`; departure-point location on triangles deferred) | `rki_config.F90`, SI advect guard | The unstructured×semi-implicit cell of the four-path matrix is the weakest of the four |
| `wetdry_bsub_fill='linear'` unsupported with multilayer | `:1148` | Minor; use `'flat'` |

**Planning consequence.** The `docs/plans/README.md` convention — "all four paths:
{structured, unstructured} × {explicit, semi-implicit}" — currently over-promises. Three of the
four cells are production-capable; the semi-implicit column is gated behind single-rank
execution and a three-BC menu. Two honest options, and the choice should be explicit:

1. **Harden SI first** — lift the single-rank restriction and widen the BC menu as a
   prerequisite PR series, *then* apply the four-path rule literally to every physics gap.
2. **Stage the parity work** — implement each physics gap explicit-first, with the SI variant
   as a declared follow-on phase per gap, and record the SI cell as "deferred" rather than
   silently shipping an untestable implementation.

Option 1 is the smaller total effort if the parity target is genuinely both-schemes, because
otherwise every gap in §2 pays the same tax independently.

---

## 6. "Unverified" — checked shallowly or not at all (do not read silence as absence)

- **Sediment-hydrodynamic feedback depth in SCHISM** (bed roughness → drag → flow) was read at the
  module level, not traced end-to-end through `schism_step`.
- **Exact WWMIII↔SCHISM two-way coupling cadence / stability** (how often, sub-cycled?) not traced.
- **FABM licensing mechanics** — I confirmed FABM is a coupler and dual-licensed; I did not analyse
  what a GPL-clean re-coupling would require. Flagged for a licensing review before any port.
- **Rakali's `bed_sediment` conservation under wet/dry** — I read the erosion/settling kernels and
  their solver wiring but did not run them; the mass-conservation claims are from docstrings + code
  structure, not a test run.
- **SCHISM particle-tracking online vs offline** — I confirmed `ptrack` is an offline utility and
  OilSpill builds on it; I did not check whether any online drifter path exists.
- **Non-hydrostatic parity.** Rakali has structured NH (CG Poisson). SCHISM's non-hydrostatic
  pressure solve `solve_jcg_qnon` is **confirmed present** (`src/Hydro/solver_subs.F90`) — but
  existence only; the formulation, mode-splitting, and boundary treatment were not traced, so
  whether this is a gap, a parity, or a divergence is still open.
- **Wave setup / infragravity, tsunami source modules, PaHM parametric hurricane wind** (`sample_inputs/PaHM_inputs`
  seen in the tree) — not investigated.
- **Data assimilation / analysis modes** in SCHISM — not investigated.
- **Multi_ice/Icepack BGC** (the ice-algae biogeochemistry inside Icepack) — noted its existence,
  not characterised. Moot under the current scope decision (no ice).
- **Rakali's semi-implicit physics coverage per closure.** §5 lists the *config-level* guards,
  which are authoritative for what will run. What each SI kernel actually integrates (which
  closures fill `si%nu_iface`, which source terms are in the implicit column vs the explicit
  predictor) was read only where the guards pointed — a kernel-level audit is still owed before
  any physics gap is implemented on the SI path.

---

*End of report. No SCHISM source was reproduced; all algorithms are described for re-derivation from
the cited primary literature.*
