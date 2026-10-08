# Rakali — roadmap

Long-horizon backlog for the solver plus the
**cross-cutting** engineering themes (§13–18: performance, infra, DA, ML, HPC / CS
research). Near-term release work →
[`v0.1.0_checklist.md`](v0.1.0_checklist.md); what works today →
[`CAPABILITIES_AND_LIMITATIONS.md`](CAPABILITIES_AND_LIMITATIONS.md).

Scope tags for collaborator / student projects:

- **`[E]`** Easy — single well-scoped feature, mostly engineering + validation against a known case. ~3–6 months. *Honours / summer.*
- **`[M]`** Medium — multi-component, a genuine research call (scheme choice, parameter study). ~1–2 years. *Masters / first PhD chapter.*
- **`[H]`** Hard — novel methodological / physical contribution; multi-paper material. ~3–4 years. *Full PhD.*
- **`[V]`** Very hard — research programme, often a team. *Postdoc / faculty-led.*

Each item names a *Validate:* target (the examiner-defensible end state). `(tentative)` = judgement call whether Rakali should do it at all. This is a **menu, not a queue** — see the prioritization tiers at the end.

## 1. Vertical mixing & boundary layer

- [~] `[H]` **GLS turbulence closure** — k-ε shipped (structured, clean-room; see the checklist).
  Remaining: k-ω + Mellor-Yamada 2.5 via the unified Umlauf-Burchard (2003)
  generic-length-scale switch; stratification-aware stability functions
  (Canuto / Galperin / Kantha-Clayson); unstructured port. Implicit tridiagonal
  vertical diffusion now exists on the k-ε path — generalize to PP81/KPP if needed.
  *Validate:* GOTM k-ε/k-ω/MY cross-comparison; Liverpool Bay ROFI.
- [ ] `[M]` **Matched cubic shape function (KPP v2)** — continuous K at OBL base.
  *Validate:* smooth K(z), no jump in d²K/dz²; sharper diurnal thermocline.
- [ ] `[M]` **Double-diffusive mixing** — salt fingers / diffusive convection,
  optional KPP extension. *Validate:* Schmitt 1981 finger regime; C-SALT climatology.
- [ ] `[H]` **KPP-Langmuir** — Stokes-drift enhancement to w_m; needs a wave model.
  *Validate:* McWilliams 1997 LES Langmuir-cell mixing.
- [ ] `[E]` **Shortwave penetration** — `heat_flux_sw` + exponential decay Q(z);
  currently all heat absorbed at the surface. *Validate:* Paulson-Simpson 1977 Jerlov types.
- [ ] `[E]` **Convective adjustment** — instantaneous mixing of unstable columns
  as a pre/alongside-KPP safety net. *Validate:* unstable column → neutral in one
  step, column heat/salt conserved.

## 2. Equation of state and density

- [ ] `[E]` **Nonlinear EOS — UNESCO 1980** (Jackett-McDougall polynomial).
  *Validate:* ≤1e-3 kg/m³ vs UNESCO tables; linear EOS recovered when linearised.
- [ ] `[M]` **Nonlinear EOS — TEOS-10** (GSW / IAPWS-09); adds the GSW-Fortran dep.
  *Validate:* GSW check-value table to machine precision.
- [ ] `[E]` **Pressure dependence** in the EOS (not just z-depth proxy); deep-water only.
  *Validate:* 4000 m static column vs WOA13; bit-identical to depth-only in <200 m.
- [ ] `[H]` **Non-Boussinesq option** (density in continuity) *(tentative — Boussinesq is fine for coastal)*.
  *Validate:* mass (not just volume) budget closes; steric height vs 1D theory.

## 3. Multilayer extensions

- [ ] `[E]` **Dynamic vertical workspace** — `BPG_NZ_MAX=20` fixed arrays → automatic
  arrays sized by nz (audit per-thread GPU stack first). *Validate:* bit-identical at
  nz≤20; high-nz (40/80) runs within stack budget, no allocator churn.

The geometric vcoord families (sigma / z\* / hybrids) are shipped; true
isopycnal (RHO) coordinates and partial cells are not implemented.

## 4. Unstructured mesh solver — parity with structured

- [ ] `[M]` **Per-layer Coriolis** (both backends still depth-averaged; pairs with vector-invariant refactor).
  *Validate:* geostrophic-adjustment / inertial-oscillation analytical test; energy conserved.
- [ ] `[E]` **CFL-based velocity truncation** (CFL>0.8 clipping).
  *Validate:* 1 m-resolution lock-exchange doesn't blow up; clip telemetry matches structured.
- [ ] `[E]` **Per-layer bottom friction.** *Validate:* Manning bottom-layer friction matches structured on a shared shelf.
- [ ] `[M]` **Fused multilayer kernel.** *Validate:* bit-identical vs unfused; reduced launch overhead.
- [ ] `[M]` **`BC_NESTED`** (read parent NetCDF, apply at child boundary).
  *Validate:* fine triangular nest in a structured Sydney Harbour parent; tidal phase agrees at boundary.
- [ ] `[M]` **`BC_SPONGE`** (edge-distance-decaying relaxation).
  *Validate:* outgoing-wave reflection coefficient vs analytical sponge theory.
- [ ] `[H]` **Non-hydrostatic on triangles** — CG Poisson on unstructured.
  *Validate:* solitary-wave propagation matches the structured NH solver at equal resolution.

## 5. Grid coupling and nesting

- [ ] `[H]` **Two-way structured nesting** — child→parent feedback (de-aliasing, conservation).
  *Validate:* mass/momentum conserved across the flux exchange; no aliasing pile-up.
- [ ] `[H]` **Structured ↔ unstructured coupling** — triangular child inside a structured parent. **Major PR.**
  *Validate:* structured-only and unstructured-only refs agree with the coupled run on the shared region.
- [ ] `[H]` **Local sub-stepping on heterogeneous meshes** — removes the "sliver triangle dictates global dt" pathology (cross-listed §18).
  *Validate:* river-plume at 5–10× larger global dt, identical conservation budget.
- [ ] `[M]` **Telescoping grids** — >2 nested levels. *Validate:* 3-level telescope reproduces 2-level + single-domain refs.
- [ ] `[V]` **Overset / Chimera grids** *(tentative — uncommon in this domain)*.

## 6. Wetting and drying

- [ ] `[E]` **Adaptive wet/dry thresholds** — 2D mask / bathymetry-dependent (rocky shore + mud flat). *Validate:* Sydney Harbour low-water, no spurious dry-front oscillation.
- [ ] `[M]` **Wet/dry in the multilayer solver** — partial-column transitions (some layers wet, some dry). *Validate:* Parramatta spring-low-tide, conservation through transitions, no blow-up.
- [ ] `[H]` **Subgrid wetting/drying** — fractional cell wetness for coarse storm-surge grids. *Validate:* Ike/Katrina-class surge vs fine-grid ref.
- [ ] `[H]` **Tidal-flat sediment drying** — ties into §9 sediment. *Validate:* Wadden Sea consolidation / erodibility.

## 7. Time integration

- [ ] `[E]` **SSP-RK3** (Shu-Osher) — better stability margin, one more RHS eval. *Validate:* convergence-order test; stable above current RK2 CFL.
- [ ] `[E]` **SSP-RK4** — diminishing returns until higher-order spatial lands.
- [ ] `[M]` **IMEX** (implicit vertical, explicit horizontal) — large dt in well-mixed columns. *Validate:* strong-vmix 1D column at dt≫explicit limit, bounded error.
- [ ] `[E]` **Tracer-CFL adaptive dt** (not just barotropic). *Validate:* river plume, no tracer overshoot.

## 8. Flux schemes and spatial order

- [ ] `[H]` **WENO reconstruction** (3rd/5th order) — pays off at high resolution; long-term (cross-links the checklist's shared WENO item). *Validate:* order-of-accuracy on smooth solution, non-oscillatory at a dam break.
- [ ] `[M]` **MUSCL-Hancock** (predictor-corrector FV) — sometimes more robust at fronts than minmod+RK. *Validate:* convergence + oblique-front robustness.
- [ ] `[E]` **Flux limiters beyond minmod** — van Leer / van Albada / superbee, as config options. *Validate:* second-order preserved; limiter-comparison plots.
- [ ] `[M]` **Roe solver** with entropy fix *(probably skip — HLL family is fine)*. *Validate:* transonic rarefaction (entropy fix); Toro's bank.
- [ ] `[V]` **Discontinuous Galerkin** *(tentative — a different solver paradigm)*.

## 9. Tracer transport and biogeochemistry

- [ ] `[E]` **Passive dye tracer** (Lagrangian viz) — trivial on the registry. *Validate:* dye blob advects without mass loss.
- [ ] `[E]` **Age tracer** (Deleersnijder) — accumulates at 1 s/s. *Validate:* age grows linearly in a closed basin; channel flushing-time vs analytic.
- [ ] `[M]` **Generic N-tracer module** — user-configured passive/reactive tracers with per-tracer settling, decay, surface/bottom flux (SCHISM's `GEN`). *Validate:* per-source unit tests; namelist round-trip.
- [ ] `[M]` **FABM coupling** — thin binding to the Framework for Aquatic Biogeochemical Models (ERSEM, PISCES, NPZD…); parity with SCHISM's tracer-module zoo for far less effort than reimplementing any one. *Validate:* FABM 1D NPZD/ERSEM-light vs GOTM-FABM reference.
- [ ] `[H]` **Biogeochemistry v1** — NPZD (O₂, NO₃, phyto) source terms in the tracer pipeline. *Validate:* Fasham 1990 seasonal 1D cycle.
- [ ] `[H]` **Biogeochemistry v2** — CO₂ system (TA, DIC, pH), carbonate chemistry. *Validate:* PISCES carbonate-chemistry cases; alkalinity conserved under dilution.
- [ ] `[M]` **Higher-order tracer advection** — MPDATA / FCT; current upwind is diffusive. *Validate:* Smolarkiewicz rotating cone; monotonicity at fronts.
- [ ] `[E]` **Monotonic tracer limiter** — no new extrema (important for BGC). *Validate:* step advection stays in range; conservation to machine precision.

## 10. Boundary conditions and forcing

- [ ] `[E]` **NetCDF-driven heat flux + E−P** — wire into the existing `forcing_file` pipeline. *Validate:* ERA5 1D column seasonal SST; bit-identical to the programmatic path.
- [ ] `[M]` **Bulk aerodynamic flux module** (Zeng 1998 / COARE 3.0-3.5) — compute τ, sensible/latent heat, evap from `(U10, T_air, q_air, SST, SLP, rad)`. Makes ERA5/BARRA-driven runs a one-liner. *Validate:* COARE 3.5 reference fluxes on identical inputs.
- [ ] `[E]` **Vegetation drag** — `F_veg = -½C_D·a·|u|·u`, per-cell mask + height (emergent vs submerged). Salt-marsh / mangrove. *Validate:* Nepf 1999 vegetated-channel profile.
- [ ] `[M]` **Hydraulic structures** (weirs, culverts, gates) — internal-BC kernel modifying fluxes between cell pairs via Q-h relations. *Validate:* broad-crested weir Q-h (Bos 1989); levee-breach (Soares-Frazão & Zech).
- [ ] `[E]` **River discharge schedules** — time-varying Q(t) from CSV/NetCDF. *Validate:* prescribed Q(t) integrates to correct cumulative volume.
- [ ] `[E]` **T/S open-boundary relaxation** (nudging) — common in regional setups. *Validate:* boundary T,S diffuse in at the right timescale.
- [ ] `[M]` **Tidal potential forcing** (body force + self-attraction & loading). Basin/regional tides. *Validate:* Australia-wide barotropic vs TPXO M2/S2/K1/O1.
- [ ] `[E]` **2D tidal harmonics** — spatially-varying amp/phase from file. *Validate:* boundary-stage harmonics match the input file.
- [ ] `[M]` **Wave radiation stress** from a wave model (needs §12). *Validate:* longshore current vs Longuet-Higgins 1970.
- [ ] `[H]` **Ice shelf / iceberg melt flux** *(tentative — polar only)*.

## 11. Output, diagnostics, and analysis

- [ ] `[E]` **Diagnostic fields** — vorticity, Okubo-Weiss, Ertel PV, KE density, computed on demand. *Validate:* Taylor-Green / Rossby-wave vorticity vs offline Python.
- [ ] `[E]` **Energy budget diagnostics** — PE, KE, dissipation. *Validate:* lock-exchange PE→KE budget closes.
- [ ] `[M]` **Lagrangian drifter tracking** — forward-integrate particles, output trajectories. *Validate:* solid-body rotation traces circles.
- [ ] `[E]` **Moving gauges** — for vessel / AUV tracks. *Validate:* prescribed track over an analytical field returns the right series.

## 12. Waves — in-model and coupled

Four entry points, increasing ambition. **Recommended ordering:** start 12D
file-based coupling → harden 12A validation → build 12B if a lighter option is
wanted → defer 12C unless a multi-year programme justifies it.

- **12A — Phase-resolving via the NH solver** (the SWASH-equivalent path). The
  NH Poisson solver already propagates surface gravity waves on fine grids; the
  work is validation + exposure as a "wave mode": `[M]` solitary-wave run-up
  (Synolakis 1987), `[E]` standing wave in a closed basin, `[M]` Dingemans-bar
  benchmark, `[H]` Boussinesq dispersion correction (Madsen-Sørensen / Nwogu),
  `[M]` wave-maker BC, `[M]` absorbing BC tuned for gravity waves, `[H]` breaking
  treatment. *Realistic:* ~2–4 weeks + validation; not a general wave model.
- **12B — Parametric wave tracers** — evolve `Hs`, `Tp`, `θ_m` as tracers with
  wind input / whitecapping / bottom-friction source terms; derived Stokes drift,
  radiation stress, wave-enhanced bottom stress feed back into momentum + mixing.
  `[M]` source terms, `[M]` directional refraction, `[E]` Stokes/radiation-stress/
  bottom-stress diagnostics, `[M]` KPP-Langmuir integration. *Realistic:* ~500–1000
  LOC; runs at native resolution (no phase-resolution refinement).
- **12C — Full spectral model** `[V]` — wave-action balance on a 2D spectrum per
  cell (propagation, wind input, whitecapping, quadruplet/triad interactions,
  depth breaking). A 1–2 FTE-year project (SWAN/WW3/WWM are ~100k LOC each).
  *Recommend against in favour of 12D* unless a deliberate multi-year commitment.
- **12D — External wave-model coupling** (the industry-standard compromise, the
  baseline recommendation): `[E]` file-based offline (SWAN NetCDF → Rakali),
  `[H]` OASIS/ESMF two-way, `[H]` library-linkage in-process; consume wave bottom
  stress + radiation stress + breaking-TKE. *Validate:* ROMS-SWAN-equivalent on a
  shared storm-surge case.

---

_The themes below (§13–18) are **cross-cutting** engineering themes._

## 13. Performance and scalability

(Several items cross-listed in §18 as CS/HPC research projects.)

- [ ] `[H]` **Dynamic load balancing** for heterogeneous wet/dry fractions — Hilbert repartition at checkpoints. *Validate:* Sydney Harbour spring-neap >30% dry, imbalance metric under threshold, throughput gain on 16+ GPUs.
- [ ] `[M]` **Mixed-precision tracers** *(tentative)* — fp32/bf16 where conservation allows; bandwidth win for many-tracer BGC. *Validate:* error-vs-precision on the analytical suite; long-rollout drift fp64 vs mixed.
- [ ] `[E]` **Memory-pool reuse audit** — quantify wins from a unified pool vs status quo. *Validate:* `acc enter data` count + step-time variance on a 24-h rollout.

## 14. Validation and testing

- [ ] `[M]` **Castro et al. SWE dam-break benchmark suite.** *Validate:* published comparison table, tolerance per case.
- [ ] `[M]` **Sydney Harbour observational validation** — tide-gauge + ADCP. *Validate:* M2/S2/N2/K2/K1/O1/P1/Q1 harmonics vs BoM; ADCP velocity RMS.
- [ ] `[M]` **Parramatta River** — thermal-stratification observations. *Validate:* stratification onset + salt-intrusion length vs climatology.
- [ ] `[M]` **Oresund exchange flow** — baroclinic validation. *Validate:* front speed + stratification vs Burchard-Bolding.
- [ ] `[M]` **Wave benchmarks** — Synolakis run-up, standing wave, Dingemans bar (prereq for claiming §12A). *Validate:* see §12A.
- [ ] `[E]` **Golden-output regression harness** — bitwise vs reference NetCDF in CI. *Validate:* seeded-perturbation detection; CI green on curated gold examples.
- [ ] `[M]` **Performance regression CI** — track throughput over time. *Validate:* seeded-slowdown detection on the GPU runner.

## 15. Infrastructure and tooling

- [ ] `[E]` **Python namelist mapping gap closure** — automated test that every Fortran key has a Python mapping; CI fails on a missing key.
- [ ] `[E]` **C API setters for EOS coefficients** (`alpha_T`, `beta_S`, `rho_0`, refs, clamp bounds) — for sensitivity / DA. *Validate:* set/get round-trip; mid-sim change reflects next step.
- [ ] `[E]` **Python notebook tutorials** — one per canonical workflow; execute end-to-end in CI.
- [ ] `[M]` **End-user documentation** — physical-interpretation / config / pitfalls guide (vs developer docs). *Validate:* domain-user review; namelist-key coverage.
- [ ] `[E]` **Container image** (Docker / Apptainer) — reproducible nvhpc + netcdf + python env. *Validate:* runs the test suite on a clean host.
- [ ] `[E]` **Schema-validating namelist reader** — reject unknown keys loudly at read time (currently tolerant). *Validate:* typo'd / missing / mistyped-key unit tests.

## 16. Data assimilation hooks

- [ ] `[M]` **Nudging operator** — relaxation toward observations over τ; a new kernel. *Validate:* state converges at the configured timescale.
- [ ] `[M]` **Ensemble support** — member-parallel runs (OpenMP task) for EnKF / PF. *Validate:* N-member mean/spread on a stochastic-perturbation test.
- [ ] `[M]` **Observation-operator hooks** — bilinear interp at obs locations, write model-at-obs (half-done via gauges). *Validate:* model-at-obs vs offline interpolation.
- [ ] `[V]` **4D-Var adjoint** *(tentative)* — needs derivative-capable kernels; probably stays a forward model coupled to an external DA framework.

## 17. Machine learning integration

The natural ML target is the **overall time evolution** (MeshGraphNets, FNO,
neural operators), not individual kernels — the `set_bc → advance → read` loop is
already the autoregressive rollout pattern, and results are pure numpy above the
Fortran side. **Defer the abstraction until a concrete model exists** (premature
abstraction risks missing real requirements).

- [ ] `[E]` **`rakali.ml.mesh_features(mesh)`** — precompute edge_index, centroids, areas, normals, boundary tags (everything a GNN needs), cached per mesh.
- [ ] `[E]` **Training-data export** — `snapshot_loader` reading UGRID snapshots + BC history into PyTorch tensors; lives in a `rakali[ml]` extra so the core doesn't pull torch.
- [ ] `[M]` **`make_ml_surrogate(model, dt, …)`** — symmetric with `make_simulation`, same stepping/state API; needs a `SimulationProtocol` contract first.
- [ ] `[M]` **Hybrid stepping** (N physics : M surrogate steps in one loop) + a benchmarking harness scoring a surrogate vs the solver on rollout error / conservation / wall-time.
- *Caveat:* Fortran GPU + PyTorch in one process contend for device memory; likely flow is offline trajectory generation → physics-free inference.

## 18. Computer-science / HPC research projects

Rakali is a production-grade GPU codebase — a credible testbed for systems /
parallel-computing research that doesn't need a domain scientist to lead.

- [ ] `[H]` **Cross-vendor portability study** — stdpar / OpenACC / OpenMP-target / SYCL / HIP / Kokkos on a hot-kernel set across NVIDIA / AMD / Intel GPUs. *Validate:* publishable perf/portability matrix (perf, LOC, compile-time, failure modes) on Frontier / Aurora / Perlmutter; conservation within ULP across backends.
- [ ] `[M]` **Hand-written CUDA hot kernels** vs the `do concurrent`+OpenACC path — quantifies the stdpar portability tax. *Validate:* bit-identical output; speedup + occupancy/register attribution for any gap.
- [ ] `[M]` **NVHPC inlining study** — formalise the `_impl` cross-TU duplication pattern; tooling-detect inlining failures and auto-suggest the refactor. *Validate:* predicted-vs-measured speedup correlation; release the analyser.
- [ ] `[H]` **Local time stepping on heterogeneous meshes** (cross-listed §5) — conservative substep flux exchange. *Validate:* river-plume at 5–10× global dt, identical conservation.
- [ ] `[M]` **GPU kernel auto-tuning** — tile sizes, gang/vector dims, register limits per architecture, generated at build / first run. *Validate:* reproducible speedup over hand-picked defaults on ≥3 architectures.
- [ ] `[M]` **Reproducibility study** — bit-for-bit across NVHPC versions, GPU architectures: what's preserved, what would it take. *Validate:* publishable matrix of which config dimensions break reproducibility.
- [ ] `[E]` **CI: real GPU smoke run** — self-hosted GPU runner on every PR (GPU tests are local today). *Validate:* green CI; flake rate <1% over 100 builds.
- [ ] `[V]` **Differentiable Rakali via Enzyme** — LLVM-level AD through the Fortran kernels; tangent-linear + adjoint without a rewrite. Unblocks 4D-Var (§16). *Validate:* finite-difference + tangent-linear consistency on a 1D test.
- [ ] `[H]` **ML-guided adaptive mesh refinement** — model proposes refinement zones from coarse-run features. *Validate:* error-vs-cost Pareto vs uniform refinement.

---

## Prioritization notes (coastal)

A menu, not a queue. A rough ordering biased toward "what Rakali needs most":

**Tier 1 (finish what's started, high value):**
- NetCDF-driven heat flux / E−P
- Bulk aerodynamic flux module (Zeng / COARE) — closes the SCHISM atmospheric-coupling gap; makes ERA5/BARRA runs a one-liner
- Nonlinear EOS (UNESCO) — small scope, high value once thermal forcing is real
- GLS remaining (k-ω + MY via the unified switch; stratification-aware stability functions) — k-ε already shipped

**Tier 2 (major coastal features):**
- Unstructured multilayer parity: per-layer Coriolis, velocity truncation, per-layer bottom friction, fused kernel (§4)
- BC_NESTED + BC_SPONGE on unstructured (unblocks unstructured nesting)
- Structured ↔ unstructured nesting
- Sediment: unstructured parity, then multi-class + bedload (§9)
- Generic N-tracer module — the path to BGC / sediment classes / dye / age through one abstraction
- Vegetation drag — small code, unlocks salt-marsh / mangrove
- SSP-RK3
- Waves: harden §12A (solitary + standing + Dingemans validation; wave-maker + absorbing BC)

**Tier 3 (ecosystem features):**
- FABM coupling — preferred path to BGC parity without building the modules
- Hydraulic structures (weirs / culverts / gates) — urban flood + levee-breach
- Waves: external SWAN coupling (§12D file-based first); parametric in-model (§12B)
- Lagrangian drifters; diagnostic vorticity / OW / PV

**Tier 4 (polish, validation, long tail):**
- Golden-output regression harness; Sydney Harbour observational validation
- End-user documentation
- Higher-order tracer advection (MPDATA)
- Waves: Boussinesq dispersion correction in the NH solver
- ML surrogate scaffolding (§17) — when a concrete model lands, not before

**Deferred / tentative:** anything marked `(tentative)`.
