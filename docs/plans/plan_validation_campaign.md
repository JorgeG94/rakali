# Coastal Validation Campaign — Executable Plan

**Scope:** this repo is the coastal path and nothing else — `sim_type='coastal'` (structured
A-grid HLL/HLLC + unstructured triangular KNP; 2D barotropic, 3D hydrostatic multilayer,
non-hydrostatic, opt-in Casulli semi-implicit). No ocean dyn-core, no ice.

**SI coverage caveat (2026-08-07):** semi-implicit cases are constrained to wall/tidal/clamped
boundaries, sigma vcoord, `nz_layers >= 2`, and **single-rank execution** — see
`schism_parity_report.md` §5. Validation cases that need open/Chapman boundaries, z\*, or a
multi-rank run have no SI variant today; mark them explicit-only rather than silently skipping.

**Input verified.** I re-read `coastal_physics_gaps.md` §5 and cross-checked the
source. Confirmed: no Carrier–Greenspan, no dry-bed dam-break front test, no standalone C-property
assertion, no flow-over-bump, no structured convergence-order test; lock exchange at ±50 % vs Benjamin;
Kato–Phillips as a 1.6–2.4 ratio band with its √t law living OFF-REPO in Python. Corrected: unstructured
tracer advection is MUSCL piecewise-linear (2nd order), not 1st. **New finding:** a coastal salt/heat
*drift print* already exists (`rki_coastal_console_stats.F90` → `total_salt`/`total_heat` via
`compute_layer_sums`); what is missing is a *closed budget with source accounting* that a test can assert
against zero — Phase 1 builds that, it is not greenfield.

---

## 1. Executive summary

Nine PR-sized phases. Two are **infrastructure** that later cases consume (a conservation diagnostic and
a convergence-order harness); seven are **analytical benchmarks** with closed-form oracles brought
in-repo. The campaign converts the coastal path from "per-kernel exact-arithmetic correct, conservation
strong, but analytical-accuracy + convergence thin" into "the standard shallow-water and stratified
benchmarks all pass at defensible, order-justified tolerances." Ordering is credibility-per-effort: the
cheapest high-signal reviewer checks (well-balancedness, dry-bed front, convergence order) come first;
the two expensive stratified re-tightenings (lock exchange, Kato–Phillips) come last because they need
the infrastructure and may expose a genuine hydrostatic-model ceiling.

**Rough size:** ~2 infra PRs (M each) + 7 benchmark PRs (mostly S–M) ≈ 9 PRs, ~3–4 weeks of focused work.

**Highest-value single item:** **Phase 3 — the well-balanced C-property assertion.** It is the cheapest
(a static lake-at-rest over irregular bathymetry, exact-zero oracle, both backends, ~half a day) and the
highest-signal: a domain scientist's first reflex on any FV coastal solver is "is it well-balanced?", and
today that property is only *implied* by other passing tests, never asserted standalone. It converts an
Audusse claim into a machine-epsilon proof.

---

## 2. Dependency graph

```
Phase 1  Conservation diagnostic (infra) ─────┐
Phase 2  Convergence-order harness  (infra) ─┐ │
                                             │ │
Phase 3  C-property (well-balanced)          │ │   (no dep — do first)
Phase 4  Dry-bed dam break ──────────────────┤ │   (uses Ritter oracle; feeds P7)
Phase 5  Carrier–Greenspan runup ────────────┘ │   (uses conservation diag for wet/dry mass check)
Phase 6  Flow over a bump ──────────────────────┤ (uses conservation diag: discharge closure)
Phase 7  Structured HLL convergence order ──────┘ (uses P2 harness; smooth case)
Phase 8  Lock-exchange tightening ──────── uses P1 (salt closure) + P2 (refinement)
Phase 9  Kato–Phillips in-repo oracle ──── uses P1 (heat/tracer closure)
```

**Sequence rule:** land Phase 1 and Phase 2 first (nothing else should invent a private mass-closure or
refinement loop). Phases 3, 4, 6 are independent and can run in parallel after infra. Phase 7 needs P2.
Phases 8–9 last.

---

## 3. Per-phase detail

Registration boilerplate (identical every phase, stated once): add `test_<name>.F90`; register the row
`"test_<name>|test_<name>|collect_<name>_tests|<name>|coastal"` in `tests/CMakeLists.txt` (5th field =
regime = `coastal`), add the target to the `add_executable`/target list block, and wire the
`use … collect_… ` + `new_testsuite("<name>", …)` pair in `tests/tester.F90`. New closed-form oracles are
added as sources to the `rakali_analytical` STATIC lib (alongside `ritter_dambreak.F90`,
`thacker_bowl.F90`) so every test can link them. Gates every PR: `pre-commit run --all` → `fortitude
check` → `ctest -R rakali` (never `-jN`). New knobs default OFF (bit-identity preserved).

---

### Phase 1 — Continuous salt / heat / mass conservation diagnostic (INFRASTRUCTURE) — **PARTIALLY DONE 2026-08-07**

**Delivered** as `tests/test_coastal_conservation.F90` (5 assertions). Verified on gcc 15.1.0 CPU
and NVHPC 26.5 GPU (cc70/V100): 119/119 `ctest -R rakali` green on both (after the BPG fix below).

**No new budget module was written, and the plan's premise needs correcting.**
`rki_console_stats` already defines `conservation_budget_t` carrying the exact residual identity
`leak_Q = (Q_total − Q0) + out_Q − src_Q`, and `rki_solver` already ships
`accumulate_boundary_mass_out` / `accumulate_mass_src` closing the **mass** budget to round-off.
The genuine hole is narrower than "no closure exists": the `salt_out`, `salt_src`, `heat_out`,
`heat_src` fields of that type are **never populated anywhere in `src/`**, so `salt_active` /
`heat_active` stay `.false.` and both quantities silently fall back to raw drift. The test forms
the source integral itself, which makes the property assertable today; wiring those accumulators
into the solver is a separate source-touching change and remains **open**.

> ### ✅ DEFECT FOUND AND FIXED — spurious sigma BPG over sloping bathymetry
>
> A closed box (walls, zero forcing) over a **sloping** bed did not conserve mass, while a flat
> bed was exact. Bisection isolated the trigger: the leak appeared only when the EOS returned a
> density differing from `rho_0`, and vanished when S and T were seeded exactly at
> `S_ref`/`T_ref` so the EOS returned `rho_0` exactly. A *spatially uniform* density offset from
> the reference was enough.
>
> **Cause.** `ml_baroclinic_pressure_gradient` ("sigma") differences the integrated density
> anomaly at constant **sigma index**: `dp_dx = (p_int_ip1(k) − p_int_im1(k))/2dx`. Over a slope,
> layer `k` sits at a different physical depth in each column, so for a uniform anomaly Δρ the
> difference is `Δρ·(z_mid_ip1 − z_mid_im1) ≠ 0` — a fictitious force proportional to
> (uniform anomaly) × (bed slope). It drove real flow, and the flow leaked mass.
>
> **Fix.** `bpg_method` now defaults to `"jacobian"` (Shchepetkin & McWilliams 2003), which
> differences along constant-z surfaces and satisfies the identity by construction. Both backend
> state defaults (`use_bpg_jacobian`) flipped to `.true.` to match. Kernel-level measurements,
> single step, uniform density over a sinusoidal bed:
>
> | BPG method | sloping bed | flat bed |
> |---|---|---|
> | sigma | 1.28e-03 | 0.0 |
> | jacobian | **2.09e-18** | 0.0 |
>
> End-to-end, 400 steps: relative mass change 1.65e-06 → **1.28e-14**.
>
> Pinned by `tests/test_bpg_uniform_density.F90`, which asserts the identity on the kernel with no
> solver run, and *also* asserts that the sigma method still violates it — so if someone repairs
> sigma, the test fails and prompts revisiting the default rather than silently rotting.
>
> **Bit-identity note:** this changes the default, so multilayer results move for any config that
> did not set `bpg_method` explicitly. `sigma` is retained for reproducing legacy runs. Full suite
> is 119/119 green on both toolchains after the change — no test depended on the old default.

> ### ✅ ALSO FOUND AND FIXED — tracer non-conservation under advection
>
> Independent of the BPG defect and untouched by fixing it. Two further causes, both now closed.
>
> **(1) The clamp destroyed tracer when the column mean left the bounds.** `ml_clamp_tracer` clips
> per-layer concentration into `[tr_min, tr_max]` then redistributes with
> `frac = min(1, net/cap_to_add)`; the `min` saturates and the residual `net − cap_to_add` was
> silently destroyed. That saturation condition is *exactly equivalent* to the column mean leaving
> the bounds — `mean ≤ tr_max ⟺ net ≤ cap_to_add` — so the two-pass scheme is provably conservative
> when the mean is in range and provably cannot be when it is not. The kernel now skips
> out-of-range columns, making it conservative for every input. Both backends fixed.
>
> It fired constantly in practice because the default `S_max = 40 PSU` sits at the top of the real
> seawater range, while the mode-split tracer/thickness inconsistency (`hTr` advected with per-layer
> CWC fluxes but `h_layer` overwritten as `dsig(k)·H_new`) pushes per-layer S above the bound at any
> sharp front.
>
> **(2) Closed walls leaked tracer.** `ml_apply_layer_bcs` distributed ghost momentum
> depth-uniformly, so at a wall `v_face(k) = ½(hu_ghost(k) + hu_int(k))` was nonzero *per layer*
> whenever the interior column was sheared — a shear-free ghost cannot cancel a sheared interior
> layer by layer. Water mass still summed to zero flux (hence exact), but each layer traded tracer
> with the ghost at its own concentration. Wall edges now mirror the interior per-layer profile
> (normal negated, tangential mirrored, matching `rki_boundary`), so `v_face(k)` vanishes identically.
>
> Lock exchange, flat bed, 40×40×6, walls, zero forcing, 400 steps, max|hu| ≈ 0.1 m²/s:
>
> | state | rel. mass | rel. salt | rel. heat |
> |---|---|---|---|
> | before both fixes | 1.5e-14 | −8.85e-06 | +4.60e-08 |
> | clamp fix only | 1.5e-14 | +4.85e-08 | +4.71e-08 |
> | **both fixes** | 1.5e-14 | **0.0** | +8.5e-15 |
>
> Salt is now *exactly* conserved, not merely small. Pinned by
> `tests/test_tracer_conservation.F90` (three kernel-level clamp assertions plus the end-to-end
> lock exchange, with a fixture guard that the gravity current actually develops).
>
> **Contract change:** the clamp no longer forces a uniform out-of-range column onto the bound,
> because doing so requires deleting mass. `test_multilayer_unstr`'s clamp cases were asserting
> exactly that behaviour and have been rewritten to use stratified columns with an in-range mean —
> the regime the clamp exists for and can serve honestly.
>
> **P8 and P9 are no longer gated by this.**

> ### ✅ THIRD DEFECT FOUND AND FIXED — closed walls passed mass over a sloping bed
>
> Localised with the solver's own `mass_out` accumulator, which is built as `−Σ_interior flux_h·area`
> — i.e. the mass change the *flux path* accounts for. In every failing case
> `(M₁−M₀) + mass_out` closed to ~1e-14, so nothing outside the flux path was touching `h`: the mass
> was genuinely leaving through the **boundary**, and the model's own budget already knew it. The
> sign also oscillated with step count (−2.2e-7 → −6.3e-7 → +1.4e-6), i.e. sloshing in and out
> rather than a monotone drain.
>
> **Cause.** `bc_wall_edge` mirrored the free surface and negated normal momentum, but the ghost
> **bed** came from the bathymetry setter — and a formula setter fills ghost rows with the
> *continuation* of the topography, not its reflection. (It must fill them: leaving `b`/`h_layer`
> at zero there makes the EOS fall back to ρ₀.) A wall *is* a mirror, so an unmirrored bed makes the
> reconstruction across the wall face asymmetric. At rest this is harmless — Audusse still returns
> exactly zero flux, which is why the C-property held to 1e-12 over irregular bathymetry — but under
> flow it yields a non-zero HLL mass flux through a closed wall.
>
> **Fix.** `bc_wall_edge` now mirrors the bed as well (`b(ig,jg) = b(im,jm)`), making the whole
> reconstruction symmetric. Done in the per-step BC rather than once at setup so it survives a bed
> that *evolves* — which is exactly what the Exner update in the morphodynamics plan will do.
> Wall/sponge edges only; open/tidal/inflow leave `b` untouched, where the domain genuinely continues.
>
> | fixture, 400 steps | rel. mass change | `mass_out`/M |
> |---|---|---|
> | sloping bed, formula ghost `b` | 1.364e-06 | −1.364e-06 |
> | sloping bed, **mirrored ghost `b`** | **−2.14e-14** | −1.1e-22 |
>
> `max|hu|` is unchanged (4.35e-3 → 4.49e-3), so the physics is identical; only the spurious wall
> flux is gone. With the fix in the BC, a fixture that writes formula bathymetry into its ghosts and
> one that pre-mirrors it now give **bit-identical** results.
>
> **Phase 1's fixtures are no longer flat-only.** `test_coastal_conservation.F90` now asserts mass
> and salt conservation over a sinusoidal bed, plus a direct wall-tightness assertion on `mass_out`
> and a budget-closure assertion — the same two-part diagnostic that localised the defect.
>
> **P5, P6 and P8 are no longer gated.**

> ### Four-path coverage of the conservation work
>
> All three fixes above were *measured* on the structured explicit path only. Two of them (the
> jacobian BPG default, the clamp guard) reach the other paths by shared kernels — which is exactly
> the "by construction" reasoning that let all three defects hide until someone measured.
> `tests/test_conservation_parity.F90` now measures two more paths.
>
> | path | fixes reach it | conservation measured |
> |---|---|---|
> | structured × explicit | ✅ | ✅ `test_coastal_conservation` |
> | structured × semi-implicit | ✅ | ✅ `test_conservation_parity` |
> | unstructured × explicit | ✅ | ✅ `test_conservation_parity` |
> | unstructured × semi-implicit | ✅ | ❌ **gap** |
>
> Both new fixtures use a sloping bed, an active density field driving genuine flow (measured
> max|hu| ≈ 3.1e-2 and max|xmom| ≈ 1.1e-2), and carry guards asserting the flow actually develops so
> they cannot decay into vacuous quiescent checks. Both conserve mass at ~1e-15.
>
> **The wall bed-mirror fix is structured-only by nature** — the unstructured KNP path already
> mirrors the bed at wall edges (`zr = zl`, `hre = hle` in `compute_boundary_state_impl`), so it
> never had that defect. Structured was the odd one out.
>
> **The remaining gap** is unstructured × semi-implicit. That path is TPFA and places η at triangle
> circumcentres, so it needs a Delaunay-quality mesh; the right-triangle split used by every
> quad-derived test mesh here gives coincident circumcentres and trips the mesh quality gate. It
> needs the per-quad diagonal-choice construction from `test_si_unstr_ml_wetdry`
> (`delaunay_diag_bl_tr`) plus `state_unstr_enable_si` — **and that routine's `ok` flag must be
> asserted**, since the existing test does `if (.not. ok) return` and would go silently vacuous on a
> mesh that failed the gate. Worth auditing whether it currently does.
>
> Separately, `test_multilayer_unstr`'s `tracer_mass_conserved_in_closed_basin` is left alone but
> should not be mistaken for coverage: four triangles, flat bed, 25 steps,
> `beta_S = alpha_T = 0` (no baroclinic coupling at all), and it asserts on the tracer only, never
> on water mass.


*Original plan text follows.*

- **What it validates.** A *closed budget* on the coastal path: `residual = (∫∫∫ h·C dV)_final −
  (…)_initial − ∫ (sources) dt`, per conserved quantity (water mass, salt, heat). Today only a per-step
  *drift print* exists; there is no source-accounted closure a test can assert to round-off. This is the
  reviewer's "can I trust a multi-week salt budget?" check (gaps §4 flags it explicitly absent).
- **Reference.** Exact — it is a discrete conservation identity, not an approximation. In a closed domain
  with no forcing the residual is machine zero; with a prescribed surface flux Q it equals ∫Q dt.
- **Setup.** New module `rki_coastal_budget` (or extend `rki_init`): reuse the existing GPU-reduction
  `compute_total_mass` / `compute_layer_sums` (device-resident, only the scalar crosses PCIe). Accumulate
  boundary tracer flux (open/tidal), surface flux (heat, E–P), precip, and river sources into running
  totals. Namelist knob `&diagnostics_nml conservation_budget=.false.` default off. No new physics.
- **Backends/regime.** Both backends; 2D + 3D multilayer. NH and semi-implicit inherit it (they share the
  tracer state).
- **Metric & tolerance.** (a) Closed box (walls, zero forcing, non-trivial IC): |residual|/|total| <
  **1e-12** for mass, **1e-10** for salt/heat — justified because the flux-form update is conservative by
  construction; anything above round-off is a real leak (this is exactly how the wet/dry leak of gaps §3
  would surface). (b) Prescribed constant surface heat flux, insulated sides: closure to ∫Q dt within
  **1e-8** relative. These are *identity* tolerances — not chosen to pass, they are what conservation
  means.
- **Where it lives.** `tests/test_coastal_conservation.F90` (analytical, closed-form zero/∫Q oracle).
  Mirror the ocean twin `test_ocean_conservation_salt_heat.F90` for the `mem:separate` device-mapping
  discipline (map every state + scratch companion; the ocean test's header documents the contract).
- **Size / infra.** Medium. New reusable budget module — the foundation Phases 5, 6, 8, 9 assert against.

---

### Phase 2 — Convergence-order harness (INFRASTRUCTURE) — **DONE 2026-08-07**

**Delivered** as `tests/convergence_harness.F90` (added to the `rakali_analytical` static lib) plus
`tests/test_convergence_harness.F90` (9 self-check assertions). Verified on gcc 15.1.0 CPU and
NVHPC 26.5 GPU (cc70/V100): 121/121 `ctest -R rakali` green on both.

**API** — `error_norms_uniform` / `error_norms_weighted` (L1, L2, Linf as integrals; the weighted
form takes per-cell areas for the unstructured backend), `convergence_rates` (pairwise + overall
observed order), and `spacing_from_counts`.

**Two deviations from the plan as written, both deliberate:**

1. **No procedure-pointer callback.** The plan called for "a problem callback + resolution list".
   In Fortran that buys ceremony rather than reuse: each consumer's setup differs enough (grid, IC,
   BC, run length, device data dance) that any callback signature would be over-fitted to its first
   caller. The reusable, independently testable content is the *rate arithmetic and the norms* —
   the consumer drives its own resolutions and hands over `(h, err)`. This also makes the harness
   fully unit-testable with no solver in the loop, which is what lets the self-check calibrate it.
2. **`convergence_rates` returns an `ok` flag** and refuses degenerate input — fewer than two
   levels, a non-positive or non-finite error, or a non-strictly-decreasing `h`. An exactly-zero
   error means the scheme reproduced the oracle, so the order is *undefined*, not infinite.
   Silently returning a plausible number from degenerate input is how a convergence test starts
   certifying nothing.

**`spacing_from_counts` exists for a reason:** the two hand-rolled rate computations already in the
suite disagree on convention — `test_convergence.F90` divides by a ratio of `dx`, while
`test_dam_break_unstr.F90` divides by a ratio of `nx`. Both are correct, which is precisely why a
third consumer would have a 50 % chance of getting it wrong. The helper makes the coarse-first
spacing convention explicit at the call site.

**The self-check includes a plateau case** — a clean 2nd-order series whose finest level flattens
against an error floor. Its *overall* rate stays above 1.8 and would pass a naive assertion, while
its final *pairwise* rate collapses to ~0. That is the concrete reason both numbers are returned,
and Phase 7 and any later refinement study should assert the overall rate **and** inspect the
pairwise series.

*Original plan text follows.*

#### Original spec

- **What it validates.** Grid-refinement observed order `p = log(E_c/E_f)/log(r)` from a sequence of
  runs. Only one refinement test exists (unstructured KNP dam-break, `test_convergence.F90` is structured
  dam-break but shock-limited). No *smooth-solution* order verification exists on either backend.
- **Reference.** The scheme's design order: structured HLL + minmod-limited MUSCL + SSP-RK2 → formally
  2nd order in smooth monotone regions, dropping to 1st at shocks and at smooth extrema (minmod clips
  them — see memory `coastal-time-schemes`). Unstructured KNP + Green–Gauss → 2nd order smooth.
- **Setup.** A reusable Fortran harness `convergence_harness` (add to `rakali_analytical`): takes a
  problem callback + resolution list `[N, 2N, 4N, 8N]`, returns L1/L2 errors and pairwise + overall
  rates. No namelist; test-only infrastructure.
- **Metric & tolerance.** Overall (coarse→fine) observed rate ≥ **1.8** on a *smooth, extremum-free*
  manufactured/analytic solution (e.g. a smoothly-varying steady state or a low-amplitude standing wave
  with no interior extremum). Per-pair floor ≥ 1.5. Justification: the interior scheme is 2nd order; the
  0.2 slack absorbs boundary/limiter contamination. Do **not** use a case with a smooth crest — minmod
  caps that at ~1 and the number would be defensible-but-misleading (this is why `test_convergence`'s
  shock case only asserts ≥ 0.3).
- **Where it lives.** Harness in `rakali_analytical`; first consumer is Phase 7. A trivial self-check
  test `tests/test_convergence_harness.F90` (feeds it an analytic O(h²) error sequence, asserts p≈2).
- **Size / infra.** Medium. Pure test infra, no source change.

---

### Phase 3 — Well-balanced C-property (lake at rest over irregular bathymetry) — **DONE 2026-08-07**

**Delivered** as `tests/test_c_property.F90` (8 assertions). Verified on gcc 15.1.0 CPU and
NVHPC 26.5 GPU (`cc70`, V100): 117/117 `ctest -R rakali` green on both.

**Corrections to the plan as written, found during implementation:**

1. **"never asserted standalone" was wrong.** `tests/test_lake_at_rest.F90` already asserted the
   C-property. It was weak, not absent — smooth fully-submerged parabola, 100 steps, 1e-10
   absolute, structured 2D only. The new test strengthens it on five axes (irregular bed,
   *emergent* island so the wet/dry front sits on topography, 1000 steps, round-off tolerances,
   both backends + 3D multilayer). `test_lake_at_rest` is kept as the cheap smoke version.
2. **The metric is momentum, not velocity.** `u = hu/h` is ill-conditioned at the wet/dry edge —
   a cell at `h ~ 1e-7` turns an acceptable `hu ~ 1e-15` into `u ~ 1e-8`, reporting a violation
   whose size depends only on proximity to the dry threshold. Momentum is what the flux balance
   conserves, so the primary assertion is `|hu|,|hv| < 1e-12 * H*sqrt(gH)`. The plan's literal
   velocity metric is kept as a second assertion, depth-gated to cells deeper than `0.1H`.
3. **Sizing `dt` on the unstructured backend must come from `mesh%radius`,** the min
   centroid-to-edge-midpoint distance that `compute_timestep_unstr` actually divides by — not the
   quad edge. For the right-triangle split the radius is ~3.4x smaller; sizing off the quad
   overshoots stability and the run overflows.
4. **A fixture guard is included** (`c_property_structured_bathy_is_emergent`) asserting the bed
   really pierces the surface and really leaves deep water. Without it, a later edit to the
   bathymetry constants could silently submerge the island and the C-property tests would keep
   passing while no longer testing the wet/dry balance they exist for.

**Validated by negative control:** a 1e-9 surface tilt (100x `ETA_TOL`) fails all three structured
assertions, confirming the tolerances are live rather than vacuous.

*Original plan text follows.*

- **What it validates.** The "C-property" (Bermúdez & Vázquez 1994): a fluid at rest over arbitrary
  bathymetry must stay *exactly* at rest — the discrete pressure/topography flux balance holds to
  round-off. Verifies Audusse hydrostatic reconstruction (`rki_kernel_flux.F90`) as a standalone
  assertion, not a byproduct.
- **Reference.** Bermúdez & Vázquez (1994); Audusse et al. (2004). Exact: `η = const`, `u = v = 0`
  forever ⇒ oracle is the initial state.
- **Setup.** Domain 100×100 (structured) / equivalent triangulation (unstructured), irregular bathymetry
  including a partially-emergent bump (tests wet/dry edge balance too). IC `η = η₀` flat, `u=v=0`. Walls.
  `dt` at CFL 0.45. Run 1000 steps. Both backends; 2D barotropic AND 3D multilayer (a flat-ρ column must
  also stay at rest — this catches the formula-bathymetry ghost-fill bug class, memory `formula-bathy-ghost`).
- **Metric & tolerance.** `max|u|, max|v| < 1e-12·√(gH)` and `max|η − η₀| < 1e-11·H` after 1000 steps.
  Justification: this is an algebraic identity of a well-balanced scheme — the only defensible tolerance
  is machine epsilon accumulated over the step count. A failure is a real well-balancedness bug, never a
  tolerance to relax.
- **Where it lives.** `tests/test_c_property.F90` (analytical; oracle = IC, no external file).
- **Size.** Small (~half day). No new infra (may use Phase 1 diag as a secondary check but not required).

---

### Phase 4 — Dry-bed dam break (Ritter dry front / Stoker) — **DONE**

**Delivered** as `tests/test_dam_break_dry.F90` (7 assertions) plus the dry-bed branch
(`ritter_dry_sample`, `ritter_dry_front_position`) added to `tests/ritter_dambreak.F90`. Verified on
gcc 15.1.0 CPU and NVHPC 26.5 GPU: 124/124 `ctest -R rakali` green on both.

Passing: front position (both backends), L1 depth error (both), positivity, and mass conservation —
the first conservation check anywhere on the **wet/dry** path, since every Phase 1 fixture is fully wet.

> ### ✅ FINDING (RESOLVED 2026-08, see update below) — the sub-1 %-depth leading tongue was not represented
>
> Front position vs the consistently-thresholded exact profile, structured, 800 cells, t = 10 s,
> 62.6 m of travel:
>
> | thresh/h_l | x_num | x_ref | lag (m) | rel |
> |---|---|---|---|---|
> | 1e-1 | 132.88 | 132.92 | **0.05** | **0.0008** |
> | 1e-2 | 146.13 | 153.24 | 7.11 | 0.114 |
> | 1e-3 | 146.38 | 159.66 | 13.29 | 0.212 |
> | 1e-6 | 146.63 | 162.54 | 15.91 | 0.254 |
>
> At a resolvable depth the front is essentially exact (0.08 %). But `x_num` barely moves from
> threshold 1e-2 down to 1e-6 — the model has a near-vertical cliff at x ≈ 146.5 while the exact
> solution tapers parabolically to 162.6. The thin tongue from h ≈ 0.03 m down to zero is missing,
> so the **tip advances at ~1.49 c₀ instead of 2 c₀**.
>
> Mass is conserved to round-off, and the missing tongue is only ~0.2 % of the water — which is why
> the L1 assertion passes comfortably. **No depth-weighted norm can see this**; only a
> front-position metric does. That is the reusable lesson for P5 and P6, which also have moving
> shorelines.
>
> **Cause not identified — three candidates eliminated by experiment.** `THIN_LAYER_THRESHOLD`
> (zeroes momentum below 1 mm): rebuilt at 1e-9, tip moved 0.25 m of a 15.9 m deficit. `FROUDE_CAP`
> (the exact Froude number *diverges* at the tip, `u → 2c₀` while `c → 0`, so a cap looked
> compelling): rebuilt at 1e8, no change. Bottom friction: `manning_n` defaults to 0. Both
> momentum-clipping mechanisms in `compute_update` are therefore excluded, which points **upstream
> to the flux** — the positivity limiter and hydrostatic reconstruction at a wet/dry face. 15.9 m at
> dx = 0.25 m is ~64 cells of lag, far too structural to be ordinary 1st-order tip error.
>
> **Not hidden behind a loosened tolerance.** Front position is asserted at a depth the scheme can
> represent, and the tip deficit is asserted separately and two-sidedly: a lower bound guards against
> regression, and an **upper** bound fails if someone fixes it — the prompt to retire the guard and
> tighten onto the true front. A one-sided assertion would let a fix land silently and leave the
> file permanently understating the model.

> ### ✅ UPDATE (2026-08) — fixed: the two "eliminated" clamps were REDUNDANT
>
> The single-knob eliminations above were both correct measurements and the wrong conclusion: the
> two clamps are **redundant** — relax either alone and the other keeps binding, which is exactly
> why each looked innocent. Dose-response (tip celerity / c₀, only `rki_constants.F90` varied):
> `THIN=1e-3, CAP=10 → 1.489` (baseline); `1e-9, 10 → 1.489`; `1e-3, 1e8 → 1.577`; `1e-9, 50 →
> 1.720`; `1e-9, 1e8 → 1.888` — but the fully-relaxed configuration also has a genuine thin-film
> velocity runaway and CFL dt collapse on beaches, so the clamps were doing necessary work with the
> wrong shape.
>
> **Fix** (`fix/thin-cell-desingularization`): both blunt clamps replaced by the Kurganov & Petrova
> (2007) desingularized velocity `u = √2·h·hu/√(h⁴+max(h⁴,eps⁴))`, `hu := h·u`, with
> `eps = DESING_EPS = 1.5e-4 m` (the codebase's `H_VANISHED` scale) — shared with the CFL dt
> estimator — plus a `FROUDE_BACKSTOP = 50` blowup arrester (the measured runaway reached
> u ~ 4×10³ m/s at h ≈ 3.5 mm with a ×2000 dt collapse; the backstop bounds it at ~10 m/s with dt
> within ~4× of the wet value). Applied to `compute_update`/`compute_rk2_average`/
> `compute_ssp_combine`, both fused kernel files, and the unstructured protect/fused kernels.
> Measured after the fix (structured): tip **1.489 → 1.736 c₀** (no-clamp ceiling 1.888 at this
> resolution); front-position lag at the 1-cm contour **7.11 m → 0.36 m** (0.6 % of travel).
> `H_FRONT_STRUCT` tightened 1e-1 → **1e-2·h_l** at the same 1.5 % tolerance; the tip guard
> re-centred two-sidedly on (1.65, 1.85).
>
> **Residuals, stated honestly:** (1) the tip still lags 2c₀ — the sub-`DESING_EPS` band is damped
> by design (that damping is what keeps beaches stable) and 1st-order smearing does the rest;
> (2) the **unstructured** tip (1.586 c₀) is *unchanged* by the same protect-kernel treatment — its
> binding constraint is the KNP reconstruction's dry-cell edge-momentum zeroing
> (`centroid_pass`, `minimum_allowed_height`), now pinned two-sidedly at (1.45, 1.75) by
> `dry_tip_celerity_unstr_is_known_deficient`, with `H_FRONT_UNSTR` kept at 1e-1·h_l;
> (3) wet/dry-heavy runs take more, smaller steps (Thacker bowl: 849 → 2612 steps/period) because
> thin fast water now exists and must be CFL-resolved — the old cheapness was the freeze;
> (4) Thacker planar L2 rose 0.078 → 0.089 (retention rim + splash shift the interior slosh ~1 cm
> RMS; tolerance re-baselined 0.08 → 0.10 with the measured history recorded in the test).

*Original plan text follows.*

#### Original spec

- **What it validates.** The dry leading-edge rarefaction. Both existing dam-break tests are *wet-bed*;
  the dry front is only NaN-smoke-tested (gaps §5). This checks front position and tip celerity against
  the exact dry-bed solution.
- **Reference.** Ritter (1892) dry-bed exact: front position `x_f(t) = 2·c₀·t`, tip celerity `2c₀` with
  `c₀=√(g h₀)`; parabolic depth profile behind the front. Extend the existing `ritter_dambreak.F90`
  oracle (currently wet-bed Stoker) with the dry-bed branch `h_R=0`.
- **Setup.** Channel L=200 m, `h_left=h₀`, `h_right=0` (dry). Structured 800 cells; unstructured
  comparable. `dt` CFL 0.45, `T_end` before the front hits the wall. 2D barotropic. Both backends (dry
  front stresses the positivity/wet-dry path on each). NH: not meaningful (hydrostatic bore). Semi-
  implicit: skip — Casulli ELM smears the sharp dry front and is not the intended tool (state why).
- **Metric & tolerance.** Front-position error < **1.5 %** (kinematic, well captured); L1 depth error <
  **5 %**. Justification: HLL is 1st-order at the vacuum-adjacent rarefaction tip (the wet-bed test hits
  2.5 % L1; the dry front is strictly harder), so 5 % L1 is the order-honest bar — reachable, and a
  refinement rate check (via Phase 2 harness) confirms it is convergence, not luck. Also assert mass
  conservation < 1e-10 (Phase 1 diag).
- **Where it lives.** `tests/test_dam_break_dry.F90` (+ dry branch in `ritter_dambreak.F90`).
- **Size.** Small–Medium. Reuses the Ritter oracle and dam-break scaffolding.

---

### Phase 5 — Carrier–Greenspan periodic runup — **PARTIAL**

**Delivered** as `tests/carrier_greenspan.F90` (oracle, added to `rakali_analytical`) and
`tests/test_carrier_greenspan.F90` (5 assertions). 125/125 `ctest -R rakali` green on gcc 15.1.0
CPU and NVHPC 26.5 GPU.

**Oracle.** The linear Carrier–Greenspan solution `eta = A0 J0(xi) cos(wt)`,
`xi = 2w sqrt(X/(g s))`. C&G's central result is that the hodograph transform maps the *nonlinear*
runup problem onto this same linear equation, so it predicts runup correctly for finite amplitude
up to breaking — this is the standard result, not a small-amplitude approximation. Bessel functions
come from the F2008 intrinsics, which gfortran 15.1 and nvfortran 26.5 evaluate identically
(`J0(3.8317059702) = -0.402759396` on both), so no series and no cross-toolchain drift.

**Fixture — a self-forcing eigenmode.** The velocity goes as `-J1(xi)`, so an offshore wall placed
at a **zero of J1** is an exact velocity node and the J0 mode becomes an eigenmode of a closed
basin: seed it and it oscillates, with no boundary forcing. That avoids entangling a runup
measurement with the tidal-BC implementation, and makes total mass exactly invariant so the wet/dry
cycle becomes a conservation assertion for free.

**✅ Validated.** Oscillation at the Carrier–Greenspan frequency to within 5 %, and mass conserved
to 1e-8 through repeated wet/dry cycling — the only fixture in the suite that wets and dries the
same cells over and over, so that is a new guarantee rather than a restatement of Phase 1.

> ### 🟡 FINDING (UPDATED 2026-08) — the momentum-clamp freeze is fixed; a mm-scale
> retention film remains
>
> Original finding: locating the shoreline at 1e-3 m depth gave an excursion of **exactly zero** —
> water thinner than ~1e-2 m never receded, frozen by the same two redundant momentum clamps as the
> Phase 4 tongue (see the Phase 4 update). After the thin-cell desingularization fix the beach
> **drains down to the hydrostatic-reconstruction retention scale**: Audusse-type well-balancing on
> a plane beach retains up to `slope·dx` (= 12.5 mm/cell here; static drain probe measures ~5 mm
> mean, film volume halves in 8 s where before it was frozen solid). The freed runup jet also
> overwets the upper beach with a mm film (to x ≈ 128.4 vs the exact 120 m maximum runup), so the
> 1-mm detector is still pinned — now at the overshoot mark, by **grid-dependent discretization
> retention**, no longer by a clamp. Still guarded two-sidedly
> (`cg_thin_film_persists_on_beach`); fixing it needs a bed-slope-aware thin-film treatment in the
> reconstruction, not another knob in the update.

> ### 🟡 REMAINING (UPDATED 2026-08) — runup amplitude now guarded at 10 %; measured 5.6 %
>
> Pre-fix the excursion scattered ~15 % both ways (18.9/22.1/22.1 by detector, 18.1/21.5/16.5 by
> period) because the frozen film moved the detector, and the guard sat at a coarse 40 %. Post-fix
> the three detection depths agree to ~2 % — **18.62 / 19.00 / 18.88 m at 3e-1 / 3e-2 / 1e-2 m**
> against the exact 20 m — i.e. 5.6 % at the 1-cm detector, with the rundown edge at x = 83.9 m vs
> the exact 80 m, and the period within 1.5 %. The guard is tightened **40 % → 10 %**.
>
> The headline **5 % is still not claimed** (task #9 not fully closed): the shoreline is a
> cell-index staircase (one cell = 0.6 % of the excursion; toolchains can differ by a cell), and
> the `~slope·dx` retention film biases the detector shoreward. Closing the last ~0.6 % honestly
> needs finer dx or a retention-aware detector — not a tolerance chosen to match one measurement.
> The backstop value matters here and was chosen by this measurement: `FROUDE_BACKSTOP = 100` let
> the freed runup jet strand a cm-scale wedge through rundown (excursion 11.8 m); 50 restores it
> (18.9 m). See `FROUDE_BACKSTOP` in `rki_constants.F90`.

**Lesson recorded for later phases:** the period assertion originally read the *shoreline* track and
passed on gcc while failing on nvfortran — 3 crossings vs 2, a factor-of-two period difference
purely from how each toolchain rounded a cell-index staircase whose rundown extent is itself
defect-affected. It now reads the offshore free-surface elevation, which sits in 5 m of water and
never wets or dries. **Do not build a periodicity measurement on a wet/dry-affected signal.**

*Original plan text follows.*

#### Original spec

- **What it validates.** Nonlinear periodic wave runup on a plane sloping beach — the canonical moving-
  shoreline analytical case. Today Thacker (parabolic bowl) is the *only* moving-shoreline analytic test
  and its planar solution is loose (8 %).
- **Reference.** Carrier & Greenspan (1958): hodograph-transform exact solution for periodic runup on a
  uniform slope; shoreline position and velocity are closed-form (Bessel-function series). Encode the
  series as a Fortran oracle in `rakali_analytical`.
- **Setup.** 1-D-in-x plane beach (slope ~1:10–1:30), narrow in y, structured; sinusoidal tidal forcing
  at the offshore open boundary at the C–G frequency; run several periods to periodic steady state.
  Structured 2D barotropic is primary. Unstructured: include if the sloping-beach triangulation + open
  BC compose cleanly; otherwise defer and say why (C–G is a 1-D beach — the structured grid is the
  natural fit; the value on triangles is marginal). NH/SI: skip (hydrostatic long-wave regime).
- **Metric & tolerance.** Runup amplitude R (max shoreline excursion) within **5 %** of C–G; shoreline
  phase within **5 %** of the forcing period. Justification: the wetting-front reconstruction is 1st
  order, so 5 % on the excursion is order-consistent; tighter would require a 2nd-order shoreline tracker
  the model does not have. Mass conserved (Phase 1) through the wet/dry cycle < 1e-8 (this also directly
  exercises the gaps §3 wet/dry leak on a *clean* geometry).
- **Where it lives.** `tests/test_carrier_greenspan.F90` (+ C–G series oracle in `rakali_analytical`).
- **Size.** Medium (the oracle series is the bulk of the work).

---

### Phase 6 — Flow over a bump (steady SWE benchmark, three regimes)

- **What it validates.** Steady flow over a smooth bump — the classic SWE verification with subcritical,
  transcritical-with-shock, and supercritical branches; tests well-balancedness *with* motion, shock
  capture, and discharge conservation.
- **Reference.** Steady analytic solution from Bernoulli + constant specific discharge `q = h·u`:
  `q²/(2 g h²) + h + b(x) = const`, with a Rankine–Hugoniot jump in the transcritical case. Classic
  Delft/SWASHES benchmark set (Goutal & Maurel 1997; Delestre et al. 2013 SWASHES). Oracle = root-find on
  the cubic per column (add to `rakali_analytical`).
- **Setup.** Channel L=25 m, bump `b(x)=0.2−0.05(x−10)²` for 8<x<12. Three cases via inflow discharge +
  outflow depth BC: subcritical (`q=4.42`, `h_out=2`), transcritical-with-shock (`q=0.18`, `h_out=0.33`),
  supercritical. Structured, ~200 cells, run to steady state. Both backends. Regimes: **2D barotropic**
  all three cases; **semi-implicit** the *subcritical* case only (no shock — ELM is appropriate; state
  that transcritical needs the explicit Godunov path); **NH** the subcritical case as a consistency check
  (NH correction ≈ 0 for this hydrostatic flow — verifies NH doesn't corrupt a hydrostatic steady state).
- **Metric & tolerance.** Free-surface L2 vs analytic < **1 %** (subcritical/supercritical, smooth), <
  **3 %** (transcritical, 1st-order at the shock); discharge `q` constant along the channel to < **1 %**
  (this is the well-balanced-with-flow check). Justification tied to scheme order at the shock exactly as
  Phase 4.
- **Where it lives.** `tests/test_flow_over_bump.F90` (+ steady-cubic oracle in `rakali_analytical`).
- **Size.** Medium. Three cases, one oracle.

---

### Phase 7 — Structured HLL convergence order — **DONE**

**Delivered** as `tests/test_convergence_structured.F90` (4 assertions), consuming the Phase 2
harness. Verified on gcc 15.1.0 CPU and NVHPC 26.5 GPU (cc70/V100): 123/123 `ctest -R rakali` green
on both.

**Result — the structured HLL path reaches its design order.**

| norm | pairwise rates | overall |
|---|---|---|
| L1 | 1.936 · 1.966 · 1.981 | **1.961** |
| L2 | 1.953 · 1.974 · 1.985 | **1.971** |
| Linf | 1.942 · 1.974 · 1.984 | **1.967** |

Two details matter more than the headline. The pairwise rates **increase** monotonically toward 2,
which is asymptotic convergence rather than a lucky average — the opposite of the plateau signature
the harness self-check pins. And **Linf converges at 1.97 as well**, not at 1 as it would if minmod
were clipping a smooth extremum; that is the direct evidence the case selection worked.

**Case selection was the whole difficulty.** The plan warned that a smooth crest caps the rate near
1 and yields a "defensible-but-misleading" number. The fundamental standing mode
`eta = A cos(pi x/L) cos(wt)` avoids it structurally rather than tolerating it: `cos(pi x/L)` is
**monotone** on `[0, L]`, so its extrema sit exactly ON the walls and the interior contains no
extremum for minmod to clip — while `d(eta)/dx = 0` at both walls also makes the mode compatible
with the wall BC, so no boundary layer contaminates the interior.

**Two constraints that would silently corrupt the measurement, both handled:**

1. **Amplitude.** The oracle is the *linearised* solution, so the model's own nonlinearity is an
   error source of relative size `O(A/H0)`. Measuring discretisation error requires
   `A/H0 << C h^2`; at `A/H0 = 1e-7` there is ~2 orders of headroom at the finest level. A larger
   amplitude flattens the series against the nonlinear floor.
2. **Timestep.** SSP-RK2 is 2nd order in time, so `dt` must shrink with `h` or the measured order
   blends space and time. Each level takes exactly 2x the steps of the previous one to the *same*
   final time, making `dt ~ h` by construction and landing every level on an identical `t`.

**The 1.8 floor was fixed by this plan before the numbers were measured** and is left unchanged —
it is not a threshold fitted to the result.

*Original plan text follows.*

#### Original spec

- **What it validates.** The missing structured refinement test (gaps §5: "no refinement test at all,
  only the unstructured KNP one"). Confirms the structured HLL+MUSCL+SSP-RK2 path achieves design order
  on a smooth solution.
- **Reference.** Formal 2nd order (see Phase 2). Use a *smooth extremum-free* analytic solution — a
  low-amplitude linear standing wave in a closed basin (analytic `η=a cos(kx)cos(ωt)`, `ω²=gk²H` in the
  linear limit) run at small amplitude so nonlinearity and limiter clipping are negligible, or the smooth
  region of the dry-bed profile away from front/extremum.
- **Setup.** Closed basin, `[N,2N,4N,8N] = [50,100,200,400]`, fixed small `dt` (or `dt∝h` to hold
  temporal order constant), a fraction of a period so the wave stays smooth. Structured 2D barotropic.
  Unstructured already has its convergence test — note that and don't duplicate; optionally extend the
  same smooth case to KNP for parity.
- **Metric & tolerance.** Observed overall L2 rate ≥ **1.8** (per Phase 2 justification). If the measured
  rate is ~1 instead, that is the minmod extremum-clipping signature — the fix is to pick a monotone
  smooth profile, NOT to lower the bar.
- **Where it lives.** `tests/test_convergence_structured.F90`, consuming the Phase 2 harness.
- **Size.** Small once Phase 2 exists.

---

### Phase 8 — Lock-exchange tightening (±50 % → defensible band)

- **What it validates.** The stratified gravity-current front speed — the flagship baroclinic coupling
  test, currently asserted only at ±50 % vs Benjamin's *inviscid energy-conserving* 0.5√(g'H).
- **Reference.** The right reference is the **dissipative** half-depth front, not the inviscid 0.5:
  Benjamin (1968) upper bound `0.5√(g'H)`; the physically-realised (and lab/DNS) energy-*consuming* front
  is `≈0.44–0.47 √(g'H)` (Shin, Dalziel & Linden 2004; Härtel, Meiburg & Necker 2000). Assert against the
  dissipative value with a justified band, not the inviscid ceiling.
- **Setup.** As `test_lock_exchange.F90` but refine: raise `NZ` 7→≥12, use the least-diffusive tracer
  path, and set `dt` for CFL ~0.3. 3D multilayer, both backends. Keep the S-only and T-only variants
  (density-equivalent) and the symmetry check.
- **Metric & tolerance.** Front speed within **±20 %** of `0.45√(g'H)` (was ±50 % of `0.5√(g'H)`).
  Justification and honest ceiling: the lock-exchange front head is **genuinely non-hydrostatic**; a
  hydrostatic sigma model + HLL/CWC tracer diffusion cannot reach the ±5 % of a NH/DNS code. ±20 % is
  what vertical refinement + reduced numerical diffusion buys on the hydrostatic path; closing further
  would require the NH-on-stratified path (gaps §11 — NH is currently homogeneous-density, so this is out
  of scope) and is documented as the residual. Add a salt/heat closure assertion (Phase 1) so the
  tightened front is not bought by a conservation regression.
- **Where it lives.** Tighten `tests/test_lock_exchange.F90` in place; document the dissipative-reference
  change in the docstring (currently it cites Benjamin 0.5 — that citation becomes the *bound*, 0.45 the
  *target*).
- **Size.** Medium (mostly re-tuning + refinement sweep + convergence confirmation via Phase 2).

---

### Phase 9 — Kato–Phillips in-repo oracle + tightening

- **What it validates.** Wind-driven entrainment deepening — the flagship mixed-layer test, currently a
  1.6–2.4 ratio band whose clean √t law lives OFF-REPO in Python (gaps §5; test docstring lines 11–14).
  Bring the oracle in-repo and tighten to an exponent + coefficient check.
- **Reference.** Kato & Phillips (1969): mixed-layer deepening `D(t) = c·u_*·√(t/N₀)` under constant
  wind stress into a linear stratification, coefficient `c ≈ 1.05` (Price 1979 re-analysis). Encode the
  closed-form law as a Fortran oracle in `rakali_analytical` (replaces the off-repo Python).
- **Setup.** As `test_kato_phillips.F90`: 1-D water column, linear initial stratification `N₀`, constant
  surface wind stress. **Critically** set the wind to the surface layer, not the depth-uniform barotropic
  mode — `&… wind_surface_layer=.true.` (memory `coastal-wind-bt`; the test docstring names this exact
  blocker to a clean √t). 3D multilayer, `NZ≥20` for a resolvable ML base. Structured primary;
  unstructured has its own `test_keps_kato_unstr.F90` — extend there for parity if the surface-wind knob
  is wired on that backend, else note the limitation.
- **Metric & tolerance.** (a) log–log slope of `D(t)` = **0.5 ± 0.05** (the √t exponent — the real
  physics content, was never asserted). (b) coefficient `c` within **±15 %** of 1.05. Justification: the
  exponent is dimensionally exact and should be tight; ±15 % on `c` absorbs the KPP/k-ε closure constant
  spread and finite vertical resolution. This replaces the order-of-magnitude 1.6–2.4 *ratio* band with a
  law-shape + coefficient assertion. If the exponent is off, the wind is still reaching the barotropic
  mode — fix the forcing routing, don't widen the band. Assert heat/tracer closure (Phase 1).
- **Where it lives.** Tighten `tests/test_kato_phillips.F90`; add the K–P law oracle to
  `rakali_analytical` and delete the off-repo Python dependency from the docstring.
- **Size.** Medium.

---

## 4. Tolerance philosophy

**A tolerance is a prediction about discretisation error, not a knob to make CI green.** Every number in
§3 is derived from one of three sources, in priority order:

1. **Algebraic identities → machine epsilon.** Well-balancedness (Phase 3), conservation closure
   (Phase 1), the √t *exponent* (Phase 9). These are properties the scheme either has exactly or does
   not; the only defensible tolerance is round-off accumulated over the step count. A failure is a bug in
   the property, full stop.
2. **Scheme order at the hardest feature → order-consistent bound.** Dry front (Phase 4, 1st-order at the
   vacuum tip → 5 % L1), transcritical shock (Phase 6, 1st-order at the jump → 3 %), runup front
   (Phase 5, 1st-order shoreline → 5 %). The tolerance is set by *where the scheme drops to 1st order*,
   and each is paired with a **refinement check** (Phase 2 harness) proving the error is *converging* —
   an absolute number without a rate is not evidence.
3. **Known model-physics ceiling → documented residual.** Lock exchange (Phase 8): a hydrostatic model
   *cannot* reach a NH front speed; the band (±20 %) is what the hydrostatic path can achieve and the gap
   to lab is attributed, in the docstring, to the missing non-hydrostatic head — with the *reference
   itself* corrected to the dissipative value so we are not measuring against an unreachable ideal.

**When a case fails, the response is never "loosen it."** In priority: (a) is it a real bug? — a failing
C-property or conservation closure is a defect to fix in the kernel, not the test; (b) is the tolerance
*mis-derived*? — re-derive from the scheme order and the harness-measured rate, which may move it up or
down, but only with a written justification; (c) is it a genuine model ceiling? — then the fix is to
correct the *reference* (dissipative vs inviscid front) and document the residual physics, not to inflate
the band around the wrong reference. A number changed without one of these three written reasons is a
regression in the campaign's credibility, which is the one thing it exists to protect.

---

## 5. What this campaign does NOT cover

- **Missing physics (gaps §2), not numerics.** No wave coupling (radiation stress / Stokes / wave setup),
  no morphodynamics / Exner, no wave-current bottom stress, no bulk air–sea flux, no tidal potential/SAL,
  no energy-extraction sink, no beta-plane. Those are *capability* PRs; you cannot validate physics that
  is not there. This campaign validates the numerics that *are* present.
- **Non-hydrostatic stratified benchmarks** (internal solitary waves, NH plume fronts) — the NH path is
  homogeneous-density (gaps §11), so a stratified-NH validation has no code to run against. Lock exchange
  (Phase 8) documents this as the residual but does not attempt to close it.
- **Observational / field cases** — no Sydney Harbour tide-gauge comparison, no storm-surge hindcast
  (e.g. named-event surge), no GRS-type field validation. Those need forcing data and a data-ingest path
  (gap §5 forcing breadth) and are a separate, later validation tier against *measurements* rather than
  analytical solutions.
- **Unstructured coverage is opportunistic.** Where a benchmark is intrinsically 1-D-on-a-slope
  (Carrier–Greenspan) or a scheme choice makes it marginal, the plan says so and stays structured rather
  than forcing a low-value triangulation.
- **Semi-implicit shock cases** — Casulli ELM is deliberately not the tool for dry fronts / transcritical
  shocks (it smears them); those regimes are validated on the explicit Godunov path only, by design.
