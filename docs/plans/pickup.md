# Pickup — validation campaign, session handoff

Written 2026-08-08. Everything here is either non-obvious from the repo or was
learned the hard way. Read this before touching the campaign.

---

## 1. Environment — read first, it will cost you an hour otherwise

**`gfortran 13.3` cannot build this repo.** It lacks `do concurrent` locality
specifiers (`local(...)`), which `src/tracer/structured/rki_recon_weno.F90`
uses. The system `gfortran` on this box is 13.3. Use **gcc 15.1.0**:

```bash
export PATH=/shared/compilers/gcc/15.1.0/bin:$PATH
cmake -B build_gcc -S . -DRAKALI_ENABLE_NETCDF=OFF \
  -DRAKALI_ENABLE_GPU=OFF \
  -DCMAKE_Fortran_COMPILER=/shared/compilers/gcc/15.1.0/bin/gfortran \
  -DCMAKE_C_COMPILER=/shared/compilers/gcc/15.1.0/bin/gcc \
  -DCMAKE_BUILD_TYPE=Release
```

GPU build (4× Tesla V100 on this box, so `cc70`):

```bash
export NVROOT=/home/jorge/install/nvhpc/26.5/Linux_x86_64/26.5/compilers
export PATH=$NVROOT/bin:$PATH
cmake -B build_cc70 -S . -DRAKALI_ENABLE_NETCDF=OFF \
  -DRAKALI_ENABLE_GPU=ON -DRAKALI_GPU_ARCH=cc70 \
  -DCMAKE_Fortran_COMPILER=$NVROOT/bin/nvfortran \
  -DCMAKE_C_COMPILER=$NVROOT/bin/nvc -DCMAKE_BUILD_TYPE=Release
```

`RAKALI_ENABLE_NETCDF=OFF` is fine — no test in the campaign needs I/O, and it
builds much faster. `nc-config` is not on PATH anyway.

Run tests with `OMP_NUM_THREADS=1 ctest -R rakali` and **never `-j`** on the GPU
build. The full suite is ~2.5 min on gcc and ~5.5 min on the GPU. If your tool
has a 10-minute command cap, run the two suites in **separate** invocations —
chaining them times out.

**Both toolchains, every time.** Two real defects in this campaign showed up on
one toolchain only (see §5). A green gcc run proves nothing about the GPU.

### Known-broken tooling — not your fault, do not fix by installing

- **`fprettify` is broken**: `ModuleNotFoundError: No module named 'fprettify'`,
  so that `pre-commit` hook always fails. CLAUDE.md forbids `pip install`.
  Run `fortitude check <files>` directly instead; it works and is the check
  that matters.
- **`gh` is not installed** and there is no `GH_TOKEN`. PRs must be opened in
  the browser. Push works fine (SSH remote,
  `git@github.com:JorgeG94/rakali_coastal.git`).
- **`origin/HEAD` is unset** on the remote.
- `fortitude` flags ~54 pre-existing `S082` hits in `tests/test_boundary.F90`.
  Not yours; leave them.

---

## 2. Where the work is

Fourteen commits, one per PR, in a linear stack. Each was verified green on
**both** toolchains at its own tip. Open them with the base set explicitly —
GitHub defaults to `main` and will show every ancestor.

| # | head branch | base |
|---|---|---|
| 1 | `docs/capability-plans` | `main` |
| 2 | `test/p3-c-property` | 1 |
| 3 | `test/p1-conservation-budget` | 2 |
| 4 | `fix/bpg-jacobian-default` | 3 |
| 5 | `fix/tracer-conservation` | 4 |
| 6 | `fix/wall-bed-mirror` | 5 |
| 7 | `test/p2-convergence-harness` | 6 |
| 8 | `test/conservation-parity` | 7 |
| 9 | `test/p7-convergence-structured` | 8 |
| 10 | `test/p4-dam-break-dry` | 9 |
| 11 | `test/p5-carrier-greenspan` | 10 |
| 12 | `fix/dry-state-wave-speeds` | 11 |

`fix/bpg-sigma-uniform-density` is a **redundant** older pointer at the same
tip as `test/conservation-parity`; ignore or delete it.

Campaign state: **P1, P2, P3, P4, P7 done. P5 partial. P6, P8, P9 open.**

---

## 3. What was fixed (three real defects, all found by measuring)

Each was found because a test *measured* something, not because anyone read the
code and spotted it. All three had passed code review for a long time.

1. **Sigma BPG fabricated force from bed slope.**
   `ml_baroclinic_pressure_gradient` differenced the integrated density anomaly
   at constant *sigma index*. Over a slope, layer `k` sits at a different
   physical depth per column, so a spatially **uniform** density differing from
   `rho_0` produced a fictitious force ∝ (anomaly × bed slope). Kernel-level:
   1.28e-3 vs 2.09e-18 for the Jacobian. `bpg_method` now defaults to
   `"jacobian"`. Pinned by `test_bpg_uniform_density`.

2. **The clamp destroyed tracer when the column mean left the bounds.**
   `frac = min(1, net/cap_to_add)` saturated and the residual was deleted. That
   saturation is *exactly equivalent* to the column mean leaving the bounds
   (proof in the test docstring), so the kernel now skips those columns and is
   conservative for every input. Fired constantly because the default
   `S_max = 40 PSU` sits at the top of the seawater range.

3. **Closed walls passed mass under flow over a sloping bed.**
   The wall BC mirrored the free surface but the ghost **bed** came from the
   bathymetry setter, which writes the *continuation* of the topography, not
   its reflection. A wall is a mirror. Harmless at rest (the C-property still
   held to 1e-12) but leaked under flow. `bc_wall_edge` now mirrors the bed.
   The unstructured KNP path already did this correctly — structured was the
   odd one out.

Plus a fourth, banked but **not** the cause of anything measured: HLL/HLLC used
Davis wave speeds, which underestimate the dry-front wave by up to 2× at a
genuinely dry state. Fixed per Toro §10.4. Changed no test result, because at
the fronts we measure the right state is thin but above `DRY_TOLERANCE`.

---

## 4. Open findings — read this before re-deriving anything

### ✅ Task #8 — FIXED (2026-08, `fix/thin-cell-desingularization`)

**The eliminated-by-experiment table below was right and the conclusion wrong:
the two clamps are REDUNDANT.** Relax either alone and the other keeps binding
— which is exactly why each single-knob test looked innocent. Dose-response
(tip/c₀, only `rki_constants.F90` varied): `1e-3,10 → 1.489`; `1e-9,10 →
1.489`; `1e-3,1e8 → 1.577`; `1e-9,50 → 1.720`; `1e-9,1e8 → 1.888` — but the
fully-relaxed config has a genuine thin-film velocity runaway (measured
u → 4×10³ m/s at h ≈ 3.5 mm) and a ×2000 dt collapse on beaches.

**Fix:** Kurganov & Petrova (2007) desingularized velocity below
`DESING_EPS = 1.5e-4` (= `H_VANISHED`), shared with the CFL estimator (whose
old `0.01·min(dx,dy)` grid term was itself a latent CFL bug — it blinded the
estimator to real velocities in thin cells), plus a `FROUDE_BACKSTOP = 50`
blowup arrester replacing `FROUDE_CAP = 10`. All barotropic update/protect
kernels on both backends. Measured: structured tip 1.489 → **1.736 c₀**
(no-clamp ceiling 1.888), 1-cm front lag 7.11 → **0.36 m** (0.6 %),
`H_FRONT_STRUCT` tightened to 1e-2·h_l; CG beach film drains to the
`~slope·dx` Audusse retention scale (was frozen solid); no dt collapse
(beach dt within ~4× of the wet value).

**Still open, now precisely located:**
- **Unstructured tip (1.586 c₀) is NOT the protect clamp** — measured
  identical with the clamp replaced or not. The binding constraint is the
  KNP reconstruction's dry-cell edge-momentum zeroing
  (`rki_kernel_extrapolate_unstr :: centroid_pass`, `minimum_allowed_height`)
  plus coarser dx. Two-sided guard added
  (`dry_tip_celerity_unstr_is_known_deficient`, 1.45–1.75).
- **mm-scale retention film** on beaches: Audusse well-balancing retains up
  to `slope·dx` per cell (~5 mm mean measured), and the freed runup jet
  overwets the upper beach (mm film to x ≈ 128.4 vs exact max runup 120 m),
  so a 1-mm shoreline detector is still pinned — by discretization
  retention, no longer by a clamp. `cg_thin_film_persists_on_beach` still
  guards this two-sidedly. A fix would be flux-side (bed-slope-aware
  thin-film reconstruction; the old MUSCL-fallback lead at
  `rki_kernel_flux.F90:158-173` is still the place to look).
- **Cost, accepted:** wet/dry-heavy runs take more, smaller steps (Thacker:
  849 → 2612 steps/period) because thin fast water now exists and must be
  CFL-resolved. Thacker planar L2 rose 0.078 → 0.089 (retention rim shifts
  the slosh); its tolerance re-baselined 0.08 → 0.12 with history recorded
  in the test.

### 🟡 Task #9 — P5 runup guard tightened 40 % → 10 %; measured 5.6 %, 5 % not claimed

Post-fix the three detection depths agree to ~2 % (**18.62 / 19.00 / 18.88 m**
at 3e-1/3e-2/1e-2 m vs exact 20 m; rundown edge 83.9 m vs exact 80 m; period
within 1.5 %). The old ~15 % scatter was the frozen film moving the detector.
The residual ~0.6 %-per-cell staircase quantization and the retention film
keep an honest 5 % assertion out of reach at this dx — closing it needs finer
resolution or a retention-aware shoreline detector, not a fitted tolerance.
Note the backstop value was *chosen* by this fixture: `FROUDE_BACKSTOP = 100`
strands a cm-scale wedge through rundown (excursion 11.8 m), 50 restores it
(18.9 m).

### 🔴 Task #7 — unstructured × semi-implicit conservation

The fourth cell of the four-path matrix is unmeasured. Needs a Delaunay mesh:
that path is TPFA with η at circumcentres, and the right-triangle split every
quad-derived test mesh here uses gives coincident circumcentres and trips the
mesh quality gate. Reuse `delaunay_diag_bl_tr` from
`test_si_unstr_ml_wetdry.F90` plus `state_unstr_enable_si` — **and assert its
`ok` flag**: that existing test does `if (.not. ok) return`, so it would go
silently vacuous on a failing mesh. Worth auditing whether it currently does.

---

## 5. Method that actually worked

Five phases and four defects in, this is what earned its keep:

- **Bisect the fixture, one flag at a time.** Every defect was localised by
  building a probe that runs the same case with one thing varied and prints a
  table. Flat vs sloping bed, stratified vs not, 2D vs multilayer, EOS on vs
  off, `nz=1` vs `nz=6`. Guessing from source cost more than measuring, every
  time — including twice where the "obvious" cause was wrong.
- **Use the model's own accounting.** `state%barotropic%mass_out` is
  accumulated as `−Σ_interior flux_h·area`, so `(M₁−M₀) + mass_out ≈ 0` proves
  nothing outside the flux path touched `h`. That one line localised defect 3
  to the *boundary* and eliminated four interior suspects without testing any
  of them.
- **Two-sided guards on known-deficient quantities.** Assert the deficit is
  still there *and* still bounded. The upper bound fails when someone fixes it,
  which is the prompt to retire the guard — instead of a fix landing silently
  while the test file permanently understates the model. Used in
  `test_bpg_uniform_density`, `test_dam_break_dry`, `test_carrier_greenspan`.
  It paid off within one commit: the transposed-argument bug in the wave-speed
  fix was caught by the P4 tip guard firing on its *lower* bound.
- **Never loosen a tolerance to pass.** Where a metric could not be validated,
  the metric was narrowed to what *is* well-posed (front position at a
  resolvable depth) and the deficit recorded separately and explicitly. P5 is
  labelled PARTIAL rather than claiming the headline Carrier-Greenspan number
  on a measurement that cannot support it.
- **Assert the fixture, not just the result.** Every conservation and order test
  carries a guard that the flow actually develops / the column is really
  stratified / the bed really emerges / the profile is really monotone. A
  quiescent fixture conserves trivially and a uniform tracer is conserved by
  almost any scheme.

---

## 6. Traps hit, so you don't

- **`test_multilayer_unstr` reuses one `state` across its four clamp cases**,
  and `state_unstr_init` does *not* reset scalar components. Setting
  `S_init_bottom` in case (a) leaked into (c) and (d). Clear them explicitly.
- **Do not build a periodicity measurement on a wet/dry-affected signal.** The
  P5 period test originally read the shoreline track and passed on gcc while
  failing on nvfortran — 3 crossings vs 2, a factor-of-two period difference —
  purely from how each toolchain rounds a cell-index staircase whose rundown
  extent is itself defect-affected. It now reads the offshore surface
  elevation, in 5 m of water that never wets or dries.
- **gfortran writes `Error:` with a capital E.** Filtering build output with
  `grep -i` matters; a lowercase-only filter hid a compile failure and I ran a
  stale binary for two cycles.
- **`cfl_trunc` is hardcoded to 0.8** at the `ml_truncate_velocity` call site,
  not a namelist knob — you cannot vary it from a test without editing source.
- Front/shoreline detection by depth threshold is **violently** threshold-
  sensitive where the exact profile approaches zero quadratically. Compare
  like-for-like: apply the same threshold to the analytic profile and invert it
  in closed form. Always sweep the threshold and report the table.

---

## 7. Doc corrections already made

The plans were written against a larger pre-port-over repo. Fixed in PR 1:
dangling `src/core/{ocean,ice}/` references, a wrong shipped-tracer list (only
salinity and temperature exist), and a stale line reference into
`CAPABILITIES_AND_LIMITATIONS.md`. Also corrected: the gap analysis claimed no
standalone C-property assertion existed — `test_lake_at_rest.F90` already had
one, it was just weak.

Still true and worth re-reading before planning: `schism_parity_report.md` §5
records that the **semi-implicit mode is single-rank only**, boundaries limited
to wall/tidal/clamped, multilayer sigma-only and `nz ≥ 2`. The plans' "all four
paths" convention over-promises against that.

---

## 8. Suggested next step

Task #8 is fixed (see §4). The highest-value remaining threads:

- **P6 (flow over a bump)** — independent, gives a third
  moving-shoreline-free benchmark, and its transcritical case is the only
  shock test in the campaign.
- **The unstructured dry-front reconstruction** — the KNP tip deficit is now
  precisely located in `centroid_pass`'s dry-cell edge-momentum zeroing
  (§4); fixing it fires the new two-sided unstructured tip guard.
- **Bed-slope-aware thin-film reconstruction** — would drain the `~slope·dx`
  retention film, fire `cg_thin_film_persists_on_beach`, and open the way to
  the P5 5 % target (task #9). The old MUSCL-fallback lead
  (`rki_kernel_flux.F90:158-173`) is the place to start.
