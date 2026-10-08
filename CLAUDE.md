# Rakali — GPU-Native Coastal Solver

Fortran solver for coastal, estuarine, and tidal applications. Structured
Cartesian grids (2D barotropic, 3D hydrostatic multilayer, non-hydrostatic) and
unstructured triangular grids (2D barotropic, 3D hydrostatic multilayer, KNP
central-upwind). GPU parallelism via `do concurrent` + OpenACC; portable across
NVHPC (GPU/multicore), gfortran, ifx.

> **Where things are.** Acronyms → [`docs/codebase/CONCEPTS.md`](docs/codebase/CONCEPTS.md) · Which closure/scheme is enabled in which regime + the tunable knobs → [`docs/CLOSURE_MATRIX.md`](docs/CLOSURE_MATRIX.md) · Design contract (vertical-layer convention, RK2 pipeline, closure dispatch) → [`src/core/coastal/README.md`](src/core/coastal/README.md) · Capabilities + limits → [`docs/CAPABILITIES_AND_LIMITATIONS.md`](docs/CAPABILITIES_AND_LIMITATIONS.md) · Physics + namelist → [`docs/REFERENCE.md`](docs/REFERENCE.md) · Style + GPU (`do concurrent` + OpenACC) programming → [`FORTRAN_STYLE.md`](FORTRAN_STYLE.md) · Repo map + run lifecycle → [`docs/codebase/INDEX.md`](docs/codebase/INDEX.md).
>
> **Per-module / per-procedure docs come from the source `!!` docstrings — rendered by FORD (built on PR to `main`). Don't hand-write parallel markdown describing what a module/kernel does; that's the copy that rots.**

## Build

```bash
module load cmake nvhpc
cmake -B build -S . && cmake --build build
cd build && ctest --output-on-failure
```

Key CMake options: `RAKALI_ENABLE_GPU` (default **OFF** — opt-in; enable with `-DRAKALI_ENABLE_GPU=ON` on the NVHPC toolchain. An explicit ON with a compiler lacking a GPU path is a configure-time error, never a silent CPU downgrade; the value prints on the `GPU offload:` configure line), `RAKALI_ENABLE_PYTHON` (default OFF), `RAKALI_ENABLE_NETCDF` (default ON; disable for kernel/solver-only portability builds), `RAKALI_GPU_ARCH` (cc70/80/90), `RAKALI_ENABLE_DOUBLE` (default ON). Full table in `CMakeLists.txt` and `README.md`.

## Local development conventions (agents + contributors)

Two standing rules for anything run against this repo — humans and AI agents alike:

- **Never `pip install` (or otherwise install packages).** Use only what the
  toolchain already provides (the loaded module environment). If a Python helper
  needs a package that isn't present, rewrite it against the standard library — a
  missing dependency is a signal to simplify, not to install.
- **All local/throwaway artifacts go in `tmp_local_artifacts/` at the repo root** —
  prototypes, scratch scripts, generated data, plots, logs, profiler/benchmark
  dumps, intermediate command output. It is git-ignored, so nothing leaks into a
  commit, and every worktree + every agent writes to the same known place. Do
  **NOT** write to `/tmp`, `/private/tmp`, or anywhere else in the tree —
  `tmp_local_artifacts/` is the *only* place. (Committed tooling under `tools/`
  is the standing exception — that's real code, not an artifact.)
- **Build directories are named `build_<toolchain>`** (`build_serial`,
  `build_cc90`, ...) — `.gitignore` covers `build*/` and `*_build/`; a build
  tree outside those patterns is how a stray `git add -A` once poisoned the
  object store. Never name a build dir anything else.

## Pre-commit workflow

Before each commit, in order:

1. **`pre-commit run --all`** — formatting (whitespace, trailing newlines, end-of-file).
2. **Fortitude lint** — Fortran static checks (`fortitude check` or via the pre-commit hook).
3. **`ctest`** — full test suite. **NEVER pass `-j N` when the GPU build is active** — every worker shares one GPU, so parallel test execution causes spurious failures + hangs. Building with `-j` is fine; only `ctest -j N` is the problem. Use `ctest -R rakali` (the rakali-only subset) for fast regression checks.
4. **Validate the docs — always, before merging anything.** FORD owns per-module / per-procedure docs from the source `!!` docstrings (don't hand-write those). But the hand-maintained **synthesis** docs rot silently and a stale cell is a bug: re-read [`docs/CLOSURE_MATRIX.md`](docs/CLOSURE_MATRIX.md) (drift-checked by `tools/check_closure_matrix.py`), the design-contract README (`src/core/coastal/README.md`), the `docs/howto/` extension guides, and [`docs/CAPABILITIES_AND_LIMITATIONS.md`](docs/CAPABILITIES_AND_LIMITATIONS.md). If your change adds, ports, removes, or renames a closure / scheme / knob / extension seam, update the relevant doc **in the same PR** — and when a doc and the code disagree, the code is authority and the doc is the bug to fix.

When adding a new capability:

1. **Namelist knob** with default = off (preserves bit-identity for existing nmls + tests).
2. **Kernel implementation** gated on the knob.
3. **Unit test** covering the new code path (analytical tests are high-leverage — every one added to the core has caught a real bug).
4. **Closure-matrix row** in `docs/CLOSURE_MATRIX.md`.
5. **`pre-commit run --all` + fortitude** — pass before committing.
6. **`ctest`** (no `-j`) green on the GPU build.
7. **One commit per capability; one PR per capability.** Smaller PRs review faster and bisect cleanly.

## Project Layout

```
src/
  core/
    rki_{constants,config,config_schema,init,forcing,nesting,profiler,timestepper,array_utils}.F90  infra
    rki_{console_stats,nml_schema}.F90   status reporting + strict namelist schema
    structured/    grid, *_state (barotropic, multilayer, nonhydrostatic),
                   boundary_types, tracer, surface_forcing
    unstructured/  *_state_unstr, mesh_unstr, mesh_gen_unstr, boundary_unstr,
                   tracer_unstr, surface_forcing_unstr, si_unstr_state
    multilayer/{structured,unstructured}/  shared ML state/helpers (ml_common, ml_truncate)
    coastal/       the A-grid path
      solver/        rki_solver, rki_solver_unstr
      kernels/structured/{barotropic,multilayer,nonhydrostatic}/
      kernels/unstructured/{barotropic,multilayer}/
  # shared physics — each split structured/ + unstructured/
  equation_of_state/   rki_ml_eos (coastal linear) + rki_eos (shared Wright/Roquet nonlinear)
  pressure_force/      rki_ml_pressure (baroclinic pressure gradient)
  tracer/              rki_ml_tracers (advection, CWC, clamp) + rki_recon_weno
  ALE/                 rki_remap_column (PCM/PLM/PPM/PQM) + rki_vcoord (vcoord + remap)
  parameterizations/
    vertical/   PP81 vertical_viscosity, friction, keps (+stability), kpp
    lateral/    ml_horizontal_viscosity (Smagorinsky + biharmonic)
  framework/     rki_scratch_3d, rki_safe_math, rki_efp (reproducing sums), rki_mem_report
  driver/        rki_driver — driver_run (testable, see test_driver_smoke)
  api/           rki_handle, rki_api, rki_api_unstr (opaque-handle FFI)
  io/            banner, bathymetry, io_netcdf + structured/ + unstructured/ (output, gauges)
app/main.F90     Thin entry — calls driver_run
python/rakali/   _core.py (structured), unstr.py, yaml_loader.py, mesh/ (builder, Hilbert renumber)
tests/           test_*.F90 (unit + analytical, flat)
benchmarks/      bench_single_gpu, bench_multilayer, bench_eos, … (+ str/unstr variants)
validation_examples/    canonical / analytical benchmarks (coastal/)
python/examples/       Python examples (self-contained) + real_bathymetry/ case studies (bring your own data)
```

The multilayer kernels are split by **physics concern** across the shared tree,
each mirrored `structured/` + `unstructured/`: `rki_ml_common`/`rki_ml_truncate`
(`core/multilayer/`), `rki_ml_dynamics` (`core/coastal/kernels/.../multilayer/`),
`rki_ml_eos` (`equation_of_state/`), `rki_ml_pressure` (`pressure_force/`),
`rki_ml_tracers` (`tracer/`), `rki_ml_{vertical_viscosity,friction,keps,kpp}`
(`parameterizations/vertical/`), `rki_ml_horizontal_viscosity`
(`parameterizations/lateral/`).

## Architecture

`sim_type='coastal'` (the only regime): HLL/HLLC (structured) + KNP
central-upwind (unstructured) on an A-grid. Everything lives under
`src/core/{coastal,structured,unstructured,multilayer}/` plus the shared
`src/{tracer, ALE, pressure_force, equation_of_state, parameterizations}/`.

An opt-in `time_scheme="semi_implicit"` mode (Casulli CG Helmholtz η-solve, no
gravity-wave CFL) is available on the **unstructured** backend only (a structured
grid with `semi_implicit` is rejected at configure time), 2D barotropic and 3D
sigma multilayer — a TPFA/circumcentre scheme:
`rki_si_unstr_{step,helmholtz,mg,wetdry}` in
  `kernels/unstructured/barotropic/`, `rki_si_unstr_{step_ml,column,closures}`
  in `kernels/unstructured/multilayer/`, and `rki_si_unstr_state` in
  `src/core/unstructured/`.

See `docs/REFERENCE.md` and `docs/CAPABILITIES_AND_LIMITATIONS.md` for the
honest limits.

### Two Backends

**Structured** (Cartesian): HLL/HLLC Riemann + hydrostatic reconstruction,
minmod limiter, SSP-RK2. SoA `(nx,ny[,nz])` for coalesced GPU access. Supports 2D
barotropic, 3D hydrostatic multilayer, and non-hydrostatic (CG Poisson,
structured only).

**Unstructured** (triangular): Kurganov-Noelle-Petrova central-upwind, Green-Gauss
gradients, SSP-RK2. 2D barotropic + 3D hydrostatic multilayer. Mesh from UGRID
NetCDF.

Both backends share the same time integrator, vertical-coordinate types, EOS, and
tracer registry.

### State Composition

`state_t` is a thin composition — each sub-type owns its allocatables and bound `init`/`destroy`. Access paths:

```
state_t
  ├── grid        : hgrid_t              — nx_total, ny_total, nghost, dx, dy
  ├── barotropic  : barotropic_state_t   — h, hu, hv, b, RK saves, flux_*, manning_n, coriolis_f
  ├── forcing     : surface_forcing_t    — wind_x/y, atm_pressure, precip, heat_flux, evap_minus_precip
  ├── multilayer  : multilayer_state_t   — h_layer, hu/hv_layer, w_ml, rho_layer, rho_ref,
  │   └── tracers : tracer_t (per-tracer)  vcoord, dsig_ml, tracers(:), idx_salinity, idx_temperature, ...
  ├── nh          : nonhydrostatic_state_t — u/v/w_layer, q_nh, CG workspaces, mode-split arrays
  └── bc          : bc_params_t          — BC config + Chapman persistent state
```

No flat `state%h` / `state%h_layer` / `state%manning_n` / `state%wind_x` aliases remain. Unstructured mirrors with `unstr_state_t%tracers(:)` + `idx_salinity`/`idx_temperature` (other fields still flat on `unstr_state_t`).

### GPU Parallelism

All data-parallel loops use `do concurrent`. Reductions use `!$acc parallel loop reduction(...)` (inert comment on non-OpenACC compilers). `!$acc enter/exit data` / `!$acc update` is used directly — no `#ifdef` wrapping. State arrays stay device-resident; only diagnostic output triggers D→H transfers. `RAKALI_GPU_OFFLOAD` is the preprocessor define set by `RAKALI_ENABLE_GPU=ON`.

### Vertical Layer Convention (load-bearing)

Both backends are **bottom-up (ROMS-style)**: `k=1` is the bed, `k=nz` is the surface layer. Surface forcings land at `(:,:,nz)` / `(nz, cell)`; bed forcings at `(:,:,1)` / `(1, cell)`. Any new multilayer kernel must follow this — BPG, rho_ref, NH init, per-layer friction, vcoord output, surface tracer fluxes, vertical exchange all assume it. Three regression gates: `test_lock_exchange` (dense current at `k=1`), `test_river_plume` (buoyant plume at `k=NZ`), `test_surface_flux_strat` (positive Q warms `k=NZ` and preserves stratification).

### Multilayer (both backends)

Mode-split barotropic-then-distribute. Each SSP-RK2 stage: barotropic flux → distribute to layers → tracer advection + diffusion → EOS → vertical mixing → BPG → friction → Smagorinsky → truncation → per-layer Coriolis → recompute barotropic. Stage 2 also RK2-averages and clamps tracers.

Two prognostic tracers registered by default: salinity (`tracers(idx_salinity)%hTr`, PSU·m) and temperature (`tracers(idx_temperature)%hTr`, °C·m). On **both** backends every passive-transport kernel (advection, ALE remap via `ml_remap_tracer{,_unstr}`, vertical exchange, vertical/horizontal diffusion, clamp, RK2) iterates the registry — append a passive tracer with `multilayer%register_passive_tracer(...)` and it "just works" (EOS + surface flux stay S/T-specific by design). The split is the same on each side: a momentum/thickness `ml_distribute_with_remap{,_unstr}` + a per-tracer `ml_remap_tracer{,_unstr}` looped over `tracers(:)`. Structured shape `(nx,ny,nz)`, unstructured `(nz, ncells)` (opposite axis order). Caveat: the structured nz=1 fused fast path still carries S/T only (extra tracers need nz>1).

Advection uses **vanishing-layer-gated centred CWC** per-face mass-flux split with per-face symmetric redistribution — guarantees per-layer conservation under deep-shallow transitions (ZSTAR_FULL); reduces to centred CWC for benign bathymetry. See `docs/REFERENCE.md` for the rest of the physics (BPG, PP81, KPP, k-ε, Smagorinsky, EOS, Coriolis, Manning, surface fluxes, clamp).

Unstructured ports the per-layer stabilisers: **velocity truncation**
(`ml_truncate_velocity_unstr`, CFL clip to `0.9*cfl_trunc*radius/dt` using the
inscribed-circle CFL radius as the isotropic length scale) and **per-layer bottom
friction** (`ml_apply_bottom_friction_unstr`, semi-implicit Manning quadratic on
the lowest active layer; the unstructured barotropic state carries no linear-drag
rate, so this is Manning-only — the structured twin's `use_linear_drag` branch
has no unstructured counterpart). Both slot in before the per-layer Coriolis in
each RK2 stage, mirroring the structured order (BPG → friction → Smagorinsky →
truncation → Coriolis → recompute). PP81, KPP, k-ε, Smagorinsky, Laplacian
h-diffusion, conservative clamp, per-layer Coriolis are also ported. **The fused
ML kernel is deliberately NOT ported**: an nsys profile of the unstructured ML
bench (20k tris, nz=5, V100) shows the step is dominated by tracer advection +
vertical/horizontal diffusion (~68% of GPU time); the fusible post-update chain
(distribute + truncate + recompute + Coriolis + friction) is only ~6.8%, and since
those kernels are launch/bandwidth-bound the realisable fusion win is a fraction of
that — well under the 5% threshold that would justify the divergence cost.

### Vertical Coordinates

Dispatched via `vcoord_type` (namelist) → `parse_vcoord_type` → `VCOORD_*` enum.

- `VCOORD_SIGMA` — terrain-following, no remap. Default.
- `VCOORD_ZSIGMA` — smoothstep blend sigma→z. Uses conservative remap.
- `VCOORD_ZSTAR` (z*-lite) — single global `z_ref` stretched per column. SSH-tracking. Conservative remap.
- `VCOORD_ZSTAR_SIGMA` — sigma in shallow water, z*-lite in deep. Conserves by construction.
- `VCOORD_ZSTAR_FULL` — per-column `z_ref(0:nz)` from local bathymetry; surface layer anchored at `zstar_h_surf_target` regardless of H. Bed-side layers can vanish (`zstar_h_min`). Operators that divide by `h_layer` gate on `VANISHING_LAYER_TOL = 1.5e-4 m`. Pairs with gated CWC + conservative clamp. **Caveat:** intertidal estuaries still leak 1-2% salt/cycle from wet/dry destruction. Recommend sigma or zstar-lite for those.

The conservative integrate-and-redistribute remap (`rki_remap_column`) is shared
across all coords; only the per-column reconstruction order changes
(`remap_method` = pcm/plm/ppm/ppm_h4/pqm, default `"ppm"`).

### Non-Hydrostatic (Structured Only)

Sigma coords + CG Poisson. Uses `NH_MIN_DEPTH` (not `DRY_TOLERANCE`) for depth guards — separate thresholds. Mode splitting via `split_ratio` (1=unsplit, 0=auto, N=fixed). All state at `state%nh`.

### Boundary Conditions

Nine types: BC_WALL(1), BC_OPEN(2/Flather), BC_TIDAL(3), BC_NESTED(4), BC_INFLOW(5), BC_DISCHARGE(6), BC_CLAMPED(7), BC_SPONGE(8), BC_CHAPMAN(9/Orlanski). Sponge ghost falls through to wall; sponge relaxation is a separate kernel. Chapman uses persistent `eta_old_*` (structured: on `bc_params_t`; unstructured: on `unstr_state_t%chapman_eta_old`). Module split: `rki_boundary_types` (backend-agnostic types/constants) vs `rki_boundary` (kernel). Unstructured supports wall, open, inflow, discharge, outflow, tidal/Flather, clamped, and Chapman (inlined into `rki_kernel_flux_unstr`), plus `BC_SPONGE` (separate `apply_sponge_unstr` relaxation kernel — sponge ghost falls through to wall, called after each RK2-stage barotropic update; v1 single-cell band, no graded width on triangles) and `BC_NESTED` (smoke-level: routed through the Flather/`BC_TIDAL` algebra with a per-group prescribed elevation; live parent-NetCDF coupling still deferred).

### I/O

Output files (`output_rank_000000.nc`; the rank suffix is a leftover of the removed MPI support). Optional deflate (3-5x). Gauges: bilinear, CF-1.8, per-gauge D→H copy.

### Python/Julia API

`librakali.so` via `iso_c_binding`. Opaque-handle FFI — multiple handles (incl. mixed structured + unstructured) coexist in one process. Zero-copy field access via `c_loc` → numpy views. Two init modes: one-shot `rki_create_from_namelist`, or two-phase `rki_create` → `set_*` → `rki_finalize`. Mid-sim callbacks use lightweight `!$acc update self` (data stays on device).

## Gotchas

- **Unstructured flux pattern**: cell-centric double-visit. Each triangle accumulates only into its own arrays — interior edges computed twice. No scatter, no atomics. **Any new unstructured kernel must follow this** — writing to neighbours = GPU race.
- **Vertical indexing**: see "Vertical Layer Convention" above. Surface = `nz`, bed = `1`. Don't flip.
- **Vanishing-layer constant**: coastal `VANISHING_LAYER_TOL`(_UNSTR) (1.5e-4 m) lives in `rki_ml_common`(`_unstr`); operators that divide by `h_layer` gate on it, and the massless-merge helper is `rki_massless`. `DRY_TOLERANCE` / `NH_MIN_DEPTH` (in `rki_constants`) are the separate wet/dry and NH depth guards — don't conflate them.
- **Tests with NetCDF/HDF5**: run with `OMP_NUM_THREADS=1` (via CTest `ENVIRONMENT`) — non-thread-safe HDF5 + `-stdpar=multicore` crashes otherwise.
- **Writing GPU tests (even with no GPU on the dev box)**: the GPU build is `-gpu=...,mem:separate` (**no** managed/unified memory), so a `do concurrent` / `!$acc` kernel gets **NO** implicit host↔device copies — *every* array a kernel touches must be device-present, or it silently reads/writes stale host memory (symptoms: budget/accumulator stays `0`, or division by an unmapped scratch array yields `NaN`, with no crash). A test that only ever runs on the multicore/host build will pass while being GPU-broken. Rules of thumb: **(1)** before calling any kernel, map every state object AND its scratch companion. **(2)** Arrays mapped `create` by `enter_data` (fluxes the production step recomputes on-device) do **not** carry host values you set before the map — push them with `!$acc update device(...)` after `enter_data`. **(3)** Host-set inputs read back after a kernel need `!$acc update self(...)`; local scratch passed straight to a `*_impl` needs its own `!$acc enter data copyin/create ... exit data`. All these directives are inert no-ops on host/multicore builds, so add them unconditionally. **Verify on the actual GPU build** — a green multicore run proves nothing about device data motion.
- **Never `!$acc update self` / `copyin` a WHOLE derived type with allocatable components**: the aggregate D→H copy overwrites the host component descriptors with DEVICE addresses, so the next host read of a component (`obj%arr(...)`) segfaults. Always update/map the COMPONENT arrays (`!$acc update self(obj%arr_x, obj%arr_y)`), never the aggregate `obj`.
- **Persistent kernel workspaces**: never local-allocate scratch in per-step kernels on `-stdpar=gpu` (stalls device on entry/exit). Structured: module-level allocatables, lazy-allocate via `*_workspace_ensure`, release via `*_cleanup` in `solver_exit_data`. Unstructured: fields on `unstr_state_t`, mapped in `state_unstr_enter_data`. Use `associate` to keep in-loop names short.
- **Outer-shim + flat-impl call sites** (for the array-of-derived-types device indirection rule): `compute_layer_sums`, `ml_equation_of_state_unstr`, `ml_remap_tracer_unstr`. Top-level allocatables (`state%barotropic%h`, `state%multilayer%h_layer`) reach the device directly; only the registry indirection needs the shim (the momentum-only `ml_distribute_with_remap{,_unstr}` no longer does).
- **Cross-TU helper inlining** (`_impl` performance variant): NVHPC's device codegen does not inline `pure !$acc routine seq` helpers across module boundaries. For hot inner loops profiled >5% GPU time, duplicate helper bodies into the calling kernel module as `*_impl` copies. Already applied to `rki_kernel_flux_unstr` / `rki_kernel_extrapolate_unstr`. **Don't preemptively duplicate** — constants and bandwidth-bound kernels don't benefit. Pay divergence cost only when profile justifies.
- **Hoist neighbour-column builds out of per-layer loops**: any kernel shaped `for cell { for k { for nbr { build_column(nbr); query_at_z(k); }}}` should build the neighbour column once per cell, then query per-k. Done for Jacobian BPG kernels (2-2.4× per kernel).
- **`if/else` clamps LAUNDER NaN under `-fast` (host AND device)**: nvfortran's
  relaxed-FP default (+`-fast`, no `-Kieee`) lowers `if (u > hi) u = hi; else if
  (u < lo) u = lo` to a NaN-blind min/max select — a NaN input comes OUT as
  `lo`/`hi` (`-Kieee` preserves NaN but costs ~10-20%). Consequence: a clamp can
  silently convert corruption into plausible extreme values. Any clamp that can
  ever see non-finite data must be `ieee_is_finite`-guarded; comparisons with NaN
  are FALSE, so `if (abs(u) > thresh)` guards also silently skip NaN — pair them
  with the catch, never rely on them to bound corrupted data.
- **NO managed/unified memory** (`-gpu=...,mem:separate`): a host-allocated array written by a `do concurrent` and then read via `sum()` on the host silently reads zeros. Use `!$acc parallel loop reduction(...)` for on-device reductions.
- **Formula bathymetry setters must fill ghost rows**: `set_bathymetry_*` formula paths that leave `h_layer=0` at ghost rows make EOS fall back to ρ=ρ₀ → a spurious density jump at wall-adjacent faces → blowup. The file-loader path fills ghosts; formula paths must too.
- **Assumed-shape dummies in `do concurrent` kernels**: NVHPC walks the descriptor with per-launch memcpys; use explicit-shape args (`arr(nx,ny,nz)`); cadence-bounded code may waive with `! assumed-shape-ok: <reason>`; enforced diff-aware by pre-commit (`dc-assumed-shape`).
- **Declare integer dims before the explicit-shape arrays that use them**: ifx warns #8586 ("non-standard out-of-order declaration") when an array's bound references a dummy integer declared on a later line. gfortran/nvfortran accept it silently. Fix: move `integer, intent(in) :: nx, ny, nz` above the first `real(...) :: arr(nx, ny, nz)` in every procedure spec part; enforced by the `decl-order` pre-commit hook.

## Coding Conventions

- File prefix: `rki_` for modules, `test_` / `test_rki_` for tests. Extensions: `.F90`.
- Types: `_t` suffix. Constants: UPPERCASE. Everything else: snake_case.
- `use` always with `only:`. `implicit none`, `private` by default, explicit `public`. FORD docs: `!!` after declarations.
- Intrinsic modules need the modifier: `use, intrinsic :: iso_fortran_env, only: ...` (fortitude rule C122). Same for `iso_c_binding`.
- No `print *` — use `pic_logger` (`global_logger`).
- **Default new procedures to `pure`** — make a new function/subroutine `pure` unless it genuinely needs a side effect (I/O, logging, mutating module state, lazy `allocate`). Pure is the norm here; non-pure is the justified exception. See `FORTRAN_STYLE.md` §Pure and Elemental Procedures.
- Working precision: `wp` from `rki_constants` (currently `real64`).
- Dependencies: `pic` (types, strings, timers, logger), `test-drive`, NetCDF-Fortran.

## Performance (V100)

| Config | Grid | GPUs | Throughput |
|---|---|---|---|
| 2D barotropic | 1652×1396 | 1 | ~707 Mcells/s |
| 3D multilayer (5 σ) | 1652×1396 | 1 | ~130 Mcells/s |
| 3D multilayer bench | 2000×2000 | 1 | ~340 Mcells/s |
