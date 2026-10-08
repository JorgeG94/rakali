# Rakali Codebase Map

GPU-native 2D/3D coastal solver. The canonical source is Fortran
`do concurrent` + OpenACC (portable across NVHPC GPU/multicore, gfortran, ifx);
an OpenMP-target variant for Intel/AMD GPUs is generated from that same source by
the acc→omp transformer (`tools/acc_to_omp.py` / `tools/dc_to_omp.py`, pushed to
the `auto/*` branches).

> **Where to read first** depends on what you're working on:
>
> - **New to the codebase / acronyms**: [CONCEPTS.md](CONCEPTS.md)
> - **Design contract (HLL/HLLC + KNP, production)**: [`src/core/coastal/README.md`](../../src/core/coastal/README.md) — conventions + the RK2 dispatch pipeline
> - **Which closure/scheme is enabled where + tunable knobs**: [`docs/CLOSURE_MATRIX.md`](../CLOSURE_MATRIX.md)
> - **Capabilities + limitations**: [`docs/CAPABILITIES_AND_LIMITATIONS.md`](../CAPABILITIES_AND_LIMITATIONS.md)
> - **Physics + namelist reference**: [`docs/REFERENCE.md`](../REFERENCE.md)
> - **Per-module / procedure detail**: the source `!!` docstrings, rendered by **FORD** (built on PR to `main`)

## Directory Structure

```
rakali/
├── app/
│   └── main.F90                        Thin entry — calls driver_run
├── src/
│   ├── core/                           Types, config, state, init, profiler
│   │   ├── structured/                  Shared sub-state types:
│   │   │                                  rki_grid (hgrid_t), rki_surface_forcing,
│   │   │                                  rki_barotropic_state, rki_multilayer_state,
│   │   │                                  rki_nonhydrostatic_state, rki_state,
│   │   │                                  rki_tracer, rki_boundary_types
│   │   ├── unstructured/                Unstructured analogues:
│   │   │                                  rki_mesh_unstr, rki_state_unstr, rki_tracer_unstr,
│   │   │                                  rki_boundary_unstr, rki_surface_forcing_unstr,
│   │   │                                  rki_barotropic_state_unstr, rki_multilayer_state_unstr
│   │   ├── coastal/                    the sim_type='coastal' solver + kernels
│   │   │   ├── solver/                  rki_solver, rki_solver_unstr
│   │   │   └── kernels/
│   │   │       ├── structured/
│   │   │       │   ├── barotropic/      HLL/HLLC + boundary + sponge + fused
│   │   │       │   ├── multilayer/      ml_dynamics, ml_truncate, kernel_fused_ml
│   │   │       │   └── nonhydrostatic/  CG Poisson solver
│   │   │       └── unstructured/
│   │   │           ├── barotropic/      KNP central-upwind (+ SI TPFA/MG)
│   │   │           └── multilayer/      ml_dynamics_unstr (+ SI column)
│   │   └── multilayer/{structured,unstructured}/  rki_ml_common, rki_ml_truncate (shared ML helpers)
│   │
│   ├── ALE/                            Shared: rki_remap_column + rki_vcoord (vertical coords + remap)
│   ├── tracer/                         Shared: rki_ml_tracers + unstr (advection, CWC, clamp) + rki_recon_weno
│   ├── equation_of_state/              rki_ml_eos (hardwired linear) + rki_eos (shared Wright/Roquet nonlinear)
│   ├── pressure_force/                 rki_ml_pressure (baroclinic pressure gradient)
│   ├── parameterizations/
│   │   ├── vertical/                    rki_ml_vertical_viscosity (PP81), rki_ml_friction,
│   │   │                                rki_ml_keps (+ stability), rki_ml_kpp
│   │   └── lateral/                     rki_ml_horizontal_viscosity (Smagorinsky + biharmonic)
│   ├── framework/                      rki_scratch_3d (scratch_3d_buffer_t), rki_safe_math,
│   │                                     rki_efp (reproducing sums), rki_mem_report
│   ├── driver/                         rki_driver — full run lifecycle
│   ├── api/                            C-binding: rki_handle, rki_api, rki_api_unstr
│   ├── io/                             NetCDF wrapper + bathymetry; structured/ + unstructured/ subdirs
├── tests/                             Flat test_*.F90 — testdrive unit + analytical
├── benchmarks/                         Single-GPU benchmarks
├── python/                             rakali._core, rakali.unstr, rakali.mesh (Hilbert SFC), examples/
├── tools/                              Lint scripts, post-processing, OpenMP↔OpenACC translators
├── docs/                               REFERENCE.md (physics + namelist), codebase/ map,
│                                       ROADMAP_COASTAL.md, v0.1.0_checklist.md, …
├── cmake/                              Compiler flags + helper modules
├── CLAUDE.md                           Agent-facing project instructions
├── FORTRAN_STYLE.md                    Coding style guide
└── CMakeLists.txt                      Build system entry point
```

## Code-level reference → FORD

Per-module and per-procedure documentation comes from the `!!` docstrings in the
source, rendered by **FORD** (built on PR to `main`). We deliberately do **not**
hand-maintain a parallel set of markdowns describing what each module/kernel
does — that's the copy that rots (the previous `KERNELS.md` / `CORE.md` /
`API.md` / … did exactly that). The docstring next to the code stays current.

The surviving **curated / cross-cutting** docs — the synthesis FORD can't
generate — are:

| Document | Covers |
|----------|--------|
| [CONCEPTS.md](CONCEPTS.md) | Acronyms (HLL, KNP, CWC, BPG, MUSCL, PP81, ALE, …) + one-line concept definitions. Start here. |
| [`src/core/coastal/README.md`](../../src/core/coastal/README.md) | Design contract: vertical-layer convention, RK2 pipeline, closure dispatch, conventions. |
| [`docs/CLOSURE_MATRIX.md`](../CLOSURE_MATRIX.md) | Which closure/scheme is enabled on which backend + the tunable knobs. |
| [`docs/CAPABILITIES_AND_LIMITATIONS.md`](../CAPABILITIES_AND_LIMITATIONS.md) | What the solver can / can't do today. |
| [`docs/REFERENCE.md`](../REFERENCE.md) | Physics, numerics, namelist keys. |
| [`FORTRAN_STYLE.md`](../../FORTRAN_STYLE.md) | Style + the GPU (`do concurrent` + OpenACC) programming guide. |

## Run lifecycle (`app/main.F90` → `driver_run`)

```
read_config(namelist)
state%init_from_config           allocate
load_bathymetry → fill ghosts
init_from_config                 initial conditions
forcing / nesting / gauges init
solver_enter_data                map state to the GPU (once)

  TIME LOOP:
    solver_step                  one SSP-RK2 step (pipeline in the design-contract README)
    output / gauges / restart     at their cadences

solver_exit_data                 map back from the GPU (once)
output_close / finalize
```

## Notes on the layout

The mental model:

- **The solver**: `src/core/coastal/` — the A-grid solver + kernels
  (structured + unstructured, barotropic / multilayer / non-hydrostatic /
  semi-implicit).
- **Shared physics + numerics**: `src/{ALE, tracer, equation_of_state, pressure_force, parameterizations, framework}/` — infrastructure the structured and unstructured backends both use.

The design contract (vertical-layer convention, the RK2 dispatch pipeline, the
closure selector) lives in `src/core/coastal/README.md`; the long-horizon
backlog is in `docs/ROADMAP_COASTAL.md`.
