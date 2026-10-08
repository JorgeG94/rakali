# Extending Rakali — contributor how-to guides

Recipe-style walkthroughs for the common extension points, written so they work
whether you're editing by hand or driving an LLM coding assistant (Claude Code,
etc.). Each guide traces one *seam*: the dispatch mechanism, the ordered list of
files to touch, the test to copy, and the GPU gotchas that bite.

These are **process docs**. They deliberately do *not* restate what each module
or procedure does — that lives in the source `!!` docstrings and is rendered by
FORD on every PR to `main`. When a guide says "the kernel computes X", go read
the kernel; the guide's job is to tell you *where the seam is and what to edit*.

> Symbol names (modules, procedures, enums, namelist groups) are stable and are
> what these guides anchor on. Line numbers drift — treat any you see in git
> history or agent output as "approximately here", and grep for the symbol.

## Read these first (orientation, not duplicated here)

| You want… | Read |
|---|---|
| Repo map + run lifecycle | [`docs/codebase/INDEX.md`](../codebase/INDEX.md) |
| Acronyms / concepts | [`docs/codebase/CONCEPTS.md`](../codebase/CONCEPTS.md) |
| Which closure/scheme is on where + the knobs | [`docs/CLOSURE_MATRIX.md`](../CLOSURE_MATRIX.md) |
| Design contract (layer convention, RK2 pipeline, dispatch) | [`src/core/coastal/README.md`](../../src/core/coastal/README.md) |
| Physics + namelist reference | [`docs/REFERENCE.md`](../REFERENCE.md) |
| Style + GPU programming (`do concurrent` + OpenACC) | [`FORTRAN_STYLE.md`](../../FORTRAN_STYLE.md) |
| Capabilities + limits | [`docs/CAPABILITIES_AND_LIMITATIONS.md`](../CAPABILITIES_AND_LIMITATIONS.md) |

## The guides

| Guide | Adds to | Has a dispatch seam already? |
|---|---|---|
| [Add a namelist knob](add_namelist_knob.md) | config plumbing | n/a — the config plumbing itself |
| [Add a passive tracer](add_passive_tracer.md) | multilayer | yes — the tracer registry "just works" |
| [Add a vertical coordinate](add_vertical_coordinate.md) | multilayer | yes — `VCOORD_*` enum + target-thickness generators |
| [Add an equation of state](add_equation_of_state.md) | multilayer | no — you build the scaffold |
| [Add a turbulence closure](add_closure.md) | multilayer | yes — the vmix `if/else if` switch |
| [Add a boundary condition](add_boundary_condition.md) | both backends | yes — the `BC_*` enum + per-edge dispatch |

## Two-minute orientation

- **The solver core** lives under `src/core/coastal/`,
  with shared physics in `src/{equation_of_state,pressure_force,tracer,ALE,parameterizations}/`.
  Many physics modules ship a `structured/` and an `unstructured/` twin.
- **Two grid backends**: structured (Cartesian, SoA arrays)
  and unstructured (triangular, cell-centric). Most extensions need both.
- **State composition**: `state_t` is a thin composition of sub-state types, each
  owning its allocatables + bound `init`/`destroy`. No flat aliases. See the
  state-composition diagram in the top-level [`CLAUDE.md`](../../CLAUDE.md).

## Anatomy of a capability (the standard workflow)

Every guide is a specialisation of this. From the top-level `CLAUDE.md`:

1. **Namelist knob, default = off.** Preserves bit-identity for existing namelists
   and tests. This is non-negotiable — every existing run must stay byte-for-byte.
2. **Kernel implementation, gated on the knob.**
3. **Unit test covering the new path.** Analytical tests are the highest-leverage
   thing you can write here — *every* analytical test added to the core has
   caught a real bug. Default to "add a test" over "audit the code".
4. **`pre-commit run --all` + fortitude** — pass before committing.
5. **`ctest` green on the GPU build** (never `-j N` — every worker shares one GPU;
   `ctest -R rakali` is the ~45 s regression subset).
6. **One capability per commit, one PR per capability.** Smaller PRs review faster
   and bisect cleanly. (Exception: several branches that each just add a namelist +
   IC + test → bundle into one `feat/<area>` PR.)

## Cross-cutting GPU rules (you *will* hit these)

The solver is GPU-native via `do concurrent` + OpenACC. State arrays stay
device-resident; only diagnostic output triggers D→H copies. These rules recur in
every guide — they're collected here so the guides can just point at them. Full
detail in [`FORTRAN_STYLE.md`](../../FORTRAN_STYLE.md) and the "Gotchas" section of
[`CLAUDE.md`](../../CLAUDE.md).

- **`do concurrent` loop order is `(k, j, i)` / `(j, i)` — j outer, i inner.**
  The reverse kills NVHPC GPU coalescing.
- **Flat-impl kernels take explicit-shape arrays** `arr(nx, ny, nz)`, never
  assumed-shape `arr(:,:,:)` — assumed-shape makes NVHPC emit a descriptor-walk
  memcpy per launch (millions of them in a real run).
- **Never `allocate` scratch inside a per-step kernel on `-stdpar=gpu`.** It stalls
  the device on entry/exit and bridges device inputs with implicit memcpys. Two
  sanctioned patterns: module-level allocatable + lazy `*_workspace_ensure` +
  `*_cleanup` (released in `solver_exit_data`), or a field on the state mapped once
  in `enter_data`. Per-column scratch uses fixed-size `NZ_STACK_MAX` stack arrays
  (a `local()` clause requires fixed size; a dummy-sized automatic array crashes
  with `CUDA_ERROR_ILLEGAL_ADDRESS`).
- **A loop-invariant `select case` goes *inside* one `do concurrent`, not split
  into one loop per case.** Splitting roughly doubles GPU launches (~1.8× slower);
  a uniform branch inside one loop is ~free.
- **A new allocatable on a mapped state type must be wired into the
  `enter_data` orchestration.** A silent omission is a 150–1500× memcpy
  explosion, not a crash — so it's easy to miss.
- **`pure` is the default** for new procedures. A `pure` helper called from inside
  a `do concurrent` GPU loop needs `!$acc routine seq`. NVHPC won't inline such
  helpers across module boundaries in hot kernels — duplicate the body as an
  `_impl` copy only when a profile justifies it (don't preemptively).
- **Bottom-up layer convention: `k=1` is the bed, `k=nz` is the surface.** Surface
  forcings land at `(:,:,nz)`; bed forcings at `(:,:,1)`. Any new per-layer kernel
  must follow it. Three regression gates guard it: `test_lock_exchange`,
  `test_river_plume`, `test_surface_flux_strat`.
- **Unstructured kernels are cell-centric double-visit**: each triangle accumulates
  only into *its own* arrays. Writing to a neighbour is a GPU race. Any new
  unstructured kernel must follow this.
- **Default-off = bit-identity** (rule 1 above, restated because it's the one most
  often forgotten): the default value must reproduce existing behaviour
  (`0.0` / `.false.` / a sentinel / a `1.0` multiplier), and the kernel must be
  gated so an un-set knob changes nothing.

## Working with an LLM coding assistant here

This repo is set up to be navigated by an agent. A few things make that work well:

- **The top-level [`CLAUDE.md`](../../CLAUDE.md) is the agent's primary brief** — it
  carries the architecture, the layer convention, the GPU rules, and a "Gotchas"
  section that is essentially a list of mistakes an LLM would otherwise make. Point
  the assistant at it first.
- **Give the assistant the relevant guide from this directory** as the task seam,
  plus the design contract README (`src/core/coastal/README.md`).
- **The assistant cannot see GPU correctness from the source.** NVHPC stdpar +
  OpenACC has sharp edges (descriptor-walk memcpys, illegal-address crashes from
  automatic arrays, silent approximate-math from intrinsic-shadowing locals). The
  ground truth is: build on the NVHPC GPU target and run `ctest` (no `-j`). Treat a
  green CPU/gfortran build as necessary-not-sufficient.
- **Prototype algorithm questions in Python first.** For a per-layer / per-cell
  numerical question, a tiny Python prototype on a representative setup answers it
  in seconds versus a ~9-minute GPU rebuild. Keep such Python prototypes and
  minimal `.F90` reproducers local — they're a development aid, not shipped
  repo artefacts.
- **Verify, then claim.** If `ctest` fails, say so with the output; don't report a
  refactor as done until it's green on the GPU build. The analytical tests are the
  arbiter.
