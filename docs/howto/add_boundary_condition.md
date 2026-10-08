# How to add a new boundary condition

There are nine BC types: `BC_WALL=1`, `BC_OPEN=2` (Flather), `BC_TIDAL=3`,
`BC_NESTED=4`, `BC_INFLOW=5`, `BC_DISCHARGE=6`, `BC_CLAMPED=7`, `BC_SPONGE=8`,
`BC_CHAPMAN=9` (Orlanski). Adding one touches the structured kernel, the
unstructured flux kernel, and the config.

## The seam — module split + dispatch

- **Backend-agnostic types/constants**: `src/core/structured/rki_boundary_types.F90`
  defines `BC_*`, `bc_params_t` (per-edge tags, tidal params, per-BC scalars, the
  Chapman persistent `eta_old_*` arrays), `bc_type_from_string`,
  `bc_init_from_config`, `bc_init_chapman`.
- **Structured kernel**: `src/core/coastal/kernels/structured/barotropic/rki_boundary.F90`.
  `apply_boundary_conditions` builds a per-edge descriptor and runs one
  `select case (bc%<side>)` per edge (west/east/south/north), each calling a shared
  per-edge helper (`bc_wall_edge`, `bc_flather_edge`, `bc_inflow_edge`,
  `bc_disch_edge`, `bc_clamped_edge`, `bc_chapman_edge`). Each helper runs a
  `do concurrent(t, k)` over (tangential index, ghost layer) and writes **only ghost
  cells**.
- **Unstructured**: there is *no* `bc_params_t`; config rides in fixed-size arrays on
  the state (`bc_type(MAX_BC_GROUPS)`, `bc_*_ext`). The BC dispatch is **inlined into
  the flux kernel** `src/core/coastal/kernels/unstructured/barotropic/rki_kernel_flux_unstr.F90`:
  a boundary edge is `neighbour < 0`, the group → `bc_type`, and
  `compute_boundary_state_impl` does the `select case`. Chapman is dispatched one
  level up (it needs per-edge mesh access) then rewrites itself to `BC_TIDAL` so the
  Flather algebra runs.

Supported on unstructured: wall, open, inflow, discharge, outflow (unstr-only,
`BC_OUTFLOW=10`), tidal/Flather, clamped, Chapman. **Not ported**: `BC_NESTED`,
`BC_SPONGE`.

## Recipe (ordered) — add `BC_FOO`

**Structured + shared + config:**

1. `rki_boundary_types.F90` — add `integer, parameter :: BC_FOO = <n>`; add to the
   `public` list; add a `case ("foo")` in `bc_type_from_string`. If FOO needs config
   scalars, add fields to `bc_params_t` and copy them in `bc_init_from_config`. If it
   needs persistent ghost state, mirror the Chapman `eta_old_*` pattern (allocated
   unconditionally so it GPU-maps) + extend `bc_init_chapman`-style seeding.
2. `rki_boundary.F90` — add `BC_FOO` to the `use` + `public`; add a `case (BC_FOO)`
   to **all four** edge `select case` blocks; write a `pure subroutine bc_foo_edge`
   following the `do concurrent(t, k)` ghost-only pattern (model on `bc_wall_edge`).
3. `rki_config.F90` — add any new `&boundary_nml` keys (declare on `config_t`, add to
   the namelist statement, read-back) per [add_namelist_knob.md](add_namelist_knob.md);
   add `"foo"` to the `warn_unknown_bc` known list; validate if needed.
4. **GPU map** any new persistent ghost array in the structured solver's
   `acc enter data` / `solver_exit_data` (the `state%bc%eta_old_*` pattern).

**Unstructured:**

5. `rki_boundary_unstr.F90` — add `integer, parameter :: BC_FOO = <n>` (match the
   structured value; use ≥10 for an unstr-only type) + `public`.
6. `rki_kernel_flux_unstr.F90` — add `BC_FOO` to the `use`; add a `case (BC_FOO)`
   inside `compute_boundary_state_impl`'s `select case`. If FOO needs per-edge mesh
   or persistent state, dispatch it at the higher level like Chapman.
7. If it needs a per-group param or persistent state: add a fixed-size array on
   `unstr_state_t` (rides along via `map(to: state)`) and/or a persistent array on
   the unstructured barotropic sub-state with its `init`/`destroy`/`acc enter/exit
   data` (the `chapman_eta_old` pattern).

**Test:** add cases to `tests/test_boundary.F90` (structured) +
`tests/test_state_unstr.F90` (unstructured default-table coverage).

## Tests — the pattern to copy

`tests/test_boundary.F90` asserts each BC's ghost-cell effect on flat and sloped
bottoms: wall negates normal momentum + mirrors tangential + mirrors η; open is
zero-gradient; tidal prescribes `A·sin(2πt/T)` (+ multi-constituent superposition);
inflow prescribes `hu = h·v`; clamped sets η; **sponge falls through to wall** (its
west ghost `hu` is negated like wall); Chapman with `dt=0` extrapolates interior η,
with `dt>0` covers two edges. `tests/test_state_unstr.F90` asserts the all-walls
default (`all(bc_type == BC_WALL)`). Pattern: a tiny grid, set one edge to your BC,
step the ghost fill, assert the closed-form ghost values.

## Gotchas

- **Unstructured cell-centric double-visit**: write only to the current cell.
  Neighbour writes are a GPU race. Chapman's per-edge `chapman_eta_old(i,k)` write is
  race-free *because* each `(i,k)` is touched only by its own cell.
- **`BC_SPONGE` shares the wall case** in the structured kernel; the actual
  relaxation is a *separate* kernel (`apply_sponge`) the solver runs after the BC
  pass. Sponge is not ported to unstructured.
- **`BC_NESTED` is a no-op case** in the structured BC kernel — ghosts come from the
  nesting module, not the BC dispatch.
- **Chapman persistent state location differs**: structured keeps `eta_old_*` on
  `bc_params_t` (mapped unconditionally); unstructured keeps `chapman_eta_old` on the
  barotropic sub-state per-edge.
- **GPU mapping**: a new *allocatable* persistent ghost array must be added to the
  relevant `enter_data` / `exit_data`. A *fixed-size* array on a mapped struct (the
  unstr `bc_*_ext(MAX_BC_GROUPS)` pattern) rides along automatically.
- **Default falls through to wall** on both backends (an unhandled / group-≤0 edge is
  a wall) — keep that invariant for FOO.

## Source pointers

- `src/core/structured/rki_boundary_types.F90` — constants, `bc_params_t`,
  `bc_type_from_string`, `bc_init_from_config`, `bc_init_chapman`.
- `src/core/coastal/kernels/structured/barotropic/rki_boundary.F90` — the per-edge
  dispatch + helpers.
- `src/core/unstructured/rki_boundary_unstr.F90` — the unstructured constant mirror.
- `src/core/coastal/kernels/unstructured/barotropic/rki_kernel_flux_unstr.F90` —
  `compute_boundary_state_impl` + the higher-level Chapman dispatch.
- `src/core/rki_config.F90` — `&boundary_nml`, `warn_unknown_bc`.
- The "Boundary Conditions" section of the top-level [`CLAUDE.md`](../../CLAUDE.md).
