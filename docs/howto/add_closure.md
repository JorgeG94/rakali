# How to add a new turbulence closure / parameterization

Covers vertical mixing (PP81, KPP, k-ε) and a note on lateral closures
(Smagorinsky, Leith). Modules live under `src/parameterizations/`.

Worked references: KPP and k-ε were both added relatively recently and are the
templates to copy.

## The seam — vertical mixing dispatch

The solver picks the closure each step with a single `if / else if / else`
keyed on two flags on `multilayer_state_t`:

```fortran
if (state%multilayer%vmix_use_keps) then       ! k-ε takes precedence
   ... ml_keps_friction_velocities / ml_keps_update / ml_keps_apply_momentum
       / per-tracer ml_keps_apply_diffusion_tracer
else if (state%multilayer%vmix_use_kpp) then    ! KPP overlay on PP81
   ... ml_kpp_boundary_layer / ml_kpp_viscosity / ml_kpp_diffusion_tracers
       (extra passive tracers fall back to PP81)
else                                            ! default
   ... ml_vertical_viscosity (PP81) / per-tracer ml_vertical_diffusion_tracer
end if
```

- Structured: in `solver_step` (the `if (nz > 1)` vmix block) of
  `src/core/coastal/solver/rki_solver.F90`.
- Unstructured: the **same three-way switch appears twice** (once per RK2 stage) in
  `src/core/coastal/solver/rki_solver_unstr.F90`.

### Interface contract a vertical closure must satisfy

- **Inputs** (read from state): per-layer `h_layer, hu_layer, hv_layer, rho_layer`
  (shear + stratification/Ri), `rho_0`; tunables `vmix_nu_max/nu_bg/ri_alpha`;
  surface forcing; for KPP/k-ε the friction velocity `u_*` and surface buoyancy
  flux (KPP computes them in its BL sweep; k-ε via `ml_keps_friction_velocities`).
- **Outputs**: `nu_t` / `kappa_t` at **interfaces**, applied as a diffusive
  momentum/tracer flux. Coastal PP81/KPP apply an **explicit** CFL-clamped flux
  in-kernel; k-ε is **implicit** (Thomas sweep) because its `nu_t` can reach values
  an explicit clamp would gut.
- **Persistent state** on `multilayer_state_t`: the flags `vmix_use_kpp` /
  `vmix_use_keps`; KPP 2D workspace (`kpp_bl_depth` persists across steps, seeds
  V_t²); k-ε interface-shaped (`nz+1`) prognostics `keps_tke/eps/nut/kt`
  (`keps_nut` persists — the implicit step reads the previous value) + 2D
  `keps_ustar_b/s`.

## Recipe (ordered) — a new vertical closure `<name>`

1. **Config knob** — `logical :: vmix_use_<name> = .false.` per
   [add_namelist_knob.md](add_namelist_knob.md) (and copy onto state).
2. **State flag + fields** — add the flag + any interface-shaped (`nz+1`) persistent
   fields to `multilayer_state_t`.
3. **Allocate** in `multilayer_state_init`, gated `if (this%vmix_use_<name>)`, with
   positive-floor seeds if step-1 division needs them; deallocate in `_destroy`.
4. **GPU map** — gated `!$acc enter data` (associate-leaf, one-deep attach) in the
   solver + the matching `exit data`.
5. **Kernel module** — `src/parameterizations/vertical/structured/rki_ml_<name>.F90`:
   `pure` subroutines, explicit-shape args `arr(nx,ny,nz)`, `do concurrent(j, i)`
   over columns (j outer, i inner), bottom-up `k`. Register in `src/CMakeLists.txt`.
6. **Wire dispatch** — add an `else if (state%multilayer%vmix_use_<name>)` branch in
   the solver vmix block; import the new publics. Decide precedence vs k-ε/KPP/PP81.
7. **Default-off** — flag `.false.`, fields unallocated → bit-identical.
8. **Both backends** — repeat for unstructured: fields on
   `unstr_multilayer_state_t` (**allocated unconditionally** there, because the
   API-driven path can flip the flag after init), map in
   `unstr_multilayer_enter_data`, add the branch **twice** in `rki_solver_unstr.F90`,
   kernel in `unstructured/rki_ml_<name>_unstr.F90` following the cell-centric
   double-visit rule (no neighbour writes).
9. **Test + matrix row** — add `tests/test_<name>.F90` to `tests/CMakeLists.txt`
   and a row to `docs/CLOSURE_MATRIX.md` (the drift hook fails the build if a
   `vmix_use_*` knob has no matrix row or a cited test is unregistered).

## Lateral closures (briefly)

Under `src/parameterizations/lateral/`:

- `rki_ml_horizontal_viscosity.F90`: Smagorinsky + biharmonic,
  dispatched by **threshold, not a switch** — `if (smagorinsky_cs > 0.0)` and
  `if (biharmonic_nu4 > 0.0)`; both can be active. Uses the lazy-allocate workspace
  pattern (`*_workspace_ensure` / `*_cleanup`).

A new lateral closure follows the same shape: add the threshold/knob, the kernel,
and gate it.

## Tests — the pattern to copy

- `tests/test_coastal_kpp.F90` — analytic shape-function checks; **`no_forcing_equals_pp81`**
  (the default-off / off-forcing bit-identity invariant — the canonical proof that a
  new overlay collapses to the base closure when inactive); BL-depth in a stratified
  vs unstratified column. Mirror: `tests/test_coastal_kpp_unstr.F90`.
- `tests/test_keps.F90` — `homogeneous_decay` against the analytic power law,
  positivity, momentum conservation, log-layer wall function, tracer mixing.
  `tests/test_kato_phillips.F90` — solver-level √t mixed-layer entrainment.
- Layer-convention gates every closure must keep green: `test_lock_exchange`
  (dense at `k=1`), `test_river_plume` (buoyant at `k=nz`),
  `test_surface_flux_strat`.

The k-ε implementation is clean-room from the published equations (Launder &
Spalding 1974; constants per Rodi 1987 / Warner et al. 2005) — **not** derived from
GOTM or any GPL source. Keep new closures the same way: cite the paper, not another
codebase. Validate against an analytic limit (decay law, log layer) where one
exists.

## Gotchas

- **No local-allocate scratch in a per-step kernel on `-stdpar=gpu`.** Use either a
  module-level allocatable + lazy `*_workspace_ensure` + `*_cleanup` (released in
  `solver_exit_data`), or fields on the state mapped once in `enter_data`. Per-column
  scratch is fixed-size `NZ_STACK_MAX` (a `local()` clause requires fixed size; a
  dummy-sized automatic array crashes with `CUDA_ERROR_ILLEGAL_ADDRESS`).
- **`do concurrent(j, i)` ordering** (j outer) for coalescing; **explicit-shape**
  kernel args, never assumed-shape.
- **`!$acc routine seq`** on cross-module `pure` helpers called from the device loop;
  for hot inner loops NVHPC won't inline them — keep helpers in-module or duplicate
  as `_impl`.
- **Bottom-up `k`** (`k=1` bed, `k=nz` surface); interface `k` sits between layers
  `k` and `k+1`; surface forcing at the top interface.
- **Unconditional allocation on the unstructured side** — the API path can flip the
  closure flag after init, so conditional allocation would leave the fields unmapped.
- **Which path has what**: default = PP81, optionally KPP or k-ε (k-ε
  overrides). There is no standalone implicit vdiff solver on the explicit path
  (KPP is fused with the explicit flux apply; k-ε's apply is implicit). Don't
  assume a closure exists on a path without checking the dispatch.
- **The CLOSURE_MATRIX can lag the code.** Treat the dispatch in the solver as the
  source of truth for what's wired; update the matrix row when you add a closure.

## Source pointers

- `src/core/coastal/solver/rki_solver.F90` (vmix dispatch + GPU map),
  `src/core/coastal/solver/rki_solver_unstr.F90` (dispatch ×2).
- `src/core/structured/rki_multilayer_state.F90` (the state contract: flags + fields
  + allocation).
- `src/parameterizations/vertical/structured/{rki_ml_vertical_viscosity,rki_ml_kpp,rki_ml_keps}.F90`
  + unstructured twins; `rki_keps_stability.F90` (shared).
- `src/parameterizations/lateral/structured/rki_ml_horizontal_viscosity.F90`.
- `docs/CLOSURE_MATRIX.md` — which closure is on where + the knobs.
