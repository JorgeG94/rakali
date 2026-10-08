# How to add a new passive tracer

Scope: the multilayer path (`nz_ml > 1`), both grid backends.

This is the easiest extension in the codebase, by design. The multilayer
tracer pipeline iterates a *registry*: append a tracer with
`register_passive_tracer(...)` and every transport kernel — advection, ALE remap,
vertical exchange, vertical/horizontal diffusion, conservative clamp, RK2 averaging
— picks it up automatically. You write zero kernel code.

## The seam

`multilayer_state_t` (structured) / `unstr_multilayer_state_t` (unstructured) own
`tracers(:)`, an allocatable array of `tracer_t` / `unstr_tracer_t`. Salinity and
temperature are the first two entries (`idx_salinity`, `idx_temperature`); a passive
tracer appends after them. Every passive-transport kernel loops
`do it = 1, size(tracers)`, so a new entry "just works". What stays S/T-specific by
design: the EOS and surface fluxes (a passive tracer has `eos_coeff = 0`, so it has
zero density feedback by construction).

The two backends differ only in array axis order (`tracer_t%hTr(nx,ny,nz)` vs
`unstr_tracer_t%hTr(nz,ncells)`) and in seeding (below).

## Recipe (ordered)

1. **Set up multilayer before `state%init`**: `use_multilayer = .true.`,
   `nz_ml > 1`. (At `nz_ml = 1` only transport generalises, not surface flux — see
   gotchas. Use `nz_ml > 1` for a fully-pipelined tracer.)
2. **After init, before the solver's `enter_data`**, call the registry method:
   ```fortran
   call state%multilayer%register_passive_tracer(grid_or_mesh, "dye", &
                                                  tr_init, tr_min, tr_max, &
                                                  kappa_bg, hdiff_kappa)
   ```
   It deep-copies the existing S/T entries (reallocating their `hTr`/`hTr0`),
   appends a new slot, `init`s it, copies your scalars, and forces
   `eos_coeff = 0` / `eos_ref = 0` (that's what "passive" means). `idx_salinity` /
   `idx_temperature` are unchanged.
3. **Capture the index**: `idx = size(state%multilayer%tracers)` — it's always the
   last slot.
4. **Seed `hTr`** — this differs by backend:
   - **Unstructured: automatic.** `solver_unstr_init_multilayer` loops the registry
     and seeds every non-S/T tracer from its `tr_init`. You do nothing.
   - **Structured: manual.** The init path seeds only S/T. After layer thicknesses
     exist, set it yourself:
     `state%multilayer%tracers(idx)%hTr = state%multilayer%h_layer * tr_init`
     (this is exactly what the tracer tests do).
5. **Done.** No kernel edits, no solver edits, no output edits — the registry loops
   carry it through transport, and the output writer loops the registry and writes
   `name`/`units`/`standard_name` verbatim.

## Why the order matters

`register_passive_tracer` grows `tracers(:)` via `move_alloc` + reallocation. The
solver's `enter_data` snapshots the array and maps each `hTr`/`hTr0` to the device
in a two-step map (`copyin` the array descriptor, then a per-element loop mapping
each `hTr`). Register *after* `enter_data` and the new slot's `hTr` is never mapped
→ `CUDA_ERROR_ILLEGAL_ADDRESS` on the first kernel.

## What you must NOT do

- Don't write `eos_coeff` — it must stay `0`, or the EOS would read density off an
  unbacked tracer.
- Don't assume `idx_salinity` / `idx_temperature` shift — they don't.
- Don't register mid-run after device mapping.
- Don't rely on automatic seeding on the **structured** backend.

## Test — the pattern to copy

The canonical pattern is a **tracking** test, not a sum check: register a "dye" with
salinity's bounds and zero density coupling, set `beta_S = 0` so salinity is
effectively passive too, run a short seiche, and assert
`maxval(|hTr_dye - hTr_S|) < 1e-9`. The dye must track salinity bit-for-bit through
advection + remap + vertical exchange + clamp + RK2. A sum check would falsely pass
if a kernel silently skipped the tracer; the tracking check won't.

- `tests/test_zstar_full_tracer_conservation.F90` — `passive_tracer_tracks_salinity`
  (nz>1) and `nz1_passive_tracer_tracks_salinity` (the nz=1 fused fast path).
- `tests/test_zstar_full_shelf_break_unstr.F90` —
  `unstr_passive_tracer_tracks_salinity` (registers before `enter_data`, relies on
  automatic unstr seeding).
- `tests/test_sediment.F90` — a real consumer: registers "sediment" as a passive
  tracer, then layers settling/erosion kernels on top (the non-generalisable physics
  stays *outside* the registry).

## Gotchas

- **`nz_ml = 1` fused fast path carries S/T only for surface flux.** Extra tracers
  get transport + RK2 + clamp on the fused path, but the surface-flux + S/T finalise
  steps are fused S/T-only. Use `nz_ml > 1` for a fully-pipelined passive tracer.
- **`enter_data` two-step map** is why registration must precede it (above).
- **Outer-shim + flat-impl device indirection.** A kernel can't dereference
  `state%multilayer%tracers(it)%hTr(...)` from inside a device loop (the
  array-of-derived-types descriptor isn't reachable on the device). The pattern: an
  outer shim dereferences the registry on the *host* and hands a flat array to a
  `pure *_impl` worker whose dummy is a plain `arr(nz,n)` / `arr(nx,ny,nz)`. Every
  registry kernel already does this — you only need to know it if you write a *new*
  registry-iterating kernel (you don't, for a plain passive tracer).

## Source pointers

- `src/core/structured/rki_tracer.F90`, `src/core/unstructured/rki_tracer_unstr.F90`
  — the tracer types + the pipeline opt-out flags (`do_horizontal_advection`,
  `do_vertical_diffusion`, …).
- `src/core/structured/rki_multilayer_state.F90`,
  `src/core/unstructured/rki_multilayer_state_unstr.F90` —
  `register_passive_tracer`, the registry, the index handles.
- `src/core/coastal/solver/rki_solver.F90`,
  `src/core/coastal/solver/rki_solver_unstr.F90` — the registry-iterating transport
  loops (search `do it = 1`).
- The Multilayer section of the top-level [`CLAUDE.md`](../../CLAUDE.md) describes
  the registry contract in prose.
