# How to add a new vertical coordinate

Adds a new `VCOORD_*` family to the coastal A-grid path (structured +
unstructured), which share one enum, one string parser, one coordinate type, and
one remap driver:

- `src/core/rki_constants.F90` — the `VCOORD_*` enum.
- `src/ALE/rki_vcoord.F90` — `vcoord_t`, `parse_vcoord_type`, the per-column
  target-thickness generator `vcoord_target_dz_column`, and the ZSTAR_FULL
  builders.
- `src/ALE/rki_remap_column.F90` — the conservative remap driver.

## The seam — key insight

The remap *reconstruction* (`remap_column` in `rki_remap_column.F90`) does **not**
switch on the coordinate — it switches only on the reconstruction method
(PCM/PLM/PPM/PPM_H4/PQM). A coordinate is defined *entirely* by its
**target-thickness generator**: the routine that, given total depth `H` (and
`η`), produces the target layer thicknesses `dz`/`target_h` that the conservative
remap maps onto. That's where the `select case (coord_type)` lives:

- `vcoord_target_dz_column(coord_type, …, dz)` in `rki_vcoord.F90` — a
  `!$acc routine seq` with a `select case`, output ROMS-ordered (`k=1` bed →
  `k=nz` surface), `sum(dz) == H`.

The per-column-reference coordinates (`VCOORD_ZSTAR_FULL`) are special-cased
*outside* the generic generator with their own builder
(`zstar_full_build_column` / `ml_build_vcoord_zstar_full{,_unstr}`) run once at
init from local bathymetry.

## Recipe (ordered) — add `VCOORD_FOO`

1. **Enum** — add `integer, parameter :: VCOORD_FOO = <n>` in
   `src/core/rki_constants.F90` and to its `public ::` list.
2. **Parser** — add a `case ("foo", "FOO", …)` → `VCOORD_FOO` in
   `parse_vcoord_type` (`src/ALE/rki_vcoord.F90`). The coastal default is
   `VCOORD_SIGMA`.
3. **`vcoord_t`** (if FOO needs new params): add knob fields to the type; handle
   template allocation in `vcoord_init` (allocate `z_ref` if FOO uses z-levels);
   update `vcoord_needs_remap` if FOO should *not* remap (currently anything
   `/= VCOORD_SIGMA` remaps, so pure-sigma-like coords skip the remap).
4. **Target generator** — add `case (VCOORD_FOO)` to `vcoord_target_dz_column`,
   producing ROMS-ordered `dz` summing to `H`. If FOO needs per-column reference
   tables (like ZSTAR_FULL), mirror the `zstar_full_build_column` + step-time
   `_zstar_full` + `ml_build_vcoord_zstar_full{,_unstr}` + state-allocation
   pattern instead of the generic generator.
5. **Config + namelist** — `vcoord_type` already exists in `&vcoord_nml`; add any
   new knob fields per [add_namelist_knob.md](add_namelist_knob.md).
6. **Driver wiring** — `src/driver/rki_driver.F90` already routes
   `parse_vcoord_type(cfg%vcoord_type)`; copy any new knob from `cfg%` onto
   `state%multilayer%vcoord%…`. If FOO needs a per-column build at setup, call the
   builder there (also built in `rki_init.F90`).
7. **State allocation** — if FOO needs a per-column z-ref table, mirror the
   `vcoord_z_ref_col` allocation (currently gated on `== VCOORD_ZSTAR_FULL`) in
   `rki_multilayer_state.F90` and the unstructured twin.
8. **Docs + tests** — add a row to `docs/CLOSURE_MATRIX.md` (vcoord table) and the
   tests below.

## Tests — the pattern to copy

- `tests/test_vcoord.F90` — string/alias → enum. **Add your aliases here.**
- `tests/test_vcoord_target.F90` — drives `vcoord_target_dz_column` with synthetic
  `dsig`/`z_ref`/`H`; asserts per-layer values + `sum(dz) == H`. **The canonical
  unit test for a new coordinate.**
- `tests/test_vcoord_zstar_full.F90` — the per-column-reference builder.
- `tests/test_zstar_ssh_tracking.F90`, `tests/test_zstar_full_shelf_break{,_unstr}.F90`,
  `tests/test_zstar_full_tracer_conservation.F90` — solver-level mass/salt/heat
  conservation + layer-ratio preservation.

Invariants to assert: (a) `sum_k dz == H` (or `H + η`) per column; (b) layers stay
≥ floor; (c) at solver level, tracer/mass conservation across steps. **A
deep-shallow bathymetry + flow + active-clamp test is required** to catch
over-allocation bugs — flat-bottom / static-stratification tests are blind to them.

## Gotchas

- **Conservation is load-bearing.** Every branch must guarantee `sum(dz) == H`
  (or `H + η`); the conservative remap depends on it. Decide where any clip deficit
  goes (ZSIGMA deposits it in the bed layer; ZSTAR_FULL trims the surface).
- **Bottom-up ROMS ordering** (`k=1` bed). The builders construct z-levels
  top-down then *reverse* into ROMS order — don't flip.
- **Vanishing layers.** ZSTAR_FULL clips bed-side layers to `zstar_h_min`; operators
  that divide by `h_layer` gate on `VANISHING_LAYER_TOL`. A negative interface NaNs
  the CFL — defend degenerate columns (`h_bed <= 0` → all-zero output).
- **Ghost fill.** Per-column z-ref must be built in ghost rows too. Leaving
  `h_layer = 0` at ghosts makes the EOS fall back to ρ₀ → a wall-adjacent density
  jump → blowup.
- **NVHPC cross-module-call codegen.** The remap kernel uses
  `!$acc parallel loop` (not `do concurrent`) precisely because it calls
  `vcoord_target_dz_column` from another module and stdpar GPU won't lower that from
  a DC even with `!$acc routine seq`. A new coordinate added to that generator
  inherits the working path.
- **`needs_remap()` gate.** Pure SIGMA skips the remap entirely; everything else
  routes through it. Choose deliberately for FOO.

## Source pointers

- `src/core/rki_constants.F90` — `VCOORD_*` enum.
- `src/ALE/rki_vcoord.F90` — `vcoord_t`, `parse_vcoord_type`,
  `vcoord_target_dz_column`, the ZSTAR_FULL builders.
- `src/ALE/rki_remap_column.F90` — the remap driver.
- The "Vertical Coordinates" section of the top-level [`CLAUDE.md`](../../CLAUDE.md).
