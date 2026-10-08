# How to add a new namelist knob / config option

A knob is the entry point for almost every capability (workflow step 1). Get the
plumbing right once and the rest of the guides reuse it.

All config lives in `src/core/rki_config.F90`, on the single public **`config_t`**:
grid, time, physics, tracer, vcoord, NH, output. Each field is declared
*with its default inline* (`real(wp) :: smagorinsky_cs = 0.0_wp`).

## The seam

Parsing goes through the strict schema engine (`rki_nml_schema`, built by
`build_rakali_schema` in `rki_config.F90`): each namelist group has a
`register_<group>` routine that binds every key to its `config_t` field by
pointer (`nml_real` / `nml_int` / `nml_logical` / `nml_enum`). The default lives
on the *type*, so an absent namelist block (or key) is a no-op — that's what
gives you bit-identity for free.

```
type field (default)  ->  g%add(nml_*("key", ptr, doc)) in register_<group>
   ->  strict parse fills the field  ->  propagate cfg -> state
   ->  consume in kernel (gated)
```

## Recipe — a new knob (e.g. a new `&physics_nml` key `foo`)

1. **Field + inline default** on `config_t` in `rki_config.F90`:
   `real(wp) :: foo = 0.0_wp  !! <FORD doc>`. The default must reproduce
   current behaviour.
2. **Register it.** In the group's `register_<group>` routine (here
   `register_physics`), bind a pointer and add the entry:
   `pr => cfg%foo` then `call g%add(nml_real("foo", pr, "<one-line doc>"))`.
   Enum-valued keys use `nml_enum` with an `allowed=` list (see
   `vcoord_type` in `register_vcoord` for the pattern); ranged ints take
   `min=`/`max=`. The strict parser rejects unknown keys, so a typo in a user
   namelist fails loud instead of being silently dropped.
3. **Validate (optional).** Range/consistency checks go in `validate_config`.
   Set `has_error` / log and let the routine `error stop`.
4. **Propagate in *two* places** — `src/core/structured/rki_state.F90` (state-init
   path) **and** `src/driver/rki_driver.F90` (driver path). Miss one and that path
   silently keeps the default. This is the most common knob bug.
5. **Consume in the kernel, gated on the knob.**
6. **FORD `!!` docstring on the field** — that *is* the documentation. Don't write
   parallel markdown describing the knob.

Adding a whole new group is the same shape one level up: a new
`register_<group>` routine + a `call register_<group>(cfg, schema)` line in
`build_rakali_schema` — see
[migrate_namelist_to_schema.md](migrate_namelist_to_schema.md).

## Default-off / bit-identity

The default must be the value that reproduces existing behaviour. Established
conventions in the tree: `0.0` (e.g. `smagorinsky_cs`), `.false.` (every
`vmix_use_*`), a `-1.0` sentinel ("unset, derive it"), a `1.0` multiplier
(`bed_factor`). The kernel is gated so that the default changes nothing — verify by
running `ctest` and confirming existing tests stay byte-for-byte.

## The CLOSURE_MATRIX drift hook

`tools/check_closure_matrix.py` (pre-commit + CI) fires when you touch
`rki_config.F90`, `docs/CLOSURE_MATRIX.md`, or `tests/CMakeLists.txt`. It enforces:

1. **Knob coverage** — every auto-discovered `logical :: vmix_use_* =` toggle, plus
   a curated list (`smagorinsky_cs`, `biharmonic_nu4`, `use_linear_bottom_drag`,
   `wind_surface_layer`, `bpg_method`), must appear as a string in
   `docs/CLOSURE_MATRIX.md`. Curated names must still exist in the config (rename
   guard).
2. **Test coverage** — every `test_*` token cited in the matrix must be registered
   in `tests/CMakeLists.txt`.

So: if you add a new `vmix_use_*` toggle (or another closure selector you add to
the curated list), you **must** add a matching row/cell to `docs/CLOSURE_MATRIX.md`
in the same change, or the hook fails the commit. A plain numeric knob that is not
a closure selector isn't enforced by the hook — but the matrix is still the
single source of truth for "how do I turn this on", so keep it current by
convention.

## Tests

There isn't a dedicated "config" test per capability — the knob is exercised by the
feature's own analytical test (the kernel gated on it). What you *do* verify:
existing tests stay bit-identical with the default, and the feature test flips the
knob on. If you add a sub-namelist group, `tests/` parsing coverage lives alongside
the feature.

## Gotchas

- **Knobs propagate in two places** (`rki_state.F90` + `rki_driver.F90`) — miss
  one and that path silently keeps the default.
- **The closure-matrix hook can block a commit that compiles fine** — update
  `docs/CLOSURE_MATRIX.md` in the same change for any `vmix_use_*` / curated knob.
- **Knobs are deliberately not in the Python/C FFI.** The API surface
  (`src/api/`) is field-pointer getters + bathymetry/IC setters; knobs reach the
  solver only through the namelist file.

## Source pointers

- `src/core/rki_config.F90` — `config_t`, `build_rakali_schema`, the
  `register_<group>` routines, `validate_config`.
- `src/core/rki_nml_schema.F90` — the strict schema engine (`nml_schema_t`).
- `src/driver/rki_driver.F90` — driver propagation.
- `src/core/structured/rki_state.F90` — state propagation.
- `tools/check_closure_matrix.py`, `docs/CLOSURE_MATRIX.md`.
