# The strict namelist schema (`rki_nml_schema`) — adding a group

The schema engine (`src/core/rki_nml_schema.F90`, type `nml_schema_t`) is the
strict, self-documenting configuration reader — every namelist group is on it
(the soft native `read(nml=)` readers are gone). One registration line per
knob carries the name, the pointer into the config, optional bounds/enum set,
units, and the doc string — and is the single list the parser, the parameter
dumps (`rakali_parameter_doc.{all,short}`), and `render_markdown`
all read from. Groups fail loudly (`error stop`, all errors collected,
did-you-mean suggestions) on unknown keys, type/range/enum violations, and
unsupported syntax. Destined for `pic` once trialed; the engine depends only
on `iso_fortran_env` + `pic_logger`.

Worked examples to copy: any `register_<group>` routine in
`build_rakali_schema` (`src/core/rki_config.F90`), e.g. `register_physics`
or `register_vcoord`.

## Contract (load-bearing)

1. **Defaults come from the config-type field initializers, never from the
   registration.** Constructors capture `default = current target value` at
   registration time — so register after the `config_t` exists with its
   defaults, before parsing. There is no `default=` argument on purpose.
2. **Group names are registered bare** (`"vcoord"`); files spell
   `&vcoord_nml`.
3. Validation that the engine can't express (exclusive bounds, cross-knob
   requirements involving non-config state) stays where it is
   (`validate_config` in `rki_config.F90`) — schema registration must never
   weaken an existing check.

## Recipe (one group per commit)

1. In `build_rakali_schema` (`rki_config.F90`): add a `register_<group>`
   routine registering **every** field of the group —
   `nml_real/int/logical/string/enum/real_array(name,
   cfg%field, doc [, units] [, min/max] [, allowed])` — and a
   `call register_<group>(cfg, schema)` line. String knobs with
   a fixed accepted set become `nml_enum` (canonicalizes case; extract the
   allowed list from the corresponding `parse_*` routine).
2. Run the **full** suite (gfortran + GPU builds): every repo `.nml` that
   reaches the driver is strictly validated for this group — fix stale
   keys in the files (or the registration, if the config type says the key
   is real).
3. Extend `tests/test_config_schema.F90` if the group introduces a new
   shape (first array knob, etc.).

Parse + dumps are already wired at both production entry points
(`app/main.F90` after the output-dir bootstrap; `rki_api.F90` handle
creation, dumps next to the input nml) — new groups need no wiring there.

## Regenerating the knob reference

`docs/generated_nml_knobs.md` is the full, code-derived table of every
namelist group + knob (name, default, units, description), rendered by the
`rakali_nml_doc` CLI tool (`app/nml_doc.F90`). It builds a default
`config_t`, registers the schema, and calls `render_markdown`. Because the
defaults come straight from the `config_t` field initialisers, the file
never drifts from the code as long as you regenerate it after touching a
group. From a build directory:

```
./rakali_nml_doc ../docs/generated_nml_knobs.md   # arg 1 = output path (default ./nml_knobs.md)
```

The tool is host-only (no NetCDF / GPU) and is built alongside `rakali`.
