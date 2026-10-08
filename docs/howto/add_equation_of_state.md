# How to add a new equation of state

There are **two EOS code paths**, and they're at different maturity:

- **Hardwired linear** — `src/equation_of_state/structured/rki_ml_eos.F90` (+ the
  unstructured twin `.../unstructured/rki_ml_eos_unstr.F90`). The legacy default:
  `rho = rho_0 + beta_S*(S - S_ref) - alpha_T*(T - T_ref)`, no enum, no dispatch.
  It also owns the `rho_ref` builder (`ml_compute_reference_density`) the Jacobian
  BPG subtracts.
- **Shared multi-variant** — `src/equation_of_state/structured/rki_eos.F90`. A
  proper variant enum + `select case` dispatch (`EOS_VARIANT_LINEAR`,
  `EOS_VARIANT_WRIGHT_97`, `EOS_VARIANT_ROQUET_SPV`; `EOS_VARIANT_TEOS10`
  reserved). This is where a new nonlinear EOS goes.

The coastal multilayer path reaches the shared dispatch when
`&tracer_nml eos_form = "wright" | "roquet_spv"`: `ml_parse_eos_form(name,
use_shared, variant)` in `rki_ml_eos.F90` translates the string into
`(use_shared, EOS_VARIANT_*)` held on the multilayer state, and
`ml_equation_of_state_shared` runs the shared kernel. `eos_form = "linear"` (the
default) stays on the hardwired path ⇒ bit-identical.

## The seam

`eos_compute_arrays(eos, …)` in `rki_eos.F90` does `select case (eos%variant)`
over the `EOS_VARIANT_*` codes, each calling a `pure *_impl` that is a single
`do concurrent(k,j,i)` with a vanishing-layer fallback to `rho_0`. The variant
lives on `eos_t%variant` (default `EOS_VARIANT_LINEAR`). `eos_validate` is
fail-loud on an unimplemented variant (TEOS-10 aborts). `parse_eos_variant(name)`
maps a string to the code.

All shared variants here are evaluated at **p = 0** (surface / potential
density). A pressure-dependent (in-situ) evaluation would couple to the
pressure-gradient integration path, which the coastal BPG does not currently do —
that is a larger change than a new arm.

## Recipe — a new variant on the shared path

1. **Enum** — add `EOS_VARIANT_<NAME>` as a public `integer, parameter` in
   `rki_eos.F90`.
2. **Kernel impl** — add `pure subroutine eos_<name>_impl(...)` modelled on the
   linear impl: explicit-shape 3D arrays, one `do concurrent(k,j,i)`,
   vanishing-layer fallback to `rho_0`.
3. **Dispatch arm** — add `case (EOS_VARIANT_<NAME>)` inside the *existing*
   `select case (eos%variant)` (keep it one select-case — see gotchas).
4. **Validate + parse** — add the code to the accepted list in `eos_validate`,
   and a `case ("<name>", ...)` in `parse_eos_variant`.
5. **New coefficients** on `eos_t` if needed (scalar defaults; any allocatable
   scratch must be released in `eos_destroy`).
6. **Route the coastal knob** — add a `case ("<name>")` in `ml_parse_eos_form`
   (`rki_ml_eos.F90`) setting `use_shared = .true.` + `variant =
   EOS_VARIANT_<NAME>`, and extend the `eos_form` enum's `allowed=` list in
   `src/core/rki_config.F90` (`&tracer_nml`). See
   [add_namelist_knob.md](add_namelist_knob.md).
7. **Test** — copy `tests/test_coastal_shared_eos.F90`.

## Recipe — extending the hardwired linear path

The hardwired `rki_ml_eos.F90` / `rki_ml_eos_unstr.F90` path is for the legacy
linear default only. If a change belongs there (e.g. a new reference-density
form), keep the `select case`-free structure and update **all call sites**: the
two structured solver stages, `rki_init.F90`, `rki_api.F90`, and the two
unstructured shim calls. But any *new EOS variant* belongs on the shared path
above, not here — that is what the `eos_form` routing exists for.

A nonlinear coastal EOS must stay compatible with the
`ml_compute_reference_density` path (the Jacobian BPG subtracts a
horizontally-uniform `rho_ref(z)` so the integrand is the small anomaly
`ρ - rho_ref(z)`).

## Tests — the pattern to copy

Analytical point-value asserts, one per branch of the formula:

- `tests/test_ml_eos.F90` (hardwired coastal linear): reference-state,
  linear-in-S, linear-in-T sign, `alpha_T = 0` makes T passive, dry-cell
  fallback, plus the `rho_ref` builder invariants.
- `tests/test_coastal_shared_eos.F90`: the coastal path routed through the shared
  `rki_eos` dispatch — a known reference seawater value within a loose tolerance;
  monotonicity (warmer → lighter, saltier → denser); a **cabbeling** case (blend
  of two equal-density parcels is denser than the linear mean — the discriminator
  a linear EOS *fails*); the vanishing-layer guard (`h ≤ 0 → rho_0`). The harness
  runs `enter_data → compute → exit_data` so the GPU path is exercised.

## Gotchas

- **One `select case` inside one `do concurrent`** — splitting variant dispatch
  into one DC-launch per case is ~1.8× slower; a uniform branch is ~free.
- **`pure` + `!$acc routine seq`** for any helper called from the device loop;
  NVHPC won't inline cross-module `acc routine seq` helpers in hot kernels (the
  Wright coefficients are inlined as module `parameter`s for exactly this reason).
- **Reference-density consistency** — the coastal BPG
  (`ml_compute_reference_density`) subtracts a reference density; a new EOS must
  remain consistent with the builder it feeds.
- **Vanishing-layer fallback is mandatory** — every impl guards `h ≤ 0` (or
  `> DRY_TOLERANCE`) and falls back to `rho_0`/reference T,S (required for
  ZSTAR_FULL where bed-side layers pinch out).

## Source pointers

- `src/equation_of_state/structured/rki_eos.F90` — shared variant enum, impls,
  dispatch, `parse_eos_variant`, `eos_validate`.
- `src/equation_of_state/structured/rki_ml_eos.F90` +
  `.../unstructured/rki_ml_eos_unstr.F90` — hardwired coastal linear EOS,
  `ml_parse_eos_form` (the `eos_form` → variant router), and the `rho_ref`
  builder.
- `src/core/rki_config.F90` — the `&tracer_nml eos_form` knob (`linear` /
  `wright` / `roquet_spv`).
- The "Equation of state" row of [`docs/CLOSURE_MATRIX.md`](../CLOSURE_MATRIX.md).
