# Rakali — Closure & Scheme Matrix

**Single source of truth for "which closure / scheme is enabled on which
backend."** When you add, port, or remove a closure, **update the relevant cell
in the same PR** — `tools/check_closure_matrix.py` (pre-commit/CI) fails the
build if a namelist closure knob has no row here, or a listed test is missing.

Cell convention: a cell holds the **namelist knob** that enables the scheme on
that backend (so the matrix doubles as "how to turn it on"), or `default` (active
unless another is chosen), or `—` (not available on that backend).

Backends: **coastal-S** = structured (Cartesian A-grid); **coastal-U** =
unstructured (triangular A-grid). Knobs live in `&nonhydrostatic_nml`
(3D/multilayer config), `&physics_nml` (forcing/drag), or `&tracer_nml`
(salinity/temperature IC/EOS/bounds).

Detailed prose + limitations: [`CAPABILITIES_AND_LIMITATIONS.md`](CAPABILITIES_AND_LIMITATIONS.md).
Physics/namelist reference: [`REFERENCE.md`](REFERENCE.md). Acronyms:
[`codebase/CONCEPTS.md`](codebase/CONCEPTS.md). When this file and the code
disagree, **the code is authority and this file is a bug** — fix it.

---

## Vertical mixing (3D / multilayer only — no-op at nz=1 and in 2D)

| Closure | coastal-S | coastal-U | Numerics | Primary test |
|---|---|---|---|---|
| PP81 (Pacanowski–Philander Ri) | `default` | `default` | explicit clamped flux | `test_lock_exchange`, `test_river_plume` |
| KPP (Large et al. 1994) | `vmix_use_kpp` | `vmix_use_kpp` | overlay on the PP81 interior | `test_coastal_kpp`, `test_coastal_kpp_unstr` |
| k-ε (two-equation, clean-room) | `vmix_use_keps` | `vmix_use_keps` | implicit Thomas | `test_keps`, `test_keps_unstr`, `test_kato_phillips`, `test_keps_kato_unstr` |

Notes: k-ε overrides KPP/PP81 when on; KPP overlays PP81. PP81/KPP apply an
**explicit** vertical-diffusion flux; the k-ε path is **implicit** (its ν_t can
reach `NUT_MAX`, which an explicit clamp would gut). Stratification-aware
stability functions are **built** — Galperin (1988) quasi-equilibrium and Canuto
(2001) Model A, selected by the `keps_stability` knob (default `constant` c_μ
baseline). The GLS-generic family (k-ω / MY2.5) is **not yet built**.

Under `time_scheme="semi_implicit"` (unstructured backend only) all three closures feed the
SI column tridiagonal through `nu_iface` and every vertical apply is implicit:
k-ε runs the same Patankar integrator once per SI step (`keps_si_picard`
optional 2-pass refinement); KPP's overlay is MAX'd into `nu_iface` **without**
the explicit-flux CFL clamp and its S/T diffusion + non-local transport run
implicitly (`rki_si_unstr_closures`; `test_si_closures_unstr`).

## Horizontal viscosity

| Closure | coastal-S | coastal-U | Test |
|---|---|---|---|
| Smagorinsky Laplacian | `smagorinsky_cs` | `smagorinsky_cs` | — (exercised by ML tests) |
| Biharmonic ∇⁴ (constant ν₄) | `biharmonic_nu4` | — | `test_biharmonic` |

Notes: explicit biharmonic is CFL-capped (`ν₄·dt·((π/dx)²+(π/dy)²)² ≤ 2`). The
Laplacian Smagorinsky has effective no-slip walls (ghost rows zero).

## Bottom drag

| Form | coastal-S | coastal-U | Notes |
|---|---|---|---|
| Manning quadratic (semi-implicit) | `default` | `default` | per-layer bed layer; the default |
| Linear (Rayleigh) | `use_linear_bottom_drag` + `bottom_drag_r` | — | structured only; the unstructured barotropic state carries no linear-drag rate |

## Wind-stress application (multilayer)

| Mode | coastal-S | coastal-U | Notes |
|---|---|---|---|
| Depth-uniform barotropic | `default` | `default` | added to the barotropic momentum |
| Surface-concentrated (top layer) | `wind_surface_layer` | — | generates Ekman/ML shear; required for the k-ε Kato-Phillips √t deepening |

## Baroclinic pressure gradient / PGF

| Variant | coastal-S | coastal-U | Knob |
|---|---|---|---|
| Jacobian (Shchepetkin–McWilliams) | `bpg_method="jacobian"` | `bpg_method="jacobian"` | the default |
| Sigma | `bpg_method` (other) | `bpg_method` (other) | |

## Vertical coordinates (dispatch via `vcoord_type`)

| Coord | coastal-S | coastal-U |
|---|---|---|
| Sigma | `default` | `default` |
| Z-sigma hybrid | `"zsigma"` | `"zsigma"` |
| Z-star (lite) | `"zstar"` | `"zstar"` |
| Z-star/sigma hybrid | `"zstar_sigma"` | `"zstar_sigma"` |
| Z-star full (per-column) | `"zstar_full"` | `"zstar_full"` |
| S-coord (ROMS) | declared, not implemented | — |

Vcoord parse/target tests: `test_vcoord`, `test_vcoord_target`. Solver-level
conservation under z*: `test_zstar_full_shelf_break`,
`test_zstar_full_shelf_break_unstr`, `test_zstar_full_tracer_conservation`.

### ALE remap reconstruction (dispatch via `remap_method`)

The conservative integrate-and-redistribute is shared across all coords; only
the per-column reconstruction order changes. All methods conserve column
integrals to roundoff and are monotone (no new extrema).

| Method | `remap_method` | Order | Notes |
|---|---|---|---|
| Piecewise constant | `"pcm"` | 0th | donor cell, diffusive |
| Piecewise linear | `"plm"` | 1st | minmod-limited |
| Piecewise parabolic | `"ppm"` (`default`) | 2nd | Colella & Woodward 1984; uniform-grid `(7/12,−1/12)` edge estimate |
| PPM, non-uniform H4 edges | `"ppm_h4"` | 2nd (parabola) + 4th edges | thickness-weighted (White & Adcroft 2008) edge estimate; cuts spurious diapycnal mixing (Ilicak 2012). Reuses the PPM limiter + parabola + redistribute verbatim. |
| Piecewise quartic (PQM) | `"pqm"` | 4th–5th | White & Adcroft (2008) PQM_IH4IH3: implicit-h4 edge values + implicit-h3 edge slopes → degree-4 reconstruction + W&A monotonicity limiter; `N<5` falls back to PPM. Opt-in (PPM stays default). |

The default is `remap_method = "ppm"` — stated as `"ppm"`. Remap tests:
`test_remap`, `test_remap_pqm`.

## Tracer advection & reconstruction

Horizontal tracer advection selects between first-order upwind and a high-order
reconstruction via `use_ppm_tracer` (default `.true.`). The scheme differs by
backend:

| Backend | Horizontal tracer advect | ALE vertical remap |
|---|---|---|
| coastal-S | upwind **or** true C-W PPM via `use_ppm_tracer` (default `.true.`); optional WENO face reconstruction via `tracer_recon` (ppm/pqm/weno5/weno7/weno9) | PCM/PLM/PPM/PPM_H4/PQM via `remap_method` (default `"ppm"`) |
| coastal-U | upwind **or** MUSCL piecewise-linear via the same `use_ppm_tracer` (default `.true.`) | PCM/PLM/PPM via the same `remap_method` (default `"ppm"`) |

**Honesty note.** `use_ppm_tracer` overstates the unstructured scheme: on
triangles it selects **MUSCL (piecewise-linear, Barth-Jespersen)**, not
parabolic PPM. The config docstring is honest ("PPM (structured) or MUSCL
(unstructured)"); the knob name is not — there is no parabolic tracer advection
on triangles. Structured face-reconstruction tests: `test_tracer_recon`.

## Equation of state

| EOS | coastal-S | coastal-U | Knob |
|---|---|---|---|
| Linear two-tracer | `default` | `default` | `&tracer_nml eos_form="linear"` |
| Wright (1997) nonlinear | `eos_form="wright"` | `eos_form="wright"` | shared `rki_eos`, evaluated at p = 0 (surface/potential density) |
| Roquet et al. (2015) SpV (TEOS-10-class) | `eos_form="roquet_spv"` | `eos_form="roquet_spv"` | shared `rki_eos`, p = 0 |

The nonlinear coastal EOS routes density through the same shared `rki_eos`
dispatch (`ml_parse_eos_form` → `use_shared` + `EOS_VARIANT_*`). Default
`"linear"` ⇒ bit-identical. EOS test: `test_ml_eos`.

## Boundary tracer fluxes

| Flux | coastal-S | coastal-U | Knob | Test |
|---|---|---|---|---|
| Surface heat / E−P (top layer `k=nz`) | scalar / field | scalar / field | `heat_flux`, `evap_minus_precip` (programmatic setters) | `test_surface_flux_strat` |

---

## Tunable parameters (key dials per scheme)

The *curated* dial set — defaults are the namelist defaults (`0` / `.false.` =
off/inert). The exhaustive per-field list (units, CF metadata) is the `config_t`
docstrings (FORD) + [`REFERENCE.md`](REFERENCE.md).

**PP81 vertical mixing** — `&nonhydrostatic_nml`
| Knob | Default | Meaning |
|---|---|---|
| `vmix_nu_max` | `0.01` | Max turbulent viscosity (m²/s). |
| `vmix_nu_bg` | `1e-4` | Background viscosity floor (m²/s). |
| `vmix_ri_alpha` | `5.0` | Ri scaling in `ν = ν_max/(1+α·Ri)² + ν_bg`. |

**KPP** — `&nonhydrostatic_nml` (overlay on PP81)
| Knob | Default | Meaning |
|---|---|---|
| `vmix_use_kpp` | `.false.` | Enable the KPP boundary-layer overlay. |
| `kpp_ri_crit` | `0.3` | Critical bulk Richardson number for the BL-depth sweep. |
| `kpp_cs_nonlocal` | `6.3` | Non-local (counter-gradient) transport coefficient `C_s`. |
| `kpp_c_vt2` | `1.8` | Unresolved-shear `V_t²` coefficient. |

**k-ε** — `&nonhydrostatic_nml` (both backends; overrides KPP/PP81)
| Knob | Default | Meaning |
|---|---|---|
| `vmix_use_keps` | `.false.` | Enable the two-equation k-ε closure. |
| `keps_stability` | `"constant"` | Stability-function scheme: `"constant"` (c_μ baseline), `"galperin"` (Galperin 1988 quasi-equilibrium), or `"canuto"` (Canuto 2001 Model A). |
| `keps_si_picard` | `1` | Picard passes of the k↔ε coupling per SEMI-IMPLICIT step (SI + k-ε only; explicit path ignores it). |

The k-ε model **constants** (`CMU=0.09`, `CE1=1.44`, `CE2=1.92`, `CE3=−0.4`,
`SIG_K=1.0`, `SIG_E=1.3`, `PR_T=1.0`, `KAPPA_EFF=0.4327`, roughness `Z0=0.01 m`,
`NUT_MAX=10`) are compile-time parameters in `rki_ml_keps`, **not yet namelist-
exposed** — the runtime dials are the on/off switch and `keps_stability`.

**Horizontal viscosity** — `&nonhydrostatic_nml`
| Knob | Default | Meaning |
|---|---|---|
| `smagorinsky_cs` | `0.0` | Smagorinsky coefficient `C_s` (0 = off; typical 0.1–0.3). |
| `biharmonic_nu4` | `0.0` | Biharmonic ∇⁴ viscosity (m⁴/s; 0 = off). Explicit-CFL capped: `ν₄·dt·((π/dx)²+(π/dy)²)² ≤ 2`. |

**Bottom drag** — `&physics_nml`
| Knob | Default | Meaning |
|---|---|---|
| `manning_n` | `0.0` | Manning roughness (s·m⁻¹ᐟ³); quadratic `C_d = g·n²/H^⅓`. Typical 0.02–0.04. |
| `use_linear_bottom_drag` | `.false.` | Switch to linear (Rayleigh) drag instead of Manning (structured only). |
| `bottom_drag_r` | `0.0` | Linear drag rate (1/s); used when the linear switch is on. |

**Wind** — `&physics_nml`
| Knob | Default | Meaning |
|---|---|---|
| `wind_surface_layer` | `.false.` | Apply wind stress to the top layer (shear) vs the depth-uniform barotropic mode (structured only). |

**Baroclinic pressure gradient** — `&nonhydrostatic_nml`
| Knob | Default | Meaning |
|---|---|---|
| `bpg_method` | `"jacobian"` | `"jacobian"` (Shchepetkin–McWilliams) or the sigma variant. |

**Tracer advection / EOS** — `&tracer_nml`
| Knob | Default | Meaning |
|---|---|---|
| `use_ppm_tracer` | `.true.` | High-order horizontal tracer advection (PPM structured / MUSCL unstructured) vs first-order upwind. |
| `remap_method` | `"ppm"` | Vertical ALE remap order (`pcm`/`plm`/`ppm`/`ppm_h4`/`pqm`). |
| `eos_form` | `"linear"` | Coastal EOS: `linear` (legacy default) / `wright` / `roquet_spv` (shared nonlinear). |

---

## Keeping this current

1. **Per-capability PR step.** The CLAUDE.md "adding a new capability" workflow
   includes: *namelist knob → kernel → unit test → **closure-matrix row** →
   pre-commit + fortitude → ctest → commit.* One row per scheme is cheap.
2. **Automated gate.** `tools/check_closure_matrix.py` parses the closure knobs
   from `src/core/rki_config.F90`, asserts each appears somewhere in this file,
   and checks that every test named here exists in `tests/CMakeLists.txt`. Wired
   into pre-commit + CI so drift fails the build instead of rotting silently.
