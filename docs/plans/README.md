# Capability plans

Executable plans for closing the coastal path's physics and validation gaps,
plus the analyses they came from. Written 2026-08-06.

> **Resuming work? Start with [`pickup.md`](pickup.md).** It carries the build
> environment (the system gfortran cannot compile this repo), the PR stack, the
> defects already fixed, the open findings with every candidate already
> eliminated, and the traps hit along the way.

## Read in this order

1. **`coastal_physics_gaps.md`** — the source-verified gap analysis. Ranked by
   *scientific consequence*: what would make a domain scientist reject this
   model for real work, and for which applications it is already credible.
   Start here; the three plans below are its answer.
2. **`plan_validation_campaign.md`** — 9 phases. Recommended first, because
   validation is the cheapest credibility per unit effort and every later
   capability needs the conservation diagnostic and convergence harness it
   builds. Highest-value single item is the well-balanced C-property
   assertion (~half a day, exact-zero oracle, both backends).
3. **`plan_wave_coupling.md`** — 6 phases. Architecture decision: do NOT port a
   spectral model; ingest an external wave field via NetCDF and build the
   per-cell forcing terms here.
4. **`plan_morphodynamics.md`** — 9 phases. Keystone is P3, Exner bed
   evolution — the first time the seabed `b` moves.

`schism_parity_report.md` is background: a feature comparison against SCHISM
(Apache-2.0, freely readable). **Revised 2026-08-07** — re-verified line-by-line
against both trees. The SCHISM-side claims (module paths, sizes, identifiers)
all held up; the Rakali side had been written against the pre-port-over repo and
carried dangling `src/core/{ocean,ice}/` references, a wrong shipped-tracer list,
and an over-generous reading of semi-implicit maturity. Those are fixed. It is
still *input*, not truth — verify against source before acting on it, and when a
doc and the code disagree, the code is authority.

## Scope of this repo

Coastal only. `sim_type` accepts `'coastal'` and nothing else
(`src/core/rki_config.F90:847`); there is no ocean dyn-core and no ice. The
target is **SCHISM physics parity on both time schemes** — explicit Godunov and
Casulli semi-implicit. Numerical-method divergence (Godunov vs semi-implicit
ELM, triangles vs mixed tri/quad, z\*/ALE vs LSC2) stays deliberate. Sea ice is
out of scope by decision.

## Conventions these plans assume

- **Every phase is one PR**: namelist knob defaulting OFF (bit-identity for
  existing configs), kernel, and a test — preferably analytical with a
  closed-form oracle.
- **All four paths**: {structured, unstructured} x {explicit, semi-implicit}.
  The Casulli semi-implicit mode runs on both backends at both 2D barotropic
  and 3D sigma multilayer; wave and bed forcing enter the explicit predictor
  and the implicit solve differently. Each plan has a four-path table.

  **Caveat added 2026-08-07 — read before applying this rule literally.** The
  semi-implicit column is *not* yet a peer of the explicit column. It is
  **single-rank only** (`px*py` must be 1 — no halo exchange in the eta-solve
  CG), its boundary menu is **wall / tidal / clamped only**, multilayer
  requires **`vcoord_type='sigma'` and `nz_layers >= 2`**, `use_sediment` is a
  **hard config error**, and the unstructured SI predictor has no
  semi-Lagrangian advection. All are `logger%error` guards in
  `src/core/rki_config.F90` (see the table in `schism_parity_report.md` §5),
  so nothing degrades silently — but it means SI physics written today cannot
  be validated at scale or in a realistic open-boundary config.

  Consequence: either harden SI first as a prerequisite PR series, or declare
  the SI cell **deferred** per gap and say so in that plan's four-path table.
  Do not write an SI variant you cannot test and record it as done.
- **Cite the paper, never another codebase.** SCHISM may be read as a
  behavioural oracle; implementations are re-derived from primary literature.
- **Tolerances are justified, never tuned to pass.** See the tolerance
  philosophy section in the validation plan — the response to a failing case
  is never to loosen it.

## Known cross-plan dependency

Morphodynamics P8 (wave-current bed stress) needs wave-coupling Phase 5. Until
then P8 is scoped to a prescribed-wave-field seam.
