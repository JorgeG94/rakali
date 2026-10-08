# Coastal path — design contract

This directory hosts the **coastal regime** (`sim_type='coastal'`, the default) —
Rakali's production-mature path: HLL/HLLC (structured) + KNP central-upwind
(unstructured) on an A-grid; 2D barotropic, 3D hydrostatic multilayer (both
backends), structured non-hydrostatic. Tidal / estuarine / harbour / dam-break /
inundation work.

The conventions + dispatch structure you need *before* touching this path, so
you don't reconstruct them from the kernels each session.

- **Which closure/scheme is enabled where + the tunable knobs:**
  [`docs/CLOSURE_MATRIX.md`](../../../docs/CLOSURE_MATRIX.md)
- **Physics + namelist:** [`docs/REFERENCE.md`](../../../docs/REFERENCE.md)
- **Capabilities + limits:** [`docs/CAPABILITIES_AND_LIMITATIONS.md`](../../../docs/CAPABILITIES_AND_LIMITATIONS.md)
- **Build / GPU model / gotchas:** [`CLAUDE.md`](../../../CLAUDE.md)

## Vertical layer convention (load-bearing — read first)

**Bottom-up, ROMS-style: `k=1` is the bed, `k=nz` is the surface layer.**

- **Surface forcings land at the top (`k=nz`)** — wind stress, surface tracer
  fluxes (heat, E−P), the k-ε / KPP surface wall functions all act at
  `hu_layer(:,:,nz)` / `tracers(it)%hTr(:,:,nz)`.
- **Bed forcings land at the bottom (`k=1`)** — bottom drag and the bed wall
  function act at `k=1`.
- **Every multilayer kernel assumes this** — BPG, `rho_ref`, NH init, per-layer
  friction, vcoord output, surface tracer fluxes, vertical exchange. **Any new
  multilayer kernel must follow it. Don't flip.**
- Interface-staggered closures (k-ε) place `tke/eps/nut/kt` at interfaces
  `1..nz+1` (1 = bed, nz+1 = surface); interface `m` is between layers `m-1`
  (below) and `m` (above).
- Three regression gates protect the orientation: `test_lock_exchange` (dense
  current at `k=1`), `test_river_plume` (buoyant plume at `k=nz`),
  `test_surface_flux_strat` (positive Q warms `k=nz`, preserves stratification).

## State composition

`state_t` is a thin composition; the unstructured `unstr_state_t` mirrors it.

- `state%barotropic` — depth-integrated `h, hu, hv` + bathymetry + RK saves +
  the barotropic physics scalars (`manning_n`, `coriolis_f`, the linear-drag and
  `wind_surface_layer` flags).
- `state%multilayer` — per-layer `h_layer, hu/hv_layer`, the vcoord, the tracer
  registry `tracers(:)` (salinity + temperature + appended passive tracers), the
  EOS / vmix / closure scalars, and the k-ε / KPP workspace fields (allocated
  only when the closure is on).
- `state%forcing` — wind, atm pressure, precip, heat flux, E−P.
- `state%nh` — non-hydrostatic state (structured only).

No flat `state%h` / `state%h_layer` aliases — go through the sub-types.

## SSP-RK2 multilayer stage pipeline

Mode-split barotropic-then-distribute. Each RK2 stage runs, in order:

1. Layer + barotropic BCs.
2. Barotropic HLL/KNP flux (mass + momentum).
3. `ml_distribute_with_remap` — push the barotropic acceleration to the layers
   (**increment-based**: `u_k + du`, so existing per-layer shear survives) and
   reshape thicknesses for the vcoord.
4. Tracer advection (PPM structured / MUSCL unstructured) + horizontal diffusion,
   looped over the registry.
5. Vertical exchange (diagnostic `w`).
6. EOS (`rho_layer` from S, T) — a no-op if `alpha_T = beta_S = 0`.
7. **Vertical mixing** — the closure dispatch (see below).
8. Surface tracer fluxes (heat, E−P) at `k=nz`; conservative clamp.
9. Baroclinic pressure gradient (Jacobian / sigma).
10. Per-layer bottom friction (Manning or linear) at `k=1`.
11. Smagorinsky + biharmonic horizontal viscosity.
12. CFL velocity truncation (both backends).
13. Per-layer Coriolis (closed-form rotation).
14. Surface wind stress to the top layer, if `wind_surface_layer`.
15. Recompute barotropic from the layers; barotropic wind / pressure / precip
    (the wind is suppressed here when step 14 already routed it to the top layer).

**Steps 5–14 are gated on `nz > 1`.** At `nz=1` (and in pure 2D) the whole
vertical block is a no-op — there is no interior interface to mix across.

## Vertical-mixing dispatch (the closure selector)

Inside the `nz > 1` block the vmix selector is a precedence chain:

```
if (vmix_use_keps) then       ! k-ε: friction velocities → ml_keps_update →
                              !   apply ν_t to momentum + κ_t to tracers (implicit)
else if (vmix_use_kpp) then   ! KPP overlay on the PP81 interior
else                          ! PP81 (default)
```

All three closures run on **both backends**. KPP / k-ε source
turbulence from surface wind + bed stress via u\* / wall functions. The k-ε apply
kernels are **implicit** (Thomas) because ν_t can reach `NUT_MAX`, which an
explicit flux clamp would gut.

Under `time_scheme="semi_implicit"` the same precedence chain runs inside
`si_unstr_step_ml`, but every closure is a pure FILL of the
SI column solve's `nu_iface` slot (the momentum apply IS the implicit
column tridiagonal; k-ε reuses `ml_keps_update_unstr` at SI cadence +
optional `keps_si_picard`; KPP's overlay is MAX'd in WITHOUT the explicit
CFL clamp and its S/T diffusion + non-local run implicitly —
`rki_si_unstr_closures`).

## Conventions every kernel follows

- **Bit-identity by default.** Every new closure / scheme is a namelist knob,
  default off, so existing nmls + tests stay byte-identical. Allocations + GPU
  mappings gate on the knob; the dispatch falls through to the prior behaviour.
  PR recipe: knob → kernel → unit test → **closure-matrix row** → commit.
- **GPU-resident.** State stays on device; only diagnostic output triggers D→H.
  Data-parallel loops are `do concurrent` (`j,i` order for NVHPC coalescing);
  per-step kernels never local-allocate scratch (module-level workspaces, lazy
  `*_ensure` / `*_cleanup`).
- **Tracer registry.** Passive transport (advection, remap, vertical exchange,
  vert/horiz diffusion, clamp, RK2) loops `tracers(:)` — append a passive tracer
  with `multilayer%register_passive_tracer(...)` and it "just works" (EOS +
  surface flux stay S/T-specific). Registry tracers passed into a `do concurrent`
  use the outer-shim + flat-impl pattern (array-of-derived-type device
  indirection).
- **Structured ⟷ unstructured.** Array shapes are opposite: structured
  `(nx,ny,nz)`, unstructured `(nz, ncells)`. The unstructured flux is a
  cell-centric double-visit (each triangle accumulates only into its own arrays —
  writing to a neighbour is a GPU race). **Not yet ported to unstructured:**
  velocity truncation, per-layer bottom friction, fused ML kernel, biharmonic
  viscosity, k-ε, linear bottom drag, surface-layer wind stress.

Per-module / per-procedure detail lives in the source docstrings (FORD). This
file is the cross-cutting contract those docstrings can't express.
