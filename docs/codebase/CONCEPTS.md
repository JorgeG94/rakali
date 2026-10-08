# Concepts and Acronyms

Quick reference for the acronyms and short conceptual labels used
throughout the Rakali codebase docs, source comments, and PR
discussion. Grouped by topic. If you find an acronym in the docs that
isn't here, please add it.

## Numerics — schemes and discretisation

| Acronym | Expansion | What it means in Rakali |
|---------|-----------|--------------------------|
| **HLL** | Harten-Lax-van Leer | Approximate Riemann solver used for the structured barotropic flux. Cheap, robust, well-balanced with hydrostatic reconstruction. |
| **HLLC** | HLL with Contact | HLL extended with a middle contact-discontinuity wave (better for stratified / contact problems). Optional in the structured solver. |
| **KNP** | Kurganov-Noelle-Petrova | Central-upwind scheme used for the unstructured barotropic flux on triangles. |
| **CU** | Central-Upwind | Same family as KNP (central-upwind = KNP-style). |
| **CWC** | Consistency With Continuity | The constraint that tracer-flux must agree with the same depth-integrated mass flux the barotropic continuity step used. Equivalently: `sum_k f_layer(k) = mass_flux_face`. Without CWC, uniform tracers don't stay uniform. |
| **BPG** | Baroclinic Pressure Gradient | Density-gradient force in the multilayer momentum equation. Two flavours implemented: sigma method and Jacobian (Shchepetkin-McWilliams 2003). |
| **BJ** | Barth-Jespersen | Slope limiter used in the unstructured MUSCL tracer advection. Caps the linear extrapolation so it doesn't introduce new extrema. |
| **GG** | Green-Gauss | Gradient reconstruction on a triangle: `grad(Tr) = (1/A) * sum_e Tr_edge * n * L_e`. Used by both the BPG and the MUSCL limiter. |
| **PCM / PLM / PPM** | Piecewise Constant / Linear / Parabolic Method | Increasing-order conservative remap reconstructions used by `rki_remap_column`. PCM is donor-cell, PLM is minmod-limited, PPM is parabolic with monotonicity (Colella-Woodward 1984). |
| **MUSCL** | Monotonic Upstream-centred Scheme for Conservation Laws | 2nd-order tracer advection on the unstructured grid: limited gradient + edge-midpoint extrapolation from the upwind cell. |
| **CFL** | Courant-Friedrichs-Lewy | Stability number `u·dt/dx`. Sets the explicit timestep. |
| **SSP-RK** | Strong Stability Preserving Runge-Kutta | Time integrator. Rakali uses SSP-RK2 (two-stage). Preserves monotonicity properties of the spatial discretisation. |
| **CG** | Conjugate Gradient | Iterative linear-system solver used for the non-hydrostatic Poisson pressure correction. |
| **PP81** | Pacanowski & Philander (1981) | Richardson-number-dependent vertical mixing scheme. `nu_t = nu_max / (1 + α·Ri)^p + nu_bg`, with p=2 for momentum and p=3 for tracers. |
| **EOS** | Equation Of State | Density formula. Rakali uses the linear two-tracer EOS `rho = rho_0 + beta_S·(S - S_ref) - alpha_T·(T - T_ref)`. |
| **PGE** | Pressure Gradient Error | Spurious horizontal pressure gradient that arises in sigma coordinates over steep bathymetry, because constant-density surfaces don't align with sigma surfaces. Mitigated by the Jacobian BPG and by ZSTAR / ZSIGMA vertical coordinates. |

## Numerics — Rakali-specific patterns

| Term | What it means |
|------|---------------|
| **Mode split** | Separating the depth-integrated barotropic step (cheap, fast wave speed) from the per-layer baroclinic distribution (slow). Each SSP-RK stage runs the barotropic update first, then distributes per-layer. |
| **Vanishing layer** | A layer with thickness ≤ `VANISHING_LAYER_TOL` (≈ `h_min`). Occurs in ZSTAR_FULL when the local depth `H < h_bed_ref`. Carries negligible mass; treated specially by the vanishing-layer-gated centred CWC advection (excluded from the per-layer flux split) and by per-layer kernels that gate on `VANISHING_LAYER_TOL` to avoid `1/h` blow-ups. |
| **Vanishing-layer-gated centred CWC** | The per-face per-layer mass-flux split used by `ml_advect_tracer*` and `ml_advect_tracer_unstr*`. For each face, layers where either side is vanishing are excluded from the per-layer flux (`v_xL(k) = 0`); the residual is redistributed across active layers via a per-face symmetric `frac_face_active = h_face_active(k) / sum h_face_active`. Reduces to the original centred CWC for sigma / zsigma / zstar-lite. |
| **Conservative soft clamp** | Two-pass column-conservative tracer limiter in `ml_clamp_tracer`. Pass 1 clips violators to `[tr_min, tr_max]` and accumulates per-cell `excess_above` and `deficit_below`; pass 2 redistributes the net imbalance into spare capacity within the same column. Per-cell tracer total preserved exactly when the column-mean is in bounds. |
| **Outer-shim + flat-impl** | Pattern for kernels that need to read array-of-derived-types fields (e.g. `tracers(it)%hTr`) from inside a `do concurrent`. The outer shim dereferences once on the host and passes flat 3D arrays into the inner impl. NVHPC stdpar can't dereference array-of-derived-types indirection from device code; the shim works around that. |
| **`_impl` duplication** | NVHPC won't inline cross-module `acc routine seq` helpers across translation units. For hot kernels we *duplicate* the helper bodies as `*_impl` copies in the calling kernel module so the compiler sees the body in the same TU and can inline it. Trade-off: parallel implementation drift risk vs. measured 3-4× speedups on flux/extrapolate. |
| **ALE** | Arbitrary Lagrangian-Eulerian | A vertical-grid philosophy where layers can move with the flow during dynamics, with periodic remapping to a fixed target grid. Rakali uses a "mode-split + remap" half-step that has some of the benefits of full ALE but rescales `dz_old` to fit `h_val`. The rescaling has known interactions with hard tracer clamps; Rakali's conservative soft clamp absorbs them. |

## Computing — hardware and toolchain

| Acronym | Expansion | What it means |
|---------|-----------|---------------|
| **GPU / CPU** | Graphics / Central Processing Unit | Rakali targets GPUs as primary (via `do concurrent` + OpenACC). CPU multicore is supported via `-stdpar=multicore`. |
| **CUDA** | Compute Unified Device Architecture | NVIDIA's GPU programming stack. Used implicitly via NVHPC's stdpar / OpenACC backends. |
| **OpenACC** | Open Accelerators | Directive-based GPU programming. Rakali uses bare directives (`!$acc enter data`, `!$acc parallel loop`) for data movement and reductions; `do concurrent` covers the rest. Inert as comments to non-OpenACC compilers. |
| **NVHPC** | NVIDIA HPC SDK | NVIDIA's Fortran compiler suite (`nvfortran`). Required for GPU offload; Rakali also supports gfortran and ifx for CPU builds. |
| **DGX** | Deep-learning GPU System | NVIDIA's multi-GPU server hardware. Most Rakali development happens on a 4-V100 DGX. |
| **SoA / AoS** | Structure of Arrays / Array of Structures | Memory layout. Rakali uses SoA (one allocatable per field, e.g. `h(:,:)`, `hu(:,:)`) for coalesced GPU access. |
| **TU** | Translation Unit | A single source file as seen by the compiler. Cross-TU inlining is the issue that motivates the `_impl` duplication pattern. |
| **IPA** | Interprocedural Analysis | Whole-program analysis that *would* enable cross-TU inlining (`-Mipa=fast,inline` in NVHPC). Slow and brittle; we don't use it — `_impl` duplication is the deliberate alternative. |
| **FP** | Floating Point | IEEE 754 floating-point arithmetic. Working precision is `wp = real64` (double). |
| **ULP** | Units in the Last Place | The smallest representable difference between two FP numbers. Used to describe rounding-level differences ("ULP-level drift"). |
| **FFI** | Foreign Function Interface | The C-compatible interface from Python / Julia into the Fortran solver, via `iso_c_binding`. |

## Domain — physics, geography, reference codes

| Acronym | Expansion | What it means |
|---------|-----------|---------------|
| **NH** | Non-Hydrostatic | The full vertical-momentum equation including `dw/dt` (not just hydrostatic balance). Rakali's NH extension uses a CG Poisson solver. Structured-only. |
| **PSU** | Practical Salinity Units | Salinity unit. Open seawater is ~34-35 PSU; Rakali uses `S_min = 0`, `S_max = 40` as default conservative bounds. |
| **SSH** | Sea Surface Height | The free-surface elevation `η` above a reference (z=0). `H = h_bed + η`. |
| **vcoord** | Vertical coordinate | Rakali-specific shorthand. The choice of how to discretise the vertical (sigma, zsigma, zstar, zstar-sigma, zstar-full). |
| **LE** | Lakes Entrance | A coastal estuary in Victoria, Australia. The canonical real-bathymetry test case for tracer / vcoord work. Has intertidal flats — relevant for the wet/dry discussion. |
| **ROMS** | Regional Ocean Modeling System | Reference for the bottom-up `k=1=bed, k=nz=surface` indexing convention (cited throughout as "ROMS-style") and the s-coordinate (Song & Haidvogel) stretching. The `VCOORD_ZSTAR_FULL` z\* design follows Adcroft & Campin 2004. |
| **ANUGA** | Australian National University Geographic Hazards Assessment tool | Reference open-source coastal-inundation solver. Rakali's unstructured KNP central-upwind on triangles mimics ANUGA's design; ANUGA is also the reference for partial-wetting schemes. |
| **SCHISM** | Semi-implicit Cross-scale Hydroscience Integrated System Model | Reference semi-implicit coastal model. Mentioned in capability discussions as the alternative regime to Rakali's explicit time-stepping. |

## I/O — file formats and conventions

| Acronym | Expansion | What it means |
|---------|-----------|---------------|
| **NetCDF** | Network Common Data Form | Self-describing binary format used for all Rakali I/O (output, restart, bathymetry, forcing, gauges). Built on HDF5. |
| **HDF5** | Hierarchical Data Format 5 | The underlying storage layer of NetCDF-4. Has a known thread-safety issue under multi-threaded NVHPC builds — tests that touch NetCDF/HDF5 run with `OMP_NUM_THREADS=1` (set via the CTest `ENVIRONMENT` property; see the gotcha in `CLAUDE.md`). |
| **UGRID** | Unstructured GRID conventions | NetCDF metadata convention for unstructured-mesh data (faces, edges, nodes, connectivity). Used by Rakali's unstructured I/O so output can be read by xarray, ParaView, etc. |
| **CF** | Climate and Forecast (metadata conventions) | Standard NetCDF metadata convention. Rakali's gauge output is CF-1.8-compliant. |

## Documentation, build, and tooling

| Acronym | Expansion | What it means |
|---------|-----------|---------------|
| **FORD** | FORtran Documenter | Tool that auto-generates HTML API docs from `!!`-prefixed Fortran comments. Rakali uses `!!` after declarations for FORD discoverability. |
| **API** | Application Programming Interface | The C-compatible callable surface in `src/api/rki_api.F90` and `rki_api_unstr.F90`. |
| **PR** | Pull Request | GitHub change-review unit. Most Rakali features ship in one or more PRs. |
| **WIP** | Work In Progress | A branch / commit that isn't ready to merge yet. |
| **TODO** | "to do" | Inline marker for known follow-up work. Larger TODOs are tracked in the roadmap / release checklist (`docs/ROADMAP_COASTAL.md`, `docs/v0.1.0_checklist.md`). |
| **CMake** | Cross-platform Make | Rakali's build system. See the CMake-options table in `CLAUDE.md`. |
| **CI** | Continuous Integration | Automated build + test on every PR. Rakali has GitHub Actions for the structured / unstructured builds. |
| **MWE / MRE** | Minimal Working / Reproducing Example | Smallest possible code that demonstrates an issue — e.g. a throwaway `do concurrent` compiler-codegen reproducer used to isolate an NVHPC/flang bug. |

## Project glossary

| Term | What it means |
|------|---------------|
| **Rakali** | The solver. Named after the Australian water rat. |
| **rki_** | Source-file prefix for all Rakali Fortran modules. |
| **`wp`** | Working precision — currently `real64` (`double` in C). All physical state and arithmetic are in `wp`. |
| **`hgrid_t`** | Horizontal-grid metadata type (cell counts, ghost width, dx/dy). |
| **`state_t`** | Top-level solution state (composition of barotropic, multilayer, NH, forcing and BC sub-objects). |
| **`mesh_unstr_t`** | Triangular-mesh metadata (nodes, edges, neighbours, areas) for the unstructured backend. |
| **Mode-split** | (See Numerics section above.) |
| **Outer-shim + flat-impl** | (See Numerics section above.) |
| **`_impl` duplication** | (See Numerics section above.) |

---

Cross-references:
- Codebase map + run lifecycle → [INDEX.md](INDEX.md)
- Per-module / per-kernel / solver call-chain detail → the source `!!` docstrings, rendered by FORD (the old `KERNELS.md` / `SOLVER.md` / `DATA_FLOW.md` were retired in favour of this)
- Coastal design contract → [`src/core/coastal/README.md`](../../src/core/coastal/README.md)
- GPU patterns and gotchas → [`FORTRAN_STYLE.md`](../../FORTRAN_STYLE.md), [CLAUDE.md](../../CLAUDE.md) → "Gotchas"
- Capability matrix and known limitations → [CAPABILITIES_AND_LIMITATIONS.md](../CAPABILITIES_AND_LIMITATIONS.md)
