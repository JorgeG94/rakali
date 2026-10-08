# Rakali

A coastal, estuarine and tidal ocean solver in Fortran. It runs on GPUs and CPUs from one
source: `do concurrent` does the parallelism, with a thin layer of OpenACC for reductions and
device data. It builds with NVHPC (GPU and multicore), gfortran 15 and ifx, and has a Python
interface.

Rakali is under active development and will not be stable until v1.0.0: expect breaking changes
to the namelist, the Python API and the output files. It is an experiment, not a production
model. See [Status](#status) before relying on it.

## What it does

- **Grids:** structured Cartesian and unstructured triangular, both on an A-grid.
- **Equations:** 2D shallow water, 3D hydrostatic multilayer, and non-hydrostatic (structured only).
- **Numerics:** HLL/HLLC on structured grids, Kurganov-Noelle-Petrova central-upwind on triangles,
  SSP-RK2, and a semi-implicit mode on triangles. Wetting and drying are well-balanced.
- **Vertical:** sigma, z-sigma and the z-star family, with conservative ALE remap.
- **Physics:** salinity and temperature plus passive tracers, vertical mixing (PP81, KPP, k-epsilon),
  Smagorinsky and biharmonic viscosity, wind, pressure and surface fluxes.
- **Boundaries:** wall, open, tidal, nested, inflow, discharge, clamped, sponge, Chapman.
- **I/O:** NetCDF (CF-1.8) output, restart, gauges.

[`docs/CAPABILITIES_AND_LIMITATIONS.md`](docs/CAPABILITIES_AND_LIMITATIONS.md) lists what works and
what doesn't. [`docs/CLOSURE_MATRIX.md`](docs/CLOSURE_MATRIX.md) says which closure runs in which regime.

## Build

```bash
cmake -B build -S . && cmake --build build
cd build && ctest --output-on-failure      # never pass -j to ctest on a GPU build
```

| Option | Default | |
|---|---|---|
| `RAKALI_ENABLE_GPU` | `OFF` | GPU offload (NVHPC). `ON` with a compiler that can't offload is an error. |
| `RAKALI_ENABLE_PYTHON` | `OFF` | Also build `librakali.so` for the Python interface. |
| `RAKALI_ENABLE_NETCDF` | `ON` | NetCDF I/O. Turn off for kernel-only builds. |
| `RAKALI_GPU_ARCH` | `cc70` | `cc70`, `cc80`, `cc90`. |
| `RAKALI_EXTRA_FORTRAN_FLAGS` | empty | Appended last, e.g. `-tp=haswell` or `-march=x86-64-v3`. |

Requires CMake, a Fortran compiler (gfortran 15 or newer, ifx, or NVHPC) and NetCDF-Fortran.
[pic](https://github.com/JorgeG94/pic) and test-drive are fetched by CMake. `tools/build_netcdf_fortran.sh`
builds NetCDF-Fortran for your compiler on top of an existing NetCDF-C.

## Run

```bash
./build/rakali validation_examples/coastal/castro_dam_break/castro_dam_break.nml
```

## Python

```bash
cmake -B build -S . -DRAKALI_ENABLE_PYTHON=ON && cmake --build build --target rakali_shared
pip install -e ".[all]"
```

```python
from rakali import Rakali
with Rakali(nx=200, ny=100, dx=10.0, dy=10.0, t_end=3600.0,
            bc_west="wall", bc_east="wall") as sim:
    sim.advance(60.0)
    print(sim.h_phys.max())      # zero-copy numpy view
```

See [`python/README.md`](python/README.md), including the unstructured API and known limitations.

## Status

Most of Rakali's code was written by Claude, Anthropic's AI model, working through Claude Code.
The design is Jorge Luis Galvez Vallejo's: the architecture, the API, the memory model and the
parallelism model. It claims no novelty.

There is no reference model to check against. The unstructured scheme follows the approach of
[ANUGA](https://github.com/anuga-community/anuga_core); everything else (the structured solver,
multilayer, vertical coordinates, closures) was designed here. Correctness rests on the tests in
`tests/`: analytical cases (lake at rest, Ritter and Thacker, Carrier-Greenspan, seiche
convergence, Kato-Phillips, lock exchange) and conservation checks.

Treat these as experimental: k-epsilon, `zstar_full` (it leaks 1-2% of salt per tidal cycle in
intertidal estuaries), the non-hydrostatic solver, and nested boundaries on triangles.
Fatal errors inside the Fortran library end in `ERROR STOP`, which kills a Python session.

## Documentation

- [`docs/REFERENCE.md`](docs/REFERENCE.md): physics and namelist reference
- [`src/core/coastal/README.md`](src/core/coastal/README.md): design contract
- [`docs/codebase/INDEX.md`](docs/codebase/INDEX.md): repo map and run lifecycle
- [`FORTRAN_STYLE.md`](FORTRAN_STYLE.md): style and the GPU programming guide
- API docs are generated from the source docstrings with FORD and published by CI.

## License

Apache 2.0. See [`LICENSE`](LICENSE).
