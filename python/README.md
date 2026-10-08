# rakali — Python interface to the Rakali shallow water solver

`rakali` is a thin Python wrapper around the Rakali GPU-native 2D
shallow water / multi-layer coastal solver written in Fortran with
`do concurrent` and OpenACC.  It exposes the solver
through a `ctypes` binding to a C-compatible API in `librakali.so`,
giving Python users zero-copy numpy views of the prognostic state and
mid-simulation hooks for forcing, boundary conditions, and data
assimilation.

## Features

- **One-shot or two-phase initialisation** — load a Fortran namelist,
  pass keyword arguments, or build the bathymetry and initial
  conditions from numpy arrays.
- **Stepping with callbacks** — `sim.advance(target, callback=fn,
  callback_interval=dt)` drives the time loop with a Python function
  invoked at each interval.  Lightweight GPU↔host sync via
  `!$acc update self` / `!$acc update device` keeps GPU data
  resident across callback iterations.
- **Mid-simulation state injection** — `set_h`, `set_hu`, `set_hv`,
  `set_salinity` for DA / observation insertion.
- **Mid-simulation forcing** — `set_wind`, `set_atm_pressure`,
  `set_precipitation` for time-varying atmospheric drivers.
- **Mid-simulation boundary conditions** — `set_inflow_velocity`,
  `set_discharge_flux`, `set_clamped_elevation`, and
  `set_tidal_constituents` for time-varying river discharge,
  tide-gauge-driven open boundaries, etc.
- **Restart save / load** — `sim.save_restart("checkpoint.nc")` and
  `Rakali.from_restart("checkpoint.nc", "config.nml")`.  Multilayer
  arrays (h_layer, hu_layer, hv_layer, hS_layer) are preserved
  across the round-trip.
- **xarray export** — `sim.to_xarray()` returns a CF-compliant
  `xarray.Dataset` snapshot with cell-centred coordinates and
  per-variable units.
- **Jupyter-friendly** — rich `_repr_html_` and an optional `tqdm`
  progress bar via `sim.advance(progress=True)`.

## Installation

`rakali` is a *Python wrapper around a Fortran library*.  You need to
build the Fortran shared library first, then install the Python
package on top of it.

### 1. Build `librakali.so`

The full build documentation is in the project's top-level README.
The minimum needed for the Python interface is:

```bash
# NVHPC GPU build
cmake -B build_nvhpc -S . -DRAKALI_ENABLE_PYTHON=ON
cmake --build build_nvhpc

# Or gfortran CPU build
cmake -B build_gcc -S . -DCMAKE_Fortran_COMPILER=gfortran \
                        -DRAKALI_ENABLE_GPU=OFF \
                        -DRAKALI_ENABLE_PYTHON=ON
cmake --build build_gcc
```

System prerequisites: a Fortran compiler (NVHPC `nvfortran` or
`gfortran`), NetCDF-Fortran, and HDF5.  These are *not*
managed by pip — install them via your system package manager (or
your HPC's module system).

### 2. Install the Python wrapper

From the repository root:

```bash
pip install -e .                    # editable install (recommended for development)
pip install -e ".[all]"             # + xarray, matplotlib, tqdm, meshpy, pyyaml
```

The pip package does **not** contain the Fortran library: it is built
separately (step 1) and found at import time. The editable install picks
up `librakali.so` from any `build*/` directory next to `pyproject.toml`.
You can point at a specific build with the `RAKALI_LIB` environment variable:

```bash
export RAKALI_LIB=/abs/path/to/librakali.so
```

### Optional dependencies

| Extra | Pulls in | Used by |
|---|---|---|
| `analysis` | `xarray` | `sim.to_xarray()` |
| `plotting` | `matplotlib` | `python/examples/live_*.py` |
| `progress` | `tqdm` | `sim.advance(progress=True)` |
| `mesh` | `meshpy` (Triangle) | `UnstructuredMesh.rectangle/disk`, the mesh builder |
| `yaml` | `pyyaml` | `rakali.yaml_loader` |
| `dev` | `pytest` | `python/tests` |
| `all` | all of the above except `dev` | everything |

## Quick start

```python
from rakali import Rakali

# From a Fortran namelist file (the river_plume case needs its bathymetry
# generated first: `python validation_examples/coastal/river_plume/generate_bathymetry.py`,
# which requires netCDF4)
sim = Rakali("validation_examples/coastal/river_plume/river_plume.nml")
sim.advance(3600.0)
print(sim.h_phys.max())   # zero-copy numpy view; valid while `sim` is open

# Programmatic configuration via kwargs
sim = Rakali(
    nx=200, ny=100, dx=10.0, dy=10.0,
    t_end=3600.0,
    bc_west="inflow", bc_east="open",
    inflow_velocity=0.5,
    use_multilayer=True, nz_layers=5,
    initial_salinity=35.0,
)

# Mid-simulation state injection (DA workflows)
sim.set_h(my_observed_h)
sim.set_salinity(my_observed_salinity)

# Time-varying forcing via callback
def update_river(s):
    Q = my_discharge_curve(s.time)
    s.set_discharge_flux(Q)

sim.advance(86400.0, callback=update_river,
            callback_interval=300.0, progress=True)

# Checkpoint and resume
sim.save_restart("checkpoint.nc")
sim2 = Rakali.from_restart("checkpoint.nc",
                            nx=200, ny=100, dx=10.0, dy=10.0,
                            t_end=86400.0,
                            bc_west="inflow", bc_east="open",
                            use_multilayer=True, nz_layers=5)

# xarray for analysis
ds = sim.to_xarray()
ds.salinity.isel(z=-1).plot()   # surface layer (layers are bottom-up: z=0 is the bed)
```

## Unstructured (triangular) grids

```python
import numpy as np
from rakali import UnstructuredMesh, UnstructuredSimulation   # needs meshpy

mesh = UnstructuredMesh.rectangle(
    x0=0.0, y0=0.0, x1=100.0, y1=2.0, max_area=0.5,
    bc_tags={"left": 2, "right": 2, "bottom": 1, "top": 1},
)
with UnstructuredSimulation(mesh) as sim:
    sim.set_bc_type(group=2, bc_type="open")        # group 1 = wall (default)
    stage = np.where(sim.centroid_x < 50.0, 2.0, 1.0)
    sim.set_bathymetry(np.zeros_like(stage))
    sim.set_initial_stage(stage)
    sim.advance(10.0)                               # prepares on first use
    h = sim.height                                  # (ncells,) zero-copy view
```

`tri_v` connectivity is **1-based**. Layers in multilayer runs are
bottom-up: index 0 is the bed layer, the last index the surface.

Views returned by `h`, `hu`, ... point into Fortran memory and keep the
simulation object alive, but become invalid after an explicit `close()`
(and after `enable_multilayer()`, which reallocates state). Copy
(`arr.copy()`) anything you need past that point. On GPU builds a held
view only refreshes when the property is read again.

Do not import `netCDF4`/`h5py` before creating a simulation in the same
process on NVHPC builds: they can shadow the HDF5 the solver links against.
Post-process after the solver has been closed.

The library is not thread-safe: use one thread per process.

## Documentation

Each public method on `rakali.Rakali` has a
detailed docstring; use `help(Rakali)` or your IDE's introspection
to browse them.

Worked examples live in `python/examples/`:

- `live_dam_break.py` — 2D dam break with matplotlib live update
- `live_river_plume.py` — 4-layer salinity plume with tqdm progress

Run the tests with the library built and `RAKALI_LIB` set (or a `build*/`
directory present):

```bash
pip install -e ".[all,dev]"
OMP_NUM_THREADS=1 pytest python/tests
```

Tests must not run in parallel (no `pytest -n`): they share one process and
the library's global logger/NetCDF state.

### Known limitation: fatal errors end the process

Unrecoverable errors inside the Fortran library (for example a malformed
namelist value or an I/O failure deep in a run) end in `ERROR STOP`, which
terminates the Python interpreter instead of raising an exception. The
wrapper checks for missing namelist, bathymetry, forcing, gauge and nesting
files up front and raises `FileNotFoundError`, but other failures can still
exit the process. Run long jobs as scripts, not inside a notebook kernel you
cannot afford to lose.

## License

Apache License 2.0.  See [LICENSE](../LICENSE) at the repository root.
