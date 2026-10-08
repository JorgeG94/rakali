# Validation examples

Canonical / analytical-solution benchmarks for the Rakali dyn-core.
Each subdir documents a setup that has a known reference behaviour
(analytical solution, published benchmark, or canonical regression
config) — runnable end-to-end with the in-tree binary and small
enough to spin up on a single GPU in minutes.

Real-bathymetry case studies (data not distributed) live under
[`python/examples/real_bathymetry/`](../python/examples/real_bathymetry/) instead.

## Layout

```
validation_examples/
  coastal/    Canonical configs for the coastal (HLL/HLLC/KNP) path
    castro_dam_break/          Castro et al. (2011) circular dam break
    river_plume/               fresh-water plume into a saline basin
    multilayer_benchmark/      throughput / mass-conservation regression
```

## Coastal benchmarks

| Subdir | What it tests | Reference |
| ------ | ------------- | --------- |
| `coastal/castro_dam_break/` | Shock-capturing on a 2D circular dam break over a Gaussian bump; convergence under refinement. | Castro et al. (2011). |
| `coastal/river_plume/` | Fresh-water buoyant plume + salt-wedge intrusion; structured (Cartesian) and unstructured (triangular) drivers. | No analytical solution; visual + Python sanity check via `check_river_plume.py`. |
| `coastal/multilayer_benchmark/` | 3D baroclinic-dam-break throughput + closed-basin mass conservation regression. | Internal regression target; no published reference. |

## Running anything here

From the repo root:

```bash
./build/rakali validation_examples/coastal/<subdir>/<file>.nml
```

Each subdir has its own README with expected wallclock, what to
look for in the output, and the relevant analytical scales / tuning
knobs.

## Adding a new benchmark

Goes here if it has a published / textbook / analytical reference,
or is a canonical regression target. Real-bathymetry / case-study
configs go in [`python/examples/real_bathymetry/`](../python/examples/real_bathymetry/) instead.
