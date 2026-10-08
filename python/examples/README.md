# Rakali Python examples

Self-contained Python scripts demonstrating the Rakali Python API.

## Setup

Build Rakali with the shared library enabled:

```bash
cmake -B build -S . -DRAKALI_ENABLE_PYTHON=ON && cmake --build build
```

The example scripts auto-discover `librakali.so` in `build_nvhpc/`,
`build_gcc/`, or `build/` (in that order).  You can override with the
`RAKALI_LIB` environment variable:

```bash
export RAKALI_LIB=/path/to/build/librakali.so
```

## Dependencies

| Package | Required | Used by |
|---|---|---|
| `numpy` | yes | all examples |
| `matplotlib` | live plotting | all `live_*.py` and `interactive_*.py` |
| `tqdm` | optional | any example using `progress=True` |
| `xarray` | optional | snapshot collection / analysis |

Install the optional bits with:

```bash
pip install matplotlib tqdm xarray
```

## Examples

### `live_dam_break.py`

A 200x200 square basin (1 km on a side) with a higher-water region in
the south-west quadrant.  Releases the "dam" at t=0 and animates the
surface depth evolving for 60 s of simulation time.  Demonstrates:

- Programmatic two-phase init (`defer_init=True`)
- `set_bathymetry` / `set_initial_h` from numpy arrays
- The `callback` parameter on `sim.advance()` driving live matplotlib
  updates via `imshow.set_data()`

Run:

```bash
python python/examples/live_dam_break.py
```

### `live_river_plume.py`

A 4-layer 12 km x 8 km basin with a fresh-water inflow on the western
edge (5 PSU vs 35 PSU background).  Watches the surface and bottom
salinity plumes develop side-by-side over 1 hour of simulation time.
Demonstrates:

- Multilayer setup with the Jacobian BPG and Smagorinsky viscosity
- `bc_west="inflow"` with `inflow_velocity` and `inflow_salinity`
- The `progress=True` flag on `sim.advance()` (uses `tqdm` if available)
- Reading per-layer salinity via `sim.salinity[:, :, k]`

Run:

```bash
python python/examples/live_river_plume.py
```

### `live_vorticity_2d.py`

Uniform barotropic inflow past a circular island in a 6 km x 3 km
channel.  Watches a Kármán vortex street shed from the island via
a live vorticity field ζ = ∂v/∂x − ∂u/∂y computed with `np.gradient`.
Demonstrates:

- Pure 2D barotropic mode (`use_multilayer=False`)
- Smagorinsky viscosity keeping the shear layer stable
- Computing a diagnostic from `hu_phys`/`hv_phys` each callback

Run:

```bash
python python/examples/live_vorticity_2d.py
```

### `live_lock_exchange.py`

Closed flat-bottomed channel with a vertical salinity gate at the
middle — dense water on the left, light on the right.  Releasing
the gate launches two counter-propagating gravity-current noses
with Kelvin-Helmholtz billows along the interface.  Rendered as an
(x, z) cross-section rebuilt from `h_layer`.  Demonstrates:

- Multilayer Jacobian BPG with Smagorinsky + horizontal salinity diffusion
- `set_salinity()` for per-layer initial stratification
- `pcolormesh` cross-section from sigma layer thicknesses

Run:

```bash
python python/examples/live_lock_exchange.py
```

### `live_tidal_estuary.py`

Idealised funnel-shaped estuary (narrow head, wide mouth) with an
M2 tidal elevation boundary on the east end and a steady freshwater
river inflow at the west head.  Twin panels show surface salinity
(the salt wedge breathing with the tide) and η.  Demonstrates:

- `bc_east="tidal"` + `set_tidal_constituents()` for M2 forcing
- `bc_west="inflow"` as a river source with `inflow_salinity`
- Multilayer stratification initialised via `set_salinity()`

Run:

```bash
python python/examples/live_tidal_estuary.py
```

### `live_solitary_wave.py`

Long narrow channel with a Gaussian free-surface bump released at
t=0.  With non-hydrostatic pressure enabled, frequency dispersion
balances nonlinear steepening and the bump unwinds into a train
of solitary-like waves.  Shown as a 1D η profile along the channel
centreline.  Demonstrates:

- `use_nonhydrostatic=True` + multilayer sigma coordinates
- Barotropic/baroclinic mode splitting (`split_ratio=0` auto)
- Line-plot animation instead of `imshow`

Run:

```bash
python python/examples/live_solitary_wave.py
```

### `live_kelvin_wave.py`

Rotating rectangular basin (f > 0) with a Gaussian bump pressed
against the south wall, decaying away from the coast with the
Rossby radius of deformation.  The bump sheds a coastally trapped
Kelvin wave that propagates eastward along the southern boundary.
Demonstrates:

- Coriolis forcing (`coriolis_f`)
- Barotropic SWE producing trapped waves from a simple IC
- Analytical Rossby-radius scaling used in the initial condition

Run:

```bash
python python/examples/live_kelvin_wave.py
```

### `interactive_obstacle.py`

Steady westerly inflow through a rectangular channel.  **Left-click**
anywhere on the plot to drop a Gaussian atmospheric pressure low at
that location — the inverse-barometer response lifts the free
surface and radiates gravity waves.  **Right-click** (or press `c`)
to clear all perturbations.  Demonstrates:

- Mid-simulation forcing updates via `set_atm_pressure()`
  (one of the few fields explicitly safe to update while the solver
  is running; bathymetry is baked in at `finalize()`)
- Matplotlib mouse/key events wired into the `callback=` hook

Run:

```bash
python python/examples/interactive_obstacle.py
```

## Tips

- **Headless runs**: matplotlib's interactive mode (`plt.ion()`) works
  with the `Agg` backend too — set `MPLBACKEND=Agg` to skip the
  on-screen window if you only want the simulation to run.
- **Saving an animation**: replace `plt.show()` calls with
  `fig.savefig(f"frame_{step:04d}.png")` inside the callback to dump
  frames, then stitch with `ffmpeg`.
- **Jupyter**: the same scripts work inline with `%matplotlib inline`
  (one final frame) or `%matplotlib widget` (live updates) — no code
  changes needed.  The `Rakali` class also has a rich `_repr_html_`
  so just typing `sim` in a notebook cell shows a formatted summary
  of grid/time/field ranges.
