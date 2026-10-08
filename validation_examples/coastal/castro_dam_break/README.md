# Castro circular dam break

2D shallow-water dam break from Castro et al. (2011), Section 5.2.
A circular column of high water collapses radially over a Gaussian
bump bathymetry inside a [-5, 5] × [-5, 5] domain. Closed-wall BCs,
no friction, t_end = 1 s.

Used as the canonical HLL/HLLC convergence + shock-capturing test
for the coastal Cartesian path.

## Run

Single resolution:

```bash
./build/rakali validation_examples/coastal/castro_dam_break/castro_dam_break.nml
```

Default 200 × 200 (`dx = 0.05`). Override `nx`/`ny`/`dx` in the
namelist for other resolutions.

## Convergence sweep

`run_sweep.py` runs the binary at multiple grid resolutions, parses
performance + mass-error metrics from stdout, and writes a CSV.

```bash
# Single resolution with timing
python validation_examples/coastal/castro_dam_break/run_sweep.py --grid 4096

# Single resolution with NVIDIA Nsight Systems profiling
python validation_examples/coastal/castro_dam_break/run_sweep.py --grid 4096 --profile

# Single resolution with custom profiler (e.g. NCU)
python validation_examples/coastal/castro_dam_break/run_sweep.py --grid 4096 \
    --profile-line 'ncu -k "my_kernel" -c 1'

# Full sweep (100 → 4096)
python validation_examples/coastal/castro_dam_break/run_sweep.py --sweep

# Custom binary path
python validation_examples/coastal/castro_dam_break/run_sweep.py \
    --sweep --binary ./build/rakali
```

Sweep sizes: 100, 200, 400, 800, 1600, 2000, 4096.

### Output

- Console summary table: Grid, Cells, dx, Steps, Wall time, Mcells/s,
  Mass error.
- CSV: `castro_results.csv` (in the dam-break dir).

## What to look for

- **Mass error → 0** under refinement (1st-order under HLL, 2nd-order
  under HLLC with MUSCL).
- **Radial shock front** stays circular (no grid-axis preference) —
  catches limiter asymmetries.
- **No spurious flow on the bump** when the column hasn't reached it
  yet — catches BPG sign bugs and unbalanced reconstructions.

## Reference

Castro, M.J., Pardo, A., Parés, C., Toro, E.F. (2011). *On some fast
well-balanced first order solvers for nonconservative systems*.
Mathematics of Computation 79, 1427-1472.
