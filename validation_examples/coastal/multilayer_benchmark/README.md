# Multilayer benchmark

3D baroclinic dam break with Coriolis and Pacanowski-Philander
vertical mixing. Exercises the full coupled coastal multilayer
solver on GPU — barotropic-then-distribute mode-split, BPG, S/T
advection + diffusion, EOS, PP81 vmix, friction.

Used as a throughput + mass-conservation regression target. Also
referenced by the ALE-remap notebook
(`python/notebooks/11_ale_remap/`) as the closed-basin
mass-conservation reproducer.

## Run

```bash
./build/rakali validation_examples/coastal/multilayer_benchmark/multilayer_benchmark.nml
```

Default 2000 × 2000 × 5 = 20 M layer-cells, dx = 10 m, 20 × 20 km
domain, 200 s integration (~200 steps at dt ~ 1 s).

Single V100 wallclock: ~30–60 s.

### Smaller / larger variants

Override in the namelist:

| Setting | Use |
|---|---|
| `nx = ny = 500` | Quick smoke (~60 k cells, finishes in seconds) |
| `nx = ny = 2000` (default) | Regression baseline |
| `nx = ny = 4000` | Heavy profiling (16 M cells × 5 layers) |

## What to look for

- **Throughput** (Mcells/s) — should track the bench numbers in
  `CLAUDE.md` "Performance" section.
- **Mass conservation** — the closed basin should hold mass to
  ~1e-12 relative under conservative remap. The
  `tests/analytical/test_zstar_full_tracer_conservation.F90`
  gate captures the regression form of this; the benchmark here is
  the long-run, full-physics version.
- **No spurious flow** before the dam release — sanity check on
  initial-condition + BPG setup.

## Profiling output

Pre-existing `report*.nsys-rep` / `report*.sqlite` files in this
dir are Nsight-Systems traces from prior profiling runs (April).
They're not tracked — kept locally for historical comparison.
