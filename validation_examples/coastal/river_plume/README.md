# River plume into a saline basin

A fresh-water river channel (200 m wide, 4 m deep) empties into a
saline harbour basin (12 m deep) via a smooth bathymetric transition.
The buoyant plume spreads across the surface while a salt wedge
intrudes along the bottom.

No analytical solution — this is a visual + integral sanity check
for the multilayer S/T transport, EOS, and PP81 / KPP vertical
mixing under a strong density contrast. Both structured (Cartesian)
and unstructured (triangular) drivers live in this dir for
cross-backend comparison.

## Generate bathymetry

```bash
cd validation_examples/coastal/river_plume
python generate_bathymetry.py        # produces bathymetry.nc
```

300 × 150 cells at 10 m resolution = 3000 × 1500 m domain. The
generator builds a narrow channel on the west side, a smooth ramp
into the wider basin, and land flanking the channel.

## Run — Fortran driver

```bash
./build/rakali validation_examples/coastal/river_plume/river_plume.nml
```

2-hour simulation, dt ~ 0.5 s (CFL bounded by channel inflow).
Single V100 wallclock ~ 1 min.

Output: `output/` with NetCDF snapshots.

## Run — Python drivers

Structured:
```bash
python validation_examples/coastal/river_plume/river_plume_struct.py
```

Unstructured:
```bash
RAKALI_LIB=build/librakali.so PYTHONPATH=python \
    python validation_examples/coastal/river_plume/river_plume_unstr.py
```

Both write surface salinity, bottom salinity, x-section, and surface
velocity PNGs into the example dir.

## Sanity-check script

```bash
python python/examples/check_river_plume.py
```

Loads the canonical namelist + bathymetry, runs ~2 h, prints one
diagnostic line per simulated minute. Use to debug "is the plume
forming at all?" without a display dependency.

## What to look for

- **Surface salinity** — fresh tongue (S ≈ 0) spreading out from the
  channel mouth, deflecting (Coriolis-dependent).
- **Bottom salinity** — saline wedge (S ≈ 35) intruding upstream
  into the channel.
- **Mixed layer** at the plume base — PP81 + KPP combined produce a
  2–4-layer gradient region. With both closures off, the interface
  is one-layer thick (numerical diffusion only).

## Tuning

| Knob | Effect |
|---|---|
| `vmix_nu_max` | PP81 strength; 0 → step interface, 1e-2 → broad gradient |
| `smagorinsky_cs` | Lateral mixing; 0 → sharp plume tongue, 0.15 → diffuse |
| `coriolis_f` | Plume deflection direction (right in NH, left in SH) |
| Channel inflow (in the nml IC) | Plume size, salt-wedge penetration |

## Tracer reconstruction — PPM vs WENO5/7

The salt front is exactly where the tracer reconstruction matters: PPM is the
default parabolic limiter; the WENO ladder (`&tracer_nml tracer_recon =
weno5|weno7|weno9`, coastal multilayer only) trades a wider halo for a sharper,
still-monotone front. Three paired nmls run the same case under different
reconstructions (all with file output on, to distinct dirs):

| nml | `tracer_recon` | `nghost` | output dir |
|-----|----------------|----------|------------|
| `river_plume_ppm.nml`   | `ppm`   | 2 | `output_ppm/`   |
| `river_plume_weno5.nml` | `weno5` | 3 | `output_weno5/` |
| `river_plume_weno7.nml` | `weno7` | 4 | `output_weno7/` |

The per-rung `nghost` (3 / 4 / 5 for weno5 / 7 / 9) is a **fail-loud** guard —
under-sizing the halo aborts at configure, not silently. `weno9` needs
`nghost = 5`; add it by analogy for the full ladder. The weno rungs also
require `use_ppm_tracer=.true.` (the high-order transport branch).

```bash
cd validation_examples/coastal/river_plume
python3 generate_bathymetry.py
for r in ppm weno5 weno7; do
  CUDA_VISIBLE_DEVICES=1 ../../../build/rakali river_plume_$r.nml
done
python3 compare_recon.py            # tabulates the three runs
```

**What to expect** — the plume front should be *sharper* under WENO (fewer
transitional surface cells across the S=0↔35 interface) with **no new S
extrema** beyond the still-water IC range and salt conservation unchanged from
PPM. `compare_recon.py` reports, per run: total-salt drift, extrema growth
(overshoot/undershoot *beyond the IC bounds* — the IC already carries a saline
patch above 35 PSU, so absolute bounds are not the right check), and the
surface front-cell count.

**Reference numbers** (single V100, nvhpc 25.5, `t_end=1000 s`):

| run | wall | salt drift (rel) | overshoot | front cells |
|-----|------|------------------|-----------|-------------|
| ppm   |  ~9 s | −1.15e-3 | 0        | 706 |
| weno5 | ~11 s | −1.15e-3 | ~1e-15   | 695 (−11, sharper) |
| weno7 | ~12 s | −1.15e-3 | 0        | 704 (−2, sharper)  |

Honest limits: the front-sharpening signal here is **modest**. This case is
short (1000 s, three snapshots) and runs with a strong lateral closure
(`smagorinsky_cs=0.15`) plus explicit horizontal tracer diffusion
(`hdiff_kappa=1.0`), so physical mixing partly masks the reconstruction
difference — yet the underlying S fields differ by up to ~12 PSU between PPM
and WENO5, so the scheme *is* engaging. WENO5 gives the narrowest front here,
consistent with the module note that WENO5-Z has the best front L1 (WENO7/9 can
overshoot marginally more). The `−1.15e-3` salt "drift" is scheme-independent
outflow through the Flather open-east boundary, not a conservation defect. To
amplify the front contrast, drop `smagorinsky_cs`/`hdiff_kappa` toward 0 and
lengthen `t_end`.
