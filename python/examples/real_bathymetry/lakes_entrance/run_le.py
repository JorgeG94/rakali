#!/usr/bin/env python
"""72-hour storm-surge simulation of Lakes Entrance (structured solver).

Driven from Python so the south (Bass Strait) boundary forcing can stack
three things explicitly:

    eta(t) = ramp(t) * tide(t) + surge(t)

  * tide(t)  — 8 EOT20 constituents at (148.0 E, -38.0 S) for epoch
               2026-04-29 00:00:00 UTC.
  * ramp(t)  — cosine spin-up 0->1 over the first TIDE_RAMP seconds, so
               the cold-start transient (still water → real tide) is
               smeared out instead of slammed in.
  * surge(t) — 1.5 m cosine bump centred at hour 36, ramping up from
               hour 24 and back down by hour 60.

The south boundary uses BC_CLAMPED + set_clamped_elevation per timestep,
so the kernel sees exactly the eta we compute here (no BC_TIDAL phase
sign convention to worry about).

Console print:  every 10 min of sim time
NetCDF snapshot: every 30 min of sim time
Total run     :  72 h

Usage:
    PYTHONPATH=/home/jorgenci/dev/rakali/python \\
        python run_lakes_entrance_surge.py
"""
from __future__ import annotations

import argparse
import os
import sys
import time as walltime
from pathlib import Path

# NetCDF/HDF5 isn't thread-safe under -stdpar; pin before any threaded import.
os.environ.setdefault('OMP_NUM_THREADS', '1')

import numpy as np
import xarray as xr

# Make the in-repo rakali package importable without setting PYTHONPATH:
# this file is <repo>/python/examples/real_bathymetry/lakes_entrance/run_le.py,
# so parents[4] = <repo>,
# and <repo>/python holds the package. librakali.so is then auto-discovered
# from <repo>/build*/ by rakali._core (RAKALI_LIB still overrides if set).
_RAKALI = Path(__file__).resolve().parents[4]
for p in [str(_RAKALI / 'python'), str(_RAKALI)]:
    if p not in sys.path:
        sys.path.insert(0, p)

from rakali import Rakali

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
BATHY_FILE = HERE / 'lakes_entrance_30m.nc'
OUT_NC = HERE / 'lakes_entrance_surge_72h.nc'
OUT_NC_HURRICANE = HERE / 'lakes_entrance_surge_72h_hurricane.nc'

# Grid (matches lakes_entrance_30m.nc)
NX, NY = 431, 326
DX = 30.0

# Tide: EOT20 @ (148.0 E, -38.0 S), epoch 2026-04-29 00:00:00 UTC, no nodal.
# Convention: eta_tide(t) = sum_k A_k * cos(omega_k * t + phi_k)
# Order: M2, S2, K1, O1, N2, K2, P1, Q1
TIDAL_AMP = np.array([0.320132, 0.091847, 0.175943, 0.111556,
                      0.087435, 0.029379, 0.054063, 0.028405])
TIDAL_PHASE = np.array([0.486437, -0.974482, 1.550926, -2.213151,
                        -1.575488, 0.099292, -2.991708, 1.900360])
TIDAL_OMEGA = np.array([1.405189e-04, 1.454441e-04, 7.292117e-05, 6.759774e-05,
                        1.378797e-04, 1.458423e-04, 7.252295e-05, 6.495854e-05])

# Run length and forcing schedule
T_END = 72 * 3600.0          # 72 hours
TIDE_RAMP = 3 * 3600.0       # cosine ramp 0->1 over first 3 h

SURGE_PEAK = 2.5             # m
SURGE_T_START = 24 * 3600.0  # surge appears at hour 24
SURGE_T_PEAK  = 36 * 3600.0  # peaks at hour 36
SURGE_T_END   = 60 * 3600.0  # back to zero at hour 60

# Time-step intervals (sim time, not wall time)
DT_BC = 60.0                 # boundary update + callback granularity
DT_PRINT = 1800.0             # console print every 10 min
DT_SNAP = 1800.0             # snapshot every 30 min

CFL = 0.45
MANNING = 0.025

# Hurricane-mode wind forcing (uniform field, ramped in around the surge).
# Stress ~ rho_air * Cd * |U10|^2; tau = 8 Pa corresponds to roughly Cat-4
# winds (~55 m/s) with Cd ~ 2.5e-3.  Direction is +y (southerly blowing
# onshore), i.e. stacks with the surge instead of fighting it.
HURRICANE_TAU_PEAK = 8.0          # Pa
HURRICANE_DIR_DEG = 0.0           # 0 = +y (north-bound), 90 = +x (east-bound)
HURRICANE_T_START = 18 * 3600.0
HURRICANE_T_PEAK  = 36 * 3600.0
HURRICANE_T_END   = 54 * 3600.0


# ---------------------------------------------------------------------------
# Boundary forcing
# ---------------------------------------------------------------------------
def tide(t: float) -> float:
    return float(np.sum(TIDAL_AMP * np.cos(TIDAL_OMEGA * t + TIDAL_PHASE)))


def tide_ramp(t: float) -> float:
    """Cosine ramp 0 -> 1 over [0, TIDE_RAMP]."""
    if t >= TIDE_RAMP:
        return 1.0
    return 0.5 * (1.0 - np.cos(np.pi * t / TIDE_RAMP))


def surge(t: float) -> float:
    """1.5 m cosine bump from SURGE_T_START -> SURGE_T_END, peak at PEAK."""
    if t <= SURGE_T_START or t >= SURGE_T_END:
        return 0.0
    if t <= SURGE_T_PEAK:
        x = (t - SURGE_T_START) / (SURGE_T_PEAK - SURGE_T_START)
    else:
        x = (SURGE_T_END - t) / (SURGE_T_END - SURGE_T_PEAK)
    return SURGE_PEAK * 0.5 * (1.0 - np.cos(np.pi * x))


def eta_bc_components(t: float):
    """Return (ramped_tide, surge, total) at sim time t."""
    rt = tide_ramp(t) * tide(t)
    sg = surge(t)
    return rt, sg, rt + sg


def hurricane_tau(t: float) -> float:
    """Cosine-bump wind stress magnitude (Pa). Zero outside [START, END]."""
    if t <= HURRICANE_T_START or t >= HURRICANE_T_END:
        return 0.0
    if t <= HURRICANE_T_PEAK:
        x = (t - HURRICANE_T_START) / (HURRICANE_T_PEAK - HURRICANE_T_START)
    else:
        x = (HURRICANE_T_END - t) / (HURRICANE_T_END - HURRICANE_T_PEAK)
    return HURRICANE_TAU_PEAK * 0.5 * (1.0 - np.cos(np.pi * x))


# ---------------------------------------------------------------------------
# Build simulation
# ---------------------------------------------------------------------------
def main(hurricane: bool = False, checkpoint_dir: 'Path | None' = None,
         checkpoint_stride: int = 1):
    if checkpoint_dir is not None:
        checkpoint_dir = Path(checkpoint_dir)
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
        if checkpoint_stride < 1:
            raise ValueError('--checkpoint-stride must be >= 1')
        print(f'*** PER-STEP CHECKPOINTS enabled: writing to {checkpoint_dir} ***')
        print(f'    stride = {checkpoint_stride} (saving every '
              f'{checkpoint_stride}{"st" if checkpoint_stride == 1 else "th"} '
              'solver step)')

    if hurricane:
        print(f'*** HURRICANE MODE: peak tau = {HURRICANE_TAU_PEAK:.1f} Pa '
              f'(~Cat-4 winds, ~55 m/s) from hour '
              f'{HURRICANE_T_START/3600:.0f}-{HURRICANE_T_END/3600:.0f}, '
              f'direction = {HURRICANE_DIR_DEG:.0f} deg ***')
    print(f'Loading bathy from {BATHY_FILE}')
    b = xr.open_dataset(BATHY_FILE)['b'].values  # (nx, ny)
    assert b.shape == (NX, NY), f'expected ({NX},{NY}), got {b.shape}'
    b_F = np.asfortranarray(b.astype(np.float64))
    print(f'  shape={b.shape}, range=[{b.min():.1f}, {b.max():.1f}] m')

    print('Building Rakali simulation (structured)...')
    sim = Rakali(
        nx=NX, ny=NY, dx=DX, dy=DX, t_end=T_END,
        cfl=CFL, manning_n=MANNING, dry_threshold=0.01,
        # All sides walled except south, which we clamp from Python.
        bc_west='wall', bc_east='wall', bc_north='wall',
        bc_south='clamped',
        log_level='warning',
        defer_init=True,
    )
    sim.set_bathymetry(b_F)
    sim.set_initial_h(np.asfortranarray(np.maximum(-b_F, 0.0)))
    sim.finalize()

    # Print forcing checkpoints so we can sanity-check the schedule.
    print('\nForcing schedule (eta at south boundary, m):')
    for h in (0, 1, 3, 12, 24, 30, 36, 42, 48, 60, 72):
        rt, sg, tot = eta_bc_components(h * 3600.0)
        print(f'  t = {h:2d} h:  tide_ramped = {rt:+.3f}  '
              f'surge = {sg:+.3f}  total = {tot:+.3f}')

    # ---- snapshot + print state held in a closure-shared dict ----
    state = {
        'next_print': DT_PRINT,
        'next_snap':  DT_SNAP,
        'snaps':      [],
        'eta_total':  [],
        'eta_tide':   [],
        'eta_surge':  [],
        't_snap':     [],
        'wall_start': walltime.perf_counter(),
        'wall_window': walltime.perf_counter(),
        'steps_window_start': 0,
        'tau_now':    0.0,
    }

    # Pre-allocate the wind buffers so we're not malloc'ing every callback.
    # Solver expects (nx_phys, ny_phys) Fortran-ordered.
    if hurricane:
        wx_buf = np.zeros((NX, NY), dtype=np.float64, order='F')
        wy_buf = np.zeros((NX, NY), dtype=np.float64, order='F')
        dir_rad = np.deg2rad(HURRICANE_DIR_DEG)
        wx_unit, wy_unit = np.sin(dir_rad), np.cos(dir_rad)
    else:
        wx_buf = wy_buf = None
        wx_unit = wy_unit = 0.0

    def step_cb(s):
        rt, sg, tot = eta_bc_components(s.time)
        s.set_clamped_elevation(tot)

        if hurricane:
            tau = hurricane_tau(s.time)
            state['tau_now'] = tau
            wx_buf.fill(tau * wx_unit)
            wy_buf.fill(tau * wy_unit)
            s.set_wind(wx_buf, wy_buf)

        # 10-minute console print
        if s.time + 1e-6 >= state['next_print']:
            wall_now = walltime.perf_counter()
            wall_total = wall_now - state['wall_start']
            wall_dt = wall_now - state['wall_window']
            steps_dt = s.step_count - state['steps_window_start']
            sps = steps_dt / wall_dt if wall_dt > 0 else 0.0
            speedup = DT_PRINT / wall_dt if wall_dt > 0 else 0.0
            tau_str = (f'  tau = {state["tau_now"]:4.1f} Pa'
                       if hurricane else '')
            print(f'  t = {s.time/3600:6.2f} h  step = {s.step_count:7d}  '
                  f'dt = {s.dt:.3e} s  '
                  f'eta_tide = {rt:+.3f}  surge = {sg:+.3f}  '
                  f'tot = {tot:+.3f}{tau_str}  |  wall = {wall_total:7.1f} s  '
                  f'({sps:5.0f} steps/s, {speedup:5.1f}x rt)',
                  flush=True)
            state['wall_window'] = wall_now
            state['steps_window_start'] = s.step_count
            state['next_print'] += DT_PRINT

        # 30-minute snapshot
        if s.time + 1e-6 >= state['next_snap']:
            ds = s.to_xarray()
            state['snaps'].append(ds)
            state['t_snap'].append(s.time)
            state['eta_tide'].append(rt)
            state['eta_surge'].append(sg)
            state['eta_total'].append(tot)
            state['next_snap'] += DT_SNAP

    ncells = NX * NY
    print(f'\nRunning {T_END/3600:.0f} h on {ncells} cells '
          f'({NX} x {NY} @ dx={DX:.0f} m) '
          f'(boundary update every {DT_BC:.0f}s, '
          f'print every {DT_PRINT:.0f}s, snapshot every {DT_SNAP:.0f}s)...\n')

    t0 = walltime.perf_counter()
    if checkpoint_dir is None:
        sim.advance(T_END, callback=step_cb, callback_interval=DT_BC)
    else:
        # Per-solver-step loop: apply BC + (optional) wind, step exactly one
        # dt, write the mid-RK checkpoint, then run the normal callback so
        # the existing 30-min snapshot + 10-min print logic still fires.
        import ctypes
        rki_step = sim._lib.rki_step
        rki_sync = sim._lib.rki_sync_host
        handle = sim._handle
        ckpt_pad = 7  # supports up to 9,999,999 steps (~10M)
        next_bc = 0.0
        while sim.time < T_END - 1e-9:
            # Refresh BC + wind every DT_BC of sim time (same cadence as
            # the old callback path) so the south-edge eta stays current.
            if sim.time + 1e-9 >= next_bc:
                rt, sg, tot = eta_bc_components(sim.time)
                sim.set_clamped_elevation(tot)
                if hurricane:
                    tau = hurricane_tau(sim.time)
                    state['tau_now'] = tau
                    wx_buf.fill(tau * wx_unit)
                    wy_buf.fill(tau * wy_unit)
                    sim.set_wind(wx_buf, wy_buf)
                next_bc += DT_BC

            rc = rki_step(handle, ctypes.c_int(1))
            if rc != 0:
                raise RuntimeError(f'rki_step failed (rc={rc})')
            rki_sync(handle)

            if sim.step_count % checkpoint_stride == 0:
                fname = checkpoint_dir / f'step_{sim.step_count:0{ckpt_pad}d}.nc'
                sim.save_checkpoint(fname)

            # Reuse the existing 10-min print / 30-min snapshot logic.
            step_cb(sim)
    elapsed = walltime.perf_counter() - t0
    print(f'\nDone: {sim.step_count} steps in {elapsed:.1f} s '
          f'({sim.step_count / elapsed:.0f} steps/s, '
          f'{T_END / elapsed:.1f}x real-time)')

    # ---- Stitch snapshots into one Dataset and write to NetCDF ----
    out_nc = OUT_NC_HURRICANE if hurricane else OUT_NC
    print(f'\nWriting {len(state["snaps"])} snapshots to {out_nc} ...')
    ds = xr.concat(state['snaps'], dim='time')
    ds = ds.assign_coords(
        eta_tide=('time', np.asarray(state['eta_tide'], dtype=np.float32)),
        eta_surge=('time', np.asarray(state['eta_surge'], dtype=np.float32)),
        eta_bc=('time', np.asarray(state['eta_total'], dtype=np.float32)),
    )
    ds.attrs['description'] = (
        'Lakes Entrance 72-hour storm-surge simulation. South boundary '
        'driven via BC_CLAMPED with eta(t) = ramped_tide(t) + surge(t).'
    )
    ds.attrs['surge_peak_m'] = SURGE_PEAK
    ds.attrs['surge_t_start_h'] = SURGE_T_START / 3600.0
    ds.attrs['surge_t_peak_h']  = SURGE_T_PEAK / 3600.0
    ds.attrs['surge_t_end_h']   = SURGE_T_END / 3600.0
    ds.attrs['tide_ramp_h']     = TIDE_RAMP / 3600.0

    if hurricane:
        ds.attrs['hurricane_tau_peak_pa'] = HURRICANE_TAU_PEAK
        ds.attrs['hurricane_dir_deg']     = HURRICANE_DIR_DEG
        ds.attrs['hurricane_t_start_h']   = HURRICANE_T_START / 3600.0
        ds.attrs['hurricane_t_peak_h']    = HURRICANE_T_PEAK / 3600.0
        ds.attrs['hurricane_t_end_h']     = HURRICANE_T_END / 3600.0

    enc = {v: {'zlib': True, 'complevel': 1}
           for v in ds.data_vars if ds[v].dtype.kind == 'f'}
    ds.to_netcdf(out_nc, encoding=enc)
    print(f'  wrote {out_nc}')

    sim.destroy()


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--hurricane', action='store_true',
                   help='Add a Cat-4 wind-stress field ramped around the '
                        'surge peak. Output goes to '
                        f'{OUT_NC_HURRICANE.name} instead of {OUT_NC.name}.')
    p.add_argument('--checkpoint', metavar='DIR', default=None,
                   help='Write a NetCDF mid-RK checkpoint per solver step '
                        'into DIR (step_NNNNNNN.nc). Each file holds '
                        'h, hu, hv, b, h0, hu0, hv0 (full ghost-included '
                        'grid) plus scalars (time, dt, step_count, '
                        'manning_n, coriolis_f, clamped_elevation). '
                        'Disables the fast advance path and switches to a '
                        'per-step Python loop. WARNING: a 72h run produces '
                        '~600k files at ~few MB each — use --checkpoint-stride '
                        'to thin.')
    p.add_argument('--checkpoint-stride', metavar='N', type=int, default=1,
                   help='Save a checkpoint every N solver steps '
                        '(default 1 = every step). N=25 produces ~22k files '
                        '(~100 GB compressed) for the default 72h Lakes '
                        'Entrance run; N=100 produces ~5.5k files (~25 GB). '
                        'Ignored unless --checkpoint is set.')
    args = p.parse_args()
    main(hurricane=args.hurricane, checkpoint_dir=args.checkpoint,
         checkpoint_stride=args.checkpoint_stride)
