#!/usr/bin/env python
"""Tidal simulation on the Parramatta River using the unstructured solver.

Cuts the river + upper harbour from the 25m Sydney Harbour bathymetry,
builds a triangulated mesh via :func:`rakali.mesh.build_coastal_mesh`,
drives a tidal simulation for a few hours, caches snapshots to a
NetCDF file via xarray, and renders a GIF.

Usage:
    PYTHONPATH=python python python/examples/real_bathymetry/sydney_basin/parramatta/run_river.py
"""
from __future__ import annotations

import sys
import time as walltime
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[5]
for p in [str(_REPO / "python"), str(_REPO)]:
    if p not in sys.path:
        sys.path.insert(0, p)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import xarray as xr

from scipy.ndimage import gaussian_filter

from rakali.mesh import build_coastal_mesh
from rakali.mesh.bed import load_bathymetry_nc
from rakali.unstr import set_log_level, LOG_INFO

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

BATHY_FILE = Path(__file__).resolve().parents[1] / "bathymetry_25m_merged.nc"
OUT_DIR = Path(__file__).resolve().parent
NC_FILE = OUT_DIR / "parramatta_snapshots_flood.nc"
GIF_FILE = OUT_DIR / "parramatta_tidal_flood.gif"
DX = 25.0  # m

WINDOW = (350 * DX, 1350 * DX, 400 * DX, 1000 * DX)
MESH_RESOLUTION = 10 * DX  # ~50 m triangle edge

# Soften the bathymetry so the river banks slope gently into a low
# flood plain instead of into cliffs:
#   1. Gaussian-smooth the raster (sigma = ~150 m) to round off steep
#      walls so the contour at FLOOD_SEA_LEVEL produces a graded bank.
#   2. Cap the maximum land elevation so the highest land near the
#      coast is still within reach of the surge.
# This expands the mesh to include initially-dry coastal terrain, which
# the +3 m surge then overtops.
BATHY_SMOOTH_SIGMA = 6.0   # raster cells (dx=25 m → ~150 m)
BATHY_CAP_MAX      = 3.0   # m, max land elevation after smoothing
FLOOD_SEA_LEVEL    = 1.5   # m, contour used for the coastline

# View window for the GIF (km). Zoom to the upper river so the flooding
# fills the frame.
VIEW_XLIM = (10.0, 25.0)
VIEW_YLIM = (10.0, 25.0)

# Tidal constituents (M2, S2, K1, O1) — amplified ~3x relative to the
# real Sydney Harbour signal to drive flooding.
TIDAL_AMP   = [1.500, 0.400, 0.500, 0.350]
TIDAL_PHASE = [-2.681322, -0.676741, -2.162299, 1.963210]
TIDAL_OMEGA = [1.405189e-04, 1.454441e-04, 7.292117e-05, 6.759774e-05]

# Storm surge that ramps in (cosine taper) over SURGE_RAMP_S then stays
# elevated. Stacks on top of the tide at the open boundary; gives the
# current the extra oomph that pushes water onto the flood plain.
SURGE_PEAK   = 1.5            # m
SURGE_RAMP_S = 6 * 3600.0     # 6-hour ramp-up

T_END = 12 * 3600.0   # 24 hours
DT_BC = 60.0          # update tidal BC every 60 s
SNAP_INTERVAL = 900.0 # snapshot every 15 minutes
MANNING_N = 0.025
CFL = 0.45


# ---------------------------------------------------------------------------
# Tidal elevation
# ---------------------------------------------------------------------------

def tidal_eta(t: float) -> float:
    eta = 0.0
    for a, p, w in zip(TIDAL_AMP, TIDAL_PHASE, TIDAL_OMEGA):
        eta += a * np.cos(w * t + p)
    ramp = 0.5 * (1.0 - np.cos(np.pi * min(t / SURGE_RAMP_S, 1.0)))
    eta += SURGE_PEAK * ramp
    return float(eta)


# ---------------------------------------------------------------------------
# Run simulation → xarray Dataset
# ---------------------------------------------------------------------------

def run_simulation(result):
    """Run the tidal simulation. Returns an xarray Dataset of snapshots."""
    mesh = result.mesh
    bed = result.bed
    ncells = mesh.ntri

    cx, cy = mesh.centroids()
    # Mesh vertices / triangles (0-based) for plotting
    new_vx = mesh.vert_x
    new_vy = mesh.vert_y
    new_tri = (mesh.tri_v - 1).astype(np.int32)

    snap_times = np.arange(SNAP_INTERVAL, T_END + 1, SNAP_INTERVAL)
    n_snaps = len(snap_times)
    eta_all = np.full((n_snaps, ncells), np.nan, dtype=np.float32)
    eta_tide_all = np.zeros(n_snaps, dtype=np.float32)

    with result.make_simulation(cfl=CFL, manning=MANNING_N) as sim:
        print(f"Running {T_END/3600:.0f}-hour tidal simulation...")
        t0 = walltime.time()
        t = 0.0
        snap_idx = 0

        while t < T_END:
            t_next = min(t + DT_BC, T_END)
            eta_bc = tidal_eta(t_next)
            sim.set_bc_eta_ext(group=result.open_marker, eta=eta_bc)
            sim.advance(t_next)
            t = t_next

            if int(t) % 600 < DT_BC:
                print(f"  t={t/3600:.1f}h  eta_tide={eta_bc:+.3f}m  "
                      f"dt={sim.dt:.3f}s  steps={sim.step_count}")

            if snap_idx < n_snaps and t >= snap_times[snap_idx]:
                stage = sim.stage.copy()
                height = sim.height.copy()

                eta = stage.copy()
                eta[height < 0.01] = np.nan

                eta_all[snap_idx] = eta.astype(np.float32)
                eta_tide_all[snap_idx] = eta_bc
                snap_idx += 1

        elapsed = walltime.time() - t0
        print(f"\nDone: {sim.step_count} steps in {elapsed:.1f}s "
              f"({sim.step_count * ncells / elapsed / 1e6:.1f} Mcells/s)")

    ds = xr.Dataset(
        {
            "eta":   (["time", "cell"], eta_all),
            "bed":   (["cell"], bed.astype(np.float32)),
            "cx":    (["cell"], cx.astype(np.float32)),
            "cy":    (["cell"], cy.astype(np.float32)),
            "tri":   (["tri_idx", "vertex"], new_tri.astype(np.int32)),
            "vx":    (["vert"], new_vx.astype(np.float32)),
            "vy":    (["vert"], new_vy.astype(np.float32)),
        },
        coords={
            "time": snap_times / 3600.0,  # hours
            "eta_tide": ("time", eta_tide_all),
        },
        attrs={
            "description": "Parramatta River tidal simulation snapshots",
            "resolution": MESH_RESOLUTION,
            "dx": DX,
            "manning_n": MANNING_N,
        },
    )
    ds["time"].attrs["units"] = "hours"
    ds["eta"].attrs["units"] = "m"
    ds["bed"].attrs["units"] = "m"
    return ds


# ---------------------------------------------------------------------------
# GIF rendering
# ---------------------------------------------------------------------------

def render_gif(ds: xr.Dataset):
    """Render a GIF from the cached xarray Dataset."""
    from matplotlib.animation import FuncAnimation, PillowWriter

    vx_km = ds["vx"].values / 1000.0
    vy_km = ds["vy"].values / 1000.0
    triangles = ds["tri"].values
    tri_plot = mtri.Triangulation(vx_km, vy_km, triangles)

    times = ds["time"].values
    eta_tide = ds["eta_tide"].values
    eta_all = ds["eta"].values

    fig, ax_eta = plt.subplots(figsize=(14, 12))
    # Light-gray background so initially-dry cells (eta = NaN) are
    # visible as land — flooding then "fills in" the gray.
    ax_eta.set_facecolor("0.82")

    tcf_eta = ax_eta.tripcolor(tri_plot, eta_all[0], cmap="RdBu_r",
                                vmin=-3.0, vmax=3.0, shading="flat")
    plt.colorbar(tcf_eta, ax=ax_eta, label="Surface elevation (m)", shrink=0.85)
    ax_eta.set_xlabel("x (km)")
    ax_eta.set_ylabel("y (km)")
    ax_eta.set_xlim(*VIEW_XLIM)
    ax_eta.set_ylim(*VIEW_YLIM)
    ax_eta.set_aspect("equal")

    title = fig.suptitle("", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    def update(frame):
        tcf_eta.set_array(eta_all[frame])
        title.set_text(
            f"Parramatta River — t = {times[frame]:.1f}h, "
            f"tide+surge = {eta_tide[frame]:+.3f}m"
        )
        return tcf_eta, title

    anim = FuncAnimation(fig, update, frames=len(times), blit=False)
    anim.save(str(GIF_FILE), writer=PillowWriter(fps=8), dpi=90)
    plt.close(fig)
    print(f"  Saved: {GIF_FILE}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    set_log_level(LOG_INFO)

    bathy_raw = load_bathymetry_nc(BATHY_FILE)
    bathy_soft = np.minimum(
        gaussian_filter(bathy_raw, sigma=BATHY_SMOOTH_SIGMA),
        BATHY_CAP_MAX,
    )
    print(f"Bathy softened: sigma={BATHY_SMOOTH_SIGMA} cells, "
          f"cap={BATHY_CAP_MAX} m, range now "
          f"[{bathy_soft.min():.1f}, {bathy_soft.max():.1f}]")

    result = build_coastal_mesh(
        bathy_soft,
        resolution=MESH_RESOLUTION,
        dx=DX, window=WINDOW,
        sea_level=FLOOD_SEA_LEVEL,
        open_boundary="east", force_wet_edges=("east",),
    )
    print(f"\nMesh ready: {result.mesh.ntri} triangles")
    print(f"Bed range: [{result.bed.min():.1f}, {result.bed.max():.1f}] m")
    print(f"Tidal eta at t=0: {tidal_eta(0.0):.3f} m")
    print(f"Tidal eta at t={T_END/3600:.0f}h: {tidal_eta(T_END):.3f} m\n")

    wall_total = walltime.time()

    if NC_FILE.exists():
        print(f"Loading cached snapshots from {NC_FILE}")
        ds = xr.open_dataset(NC_FILE)
    else:
        ds = run_simulation(result)
        ds.to_netcdf(NC_FILE)
        print(f"  Cached snapshots to {NC_FILE}")

    wall_sim = walltime.time() - wall_total

    if "--no-gif" not in sys.argv:
        print("Rendering GIF...")
        wall_gif = walltime.time()
        render_gif(ds)
        wall_gif = walltime.time() - wall_gif
        wall_total = walltime.time() - wall_total
        print(f"\nWall time: sim={wall_sim:.1f}s  gif={wall_gif:.1f}s  total={wall_total:.1f}s")
    else:
        wall_total = walltime.time() - wall_total
        print(f"\nWall time: sim={wall_sim:.1f}s  total={wall_total:.1f}s (gif skipped)")


if __name__ == "__main__":
    main()
