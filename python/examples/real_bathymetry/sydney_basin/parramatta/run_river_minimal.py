#!/usr/bin/env python
"""Minimal Parramatta River tidal simulation.

This is the shortest end-to-end coastal sim you can write with Rakali:
read bathymetry → build mesh → step.  Defaults (CFL, Manning, wall/open
BCs, still-water IC) all come from :meth:`CoastalMeshResult.make_simulation`.

Usage:
    PYTHONPATH=python python python/examples/real_bathymetry/sydney_basin/parramatta/run_river_minimal.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(_REPO / "python"))

import matplotlib
matplotlib.use("Agg")

from rakali.mesh import build_coastal_mesh

BATHY_FILE = Path(__file__).resolve().parents[1] / "bathymetry_25m_merged.nc"
DX = 25.0

# M2 tide — the dominant semidiurnal component on the NSW coast.
TIDAL_AMP, TIDAL_PHASE, TIDAL_OMEGA = 0.49, -2.68, 1.405e-04


def main():
    result = build_coastal_mesh(
        BATHY_FILE,
        resolution=2 * DX,
        dx=DX,
        window=(350 * DX, 1350 * DX, 400 * DX, 1000 * DX),
        open_boundary="east",
        force_wet_edges=("east",),
    )
    print(f"Mesh: {result.mesh.ntri} triangles, {len(result.islands)} islands")

    with result.make_simulation() as sim:
        t_end = 3600.0
        dt_bc = 60.0
        t = 0.0
        while t < t_end:
            t_next = min(t + dt_bc, t_end)
            eta = TIDAL_AMP * np.cos(TIDAL_OMEGA * t_next + TIDAL_PHASE)
            sim.set_bc_eta_ext(group=result.open_marker, eta=eta)
            sim.advance(t_next)
            t = t_next

        print(f"Done: t={sim.time/3600:.2f}h, steps={sim.step_count}")


if __name__ == "__main__":
    main()
