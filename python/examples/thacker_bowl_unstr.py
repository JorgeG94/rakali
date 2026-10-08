"""Sloshing in a parabolic bowl on an unstructured disk mesh.

Demonstrates the wet/dry treatment: a tilted lake in a parabolic bowl
oscillates back and forth, with cells near the moving shoreline
alternately drying and re-wetting.

Inputs are loaded from ``thacker_bowl_unstr.yaml`` next to this script.

Run from the repo root::

    PYTHONPATH=python python python/examples/thacker_bowl_unstr.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from rakali.unstr import UnstructuredMesh, UnstructuredSimulation
from rakali.yaml_loader import load_config


GRAVITY = 9.81


def _bed(r: np.ndarray, L: float, h0: float) -> np.ndarray:
    return -h0 * (1.0 - (r * r) / (L * L))


def main():
    cfg = load_config(Path(__file__).with_suffix(".yaml"))

    L = cfg.bathymetry.bowl.L
    h0 = cfg.bathymetry.bowl.h0
    z_lake = cfg.initial_condition.tilt.z_lake
    tilt = cfg.initial_condition.tilt.amplitude

    R_dom = cfg.domain.radius
    diameter = 2.0 * R_dom
    dx = diameter / cfg.domain.n_across
    max_area = 0.5 * dx * dx

    omega = math.sqrt(2.0 * GRAVITY * h0) / L
    period = 2.0 * math.pi / omega
    print(f"Oscillation period (linear theory): {period:.3f} s")

    print(f"Building disk mesh: R={R_dom} m")
    mesh = UnstructuredMesh.disk(
        cx=0.0, cy=0.0, radius=R_dom,
        max_area=max_area,
        n_boundary=cfg.domain.n_boundary,
        boundary_tag=1,
    )
    print(f"  -> {mesh.nverts} vertices, {mesh.ntri} triangles")

    sim = UnstructuredSimulation(mesh)
    for tag, bc in cfg.boundary_conditions.groups.items():
        sim.set_bc_type(group=tag, bc_type=bc)

    cx = sim.centroid_x.copy()
    cy = sim.centroid_y.copy()
    r = np.sqrt(cx * cx + cy * cy)
    bed = _bed(r, L, h0)
    stage = np.maximum(z_lake + tilt * (cx / L), bed)

    sim.set_bathymetry(bed.astype(np.float64))
    sim.set_initial_stage(stage.astype(np.float64))
    sim.prepare()

    area = sim.area.copy()

    snapshots = []
    snapshots.append(("t=0", sim.stage.copy()))
    for label, frac in (("T/4", 0.25), ("T/2", 0.5)):
        sim.advance(frac * period)
        snapshots.append((f"t={label}", sim.stage.copy()))
        h = np.maximum(sim.stage - bed, 0.0)
        m = float((h * area).sum())
        n_dry = int((h < 1.0e-6).sum())
        print(f"  {label}: t={sim.time:.3f}, mass={m:.4f}, "
              f"dry={n_dry}/{sim.ncells}")

    try:
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri
    except ImportError:
        print("(matplotlib not available)")
        return

    triang = mtri.Triangulation(mesh.vert_x, mesh.vert_y, mesh.tri_v - 1)

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    axes = axes.ravel()

    bath = axes[0].tripcolor(triang, facecolors=bed, shading="flat",
                             cmap="terrain")
    axes[0].set_title("Bathymetry")
    axes[0].set_aspect("equal")
    plt.colorbar(bath, ax=axes[0], label="bed (m)")

    for ax, (label, s) in zip(axes[1:], snapshots):
        # Show stage but mask dry cells (bed >= stage) for clarity
        h = np.maximum(s - bed, 0.0)
        s_show = np.where(h > 1.0e-6, s, np.nan)
        tpc = ax.tripcolor(triang, facecolors=s_show, shading="flat",
                           cmap="viridis")
        ax.set_title(label)
        ax.set_aspect("equal")
        plt.colorbar(tpc, ax=ax, label="stage (m)")

    plt.tight_layout()
    plt.savefig("thacker_bowl_unstr.png", dpi=120)
    print("Saved thacker_bowl_unstr.png")


if __name__ == "__main__":
    main()
