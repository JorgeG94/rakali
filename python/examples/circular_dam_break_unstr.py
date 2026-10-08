"""Circular dam break on an unstructured disk mesh.

A cylindrical column of water of radius R0 and depth h_in sits inside
a flat-bed disk domain of radius R_outer at depth h_out. At t = 0 the
column is released; the resulting wave is a radially-symmetric outward
shock plus an inward-going rarefaction.

The structured-grid solver fundamentally cannot reproduce radial
symmetry exactly because its grid is anisotropic; the unstructured
solver can, to within the limiter's accuracy. This example is a
visual proof: a clean ring of water moving outward.

Inputs are loaded from ``circular_dam_break_unstr.yaml`` next to this script.

Run from the repo root::

    PYTHONPATH=python python python/examples/circular_dam_break_unstr.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from rakali.unstr import UnstructuredMesh, UnstructuredSimulation
from rakali.yaml_loader import load_config


def main():
    cfg = load_config(Path(__file__).with_suffix(".yaml"))

    R_out = cfg.domain.radius
    R_dam = cfg.initial_condition.dam_break.radius
    diameter = 2.0 * R_out
    dx_target = diameter / cfg.domain.n_across
    max_area = 0.5 * dx_target * dx_target

    print(f"Building disk mesh: R={R_out} m, target dx ~ {dx_target:.3f} m")
    mesh = UnstructuredMesh.disk(
        cx=0.0, cy=0.0, radius=R_out,
        max_area=max_area,
        n_boundary=cfg.domain.n_boundary,
        boundary_tag=1,
    )
    print(f"  -> {mesh.nverts} vertices, {mesh.ntri} triangles")

    sim = UnstructuredSimulation(mesh)
    for tag, bc in cfg.boundary_conditions.groups.items():
        sim.set_bc_type(group=tag, bc_type=bc)
    print(f"  -> solver allocated {sim.ncells} owned cells")

    cx = sim.centroid_x.copy()
    cy = sim.centroid_y.copy()
    r = np.sqrt(cx * cx + cy * cy)

    h_in = cfg.initial_condition.dam_break.h_high
    h_out = cfg.initial_condition.dam_break.h_low
    bed = np.full(sim.ncells, cfg.bathymetry.depth, dtype=np.float64)
    stage = np.where(r < R_dam, h_in, h_out).astype(np.float64)

    sim.set_bathymetry(bed)
    sim.set_initial_stage(stage)

    sim.prepare()
    print(f"Advancing to t = {cfg.domain.t_end} s ...")
    sim.advance(cfg.domain.t_end)

    print(f"Done: t = {sim.time:.4f} s, steps = {sim.step_count}, "
          f"final dt = {sim.dt:.4e} s")
    print(f"Mass error: {sim.mass_error:.3e} m^3")

    s = sim.stage.copy()
    h = sim.height
    print(f"Stage range: [{s.min():.4f}, {s.max():.4f}] m")
    print(f"Height range: [{h.min():.4f}, {h.max():.4f}] m")
    print(f"|xmom| max: {np.abs(sim.xmom).max():.4f}")
    print(f"|ymom| max: {np.abs(sim.ymom).max():.4f}")

    # Quick radial-symmetry diagnostic: bin by r, report std/mean per bin.
    r_now = np.sqrt(cx * cx + cy * cy)
    n_bins = 40
    edges = np.linspace(0.0, R_out, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    idx = np.clip(np.digitize(r_now, edges) - 1, 0, n_bins - 1)
    means = np.zeros(n_bins)
    stds = np.zeros(n_bins)
    counts = np.bincount(idx, minlength=n_bins)
    for i in range(n_bins):
        mask = idx == i
        if counts[i] >= 3:
            means[i] = float(s[mask].mean())
            stds[i] = float(s[mask].std())
    nonempty = counts >= 3
    max_relstd = float(np.max(stds[nonempty] /
                              np.maximum(np.abs(means[nonempty] - h_out),
                                         1.0e-6)))
    print(f"max per-bin relative std (radial-symmetry proxy): {max_relstd:.4f}")

    try:
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri
    except ImportError:
        print("(matplotlib not available; skipping plot)")
        return

    # 2D contour of stage at t_end
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    triang = mtri.Triangulation(mesh.vert_x, mesh.vert_y, mesh.tri_v - 1)
    tpc = axes[0].tripcolor(triang, facecolors=s, shading="flat",
                            cmap="viridis")
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel("y (m)")
    axes[0].set_title(f"Stage at t = {sim.time:.3f} s")
    plt.colorbar(tpc, ax=axes[0], label="stage (m)")

    # Radial profile: stage vs r
    axes[1].plot(r_now, s, "k.", ms=1, alpha=0.4, label="all cells")
    axes[1].plot(centers[nonempty], means[nonempty], "r-", lw=1.5,
                 label="bin mean")
    axes[1].axvline(R_dam, color="b", lw=0.5, ls="--", label="initial dam")
    axes[1].set_xlabel("r (m)")
    axes[1].set_ylabel("stage (m)")
    axes[1].set_title("Radial profile")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig("circular_dam_break_unstr.png", dpi=120)
    print("Saved circular_dam_break_unstr.png")


if __name__ == "__main__":
    main()
