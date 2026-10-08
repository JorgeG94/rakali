"""End-to-end Python smoke test of the unstructured Rakali solver.

Builds a triangulated channel with meshpy.triangle, sets up a Ritter
wet/wet dam break, runs the solver, and prints summary statistics.
Optionally plots the free-surface profile if matplotlib is available.

Inputs are loaded from ``dam_break_unstr.yaml`` next to this script.

Run from the repo root after building the shared library::

    cmake -B build_nvhpc -S . -DRAKALI_ENABLE_PYTHON=ON
    cmake --build build_nvhpc
    PYTHONPATH=python python python/examples/dam_break_unstr.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from rakali.unstr import UnstructuredMesh, UnstructuredSimulation
from rakali.yaml_loader import load_config


def main():
    cfg = load_config(Path(__file__).with_suffix(".yaml"))

    L = cfg.domain.length
    W = cfg.domain.width
    dx_target = L / cfg.domain.n_across
    max_area = 0.5 * dx_target * dx_target  # ~ right-isoceles triangle area

    print(f"Building channel mesh: {L} x {W} m, target dx ~ {dx_target:.3f} m")
    mesh = UnstructuredMesh.rectangle(
        x0=0.0, y0=0.0, x1=L, y1=W,
        max_area=max_area,
        bc_tags={"left": 2, "right": 2, "bottom": 1, "top": 1},
    )
    print(f"  -> {mesh.nverts} vertices, {mesh.ntri} triangles")

    sim = UnstructuredSimulation(mesh)
    print(f"  -> solver allocated {sim.ncells} owned cells")

    for tag, bc in cfg.boundary_conditions.groups.items():
        sim.set_bc_type(group=tag, bc_type=bc)

    cx = sim.centroid_x.copy()  # copy because we're about to call setters
    bed = np.full(sim.ncells, cfg.bathymetry.depth, dtype=np.float64)
    h_left = cfg.initial_condition.dam_break.h_high
    h_right = cfg.initial_condition.dam_break.h_low
    split = cfg.initial_condition.dam_break.split_fraction * L
    stage = np.where(cx < split, h_left, h_right)

    sim.set_bathymetry(bed)
    sim.set_initial_stage(stage)

    sim.prepare()
    print(f"Advancing to t = {cfg.domain.t_end} s ...")
    sim.advance(cfg.domain.t_end)

    print(f"Done: t = {sim.time:.6f} s, steps = {sim.step_count}, "
          f"final dt = {sim.dt:.4e} s")
    print(f"Mass error (cumulative protect clip): {sim.mass_error:.3e} m^3")

    h = sim.height
    print(f"Stage range: [{sim.stage.min():.4f}, {sim.stage.max():.4f}] m")
    print(f"Height range: [{h.min():.4f}, {h.max():.4f}] m")
    print(f"|xmom| max: {np.abs(sim.xmom).max():.4f} m^2/s")
    print(f"|ymom| max: {np.abs(sim.ymom).max():.4f} m^2/s")

    # Sanity checks
    assert h.min() > 0.0, "negative depth produced"
    assert sim.stage.max() <= h_left + 1e-6, "stage exceeded h_left"
    assert sim.stage.min() >= h_right - 0.5, "stage dropped well below h_right"

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("(matplotlib not available; skipping plot)")
        return
    order = np.argsort(cx)
    plt.figure(figsize=(8, 4))
    plt.plot(cx[order], sim.stage[order], "k.", ms=2, label="stage (Rakali)")
    plt.axvline(split, color="r", lw=0.5, ls="--", label="initial dam")
    plt.xlabel("x (m)")
    plt.ylabel("stage (m)")
    plt.title(f"Unstructured Ritter dam break, t = {sim.time:.3f} s")
    plt.legend()
    plt.tight_layout()
    plt.savefig("dam_break_unstr.png", dpi=120)
    print("Saved dam_break_unstr.png")


if __name__ == "__main__":
    main()
