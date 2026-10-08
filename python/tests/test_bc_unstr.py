"""Sanity tests for the new BC types: BC_INFLOW, BC_OUTFLOW, BC_TIDAL.

These don't compare against an analytical solution -- they verify that
each BC has the qualitatively correct effect on the flow:

* BC_INFLOW pumps mass and momentum into a closed-on-three-sides domain.
* BC_OUTFLOW lets mass leave at a rate consistent with the prescribed depth.
* BC_INFLOW + BC_OUTFLOW together drive a steady-state channel flow.
* BC_TIDAL pulls the free surface near the boundary toward the
  prescribed external stage.

Steady-state validation against an exact 1D channel solution is left
for the convergence-test file -- this is the smoke / regression set.
"""
from __future__ import annotations

import numpy as np
import pytest

from rakali.unstr import UnstructuredMesh, UnstructuredSimulation


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_channel(nx: int = 50, length: float = 50.0, width: float = 2.0):
    dx = length / nx
    max_area = 0.5 * dx * dx
    return UnstructuredMesh.rectangle(
        x0=0.0, y0=0.0, x1=length, y1=width,
        max_area=max_area,
        bc_tags={"left": 2, "right": 3, "bottom": 1, "top": 1},
    )


def _flat_lake(sim, depth: float):
    bed = np.zeros(sim.ncells, dtype=np.float64)
    stage = np.full(sim.ncells, depth, dtype=np.float64)
    sim.set_bathymetry(bed)
    sim.set_initial_stage(stage)


def _mass(sim) -> float:
    """Total water volume (stage - bed). Using stage avoids the issue
    that ``sim.height`` is zero until the first kernel pass populates
    ``height_cv``.
    """
    return float((sim.stage - sim.bed).sum())


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_inflow_pumps_mass_in():
    """Inflow on the left, walls everywhere else (no outlet) -- mass
    must increase monotonically and xmom near the inflow must be
    positive."""
    mesh = _build_channel(nx=50, length=50.0)
    with UnstructuredSimulation(mesh) as sim:
        sim.set_bc_type(group=1, bc_type="wall")    # long sides
        sim.set_bc_type(group=2, bc_type="inflow")  # left
        sim.set_bc_type(group=3, bc_type="wall")    # right (sealed)
        sim.set_bc_inflow_speed(group=2, speed=0.5)

        _flat_lake(sim, depth=1.0)

        sim.prepare()

        # Track mass over a few intervals
        masses = []
        for t_target in (0.0, 1.0, 2.0):
            if t_target > 0.0:
                sim.advance(t_target)
            masses.append(_mass(sim))

        print(f"\n  masses: {masses}")
        assert masses[1] > masses[0] + 1e-6, "mass did not grow at all"
        assert masses[2] > masses[1] + 1e-6, "mass growth did not continue"

        # Flow direction near the inflow (left tenth of the channel)
        cx = sim.centroid_x.copy()
        left_band = cx < 5.0
        assert np.any(left_band)
        mean_xmom_left = float(sim.xmom[left_band].mean())
        print(f"  mean xmom near inflow: {mean_xmom_left:.4f}")
        assert mean_xmom_left > 0.0, "inflow did not produce positive xmom"


def test_outflow_drains_when_higher_depth_inside():
    """Closed on three sides, outflow on the right with h_ext < interior
    depth -- mass should monotonically decrease."""
    mesh = _build_channel(nx=50, length=50.0)
    with UnstructuredSimulation(mesh) as sim:
        sim.set_bc_type(group=1, bc_type="wall")
        sim.set_bc_type(group=2, bc_type="wall")     # left sealed
        sim.set_bc_type(group=3, bc_type="outflow")  # right
        sim.set_bc_h_ext(group=3, h=0.5)             # downstream depth

        _flat_lake(sim, depth=1.5)
        sim.prepare()

        m0 = _mass(sim)
        sim.advance(2.0)
        m1 = _mass(sim)
        sim.advance(4.0)
        m2 = _mass(sim)

        print(f"\n  masses: {m0:.3f} -> {m1:.3f} -> {m2:.3f}")
        assert m1 < m0 - 1e-6, "outflow did not drain any mass"
        assert m2 < m1 - 1e-6, "drain stalled"
        assert sim.height.min() > 0.0, "produced negative depth"


def test_inflow_outflow_steady_channel():
    """Inflow on the left, outflow on the right, walls on the long
    sides. After enough time the flow should look like a near-steady
    channel: xmom positive everywhere, depth bounded, no NaNs."""
    mesh = _build_channel(nx=80, length=50.0)
    with UnstructuredSimulation(mesh) as sim:
        sim.set_bc_type(group=1, bc_type="wall")
        sim.set_bc_type(group=2, bc_type="inflow")
        sim.set_bc_type(group=3, bc_type="outflow")
        sim.set_bc_inflow_speed(group=2, speed=0.3)
        sim.set_bc_h_ext(group=3, h=1.0)

        _flat_lake(sim, depth=1.0)
        sim.prepare()
        sim.advance(20.0)

        h = sim.height.copy()
        xm = sim.xmom.copy()

        print(f"\n  steps={sim.step_count}, t={sim.time:.2f}")
        print(f"  h range: [{h.min():.4f}, {h.max():.4f}]")
        print(f"  xmom mean: {xm.mean():.4f}, max: {xm.max():.4f}")

        assert not np.any(np.isnan(h)), "NaN in depth"
        assert h.min() > 0.0, "negative depth"
        assert h.max() < 2.0, "depth blew up"
        # Mean discharge should be positive (left -> right)
        assert xm.mean() > 0.0, "no net flow developed"
        # The mean xmom should be roughly h * u_inflow ~ 0.3
        assert 0.1 < xm.mean() < 0.6, (
            f"discharge {xm.mean():.3f} far from expected ~0.3"
        )


def test_tidal_clamps_stage_near_boundary():
    """One end of a closed basin set to BC_TIDAL with eta_ext = 1.5,
    interior at h = 1.0. The free surface near the boundary should
    rise toward 1.5 and the basin mass should grow."""
    mesh = _build_channel(nx=50, length=50.0)
    with UnstructuredSimulation(mesh) as sim:
        sim.set_bc_type(group=1, bc_type="wall")
        sim.set_bc_type(group=2, bc_type="tidal")    # left
        sim.set_bc_type(group=3, bc_type="wall")     # right sealed
        sim.set_bc_eta_ext(group=2, eta=1.5)

        _flat_lake(sim, depth=1.0)
        sim.prepare()

        m0 = _mass(sim)
        sim.advance(3.0)
        m1 = _mass(sim)

        cx = sim.centroid_x.copy()
        near_band = cx < 5.0
        far_band = cx > 45.0
        h = sim.height
        assert np.any(near_band) and np.any(far_band)

        print(f"\n  mass {m0:.3f} -> {m1:.3f}")
        print(f"  h near tidal: mean={h[near_band].mean():.4f}")
        print(f"  h far from tidal: mean={h[far_band].mean():.4f}")

        assert m1 > m0 + 1e-6, "tidal forcing did not raise mass"
        # Cells near the boundary should have moved toward eta_ext
        assert h[near_band].mean() > 1.05, (
            "stage near tidal boundary did not rise"
        )
