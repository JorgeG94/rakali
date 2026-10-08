"""Round-trip test for unstructured solver restart files.

Runs a dam break, saves a restart, creates a new simulation from the
restart, and verifies fields match.
"""
from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

from rakali.unstr import UnstructuredMesh, UnstructuredSimulation


def _build_sim():
    mesh = UnstructuredMesh.rectangle(
        x0=0.0, y0=0.0, x1=100.0, y1=2.0,
        max_area=0.5,
        bc_tags={"left": 2, "right": 2, "bottom": 1, "top": 1},
    )
    sim = UnstructuredSimulation(mesh)
    sim.set_bc_type(group=1, bc_type="wall")
    sim.set_bc_type(group=2, bc_type="open")

    cx = sim.centroid_x.copy()
    bed = np.zeros(sim.ncells, dtype=np.float64)
    stage = np.where(cx < 50.0, 2.0, 1.0).astype(np.float64)
    sim.set_bathymetry(bed)
    sim.set_initial_stage(stage)
    return sim, mesh


def test_restart_round_trip():
    """Save/load preserves stage, xmom, ymom, bed, and time."""
    with tempfile.TemporaryDirectory() as tmpdir:
        restart_file = os.path.join(tmpdir, "restart.nc")

        # Run and save
        sim1, mesh = _build_sim()
        sim1.prepare()
        sim1.advance(0.1)

        stage_ref = sim1.stage.copy()
        xmom_ref = sim1.xmom.copy()
        ymom_ref = sim1.ymom.copy()
        bed_ref = sim1.bed.copy()
        time_ref = sim1.time

        sim1.save_restart(restart_file)
        sim1.destroy()

        # Load into a new simulation with the same mesh
        sim2 = UnstructuredSimulation(mesh)
        sim2.load_restart(restart_file)
        sim2.prepare()

        np.testing.assert_allclose(sim2.stage, stage_ref, atol=1e-14)
        np.testing.assert_allclose(sim2.xmom, xmom_ref, atol=1e-14)
        np.testing.assert_allclose(sim2.ymom, ymom_ref, atol=1e-14)
        np.testing.assert_allclose(sim2.bed, bed_ref, atol=1e-14)
        assert sim2.time == pytest.approx(time_ref, abs=1e-12)

        sim2.destroy()


def test_restart_continue_simulation():
    """A restarted simulation continues and produces the same result
    as running straight through."""
    with tempfile.TemporaryDirectory() as tmpdir:
        restart_file = os.path.join(tmpdir, "restart.nc")

        # Run straight through
        sim_ref, mesh = _build_sim()
        sim_ref.prepare()
        sim_ref.advance(0.2)
        stage_full = sim_ref.stage.copy()
        sim_ref.destroy()

        # Run to midpoint, restart, continue
        sim1, mesh = _build_sim()
        sim1.prepare()
        sim1.advance(0.1)
        sim1.save_restart(restart_file)
        sim1.destroy()

        sim2 = UnstructuredSimulation(mesh)
        sim2.load_restart(restart_file)
        sim2.set_bc_type(group=1, bc_type="wall")
        sim2.set_bc_type(group=2, bc_type="open")
        sim2.prepare()
        sim2.advance(0.2)
        stage_restarted = sim2.stage.copy()
        sim2.destroy()

        # Not bitwise — the restarted run re-bootstraps the CFL (one
        # extra extrapolate+flux pass) so the first dt differs slightly
        # and errors compound. 1e-3 is generous but correct.
        np.testing.assert_allclose(stage_restarted, stage_full, atol=1e-3)
