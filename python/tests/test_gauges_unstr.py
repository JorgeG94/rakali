"""Test gauge output for the unstructured solver.

Places gauges at known locations on a simple mesh, runs a short
simulation, and verifies the output file contains sensible data.
"""
from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

from rakali.unstr import UnstructuredMesh, UnstructuredSimulation

netCDF4 = pytest.importorskip("netCDF4")


def _build_sim():
    mesh = UnstructuredMesh.rectangle(
        x0=0.0, y0=0.0, x1=100.0, y1=10.0,
        max_area=1.0,
        bc_tags={"left": 1, "right": 1, "bottom": 1, "top": 1},
    )
    sim = UnstructuredSimulation(mesh)
    sim.set_bc_type(group=1, bc_type="wall")

    bed = np.zeros(sim.ncells, dtype=np.float64)
    stage = np.full(sim.ncells, 5.0, dtype=np.float64)
    sim.set_bathymetry(bed)
    sim.set_initial_stage(stage)
    return sim, mesh


def test_gauge_file_structure():
    """Gauge output file has correct dimensions and variables."""
    with tempfile.TemporaryDirectory() as tmpdir:
        gauge_file = os.path.join(tmpdir, "gauges.nc")

        sim, _ = _build_sim()
        sim.prepare()

        # Place 3 gauges inside the domain
        gx = np.array([25.0, 50.0, 75.0])
        gy = np.array([5.0, 5.0, 5.0])
        sim.init_gauges(gauge_file, gx, gy)

        sim.advance(0.01)
        sim.write_gauges()
        sim.advance(0.02)
        sim.write_gauges()

        sim.close_gauges()
        sim.destroy()

        ds = netCDF4.Dataset(gauge_file, "r")

        assert "station" in ds.dimensions
        assert "time" in ds.dimensions
        assert ds.dimensions["station"].size == 3
        assert ds.dimensions["time"].size == 2

        assert "station_x" in ds.variables
        assert "station_y" in ds.variables
        assert "eta" in ds.variables
        assert "u" in ds.variables
        assert "v" in ds.variables
        assert "time" in ds.variables

        np.testing.assert_allclose(ds.variables["station_x"][:], gx, atol=1e-10)
        np.testing.assert_allclose(ds.variables["station_y"][:], gy, atol=1e-10)

        ds.close()


def test_gauge_values_lake_at_rest():
    """In a lake at rest, eta == initial stage and velocities == 0."""
    with tempfile.TemporaryDirectory() as tmpdir:
        gauge_file = os.path.join(tmpdir, "gauges.nc")

        sim, _ = _build_sim()
        sim.prepare()

        gx = np.array([30.0, 60.0])
        gy = np.array([5.0, 5.0])
        sim.init_gauges(gauge_file, gx, gy)

        # Lake at rest — advance a bit
        sim.advance(0.5)
        sim.write_gauges()
        sim.close_gauges()
        sim.destroy()

        ds = netCDF4.Dataset(gauge_file, "r")

        eta = ds.variables["eta"][0, :]
        u = ds.variables["u"][0, :]
        v = ds.variables["v"][0, :]

        # Stage should be exactly the initial 5.0 (lake at rest)
        np.testing.assert_allclose(eta, 5.0, atol=1e-10)
        # Velocities should be zero
        np.testing.assert_allclose(u, 0.0, atol=1e-10)
        np.testing.assert_allclose(v, 0.0, atol=1e-10)

        ds.close()


def test_gauge_outside_domain_skipped():
    """Gauges outside the domain don't crash; they're silently skipped."""
    with tempfile.TemporaryDirectory() as tmpdir:
        gauge_file = os.path.join(tmpdir, "gauges.nc")

        sim, _ = _build_sim()
        sim.prepare()

        # One inside, one way outside
        gx = np.array([50.0, 999.0])
        gy = np.array([5.0, 999.0])
        sim.init_gauges(gauge_file, gx, gy)

        sim.advance(0.01)
        sim.write_gauges()
        sim.close_gauges()
        sim.destroy()

        ds = netCDF4.Dataset(gauge_file, "r")
        # Only 1 station should be in the file
        assert ds.dimensions["station"].size == 1
        np.testing.assert_allclose(ds.variables["station_x"][0], 50.0, atol=1e-10)
        ds.close()


def test_destroy_closes_gauges():
    """Calling destroy auto-closes the gauge file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        gauge_file = os.path.join(tmpdir, "gauges.nc")

        sim, _ = _build_sim()
        sim.prepare()

        gx = np.array([50.0])
        gy = np.array([5.0])
        sim.init_gauges(gauge_file, gx, gy)

        sim.advance(0.01)
        sim.write_gauges()
        sim.destroy()  # should close gauges too

        ds = netCDF4.Dataset(gauge_file, "r")
        assert ds.dimensions["time"].size == 1
        ds.close()
