"""Round-trip test for UGRID NetCDF output from the unstructured solver.

Runs a short dam break, writes snapshots via ``sim.write_output()``,
then reads the file back with netCDF4 and checks:

1. UGRID topology: mesh2d attributes, node coords, face_node_connectivity
2. Static fields: bathymetry, centroids, areas
3. Time-varying fields: stage/xmom/ymom match the solver state
4. Multiple snapshots accumulate correctly
5. Deflate compression option works
"""
from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

from rakali.unstr import UnstructuredMesh, UnstructuredSimulation

netCDF4 = pytest.importorskip("netCDF4")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_sim(nx: int = 50, length: float = 100.0, width: float = 2.0):
    dx = length / nx
    max_area = 0.5 * dx * dx
    mesh = UnstructuredMesh.rectangle(
        x0=0.0, y0=0.0, x1=length, y1=width,
        max_area=max_area,
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


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_ugrid_topology():
    """mesh2d attributes, node coords, and face_node_connectivity are correct."""
    with tempfile.TemporaryDirectory() as tmpdir:
        outfile = os.path.join(tmpdir, "topo.nc")

        sim, mesh = _build_sim()
        sim.prepare()
        sim.init_output(outfile)
        sim.write_output()  # one snapshot so the file is valid
        sim.close_output()

        ncells = sim.ncells
        nverts = mesh.nverts
        sim.destroy()

        ds = netCDF4.Dataset(outfile, "r")

        # mesh2d topology variable
        assert "mesh2d" in ds.variables
        m = ds.variables["mesh2d"]
        assert m.getncattr("cf_role") == "mesh_topology"
        assert m.getncattr("topology_dimension") == 2

        # node coordinates
        assert "mesh_node_x" in ds.variables
        assert "mesh_node_y" in ds.variables
        nx = ds.variables["mesh_node_x"][:]
        ny = ds.variables["mesh_node_y"][:]
        assert nx.shape == ny.shape
        assert nx.shape[0] == nverts

        # face-node connectivity
        # Fortran writes (three, face); Python/C reads as (face, three)
        assert "mesh_face_nodes" in ds.variables
        conn = ds.variables["mesh_face_nodes"][:]
        assert conn.shape == (ncells, 3)
        assert ds.variables["mesh_face_nodes"].getncattr("start_index") == 1
        # all vertex indices should be in [1, nverts]
        assert conn.min() >= 1
        assert conn.max() <= nverts

        ds.close()


def test_static_fields():
    """Bathymetry, centroids, and areas match the solver values."""
    with tempfile.TemporaryDirectory() as tmpdir:
        outfile = os.path.join(tmpdir, "static.nc")

        sim, _mesh = _build_sim()
        sim.prepare()
        sim.init_output(outfile)
        sim.write_output()
        sim.close_output()

        bed_ref = sim.bed.copy()
        cx_ref = sim.centroid_x.copy()
        cy_ref = sim.centroid_y.copy()
        area_ref = sim.area.copy()
        sim.destroy()

        ds = netCDF4.Dataset(outfile, "r")

        np.testing.assert_allclose(ds.variables["bed"][:], bed_ref, atol=1e-14)
        np.testing.assert_allclose(ds.variables["mesh_face_x"][:], cx_ref, atol=1e-14)
        np.testing.assert_allclose(ds.variables["mesh_face_y"][:], cy_ref, atol=1e-14)
        np.testing.assert_allclose(ds.variables["mesh_face_area"][:], area_ref,
                                   atol=1e-14)

        ds.close()


def test_snapshot_values():
    """Stage, xmom, ymom at a single snapshot match the solver state."""
    with tempfile.TemporaryDirectory() as tmpdir:
        outfile = os.path.join(tmpdir, "snap.nc")

        sim, _mesh = _build_sim()
        sim.prepare()
        sim.advance(0.05)

        sim.init_output(outfile)
        sim.write_output()
        sim.close_output()

        # After write_output the GPU mapping is still intact,
        # but reading field views calls ensure_host. Copy now.
        stage_ref = sim.stage.copy()
        xmom_ref = sim.xmom.copy()
        ymom_ref = sim.ymom.copy()
        t_ref = sim.time
        sim.destroy()

        ds = netCDF4.Dataset(outfile, "r")

        assert ds.dimensions["time"].size == 1
        np.testing.assert_allclose(ds.variables["time"][0], t_ref, atol=1e-12)
        np.testing.assert_allclose(ds.variables["stage"][0, :], stage_ref, atol=1e-14)
        np.testing.assert_allclose(ds.variables["xmom"][0, :], xmom_ref, atol=1e-14)
        np.testing.assert_allclose(ds.variables["ymom"][0, :], ymom_ref, atol=1e-14)

        ds.close()


def test_multiple_snapshots():
    """Writing multiple snapshots accumulates the time dimension."""
    n_snaps = 3
    dt_out = 0.02

    with tempfile.TemporaryDirectory() as tmpdir:
        outfile = os.path.join(tmpdir, "multi.nc")

        sim, _mesh = _build_sim()
        sim.prepare()
        sim.init_output(outfile)
        times = []
        for i in range(n_snaps):
            sim.advance((i + 1) * dt_out)
            sim.write_output()
            times.append(sim.time)
        sim.close_output()
        sim.destroy()

        ds = netCDF4.Dataset(outfile, "r")

        assert ds.dimensions["time"].size == n_snaps
        for i, t_expected in enumerate(times):
            np.testing.assert_allclose(
                ds.variables["time"][i], t_expected, atol=1e-12,
                err_msg=f"time mismatch at snapshot {i}",
            )
        # stage should vary between snapshots (dam break is evolving)
        s0 = ds.variables["stage"][0, :]
        s_last = ds.variables["stage"][n_snaps - 1, :]
        assert not np.allclose(s0, s_last), "snapshots are identical — not evolving"

        ds.close()


def test_stepping_continues_after_write():
    """Solver can keep stepping after write_output without re-finalize."""
    with tempfile.TemporaryDirectory() as tmpdir:
        outfile = os.path.join(tmpdir, "continue.nc")

        sim, _mesh = _build_sim()
        sim.prepare()
        sim.init_output(outfile)

        sim.advance(0.01)
        sim.write_output()
        t1 = sim.time

        # This should NOT require a manual re-finalize
        sim.advance(0.02)
        t2 = sim.time
        sim.write_output()
        sim.close_output()
        sim.destroy()

        assert t2 > t1, "time did not advance after write_output"

        ds = netCDF4.Dataset(outfile, "r")
        assert ds.dimensions["time"].size == 2
        assert ds.variables["time"][1] > ds.variables["time"][0]
        ds.close()


def test_deflate_compression():
    """Output with deflate_level > 0 produces a valid, compressed file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        f_deflate = os.path.join(tmpdir, "deflate.nc")

        sim, _mesh = _build_sim()
        sim.prepare()
        sim.init_output(f_deflate, deflate_level=1)
        sim.advance(0.05)
        sim.write_output()
        sim.close_output()
        sim.destroy()

        # Verify deflated file is readable and variables have compression
        ds = netCDF4.Dataset(f_deflate, "r")
        assert ds.dimensions["time"].size == 1
        assert ds.variables["stage"][:].shape[1] > 0

        # Check that at least one variable has deflate enabled
        filt = ds.variables["stage"].filters()
        assert filt["zlib"] is True, "deflate not enabled on stage variable"

        ds.close()


def test_destroy_closes_output():
    """Calling destroy() on the sim auto-closes the output file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        outfile = os.path.join(tmpdir, "autoclose.nc")

        sim, _mesh = _build_sim()
        sim.prepare()
        sim.init_output(outfile)
        sim.advance(0.01)
        sim.write_output()
        sim.destroy()  # should close output too

        # File should be readable (properly closed)
        ds = netCDF4.Dataset(outfile, "r")
        assert ds.dimensions["time"].size == 1
        ds.close()
