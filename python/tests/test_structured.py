"""Structured-grid (``Rakali``) wrapper tests.

Exercises the documented workflow end to end: kwargs / namelist creation,
two-phase init with numpy bathymetry, state injection, callbacks, restart
round-trip, multilayer tracers, and xarray export.  Physics checks are
analytical (lake at rest, closed-basin mass conservation, uniform tracer).
"""
from pathlib import Path

import numpy as np
import pytest

from rakali import Rakali

REPO = Path(__file__).resolve().parents[2]
NX, NY = 32, 16


def _basin(**kw):
    """Closed (all-wall) flat basin, two-phase init left to the caller."""
    cfg = dict(nx=NX, ny=NY, dx=10.0, dy=10.0, t_end=60.0,
               bc_west="wall", bc_east="wall", bc_south="wall", bc_north="wall")
    cfg.update(kw)
    return Rakali(defer_init=True, **cfg)


def _finalized_basin(depth=2.0, bump=None, **kw):
    sim = _basin(**kw)
    b = np.full((NX, NY), -depth)
    sim.set_bathymetry(b)
    h = np.full((NX, NY), depth)
    if bump is not None:
        h = h + bump
    sim.set_initial_h(h)
    sim.finalize()
    return sim


def test_kwargs_create_and_advance():
    with Rakali(nx=NX, ny=NY, dx=10.0, dy=10.0, t_end=10.0,
                bc_west="wall", bc_east="wall") as sim:
        sim.advance(2.0)
        assert sim.time == pytest.approx(2.0)
        assert sim.h_phys.shape == (NX, NY)
        assert np.all(np.isfinite(sim.h_phys))


def test_lake_at_rest_stays_at_rest():
    with _finalized_basin() as sim:
        h0 = sim.h_phys.copy()
        sim.advance(30.0)
        assert np.max(np.abs(sim.h_phys - h0)) < 1e-10
        assert np.max(np.abs(sim.hu_phys)) < 1e-10


def test_closed_basin_conserves_mass():
    bump = np.zeros((NX, NY))
    bump[NX // 2 - 3:NX // 2 + 3, NY // 2 - 3:NY // 2 + 3] = 0.3
    with _finalized_basin(bump=bump) as sim:
        m0 = sim.h_phys.sum()
        sim.advance(30.0)
        assert sim.h_phys.sum() == pytest.approx(m0, rel=1e-12)
        # ...and the bump actually evolved, so the check is not vacuous.
        assert sim.h_phys.max() < 2.3 - 1e-3


def test_set_h_round_trip():
    with _finalized_basin() as sim:
        new_h = np.full((NX, NY), 2.5)
        sim.set_h(new_h)
        sim.sync_host()
        assert np.allclose(sim.h_phys, 2.5)


def test_callback_fires_at_interval():
    times = []
    with _finalized_basin() as sim:
        sim.advance(20.0, callback=lambda s: times.append(s.time),
                    callback_interval=5.0)
    assert len(times) >= 3
    assert all(b > a for a, b in zip(times, times[1:]))


def test_wrong_shape_rejected():
    with _finalized_basin() as sim:
        with pytest.raises(ValueError):
            sim.set_h(np.ones((NX + 1, NY)))


def test_namelist_file_creation(tmp_path):
    nml = tmp_path / "basin.nml"
    nml.write_text(Rakali._kwargs_to_namelist(
        dict(nx=NX, ny=NY, dx=10.0, dy=10.0, t_end=10.0,
             bc_west="wall", bc_east="wall")))
    with Rakali(str(nml)) as sim:
        sim.advance(1.0)
        assert np.all(np.isfinite(sim.h_phys))


def test_missing_namelist_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        Rakali(str(tmp_path / "nope.nml"))


def test_missing_bathymetry_raises_not_kills_interpreter(tmp_path):
    # Regression: a missing input file used to reach a Fortran ERROR STOP
    # that terminated the whole Python process.
    with pytest.raises(FileNotFoundError, match="bathymetry_file"):
        Rakali(nx=NX, ny=NY, dx=10.0, dy=10.0, t_end=10.0,
               bathymetry_file=str(tmp_path / "nope.nc"))
    nml = tmp_path / "bad.nml"
    nml.write_text('&output_nml\n  bathymetry_file = "nope.nc"\n/\n')
    with pytest.raises(FileNotFoundError, match="bathymetry_file"):
        Rakali(str(nml))


def test_multilayer_uniform_salinity_stays_uniform():
    with Rakali(nx=NX, ny=NY, dx=10.0, dy=10.0, t_end=20.0,
                bc_west="wall", bc_east="wall", bc_south="wall",
                bc_north="wall", use_multilayer=True, nz_layers=4,
                initial_salinity=30.0) as sim:
        sim.advance(10.0)
        S = sim.salinity
        assert S is not None
        live = S[np.isfinite(S) & (sim.h_layer > 1e-6)]
        assert np.allclose(live, 30.0, atol=1e-8)


def test_restart_round_trip(tmp_path):
    cfg = dict(nx=NX, ny=NY, dx=10.0, dy=10.0, t_end=60.0,
               bc_west="wall", bc_east="wall", bc_south="wall",
               bc_north="wall")
    ckpt = str(tmp_path / "ckpt.nc")
    with Rakali(**cfg) as a:
        a.advance(5.0)
        t_saved = a.time
        h_saved = a.h_phys.copy()
        a.save_restart(ckpt)
    with Rakali.from_restart(ckpt, **cfg) as b:
        assert b.time == pytest.approx(t_saved)
        assert np.allclose(b.h_phys, h_saved, atol=1e-12)


def test_to_xarray():
    xr = pytest.importorskip("xarray")
    with _finalized_basin() as sim:
        sim.advance(2.0)
        ds = sim.to_xarray()
        assert isinstance(ds, xr.Dataset)
        assert {"h", "u", "v"} <= set(ds.data_vars)
        assert ds["h"].shape[-2:] == (NX, NY) or ds["h"].shape[-2:] == (NY, NX)
