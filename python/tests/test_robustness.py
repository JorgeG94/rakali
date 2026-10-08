"""Guards that keep user mistakes from killing the interpreter.

Many invalid inputs would otherwise reach a Fortran ``ERROR STOP`` inside the
shared library, which terminates the whole Python process.  Each test here
would take pytest down with it if the guard regressed, so a green run is the
evidence that the guard works.
"""
import gc
import re
from pathlib import Path

import numpy as np
import pytest

from rakali import Rakali, UnstructuredMesh, UnstructuredSimulation

REPO = Path(__file__).resolve().parents[2]
BASE = dict(nx=8, ny=8, dx=1.0, dy=1.0, t_end=5.0)


# --- kwarg <-> namelist schema consistency ------------------------------


def _fortran_schema():
    src = (REPO / "src/core/rki_config.F90").read_text()
    body = src[src.index("subroutine build_rakali_schema"):]
    schema, grp = {}, None
    for line in body.splitlines():
        m = re.search(r'g%name\s*=\s*"([^"]+)"', line)
        if m:
            grp = m.group(1)
            schema.setdefault(grp, set())
            continue
        for m in re.finditer(r'nml_\w+\(\s*"([^"]+)"', line):
            if grp:
                schema[grp].add(m.group(1).lower())
    return schema


def test_every_kwarg_is_in_the_namelist_schema():
    # An unregistered key is a fatal strict-parse stop, not an exception.
    schema = _fortran_schema()
    assert schema, "could not parse the Fortran schema"
    bad = []
    for key, (group, fname) in Rakali._PARAM_MAP.items():
        g = group[:-4] if group.endswith("_nml") else group
        if fname.lower() not in schema.get(g, ()):
            bad.append(key)
    assert not bad, f"kwargs missing from the namelist schema: {bad}"


# --- kwarg validation ----------------------------------------------------


@pytest.mark.parametrize("bad", [
    dict(bc_west="foo"),
    dict(nx=0),
    dict(t_end=-1.0),
    dict(dx=float("nan")),
    dict(dx=None),
    dict(time_scheme="semi_implicit"),
    dict(time_scheme="rk4"),
    dict(vcoord_type="nope"),
    dict(remap_method="nope"),
    dict(vcoord_hc=5.0),
    dict(output_to_file=True),
    dict(forcing_file="f.nc"),
])
def test_bad_kwargs_raise_instead_of_stopping(bad):
    cfg = {**BASE, **bad}
    with pytest.raises(ValueError):
        Rakali(**cfg)


def test_numpy_scalars_are_accepted():
    cfg = dict(nx=np.int64(8), ny=np.int32(8), dx=np.float64(1.0),
               dy=np.float64(1.0), t_end=np.float64(5.0),
               use_multilayer=np.bool_(False))
    with Rakali(**cfg) as sim:
        sim.advance(0.5)


# --- file guards ---------------------------------------------------------


def test_from_restart_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        Rakali.from_restart(str(tmp_path / "nope.nc"), **BASE)


def test_save_restart_into_missing_directory_raises(tmp_path):
    with Rakali(**BASE) as sim:
        with pytest.raises(FileNotFoundError):
            sim.save_restart(str(tmp_path / "no_such_dir" / "x.nc"))


# --- structured run control ---------------------------------------------


def test_step_past_t_end_raises():
    with Rakali(**{**BASE, "t_end": 0.2}) as sim:
        sim.step(2000)  # the library clamps the last step to t_end
        assert sim.time == pytest.approx(0.2)
        with pytest.raises(RuntimeError, match="t_end"):
            sim.step(1)


def test_use_after_close_raises():
    sim = Rakali(**BASE)
    sim.close()
    with pytest.raises(Exception):
        sim.time
    sim.close()  # idempotent


def test_views_keep_the_simulation_alive():
    # Regression: a view returned past the owner's lifetime used to point at
    # freed Fortran memory (the wrapper's __del__ destroyed the handle).
    def make():
        sim = Rakali(**BASE)
        sim.advance(0.2)
        return sim.h_phys

    view = make()
    gc.collect()
    expected = view.copy()
    gc.collect()
    assert np.all(np.isfinite(view))
    assert np.array_equal(view, expected)


# --- unstructured input validation --------------------------------------

TRI = np.array([[1, 2, 3]], dtype=np.int32)
VX = np.array([0.0, 1.0, 0.0])
VY = np.array([0.0, 0.0, 1.0])


def test_mesh_rejects_zero_based_connectivity():
    with pytest.raises(ValueError, match="1-based"):
        UnstructuredMesh(VX, VY, TRI - 1)


def test_mesh_rejects_out_of_range_vertex():
    with pytest.raises(ValueError, match="1-based"):
        UnstructuredMesh(VX, VY, TRI + 1)


def test_mesh_rejects_oversized_bc_tag():
    with pytest.raises(ValueError, match="bc_tags"):
        UnstructuredMesh(VX, VY, TRI, np.array([[0, 0, 33]], dtype=np.int32))


def _usim():
    mesh = UnstructuredMesh.rectangle(x0=0, y0=0, x1=10, y1=10, max_area=5.0)
    return UnstructuredSimulation(mesh)


def test_unstr_rejects_bad_time_scheme_and_vcoord():
    with _usim() as sim:
        with pytest.raises(ValueError):
            sim.set_time_scheme("semi_implicit")
        with pytest.raises(ValueError):
            sim.set_time_scheme("rk4")
        with pytest.raises(ValueError):
            sim.set_vcoord(coord_type="nope")
        with pytest.raises(ValueError):
            sim.set_vcoord(remap_method="nope")
        sim.enable_multilayer(3)
        sim.set_vcoord(coord_type="sigma", remap_method="pqm")


def test_unstr_load_restart_missing_file_raises(tmp_path):
    with _usim() as sim:
        with pytest.raises(FileNotFoundError):
            sim.load_restart(str(tmp_path / "nope.nc"))


def test_missing_library_env_var_raises(monkeypatch):
    from rakali._core import _find_library
    monkeypatch.setenv("RAKALI_LIB", "/nonexistent/librakali.so")
    with pytest.raises(FileNotFoundError):
        _find_library()
