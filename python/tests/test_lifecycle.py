"""Lifecycle + create-path regression tests for the Python wrappers.

Covers the structured-API unification work:
  * P1 — the kwargs create path is fully in-memory (no temp file).
  * P2 — Rakali + UnstructuredSimulation share the _SimBase lifecycle
    (prepare / close / destroy / context-manager) and read the same.

The unstructured numerical path is exercised by the other test_*_unstr
files; here we focus on the shared front-end shape (structured side,
plus the cross-class invariants).
"""
import tempfile

import numpy as np
import pytest

from rakali import Rakali
from rakali.unstr import UnstructuredSimulation
from rakali._lifecycle import _SimBase


def test_both_wrappers_share_simbase_lifecycle():
    # The whole point of P2: the two front ends feel identical.
    assert issubclass(Rakali, _SimBase)
    assert issubclass(UnstructuredSimulation, _SimBase)
    for cls in (Rakali, UnstructuredSimulation):
        for name in ("prepare", "close", "destroy", "__enter__", "__exit__",
                     "step", "advance"):
            assert hasattr(cls, name), f"{cls.__name__} missing {name}"


def test_kwargs_create_uses_no_temp_file(monkeypatch):
    # P1: building a sim from kwargs must NOT round-trip through a temp
    # namelist file. Make any temp-file creation an immediate failure.
    def _boom(*args, **kwargs):
        raise AssertionError("Rakali(**kwargs) must not create a temp file")

    monkeypatch.setattr(tempfile, "NamedTemporaryFile", _boom)

    with Rakali(nx=24, ny=8, dx=1.0, dy=1.0, t_end=5.0,
                init_condition="dam_break", log_level="warning") as sim:
        sim.step(3)
        assert sim.step_count == 3
        assert sim.h_phys.max() > 1.0  # dam break has water

    # Context-manager exit tore the handle down.
    assert sim._alive is False


def test_prepare_is_canonical_and_close_is_idempotent():
    # P2: prepare() is the canonical "ready to run" verb on both backends
    # (alias of the legacy finalize()); close()/destroy() are idempotent.
    sim = Rakali(nx=16, ny=8, dx=1.0, dy=1.0, t_end=5.0,
                 defer_init=True, log_level="warning")
    sim.set_bathymetry(np.zeros((16, 8), order="F"))
    sim.set_initial_h(np.ones((16, 8), order="F"))
    sim.prepare()
    sim.step(2)
    assert sim.step_count == 2

    sim.close()
    sim.close()    # idempotent — no double-free
    sim.destroy()  # alias, also idempotent


def test_context_manager_frees_on_exception():
    # An exception inside the `with` block must still release the handle.
    sim_ref = {}
    with pytest.raises(ValueError):
        with Rakali(nx=12, ny=6, dx=1.0, dy=1.0, t_end=5.0,
                    init_condition="dam_break", log_level="warning") as sim:
            sim_ref["s"] = sim
            raise ValueError("boom")
    assert sim_ref["s"]._alive is False
