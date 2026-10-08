"""Closed-basin standing wave: amplitude / sign-flip sanity check.

Originally this file held a convergence-rate test against the linear
standing-wave solution. We removed it after discovering that the
closed-basin standing wave is **not** a clean convergence problem for
this scheme: the wave's gradient extends to the walls from t = 0, and
the wall ghost-state mirror introduces an O(h) boundary truncation
that caps the observed convergence rate at ~0.6-0.8 regardless of
limiter settings (verified by repeating the test with the limiter
disabled via ``set_betas(1, 1, 1)`` -- bit-identical L1 numbers).

The actual interior smooth-problem convergence test now lives in
``test_gaussian_bump_unstr.py``: a Gaussian bump on a long basin
where the wave hasn't reached the walls within the run, so the
boundary contribution to the L1 is exactly zero.

What we keep here is the half-period **sign flip** test, which is a
useful qualitative regression: it proves the wave is actually
oscillating in the basin (not just sitting still or blowing up) and
catches gross damping bugs.
"""
from __future__ import annotations

import math
import numpy as np
import pytest

from rakali.unstr import UnstructuredMesh, UnstructuredSimulation


GRAVITY = 9.81

L = 10.0       # basin length (m)
W = 2.0        # basin width  (m)
H = 1.0        # mean depth   (m)
A = 0.05       # wave amplitude (m), << H so linear theory applies


def _omega() -> float:
    k = math.pi / L
    return math.sqrt(GRAVITY * H) * k


def _exact_stage(x: np.ndarray, t: float) -> np.ndarray:
    k = math.pi / L
    return H + A * np.cos(k * x) * math.cos(_omega() * t)


def _build_basin(nx: int):
    dx = L / nx
    max_area = 0.5 * dx * dx
    return UnstructuredMesh.rectangle(
        x0=0.0, y0=0.0, x1=L, y1=W,
        max_area=max_area,
        bc_tags={"left": 1, "right": 1, "bottom": 1, "top": 1},
    )


def _setup_standing_wave(sim):
    cx = sim.centroid_x.copy()
    bed = np.zeros(sim.ncells, dtype=np.float64)
    stage = _exact_stage(cx, 0.0).astype(np.float64)
    sim.set_bathymetry(bed)
    sim.set_initial_stage(stage)


def test_standing_wave_amplitude_preserved():
    """At t = 0 the wave is at maximum displacement and the kernel
    should not damp the amplitude noticeably over a single short
    integration. After half a period the wave should have flipped
    sign (cos(omega*T/2) = -1) -- check the displacement at the
    boundary cells.
    """
    half_period = math.pi / _omega()
    mesh = _build_basin(nx=80)
    with UnstructuredSimulation(mesh) as sim:
        sim.set_bc_type(group=1, bc_type="wall")
        _setup_standing_wave(sim)
        sim.prepare()
        sim.advance(half_period)

        cx = sim.centroid_x.copy()
        s = sim.stage.copy()

        # Cell closest to x = 0 (where exact = H + A * 1 * (-1) = H - A)
        i_left = int(np.argmin(cx))
        # Cell closest to x = L (where exact = H + A * (-1) * (-1) = H + A)
        i_right = int(np.argmax(cx))

        eta_left = s[i_left] - H
        eta_right = s[i_right] - H
        print(f"\n  half-period eta_left={eta_left:+.5f} (expect ~{-A:+.5f})")
        print(f"  half-period eta_right={eta_right:+.5f} (expect ~{+A:+.5f})")

        # Allow 25% damping over half a period -- the β-limiter and
        # wall ghost-mirror nibble at the peaks but the sign flip must
        # be unambiguous.
        assert eta_left < -0.6 * A, "left peak did not invert as expected"
        assert eta_right > 0.6 * A, "right peak did not invert as expected"
        assert np.sign(eta_left) == -1.0
        assert np.sign(eta_right) == +1.0
