"""Pure-Python validation of the unstructured solver against the
exact Ritter / Toro Riemann solution for a wet/wet dam break.

Mirrors `tests/analytical/test_dam_break_unstr.F90`, but driven entirely
from Python via the C ABI added in `rki_api_unstr.F90`. This is the
regression net for every future change to the Python path: if Python
and Fortran agree on the same numerical answer to within ~5% L1 of
the analytical solution, the whole stack is consistent.

Tolerances mirror the Fortran test: L1 / jump < 5% on a ~100 cell
channel, no negative depth, no NaN, undisturbed endpoints far from
the wave system.
"""
from __future__ import annotations

import math
import numpy as np
import pytest

from rakali.unstr import UnstructuredMesh, UnstructuredSimulation


GRAVITY = 9.81


# ---------------------------------------------------------------------------
# Pure-Python Toro Riemann solver -- ports tests/analytical/ritter_dambreak.F90
# ---------------------------------------------------------------------------


def _exact_riemann_dambreak(h_l: float, h_r: float) -> tuple[float, float]:
    """Solve the wet/wet Riemann problem with u_L = u_R = 0.

    Returns (h_star, u_star). Uses Newton iteration on
    f(h) = f_L(h) + f_R(h) = 0, following Toro 3rd ed. Ch. 5.
    """
    a_l = math.sqrt(GRAVITY * h_l)
    a_r = math.sqrt(GRAVITY * h_r)
    # Two-rarefaction approximation as initial guess
    h_guess = ((a_l + a_r) / 2.0) ** 2 / GRAVITY

    for _ in range(50):
        # Left wave
        if h_guess <= h_l:
            f_l = 2.0 * (math.sqrt(GRAVITY * h_guess) - a_l)
            df_l = math.sqrt(GRAVITY / h_guess)
        else:
            q_k = math.sqrt(0.5 * GRAVITY * (h_guess + h_l) / (h_guess * h_l))
            f_l = (h_guess - h_l) * q_k
            df_l = q_k - 0.25 * GRAVITY * (h_guess - h_l) / (q_k * h_guess * h_guess)

        # Right wave
        if h_guess <= h_r:
            f_r = 2.0 * (math.sqrt(GRAVITY * h_guess) - a_r)
            df_r = math.sqrt(GRAVITY / h_guess)
        else:
            q_k = math.sqrt(0.5 * GRAVITY * (h_guess + h_r) / (h_guess * h_r))
            f_r = (h_guess - h_r) * q_k
            df_r = q_k - 0.25 * GRAVITY * (h_guess - h_r) / (q_k * h_guess * h_guess)

        delta = (f_l + f_r) / (df_l + df_r)
        h_guess = max(h_guess - delta, 1.0e-15)
        if abs(delta) < 1.0e-12 * h_guess:
            break

    h_star = h_guess
    u_star = 0.5 * (f_r - f_l)
    return h_star, u_star


def _exact_riemann_sample(
    x: np.ndarray, t: float, h_l: float, h_r: float, x_dam: float
) -> np.ndarray:
    """Vectorised sample of the exact h(x, t) profile.

    Standard wet/wet dam break with h_L > h_R: left rarefaction,
    right shock, star region between.
    """
    h_star, u_star = _exact_riemann_dambreak(h_l, h_r)
    a_l = math.sqrt(GRAVITY * h_l)
    a_r = math.sqrt(GRAVITY * h_r)
    a_star = math.sqrt(GRAVITY * h_star)

    h_out = np.empty_like(x)

    if t < 1.0e-15:
        h_out[:] = np.where(x < x_dam, h_l, h_r)
        return h_out

    xi = (x - x_dam) / t

    # Left rarefaction (h_star < h_l for the standard case)
    s_head_l = -a_l
    s_tail_l = u_star - a_star
    # Right shock (h_star > h_r)
    s_shock_r = a_r * math.sqrt(h_star * (h_star + h_r) / (2.0 * h_r * h_r))

    h_fan = ((2.0 * a_l - xi) / 3.0) ** 2 / GRAVITY

    h_out[:] = h_l
    in_fan = (xi > s_head_l) & (xi <= s_tail_l)
    in_star = (xi > s_tail_l) & (xi < s_shock_r)
    in_right = xi >= s_shock_r

    h_out = np.where(in_fan, h_fan, h_out)
    h_out = np.where(in_star, h_star, h_out)
    h_out = np.where(in_right, h_r, h_out)
    return h_out


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_channel_sim(nx: int, length: float = 100.0, width: float = 2.0):
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
    return sim


def _setup_dam_break(sim, h_l: float, h_r: float, x_dam: float):
    cx = sim.centroid_x.copy()
    bed = np.zeros(sim.ncells, dtype=np.float64)
    stage = np.where(cx < x_dam, h_l, h_r).astype(np.float64)
    sim.set_bathymetry(bed)
    sim.set_initial_stage(stage)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_ritter_l1_under_5_percent():
    """Area-weighted L1 error against the exact solution stays under 5%
    of the dam jump on a ~100-cell channel."""
    h_l, h_r = 2.0, 1.0
    x_dam = 50.0
    t_end = 0.1
    length = 100.0

    with _build_channel_sim(nx=100, length=length) as sim:
        _setup_dam_break(sim, h_l, h_r, x_dam)
        sim.prepare()
        sim.advance(t_end)

        cx = sim.centroid_x.copy()
        h_num = sim.height.copy()
        # Use mesh areas via the relation area = (radius * perimeter / 2);
        # we don't expose radius from Python, so weight by uniform 1.0
        # -- the channel is fairly uniform and this matches the Fortran
        # test pattern (see test_dam_break_unstr.F90).
        h_exact = _exact_riemann_sample(cx, t_end, h_l, h_r, x_dam)

        # Sanity
        assert sim.time == pytest.approx(t_end, abs=1e-9)
        assert h_num.min() > 0.0, "produced a non-positive depth"
        assert not np.any(np.isnan(h_num)), "NaN in height field"

        # L1 error normalised by the dam jump
        l1 = float(np.mean(np.abs(h_num - h_exact)))
        l1_norm = l1 / (h_l - h_r)
        print(f"\n  ncells={sim.ncells}, t={sim.time:.3f}, "
              f"steps={sim.step_count}, L1/jump={l1_norm:.4f}")
        assert l1_norm < 0.05, f"L1/jump = {l1_norm:.4f} >= 0.05"


def test_undisturbed_endpoints():
    """Cells far from the wave system stay at the original depth."""
    h_l, h_r = 2.0, 1.0
    x_dam = 50.0
    t_end = 0.1
    length = 100.0

    a_l = math.sqrt(GRAVITY * h_l)
    a_r = math.sqrt(GRAVITY * h_r)
    # Left rarefaction head moves at -a_l, right shock moves at +s_shock_r
    h_star, u_star = _exact_riemann_dambreak(h_l, h_r)
    s_shock_r = a_r * math.sqrt(h_star * (h_star + h_r) / (2.0 * h_r * h_r))
    # Wave reach with safety margin
    left_reach = x_dam - a_l * t_end - 2.0
    right_reach = x_dam + s_shock_r * t_end + 2.0

    with _build_channel_sim(nx=100, length=length) as sim:
        _setup_dam_break(sim, h_l, h_r, x_dam)
        sim.prepare()
        sim.advance(t_end)

        cx = sim.centroid_x.copy()
        h_num = sim.height.copy()

        far_left = cx < left_reach
        far_right = cx > right_reach
        assert np.any(far_left) and np.any(far_right)
        assert np.allclose(h_num[far_left], h_l, atol=1.0e-8)
        assert np.allclose(h_num[far_right], h_r, atol=1.0e-8)


def test_no_y_drift():
    """In a channel symmetric about y, |ymom| should stay much smaller
    than |xmom|. Tolerance is loose because the triangulation is not
    y-symmetric and ymom can pick up a small bias near walls."""
    with _build_channel_sim(nx=100) as sim:
        _setup_dam_break(sim, h_l=2.0, h_r=1.0, x_dam=50.0)
        sim.prepare()
        sim.advance(0.1)

        max_xmom = float(np.abs(sim.xmom).max())
        max_ymom = float(np.abs(sim.ymom).max())
        assert max_xmom > 0.5, "no flow developed at all"
        assert max_ymom < 0.5 * max_xmom, (
            f"y-momentum {max_ymom:.3f} too large vs x-momentum {max_xmom:.3f}"
        )


def test_mesh_refinement_reduces_error():
    """L1 error decreases monotonically as the mesh refines."""
    h_l, h_r = 2.0, 1.0
    x_dam = 50.0
    t_end = 0.1
    length = 100.0

    errors = []
    for nx in (50, 100, 200):
        with _build_channel_sim(nx=nx, length=length) as sim:
            _setup_dam_break(sim, h_l, h_r, x_dam)
            sim.prepare()
            sim.advance(t_end)

            cx = sim.centroid_x.copy()
            h_num = sim.height.copy()
            h_exact = _exact_riemann_sample(cx, t_end, h_l, h_r, x_dam)
            l1 = float(np.mean(np.abs(h_num - h_exact)))
            errors.append(l1)
            print(f"\n  nx={nx} ncells={sim.ncells} L1={l1:.5f}")

    assert errors[1] < errors[0], "L1 did not decrease 50 -> 100"
    assert errors[2] < errors[1], "L1 did not decrease 100 -> 200"
