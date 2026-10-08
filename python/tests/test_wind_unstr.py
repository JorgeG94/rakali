"""Wind stress tests for the unstructured solver.

Two checks:

1. **Mechanism**: a uniform wind on a closed flat basin must inject
   x-momentum at the linear rate ``rho_water * V * (tau_x / rho_water)
   = V * tau_x`` per unit time, where ``V`` is the total water volume.
   We measure the rate-of-change of total xmom and compare against
   ``A * tau_x / rho_water`` (since hu has units m^2/s, the integral
   over area gives m^4/s, and d/dt of that is A * tau / rho).

2. **Wind setup**: a longer integration with bottom friction should
   tilt the free surface in the wind direction, with the slope
   approaching the analytical steady state
   ``deta/dx ~ tau_x / (rho * g * H)``. This is a coarse sanity check;
   reaching the actual steady state requires very long integrations,
   so we only check the *sign* and order-of-magnitude of the tilt.
"""
from __future__ import annotations

import numpy as np
import pytest

from rakali.unstr import UnstructuredMesh, UnstructuredSimulation

RHO_WATER = 1025.0
GRAVITY = 9.80665


def _build_basin(length=200.0, width=20.0, nx=40):
    dx = length / nx
    max_area = 0.5 * dx * dx
    return UnstructuredMesh.rectangle(
        x0=0.0, y0=0.0, x1=length, y1=width,
        max_area=max_area,
        bc_tags={"left": 1, "right": 1, "bottom": 1, "top": 1},  # all walls
    )


def _flat_lake(sim, depth):
    bed = np.zeros(sim.ncells, dtype=np.float64)
    stage = np.full(sim.ncells, depth, dtype=np.float64)
    sim.set_bathymetry(bed)
    sim.set_initial_stage(stage)


def test_wind_injects_momentum_at_predicted_rate():
    """Closed basin, no friction. Apply uniform wind for a short time
    and check that the integrated x-momentum grows at the analytical
    rate ``d/dt sum(hu*A) = sum(A) * tau_x / rho``.
    """
    mesh = _build_basin(length=100.0, width=10.0, nx=20)
    with UnstructuredSimulation(mesh) as sim:
        _flat_lake(sim, depth=5.0)
        tau_x = 0.5  # Pa
        sim.set_wind(tau_x, 0.0)
        sim.set_manning(0.0)
        sim.prepare()

        total_area = float(sim.area.sum())
        # advance just long enough for the linear regime to dominate
        # (depth is uniform, no pressure gradient yet -> momentum
        #  growth is essentially linear in t for small t)
        t_end = 1.0
        sim.advance(t_end)

        observed_xmom_integral = float((sim.xmom * sim.area).sum())
        predicted = total_area * tau_x / RHO_WATER * t_end

        # 20% tolerance: pressure gradient from the developing tilt
        # bites into the linear growth even at t = 1 s.
        rel = abs(observed_xmom_integral - predicted) / predicted
        assert rel < 0.20, (
            f"xmom integral {observed_xmom_integral:.3e} vs predicted "
            f"{predicted:.3e} (rel error {rel:.2%})")
        # And it must be of the right sign.
        assert observed_xmom_integral > 0.0


def test_wind_zero_no_motion():
    """tau = 0 -> still water stays still even with the wind kernel
    in the call chain (no spurious motion from the new code path)."""
    mesh = _build_basin(length=100.0, width=10.0, nx=20)
    with UnstructuredSimulation(mesh) as sim:
        _flat_lake(sim, depth=5.0)
        sim.set_wind(0.0, 0.0)
        sim.prepare()
        sim.advance(1.0)

        assert np.max(np.abs(sim.xmom)) < 1e-10
        assert np.max(np.abs(sim.ymom)) < 1e-10
        assert np.max(np.abs(sim.stage - 5.0)) < 1e-10


def test_wind_setup_sign():
    """Wind blowing in +x in a closed basin must raise the surface at
    the downwind (right) end and lower it at the upwind (left) end.
    Sign-only sanity check; reaching the analytical steady-state slope
    needs much longer integration.
    """
    length = 200.0
    depth = 5.0
    mesh = _build_basin(length=length, width=20.0, nx=40)
    with UnstructuredSimulation(mesh) as sim:
        _flat_lake(sim, depth=depth)
        sim.set_wind(1.0, 0.0)   # 1 Pa, strong wind
        sim.set_manning(0.025)
        sim.prepare()
        sim.advance(20.0)

        # Average stage in the leftmost 10% vs the rightmost 10%
        cx = sim.centroid_x
        left = sim.stage[cx < 0.1 * length].mean()
        right = sim.stage[cx > 0.9 * length].mean()

        assert right > left, (
            f"downwind stage {right:.6f} should exceed upwind {left:.6f}")
        # Order-of-magnitude: analytical steady tilt over the full
        # length is tau*L/(rho*g*H) = 1*200/(1025*9.80665*5) ~ 4e-3 m.
        # We are nowhere near steady state, so just check the tilt is
        # non-trivial (> 1e-6 m) and below 10x the steady value.
        tilt = right - left
        steady = 1.0 * length / (RHO_WATER * GRAVITY * depth)
        assert 1e-6 < tilt < 10.0 * steady, (
            f"tilt {tilt:.3e} not in (1e-6, {10*steady:.3e})")
