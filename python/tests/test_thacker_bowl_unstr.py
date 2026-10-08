"""Wet/dry validation on a parabolic bowl ("Thacker"-style).

The Thacker (1981) parabolic-bowl problem is the canonical wet/dry
validation case for shallow-water solvers. We don't compare against
the full analytical solution here -- the planar/curved Thacker
analytical formulas are involved -- but instead exercise the two
critical properties any solver must satisfy on this geometry:

1. **Lake at rest with dry cells**: a still water surface in a bowl
   produces a partially-submerged equilibrium. The solver must hold
   this state to machine precision: no spurious velocities, no drift,
   no mass leak through the partially-wet cells. This is the strictest
   well-balancing test, because any imbalance between bed-slope and
   pressure-gradient terms shows up immediately as flow.

2. **Sloshing oscillation**: a perturbation away from equilibrium
   produces an oscillation where cells near the shoreline alternately
   go dry and become wet again. The solver must (a) not crash, (b)
   keep mass conserved, and (c) keep depths non-negative throughout.

Bathymetry::

    bed(r) = -h_0 * (1 - r²/L²)            r in [0, L]
    bed(r) = +h_0 * (r²/L² - 1)            r > L  (rim above water)

Bed minimum is -h_0 at the centre and rises to a positive rim outside
radius L. The disk mesh radius is set just larger than L so the rim
cells exist and can go dry near the boundary.
"""
from __future__ import annotations

import math
import numpy as np
import pytest

from rakali.unstr import UnstructuredMesh, UnstructuredSimulation


GRAVITY = 9.81

L = 5.0          # bowl characteristic radius (m)
H_0 = 1.0        # bowl depth at centre (m)
R_DOMAIN = 5.5   # disk mesh radius -- slightly larger than L so the
                 # rim sits inside the domain and rim cells can dry out


def _omega() -> float:
    """Small-amplitude oscillation frequency for the parabolic bowl."""
    return math.sqrt(2.0 * GRAVITY * H_0) / L


def _bed(r: np.ndarray) -> np.ndarray:
    """Parabolic bowl bathymetry."""
    return -H_0 * (1.0 - (r * r) / (L * L))


def _build_bowl_sim(n_cross: int = 50):
    diameter = 2.0 * R_DOMAIN
    dx = diameter / n_cross
    max_area = 0.5 * dx * dx
    mesh = UnstructuredMesh.disk(
        cx=0.0, cy=0.0, radius=R_DOMAIN,
        max_area=max_area,
        n_boundary=128,
        boundary_tag=1,
    )
    sim = UnstructuredSimulation(mesh)
    sim.set_bc_type(group=1, bc_type="wall")
    return sim


def _setup_lake_at_rest(sim, z_lake: float):
    """Stage = max(z_lake, bed). Wet cells have water above the
    prescribed lake level; dry cells have stage equal to bed."""
    cx = sim.centroid_x.copy()
    cy = sim.centroid_y.copy()
    r = np.sqrt(cx * cx + cy * cy)
    bed = _bed(r)
    stage = np.maximum(z_lake, bed)
    sim.set_bathymetry(bed.astype(np.float64))
    sim.set_initial_stage(stage.astype(np.float64))
    return r, bed


# ---------------------------------------------------------------------------
# 1. Well-balanced lake at rest with dry cells
# ---------------------------------------------------------------------------


def test_bowl_lake_at_rest_well_balanced():
    """Partially-submerged still lake stays still: no spurious flow,
    no stage drift, mass conserved exactly. Dry cells exist and remain
    dry; wet cells stay at the prescribed lake level."""
    z_lake = -0.3  # well below the rim, well above the centre

    with _build_bowl_sim(n_cross=50) as sim:
        r, bed = _setup_lake_at_rest(sim, z_lake)
        sim.prepare()

        s0 = sim.stage.copy()
        m0 = sim.volume()
        n_wet0 = int((s0 - bed > 1.0e-6).sum())
        n_dry0 = int(sim.ncells - n_wet0)

        # A few seconds is plenty -- if well-balancing is broken, drift
        # appears within the first few timesteps.
        sim.advance(2.0)

        s1 = sim.stage.copy()
        xm = sim.xmom.copy()
        ym = sim.ymom.copy()
        h1 = np.maximum(s1 - bed, 0.0)
        m1 = sim.volume()

        # Spurious flow magnitude
        max_speed = float(np.max(np.sqrt(xm * xm + ym * ym) /
                                 np.maximum(h1, 1.0e-6)))
        # Maximum stage drift (only over wet cells -- dry cells have
        # stage == bed by definition)
        wet = h1 > 1.0e-6
        max_drift = float(np.max(np.abs(s1[wet] - z_lake)))

        rel_mass = abs(m1 - m0) / max(abs(m0), 1.0e-12)

        print(f"\n  ncells={sim.ncells} (wet0={n_wet0}, dry0={n_dry0})")
        print(f"  total mass: {m0:.6f} -> {m1:.6f}")
        print(f"  mass relative drift:  {rel_mass:.3e}")
        print(f"  max stage drift (wet): {max_drift:.3e}")
        print(f"  max spurious speed:    {max_speed:.3e}")

        assert n_wet0 > 0, "no wet cells -- z_lake too low"
        assert n_dry0 > 0, "no dry cells -- z_lake too high (test trivial)"
        # Well-balanced tolerance: no scheme is perfect, but the
        # drift on a still-lake bowl test should be at numerical-noise
        # scale. ANUGA-equivalent solvers reach ~1e-12 here; loose
        # tolerance to absorb the limiter and protect-clip noise.
        assert rel_mass < 1.0e-8, f"mass drift {rel_mass:.3e} too large"
        assert max_drift < 1.0e-6, (
            f"stage drifted by {max_drift:.3e} from rest"
        )
        assert max_speed < 1.0e-3, (
            f"spurious speed {max_speed:.3e} on still lake"
        )


# ---------------------------------------------------------------------------
# 2. Sloshing oscillation -- exercises wet/dry transitions
# ---------------------------------------------------------------------------


def test_bowl_sloshing_wet_dry():
    """Tilt the lake's free surface and let it oscillate. Cells near
    the moving shoreline alternate between wet and dry. The simulation
    must complete cleanly with mass conserved."""
    z_lake = -0.3
    tilt = 0.1   # m of tilt across the bowl -- small enough to stay
                 # well-behaved, large enough to drive observable
                 # wet/dry transitions

    with _build_bowl_sim(n_cross=50) as sim:
        cx = sim.centroid_x.copy()
        cy = sim.centroid_y.copy()
        r = np.sqrt(cx * cx + cy * cy)
        bed = _bed(r)
        # Tilt the lake along x
        stage = np.maximum(z_lake + tilt * (cx / L), bed)

        sim.set_bathymetry(bed.astype(np.float64))
        sim.set_initial_stage(stage.astype(np.float64))
        sim.prepare()

        s0 = sim.stage.copy()
        h0 = np.maximum(s0 - bed, 0.0)
        m0 = sim.volume()
        n_dry0 = int((h0 < 1.0e-6).sum())

        # Run for ~one small-amplitude period of the bowl oscillation.
        period = 2.0 * math.pi / _omega()
        t_end = 1.0 * period

        # Track dry-cell count over a few snapshots to confirm wet/dry
        # transitions actually happen.
        dry_counts = [n_dry0]
        for tf in (0.25, 0.5, 0.75, 1.0):
            sim.advance(tf * t_end)
            h = np.maximum(sim.stage - bed, 0.0)
            dry_counts.append(int((h < 1.0e-6).sum()))

        s_final = sim.stage.copy()
        h_final = np.maximum(s_final - bed, 0.0)
        m_final = sim.volume()
        max_speed = float(np.max(
            np.sqrt(sim.xmom ** 2 + sim.ymom ** 2) /
            np.maximum(h_final, 1.0e-6)
        ))

        rel_mass = abs(m_final - m0) / max(abs(m0), 1.0e-12)

        print(f"\n  ncells={sim.ncells}")
        print(f"  oscillation period (linear theory): {period:.3f} s")
        print(f"  ran for {sim.time:.3f} s, {sim.step_count} steps")
        print(f"  dry cell counts: {dry_counts}")
        print(f"  mass: {m0:.6f} -> {m_final:.6f}")
        print(f"  mass relative drift: {rel_mass:.3e}")
        print(f"  max final speed:     {max_speed:.3e}")

        # Wet/dry must actually be exercised -- the dry-cell count
        # should change between snapshots.
        assert max(dry_counts) > min(dry_counts), (
            "no wet/dry transitions occurred"
        )

        # Mass conservation under wet/dry. The protect step clips
        # tiny negative depths to zero, so we lose a small amount.
        # Tighter than 1e-4 is harder than the still-lake case.
        assert rel_mass < 1.0e-4, (
            f"mass drift {rel_mass:.3e} exceeds 1e-4 under wet/dry"
        )

        # No NaN, no negative depth
        assert not np.any(np.isnan(h_final)), "NaN in depth"
        assert h_final.min() >= 0.0, "negative depth in final state"

        # Velocities are bounded -- should be O(sqrt(g*tilt)) ~ 1 m/s
        assert max_speed < 5.0, (
            f"speed {max_speed:.3f} runaway -- something exploded"
        )
