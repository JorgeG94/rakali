"""Tests for the rakali.mesh coastal builder.

These exercise the pure-Python mesh construction primitives and the
top-level :func:`build_coastal_mesh` on a synthetic raster bathymetry.
No solver state — the library is never created, so these tests run without the
shared library.
"""
from __future__ import annotations

import numpy as np
import pytest

from rakali.mesh import (
    build_coastal_mesh,
    build_pslg,
    extract_coastlines,
    interpolate_bed,
    signed_area,
    simplify_polygon,
    triangulate,
)
from rakali.mesh.coastline import filter_islands_by_centroid
from rakali.mesh.triangulate import make_side_marker


# ---------------------------------------------------------------------------
# Synthetic bathymetries
# ---------------------------------------------------------------------------


def _disk_bathy(nx=120, ny=120, dx=50.0, depth=-10.0):
    """Circular basin: wet inside a radius of ``0.4 * nx * dx``, dry outside.

    Returns ``(b, dx, dy)`` with cell centres on the standard ``i*dx`` grid.
    """
    x = np.arange(nx) * dx
    y = np.arange(ny) * dx
    cx, cy = (nx - 1) * dx / 2.0, (ny - 1) * dx / 2.0
    r = 0.4 * nx * dx
    xx, yy = np.meshgrid(x, y, indexing="ij")
    dist = np.hypot(xx - cx, yy - cy)
    b = np.where(dist <= r, depth, 20.0)
    return b, dx, dx


def _island_bathy(nx=200, ny=200, dx=50.0):
    """Sea basin with one circular island near the centre.  Used to
    check island detection and hole insertion."""
    x = np.arange(nx) * dx
    y = np.arange(ny) * dx
    cx, cy = (nx - 1) * dx / 2.0, (ny - 1) * dx / 2.0
    xx, yy = np.meshgrid(x, y, indexing="ij")
    dist_outer = np.hypot(xx - cx, yy - cy)
    b = np.where(dist_outer <= 0.4 * nx * dx, -15.0, 50.0)
    # Island: dry disk inside the basin
    dist_isl = np.hypot(xx - cx - 1000.0, yy - cy)
    b[dist_isl <= 500.0] = 10.0
    return b, dx, dx


# ---------------------------------------------------------------------------
# Primitive tests
# ---------------------------------------------------------------------------


def test_signed_area_sign_convention():
    """CCW square has positive signed area; flipped version negative."""
    sq = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]])
    assert signed_area(sq) == pytest.approx(1.0)
    assert signed_area(sq[::-1]) == pytest.approx(-1.0)


def test_simplify_polygon_preserves_endpoints():
    """Simplification must keep the first and last point."""
    pts = np.array([[float(i), 0.0] for i in range(11)])
    out = simplify_polygon(pts, min_dist=2.5)
    assert np.allclose(out[0], pts[0])
    assert np.allclose(out[-1], pts[-1])
    # No two consecutive kept points closer than min_dist (except
    # possibly the last pair, because endpoints are always kept).
    d = np.hypot(np.diff(out[:-1, 0]), np.diff(out[:-1, 1]))
    assert np.all(d >= 2.5 - 1e-12)


def test_extract_coastlines_disk():
    """A circular wet basin extracts a single closed outer contour."""
    b, dx, _ = _disk_bathy()
    x = np.arange(b.shape[0]) * dx
    y = np.arange(b.shape[1]) * dx
    outer, islands = extract_coastlines(b, x, y, sea_level=0.0, simplify_tol=dx)
    assert outer.shape[1] == 2
    assert len(outer) >= 16  # discretised circle, not collapsed to a point
    assert islands == []
    # CCW by construction
    assert signed_area(outer) > 0.0


def test_make_side_marker():
    """make_side_marker tags east edges as open_marker and others as wall."""
    bbox = (0.0, 100.0, 0.0, 50.0)
    m = make_side_marker(bbox, ["east"], tol=1.0, open_marker=2, wall_marker=1)
    assert m(99.9, 25.0) == 2  # on east side
    assert m(0.5, 25.0) == 1   # west interior
    assert m(50.0, 0.5) == 1   # south
    assert m(50.0, 49.9) == 1  # north


# ---------------------------------------------------------------------------
# End-to-end
# ---------------------------------------------------------------------------


def test_build_coastal_mesh_disk():
    """Full pipeline: synthetic disk bathy → mesh + bed at centroids."""
    b, dx, _ = _disk_bathy()
    result = build_coastal_mesh(
        b, resolution=200.0, dx=dx, sea_level=0.0, verbose=False,
    )
    assert result.mesh.ntri > 10
    assert result.bed.shape == (result.mesh.ntri,)
    # Centroids lie inside the wet region, so bed should be negative
    # (the disk bathy uses -10 m inside).
    assert np.all(result.bed < 0.0)
    # All boundary edges on the outer polygon should be walls (no
    # open_boundary was requested).
    assert int(np.sum(result.mesh.bc_tags == result.wall_marker)) > 0
    assert int(np.sum(result.mesh.bc_tags == result.open_marker)) == 0


def test_build_coastal_mesh_open_east():
    """``open_boundary='east'`` tags at least some edges with the open
    marker on a domain that actually reaches the east edge."""
    # Half-plane bathy: wet in x >= 2000 m, dry west of that.  The wet
    # region touches the east raster boundary so the coastline has a
    # segment along x = x_max.
    nx, ny, dx = 120, 120, 50.0
    x = np.arange(nx) * dx
    y = np.arange(ny) * dx
    xx, _ = np.meshgrid(x, y, indexing="ij")
    b = np.where(xx >= 2000.0, -10.0, 20.0)
    result = build_coastal_mesh(
        b, resolution=200.0, dx=dx, sea_level=0.0,
        open_boundary="east", verbose=False,
    )
    n_open = int(np.sum(result.mesh.bc_tags == result.open_marker))
    assert n_open > 0, "expected at least one open boundary edge"


def test_build_coastal_mesh_island():
    """Island bathy produces an island hole in the PSLG."""
    b, dx, _ = _island_bathy()
    result = build_coastal_mesh(
        b, resolution=250.0, dx=dx, sea_level=0.0,
        min_island_area=1.0e5, verbose=False,
    )
    assert len(result.islands) == 1
    # No triangle should have its centroid inside the island.
    cx, cy = result.mesh.centroids()
    xx = (result.bathymetry.shape[0] - 1) * dx / 2.0 + 1000.0
    yy = (result.bathymetry.shape[1] - 1) * dx / 2.0
    d = np.hypot(cx - xx, cy - yy)
    assert np.all(d > 450.0)  # island radius is 500 m, leave margin


def test_interpolate_bed_matches_raster():
    """``interpolate_bed`` with nearest method returns the raster value
    at each centroid."""
    b, dx, _ = _disk_bathy()
    x = np.arange(b.shape[0]) * dx
    y = np.arange(b.shape[1]) * dx
    # Sample at cell centres: should return exactly b[i, j].
    cx, cy = np.meshgrid(x, y, indexing="ij")
    vals = interpolate_bed(b, x, y, cx.ravel(), cy.ravel(), method="nearest")
    assert np.allclose(vals, b.ravel())
