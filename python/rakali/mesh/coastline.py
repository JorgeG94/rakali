"""Coastline extraction from a raster bathymetry.

The coastline of a structured-grid bathymetry is, by definition, the
``b = sea_level`` contour.  We extract it with matplotlib's contouring,
then simplify (Douglas-Peucker-lite: drop vertices closer than a
tolerance) and classify each closed path as either the outer boundary
or an island.
"""
from __future__ import annotations

from typing import List, Tuple

import numpy as np


def signed_area(pts: np.ndarray) -> float:
    """Shoelace signed area of a closed 2-D polygon ``pts[N, 2]``.

    Positive for counter-clockwise, negative for clockwise.  The
    polygon is assumed closed (``pts[0] == pts[-1]``); if not, the
    caller should close it first.
    """
    x, y = pts[:, 0], pts[:, 1]
    return 0.5 * float(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]))


def simplify_polygon(pts: np.ndarray, min_dist: float) -> np.ndarray:
    """Drop intermediate vertices closer than ``min_dist`` to the last
    kept vertex.  Endpoints are always preserved.

    This is a cheap linear-time simplification, not Douglas-Peucker.
    It is sufficient for coastline polygons at mesh resolution because
    ``min_dist`` is tuned to the target triangle edge length — a finer
    algorithm would not change the final mesh.
    """
    if len(pts) <= 2:
        return np.asarray(pts)
    keep = [pts[0]]
    for p in pts[1:-1]:
        if np.hypot(p[0] - keep[-1][0], p[1] - keep[-1][1]) >= min_dist:
            keep.append(p)
    keep.append(pts[-1])
    return np.asarray(keep)


def _close_polygon(pts: np.ndarray) -> np.ndarray:
    """Return ``pts`` with its first vertex appended if the polygon is
    not already closed (first != last)."""
    if not np.allclose(pts[0], pts[-1]):
        return np.vstack([pts, pts[0:1]])
    return pts


def extract_coastlines(
    bathy: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    sea_level: float = 0.0,
    simplify_tol: float | None = None,
    min_points: int = 4,
) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Extract the coastline from a structured-grid bathymetry.

    Parameters
    ----------
    bathy : ``(nx, ny)`` float array
        Bed elevation in metres, positive-up (wet cells ``< sea_level``).
    x, y : 1-D float arrays
        Cell-centre coordinates (length ``nx`` and ``ny``).
    sea_level : float
        Contour level that defines the coastline (default 0).
    simplify_tol : float or None
        Drop vertices closer than this (metres).  ``None`` disables
        simplification.
    min_points : int
        Discard contour segments with fewer than this many points
        (they are usually numerical noise).

    Returns
    -------
    outer : ``(M, 2)`` array
        Outer boundary polygon, CCW, closed.  Selected as the contour
        with the largest number of vertices (robust proxy for
        "encloses the domain").
    islands : list of ``(M_i, 2)`` arrays
        Remaining closed contours, each closed.  No winding convention
        is enforced here — the PSLG builder orients them appropriately.
    """
    import matplotlib

    # Headless backend is fine; the user may already have set Agg.
    # matplotlib.contour does not need the backend to be interactive.
    import matplotlib.pyplot as plt

    fig_tmp, ax_tmp = plt.subplots()
    cs = ax_tmp.contour(x, y, bathy.T, levels=[float(sea_level)])
    raw_paths = [seg.copy() for seg in cs.allsegs[0] if len(seg) >= min_points + 1]
    plt.close(fig_tmp)

    if not raw_paths:
        raise RuntimeError(
            f"No coastline contours found at sea_level={sea_level}. "
            "Check the bathymetry sign convention (positive-up expected)."
        )

    # Simplify and re-close each path.
    paths = []
    if simplify_tol is None or simplify_tol <= 0.0:
        for p in raw_paths:
            if len(p) >= min_points:
                paths.append(_close_polygon(p))
    else:
        for p in raw_paths:
            s = simplify_polygon(p, simplify_tol)
            if len(s) < min_points:
                continue
            if np.hypot(*(s[0] - s[-1])) < simplify_tol * 2:
                s = _close_polygon(s)
            paths.append(s)

    if not paths:
        raise RuntimeError(
            "No coastline contours survived simplification. "
            f"Try a smaller simplify_tol (current {simplify_tol})."
        )

    # Pick outer = the longest closed polygon. This is a robust proxy
    # for "encloses the domain" for raster bathymetries where the sea
    # is topologically connected.
    main_idx = max(range(len(paths)), key=lambda i: len(paths[i]))
    outer = paths[main_idx]
    inner = [paths[i] for i in range(len(paths)) if i != main_idx]

    # Force outer to CCW
    if signed_area(outer) < 0:
        outer = outer[::-1]

    return outer, inner


def filter_islands_by_centroid(
    islands: List[np.ndarray],
    outer: np.ndarray,
    bathy: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    min_area: float = 0.0,
) -> List[np.ndarray]:
    """Keep only closed contours whose centroid (a) lies inside
    ``outer`` and (b) sits on a dry raster cell.  Optionally drop
    contours smaller than ``min_area`` (m²).

    Used to reject spurious contour fragments that matplotlib can
    produce at the raster edge or across pinched-off wet regions.
    """
    from matplotlib.path import Path as MplPath

    if not islands:
        return []

    outer_mpl = MplPath(outer)
    dx = float(x[1] - x[0]) if len(x) > 1 else 1.0
    dy = float(y[1] - y[0]) if len(y) > 1 else 1.0
    nx, ny = bathy.shape

    kept = []
    for isl in islands:
        isl_closed = _close_polygon(isl)
        n = len(isl_closed) - 1
        if n < 3:
            continue
        if min_area > 0.0 and abs(signed_area(isl_closed)) < min_area:
            continue
        centroid = isl_closed[:n].mean(axis=0)
        if not outer_mpl.contains_point(centroid):
            continue
        ix = int(np.clip(round((centroid[0] - x[0]) / dx), 0, nx - 1))
        iy = int(np.clip(round((centroid[1] - y[0]) / dy), 0, ny - 1))
        if bathy[ix, iy] > 0:
            kept.append(isl_closed)

    return kept
