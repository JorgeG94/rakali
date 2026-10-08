"""Cell renumbering for cache locality on the unstructured solver.

The triangle list coming out of meshpy has no useful spatial order — two
adjacent triangles in memory can sit on opposite sides of the domain. On a
GPU that turns every neighbour lookup into an uncoalesced load. Sorting the
triangles along a 2D Hilbert curve so cell ``i`` and its three neighbours
land in nearby cache lines is a free single-GPU optimisation: zero comms,
zero per-step cost, applied once at mesh-build time.

The Hilbert curve preserves spatial locality: cells that are close in
space land close on the curve.

Algorithm
---------
1. Normalise the cell centroids ``(cx, cy)`` to integer coordinates in
   ``[0, 2**ORDER - 1]``.
2. Map each ``(x, y)`` to its 1D distance along the Hilbert curve.
3. ``argsort`` by distance to get a spatially-coherent traversal order.
"""
from __future__ import annotations

import numpy as np

from .builder import CoastalMeshResult

# 16 bits per axis → 65 535 grid points along each direction.  Plenty
# for any realistic coastal mesh; doubling it would only matter for
# basin-scale grids with sub-metre triangles, where you would partition
# with METIS anyway.
_ORDER = 16


def _xy_to_hilbert_d(xi: np.ndarray, yi: np.ndarray, order: int) -> np.ndarray:
    """Vectorised Hilbert-curve distance for integer coordinates.

    Implements the iterative quadrant-rotation algorithm from
    Skilling (2004), "Programming the Hilbert curve", with the
    standard "U" base shape.  Distances are returned as ``int64``
    (sufficient for ``order <= 31``).
    """
    x = xi.astype(np.int64, copy=True)
    y = yi.astype(np.int64, copy=True)
    d = np.zeros_like(x)

    s = 1 << (order - 1)
    while s > 0:
        rx = ((x & s) > 0).astype(np.int64)
        ry = ((y & s) > 0).astype(np.int64)
        d += (s * s) * ((3 * rx) ^ ry)

        # Quadrant-conditional rotate / flip.  Where ry == 0 we rotate
        # the lower half so the curve enters / exits the right way.
        flip = (ry == 0) & (rx == 1)
        x[flip] = s - 1 - x[flip]
        y[flip] = s - 1 - y[flip]
        rotate = ry == 0
        # In-place swap of (x, y) for the rows that need rotating.
        tmp = x[rotate].copy()
        x[rotate] = y[rotate]
        y[rotate] = tmp

        s >>= 1
    return d


def hilbert_keys(cx: np.ndarray, cy: np.ndarray,
                 order: int = _ORDER) -> np.ndarray:
    """Hilbert-curve distance for each centroid, after [lo, hi] normalisation.

    Degenerate axes (constant coordinate) collapse to 0
    so ordering still respects the other axis.
    """
    cx = np.asarray(cx, dtype=np.float64)
    cy = np.asarray(cy, dtype=np.float64)
    n = cx.size
    if cy.size != n:
        raise ValueError("cx and cy must have the same length")
    if n == 0:
        return np.zeros(0, dtype=np.int64)

    scale = (1 << order) - 1
    x_lo, x_hi = cx.min(), cx.max()
    y_lo, y_hi = cy.min(), cy.max()
    x_span = x_hi - x_lo
    y_span = y_hi - y_lo

    xi = np.zeros(n, dtype=np.int64) if x_span == 0.0 else \
         np.clip(((cx - x_lo) / x_span * scale).round().astype(np.int64),
                 0, scale)
    yi = np.zeros(n, dtype=np.int64) if y_span == 0.0 else \
         np.clip(((cy - y_lo) / y_span * scale).round().astype(np.int64),
                 0, scale)

    return _xy_to_hilbert_d(xi, yi, order)



def hilbert_renumber(result: CoastalMeshResult,
                     order: int = 16) -> CoastalMeshResult:
    """Reorder triangles along a 2D Hilbert curve, in place.

    Permutes ``result.mesh.tri_v``, ``result.mesh.bc_tags``, and
    ``result.bed`` so cell ``i`` sits next to its spatial neighbours in
    memory. Vertex coordinates are left alone — only the triangle axis
    moves.

    Parameters
    ----------
    result : CoastalMeshResult
        Mesh result returned by :func:`build_coastal_mesh`. Mutated in
        place; also returned for chaining.
    order : int
        Hilbert curve resolution (``2^order`` cells per side). 16 is
        plenty for any mesh that fits in memory.

    Returns
    -------
    CoastalMeshResult
        The same object, with its triangles reordered.
    """
    mesh = result.mesh
    cx, cy = mesh.centroids()
    keys = hilbert_keys(cx, cy, order=order)
    perm = np.argsort(keys, kind="stable")

    mesh.tri_v = np.ascontiguousarray(mesh.tri_v[perm])
    if mesh.bc_tags is not None:
        mesh.bc_tags = np.ascontiguousarray(mesh.bc_tags[perm])
    result.bed = np.ascontiguousarray(result.bed[perm])
    return result
