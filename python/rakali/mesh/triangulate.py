"""Constrained Delaunay triangulation of coastal polygons.

Wraps ``meshpy.triangle`` with the idioms the Rakali unstructured
solver expects: CCW triangle orientation, 1-based vertex ids, per-edge
boundary tags on the ``(ntri, 3)`` layout, and an API that takes an
outer polygon plus a list of island (hole) polygons directly.
"""
from __future__ import annotations

from typing import Callable, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from ..unstr import UnstructuredMesh
from .coastline import signed_area


# Facet-marker callable signature: (midpoint_x, midpoint_y) -> int
# Default marker value (wall).  0 would work too, but using 1 keeps
# generated meshes consistent with the hand-written examples.
DEFAULT_WALL_MARKER = 1


def build_pslg(
    outer: np.ndarray,
    islands: Sequence[np.ndarray] = (),
    marker_fn: Optional[Callable[[float, float], int]] = None,
):
    """Build a Planar Straight-Line Graph for ``meshpy.triangle``.

    Parameters
    ----------
    outer : ``(N+1, 2)`` array
        Closed outer polygon.  Must be CCW (use ``signed_area`` to
        check).
    islands : sequence of closed polygons
        Hole contours.  Hole markers are generated at each island's
        centroid.
    marker_fn : callable ``(mx, my) -> int`` or None
        Boundary-marker rule applied to each outer-polygon edge.  The
        midpoint of the edge is passed.  ``None`` marks every outer
        edge as :data:`DEFAULT_WALL_MARKER`.  Island edges are always
        marked as walls.

    Returns
    -------
    info : ``meshpy.triangle.MeshInfo``
        Ready to pass to :func:`triangulate`.
    """
    try:
        from meshpy import triangle as mt
    except ImportError as exc:
        raise ImportError(
            "build_pslg requires meshpy. "
            "Install with `conda install -c conda-forge meshpy`."
        ) from exc

    if marker_fn is None:
        def marker_fn(mx, my, _default=DEFAULT_WALL_MARKER):  # noqa: E306
            return _default

    points: List[Tuple[float, float]] = []
    facets: List[Tuple[int, int]] = []
    markers: List[int] = []
    offset = 0

    # Outer
    n_outer = len(outer) - 1
    for i in range(n_outer):
        points.append((float(outer[i, 0]), float(outer[i, 1])))
    for i in range(n_outer):
        j = (i + 1) % n_outer
        mx = 0.5 * (outer[i, 0] + outer[j, 0])
        my = 0.5 * (outer[i, 1] + outer[j, 1])
        facets.append((offset + i, offset + j))
        markers.append(int(marker_fn(float(mx), float(my))))
    offset += n_outer

    # Islands
    hole_points: List[Tuple[float, float]] = []
    for isl in islands:
        n_isl = len(isl) - 1
        if n_isl < 3:
            continue
        for i in range(n_isl):
            points.append((float(isl[i, 0]), float(isl[i, 1])))
        for i in range(n_isl):
            j = (i + 1) % n_isl
            facets.append((offset + i, offset + j))
            markers.append(DEFAULT_WALL_MARKER)
        centroid = isl[:n_isl].mean(axis=0)
        hole_points.append((float(centroid[0]), float(centroid[1])))
        offset += n_isl

    info = mt.MeshInfo()
    info.set_points(points)
    info.set_facets(facets, facet_markers=markers)
    if hole_points:
        info.set_holes(hole_points)
    return info


def triangulate(
    info,
    max_area: float,
    min_angle: float = 25.0,
    refine_regions: Iterable[Tuple[np.ndarray, float]] = (),
) -> UnstructuredMesh:
    """Run Triangle and wrap the result as an :class:`UnstructuredMesh`.

    Parameters
    ----------
    info : ``meshpy.triangle.MeshInfo``
        The PSLG built by :func:`build_pslg`.
    max_area : float
        Global triangle area upper bound (m²).  For a target edge
        length ``L`` on an equilateral triangle, ``max_area ≈ 0.5*L²``
        is a reasonable default (meshpy keeps triangles well inside
        that bound after refinement).
    min_angle : float
        Minimum interior angle (degrees).  ``25`` matches the existing
        hand-written examples; ``30`` is the Triangle default and
        produces slightly better shape quality at the cost of more
        triangles.
    refine_regions : iterable of ``(polygon, max_area)``
        Optional local refinement regions.  Each polygon (closed,
        ``(N+1, 2)``) gets its own ``max_area`` constraint.  Uses
        ``MeshInfo.set_regions``.  Currently unused by the shipped
        examples — exposed for future finer-resolution studies.
    """
    try:
        from meshpy import triangle as mt
    except ImportError as exc:
        raise ImportError(
            "triangulate requires meshpy. "
            "Install with `conda install -c conda-forge meshpy`."
        ) from exc

    refine_regions = list(refine_regions)
    if refine_regions:
        # Region records: (x, y, attribute, max_area).  Attribute is a
        # region id (0, 1, 2, ...) that Triangle carries through but we
        # do not use downstream.
        regions = []
        for rid, (poly, a_max) in enumerate(refine_regions):
            centroid = np.asarray(poly)[:-1].mean(axis=0)
            regions.append(
                (float(centroid[0]), float(centroid[1]), rid, float(a_max))
            )
        info.set_regions(regions)

    mesh = mt.build(
        info,
        max_volume=float(max_area),
        min_angle=float(min_angle),
        generate_faces=True,
    )

    pts = np.asarray(mesh.points, dtype=np.float64)
    elems = np.asarray(mesh.elements, dtype=np.int64)

    vert_x = np.ascontiguousarray(pts[:, 0])
    vert_y = np.ascontiguousarray(pts[:, 1])

    # Flip CW triangles to CCW.
    x_e, y_e = vert_x[elems], vert_y[elems]
    signed = 0.5 * (
        (x_e[:, 1] - x_e[:, 0]) * (y_e[:, 2] - y_e[:, 0])
        - (x_e[:, 2] - x_e[:, 0]) * (y_e[:, 1] - y_e[:, 0])
    )
    cw = signed < 0.0
    if np.any(cw):
        elems[cw] = elems[cw][:, [0, 2, 1]]

    ntri = elems.shape[0]
    tri_v = np.ascontiguousarray((elems + 1).astype(np.int32))

    bc_tags = np.zeros((ntri, 3), dtype=np.int32)
    if hasattr(mesh, "faces") and len(mesh.faces) > 0:
        face_pairs = np.asarray(mesh.faces, dtype=np.int64)
        face_marks = np.asarray(mesh.face_markers, dtype=np.int32)
        edge_marker = {}
        for (a, b), m in zip(face_pairs, face_marks):
            if m != 0:
                edge_marker[frozenset((int(a), int(b)))] = int(m)
        if edge_marker:
            opp = [(1, 2), (2, 0), (0, 1)]
            for k in range(ntri):
                for e, (va_idx, vb_idx) in enumerate(opp):
                    va = int(elems[k, va_idx])
                    vb = int(elems[k, vb_idx])
                    m = edge_marker.get(frozenset((va, vb)), 0)
                    if m != 0:
                        bc_tags[k, e] = m

    return UnstructuredMesh(vert_x, vert_y, tri_v, bc_tags=bc_tags)


def make_side_marker(
    bbox: Tuple[float, float, float, float],
    open_sides: Iterable[str],
    tol: float,
    open_marker: int = 2,
    wall_marker: int = DEFAULT_WALL_MARKER,
) -> Callable[[float, float], int]:
    """Build a ``marker_fn`` that tags edges on the listed domain
    sides as ``open_marker`` and all others as ``wall_marker``.

    ``bbox`` is ``(x0, x1, y0, y1)``.  ``open_sides`` accepts any of
    ``"west"``, ``"east"``, ``"south"``, ``"north"``.  ``tol`` is the
    distance (metres) from a side within which an edge midpoint is
    considered to lie on that side — a cell-size of the underlying
    raster is a sensible default.
    """
    x0, x1, y0, y1 = bbox
    sides = {s.lower() for s in open_sides}
    unknown = sides - {"west", "east", "south", "north"}
    if unknown:
        raise ValueError(
            f"unknown side names: {sorted(unknown)}. "
            "Expected any of west, east, south, north."
        )

    def marker(mx: float, my: float) -> int:
        if "east" in sides and mx > x1 - tol:
            return open_marker
        if "west" in sides and mx < x0 + tol:
            return open_marker
        if "north" in sides and my > y1 - tol:
            return open_marker
        if "south" in sides and my < y0 + tol:
            return open_marker
        return wall_marker

    return marker
