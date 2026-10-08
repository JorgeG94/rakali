"""Top-level coastal mesh builder.

``build_coastal_mesh`` composes the primitives in
``rakali.mesh.{bed, coastline, triangulate}`` into a single call that
turns a bathymetry file into a mesh-plus-bed ready for the
unstructured solver.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import (
    Callable,
    Iterable,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import numpy as np

from ..unstr import UnstructuredMesh, UnstructuredSimulation
from .bed import interpolate_bed, load_bathymetry_nc, window_bathymetry
from .coastline import (
    extract_coastlines,
    filter_islands_by_centroid,
    signed_area,
)
from .triangulate import (
    DEFAULT_WALL_MARKER,
    build_pslg,
    make_side_marker,
    triangulate,
)


BathymetrySpec = Union[str, Path, np.ndarray]
OpenBoundarySpec = Union[
    None,
    str,
    Sequence[str],
    Callable[[float, float], int],
]


@dataclass
class CoastalMeshResult:
    """Output of :func:`build_coastal_mesh`.

    Attributes
    ----------
    mesh : UnstructuredMesh
        Triangulation, ready to hand to
        :class:`~rakali.unstr.UnstructuredSimulation`.
    bed : ``(mesh.ntri,)`` float64 array
        Bathymetry sampled at each triangle centroid, metres,
        positive-up.  Same ordering as ``mesh``.
    bathymetry : ``(nx, ny)`` float64 array
        The windowed raster that the mesh was cut from.  Kept so callers
        can re-interpolate onto new points (e.g. when adding gauges) or
        plot a background.
    x, y : 1-D float64 arrays
        Cell-centre coordinates of ``bathymetry``.
    outer : ``(M, 2)`` float64 array
        Outer coastline polygon that became the PSLG outer boundary
        (CCW, closed).
    islands : list of ``(M_i, 2)`` arrays
        Closed island polygons inserted as holes.
    bbox : tuple
        ``(x_min, x_max, y_min, y_max)`` of the windowed raster.
    wall_marker, open_marker : int
        Boundary-group ids used on wall and open edges respectively.
        Pass these to :meth:`UnstructuredSimulation.set_bc_type`.
    """

    mesh: UnstructuredMesh
    bed: np.ndarray
    bathymetry: np.ndarray
    x: np.ndarray
    y: np.ndarray
    outer: np.ndarray
    islands: List[np.ndarray]
    bbox: Tuple[float, float, float, float]
    wall_marker: int = DEFAULT_WALL_MARKER
    open_marker: int = 2
    metadata: dict = field(default_factory=dict)

    def make_simulation(
        self,
        *,
        bc_wall: str = "wall",
        bc_open: str = "tidal",
        cfl: float = 0.45,
        manning: float = 0.025,
        initial_stage: Union[np.ndarray, float, None] = None,
    ) -> UnstructuredSimulation:
        """Build and configure an :class:`UnstructuredSimulation` with
        sensible defaults for a coastal run.

        Wires up:

        * ``bc_wall`` (default ``"wall"``) on the wall boundary group.
        * ``bc_open`` (default ``"tidal"``) on the open boundary group.
          Pass ``"open"`` for a Flather radiation BC, ``"inflow"`` for
          a prescribed normal velocity, etc.
        * CFL number (default 0.45) and Manning friction (default
          0.025) — typical values for tidal and river applications.
        * Bathymetry interpolated at triangle centroids.
        * Initial sea-surface elevation: ``max(bed, 0)`` (still water
          at z=0, dry cells pinned to the bed) unless you pass
          ``initial_stage``.

        Parameters
        ----------
        bc_wall, bc_open : str
            BC type names for the wall / open groups.  Any value
            accepted by :meth:`UnstructuredSimulation.set_bc_type`.
        cfl : float
            CFL number.
        manning : float
            Manning roughness coefficient.
        initial_stage : array or scalar, optional
            Override the default still-water IC.  A scalar is
            broadcast across all cells; an array must match the
            owned-cell count.

        Returns
        -------
        UnstructuredSimulation
            Configured but not yet prepared — stepping will auto-call
            :meth:`~UnstructuredSimulation.prepare` on first use, or
            call it manually if you want to time setup separately.
            Use as a context manager to release the C handle on exit.
        """
        sim = UnstructuredSimulation(self.mesh)
        bed = self.bed

        sim.set_bc_type(group=self.wall_marker, bc_type=bc_wall)
        sim.set_bc_type(group=self.open_marker, bc_type=bc_open)
        sim.set_cfl(cfl)
        sim.set_manning(manning)

        sim.set_bathymetry(bed)
        if initial_stage is None:
            stage0 = np.maximum(bed, 0.0)
        elif np.isscalar(initial_stage):
            stage0 = np.full_like(bed, float(initial_stage))
        else:
            stage0 = np.asarray(initial_stage, dtype=np.float64)
        sim.set_initial_stage(stage0)

        return sim


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def build_coastal_mesh(
    bathymetry: BathymetrySpec,
    resolution: float,
    *,
    dx: Optional[float] = None,
    dy: Optional[float] = None,
    window: Optional[Tuple[float, float, float, float]] = None,
    sea_level: float = 0.0,
    open_boundary: OpenBoundarySpec = None,
    open_marker: int = 2,
    wall_marker: int = DEFAULT_WALL_MARKER,
    force_wet_edges: Sequence[str] = (),
    wet_depth: float = -50.0,
    min_island_area: float = 0.0,
    simplify_tol: Optional[float] = None,
    min_angle: float = 25.0,
    max_area: Optional[float] = None,
    refine_regions: Iterable[Tuple[np.ndarray, float]] = (),
    interp_method: str = "nearest",
    bathymetry_var: Optional[str] = None,
    verbose: bool = True,
    _log: Optional[Callable[[str], None]] = None,
) -> CoastalMeshResult:
    """Build a coastal triangular mesh from a raster bathymetry.

    Parameters
    ----------
    bathymetry : str, Path, or ``(nx, ny)`` ndarray
        Path to a NetCDF file (schema: dimensions ``x``, ``y``;
        variable ``b`` / ``elevation`` / ``depth``) or a pre-loaded
        array.  Positive-up sign convention: ``b < sea_level`` is wet.
    resolution : float
        Target triangle edge length (m).  Drives ``simplify_tol`` and
        ``max_area`` when those are not given.
    dx, dy : float, optional
        Raster cell size in metres.  Required.  ``dy`` defaults to
        ``dx`` (square cells, matches every shipped example).
    window : ``(x0, x1, y0, y1)`` or None
        Crop the raster to this rectangle (metres).  ``None`` uses the
        full raster.
    sea_level : float
        Bathymetry contour treated as the coastline.  Default 0.
    open_boundary :
        Specifies which outer-polygon edges carry the ``open_marker``
        boundary-group id.  Accepts:

        * ``None``        — every edge is a wall.
        * ``"east"`` / ``"west"`` / ``"north"`` / ``"south"`` — all
          edges whose midpoint lies within one cell of the named
          domain side are open; matches the existing examples.
        * a list of side names — any of the above.
        * a callable ``(mx, my) -> int`` — invoked per-edge midpoint;
          the return value is the marker id.  Escape hatch for
          non-rectangular domains.
    open_marker, wall_marker : int
        Integer group ids.  ``open_marker`` defaults to 2 and
        ``wall_marker`` to 1, matching the conventions in the shipped
        examples.
    force_wet_edges : sequence of side names
        Force the outermost cell row along each named side to at least
        ``wet_depth`` metres deep **before** contour extraction.  This
        exists because open BCs need wet cells at the boundary; on
        domains where the bathymetry rasters as dry at the edge, the
        contour would close the coastline along that edge and there is
        no open boundary to tag.  Applied as
        ``b[edge] = min(b[edge], wet_depth)``.  Default empty (off).
    wet_depth : float
        Depth (negative value) used by ``force_wet_edges``.
    min_island_area : float
        Drop detected islands smaller than this (m²).  Default 0 keeps
        all valid islands (those with a dry centroid inside the outer
        polygon).
    simplify_tol : float, optional
        Minimum distance between adjacent polygon vertices after
        simplification (m).  Default ``0.7 * resolution`` — tight
        enough that the mesh resolves every curve it needs to.
    min_angle : float
        Minimum interior triangle angle, degrees.  Default 25.
    max_area : float, optional
        Upper bound on triangle area (m²).  Default
        ``0.5 * resolution**2`` (equilateral triangle with edge
        ``resolution`` has area ``~0.43 * L²``).
    refine_regions : iterable of ``(polygon, max_area)``
        Local refinement regions.  Currently unused by shipped
        examples; wired through for future studies.
    interp_method : str
        Passed to ``scipy.interpolate.RegularGridInterpolator``.
        ``"nearest"`` is the default (safe for step-feature bathymetry);
        ``"linear"`` smooths.
    bathymetry_var : str, optional
        NetCDF variable name override.  Default auto-detects.
    verbose : bool
        Print progress messages.  Disable for quiet test runs.

    Returns
    -------
    CoastalMeshResult
    """
    log = _log if _log is not None else (print if verbose else (lambda *_: None))

    if dx is None:
        raise ValueError(
            "build_coastal_mesh: dx is required "
            "(the bathymetry file does not carry cell spacing)."
        )
    if dy is None:
        dy = dx

    # -- Load raster ---------------------------------------------------
    if isinstance(bathymetry, (str, Path)):
        log(f"[mesh] loading bathymetry: {bathymetry}")
        bathy_full = load_bathymetry_nc(bathymetry, varname=bathymetry_var)
    else:
        bathy_full = np.asarray(bathymetry, dtype=np.float64)
        if bathy_full.ndim != 2:
            raise ValueError("bathymetry array must be 2-D")

    # -- Window --------------------------------------------------------
    bathy, x, y = window_bathymetry(bathy_full, dx=dx, dy=dy, window=window)
    nx_sub, ny_sub = bathy.shape
    n_wet = int((bathy <= sea_level).sum())
    log(
        f"[mesh] window: {nx_sub}x{ny_sub} cells "
        f"({n_wet} wet @ sea_level={sea_level})"
    )
    bbox = (float(x[0]), float(x[-1]), float(y[0]), float(y[-1]))

    # -- Force wet edges ----------------------------------------------
    bathy_work = bathy.copy()
    if force_wet_edges:
        sides = {s.lower() for s in force_wet_edges}
        if "east" in sides:
            bathy_work[-1, :] = np.minimum(bathy_work[-1, :], wet_depth)
        if "west" in sides:
            bathy_work[0, :] = np.minimum(bathy_work[0, :], wet_depth)
        if "north" in sides:
            bathy_work[:, -1] = np.minimum(bathy_work[:, -1], wet_depth)
        if "south" in sides:
            bathy_work[:, 0] = np.minimum(bathy_work[:, 0], wet_depth)

    # Pad with a ring of land so the contour at the raster boundary
    # closes into proper polygons.  The padded cells are never in the
    # mesh — they only serve to trap the contour extraction.
    b_pad = np.full((nx_sub + 2, ny_sub + 2), 100.0)
    b_pad[1:-1, 1:-1] = bathy_work
    x_pad = np.concatenate([[x[0] - dx], x, [x[-1] + dx]])
    y_pad = np.concatenate([[y[0] - dy], y, [y[-1] + dy]])

    # -- Extract coastline --------------------------------------------
    tol = simplify_tol if simplify_tol is not None else 0.7 * resolution
    log(f"[mesh] extracting coastline (simplify_tol={tol:.1f}m)")
    outer, islands_raw = extract_coastlines(
        b_pad, x_pad, y_pad, sea_level=sea_level, simplify_tol=tol,
    )
    islands = filter_islands_by_centroid(
        islands_raw, outer, bathy_work, x, y, min_area=min_island_area,
    )
    log(
        f"[mesh] outer: {len(outer) - 1} pts, "
        f"area={abs(signed_area(outer)) / 1e6:.2f} km²; "
        f"islands: {len(islands)} (of {len(islands_raw)} raw)"
    )

    # -- Boundary-marker callable --------------------------------------
    if open_boundary is None:
        marker_fn = None
    elif callable(open_boundary):
        marker_fn = open_boundary
    else:
        if isinstance(open_boundary, str):
            sides = [open_boundary]
        else:
            sides = list(open_boundary)
        marker_fn = make_side_marker(
            bbox, sides, tol=max(dx, dy),
            open_marker=open_marker, wall_marker=wall_marker,
        )

    # -- Triangulate --------------------------------------------------
    if max_area is None:
        max_area = 0.5 * resolution ** 2
    info = build_pslg(outer, islands, marker_fn=marker_fn)
    log(
        f"[mesh] triangulating: max_area={max_area:.1f} m², "
        f"min_angle={min_angle:.1f}°"
    )
    mesh = triangulate(
        info, max_area=max_area, min_angle=min_angle,
        refine_regions=refine_regions,
    )

    n_wall = int(np.sum(mesh.bc_tags == wall_marker))
    n_open = int(np.sum(mesh.bc_tags == open_marker))
    log(
        f"[mesh] built: {mesh.ntri} triangles, {mesh.nverts} vertices; "
        f"{n_wall} wall edges, {n_open} open edges"
    )

    # -- Sample bed at centroids --------------------------------------
    cx, cy = mesh.centroids()
    bed = interpolate_bed(bathy, x, y, cx, cy, method=interp_method)

    return CoastalMeshResult(
        mesh=mesh,
        bed=bed,
        bathymetry=bathy,
        x=x, y=y,
        outer=outer,
        islands=islands,
        bbox=bbox,
        wall_marker=wall_marker,
        open_marker=open_marker,
        metadata={
            "n_wall_edges": n_wall,
            "n_open_edges": n_open,
            "resolution": resolution,
            "max_area": max_area,
            "simplify_tol": tol,
        },
    )
