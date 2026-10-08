"""Coastal mesh construction for the unstructured Rakali solver.

This subpackage turns a NetCDF bathymetry file into a Delaunay-refined
triangular mesh with boundary markers, ready to hand to
``UnstructuredSimulation``.  The intended entry point is
:func:`build_coastal_mesh`; the lower-level primitives in
``coastline``, ``triangulate`` and ``bed`` are exposed for users who
need to step outside the standard flow.

Example
-------
>>> from rakali.mesh import build_coastal_mesh
>>> from rakali.unstr import UnstructuredSimulation
>>> result = build_coastal_mesh(
...     "bathymetry.nc",
...     dx=25.0, resolution=100.0, open_boundary="east",
... )
>>> sim = UnstructuredSimulation(result.mesh)
>>> sim.set_bathymetry(result.bed)
"""
from .builder import (
    CoastalMeshResult,
    build_coastal_mesh,
)
from .coastline import (
    extract_coastlines,
    simplify_polygon,
    signed_area,
)
from .triangulate import build_pslg, triangulate
from .bed import interpolate_bed, load_bathymetry_nc
from .renumber import hilbert_renumber

__all__ = [
    "CoastalMeshResult",
    "build_coastal_mesh",
    "extract_coastlines",
    "simplify_polygon",
    "signed_area",
    "build_pslg",
    "triangulate",
    "interpolate_bed",
    "load_bathymetry_nc",
    "hilbert_renumber",
]
