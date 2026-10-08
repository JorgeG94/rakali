"""Python interface to the unstructured Rakali solver.

Two classes:

* ``UnstructuredMesh`` -- thin wrapper around ``meshpy.triangle`` that
  builds a Delaunay-refined triangulation from a polygon, hands back
  numpy arrays for vertices, triangles, and per-edge boundary tags.
  No solver state -- this is purely a preprocessing helper.

* ``UnstructuredSimulation`` -- ctypes binding to ``librakali.so``'s
  ``rki_unstr_*`` entrypoints. Hands a mesh across the C ABI, accepts
  per-cell bathymetry / initial state from numpy, runs the solver, and
  exposes zero-copy numpy views of the result fields.

Deliberate non-imports: ``netCDF4``, ``mpi4py``, ``h5py``. The Python
process that loads ``librakali.so`` must not pull in any library that
would shadow the NVHPC-built NetCDF/HDF5 stack the solver itself
links against. Post-processing (plotting, file I/O) belongs in
*separate* scripts that import ``netCDF4``/``xarray`` after the solver
has finished and unloaded.

Example -- channel dam break, walls on the long sides, open ends::

    from rakali.unstr import UnstructuredMesh, UnstructuredSimulation
    import numpy as np

    # 100m x 2m channel, ~0.5m max area triangles
    mesh = UnstructuredMesh.rectangle(
        x0=0.0, y0=0.0, x1=100.0, y1=2.0,
        max_area=0.5,
        bc_tags={"left": 2, "right": 2, "bottom": 1, "top": 1},
    )
    # group 1 -> wall (default), group 2 -> open
    sim = UnstructuredSimulation(mesh)
    sim.set_bc_type(group=2, bc_type="open")

    # Initial dam: stage = 2 on the left half, 1 on the right
    cx = sim.centroid_x
    stage = np.where(cx < 50.0, 2.0, 1.0)
    sim.set_bathymetry(np.zeros_like(stage))
    sim.set_initial_stage(stage)

    sim.advance(0.1)          # auto-calls prepare() on first use

    print(f"t = {sim.time:.4f}, steps = {sim.step_count}")
    h = sim.height           # zero-copy numpy view
"""
from __future__ import annotations

import ctypes
import numpy as np
from typing import Mapping, Optional

from rakali._core import _find_library
from rakali._lifecycle import (_SimBase, _owned_view, _require_file,
                               _require_writable_dir)


# ---------------------------------------------------------------------------
# BC type integer constants -- mirror src/core/rki_boundary_unstr.F90
# ---------------------------------------------------------------------------

BC_INTERIOR = 0
BC_WALL = 1
BC_OPEN = 2
BC_TIDAL = 3
BC_INFLOW = 5
BC_CLAMPED = 7
BC_CHAPMAN = 9
BC_OUTFLOW = 10

_BC_NAME_TO_INT = {
    "interior": BC_INTERIOR,
    "wall": BC_WALL,
    "open": BC_OPEN,
    "inflow": BC_INFLOW,
    "outflow": BC_OUTFLOW,
    "tidal": BC_TIDAL,
    "clamped": BC_CLAMPED,
    "chapman": BC_CHAPMAN,
}


def _bc_to_int(bc) -> int:
    if isinstance(bc, str):
        try:
            return _BC_NAME_TO_INT[bc.lower()]
        except KeyError:
            raise ValueError(
                f"Unknown BC name {bc!r}; expected one of "
                f"{sorted(_BC_NAME_TO_INT)}"
            )
    return int(bc)


# ---------------------------------------------------------------------------
# Library loading
# ---------------------------------------------------------------------------

_lib = None


def _load():
    """Load librakali.so on first use and bind unstructured signatures."""
    global _lib
    if _lib is not None:
        return _lib
    lib = ctypes.CDLL(_find_library())

    c_int = ctypes.c_int
    c_double = ctypes.c_double
    c_void_p = ctypes.c_void_p
    p_int = ctypes.POINTER(c_int)
    p_double = ctypes.POINTER(c_double)
    p_void = ctypes.POINTER(c_void_p)

    # Every per-simulation entrypoint takes the opaque handle as its
    # FIRST argument (by value, c_void_p).  rki_unstr_create / destroy
    # mutate the caller's handle slot via POINTER(c_void_p).  set_log_level
    # stays handle-less because the logger is process-global.
    HANDLE = c_void_p
    HANDLE_OUT = ctypes.POINTER(c_void_p)

    # Solver working precision (bytes): pick the numpy dtype + ctypes
    # scalar for the zero-copy field views so a single-precision build
    # (wp=real32) is read as float32 rather than mis-read as float64.
    # `rki_working_precision` is a process-global symbol in the shared
    # librakali (defined in the structured rki_api), so it resolves here
    # too.  Setters stay c_double (Fortran converts to wp internally).
    lib.rki_working_precision.argtypes = []
    lib.rki_working_precision.restype = c_int
    _wp_bytes = lib.rki_working_precision()
    lib._rki_wp_ctype = ctypes.c_float if _wp_bytes == 4 else c_double
    lib._rki_wp_dtype = np.float32 if _wp_bytes == 4 else np.float64
    lib.rki_unstr_set_log_level.argtypes = [c_int]
    lib.rki_unstr_set_log_level.restype = c_int

    # ----- lifecycle -----
    lib.rki_unstr_create.argtypes = [
        c_int, c_int,         # nverts, ntri
        p_double, p_double,   # vert_x, vert_y
        p_int,                # tri_v(3, ntri), 1-based, column-major
        c_void_p,             # bc_tags(3, ntri), nullable
        HANDLE_OUT,           # &handle (out)
    ]
    lib.rki_unstr_create.restype = c_int

    lib.rki_unstr_finalize.argtypes = [HANDLE]
    lib.rki_unstr_finalize.restype = c_int

    lib.rki_unstr_destroy.argtypes = [HANDLE_OUT]
    lib.rki_unstr_destroy.restype = None

    # ----- stepping -----
    lib.rki_unstr_step.argtypes = [HANDLE, c_double]
    lib.rki_unstr_step.restype = c_int
    lib.rki_unstr_advance.argtypes = [HANDLE, c_double]
    lib.rki_unstr_advance.restype = c_int

    # ----- queries -----
    lib.rki_unstr_get_time.argtypes = [HANDLE]
    lib.rki_unstr_get_time.restype = c_double
    lib.rki_unstr_get_dt.argtypes = [HANDLE]
    lib.rki_unstr_get_dt.restype = c_double
    lib.rki_unstr_get_step_count.argtypes = [HANDLE]
    lib.rki_unstr_get_step_count.restype = c_int
    lib.rki_unstr_get_ncells_owned.argtypes = [HANDLE]
    lib.rki_unstr_get_ncells_owned.restype = c_int
    lib.rki_unstr_get_ncells_total.argtypes = [HANDLE]
    lib.rki_unstr_get_ncells_total.restype = c_int
    lib.rki_unstr_get_mass_error.argtypes = [HANDLE]
    lib.rki_unstr_get_mass_error.restype = c_double

    # ----- setters -----
    # Array setters take the buffer as a void* (the Fortran side
    # rebinds it via c_f_pointer using the cell count from the handle).
    # Passing a numpy POINTER(c_double) auto-converts to c_void_p.
    for name in (
        "rki_unstr_set_bathymetry",
        "rki_unstr_set_initial_stage",
        "rki_unstr_set_initial_xmom",
        "rki_unstr_set_initial_ymom",
        "rki_unstr_set_heat_flux",
        "rki_unstr_set_evap_minus_precip",
    ):
        f = getattr(lib, name)
        f.argtypes = [HANDLE, c_void_p]
        f.restype = c_int

    lib.rki_unstr_set_bc_type.argtypes = [HANDLE, c_int, c_int]
    lib.rki_unstr_set_bc_type.restype = c_int
    for name in (
        "rki_unstr_set_bc_eta_ext",
        "rki_unstr_set_bc_un_ext",
        "rki_unstr_set_bc_h_ext",
        "rki_unstr_set_bc_qn_ext",
    ):
        f = getattr(lib, name)
        f.argtypes = [HANDLE, c_int, c_double]
        f.restype = c_int
    lib.rki_unstr_set_manning.argtypes = [HANDLE, c_double]
    lib.rki_unstr_set_manning.restype = c_int
    lib.rki_unstr_set_time_scheme.argtypes = [HANDLE, c_void_p, c_int]
    lib.rki_unstr_set_time_scheme.restype = c_int
    lib.rki_unstr_set_eos.argtypes = [HANDLE, c_double, c_double, c_double,
                                      c_double, c_double]
    lib.rki_unstr_set_eos.restype = c_int
    lib.rki_unstr_set_bpg_jacobian.argtypes = [HANDLE, c_int]
    lib.rki_unstr_set_bpg_jacobian.restype = c_int
    lib.rki_unstr_set_vcoord.argtypes = [HANDLE, c_int, c_int, c_double, c_double]
    lib.rki_unstr_set_vcoord.restype = c_int
    lib.rki_unstr_set_zstar_full_params.argtypes = [HANDLE, c_double, c_double, c_int, c_int]
    lib.rki_unstr_set_zstar_full_params.restype = c_int
    lib.rki_unstr_set_cfl.argtypes = [HANDLE, c_double]
    lib.rki_unstr_set_cfl.restype = c_int
    lib.rki_unstr_set_betas.argtypes = [HANDLE, c_double, c_double, c_double]
    lib.rki_unstr_set_betas.restype = c_int
    lib.rki_unstr_set_wind.argtypes = [HANDLE, c_double, c_double]
    lib.rki_unstr_set_wind.restype = c_int
    lib.rki_unstr_set_wind_field.argtypes = [HANDLE, c_void_p, c_void_p]
    lib.rki_unstr_set_wind_field.restype = c_int
    lib.rki_unstr_set_initial_salinity.argtypes = [HANDLE, c_double]
    lib.rki_unstr_set_initial_salinity.restype = c_int
    lib.rki_unstr_set_initial_temperature.argtypes = [HANDLE, c_double]
    lib.rki_unstr_set_initial_temperature.restype = c_int
    for name in (
        "rki_unstr_set_bc_salinity",
        "rki_unstr_set_bc_temperature",
    ):
        f = getattr(lib, name)
        f.argtypes = [HANDLE, c_int, c_double]
        f.restype = c_int
    lib.rki_unstr_set_smagorinsky_cs.argtypes = [HANDLE, c_double]
    lib.rki_unstr_set_smagorinsky_cs.restype = c_int
    lib.rki_unstr_set_hdiff_kappa.argtypes = [HANDLE, c_double]
    lib.rki_unstr_set_hdiff_kappa.restype = c_int
    lib.rki_unstr_set_use_ppm_tracer.argtypes = [HANDLE, c_int]
    lib.rki_unstr_set_use_ppm_tracer.restype = c_int
    lib.rki_unstr_set_vmix.argtypes = [HANDLE, c_double, c_double, c_double,
                                       c_double, c_double]
    lib.rki_unstr_set_vmix.restype = c_int
    lib.rki_unstr_set_kpp.argtypes = [HANDLE, c_int, c_double, c_double,
                                      c_double]
    lib.rki_unstr_set_kpp.restype = c_int

    # ----- field getters -----
    for name in (
        "rki_unstr_get_stage_ptr",
        "rki_unstr_get_xmom_ptr",
        "rki_unstr_get_ymom_ptr",
        "rki_unstr_get_bed_ptr",
        "rki_unstr_get_height_ptr",
        "rki_unstr_get_centroid_x_ptr",
        "rki_unstr_get_centroid_y_ptr",
        "rki_unstr_get_area_ptr",
    ):
        f = getattr(lib, name)
        f.argtypes = [HANDLE, p_void, p_int]
        f.restype = c_int

    # ----- multilayer (Phase 1b) -----
    lib.rki_unstr_enable_multilayer.argtypes = [HANDLE, c_int]
    lib.rki_unstr_enable_multilayer.restype = c_int
    lib.rki_unstr_get_nz_ml.argtypes = [HANDLE]
    lib.rki_unstr_get_nz_ml.restype = c_int
    for name in (
        "rki_unstr_get_h_layer_ptr",
        "rki_unstr_get_hu_layer_ptr",
        "rki_unstr_get_hv_layer_ptr",
        "rki_unstr_get_hS_layer_ptr",
        "rki_unstr_get_hT_layer_ptr",
    ):
        f = getattr(lib, name)
        f.argtypes = [HANDLE, p_void, p_int, p_int]
        f.restype = c_int

    # ----- restart -----
    lib.rki_unstr_save_restart.argtypes = [HANDLE, c_void_p, c_int]
    lib.rki_unstr_save_restart.restype = c_int
    lib.rki_unstr_load_restart.argtypes = [HANDLE, c_void_p, c_int]
    lib.rki_unstr_load_restart.restype = c_int

    # ----- gauges -----
    lib.rki_unstr_gauges_init.argtypes = [HANDLE, c_void_p, c_int, c_int,
                                          p_double, p_double]
    lib.rki_unstr_gauges_init.restype = c_int
    lib.rki_unstr_gauges_write.argtypes = [HANDLE]
    lib.rki_unstr_gauges_write.restype = c_int
    lib.rki_unstr_gauges_close.argtypes = [HANDLE]
    lib.rki_unstr_gauges_close.restype = c_int

    # ----- UGRID output -----
    lib.rki_unstr_output_init.argtypes = [HANDLE, c_void_p, c_int, c_int]
    lib.rki_unstr_output_init.restype = c_int
    lib.rki_unstr_output_write.argtypes = [HANDLE]
    lib.rki_unstr_output_write.restype = c_int
    lib.rki_unstr_output_close.argtypes = [HANDLE]
    lib.rki_unstr_output_close.restype = c_int

    _lib = lib
    return lib


def _check(ierr: int, what: str):
    if ierr != 0:
        raise RuntimeError(f"{what} failed with ierr={ierr}")


# pic_logger numeric levels (mirror of pic_logger.f90 constants).
# Lower number = quieter.
LOG_ERROR = 5
LOG_WARNING = 6
LOG_PERFORMANCE = 7
LOG_INFO = 8       # default
LOG_VERBOSE = 9
LOG_DEBUG = 10


def set_log_level(level: int) -> None:
    """Configure the global pic_logger verbosity for the unstructured
    solver. Use one of the ``LOG_*`` constants in this module
    (``LOG_WARNING`` to silence the per-test mesh/state INFO chatter).
    Safe to call before any simulation has been created.
    """
    lib = _load()
    _check(lib.rki_unstr_set_log_level(int(level)), "rki_unstr_set_log_level")


# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------


_MAX_BC_GROUPS = 32  # rki_boundary_unstr::MAX_BC_GROUPS


class UnstructuredMesh:
    """A triangular mesh ready to hand to ``UnstructuredSimulation``.

    Stores numpy arrays in the layout the C ABI expects:

    * ``vert_x``, ``vert_y`` -- shape ``(nverts,)``, ``float64``.
    * ``tri_v`` -- shape ``(ntri, 3)``, ``int32``, **1-based** vertex
      ids in CCW order. The numpy C-contiguous ``(ntri, 3)`` layout has
      the same byte order as Fortran column-major ``(3, ntri)``, which
      is what the C ABI expects. Do not transpose: ``arr[k, i]`` in
      numpy is the same memory as ``tri_v(i+1, k+1)`` in Fortran.
    * ``bc_tags`` -- shape ``(ntri, 3)``, ``int32``. Same layout
      convention. Boundary group id per edge of each cell, zero on
      interior edges. Edge ``i`` is opposite local vertex ``i``
      (ANUGA convention).

    The default constructor takes raw arrays already in this layout.
    Use the classmethod ``rectangle`` for the common case of a
    rectangular domain with named edges.
    """

    def __init__(
        self,
        vert_x: np.ndarray,
        vert_y: np.ndarray,
        tri_v: np.ndarray,
        bc_tags: Optional[np.ndarray] = None,
    ):
        vx = np.ascontiguousarray(vert_x, dtype=np.float64)
        vy = np.ascontiguousarray(vert_y, dtype=np.float64)
        if vx.ndim != 1 or vy.ndim != 1 or vx.shape != vy.shape:
            raise ValueError("vert_x, vert_y must be 1D arrays of equal length")

        tv = np.ascontiguousarray(tri_v, dtype=np.int32)
        if tv.ndim != 2 or tv.shape[1] != 3:
            raise ValueError("tri_v must have shape (ntri, 3)")
        if tv.size and (tv.min() < 1 or tv.max() > vx.size):
            raise ValueError(
                f"tri_v is 1-based and must lie in [1, {vx.size}]; got "
                f"[{tv.min()}, {tv.max()}] (0-based connectivity from "
                f"scipy/meshio needs +1)")

        if bc_tags is None:
            bt = np.zeros_like(tv, dtype=np.int32)
        else:
            bt = np.ascontiguousarray(bc_tags, dtype=np.int32)
            if bt.shape != tv.shape:
                raise ValueError("bc_tags must have the same shape as tri_v")
            if bt.size and (bt.min() < 0 or bt.max() > _MAX_BC_GROUPS):
                raise ValueError(
                    f"bc_tags must lie in [0, {_MAX_BC_GROUPS}]; got "
                    f"[{bt.min()}, {bt.max()}]")

        self.vert_x = vx
        self.vert_y = vy
        self.tri_v = tv
        self.bc_tags = bt

    @property
    def nverts(self) -> int:
        return int(self.vert_x.size)

    @property
    def ntri(self) -> int:
        return int(self.tri_v.shape[0])

    # ---------- geometry helpers ----------

    def centroids(self) -> tuple[np.ndarray, np.ndarray]:
        """Per-triangle centroid x and y coordinates as float64 arrays."""
        # tri_v is 1-based; convert to 0-based for numpy indexing.
        idx = self.tri_v.astype(np.int64) - 1
        cx = self.vert_x[idx].mean(axis=1)
        cy = self.vert_y[idx].mean(axis=1)
        return np.ascontiguousarray(cx), np.ascontiguousarray(cy)

    # ---------- builders ----------

    @classmethod
    def disk(
        cls,
        cx: float,
        cy: float,
        radius: float,
        max_area: float,
        n_boundary: int = 64,
        min_angle: float = 30.0,
        boundary_tag: int = 1,
    ) -> "UnstructuredMesh":
        """Build a Delaunay-refined triangulation of a disk.

        Approximates the circle by a regular ``n_boundary``-gon and
        hands it to ``meshpy.triangle`` for refinement. All boundary
        edges share the single ``boundary_tag`` group id (default 1)
        because the wall normal is the same physical condition all
        the way around -- callers can map that group to whatever BC
        type they like via ``set_bc_type``.

        Increase ``n_boundary`` for a smoother circle (256 gives a
        deviation of ~0.0003 * radius from a perfect circle); the
        default 64 is fine for most validation cases.
        """
        try:
            from meshpy import triangle as mt
        except ImportError as exc:
            raise ImportError(
                "UnstructuredMesh.disk requires meshpy. "
                "Install with `conda install -c conda-forge meshpy`."
            ) from exc

        if radius <= 0.0:
            raise ValueError("disk: radius must be positive")
        if n_boundary < 8:
            raise ValueError("disk: n_boundary must be >= 8")

        thetas = np.linspace(0.0, 2.0 * np.pi, n_boundary, endpoint=False)
        points = [
            (float(cx + radius * np.cos(t)), float(cy + radius * np.sin(t)))
            for t in thetas
        ]
        # Close the loop: facet i connects vertex i to vertex (i+1) mod N.
        facets = [(i, (i + 1) % n_boundary) for i in range(n_boundary)]
        markers = [int(boundary_tag)] * n_boundary

        info = mt.MeshInfo()
        info.set_points(points)
        info.set_facets(facets, facet_markers=markers)

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

        # CCW orientation -- flip any clockwise triangles meshpy emits
        x = vert_x[elems]
        y = vert_y[elems]
        signed = 0.5 * (
            (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0])
            - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
        )
        cw = signed < 0.0
        if np.any(cw):
            elems[cw] = elems[cw][:, [0, 2, 1]]

        ntri = elems.shape[0]
        tri_v = np.ascontiguousarray((elems + 1).astype(np.int32))

        # Boundary tags via the meshpy face_markers, same machinery as
        # the rectangle builder.
        bc_arr = np.zeros((ntri, 3), dtype=np.int32)
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
                    for i, (va_idx, vb_idx) in enumerate(opp):
                        va = int(elems[k, va_idx])
                        vb = int(elems[k, vb_idx])
                        m = edge_marker.get(frozenset((va, vb)), 0)
                        if m != 0:
                            bc_arr[k, i] = m

        return cls(vert_x, vert_y, tri_v, bc_tags=bc_arr)

    @classmethod
    def rectangle(
        cls,
        x0: float,
        y0: float,
        x1: float,
        y1: float,
        max_area: float,
        min_angle: float = 30.0,
        bc_tags: Optional[Mapping[str, int]] = None,
    ) -> "UnstructuredMesh":
        """Build a Delaunay-refined triangulation of an axis-aligned box.

        ``bc_tags`` maps named edges (``"left"``, ``"right"``, ``"bottom"``,
        ``"top"``) to integer group ids. Edges not listed default to 0
        (which the solver treats as the default-wall group).

        Requires ``meshpy.triangle`` (installed via the conda env).
        """
        try:
            from meshpy import triangle as mt
        except ImportError as exc:
            raise ImportError(
                "UnstructuredMesh.rectangle requires meshpy. "
                "Install with `conda install -c conda-forge meshpy`."
            ) from exc

        if x1 <= x0 or y1 <= y0:
            raise ValueError("rectangle: need x1>x0 and y1>y0")

        # CCW corner walk: bottom -> right -> top -> left
        points = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        # Edge i goes from point i to point (i+1)%4
        edge_names = ["bottom", "right", "top", "left"]
        tags = bc_tags or {}
        markers = [int(tags.get(name, 0)) for name in edge_names]
        facets = [(0, 1), (1, 2), (2, 3), (3, 0)]

        info = mt.MeshInfo()
        info.set_points(points)
        info.set_facets(facets, facet_markers=markers)

        mesh = mt.build(
            info,
            max_volume=float(max_area),
            min_angle=float(min_angle),
            generate_faces=True,
        )

        pts = np.asarray(mesh.points, dtype=np.float64)         # (nv, 2)
        elems = np.asarray(mesh.elements, dtype=np.int64)       # (ntri, 3) 0-based

        vert_x = np.ascontiguousarray(pts[:, 0])
        vert_y = np.ascontiguousarray(pts[:, 1])

        # Ensure CCW orientation. Triangle usually emits CCW already, but
        # double-check on the signed area and flip any clockwise triangles
        # so the geometry kernel doesn't compute negative areas.
        x = vert_x[elems]
        y = vert_y[elems]
        signed = 0.5 * (
            (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0])
            - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
        )
        cw = signed < 0.0
        if np.any(cw):
            elems[cw] = elems[cw][:, [0, 2, 1]]

        ntri = elems.shape[0]
        # Numpy C-contiguous (ntri, 3) == Fortran column-major (3, ntri).
        # Pass through as-is, just bump to 1-based for Fortran indexing.
        tri_v = np.ascontiguousarray((elems + 1).astype(np.int32))

        # Build per-edge boundary tags from the meshpy face markers.
        # meshpy.triangle exposes boundary edges via mesh.faces and
        # mesh.face_markers as (vert_a, vert_b) pairs. We map each
        # tagged boundary edge to the (cell, local_edge) it belongs to.
        bc_arr = np.zeros((ntri, 3), dtype=np.int32)
        if hasattr(mesh, "faces") and len(mesh.faces) > 0:
            face_pairs = np.asarray(mesh.faces, dtype=np.int64)         # (nf, 2)
            face_marks = np.asarray(mesh.face_markers, dtype=np.int32)  # (nf,)

            # Build a dict {frozenset(va, vb): marker} for tagged boundary edges
            edge_marker = {}
            for (a, b), m in zip(face_pairs, face_marks):
                if m != 0:
                    edge_marker[frozenset((int(a), int(b)))] = int(m)

            if edge_marker:
                # Edge i of cell k is opposite local vertex i; its
                # endpoints are the other two vertices.
                opp = [(1, 2), (2, 0), (0, 1)]
                for k in range(ntri):
                    for i, (va_idx, vb_idx) in enumerate(opp):
                        va = int(elems[k, va_idx])
                        vb = int(elems[k, vb_idx])
                        m = edge_marker.get(frozenset((va, vb)), 0)
                        if m != 0:
                            bc_arr[k, i] = m

        return cls(vert_x, vert_y, tri_v, bc_tags=bc_arr)


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


class UnstructuredSimulation(_SimBase):
    """Drives the unstructured Rakali solver from Python.

    Lifetime: ``__init__`` calls ``rki_unstr_create`` so the C side
    immediately allocates mesh + state. Set bathymetry / IC / BC types,
    then call :meth:`prepare` to move data onto the device (or let
    :meth:`step` / :meth:`advance` auto-call it on first use). After
    that, :meth:`step` and :meth:`advance` run the solver, and the
    field properties give zero-copy numpy views.

    Multi-instance: each ``UnstructuredSimulation`` carries its own
    opaque handle, so several can coexist in one process (and a
    structured ``Rakali`` simulation can run side-by-side without
    interference).

    """

    def __init__(self, mesh):
        self._lib = _load()
        self._destroyed = False
        self._prepared = False
        self._handle = ctypes.c_void_p()

        self._mesh = mesh
        self._create_serial(mesh)

        self._ncells = self._lib.rki_unstr_get_ncells_owned(self._handle)

    def _create_serial(self, mesh: UnstructuredMesh):
        ntri = mesh.ntri
        nverts = mesh.nverts

        vx = mesh.vert_x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        vy = mesh.vert_y.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        tv = mesh.tri_v.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

        bt_ptr = ctypes.c_void_p(0)
        if mesh.bc_tags is not None and np.any(mesh.bc_tags != 0):
            bt_ptr = ctypes.c_void_p(mesh.bc_tags.ctypes.data)

        ierr = self._lib.rki_unstr_create(nverts, ntri, vx, vy, tv, bt_ptr,
                                          ctypes.byref(self._handle))
        _check(ierr, "rki_unstr_create")

    # ---------- destruction ----------

    def _destroy_handle(self):
        """Free the C handle (idempotent via `_destroyed`).

        Teardown entry points — `close()`, `destroy()`, `with` exit and
        `__del__` — come from the shared :class:`_SimBase` mixin so the
        structured and unstructured wrappers behave identically.
        """
        if getattr(self, "_destroyed", True):
            return
        self._lib.rki_unstr_destroy(ctypes.byref(self._handle))
        self._destroyed = True
        self._prepared = False

    # ---------- properties ----------

    @property
    def ncells(self) -> int:
        return self._ncells

    @property
    def time(self) -> float:
        return float(self._lib.rki_unstr_get_time(self._handle))

    @property
    def dt(self) -> float:
        return float(self._lib.rki_unstr_get_dt(self._handle))

    @property
    def step_count(self) -> int:
        return int(self._lib.rki_unstr_get_step_count(self._handle))

    @property
    def mass_error(self) -> float:
        return float(self._lib.rki_unstr_get_mass_error(self._handle))

    # ---------- setters ----------

    def set_bathymetry(self, bed: np.ndarray):
        self._set_field("rki_unstr_set_bathymetry", bed, "bathymetry")

    def set_initial_stage(self, stage: np.ndarray):
        self._set_field("rki_unstr_set_initial_stage", stage, "initial stage")

    def set_initial_xmom(self, xmom: np.ndarray):
        self._set_field("rki_unstr_set_initial_xmom", xmom, "initial xmom")

    def set_initial_ymom(self, ymom: np.ndarray):
        self._set_field("rki_unstr_set_initial_ymom", ymom, "initial ymom")

    def set_heat_flux(self, q: np.ndarray):
        """Net surface heat flux (W/m^2, positive into ocean) per owned cell.

        Applied to the top layer's temperature tracer every RK2 stage:
        ``hT(top) += dt * Q / (rho_0 * cp)``.
        """
        self._set_field("rki_unstr_set_heat_flux", q, "heat flux")

    def set_evap_minus_precip(self, ep: np.ndarray):
        """Virtual surface freshwater flux (m/s, positive for net evaporation)
        per owned cell.

        Applied as a virtual salt flux on the top layer:
        ``hS(top) += dt * (E - P) * S_top``. Does not change ``h_layer``.
        """
        self._set_field("rki_unstr_set_evap_minus_precip", ep, "evap_minus_precip")

    def _set_field(self, sym: str, arr: np.ndarray, what: str):
        a = np.ascontiguousarray(arr, dtype=np.float64)
        if a.size != self._ncells:
            raise ValueError(
                f"{what}: expected length {self._ncells}, got {a.size}"
            )
        # Stash a reference so the buffer outlives the call.
        ptr = a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        ierr = getattr(self._lib, sym)(self._handle, ptr)
        _check(ierr, sym)

    def set_bc_type(self, group: int, bc_type):
        ierr = self._lib.rki_unstr_set_bc_type(
            self._handle, int(group), _bc_to_int(bc_type))
        _check(ierr, "rki_unstr_set_bc_type")

    def set_bc_eta_ext(self, group: int, eta: float):
        """External free-surface elevation (m) for a BC_TIDAL group."""
        _check(self._lib.rki_unstr_set_bc_eta_ext(
                   self._handle, int(group), float(eta)),
               "rki_unstr_set_bc_eta_ext")

    def set_bc_un_ext(self, group: int, un: float):
        """External outward-normal velocity (m/s) for a BC_INFLOW group.

        Sign convention: along the outward normal, so negative means
        flow into the domain. Most users want :meth:`set_bc_inflow_speed`
        instead, which takes a positive number for inflow.
        """
        _check(self._lib.rki_unstr_set_bc_un_ext(
                   self._handle, int(group), float(un)),
               "rki_unstr_set_bc_un_ext")

    def set_bc_inflow_speed(self, group: int, speed: float):
        """Friendlier wrapper for ``set_bc_un_ext``: positive ``speed``
        means flow INTO the domain, regardless of which boundary the
        group sits on. Internally negates and forwards.
        """
        self.set_bc_un_ext(group, -float(speed))

    def set_bc_h_ext(self, group: int, h: float):
        """External depth (m) for a BC_OUTFLOW group."""
        _check(self._lib.rki_unstr_set_bc_h_ext(
                   self._handle, int(group), float(h)),
               "rki_unstr_set_bc_h_ext")

    def set_bc_qn_ext(self, group: int, qn: float):
        """External specific discharge (m^2/s) for a BC_DISCHARGE group.

        Sign convention: along the outward normal, so negative means
        flow into the domain. Most users want :meth:`set_bc_discharge`
        instead, which takes a positive number for inflow.
        """
        _check(self._lib.rki_unstr_set_bc_qn_ext(
                   self._handle, int(group), float(qn)),
               "rki_unstr_set_bc_qn_ext")

    def set_bc_discharge(self, group: int, qn: float):
        """Friendlier wrapper for ``set_bc_qn_ext``: positive ``qn``
        means specific discharge INTO the domain (m^2/s), regardless
        of which boundary the group sits on. Internally negates and
        forwards.
        """
        self.set_bc_qn_ext(group, -float(qn))

    def set_manning(self, n: float):
        _check(self._lib.rki_unstr_set_manning(self._handle, float(n)),
               "rki_unstr_set_manning")

    def set_time_scheme(self, scheme: str):
        """Select the explicit time integrator (2D barotropic path).

        ``scheme`` is one of ``"euler"``, ``"rk2"`` (default), ``"rk3"``,
        ``"rk43"``. Must be called before :meth:`prepare`. Raises if a
        non-``rk2`` scheme is requested on a multilayer simulation (which
        keeps SSP-RK2).
        """
        if scheme not in ("euler", "rk2", "rk3", "rk43"):
            raise ValueError(
                f"time scheme {scheme!r} is not valid; choose one of "
                f"'euler', 'rk2', 'rk3', 'rk43' (semi-implicit is selected "
                f"through the semi-implicit API, not here)")
        name = scheme.encode("utf-8")
        _check(self._lib.rki_unstr_set_time_scheme(
            self._handle,
            ctypes.c_void_p(
                ctypes.cast(ctypes.c_char_p(name), ctypes.c_void_p).value
            ),
            len(name)),
            "rki_unstr_set_time_scheme")

    def set_eos(self, rho_0: float = 1000.0, beta_s: float = 0.78,
                s_ref: float = 0.0, alpha_t: float = 0.17,
                t_ref: float = 15.0):
        """Override the linear EOS knobs used by the multilayer solver.

        ``rho = rho_0 + beta_s*(S - s_ref) - alpha_t*(T - t_ref)``

        Setting both ``beta_s`` and ``alpha_t`` to zero disables the
        baroclinic path (EOS / BPG / vertical exchange are skipped),
        which is the right setting for barotropic-invariance tests
        that must not see sigma PGE noise from a tilted stage.

        Must be called before :meth:`prepare` so the values are
        captured in the device copy of the state struct.
        """
        _check(self._lib.rki_unstr_set_eos(
            self._handle,
            float(rho_0), float(beta_s), float(s_ref),
            float(alpha_t), float(t_ref)),
            "rki_unstr_set_eos")

    def set_cfl(self, cfl: float):
        _check(self._lib.rki_unstr_set_cfl(self._handle, float(cfl)),
               "rki_unstr_set_cfl")

    def set_bpg_method(self, method: str):
        """Select baroclinic pressure gradient method: ``"sigma"`` (default)
        or ``"jacobian"`` (Shchepetkin-McWilliams 2003, no PGE on steep
        bathymetry).

        Must be called before :meth:`prepare` so the flag is captured in
        the device snapshot of the state.
        """
        m = method.lower()
        if m not in ("sigma", "jacobian"):
            raise ValueError(
                f"set_bpg_method: expected 'sigma' or 'jacobian', got {method!r}")
        flag = 1 if m == "jacobian" else 0
        _check(self._lib.rki_unstr_set_bpg_jacobian(self._handle, flag),
               "rki_unstr_set_bpg_jacobian")

    def set_vcoord(self, coord_type: str = "sigma",
                   remap_method: str = "plm",
                   depth_transition: float = 100.0,
                   blend_width: float = 50.0,
                   zstar_h_surf_target: float = 0.0,
                   zstar_h_min: float = 1.0e-4,
                   zstar_n_surf: int = 0,
                   zstar_stretching: str = "log"):
        """Set the vertical coordinate type and remapping method.

        Must be called after :meth:`enable_multilayer` and before
        :meth:`prepare` / :meth:`advance`.

        Parameters
        ----------
        coord_type : str
            ``"sigma"`` (default), ``"zsigma"`` (z-sigma hybrid),
            ``"zstar"`` (z*-lite), ``"zstar_full"`` (full z*
            with per-cell reference profile + vanishing layers), or
            ``"zstar_sigma"`` (z*/sigma smoothstep hybrid: sigma in
            shallow water, z*-lite in deep).
        remap_method : str
            ``"pcm"``, ``"plm"`` (default), ``"ppm"``, ``"ppm_h4"`` or ``"pqm"``.
        depth_transition : float
            For zsigma: depth (m) where sigma-to-z blending starts.
        blend_width : float
            For zsigma: blending zone width (m).
        zstar_h_surf_target : float
            For zstar_full: target physical thickness of the surface
            layer (m).  0 (default) falls back to per-column uniform
            spacing.  Set > 0 to anchor the surface layer at a fixed
            thickness regardless of total depth.
        zstar_h_min : float
            For zstar_full: vanishing-layer floor (m).  Layers that
            would land below the local bed get clipped to this
            thickness.  Default 1e-4 m.
        zstar_n_surf : int
            For zstar_full: number of fine near-surface layers using
            the stretching mode.  0 (default) = auto (max(1, nz/3)).
        zstar_stretching : str
            For zstar_full: ``"log"`` (default, geometric progression
            in the surface fine zone) or ``"uniform"``.
        """
        _COORD_MAP = {"sigma": 1, "zsigma": 2, "zstar": 4,
                      "zstar_full": 5, "zstar_sigma": 6}
        _REMAP_MAP = {"pcm": 1, "plm": 2, "ppm": 3, "ppm_h4": 4, "pqm": 5}
        _STRETCH_MAP = {"uniform": 0, "log": 1}
        if coord_type.lower() not in _COORD_MAP:
            raise ValueError(f"coord_type={coord_type!r} is not valid; "
                             f"choose one of {sorted(_COORD_MAP)}")
        if remap_method.lower() not in _REMAP_MAP:
            raise ValueError(f"remap_method={remap_method!r} is not valid; "
                             f"choose one of {sorted(_REMAP_MAP)}")
        if zstar_stretching.lower() not in _STRETCH_MAP:
            raise ValueError(f"zstar_stretching={zstar_stretching!r} is not "
                             f"valid; choose one of {sorted(_STRETCH_MAP)}")
        ct = _COORD_MAP[coord_type.lower()]
        rm = _REMAP_MAP[remap_method.lower()]
        _check(self._lib.rki_unstr_set_vcoord(
            self._handle, ct, rm, float(depth_transition), float(blend_width)),
            "rki_unstr_set_vcoord")
        if ct == 5:  # VCOORD_ZSTAR_FULL
            stretch = _STRETCH_MAP[zstar_stretching.lower()]
            _check(self._lib.rki_unstr_set_zstar_full_params(
                self._handle, float(zstar_h_surf_target),
                float(zstar_h_min), int(zstar_n_surf), int(stretch)),
                "rki_unstr_set_zstar_full_params")

    def set_betas(self, beta_w: float, beta_uh: float, beta_vh: float):
        """Override the slope-limiter beta values for the wet branch.

        Set all three to ``1.0`` to disable limiting -- recovers full
        second-order convergence on smooth problems but the scheme is
        no longer TVD, so do not use this for shock-containing runs.
        Default ANUGA values are 1.0 / 0.8 / 0.8.
        """
        _check(self._lib.rki_unstr_set_betas(
            self._handle, float(beta_w), float(beta_uh), float(beta_vh)),
            "rki_unstr_set_betas")

    def set_wind(self, tau_x: float, tau_y: float):
        """Set spatially-uniform wind stress (Pa) on the free surface.

        Applied to ``hu`` / ``hv`` as ``+= dt * tau / RHO_WATER`` on
        wet cells. Pass ``0, 0`` to disable. Typical magnitudes for a
        10 m/s wind are around 0.1-0.2 Pa via the standard quadratic
        drag law (the caller is responsible for converting wind speed
        to stress; this entrypoint takes the stress directly so that
        non-quadratic drag laws can be used).
        """
        _check(self._lib.rki_unstr_set_wind(
                   self._handle, float(tau_x), float(tau_y)),
               "rki_unstr_set_wind")

    def set_wind_field(self, wind_x: np.ndarray, wind_y: np.ndarray):
        """Set per-cell wind stress (Pa), x/y components, per owned cell.

        Added to ``hu`` / ``hv`` as ``+= dt * wind / RHO_WATER`` on wet
        cells, on top of any spatially-uniform stress from :meth:`set_wind`.
        Length is the cell count. The caller
        converts wind speed to stress (e.g. a quadratic drag law).
        """
        ax = np.ascontiguousarray(wind_x, dtype=np.float64)
        ay = np.ascontiguousarray(wind_y, dtype=np.float64)
        if ax.size != self._ncells or ay.size != self._ncells:
            raise ValueError(
                f"wind field: expected length {self._ncells}, "
                f"got {ax.size} / {ay.size}"
            )
        px = ax.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        py = ay.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        _check(self._lib.rki_unstr_set_wind_field(self._handle, px, py),
               "rki_unstr_set_wind_field")

    def set_initial_salinity(self, value: float):
        """Set uniform initial salinity (PSU) for all layers.

        Must be called before :meth:`prepare`.
        """
        _check(self._lib.rki_unstr_set_initial_salinity(
                   self._handle, float(value)),
               "rki_unstr_set_initial_salinity")

    def set_initial_temperature(self, value: float):
        """Set uniform initial temperature (degC) for all layers.

        Must be called before :meth:`prepare`.
        """
        _check(self._lib.rki_unstr_set_initial_temperature(
                   self._handle, float(value)),
               "rki_unstr_set_initial_temperature")

    def set_bc_salinity(self, group: int, value: float):
        """Set external salinity (PSU) for a boundary group.

        Used by BC_INFLOW / BC_TIDAL edges: when mass flux enters the
        domain through boundary edges tagged with ``group``, the tracer
        advection kernel uses this value as the inflow salinity.
        """
        _check(self._lib.rki_unstr_set_bc_salinity(
                   self._handle, int(group), float(value)),
               "rki_unstr_set_bc_salinity")

    def set_bc_temperature(self, group: int, value: float):
        """Set external temperature (degC) for a boundary group."""
        _check(self._lib.rki_unstr_set_bc_temperature(
                   self._handle, int(group), float(value)),
               "rki_unstr_set_bc_temperature")

    def set_smagorinsky_cs(self, cs: float):
        """Set Smagorinsky horizontal viscosity coefficient.

        ``nu_h = (Cs * sqrt(A))^2 * |S|`` where A is cell area.
        Typical values: 0.1-0.2. Set to 0 to disable (default).
        """
        _check(self._lib.rki_unstr_set_smagorinsky_cs(self._handle, float(cs)),
               "rki_unstr_set_smagorinsky_cs")

    def set_hdiff_kappa(self, kappa: float):
        """Set horizontal tracer diffusion coefficient (m^2/s).

        Applied as Laplacian diffusion to salinity and temperature.
        Typical values: 1-10 m^2/s. Set to 0 to disable (default).
        """
        _check(self._lib.rki_unstr_set_hdiff_kappa(self._handle, float(kappa)),
               "rki_unstr_set_hdiff_kappa")

    def set_use_ppm_tracer(self, enable: bool = True):
        """Enable/disable high-order tracer advection.

        When True (default), uses MUSCL (Green-Gauss + Barth-Jespersen
        limiter) for second-order tracer advection.  When False, falls
        back to first-order upwind.

        Can be called before or after ``prepare()``.
        """
        _check(self._lib.rki_unstr_set_use_ppm_tracer(
                   self._handle, int(bool(enable))),
               "rki_unstr_set_use_ppm_tracer")

    def set_vmix(self, nu_max: float = 0.01, nu_bg: float = 1e-4,
                 ri_alpha: float = 5.0, kappa_s_bg: float = 1e-5,
                 kappa_t_bg: float = 1e-5):
        """Set Pacanowski-Philander vertical mixing parameters.

        Momentum:  ``nu_t  = nu_max / (1 + alpha*Ri)^2 + nu_bg``
        Tracers:   ``kappa = nu_max / (1 + alpha*Ri)^3 + kappa_bg``

        Parameters
        ----------
        nu_max : float
            Maximum turbulent viscosity (m^2/s). Default 0.01.
        nu_bg : float
            Background vertical viscosity (m^2/s). Default 1e-4.
        ri_alpha : float
            Richardson number scaling. Default 5.0.
        kappa_s_bg : float
            Background salinity diffusivity (m^2/s). Default 1e-5.
        kappa_t_bg : float
            Background temperature diffusivity (m^2/s). Default 1e-5.
        """
        _check(self._lib.rki_unstr_set_vmix(
            self._handle,
            float(nu_max), float(nu_bg), float(ri_alpha),
            float(kappa_s_bg), float(kappa_t_bg)),
            "rki_unstr_set_vmix")

    def set_kpp(self, enable: bool = True, ri_crit: float = 0.3,
                cs_nonlocal: float = 6.3, c_vt2: float = 1.8):
        """Enable/configure the KPP surface-boundary-layer overlay.

        KPP (Large et al. 1994) runs on top of the PP81 interior closure
        (see :meth:`set_vmix`): a bulk-Richardson boundary-layer depth
        sweep, a shape-function-enhanced diffusivity overlaid on PP81,
        and non-local counter-gradient tracer transport. The surface
        friction velocity ``u_*`` comes from the (uniform) wind stress
        set via :meth:`set_wind`; the buoyancy flux from the surface
        heat / virtual-salt fluxes.

        Parameters
        ----------
        enable : bool
            Turn the overlay on (default True). When False the path runs
            pure PP81, bit-identical to pre-KPP.
        ri_crit : float
            Critical bulk Richardson number for the BL-depth sweep.
            Default 0.3.
        cs_nonlocal : float
            Non-local transport coefficient C_s (LMD94 eq 20). Default 6.3.
        c_vt2 : float
            Unresolved-turbulence coefficient for V_t^2 (LMD94 eq 23).
            Default 1.8; set 0 to disable V_t^2.
        """
        _check(self._lib.rki_unstr_set_kpp(
            self._handle,
            1 if enable else 0,
            float(ri_crit), float(cs_nonlocal), float(c_vt2)),
            "rki_unstr_set_kpp")

    # ---------- multilayer (Phase 1b of the 3D unstructured plan) ----------

    def enable_multilayer(self, nz: int):
        """Turn on sigma-layer multilayer mode with ``nz`` layers.

        Call this **immediately after construction** and **before** any
        field setters — the underlying C entrypoint re-initialises the
        state struct, which wipes any bathymetry / stage / BC values
        already loaded. The expected lifecycle is:

        .. code-block:: python

            sim = UnstructuredSimulation(mesh)
            sim.enable_multilayer(nz=5)
            sim.set_bathymetry(bed)
            sim.set_initial_stage(stage)
            sim.advance(1.0)

        In Phase 1 every layer carries a uniform fraction of the
        barotropic momentum so a uniform IC stays shear-free. Future
        phases will populate salinity / temperature ICs through
        additional setters.
        """
        if self._prepared:
            raise RuntimeError(
                "enable_multilayer must be called before prepare()/advance()"
            )
        _check(self._lib.rki_unstr_enable_multilayer(self._handle, int(nz)),
               "rki_unstr_enable_multilayer")

    @property
    def nz_ml(self) -> int:
        """Active number of sigma layers, or 0 when multilayer is off."""
        return int(self._lib.rki_unstr_get_nz_ml(self._handle))

    # ---------- lifecycle ----------

    def prepare(self):
        """Move bathymetry / IC / BC data onto the device and lock in
        the setup.  Must happen before the first :meth:`step` or
        :meth:`advance` — those methods auto-call :meth:`prepare` if
        you haven't, so most users never need to invoke it explicitly.
        Call it by hand only when you want to time setup separately
        from stepping.
        """
        if self._prepared:
            return
        _check(self._lib.rki_unstr_finalize(self._handle), "rki_unstr_finalize")
        self._prepared = True

    def finalize(self):
        """Deprecated alias for :meth:`prepare`.  The old name reads as
        "shut the simulation down" but this call in fact gets the
        simulation ready to run; the rename makes that unambiguous.
        """
        import warnings
        warnings.warn(
            "UnstructuredSimulation.finalize() has been renamed to "
            "prepare(); finalize() will be removed in a future release.",
            DeprecationWarning,
            stacklevel=2,
        )
        self.prepare()

    def step(self, dt: float):
        if not self._prepared:
            self.prepare()
        _check(self._lib.rki_unstr_step(self._handle, float(dt)),
               "rki_unstr_step")

    def advance(self, t_end: float):
        if not self._prepared:
            self.prepare()
        _check(self._lib.rki_unstr_advance(self._handle, float(t_end)),
               "rki_unstr_advance")

    # ---------- gauges ----------

    def init_gauges(self, filename: str, gauge_x, gauge_y):
        """Initialise gauge output at the given (x, y) locations.

        Performs a point-in-triangle search over owned cells to locate
        each gauge. Gauges outside the mesh are silently
        skipped. Creates a CF-1.8 station NetCDF file.

        Parameters
        ----------
        filename : str
            Path to the gauge output file.
        gauge_x, gauge_y : array-like
            Coordinates of gauge stations (metres).
        """
        if not self._prepared:
            self.prepare()
        gx = np.ascontiguousarray(gauge_x, dtype=np.float64)
        gy = np.ascontiguousarray(gauge_y, dtype=np.float64)
        if gx.size != gy.size:
            raise ValueError("gauge_x and gauge_y must have the same length")
        _require_writable_dir(filename)
        fname = filename.encode("utf-8")
        ierr = self._lib.rki_unstr_gauges_init(
            self._handle,
            ctypes.c_void_p(
                ctypes.cast(ctypes.c_char_p(fname), ctypes.c_void_p).value
            ),
            len(fname),
            gx.size,
            gx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            gy.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )
        _check(ierr, "rki_unstr_gauges_init")

    def write_gauges(self):
        """Write gauge values (eta, u, v) at the current simulation time."""
        _check(self._lib.rki_unstr_gauges_write(self._handle),
               "rki_unstr_gauges_write")

    def close_gauges(self):
        """Close the gauge output file."""
        _check(self._lib.rki_unstr_gauges_close(self._handle),
               "rki_unstr_gauges_close")

    # ---------- restart ----------

    def save_restart(self, filename: str):
        """Write a restart checkpoint at the current simulation time.

        Uses lightweight host sync so the GPU data mapping stays intact.
        """
        if not self._prepared:
            self.prepare()
        _require_writable_dir(filename)
        fname = filename.encode("utf-8")
        ierr = self._lib.rki_unstr_save_restart(
            self._handle,
            ctypes.c_void_p(
                ctypes.cast(ctypes.c_char_p(fname), ctypes.c_void_p).value
            ),
            len(fname),
        )
        _check(ierr, "rki_unstr_save_restart")

    def load_restart(self, filename: str):
        """Load state from a restart file.

        Must be called before :meth:`prepare` — the restart populates
        the host-side state arrays (stage, xmom, ymom, bed) and sets
        the simulation time.  Subsequent :meth:`step` / :meth:`advance`
        will auto-prepare and move data onto the device.
        """
        if self._prepared:
            raise RuntimeError("load_restart must be called before prepare")
        _require_file(filename)
        fname = filename.encode("utf-8")
        ierr = self._lib.rki_unstr_load_restart(
            self._handle,
            ctypes.c_void_p(
                ctypes.cast(ctypes.c_char_p(fname), ctypes.c_void_p).value
            ),
            len(fname),
        )
        _check(ierr, "rki_unstr_load_restart")

    # ---------- UGRID NetCDF output ----------

    def init_output(self, filename: str, deflate_level: int = 0):
        """Create a UGRID NetCDF output file and write static topology.

        Call after :meth:`prepare` (or after setting bathymetry at
        minimum). Uses a lightweight host sync so the GPU data mapping
        stays intact — no need to re-prepare after this call.

        Parameters
        ----------
        filename : str
            Path to the output file (e.g. ``"output.nc"``).
        deflate_level : int
            NetCDF deflate compression level (0 = none, 1-9). Default 0.
        """
        if not self._prepared:
            self.prepare()
        _require_writable_dir(filename)
        fname = filename.encode("utf-8")
        ierr = self._lib.rki_unstr_output_init(
            self._handle,
            ctypes.c_void_p(ctypes.cast(ctypes.c_char_p(fname),
                                        ctypes.c_void_p).value),
            len(fname),
            int(deflate_level),
        )
        _check(ierr, "rki_unstr_output_init")

    def write_output(self):
        """Append the current (stage, xmom, ymom) to the output file.

        Uses lightweight host sync — the GPU data mapping stays intact
        and stepping can continue without re-finalizing.
        """
        _check(self._lib.rki_unstr_output_write(self._handle),
               "rki_unstr_output_write")

    def close_output(self):
        """Close the UGRID output file."""
        _check(self._lib.rki_unstr_output_close(self._handle),
               "rki_unstr_output_close")

    # ---------- field views ----------
    #
    # Each accessor calls the C getter, which triggers a host sync if
    # data is on the device. Use sparingly mid-run; for batch
    # post-processing this is the right pattern.

    def _view(self, sym: str) -> np.ndarray:
        ptr = ctypes.c_void_p()
        n = ctypes.c_int(0)
        ierr = getattr(self._lib, sym)(self._handle,
                                        ctypes.byref(ptr), ctypes.byref(n))
        _check(ierr, sym)
        # The C getter calls refresh_host() which uses `!$acc update
        # self` to sync device->host while keeping the GPU mapping
        # intact, so we leave `_prepared` alone. The numpy view does not
        # own the memory but keeps this wrapper alive (see _owned_view).
        return _owned_view(ptr.value, self._lib._rki_wp_ctype, n.value, self)

    def _view_layer(self, sym: str) -> np.ndarray:
        """Return a ``(ncells_owned, nz_ml)`` numpy view of a layer field.

        The Fortran storage is ``(nz_ml, ncells_total)`` with layer as
        the fast (leading) axis. Numpy is row-major by default, so we
        reshape the flat buffer as ``(ncells_owned, nz_ml)`` — reading
        ``arr[cell, layer]`` matches the Fortran ``arr(layer, cell)``
        lookup without a transpose.
        """
        ptr = ctypes.c_void_p()
        nz = ctypes.c_int(0)
        ncells = ctypes.c_int(0)
        ierr = getattr(self._lib, sym)(
            self._handle,
            ctypes.byref(ptr), ctypes.byref(nz), ctypes.byref(ncells)
        )
        _check(ierr, sym)
        total = int(nz.value) * int(ncells.value)
        arr = _owned_view(ptr.value, self._lib._rki_wp_ctype, total, self)
        return arr.reshape((int(ncells.value), int(nz.value)))

    @property
    def stage(self) -> np.ndarray:
        return self._view("rki_unstr_get_stage_ptr")

    @property
    def xmom(self) -> np.ndarray:
        return self._view("rki_unstr_get_xmom_ptr")

    @property
    def ymom(self) -> np.ndarray:
        return self._view("rki_unstr_get_ymom_ptr")

    @property
    def bed(self) -> np.ndarray:
        return self._view("rki_unstr_get_bed_ptr")

    @property
    def height(self) -> np.ndarray:
        return self._view("rki_unstr_get_height_ptr")

    @property
    def h_layer(self) -> np.ndarray:
        """Layer thicknesses, shape ``(ncells_owned, nz_ml)``.

        Only valid when :meth:`enable_multilayer` has been called; the
        C getter errors out otherwise.
        """
        return self._view_layer("rki_unstr_get_h_layer_ptr")

    @property
    def hu_layer(self) -> np.ndarray:
        """Layer x-momentum ``h_k * u_k``, shape ``(ncells_owned, nz_ml)``."""
        return self._view_layer("rki_unstr_get_hu_layer_ptr")

    @property
    def hv_layer(self) -> np.ndarray:
        """Layer y-momentum ``h_k * v_k``, shape ``(ncells_owned, nz_ml)``."""
        return self._view_layer("rki_unstr_get_hv_layer_ptr")

    @property
    def hS_layer(self) -> np.ndarray:
        """Thickness-weighted salinity ``h_k * S_k``, shape ``(ncells_owned, nz_ml)``.

        Divide by :attr:`h_layer` to recover salinity in PSU.
        """
        return self._view_layer("rki_unstr_get_hS_layer_ptr")

    @property
    def hT_layer(self) -> np.ndarray:
        """Thickness-weighted temperature ``h_k * T_k``, shape ``(ncells_owned, nz_ml)``."""
        return self._view_layer("rki_unstr_get_hT_layer_ptr")

    @property
    def salinity(self) -> np.ndarray:
        """Per-layer salinity (PSU), shape ``(ncells_owned, nz_ml)``.

        Computed as ``hS_layer / h_layer``. Dry cells (h_k ~ 0) return 0.
        """
        hS = self.hS_layer
        h = self.h_layer
        with np.errstate(divide="ignore", invalid="ignore"):
            S = np.where(h > 1e-10, hS / h, 0.0)
        return S

    @property
    def temperature(self) -> np.ndarray:
        """Per-layer temperature (degC), shape ``(ncells_owned, nz_ml)``."""
        hT = self.hT_layer
        h = self.h_layer
        with np.errstate(divide="ignore", invalid="ignore"):
            T = np.where(h > 1e-10, hT / h, 0.0)
        return T

    @property
    def centroid_x(self) -> np.ndarray:
        return self._view("rki_unstr_get_centroid_x_ptr")

    @property
    def centroid_y(self) -> np.ndarray:
        return self._view("rki_unstr_get_centroid_y_ptr")

    @property
    def area(self) -> np.ndarray:
        """Per-cell triangle area (m^2). Static -- safe to cache."""
        return self._view("rki_unstr_get_area_ptr")

    def volume(self) -> float:
        """Total water volume in the domain (m^3).

        Computed as ``sum((stage - bed) * area)`` -- correct from the
        moment the IC is set, regardless of whether ``height_cv`` has
        been populated by a kernel pass yet. Use this for mass
        conservation diagnostics in tests.
        """
        h = np.maximum(self.stage - self.bed, 0.0)
        return float((h * self.area).sum())
