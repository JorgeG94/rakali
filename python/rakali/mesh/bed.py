"""Bathymetry loading and interpolation helpers.

Works with the minimal Rakali bathymetry schema:

* NetCDF file with dimensions ``x`` (nx) and ``y`` (ny).
* 2-D variable named ``b``, ``elevation``, or ``depth`` (positive-up,
  so ``b < 0`` is wet).
* Cell spacing ``dx`` is not stored — it must be supplied by the
  caller.  Cell centres are assumed at ``(i*dx, j*dx)``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import numpy as np


_BATHY_VAR_CANDIDATES = ("b", "elevation", "depth")


def load_bathymetry_nc(
    filename: str | Path,
    varname: Optional[str] = None,
) -> np.ndarray:
    """Load the 2-D bathymetry variable from ``filename`` as an
    ``(nx, ny)`` ``float64`` array.

    Tries ``varname`` first if given; otherwise walks the default list
    ``("b", "elevation", "depth")``.  Matches the variable lookup
    convention in ``src/io/rki_bathymetry.F90``.
    """
    try:
        from netCDF4 import Dataset
    except ImportError as exc:
        raise ImportError(
            "load_bathymetry_nc requires netCDF4. "
            "Install with `conda install -c conda-forge netcdf4`."
        ) from exc

    candidates = (varname,) if varname else _BATHY_VAR_CANDIDATES

    ds = Dataset(str(filename), "r")
    try:
        for name in candidates:
            if name in ds.variables:
                arr = np.array(ds.variables[name][:], dtype=np.float64)
                break
        else:
            raise KeyError(
                f"no bathymetry variable in {filename}. "
                f"Tried: {candidates}"
            )
    finally:
        ds.close()

    if arr.ndim != 2:
        raise ValueError(
            f"bathymetry variable must be 2-D; got shape {arr.shape}"
        )
    return arr


def window_bathymetry(
    bathy: np.ndarray,
    dx: float,
    window: Optional[Tuple[float, float, float, float]] = None,
    dy: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract a rectangular sub-region of a raster bathymetry.

    Parameters
    ----------
    bathy : ``(nx, ny)`` float array
    dx : float
        Cell size in the x direction (metres).
    window : ``(x0, x1, y0, y1)`` or None
        Sub-region in metres.  ``None`` returns the full raster.
        Coordinates are clipped to the raster extent.
    dy : float, optional
        Cell size in the y direction.  Defaults to ``dx`` (square
        cells), which matches every shipped example.

    Returns
    -------
    sub : ``(nx_sub, ny_sub)`` array
    x, y : 1-D coordinate arrays for ``sub``
    """
    nx, ny = bathy.shape
    if dy is None:
        dy = dx

    if window is None:
        ix0, ix1 = 0, nx
        iy0, iy1 = 0, ny
    else:
        x0, x1, y0, y1 = window
        ix0 = int(np.clip(round(x0 / dx), 0, nx))
        ix1 = int(np.clip(round(x1 / dx), 0, nx))
        iy0 = int(np.clip(round(y0 / dy), 0, ny))
        iy1 = int(np.clip(round(y1 / dy), 0, ny))
        if ix1 <= ix0 or iy1 <= iy0:
            raise ValueError(
                f"window {window} is empty or outside the raster extent "
                f"(nx={nx}, ny={ny}, dx={dx}, dy={dy})."
            )

    sub = bathy[ix0:ix1, iy0:iy1].copy()
    x = np.arange(ix0, ix1) * dx
    y = np.arange(iy0, iy1) * dy
    return sub, x, y


def interpolate_bed(
    bathy: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    cx: np.ndarray,
    cy: np.ndarray,
    method: str = "nearest",
    fill_value: float = 0.0,
) -> np.ndarray:
    """Sample a raster bathymetry at ``(cx, cy)`` points.

    Uses ``scipy.interpolate.RegularGridInterpolator``.  ``"nearest"``
    matches the hand-written examples and is a conservative choice for
    coastal bathymetry with sharp step features; ``"linear"`` smooths
    but can introduce spurious shallow regions across steep grid cells.
    """
    try:
        from scipy.interpolate import RegularGridInterpolator
    except ImportError as exc:
        raise ImportError(
            "interpolate_bed requires scipy. "
            "Install with `conda install -c conda-forge scipy`."
        ) from exc

    interp = RegularGridInterpolator(
        (x, y), bathy,
        method=method, bounds_error=False, fill_value=float(fill_value),
    )
    cx = np.ascontiguousarray(cx, dtype=np.float64)
    cy = np.ascontiguousarray(cy, dtype=np.float64)
    return interp(np.column_stack([cx, cy]))
