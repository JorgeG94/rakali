#!/usr/bin/env python3
"""Live plotting example: 4-layer river plume on a coastal shelf.

A more physically realistic plume than a simple uniform west-edge
inflow:

* The west edge is a coastline (above sea level) **except** for a
  narrow river mouth ~1 km wide, positioned north of centre so the
  plume has room to develop.
* Bathymetry is uniform across the basin (sigma-layer artefacts make
  a sloping shelf misleading — see comments below).
* Coriolis is on, so the plume curves to the right (i.e. southwards
  in this domain) as it leaves the river mouth — the classic
  "anticyclonic bulge" of buoyancy-driven coastal currents.
* Smagorinsky horizontal viscosity and a small horizontal salinity
  diffusion smooth the front without erasing it.
* The east boundary is a Flather open boundary so the plume can
  radiate out instead of bouncing back.

The animation shows two panels: surface salinity (the plume head)
and the depth-averaged salinity (the depth-integrated freshness
signal — picks up the plume even where it has thinned vertically).
The coastline contour (b = 0) is overlaid on both panels in white.

Inputs are loaded from ``live_river_plume.yaml`` next to this script.

Run with::

    python live_river_plume.py
"""
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "python"))

from rakali import Rakali
from rakali.yaml_loader import load_config, to_rakali_kwargs

cfg = load_config(Path(__file__).with_suffix(".yaml"))

# ---- Domain ----
NX, NY = cfg.domain.nx, cfg.domain.ny
DX, DY = cfg.domain.dx, cfg.domain.dy
NZ = cfg.domain.nz
LX_KM = NX * DX / 1000.0
LY_KM = NY * DY / 1000.0

# ---- Coastal geometry ----
# Two important notes after debugging "no plume visible":
#
# 1. Use a UNIFORM basin depth.  With sigma layers, the surface layer
#    thickness is H/NZ — so a sloping shelf (4 m at the coast → 25 m
#    offshore) makes the surface layer 4-6x thicker offshore, which
#    *automatically* dilutes any incoming freshwater by the same
#    factor.  Plume signal disappears into the mean offshore S.
#
# 2. Match the river mouth depth to the basin depth.  Same reason:
#    if the channel is 4 m and the offshore is 15 m, the leading
#    edge of the plume crosses an instantaneous 4× layer-thickness
#    increase and ~73 % of its freshness is lost to the deeper
#    surface layer offshore.
#
# Net result: keep the geometry simple, and let buoyancy + Coriolis
# do the visual work.
BASIN_DEPTH = cfg.bathymetry.depth
LAND_HEIGHT = cfg.bathymetry.land_height
COAST_WIDTH = cfg.bathymetry.coastline.west_strip_cells
RIVER_MOUTH_J = int(NY * cfg.bathymetry.coastline.river_mouth_j_fraction)
RIVER_MOUTH_HALFWIDTH = cfg.bathymetry.coastline.river_mouth_halfwidth_cells

# ---- Forcing / physics ----
INFLOW_V = cfg.forcing.inflow_velocity
INFLOW_S = cfg.forcing.inflow_salinity
BG_S = cfg.salinity.initial
CORIOLIS_F = cfg.forcing.coriolis_f
T_END = cfg.domain.t_end

# ---- Build the Rakali simulation ----
sim = Rakali(**to_rakali_kwargs(cfg))

# ---- Bathymetry: flat basin with a coastline and a narrow mouth ----
# Sign convention: b > 0 is above sea level (dry land), b < 0 is wet.
# Uniform BASIN_DEPTH everywhere except for the western coastline strip.
bathy = np.full((NX, NY), -BASIN_DEPTH, dtype=np.float64, order='F')

# Carve a strip of "land" along the west edge.
bathy[:COAST_WIDTH, :] = LAND_HEIGHT

# Cut a narrow river mouth through the land at j = RIVER_MOUTH_J.
# Keep the mouth at full BASIN_DEPTH so there's no thickness jump.
mouth_lo = RIVER_MOUTH_J - RIVER_MOUTH_HALFWIDTH
mouth_hi = RIVER_MOUTH_J + RIVER_MOUTH_HALFWIDTH
bathy[:COAST_WIDTH, mouth_lo:mouth_hi] = -BASIN_DEPTH

sim.set_bathymetry(bathy)

# Initial water depth: still water at eta=0.  h = max(-b, 0) gives
# zero on land cells (dry) and the local depth elsewhere.
h_init = np.maximum(-bathy, 0.0).astype(np.float64)
h_init = np.asfortranarray(h_init)
sim.set_initial_h(h_init)

sim.finalize()
print(f"Created: {sim}")
print(f"  Domain: {LX_KM:.1f} km x {LY_KM:.1f} km, {NZ} layers, "
      f"uniform depth {BASIN_DEPTH:.0f} m")
print(f"  River mouth: ~{2 * RIVER_MOUTH_HALFWIDTH * DX / 1000:.1f} km wide "
      f"at y = {RIVER_MOUTH_J * DY / 1000:.1f} km")
print(f"  Inflow: v = {INFLOW_V} m/s, S = {INFLOW_S} PSU "
      f"(background {BG_S} PSU)")
print(f"  Coriolis f = {CORIOLIS_F:.2e} rad/s; "
      f"inertial radius U/f = {INFLOW_V / CORIOLIS_F / 1000:.2f} km")

# ---- Mixing diagnostics computed from the salinity field ----
#
# All three are pure numpy ops on sim.salinity / sim.h_layer — no
# new C API needed.  See log_book entries for the full reasoning.
#
#   contrast(x, y)  = S[k=0] - S[k=nz-1]
#       Signed surface-to-bottom salinity difference.  Negative
#       (fresher on top) for a buoyant plume — that's the normal
#       state.  Magnitude is the strength of the stratification.
#
#   sigma(x, y)     = std(S, axis=z)
#       Column standard deviation of salinity.  Zero in a
#       perfectly mixed water column, large where stratification
#       is sharp.  Single non-negative scalar per cell.
#
#   pycno_depth(x, y) = depth (m) of the layer interface with the
#       largest |dS/dz|, i.e. the pycnocline location.  Computed
#       as cumsum(h_layer) at the argmax interface.  With NZ small
#       this is necessarily coarse (only NZ-1 distinct values), so
#       it works best as a "where has the front sunk to" indicator.

def compute_diagnostics(s):
    """Return (surface_S, contrast, sigma, pycno_depth) as masked arrays.

    Dry land cells (where h_phys ≈ 0) are masked out so they don't
    contaminate min/max diagnostics or fake a "freshwater" signal in
    the plot via the salinity accessor's `S = 0 where dry` convention.
    """
    S = s.salinity                          # (nx, ny, nz)
    h_lay = s.h_layer                       # (nx, ny, nz)
    dry = s.h_phys < 1.0e-6                 # (nx, ny) — True on land

    surface = S[:, :, 0]
    contrast = S[:, :, 0] - S[:, :, -1]
    sigma = S.std(axis=2)

    if S.shape[2] >= 2:
        dS_dz = np.abs(np.diff(S, axis=2))  # (nx, ny, nz-1)
        pycno_k = np.argmax(dS_dz, axis=2)  # (nx, ny)
        # Depth of the bottom of layer k = depth of interface k
        depths = np.cumsum(h_lay, axis=2)
        pycno_depth = np.take_along_axis(
            depths, pycno_k[:, :, np.newaxis], axis=2
        ).squeeze(-1)
    else:
        pycno_depth = np.zeros_like(contrast)

    return (
        np.ma.masked_array(surface,     mask=dry),
        np.ma.masked_array(contrast,    mask=dry),
        np.ma.masked_array(sigma,       mask=dry),
        np.ma.masked_array(pycno_depth, mask=dry),
    )


# ---- Set up the live figure ----
plt.ion()
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
ax_surf, ax_contrast = axes[0]
ax_sigma, ax_pycno = axes[1]

extent = [0, LX_KM, 0, LY_KM]

# Bathymetry contour for the coastline (b = 0, i.e. sea level)
b_phys = sim.b_phys
xs = (np.arange(NX) + 0.5) * DX / 1000.0
ys = (np.arange(NY) + 0.5) * DY / 1000.0

# Initial state
surf0, contrast0, sigma0, pycno0 = compute_diagnostics(sim)
contrast_lim = BG_S - INFLOW_S          # symmetric range for the diverging cmap
sigma_lim = 0.5 * (BG_S - INFLOW_S)
# Pycnocline depth scale: cap at the basin depth
pycno_max = float(BASIN_DEPTH)


def _make_cmap(name):
    """Return a copy of a matplotlib colormap with land set to light grey."""
    import matplotlib as mpl
    cmap = mpl.colormaps[name].copy()
    cmap.set_bad('#cccccc')
    return cmap


def _decorate(ax, title, label, img):
    ax.contour(xs, ys, b_phys.T, levels=[0.0], colors='black',
               linewidths=0.8)
    ax.set_title(title)
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_aspect('equal')
    fig.colorbar(img, ax=ax, label=label)


# (0, 0) — Surface salinity (the plume head)
img_surf = ax_surf.imshow(
    surf0.T,
    origin='lower', extent=extent, cmap=_make_cmap('viridis_r'),
    vmin=INFLOW_S, vmax=BG_S, interpolation='nearest',
)
_decorate(ax_surf, "Surface salinity (k = 0)", "S (PSU)", img_surf)

# (0, 1) — Surface minus bottom salinity (signed stratification)
img_contrast = ax_contrast.imshow(
    contrast0.T,
    origin='lower', extent=extent, cmap=_make_cmap('RdBu'),
    vmin=-contrast_lim, vmax=+contrast_lim, interpolation='nearest',
)
_decorate(
    ax_contrast,
    "Surface − bottom salinity  (negative = stable buoyant plume)",
    "ΔS (PSU)", img_contrast,
)

# (1, 0) — Column standard deviation (mixing intensity)
img_sigma = ax_sigma.imshow(
    sigma0.T,
    origin='lower', extent=extent, cmap=_make_cmap('magma'),
    vmin=0.0, vmax=sigma_lim, interpolation='nearest',
)
_decorate(
    ax_sigma,
    "Column std of salinity  (large = strongly stratified)",
    "σ_S (PSU)", img_sigma,
)

# (1, 1) — Pycnocline depth (interface of strongest |dS/dz|)
img_pycno = ax_pycno.imshow(
    pycno0.T,
    origin='lower', extent=extent, cmap=_make_cmap('cividis_r'),
    vmin=0.0, vmax=pycno_max, interpolation='nearest',
)
_decorate(
    ax_pycno,
    f"Pycnocline depth  (NZ={NZ} → {NZ - 1} possible interfaces)",
    "depth (m)", img_pycno,
)

suptitle = fig.suptitle(f"t = 0.0 min")
fig.tight_layout()
plt.show()


def redraw(s):
    surf, contrast, sigma, pycno = compute_diagnostics(s)
    img_surf.set_data(surf.T)
    img_contrast.set_data(contrast.T)
    img_sigma.set_data(sigma.T)
    img_pycno.set_data(pycno.T)
    # min / max are taken over the masked array → wet cells only,
    # ignoring the dry land strip on the west.
    s_min = surf.min()
    s_max = surf.max()
    dS_max = np.abs(contrast).max()
    suptitle.set_text(
        f"t = {s.time / 60.0:6.1f} min   "
        f"wet surf S ∈ [{s_min:.2f}, {s_max:.2f}] PSU   "
        f"max |ΔS| = {dS_max:.2f} PSU"
    )
    fig.canvas.draw_idle()
    fig.canvas.flush_events()
    # Console debug — useful when the figure looks wrong
    print(f"  t = {s.time / 60.0:6.1f} min   "
          f"wet surf S min = {s_min:.3f} PSU   "
          f"max |ΔS| = {dS_max:.3f} PSU")


print("Animating river plume (2 hours of sim time)...")
sim.advance(
    T_END,
    callback=redraw,
    callback_interval=600.0,   # update every 60 sim-seconds (≈120 frames)
    progress=True,            # tqdm bar if available
)
print("Done.  Close the figure window to exit.")

plt.ioff()
plt.show()

sim.destroy()
