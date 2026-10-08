#!/usr/bin/env python3
"""Live plotting example: salt-wedge estuary breathing with M2 tides.

Idealised funnel-shaped estuary with:

* **M2 tidal elevation boundary** on the east (ocean) end
* **Steady freshwater river discharge** at the west (head)
* **Stratified initial condition** — fresher on top, salty on bottom

The money-shot panel is a vertical cross-section along the estuary
thalweg (y = mid-channel) that shows the **salt wedge** advancing
up-estuary on flood tide and retreating on ebb — the classic
partially-mixed estuary behaviour.  The dense ocean water intrudes
along the bottom while the buoyant river plume flows out along the
top, so over a tidal cycle the isohalines tilt and wobble back and
forth in the x-z plane.

Three panels, synchronised:

1. **Top-left**: surface salinity map (the river plume signal)
2. **Top-right**: free-surface η map (the tidal wave propagating
   up-estuary, amplifying as the funnel narrows)
3. **Bottom**: x-z salinity cross-section along the thalweg — the
   salt wedge itself

Inputs are loaded from ``live_tidal_estuary.yaml`` next to this script.

Run::

    python live_tidal_estuary.py
"""
import os
import sys
import math
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

# ---- Physics ----
M2_AMP = cfg.forcing.tide.amplitude
M2_PERIOD = cfg.forcing.tide.period
M2_OMEGA = 2.0 * math.pi / M2_PERIOD
RIVER_V = cfg.forcing.inflow_velocity
RIVER_S = cfg.forcing.inflow_salinity
BG_S = cfg.salinity.initial
T_END = cfg.domain.t_end
N_CYCLES = T_END / M2_PERIOD

# ---- Build ----
sim = Rakali(**to_rakali_kwargs(cfg))

# ---- Funnel bathymetry ----
# Narrow & shallow at the head, wide & deep at the mouth.  Carve
# land strips on the north and south as the estuary widens.
DEPTH_HEAD = cfg.bathymetry.funnel.depth_head
DEPTH_MOUTH = cfg.bathymetry.funnel.depth_mouth
WIDTH_HEAD_CELLS = cfg.bathymetry.funnel.width_head_cells
WIDTH_MOUTH_CELLS = cfg.bathymetry.funnel.width_mouth_cells
LAND_HEIGHT = cfg.bathymetry.land_height

bathy = np.full((NX, NY), LAND_HEIGHT, dtype=np.float64, order='F')
jc = NY / 2.0
for i in range(NX):
    frac = i / (NX - 1)
    half = WIDTH_HEAD_CELLS + frac * (WIDTH_MOUTH_CELLS - WIDTH_HEAD_CELLS)
    depth = DEPTH_HEAD + frac * (DEPTH_MOUTH - DEPTH_HEAD)
    j_lo = int(round(jc - half))
    j_hi = int(round(jc + half))
    bathy[i, j_lo:j_hi] = -depth

sim.set_bathymetry(bathy)
h_init = np.maximum(-bathy, 0.0).astype(np.float64)
sim.set_initial_h(np.asfortranarray(h_init))
sim.finalize()

# M2 tidal constituent on the east (ocean) boundary.
sim.set_tidal_constituents(
    amp=np.array([M2_AMP]),
    phase=np.array([0.0]),
    omega=np.array([M2_OMEGA]),
)

# Initial stratification: linear fresh→salty from surface (k=0) to
# bottom (k=NZ-1), *and* a pre-existing salinity gradient along the
# estuary so the wedge starts somewhere sensible instead of having
# to advect in over many cycles.  Surface is 2 PSU fresher than
# bottom everywhere; along-channel we ramp from RIVER_S at the head
# to BG_S at the mouth.
S_init = np.empty((NX, NY, NZ), dtype=np.float64)
for i in range(NX):
    frac = i / (NX - 1)
    S_column_mean = RIVER_S + frac * (BG_S - RIVER_S)
    for k in range(NZ):
        # surface k=0 fresher, bottom k=NZ-1 saltier
        dS = -2.0 * (1.0 - 2.0 * k / max(NZ - 1, 1))
        S_init[i, :, k] = S_column_mean + dS
S_init = np.clip(S_init, RIVER_S, BG_S)
sim.set_salinity(np.asfortranarray(S_init))

print(f"Created: {sim}")
print(f"  Estuary: {LX_KM:.1f} km x {LY_KM:.1f} km, depth {DEPTH_HEAD}-{DEPTH_MOUTH} m")
print(f"  M2 tide: amp = {M2_AMP} m, period = {M2_PERIOD/3600:.2f} h")
print(f"  River: v = {RIVER_V} m/s at head, S = {RIVER_S} PSU")
print(f"  Run: {N_CYCLES:.1f} M2 cycles ({T_END/3600:.1f} h)")

# ---- Figure: 3-panel layout ----
plt.ion()
fig = plt.figure(figsize=(13, 9))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.2], hspace=0.28, wspace=0.22)
ax_surf = fig.add_subplot(gs[0, 0])
ax_eta = fig.add_subplot(gs[0, 1])
ax_xz = fig.add_subplot(gs[1, :])

extent = [0, LX_KM, 0, LY_KM]

import matplotlib as mpl
cmap_s = mpl.colormaps['viridis_r'].copy(); cmap_s.set_bad('#aaaaaa')
cmap_e = mpl.colormaps['RdBu_r'].copy();    cmap_e.set_bad('#aaaaaa')

J_MID = NY // 2
b_mid = sim.b_phys[:, J_MID]          # negative where wet
# Static land mask: land is defined purely by the bathymetry (b > 0).
# Using `h_phys < 1e-6` as the mask would be dynamic and fail the
# instant a BC seeps a film of water onto a nominally-dry cell,
# making the grey strips flicker into view mid-simulation.
LAND_MASK = sim.b_phys > 0.0


def masked_surface_S(s):
    # Surface = layer nz-1 under k=1=bed, k=nz=surface bottom-up.
    S = s.salinity[:, :, -1]
    return np.ma.masked_array(np.where(np.isnan(S), 0.0, S),
                              mask=LAND_MASK)


def masked_eta(s):
    eta = s.h_phys + s.b_phys
    return np.ma.masked_array(eta, mask=LAND_MASK)


def thalweg_cross_section(s):
    """Return (X, Z, S) for a pcolormesh along y = J_MID.

    Z interfaces per column are rebuilt from cumulative layer
    thickness starting at the free surface.
    """
    # Reverse the layer axis so [0]=surface (matches the build-from-
    # top z_interfaces below). Native k=1=bed, k=nz=surface (bottom-up).
    h_lay = s.h_layer[:, J_MID, ::-1]              # (nx, nz), [0]=surface
    salt = s.salinity[:, J_MID, ::-1]              # (nx, nz), [0]=surface
    eta_1d = s.h_phys[:, J_MID] + s.b_phys[:, J_MID]
    z_top = eta_1d
    z_interfaces = np.empty((NX, NZ + 1))
    z_interfaces[:, 0] = z_top
    z_interfaces[:, 1:] = z_top[:, None] - np.cumsum(h_lay, axis=1)
    x_edges = np.arange(NX + 1) * DX / 1000.0
    X = np.repeat(x_edges[:, None], NZ + 1, axis=1)
    Z = np.empty_like(X)
    Z[:-1, :] = z_interfaces
    Z[-1, :] = z_interfaces[-1, :]
    return X, Z, salt


# -- surface salinity --
img_s = ax_surf.imshow(
    masked_surface_S(sim).T, origin='lower', extent=extent,
    cmap=cmap_s, vmin=RIVER_S, vmax=BG_S, interpolation='nearest',
)
fig.colorbar(img_s, ax=ax_surf, label='S (PSU)')
ax_surf.set_title("surface salinity")
ax_surf.set_xlabel('x (km)'); ax_surf.set_ylabel('y (km)')
ax_surf.set_aspect('equal')

# -- η --
img_e = ax_eta.imshow(
    masked_eta(sim).T, origin='lower', extent=extent,
    cmap=cmap_e, vmin=-M2_AMP, vmax=M2_AMP, interpolation='nearest',
)
fig.colorbar(img_e, ax=ax_eta, label='η (m)')
ax_eta.set_title("free surface η")
ax_eta.set_xlabel('x (km)'); ax_eta.set_ylabel('y (km)')
ax_eta.set_aspect('equal')

# -- x-z salinity along the thalweg --
X0, Z0, S0 = thalweg_cross_section(sim)
mesh_xz = ax_xz.pcolormesh(
    X0, Z0, S0, cmap=cmap_s, vmin=RIVER_S, vmax=BG_S, shading='flat',
)
fig.colorbar(mesh_xz, ax=ax_xz, label='S (PSU)')
# Bathymetry outline on the cross-section
xs_km = (np.arange(NX) + 0.5) * DX / 1000.0
(bline,) = ax_xz.plot(xs_km, b_mid, color='#222', lw=1.5, label='bed')
# 5/20/30 PSU isohalines are a useful "salt wedge" indicator
ax_xz.set_xlabel('along-estuary x (km)')
ax_xz.set_ylabel('z (m)')
ax_xz.set_title("salinity along the thalweg (y = mid-channel)  "
                "— watch the salt wedge breathe with the tide")
ax_xz.set_xlim(0, LX_KM)
ax_xz.set_ylim(-DEPTH_MOUTH - 1.0, 3.0)
ax_xz.axhline(0.0, color='k', lw=0.4, alpha=0.5)

suptitle = fig.suptitle("t = 0.00 h   (0.0 M2 cycles)",
                        fontsize=13, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.97))
plt.show()


def redraw(s):
    img_s.set_data(masked_surface_S(s).T)
    img_e.set_data(masked_eta(s).T)
    _, _, salt = thalweg_cross_section(s)
    mesh_xz.set_array(salt.ravel())
    t_h = s.time / 3600.0
    n_cyc = s.time / M2_PERIOD
    # M2 phase indicator: sin(ωt) → +1 at flood peak, -1 at ebb peak
    phase = math.sin(M2_OMEGA * s.time)
    phase_label = ("flood" if phase > 0.3 else
                   "ebb" if phase < -0.3 else
                   "slack")
    suptitle.set_text(
        f"t = {t_h:5.2f} h   ({n_cyc:.2f} M2 cycles)   "
        f"tide: {phase_label}"
    )
    fig.canvas.draw_idle()
    fig.canvas.flush_events()


print("Animating salt-wedge estuary...")
sim.advance(T_END, callback=redraw, callback_interval=300.0, progress=True)
print("Done.  Close the figure window to exit.")

plt.ioff()
plt.show()
sim.destroy()
