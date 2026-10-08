#!/usr/bin/env python3
"""Live plotting example: vorticity in a pure 2D barotropic shallow-water sim.

Classic "flow past an island" setup.  The 2D shallow-water equations
carry a perfectly good vorticity equation,

    D(ζ + f)/Dt = -(ζ + f) ∇·u_H ,

so even though there is no stratification and no third dimension, a
barotropic solver can generate, advect, and stretch vorticity.  The
cleanest demo is to drive uniform flow past a no-slip obstacle: the
boundary layer separates from the sides of the island, rolls up, and
sheds a Kármán vortex street downstream.  Vorticity starts at
identically zero and is generated entirely by the interaction between
the pressure field and the coast.

Diagnostic: vorticity ζ = ∂v/∂x - ∂u/∂y computed on cell centres from
the Rakali state via u = hu/h, v = hv/h.  Displayed with a diverging
colormap so positive (counter-clockwise) and negative (clockwise)
vortices show up as red and blue.

Inputs are loaded from ``live_vorticity_2d.yaml`` next to this script.

Run with::

    python live_vorticity_2d.py
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

NX, NY = cfg.domain.nx, cfg.domain.ny
DX, DY = cfg.domain.dx, cfg.domain.dy
T_END = cfg.domain.t_end
LX_KM = NX * DX / 1000.0
LY_KM = NY * DY / 1000.0

BASIN_DEPTH = cfg.bathymetry.depth
LAND_HEIGHT = cfg.bathymetry.land_height
ISLAND_CX = int(NX * cfg.bathymetry.island.cx_fraction)
ISLAND_CY = int(NY * cfg.bathymetry.island.cy_fraction)
ISLAND_R = cfg.bathymetry.island.radius_cells

U_INFLOW = cfg.forcing.inflow_velocity

sim = Rakali(**to_rakali_kwargs(cfg))

# ---- Bathymetry: flat channel with a circular island ----
# Sign convention: b > 0 above sea level (land), b < 0 wet.
bathy = np.full((NX, NY), -BASIN_DEPTH, dtype=np.float64, order='F')

ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
island_mask = (ii - ISLAND_CX) ** 2 + (jj - ISLAND_CY) ** 2 <= ISLAND_R ** 2
bathy[island_mask] = LAND_HEIGHT

sim.set_bathymetry(bathy)

# Initial depth: still water at eta = 0; dry on the island.
h_init = np.maximum(-bathy, 0.0).astype(np.float64)
h_init = np.asfortranarray(h_init)
sim.set_initial_h(h_init)

sim.finalize()
print(f"Created: {sim}")
print(f"  Domain: {LX_KM:.1f} km x {LY_KM:.1f} km, flat depth {BASIN_DEPTH:.0f} m")
print(f"  Island: radius {ISLAND_R * DX:.0f} m at "
      f"({ISLAND_CX * DX / 1000:.2f} km, {ISLAND_CY * DY / 1000:.2f} km)")
print(f"  Inflow: U = {U_INFLOW:.2f} m/s  "
      f"(advective crossing time ≈ {LX_KM * 1000 / U_INFLOW / 60:.1f} min)")

# ---- Vorticity diagnostic ----
#
# 2D barotropic vorticity ζ = ∂v/∂x - ∂u/∂y, computed on cell centres.
# u = hu / h, v = hv / h wherever h is wet; land cells are masked so they
# don't contaminate the finite-difference stencil or the colour range.
#
# np.gradient uses centred differences in the interior and one-sided at
# the edges, which is more than enough resolution for a visualisation.


def compute_vorticity(s):
    h = s.h_phys
    hu = s.hu_phys
    hv = s.hv_phys
    wet = h > 1.0e-6

    u = np.zeros_like(h)
    v = np.zeros_like(h)
    np.divide(hu, h, out=u, where=wet)
    np.divide(hv, h, out=v, where=wet)

    # Zero velocity on dry cells so the no-slip boundary shows up cleanly
    # in the gradient.
    u[~wet] = 0.0
    v[~wet] = 0.0

    dv_dx = np.gradient(v, DX, axis=0)
    du_dy = np.gradient(u, DY, axis=1)
    zeta = dv_dx - du_dy

    return np.ma.masked_array(zeta, mask=~wet)


# ---- Set up the live figure ----
plt.ion()
fig, ax = plt.subplots(figsize=(11, 6))

extent = [0, LX_KM, 0, LY_KM]

# Colour limits: scale with U/L for a typical shear.  Start generous
# and let the user see the full dynamic range.
ZETA_LIM = 4.0 * U_INFLOW / (ISLAND_R * DX)

import matplotlib as mpl
cmap = mpl.colormaps['RdBu_r'].copy()
cmap.set_bad('#333333')          # dark grey for the island

zeta0 = compute_vorticity(sim)
img = ax.imshow(
    zeta0.T,
    origin='lower', extent=extent, cmap=cmap,
    vmin=-ZETA_LIM, vmax=+ZETA_LIM, interpolation='nearest',
)
cbar = fig.colorbar(img, ax=ax, label=r'vorticity $\zeta = \partial_x v - \partial_y u$  (s$^{-1}$)')

# Coastline contour (b = 0) — the outline of the island
xs = (np.arange(NX) + 0.5) * DX / 1000.0
ys = (np.arange(NY) + 0.5) * DY / 1000.0
ax.contour(xs, ys, sim.b_phys.T, levels=[0.0], colors='black', linewidths=1.0)

ax.set_title("t = 0.0 min")
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
ax.set_aspect('equal')
fig.tight_layout()
plt.show()


def redraw(s):
    zeta = compute_vorticity(s)
    img.set_data(zeta.T)
    zmax = float(np.abs(zeta).max())
    ax.set_title(
        f"t = {s.time / 60.0:6.1f} min   "
        rf"$|\zeta|_{{\max}}$ = {zmax:.2e} s$^{{-1}}$"
    )
    fig.canvas.draw_idle()
    fig.canvas.flush_events()
    print(f"  t = {s.time / 60.0:6.1f} min   max |ζ| = {zmax:.3e} s^-1")


print("Animating 2D barotropic vorticity...")
sim.advance(
    T_END,
    callback=redraw,
    callback_interval=60.0,      # update every 60 sim-seconds
    progress=True,
)
print("Done.  Close the figure window to exit.")

plt.ioff()
plt.show()

sim.destroy()
