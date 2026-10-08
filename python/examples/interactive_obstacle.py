#!/usr/bin/env python3
"""Interactive example: click to drop splashes into a flowing channel.

A steady westerly inflow drives water through a rectangular
channel.  **Left-click** anywhere on the plot to drop an
instantaneous Gaussian free-surface bump at that location — it
radiates as gravity waves immediately, advected downstream by
the mean flow.  **Right-click** a location to drop a *negative*
bump (a hole).  Press **'r'** to reset the free surface to flat.

Unlike an atmospheric-pressure low (which would take ~L/c to
equilibrate via the inverse-barometer response), a direct free-
surface perturbation is instantly visible and radiates waves
within a few time steps.  This uses the post-finalize
``set_h`` setter which is designed for mid-simulation state
modifications.

Demonstrates:

* Mid-simulation state modification via ``set_h`` (read-modify-write
  on the zero-copy ``h_phys`` view)
* Matplotlib mouse/key events wired into the ``callback=`` hook
* Live surface elevation rendering

Inputs are loaded from ``interactive_obstacle.yaml`` next to this script.

Run::

    python interactive_obstacle.py
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
DEPTH = cfg.bathymetry.depth
T_END = cfg.domain.t_end

sim = Rakali(**to_rakali_kwargs(cfg))

bathy = np.full((NX, NY), -DEPTH, dtype=np.float64, order='F')
sim.set_bathymetry(bathy)
sim.set_initial_h(np.asfortranarray(np.full((NX, NY), DEPTH)))
sim.finalize()

print(f"Created: {sim}")
print("Left-click = splash (+η bump)   "
      "Right-click = hole (-η bump)   "
      "'r' = reset to flat")

# ---- Figure ----
plt.ion()
fig, ax = plt.subplots(figsize=(11, 6))
LX_KM = NX * DX / 1000.0
LY_KM = NY * DY / 1000.0
extent = [0, LX_KM, 0, LY_KM]

eta = sim.h_phys + sim.b_phys
ETA_LIM = 0.5
img = ax.imshow(eta.T, origin='lower', extent=extent,
                cmap='RdBu_r', vmin=-ETA_LIM, vmax=ETA_LIM,
                interpolation='nearest')
fig.colorbar(img, ax=ax, label='η (m)')
ax.set_xlabel('x (km)'); ax.set_ylabel('y (km)'); ax.set_aspect('equal')
title = ax.set_title("t = 0.0 min   "
                     "(left-click: splash   right-click: hole   r: reset)")

# ---- Mouse wiring ----
# A "splash" is an instantaneous Gaussian bump added to the water
# depth h.  It immediately radiates as gravity waves (±sqrt(gH)).
SPLASH_AMP = 0.35               # m peak perturbation
SPLASH_RADIUS_CELLS = 6
XX, YY = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')


def add_splash(i0, j0, sign=+1.0):
    """Add a Gaussian η bump to the running state at cell (i0, j0)."""
    r2 = (XX - i0) ** 2 + (YY - j0) ** 2
    bump = sign * SPLASH_AMP * np.exp(-r2 / (2.0 * SPLASH_RADIUS_CELLS ** 2))
    # Read current h, add bump, write back.  Mask out dry cells so
    # we don't create water where there is none.
    h_cur = sim.h_phys.copy()
    wet = h_cur > 1e-6
    h_new = h_cur + bump * wet
    # Clip to a safe minimum so negative splashes can't go below zero.
    np.maximum(h_new, 1e-6, out=h_new)
    sim.set_h(np.asfortranarray(h_new))


def reset_surface():
    """Restore a flat free surface at η = 0."""
    h_flat = np.maximum(-sim.b_phys, 0.0).astype(np.float64)
    sim.set_h(np.asfortranarray(h_flat))


def on_click(event):
    if event.inaxes is not ax or event.xdata is None:
        return
    i = int(event.xdata * 1000.0 / DX)
    j = int(event.ydata * 1000.0 / DY)
    if not (0 <= i < NX and 0 <= j < NY):
        return
    if event.button == 1:
        add_splash(i, j, +1.0)
        print(f"  splash at ({event.xdata:.2f}, {event.ydata:.2f}) km")
    elif event.button == 3:
        add_splash(i, j, -1.0)
        print(f"  hole at ({event.xdata:.2f}, {event.ydata:.2f}) km")


def on_key(event):
    if event.key == 'r':
        reset_surface()
        print("  reset to flat surface")


fig.canvas.mpl_connect('button_press_event', on_click)
fig.canvas.mpl_connect('key_press_event', on_key)

fig.tight_layout()
plt.show()


def redraw(s):
    eta = s.h_phys + s.b_phys
    img.set_data(eta.T)
    title.set_text(f"t = {s.time / 60.0:6.1f} min   "
                   f"η ∈ [{eta.min():+.2f}, {eta.max():+.2f}] m")
    fig.canvas.draw_idle()
    fig.canvas.flush_events()


print("Running.  Click the plot to drop lows; close the window to exit.")
sim.advance(T_END, callback=redraw, callback_interval=1.0, progress=True)

plt.ioff()
plt.show()
sim.destroy()
