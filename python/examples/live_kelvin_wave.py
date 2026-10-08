#!/usr/bin/env python3
"""Live plotting example: coastally trapped Kelvin wave.

A rotating rectangular basin with walls on all four sides.  The
initial condition is a Gaussian free-surface bump pressed against
the south wall.  In the Northern Hemisphere (f > 0) the geostrophic
balance requires that free waves propagate with the coast on their
*right* — so this bump sheds a Kelvin wave that travels eastward
along the southern boundary, decaying exponentially away from the
wall with the Rossby radius of deformation L_R = sqrt(gH) / f.

Demonstrates:

* Coriolis rotation (``coriolis_f``)
* Barotropic SWE producing coastally trapped modes from a simple IC
* Diverging colormap for free-surface elevation

Inputs are loaded from ``live_kelvin_wave.yaml`` next to this script.

Run::

    python live_kelvin_wave.py
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

NX, NY = cfg.domain.nx, cfg.domain.ny
DX, DY = cfg.domain.dx, cfg.domain.dy
DEPTH = cfg.bathymetry.depth
F = cfg.forcing.coriolis_f
T_END = cfg.domain.t_end

LX_KM = NX * DX / 1000.0
LY_KM = NY * DY / 1000.0
C = math.sqrt(9.81 * DEPTH)
ROSSBY_R = C / F
print(f"Kelvin wave speed c = sqrt(gH) = {C:.1f} m/s")
print(f"Rossby radius of deformation = {ROSSBY_R / 1000.0:.1f} km")

sim = Rakali(**to_rakali_kwargs(cfg))

bathy = np.full((NX, NY), -DEPTH, dtype=np.float64, order='F')
sim.set_bathymetry(bathy)

# Gaussian bump near the south wall, near the west end.
bump_cfg = cfg.initial_condition.bump
AMP = bump_cfg.amplitude
X0 = bump_cfg.x0_fraction * LX_KM * 1000.0
SIGMA_X = bump_cfg.sigma_x
# Exponential decay away from the south wall with the Rossby radius.
ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
xx = (ii + 0.5) * DX
yy = (jj + 0.5) * DY
bump = AMP * np.exp(-((xx - X0) ** 2) / (2.0 * SIGMA_X ** 2)) \
           * np.exp(-yy / ROSSBY_R)

h_init = DEPTH + bump
sim.set_initial_h(np.asfortranarray(h_init))
sim.finalize()
print(f"Created: {sim}")

# ---- Figure ----
plt.ion()
fig, ax = plt.subplots(figsize=(11, 6))
extent = [0, LX_KM, 0, LY_KM]

eta0 = sim.h_phys + sim.b_phys
img = ax.imshow(eta0.T, origin='lower', extent=extent,
                cmap='RdBu_r', vmin=-AMP, vmax=AMP, interpolation='nearest')
fig.colorbar(img, ax=ax, label='η (m)')
ax.set_xlabel('x (km)'); ax.set_ylabel('y (km)'); ax.set_aspect('equal')
title = ax.set_title("t = 0.00 h")
fig.tight_layout()
plt.show()


def redraw(s):
    eta = s.h_phys + s.b_phys
    img.set_data(eta.T)
    title.set_text(f"t = {s.time / 3600.0:5.2f} h   "
                   f"η ∈ [{eta.min():+.2f}, {eta.max():+.2f}] m")
    fig.canvas.draw_idle()
    fig.canvas.flush_events()


print("Animating coastal Kelvin wave...")
sim.advance(T_END, callback=redraw, callback_interval=120.0, progress=True)
print("Done.  Close the figure window to exit.")

plt.ioff()
plt.show()
sim.destroy()
