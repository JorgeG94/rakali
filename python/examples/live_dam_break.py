#!/usr/bin/env python3
"""Live plotting example: classic 1D-style dam break in a 2D basin.

Sets up a rectangular basin with a virtual dam at x = nx/2: the
western half is filled to H_HIGH, the eastern half to H_LOW.  The
dam is released at t=0 and the resulting bore propagates eastward
across the basin.  Walls on all four sides reflect the wave back.

Inputs are loaded from ``live_dam_break.yaml`` next to this script.

Run with::

    python live_dam_break.py
"""
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Make the rakali package importable when running from the source tree
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "python"))

from rakali import Rakali
from rakali.yaml_loader import load_config, to_rakali_kwargs

cfg = load_config(Path(__file__).with_suffix(".yaml"))

NX, NY = cfg.domain.nx, cfg.domain.ny
DX, DY = cfg.domain.dx, cfg.domain.dy
T_END = cfg.domain.t_end
H_HIGH = cfg.initial_condition.dam_break.h_high
H_LOW = cfg.initial_condition.dam_break.h_low

sim = Rakali(**to_rakali_kwargs(cfg))

# Flat bottom (deep enough to hold the high water)
bathy = np.full((NX, NY), -H_HIGH, dtype=np.float64, order='F')
sim.set_bathymetry(bathy)

# Initial depth: H_HIGH on the western half, H_LOW on the eastern half.
# The "dam" is a vertical line at x = NX/2; the bore travels eastward.
split = int(NX * cfg.initial_condition.dam_break.split_fraction)
h_init = np.full((NX, NY), H_LOW, dtype=np.float64, order='F')
h_init[:split, :] = H_HIGH
sim.set_initial_h(h_init)

sim.finalize()
print(f"Created: {sim}")

# ---- Set up the live figure ----
plt.ion()
fig, ax = plt.subplots(figsize=(7, 6))

# imshow expects (rows, cols) = (y, x); transpose so y is vertical.
extent = [0, NX * DX, 0, NY * DY]
img = ax.imshow(
    sim.h_phys.T,
    origin='lower',
    extent=extent,
    cmap='viridis',
    vmin=H_LOW, vmax=H_HIGH,
    interpolation='nearest',
)
cbar = fig.colorbar(img, ax=ax, label='water depth h (m)')
title = ax.set_title(f"t = 0.00 s")
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_aspect('equal')
fig.tight_layout()
plt.show()

# ---- Run with callback that updates the image ----
def redraw(s):
    img.set_data(s.h_phys.T)
    title.set_text(f"t = {s.time:.2f} s   "
                   f"h ∈ [{s.h_phys.min():.2f}, {s.h_phys.max():.2f}] m")
    fig.canvas.draw_idle()
    fig.canvas.flush_events()

print("Animating dam break...")
sim.advance(T_END, callback=redraw, callback_interval=0.5)
print("Done.  Close the figure window to exit.")

plt.ioff()
plt.show()

sim.destroy()
