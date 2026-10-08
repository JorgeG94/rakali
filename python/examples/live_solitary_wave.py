#!/usr/bin/env python3
"""Live plotting example: shape-preserving solitary wave.

A true solitary wave is a delicate balance between nonlinear
steepening and frequency dispersion.  The correct initial
condition is **not** a free-surface bump with u = 0 (that just
splits into two halves and flattens out), but a sech² profile
with the velocity field chosen to match:

    η(x)    = A sech²(κ (x - x₀))
    u(x)    = c η / (H + η)        with c = sqrt(g(H + A))
    κ       = sqrt(3 A / (4 H³))

Released at t=0 the wave should propagate rightward without
visibly changing shape — that's the whole point.  Run it with
``use_nonhydrostatic=True`` for the dispersive balance; with NH
off the same IC would noticeably deform as the hydrostatic
equations lack the dispersive term that holds the sech² shape
together.

Inputs are loaded from ``live_solitary_wave.yaml`` next to this script.

Run::

    python live_solitary_wave.py
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
NZ = cfg.domain.nz
H0 = cfg.bathymetry.depth
AMP = cfg.initial_condition.solitary.amplitude
X0 = cfg.initial_condition.solitary.x0
T_END = cfg.domain.t_end
G = 9.81

kappa = math.sqrt(3.0 * AMP / (4.0 * H0 ** 3))
c = math.sqrt(G * (H0 + AMP))
print(f"Solitary wave: A/H = {AMP / H0:.2f}, κ = {kappa:.4f} /m, "
      f"c = {c:.3f} m/s")
print(f"Crossing time ≈ {NX * DX / c:.1f} s")

sim = Rakali(**to_rakali_kwargs(cfg))

bathy = np.full((NX, NY), -H0, dtype=np.float64, order='F')
sim.set_bathymetry(bathy)

x = (np.arange(NX) + 0.5) * DX
eta_1d = AMP / np.cosh(kappa * (x - X0)) ** 2      # sech²
u_1d = c * eta_1d / (H0 + eta_1d)                  # matched velocity
h_1d = H0 + eta_1d
hu_1d = h_1d * u_1d

h_init = np.broadcast_to(h_1d[:, None], (NX, NY)).copy()
sim.set_initial_h(np.asfortranarray(h_init))
sim.finalize()

# Matched momentum — set after finalize with set_hu.
hu_init = np.broadcast_to(hu_1d[:, None], (NX, NY)).copy()
sim.set_hu(np.asfortranarray(hu_init))

print(f"Created: {sim}")

# ---- Figure: η profile along the channel ----
plt.ion()
fig, ax = plt.subplots(figsize=(12, 4.5))
xs = x
J_MID = NY // 2

eta0 = sim.h_phys[:, J_MID] + sim.b_phys[:, J_MID]
(line_init,) = ax.plot(xs, eta0, 'k--', lw=0.8, alpha=0.4,
                       label='initial η')
(line,) = ax.plot(xs, eta0, lw=1.6, color='tab:blue', label='η(t)')
ax.set_xlim(xs[0], xs[-1])
ax.set_ylim(-0.2, AMP + 0.2)
ax.set_xlabel('x (m)'); ax.set_ylabel('η (m)')
ax.axhline(0.0, color='k', lw=0.5)
ax.grid(alpha=0.3)
ax.legend(loc='upper right')
title = ax.set_title("t = 0.0 s")
fig.tight_layout()
plt.show()


def redraw(s):
    eta = s.h_phys[:, J_MID] + s.b_phys[:, J_MID]
    line.set_ydata(eta)
    title.set_text(f"t = {s.time:5.2f} s   "
                   f"peak η = {eta.max():.3f} m   "
                   f"x_peak ≈ {xs[eta.argmax()]:.1f} m")
    fig.canvas.draw_idle()
    fig.canvas.flush_events()


print("Animating shape-preserving solitary wave...")
sim.advance(T_END, callback=redraw, callback_interval=0.5, progress=True)
print("Done.  Close the figure window to exit.")

plt.ioff()
plt.show()
sim.destroy()
