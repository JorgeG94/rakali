#!/usr/bin/env python3
"""Live plotting example: baroclinic lock exchange.

Classic gravity-current benchmark.  A vertical "gate" at x = L/2
separates dense (S_LEFT) and light (S_RIGHT) water in a closed,
flat-bottomed channel.  Releasing the gate launches two
counter-propagating noses: a dense underflow moving right and a
light overflow moving left.  Kelvin-Helmholtz billows along the
interface give the characteristic "teeth" pattern.

Demonstrates:

* Multilayer sigma coordinates with the Jacobian BPG
* ``set_salinity()`` for per-layer initial stratification
* Smagorinsky horizontal viscosity + horizontal salinity diffusion
  to damp grid-scale noise without smearing the front
* Pseudo-2D (x, z) salinity cross-section rebuilt from h_layer

Inputs are loaded from ``live_lock_exchange.yaml`` next to this script.

Run::

    python live_lock_exchange.py
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
NZ = cfg.domain.nz
DEPTH = cfg.bathymetry.depth
T_END = cfg.domain.t_end
S_LEFT = cfg.salinity.s_left
S_RIGHT = cfg.salinity.s_right

sim = Rakali(**to_rakali_kwargs(cfg))

bathy = np.full((NX, NY), -DEPTH, dtype=np.float64, order='F')
sim.set_bathymetry(bathy)
h_init = np.full((NX, NY), DEPTH, dtype=np.float64, order='F')
sim.set_initial_h(h_init)
sim.finalize()

# Vertical gate at x = NX/2: dense on the left, light on the right.
S = np.full((NX, NY, NZ), 0.5 * (S_LEFT + S_RIGHT), dtype=np.float64)
S[:NX // 2, :, :] = S_LEFT
S[NX // 2:, :, :] = S_RIGHT
sim.set_salinity(np.asfortranarray(S))

print(f"Created: {sim}")
print(f"  Lock exchange: ΔS = {S_LEFT - S_RIGHT:.1f} PSU, NZ = {NZ}")

# ---- Build a (x, z) cross-section from h_layer at mid-channel j ----
J_MID = NY // 2


def cross_section(s):
    """Return (X_edges, Z_edges, S_slice) for pcolormesh."""
    # k=1=bed, k=nz=surface (bottom-up). Reverse so axis index 0 ==
    # surface layer; that way cumsum from index 0 walks the column
    # downward (surface -> bed) which is what the z_interfaces build
    # below assumes.
    h_lay = s.h_layer[:, J_MID, ::-1]       # (nx, nz), [0]=surface
    salt = s.salinity[:, J_MID, ::-1]       # (nx, nz), [0]=surface
    # Build z interfaces per column: top at eta, going down.
    eta = s.h_phys[:, J_MID] + s.b_phys[:, J_MID]   # eta = h + b
    z_top = eta
    # Layer k thickness = h_lay[:, k]; bottom of layer k at z_top - cumsum
    z_interfaces = np.empty((NX, NZ + 1))
    z_interfaces[:, 0] = z_top
    z_interfaces[:, 1:] = z_top[:, None] - np.cumsum(h_lay, axis=1)

    x_edges = np.arange(NX + 1) * DX
    # pcolormesh wants 2D X, Z edge grids.  Stagger x, replicate z.
    X = np.repeat(x_edges[:, None], NZ + 1, axis=1)          # (nx+1, nz+1)
    Z = np.empty_like(X)
    Z[:-1, :] = z_interfaces
    Z[-1, :] = z_interfaces[-1, :]
    return X, Z, salt


# ---- Live figure ----
plt.ion()
fig, ax = plt.subplots(figsize=(11, 4.5))
X0, Z0, S0 = cross_section(sim)
mesh = ax.pcolormesh(X0, Z0, S0, cmap='RdBu_r',
                     vmin=S_RIGHT, vmax=S_LEFT, shading='flat')
fig.colorbar(mesh, ax=ax, label='salinity (PSU)')
ax.set_xlabel('x (m)')
ax.set_ylabel('z (m)')
ax.set_ylim(-DEPTH - 0.5, 1.0)
title = ax.set_title("t = 0.0 s")
fig.tight_layout()
plt.show()


def redraw(s):
    X, Z, salt = cross_section(s)
    # pcolormesh with shading='flat' needs a flat 2D C of shape (nx, nz)
    mesh.set_array(salt.ravel())
    title.set_text(f"t = {s.time:6.1f} s")
    fig.canvas.draw_idle()
    fig.canvas.flush_events()


print("Animating lock exchange...")
sim.advance(T_END, callback=redraw, callback_interval=1.0, progress=True)
print("Done.  Close the figure window to exit.")

plt.ioff()
plt.show()
sim.destroy()
