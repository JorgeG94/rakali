#!/usr/bin/env python3
"""Diagnostic plots for the river plume example."""
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ds = Dataset("output/output_rank_000000.nc")
t = ds.variables["time"][:]
x = ds.variables["x"][:]
y = ds.variables["y"][:]
b = ds.variables["b"][:]  # (y, x) from C-order NetCDF
dsig = ds.variables["dsig_ml"][:]

nt = len(t)
nz = len(dsig)

# Pick time indices: early, mid, late
tidx = [0, nt // 4, nt // 2, 3 * nt // 4, nt - 1]
tidx = sorted(set(tidx))

# ---- Figure 1: Surface salinity at multiple times ----
fig, axes = plt.subplots(len(tidx), 1, figsize=(14, 3 * len(tidx)), sharex=True, sharey=True)
if len(tidx) == 1:
    axes = [axes]

for ax, ti in zip(axes, tidx):
    # Surface = last z-index under k=1=bed, k=nz=surface (bottom-up).
    S = ds.variables["salinity"][ti, -1, :, :]  # surface layer (k=nz-1), (y, x)
    pcm = ax.pcolormesh(x, y, S, cmap="RdYlBu_r", vmin=0, vmax=35)
    ax.set_title(f"Surface salinity at t = {t[ti]:.0f} s ({t[ti]/60:.0f} min)")
    ax.set_ylabel("y (m)")
    ax.set_aspect("equal")
    fig.colorbar(pcm, ax=ax, label="PSU", shrink=0.8)

axes[-1].set_xlabel("x (m)")
plt.tight_layout()
plt.savefig("surface_salinity.png", dpi=150)
print("Saved surface_salinity.png")
plt.close()

# ---- Figure 2: Along-channel salinity cross-section (y = centre) ----
jc = len(y) // 2  # centre row

fig, axes = plt.subplots(len(tidx), 1, figsize=(14, 2.5 * len(tidx)), sharex=True)
if len(tidx) == 1:
    axes = [axes]

for ax, ti in zip(axes, tidx):
    h_total = ds.variables["h"][ti, :, :]  # (y, x)
    # Reverse the z axis so [0]=surface (matches the z_interface build
    # that walks from eta downward). Native layout is k=1=bed,
    # k=nz=surface, so reading [ti, ::-1, :, :] puts surface first.
    h_lay = ds.variables["h_layer"][ti, ::-1, :, :]  # (z_surf_first, y, x)
    S_lay = ds.variables["salinity"][ti, ::-1, :, :]  # (z_surf_first, y, x)

    # Build vertical grid at centre row
    eta = h_total[jc, :] + b[jc, :]  # surface elevation
    z_interface = np.zeros((nz + 1, len(x)))
    z_interface[0, :] = eta
    for k in range(nz):
        z_interface[k + 1, :] = z_interface[k, :] - h_lay[k, jc, :]

    # Cell centres
    z_centre = 0.5 * (z_interface[:-1, :] + z_interface[1:, :])

    # Clip salinity for display
    S_section = np.clip(S_lay[:, jc, :], 0, 36)

    # Plot as pcolormesh — need to construct (nz, nx) mesh
    # Use scatter for simplicity
    for k in range(nz):
        for ix in range(0, len(x), 1):
            ax.fill_between(
                [x[ix] - 5, x[ix] + 5],
                z_interface[k, ix], z_interface[k + 1, ix],
                color=plt.cm.RdYlBu_r(S_section[k, ix] / 35.0),
                linewidth=0,
            )

    # Bottom
    ax.fill_between(x, b[jc, :], b[jc, :].min() - 2, color="0.4")
    ax.set_title(f"Along-channel salinity at t = {t[ti]:.0f} s (y = {y[jc]:.0f} m)")
    ax.set_ylabel("z (m)")
    ax.set_ylim(b[jc, :].min() - 1, 2)

axes[-1].set_xlabel("x (m)")

# Add a colorbar
sm = plt.cm.ScalarMappable(cmap="RdYlBu_r", norm=plt.Normalize(0, 35))
fig.colorbar(sm, ax=axes, label="Salinity (PSU)", shrink=0.6)
plt.tight_layout()
plt.savefig("salinity_xsection.png", dpi=150)
print("Saved salinity_xsection.png")
plt.close()

# ---- Figure 3: Bottom salinity at multiple times ----
fig, axes = plt.subplots(len(tidx), 1, figsize=(14, 3 * len(tidx)), sharex=True, sharey=True)
if len(tidx) == 1:
    axes = [axes]

for ax, ti in zip(axes, tidx):
    S = ds.variables["salinity"][ti, -1, :, :]  # bottom layer, (y, x)
    S_clip = np.clip(S, 0, 36)
    pcm = ax.pcolormesh(x, y, S_clip, cmap="RdYlBu_r", vmin=0, vmax=35)
    ax.set_title(f"Bottom salinity at t = {t[ti]:.0f} s ({t[ti]/60:.0f} min)")
    ax.set_ylabel("y (m)")
    ax.set_aspect("equal")
    fig.colorbar(pcm, ax=ax, label="PSU", shrink=0.8)

axes[-1].set_xlabel("x (m)")
plt.tight_layout()
plt.savefig("bottom_salinity.png", dpi=150)
print("Saved bottom_salinity.png")
plt.close()

# ---- Figure 4: Surface velocity at last timestep ----
fig, ax = plt.subplots(1, 1, figsize=(14, 5))
ti = tidx[-1]
u = ds.variables["u_ml"][ti, 0, :, :]  # surface layer (y, x)
v = ds.variables["v_ml"][ti, 0, :, :]
speed = np.sqrt(u**2 + v**2)

pcm = ax.pcolormesh(x, y, speed, cmap="magma", vmin=0, vmax=1.0)
# Quiver (subsample)
skip = 5
ax.quiver(x[::skip], y[::skip], u[::skip, ::skip], v[::skip, ::skip],
          color="white", scale=10, width=0.002, alpha=0.7)
ax.set_title(f"Surface velocity at t = {t[ti]:.0f} s")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_aspect("equal")
fig.colorbar(pcm, ax=ax, label="Speed (m/s)")
plt.tight_layout()
plt.savefig("surface_velocity.png", dpi=150)
print("Saved surface_velocity.png")
plt.close()

# ---- Quick stats ----
print("\n=== Diagnostics ===")
for ti in [0, nt - 1]:
    S = ds.variables["salinity"][ti, :, :, :]
    h = ds.variables["h"][ti, :, :]
    hS_total = 0.0
    for k in range(nz):
        h_lay = ds.variables["h_layer"][ti, k, :, :]
        hS_total += np.sum(h_lay * S[k, :, :])
    print(f"t={t[ti]:7.0f}s: S range=[{S.min():8.2f}, {S.max():8.2f}], "
          f"total hS={hS_total:.2e}")

ds.close()
