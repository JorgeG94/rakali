#!/usr/bin/env python3
"""River plume on a structured Cartesian grid.

Companion to ``river_plume_unstr.py``. Same channel + basin geometry,
same physics knobs as ``river_plume.nml``, but driven entirely from
Python so both solvers can be run from a single environment and their
diagnostic logs compared line-for-line.

Layout (matches ``generate_bathymetry.py``)::

    +-----------------------------------------------+   y = 1500 m
    |                       | land (b = +5 m)       |
    |  channel (200 m wide) |-----------------------|
    |  4 m deep             | transition   basin    |   y = 750 m
    |                       |  (200 m)     (12 m)   |
    |                       |-----------------------|
    |                       | land                  |
    +-----------------------------------------------+   y = 0
    x = 0          x = 1000    x = 1200   x = 3000

Run from the repo root after building the shared library::

    cmake -B build -S . -DRAKALI_ENABLE_PYTHON=ON
    cmake --build build
    RAKALI_LIB=build/librakali.so PYTHONPATH=python \\
        python validation_examples/coastal/river_plume/river_plume_struct.py

Add ``--plot`` to save diagnostic PNGs; other flags mirror
``river_plume_unstr.py``.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

from rakali import Rakali

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from generate_bathymetry import generate_bathymetry  # noqa: E402


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--t-end", type=float, default=3600.0,
                        help="simulation end time in seconds (default: 3600)")
    parser.add_argument("--nz", type=int, default=5,
                        help="number of sigma layers (default: 5)")
    parser.add_argument("--dx", type=float, default=10.0,
                        help="cell size in x and y (default: 10 m)")
    parser.add_argument("--inflow-speed", type=float, default=0.5,
                        help="river inflow speed in m/s (default: 0.5)")
    parser.add_argument("--ambient-salinity", type=float, default=35.0,
                        help="ambient basin salinity in PSU (default: 35)")
    parser.add_argument("--river-salinity", type=float, default=0.0,
                        help="river inflow salinity in PSU (default: 0)")
    parser.add_argument("--beta-s", type=float, default=0.78,
                        help="EOS salinity coefficient (default: 0.78). "
                             "Set to 0 to make salinity a passive tracer "
                             "(no buoyancy feedback).")
    parser.add_argument("--advection", choices=["ppm", "upwind"],
                        default="ppm",
                        help="tracer horizontal advection scheme "
                             "(default: ppm). 'upwind' is first-order "
                             "for A/B comparison vs PPM.")
    parser.add_argument("--vcoord", choices=["sigma", "zsigma"],
                        default="sigma",
                        help="vertical coordinate: 'sigma' (default) or "
                             "'zsigma' (z-sigma hybrid, reduces PGE on "
                             "steep bathymetry)")
    parser.add_argument("--remap", choices=["pcm", "plm", "ppm"],
                        default="ppm",
                        help="vertical remapping method for zsigma "
                             "(default: ppm)")
    parser.add_argument("--depth-transition", type=float, default=6.0,
                        help="zsigma: depth (m) where blending starts "
                             "(default: 6, between channel 4 m and basin 12 m)")
    parser.add_argument("--blend-width", type=float, default=4.0,
                        help="zsigma: blending zone width in m (default: 4)")
    parser.add_argument("--output-interval", type=float, default=300.0,
                        help="seconds between status prints (default: 300)")
    parser.add_argument("--plot", action="store_true",
                        help="save diagnostic plots at the end")
    parser.add_argument("--channel-length", type=float, default=1000.0)
    parser.add_argument("--channel-width", type=float, default=200.0)
    parser.add_argument("--basin-length", type=float, default=2000.0)
    parser.add_argument("--basin-width", type=float, default=1500.0)
    parser.add_argument("--channel-depth", type=float, default=4.0)
    parser.add_argument("--basin-depth", type=float, default=12.0)
    args = parser.parse_args()

    # ---- Grid ----
    total_length = args.channel_length + args.basin_length
    nx = int(round(total_length / args.dx))
    ny = int(round(args.basin_width / args.dx))
    ncells_2d = nx * ny
    ncells_3d = ncells_2d * args.nz
    print(f"Building grid: {nx} x {ny} cells at dx = {args.dx} m "
          f"({nx * args.dx:.0f} x {ny * args.dx:.0f} m domain)")
    print(f"  {ncells_2d} 2D cells, {ncells_3d} 3D cells "
          f"({args.nz} sigma layers)")

    # ---- Bathymetry (matches river_plume.nml via generate_bathymetry.py) ----
    bathy = generate_bathymetry(nx=nx, ny=ny, dx=args.dx, dy=args.dx)
    bathy = np.asfortranarray(bathy.astype(np.float64))
    print(f"  Bathymetry: depth range [{bathy.min():.1f}, {bathy.max():.1f}] m")

    # ---- Simulation setup ----
    # Physics knobs match validation_examples/coastal/river_plume/river_plume.nml exactly so the
    # structured and unstructured drivers are configured the same way.
    sim = Rakali(
        nx=nx, ny=ny, dx=args.dx, dy=args.dx,
        t_end=args.t_end,
        cfl=0.45,
        manning_n=0.025,
        dry_threshold=0.01,
        # BCs
        bc_west="inflow",
        bc_east="open",
        bc_south="wall",
        bc_north="wall",
        inflow_velocity=args.inflow_speed,
        inflow_salinity=args.river_salinity,
        # Multilayer + BPG
        use_multilayer=True,
        nz_layers=args.nz,
        bpg_method="jacobian",
        ml_min_depth=0.5,
        rho_0=1025.0,
        # Horizontal tracer diffusion + Smagorinsky
        hdiff_kappa=1.0,
        smagorinsky_cs=0.15,
        # Pacanowski-Philander vertical mixing
        vmix_nu_max=0.01,
        vmix_nu_bg=1.0e-4,
        vmix_ri_alpha=5.0,
        # Tracers
        initial_salinity=args.ambient_salinity,
        beta_S=args.beta_s,
        kappa_S_bg=1.0e-5,
        alpha_T=0.0,            # temperature decoupled
        use_ppm_tracer=(args.advection == "ppm"),
        # Vertical coordinate
        vcoord_type=args.vcoord,
        remap_method=args.remap,
        vcoord_depth_transition=args.depth_transition,
        vcoord_blend_width=args.blend_width,
        # Silence per-step logging; we print our own diagnostic line
        log_level="warning",
        defer_init=True,
    )

    sim.set_bathymetry(bathy)

    # Still water at eta = 0 everywhere (dry on land, local depth elsewhere)
    h_init = np.asfortranarray(np.maximum(-bathy, 0.0))
    sim.set_initial_h(h_init)

    sim.finalize()


    # ---- Run ----
    vcoord_label = args.vcoord
    if args.vcoord == "zsigma":
        vcoord_label += (f" (transition={args.depth_transition}m, "
                         f"blend={args.blend_width}m, remap={args.remap})")
    print(f"\nRunning to t = {args.t_end} s with {args.nz} layers ...")
    print(f"  Vertical coord: {vcoord_label}")
    print(f"  Inflow: {args.inflow_speed} m/s, S = {args.river_salinity} PSU")
    print(f"  Ambient: S = {args.ambient_salinity} PSU")
    print(f"  Tracer advection: {args.advection.upper()}")
    print()

    t_next_print = args.output_interval
    wall_start = time.perf_counter()
    wall_window = wall_start
    steps_window = 0

    while sim.time < args.t_end:
        t_target = min(t_next_print, args.t_end)
        steps_before = sim.step_count
        sim.advance(t_target)
        steps_window += sim.step_count - steps_before
        if (sim.time >= t_next_print or sim.time >= args.t_end):
            wall_now = time.perf_counter()
            wall_total = wall_now - wall_start
            wall_dt = wall_now - wall_window
            sim_dt = args.output_interval if sim.time >= t_next_print else (
                sim.time - (t_next_print - args.output_interval))
            sps = steps_window / wall_dt if wall_dt > 0 else 0.0
            speedup = sim_dt / wall_dt if wall_dt > 0 else 0.0

            S = sim.salinity                    # (nx, ny, nz)
            h_phys = sim.h_phys                 # (nx, ny)
            dry = h_phys < 1.0e-6
            # Surface layer is k = 0, bottom layer is k = nz-1
            S_surf = np.ma.masked_array(S[:, :, 0], mask=dry)
            S_bot = np.ma.masked_array(S[:, :, -1], mask=dry)
            print(f"  t = {sim.time:8.1f} s  |  steps = {sim.step_count:6d}  |  "
                  f"dt = {sim.dt:.4e}  |  "
                  f"wall = {wall_total:7.2f} s ({sps:5.0f} steps/s, "
                  f"{speedup:5.1f}x rt)  |  "
                  f"S_surf [{S_surf.min():.1f}, {S_surf.max():.1f}]  "
                  f"S_bot [{S_bot.min():.1f}, {S_bot.max():.1f}]")
            # Localise the per-layer max -- ignore dry land via filling the
            # masked array with -inf before argmax.
            i_surf, j_surf = np.unravel_index(
                np.argmax(S_surf.filled(-np.inf)), S_surf.shape)
            i_bot, j_bot = np.unravel_index(
                np.argmax(S_bot.filled(-np.inf)), S_bot.shape)
            x_surf = (i_surf + 0.5) * sim.dx
            y_surf = (j_surf + 0.5) * sim.dy
            x_bot = (i_bot + 0.5) * sim.dx
            y_bot = (j_bot + 0.5) * sim.dy
            print(f"      S_surf max at (i={i_surf}, j={j_surf})  "
                  f"(x={x_surf:.1f}, y={y_surf:.1f}), "
                  f"h = {h_phys[i_surf, j_surf]:.2f} m  |  "
                  f"S_bot max at (i={i_bot}, j={j_bot})  "
                  f"(x={x_bot:.1f}, y={y_bot:.1f}), "
                  f"h = {h_phys[i_bot, j_bot]:.2f} m")
            wall_window = wall_now
            steps_window = 0
        if sim.time >= t_next_print or sim.time >= args.t_end:
            t_next_print += args.output_interval

    wall_total = time.perf_counter() - wall_start
    print(f"\nDone: {sim.step_count} steps, final dt = {sim.dt:.4e} s, "
          f"wall = {wall_total:.2f} s "
          f"({sim.step_count / wall_total:.0f} steps/s, "
          f"{args.t_end / wall_total:.1f}x real-time)")

    # ---- Plot ----
    if args.plot:
        plot_results(sim, args)

    sim.destroy()


def plot_results(sim, args):
    """Save diagnostic plots of the plume (structured)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plots")
        return

    nx, ny = sim.nx, sim.ny
    dx = sim.dx
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dx
    xe = np.arange(nx + 1) * dx
    ye = np.arange(ny + 1) * dx

    S = sim.salinity          # (nx, ny, nz)
    h_lay = sim.h_layer       # (nx, ny, nz)
    b = sim.b_phys            # (nx, ny)
    h = sim.h_phys            # (nx, ny)
    dry = h < 1.0e-6
    nz = sim.nz

    # ---- Surface salinity ----
    S_surf = np.ma.masked_array(S[:, :, 0], mask=dry)
    fig, ax = plt.subplots(1, 1, figsize=(14, 5))
    pcm = ax.pcolormesh(xe, ye, S_surf.T, cmap="RdYlBu_r", vmin=0, vmax=35)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"Surface salinity at t = {sim.time:.0f} s (structured)")
    ax.set_aspect("equal")
    fig.colorbar(pcm, ax=ax, label="Salinity (PSU)")
    plt.tight_layout()
    plt.savefig("validation_examples/coastal/river_plume/surface_salinity_struct.png", dpi=150)
    print("Saved validation_examples/coastal/river_plume/surface_salinity_struct.png")
    plt.close()

    # ---- Bottom salinity ----
    S_bot = np.ma.masked_array(S[:, :, -1], mask=dry)
    fig, ax = plt.subplots(1, 1, figsize=(14, 5))
    pcm = ax.pcolormesh(xe, ye, S_bot.T, cmap="RdYlBu_r", vmin=0, vmax=35)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"Bottom salinity at t = {sim.time:.0f} s (structured)")
    ax.set_aspect("equal")
    fig.colorbar(pcm, ax=ax, label="Salinity (PSU)")
    plt.tight_layout()
    plt.savefig("validation_examples/coastal/river_plume/bottom_salinity_struct.png", dpi=150)
    print("Saved validation_examples/coastal/river_plume/bottom_salinity_struct.png")
    plt.close()

    # ---- Along-channel cross-section (centreline row) ----
    jc = ny // 2
    b_band = b[:, jc]
    h_band = h_lay[:, jc, :]           # (nx, nz)
    S_band = S[:, jc, :]

    # Interface z levels, surface down
    z_iface = np.zeros((nx, nz + 1))
    z_iface[:, 0] = b_band + h_band.sum(axis=1)
    for k in range(nz):
        z_iface[:, k + 1] = z_iface[:, k] - h_band[:, k]

    fig, ax = plt.subplots(1, 1, figsize=(14, 5))
    sm = plt.cm.ScalarMappable(cmap="RdYlBu_r",
                               norm=plt.Normalize(0, 35))
    for k in range(nz):
        for j in range(nx - 1):
            color = sm.to_rgba(S_band[j, k])
            ax.fill_between(
                [x[j], x[j + 1]],
                [z_iface[j, k], z_iface[j + 1, k]],
                [z_iface[j, k + 1], z_iface[j + 1, k + 1]],
                color=color, linewidth=0,
            )
    ax.fill_between(x, b_band, b_band.min() - 1, color="0.4")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_title(f"Along-channel salinity cross-section at t = {sim.time:.0f} s "
                 f"(structured)")
    ax.set_ylim(b_band.min() - 1, 2)
    fig.colorbar(sm, ax=ax, label="Salinity (PSU)")
    plt.tight_layout()
    plt.savefig("validation_examples/coastal/river_plume/xsection_salinity_struct.png", dpi=150)
    print("Saved validation_examples/coastal/river_plume/xsection_salinity_struct.png")
    plt.close()

    # ---- Surface velocity ----
    # Grab layer-wise momentum via the 3D accessor (no Python property yet).
    hu_lay = sim._phys_slice(sim._get_3d_field("rki_get_hu_layer_ptr"))
    hv_lay = sim._phys_slice(sim._get_3d_field("rki_get_hv_layer_ptr"))
    h_top = h_lay[:, :, 0]
    u_surf = np.where(h_top > 1e-6, hu_lay[:, :, 0] / h_top, 0.0)
    v_surf = np.where(h_top > 1e-6, hv_lay[:, :, 0] / h_top, 0.0)
    speed = np.sqrt(u_surf**2 + v_surf**2)
    speed_m = np.ma.masked_array(speed, mask=dry)

    fig, ax = plt.subplots(1, 1, figsize=(14, 5))
    pcm = ax.pcolormesh(xe, ye, speed_m.T, cmap="magma", vmin=0, vmax=1.0)
    skip = max(1, nx // 60)
    ax.quiver(x[::skip], y[::skip],
              u_surf[::skip, ::skip].T, v_surf[::skip, ::skip].T,
              color="white", scale=10, width=0.002, alpha=0.7)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"Surface velocity at t = {sim.time:.0f} s (structured)")
    ax.set_aspect("equal")
    fig.colorbar(pcm, ax=ax, label="Speed (m/s)")
    plt.tight_layout()
    plt.savefig("validation_examples/coastal/river_plume/surface_velocity_struct.png", dpi=150)
    print("Saved validation_examples/coastal/river_plume/surface_velocity_struct.png")
    plt.close()

    print("\nAll plots saved.")


if __name__ == "__main__":
    main()
