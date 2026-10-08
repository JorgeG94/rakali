#!/usr/bin/env python3
"""River plume on an unstructured triangular mesh.

A fresh-water river channel (200 m wide, 4 m deep) empties into a
saline basin (12 m deep). The buoyant plume spreads across the surface
while a salt wedge intrudes along the bottom.

Uses the Python API to drive the unstructured multilayer solver with
5 sigma layers. Exercises the full API surface: mesh generation,
multilayer, EOS, boundary salinity, inflow BC, and tracer output.

Run from the repo root after building the shared library::

    cmake -B build -S . -DRAKALI_ENABLE_GPU=OFF -DRAKALI_ENABLE_PYTHON=ON
    cmake --build build
    RAKALI_LIB=build/librakali.so PYTHONPATH=python python validation_examples/coastal/river_plume/river_plume_unstr.py

Add ``--plot`` to save diagnostic PNGs, ``--t-end 3600`` for a longer
run, ``--max-area 200`` for coarser mesh (faster), etc.
"""
from __future__ import annotations

import argparse
import time
import numpy as np
import sys
import os

from rakali.unstr import (
    UnstructuredMesh, UnstructuredSimulation,
    set_log_level, LOG_INFO,
)

GRAV = 9.81


def compute_inradii(mesh):
    """Per-triangle inscribed-circle radius, matching the Fortran kernel
    convention `r = 2A/perimeter`. Returns a (ntri,) float64 array."""
    idx = mesh.tri_v.astype(np.int64) - 1   # tri_v is 1-based
    px = mesh.vert_x[idx]                   # (ntri, 3)
    py = mesh.vert_y[idx]
    e01 = np.hypot(px[:, 1] - px[:, 0], py[:, 1] - py[:, 0])
    e12 = np.hypot(px[:, 2] - px[:, 1], py[:, 2] - py[:, 1])
    e20 = np.hypot(px[:, 0] - px[:, 2], py[:, 0] - py[:, 2])
    perim = e01 + e12 + e20
    cross = ((px[:, 1] - px[:, 0]) * (py[:, 2] - py[:, 0])
             - (px[:, 2] - px[:, 0]) * (py[:, 1] - py[:, 0]))
    area = 0.5 * np.abs(cross)
    return 2.0 * area / np.maximum(perim, 1e-30)


def build_mesh(
    channel_length: float = 1000.0,
    channel_width: float = 200.0,
    basin_length: float = 2000.0,
    basin_width: float = 1500.0,
    max_area: float = 200.0,
    min_angle: float = 33.0,
) -> UnstructuredMesh:
    """Build a channel + basin mesh with tagged boundary groups.

    Layout (not to scale)::

        wall (group 1)
        ┌──────────────────────────────────────────┐
        │                                          │
        │  channel (1000m)  │   basin (2000m)      │
        │                   │                      │
   inflow (group 2)        transition      open (group 3)
        │                   │                      │
        │                                          │
        └──────────────────────────────────────────┘
        wall (group 1)

    The channel is centred in the y-direction. The basin spans the
    full width east of the transition zone.

    BC groups:
      1 = wall (top, bottom, channel walls)
      2 = inflow (west end of channel)
      3 = open (east end of basin)
    """
    # Domain is a polygon: channel on the left, basin on the right
    # We build a polygon that traces the outer boundary
    y_mid = basin_width / 2.0
    half_ch = channel_width / 2.0
    x_basin = channel_length + basin_length

    # Walk the polygon CCW, starting at bottom-left of channel
    # Channel bottom wall
    p0 = (0.0, y_mid - half_ch)           # bottom-left of channel (inflow corner)
    p1 = (channel_length, y_mid - half_ch) # bottom-right of channel
    # Basin bottom wall
    p2 = (channel_length, 0.0)             # bottom-left of basin
    p3 = (x_basin, 0.0)                    # bottom-right of basin (open corner)
    # Basin right wall (open)
    p4 = (x_basin, basin_width)            # top-right of basin (open corner)
    # Basin top wall
    p5 = (channel_length, basin_width)     # top-left of basin
    # Channel top wall
    p6 = (channel_length, y_mid + half_ch) # top-right of channel
    p7 = (0.0, y_mid + half_ch)            # top-left of channel (inflow corner)

    points = [p0, p1, p2, p3, p4, p5, p6, p7]
    # Facets connect consecutive points; the last connects back to 0
    facets = [(i, (i + 1) % len(points)) for i in range(len(points))]

    # Edge -> group mapping:
    #   0: p0->p1  = channel bottom wall  -> group 1
    #   1: p1->p2  = channel-basin step S -> group 1
    #   2: p2->p3  = basin bottom wall    -> group 1
    #   3: p3->p4  = basin east (open)    -> group 3
    #   4: p4->p5  = basin top wall       -> group 1
    #   5: p5->p6  = channel-basin step N -> group 1
    #   6: p6->p7  = channel top wall     -> group 1
    #   7: p7->p0  = channel west (inflow)-> group 2
    markers = [1, 1, 1, 3, 1, 1, 1, 2]

    try:
        from meshpy import triangle as mt
    except ImportError:
        raise SystemExit(
            "meshpy required: conda install -c conda-forge meshpy"
        )

    info = mt.MeshInfo()
    info.set_points(points)
    info.set_facets(facets, facet_markers=markers)

    mesh_data = mt.build(
        info,
        max_volume=float(max_area),
        min_angle=float(min_angle),
        generate_faces=True,
    )

    pts = np.asarray(mesh_data.points, dtype=np.float64)
    elems = np.asarray(mesh_data.elements, dtype=np.int64)

    vert_x = np.ascontiguousarray(pts[:, 0])
    vert_y = np.ascontiguousarray(pts[:, 1])

    # Ensure CCW orientation
    x = vert_x[elems]
    y = vert_y[elems]
    signed = 0.5 * (
        (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0])
        - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    )
    cw = signed < 0.0
    if np.any(cw):
        elems[cw] = elems[cw][:, [0, 2, 1]]

    ntri = elems.shape[0]
    tri_v = np.ascontiguousarray((elems + 1).astype(np.int32))

    # Build bc_tags from meshpy face markers
    bc_arr = np.zeros((ntri, 3), dtype=np.int32)
    if hasattr(mesh_data, "faces") and len(mesh_data.faces) > 0:
        face_pairs = np.asarray(mesh_data.faces, dtype=np.int64)
        face_marks = np.asarray(mesh_data.face_markers, dtype=np.int32)
        edge_marker = {}
        for (a, b), m in zip(face_pairs, face_marks):
            if m != 0:
                edge_marker[frozenset((int(a), int(b)))] = int(m)
        if edge_marker:
            opp = [(1, 2), (2, 0), (0, 1)]
            for k in range(ntri):
                for i, (va_idx, vb_idx) in enumerate(opp):
                    va = int(elems[k, va_idx])
                    vb = int(elems[k, vb_idx])
                    m = edge_marker.get(frozenset((va, vb)), 0)
                    if m != 0:
                        bc_arr[k, i] = m

    return UnstructuredMesh(vert_x, vert_y, tri_v, bc_tags=bc_arr)


def bathymetry(cx: np.ndarray, cy: np.ndarray,
               channel_length: float = 1000.0,
               channel_width: float = 200.0,
               basin_width: float = 1500.0,
               channel_depth: float = 4.0,
               basin_depth: float = 12.0,
               transition_length: float = 200.0) -> np.ndarray:
    """Compute bed elevation for each cell centroid.

    Returns b < 0 (below sea level) everywhere in the domain.
    The channel is shallower; the basin is deeper; a smooth ramp
    connects them over `transition_length`.
    """
    n = cx.size
    b = np.full(n, -channel_depth)
    y_mid = basin_width / 2.0

    for i in range(n):
        x = cx[i]
        if x >= channel_length + transition_length:
            # Full basin
            b[i] = -basin_depth
        elif x >= channel_length:
            # Transition zone: smooth depth ramp
            frac = (x - channel_length) / transition_length
            b[i] = -(channel_depth + frac * (basin_depth - channel_depth))
        # else: channel depth (already set)

    return b


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--t-end", type=float, default=3600.0,
                        help="simulation end time in seconds (default: 3600)")
    parser.add_argument("--nz", type=int, default=5,
                        help="number of sigma layers (default: 5)")
    parser.add_argument("--max-area", type=float, default=200.0,
                        help="max triangle area for mesh (default: 200 m^2)")
    parser.add_argument("--min-angle", type=float, default=33.0,
                        help="min triangle angle in degrees (default: 33)")
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
    parser.add_argument("--bpg-method", choices=["sigma", "jacobian"],
                        default="jacobian",
                        help="Baroclinic pressure gradient method "
                             "(default: jacobian). 'jacobian' eliminates the "
                             "sigma PGE on steep bathymetry; 'sigma' is the "
                             "legacy same-σ-level differencing.")
    parser.add_argument("--advection", choices=["muscl", "upwind"],
                        default="muscl",
                        help="tracer horizontal advection scheme "
                             "(default: muscl). 'upwind' is first-order "
                             "for A/B comparison vs MUSCL.")
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
                             "(default: 6)")
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

    set_log_level(LOG_INFO)

    # ---- Mesh ----
    print(f"Building mesh: channel {args.channel_length}x{args.channel_width} m, "
          f"basin {args.basin_length}x{args.basin_width} m, "
          f"max_area = {args.max_area} m^2")

    mesh = build_mesh(
        channel_length=args.channel_length,
        channel_width=args.channel_width,
        basin_length=args.basin_length,
        basin_width=args.basin_width,
        max_area=args.max_area,
        min_angle=args.min_angle,
    )
    print(f"  {mesh.nverts} vertices, {mesh.ntri} triangles "
          f"({mesh.ntri * args.nz} 3D cells, {args.nz} sigma layers)")

    # ---- Mesh quality diagnostic ----------------------------------
    radii = compute_inradii(mesh)
    cx_mesh, cy_mesh = mesh.centroids()
    k_min = int(np.argmin(radii))
    print(f"  radius (m): min = {radii.min():.3f}, "
          f"p1 = {np.percentile(radii, 1):.3f}, "
          f"median = {np.percentile(radii, 50):.3f}, "
          f"max = {radii.max():.3f}")
    print(f"  smallest triangle: idx = {k_min} at "
          f"(x={cx_mesh[k_min]:.1f}, y={cy_mesh[k_min]:.1f}), "
          f"r = {radii[k_min]:.3f} m")

    # Diagnostic mesh helpers aligned with `sim.centroid_x` / `sim.height`.
    n_owned = mesh.ntri
    radii = compute_inradii(mesh)[:n_owned]
    cx_full, cy_full = mesh.centroids()
    cx_mesh = cx_full[:n_owned]
    cy_mesh = cy_full[:n_owned]

    # ---- Simulation setup ----
    sim = UnstructuredSimulation(mesh)
    sim.enable_multilayer(args.nz)

    cx = sim.centroid_x.copy()
    cy = sim.centroid_y.copy()

    # Bathymetry
    bed = bathymetry(
        cx, cy,
        channel_length=args.channel_length,
        channel_width=args.channel_width,
        basin_width=args.basin_width,
        channel_depth=args.channel_depth,
        basin_depth=args.basin_depth,
    )
    sim.set_bathymetry(bed)

    # Initial condition: still water at sea level
    stage = np.zeros(sim.ncells)
    sim.set_initial_stage(stage)

    # Boundary conditions
    # group 1 = wall (default)
    sim.set_bc_type(group=1, bc_type="wall")
    # group 2 = inflow (river)
    sim.set_bc_type(group=2, bc_type="inflow")
    sim.set_bc_inflow_speed(group=2, speed=args.inflow_speed)
    sim.set_bc_h_ext(group=2, h=args.channel_depth)
    # group 3 = open (ocean, Flather radiation)
    sim.set_bc_type(group=3, bc_type="open")

    # EOS: salinity-driven buoyancy, temperature passive
    sim.set_eos(rho_0=1025.0, beta_s=args.beta_s, s_ref=0.0,
                alpha_t=0.0, t_ref=15.0)
    sim.set_bpg_method(args.bpg_method)
    if args.vcoord != "sigma":
        sim.set_vcoord(coord_type=args.vcoord, remap_method=args.remap,
                        depth_transition=args.depth_transition,
                        blend_width=args.blend_width)

    # Tracers
    sim.set_initial_salinity(args.ambient_salinity)
    sim.set_bc_salinity(group=2, value=args.river_salinity)  # fresh water inflow
    sim.set_bc_salinity(group=3, value=args.ambient_salinity)  # ocean boundary

    # Physics
    sim.set_manning(0.025)
    cfl = 0.45
    sim.set_cfl(cfl)

    # Horizontal viscosity & diffusion
    sim.set_smagorinsky_cs(0.15)
    sim.set_hdiff_kappa(1.0)

    # Vertical mixing (Pacanowski-Philander)
    sim.set_vmix(nu_max=0.01, nu_bg=1e-4, ri_alpha=5.0,
                 kappa_s_bg=1e-5, kappa_t_bg=1e-5)

    # Tracer advection scheme (MUSCL by default; upwind for A/B compare)
    sim.set_use_ppm_tracer(args.advection == "muscl")

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
    sim.prepare()

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
            # Grab salinity snapshot for diagnostics.
            # Unstructured layer convention: k = 0 (Python) is the BOTTOM
            # sigma layer; k = -1 is the SURFACE/TOP layer. This is the
            # OPPOSITE of the structured solver, see
            # rki_kernel_multilayer_unstr.F90:18.
            S = sim.salinity
            S_bot = S[:, 0]    # k = 1 in Fortran -> bottom
            S_surf = S[:, -1]  # k = nz in Fortran -> surface

            # Barotropic dt floor: cfl * min(r_k / sqrt(g h_k)) over wet cells.
            # If dt_obs == dt_floor (ratio = 1) the limiter is a wet cell
            # running at its barotropic wave speed -- the mesh / depth is the
            # culprit. If dt_obs << dt_floor (ratio << 1) some cell has wave
            # speed >> sqrt(g h) -- baroclinic / BPG noise. The argmin tells
            # us where the limiter sits.
            h = sim.height
            wet = h > 1e-3
            if wet.any():
                bt_speed = np.sqrt(GRAV * np.where(wet, h, np.inf))
                cell_dt = radii / bt_speed              # per-cell barotropic CFL contribution
                k_lim = int(np.argmin(cell_dt))
                dt_floor = cfl * cell_dt[k_lim]
                ratio = sim.dt / dt_floor if dt_floor > 0 else float('nan')
            else:
                k_lim = -1
                dt_floor = float('nan')
                ratio = float('nan')

            print(f"  t = {sim.time:8.1f} s  |  steps = {sim.step_count:6d}  |  "
                  f"dt = {sim.dt:.4e}  |  "
                  f"wall = {wall_total:7.2f} s ({sps:5.0f} steps/s, "
                  f"{speedup:5.1f}x rt)  |  "
                  f"S_surf [{S_surf.min():.1f}, {S_surf.max():.1f}]  "
                  f"S_bot [{S_bot.min():.1f}, {S_bot.max():.1f}]")
            k_surf = int(np.argmax(S_surf))
            k_bot = int(np.argmax(S_bot))
            print(f"      S_surf max at cell {k_surf} "
                  f"(x={cx_mesh[k_surf]:.1f}, y={cy_mesh[k_surf]:.1f}), "
                  f"h = {h[k_surf]:.2f} m  |  "
                  f"S_bot max at cell {k_bot} "
                  f"(x={cx_mesh[k_bot]:.1f}, y={cy_mesh[k_bot]:.1f}), "
                  f"h = {h[k_bot]:.2f} m")

            # ---- Tracer-conservation + baroclinic-shear diagnostic ----
            # 1) Total hS mass: sum_{cell, layer} hS_layer * area.
            #    The CWC scheme should preserve this exactly up to BC
            #    inflow/outflow. Drift >> 1e-10 ratio -> conservation
            #    violation in a kernel; tight ratio -> overshoot is a
            #    redistribution problem (front instability), not loss.
            # 2) Max |u_surf - u_bot|: the strongest baroclinic shear
            #    in the domain. With the depth-averaged tracer flux
            #    (dsig_k * mass_flux_e), shear of O(0.1-1 m/s) means
            #    the per-layer transport is being approximated as the
            #    depth-averaged speed -- the suspected source of the
            #    column-mode overshoot.
            h_lay = sim.h_layer                   # (ncells, nz)
            hu_lay = sim.hu_layer
            area = sim.area                       # (ncells,)
            hS_lay = sim.hS_layer

            total_hS = float((hS_lay * area[:, None]).sum())
            if 'total_hS_initial' not in locals():
                total_hS_initial = total_hS
            drift_rel = (total_hS - total_hS_initial) / max(abs(total_hS_initial), 1e-30)

            with np.errstate(divide='ignore', invalid='ignore'):
                u_lay = np.where(h_lay > 1e-6, hu_lay / h_lay, 0.0)
            shear = np.abs(u_lay[:, -1] - u_lay[:, 0])  # |u_surf - u_bot|
            wet = h > 0.5
            shear_max = float(shear[wet].max()) if wet.any() else 0.0
            k_shear = int(np.argmax(np.where(wet, shear, -np.inf)))
            print(f"      total_hS = {total_hS:.6e}  drift = {drift_rel:+.3e}  "
                  f"max |u_surf - u_bot| = {shear_max:.3f} m/s "
                  f"at cell {k_shear} (x={cx_mesh[k_shear]:.1f}, "
                  f"y={cy_mesh[k_shear]:.1f}), h = {h[k_shear]:.2f} m")
            if k_lim >= 0:
                print(f"      dt_floor = {dt_floor:.4e}  ratio = {ratio:.3f}  "
                      f"limiter cell {k_lim} at "
                      f"(x={cx_mesh[k_lim]:.1f}, y={cy_mesh[k_lim]:.1f}), "
                      f"r = {radii[k_lim]:.3f} m, h = {h[k_lim]:.2f} m")
            wall_window = wall_now
            steps_window = 0
        if sim.time >= t_next_print or sim.time >= args.t_end:
            t_next_print += args.output_interval

    wall_total = time.perf_counter() - wall_start
    print(f"\nDone: {sim.step_count} steps, final dt = {sim.dt:.4e} s, "
          f"wall = {wall_total:.2f} s "
          f"({sim.step_count / wall_total:.0f} steps/s, "
          f"{args.t_end / wall_total:.1f}x real-time)")
    print(f"Mass error: {sim.mass_error:.3e} m^3")

    # ---- Plot ----
    if args.plot:
        plot_results(sim, args)


def plot_results(sim, args):
    """Save diagnostic plots of the plume."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri
    except ImportError:
        print("matplotlib not available; skipping plots")
        return

    cx = sim.centroid_x
    cy = sim.centroid_y
    S = sim.salinity
    h_lay = sim.h_layer

    nz = sim.nz_ml
    # Unstructured: k = 0 (Python) is bottom, k = -1 is surface.
    S_bot = S[:, 0]
    S_surf = S[:, -1]

    # Use a triangulation of centroids for plotting
    tri_plot = mtri.Triangulation(cx, cy)

    # ---- Surface salinity ----
    fig, ax = plt.subplots(1, 1, figsize=(14, 5))
    tcf = ax.tripcolor(tri_plot, S_surf, cmap="RdYlBu_r", vmin=0, vmax=35)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"Surface salinity at t = {sim.time:.0f} s")
    ax.set_aspect("equal")
    fig.colorbar(tcf, ax=ax, label="Salinity (PSU)")
    plt.tight_layout()
    plt.savefig("validation_examples/coastal/river_plume/surface_salinity_unstr.png", dpi=150)
    print("Saved validation_examples/coastal/river_plume/surface_salinity_unstr.png")
    plt.close()

    # ---- Bottom salinity ----
    fig, ax = plt.subplots(1, 1, figsize=(14, 5))
    tcf = ax.tripcolor(tri_plot, S_bot, cmap="RdYlBu_r", vmin=0, vmax=35)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"Bottom salinity at t = {sim.time:.0f} s")
    ax.set_aspect("equal")
    fig.colorbar(tcf, ax=ax, label="Salinity (PSU)")
    plt.tight_layout()
    plt.savefig("validation_examples/coastal/river_plume/bottom_salinity_unstr.png", dpi=150)
    print("Saved validation_examples/coastal/river_plume/bottom_salinity_unstr.png")
    plt.close()

    # ---- Along-channel cross-section ----
    # Sample cells near the centreline (y ~ basin_width/2)
    y_mid = args.basin_width / 2.0
    band = args.channel_width / 4.0
    mask = np.abs(cy - y_mid) < band
    if mask.sum() > 10:
        cx_band = cx[mask]
        order = np.argsort(cx_band)
        cx_sorted = cx_band[order]

        # Reverse the layer axis so layer 0 = top (matches the plot's
        # "build z down from surface" convention). The unstructured
        # solver stores k = 0 as bottom, k = nz-1 as top.
        h_band = h_lay[mask][order][:, ::-1]   # (n_band, nz), top-first
        S_band = S[mask][order][:, ::-1]
        bed_band = sim.bed[mask][order]

        fig, ax = plt.subplots(1, 1, figsize=(14, 5))

        # Build z interfaces
        n_band = cx_sorted.size
        z_iface = np.zeros((n_band, nz + 1))
        z_iface[:, 0] = bed_band + h_band.sum(axis=1)  # surface
        for k in range(nz):
            z_iface[:, k + 1] = z_iface[:, k] - h_band[:, k]

        # Plot each layer as colored strips
        sm = plt.cm.ScalarMappable(cmap="RdYlBu_r",
                                   norm=plt.Normalize(0, 35))
        for k in range(nz):
            for j in range(n_band - 1):
                color = sm.to_rgba(S_band[j, k])
                ax.fill_between(
                    [cx_sorted[j], cx_sorted[j + 1]],
                    [z_iface[j, k], z_iface[j + 1, k]],
                    [z_iface[j, k + 1], z_iface[j + 1, k + 1]],
                    color=color, linewidth=0,
                )

        # Bottom
        ax.fill_between(cx_sorted, bed_band,
                        bed_band.min() - 1, color="0.4")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("z (m)")
        ax.set_title(f"Along-channel salinity cross-section at t = {sim.time:.0f} s")
        ax.set_ylim(bed_band.min() - 1, 2)
        fig.colorbar(sm, ax=ax, label="Salinity (PSU)")
        plt.tight_layout()
        plt.savefig("validation_examples/coastal/river_plume/xsection_salinity_unstr.png", dpi=150)
        print("Saved validation_examples/coastal/river_plume/xsection_salinity_unstr.png")
        plt.close()

    # ---- Surface velocity ----
    fig, ax = plt.subplots(1, 1, figsize=(14, 5))
    hu = sim.hu_layer
    hv = sim.hv_layer
    # Unstructured: surface is k = -1, not k = 0.
    u_surf = np.where(h_lay[:, -1] > 1e-6, hu[:, -1] / h_lay[:, -1], 0.0)
    v_surf = np.where(h_lay[:, -1] > 1e-6, hv[:, -1] / h_lay[:, -1], 0.0)
    speed = np.sqrt(u_surf**2 + v_surf**2)

    tcf = ax.tripcolor(tri_plot, speed, cmap="magma", vmin=0, vmax=1.0)
    # Quiver subsample
    skip = max(1, sim.ncells // 500)
    idx = np.arange(0, sim.ncells, skip)
    ax.quiver(cx[idx], cy[idx], u_surf[idx], v_surf[idx],
              color="white", scale=10, width=0.002, alpha=0.7)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"Surface velocity at t = {sim.time:.0f} s")
    ax.set_aspect("equal")
    fig.colorbar(tcf, ax=ax, label="Speed (m/s)")
    plt.tight_layout()
    plt.savefig("validation_examples/coastal/river_plume/surface_velocity_unstr.png", dpi=150)
    print("Saved validation_examples/coastal/river_plume/surface_velocity_unstr.png")
    plt.close()

    print("\nAll plots saved.")


if __name__ == "__main__":
    main()
