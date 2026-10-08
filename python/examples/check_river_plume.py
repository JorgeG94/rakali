#!/usr/bin/env python3
"""Headless sanity check using the canonical river_plume.nml.

Loads the well-tested ``validation_examples/coastal/river_plume/river_plume.nml``
and prints one diagnostic line per minute of sim time.  No matplotlib, no
display dependency — paste the console output to debug "is the plume
forming at all?" questions.

Run with::

    python python/examples/check_river_plume.py

Prerequisites:
  * ``validation_examples/coastal/river_plume/bathymetry.nc`` must exist.
    If not::

        cd validation_examples/coastal/river_plume && python generate_bathymetry.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "python"))
os.environ.setdefault("RAKALI_LIB",
                      os.path.join(REPO, "build_nvhpc", "librakali.so"))

from rakali import Rakali

NML = os.path.join(REPO, "validation_examples", "coastal", "river_plume", "river_plume.nml")
BATHY_NC = os.path.join(REPO, "validation_examples", "coastal", "river_plume", "bathymetry.nc")

if not os.path.isfile(NML):
    sys.exit(f"ERROR: namelist not found: {NML}")

if not os.path.isfile(BATHY_NC):
    sys.exit(
        f"ERROR: bathymetry file not found: {BATHY_NC}\n\n"
        "Generate it first with:\n"
        "    cd validation_examples/coastal/river_plume && python generate_bathymetry.py\n"
    )

print(f"Loading: {NML}")
sim = Rakali(NML)
print(f"Created: {sim}")
print(f"  nx = {sim.nx}, ny = {sim.ny}, nz = {sim.nz}, "
      f"dx = {sim.dx} m, dy = {sim.dy} m")
print()

# Initial state — already on host since we just created the sim
sim.sync_host()
S0 = sim.salinity
h0 = sim.h_phys
dry0 = h0 < 1.0e-6
print(f"Initial state:")
print(f"  h_phys [min..max] = [{h0.min():.3f}, {h0.max():.3f}] m")
print(f"  dry cells = {int(dry0.sum())} / {dry0.size}")
S0_wet = np.ma.masked_array(S0[:, :, 0], mask=dry0)
print(f"  surface S [wet min..max] = "
      f"[{float(S0_wet.min()):.3f}, {float(S0_wet.max()):.3f}] PSU")
print()

# Per-callback diagnostic
print("Diagnostics every 60 sim-seconds (wet cells only):")
print()
print("  t (min) | surf S [min .. max] | bot S [min .. max] | "
      "max |dS| | wet u [min..max] | argmin (i,j)")
print("  " + "-" * 105)


def report(s):
    h = s.h_phys
    dry = h < 1.0e-6
    S = s.salinity
    surf = np.ma.masked_array(S[:, :, 0], mask=dry)
    bot  = np.ma.masked_array(S[:, :, -1], mask=dry)
    contrast = np.ma.masked_array(S[:, :, 0] - S[:, :, -1], mask=dry)

    if surf.count() > 0:
        flat_argmin = int(surf.argmin(fill_value=np.inf))
        i_min, j_min = np.unravel_index(flat_argmin, surf.shape)
    else:
        i_min, j_min = -1, -1

    # u_layer at the surface (m/s) — Coriolis sanity check
    hu_lay = s._get_3d_field('rki_get_hu_layer_ptr')
    h_lay = s._get_3d_field('rki_get_h_layer_ptr')
    ng = s._nghost
    hu_phys = hu_lay[ng:-ng, ng:-ng, 0]
    h_lphys = h_lay[ng:-ng, ng:-ng, 0]
    with np.errstate(invalid='ignore', divide='ignore'):
        u_surf = np.where(h_lphys > 1e-6, hu_phys / h_lphys, 0.0)
    u_surf_m = np.ma.masked_array(u_surf, mask=dry)

    print(f"  {s.time/60:6.2f}  | "
          f"{float(surf.min()):6.2f} .. {float(surf.max()):6.2f}   | "
          f"{float(bot.min()):6.2f} .. {float(bot.max()):6.2f}  | "
          f"{float(np.abs(contrast).max()):7.3f} | "
          f"{float(u_surf_m.min()):+5.2f} .. {float(u_surf_m.max()):+5.2f} | "
          f"({i_min:3d}, {j_min:3d})")


# Run for 10 minutes of sim time, callback every 60 s
T_REPORT_END = min(600.0, sim.time + 600.0)
sim.advance(T_REPORT_END, callback=report, callback_interval=60.0)

print()
print(f"Final t = {sim.time:.1f} s   step_count = {sim.step_count}")

# Where is the freshest surface cell at the end?
# Layer index nz-1 (i.e. [..., -1]) is the surface under k=1=bed,
# k=nz=surface bottom-up convention. See CLAUDE.md vertical layer notes.
sim.sync_host()
S_end = sim.salinity[:, :, -1]
dry_end = sim.h_phys < 1.0e-6
S_end_wet = np.ma.masked_array(S_end, mask=dry_end)
print(f"Final wet surface S min = {float(S_end_wet.min()):.3f} PSU "
      f"(should be << {35.0:.0f} if the plume is working)")

sim.destroy()
