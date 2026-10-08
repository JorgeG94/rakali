#!/usr/bin/env python3
"""Run Castro dam break benchmark.

Usage:
    python run_sweep.py --grid 4096                     # single resolution
    python run_sweep.py --grid 4096 --profile           # single + nsys profile
    python run_sweep.py --sweep                         # full sweep G1-G8
    python run_sweep.py --sweep --binary ./build/rakali # custom binary path
"""

import argparse
import subprocess
import sys
import re
import tempfile
import os
from datetime import datetime
from pathlib import Path

# Grid resolutions from Castro plan: G1-G8
SWEEP_SIZES = [100, 200, 400, 800, 1600, 2000, 4096]
#SWEEP_SIZES = [100, 200, 400]

DOMAIN_LENGTH = 10.0  # [-5, 5] in each direction

NAMELIST_TEMPLATE = """\
&grid_nml
  nx = {nx}
  ny = {ny}
  dx = {dx}
  dy = {dy}
  nghost = 2
/

&time_nml
  t_end = 0.2
  cfl = 0.9
  dt_max = 1.0
/

&physics_nml
  manning_n = 0.0
  dry_threshold = 1.0e-6
  wind_stress_x = 0.0
  wind_stress_y = 0.0
/

&output_nml
  output_to_file = .false.
  output_dir = "./output"
  output_interval = 0.0
/

&boundary_nml
  bc_west  = "wall"
  bc_east  = "wall"
  bc_south = "wall"
  bc_north = "wall"
/

&initial_condition_nml
  init_condition = "castro_dam_break"
/

&logging_nml
  log_level = "info"
  status_interval = 0.0
/
"""


def make_namelist(n):
    """Write a temporary namelist for an NxN grid and return its path."""
    dx = DOMAIN_LENGTH / n
    content = NAMELIST_TEMPLATE.format(nx=n, ny=n, dx=dx, dy=dx)
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".nml", delete=False)
    f.write(content)
    f.close()
    return f.name


def parse_output(output):
    """Extract metrics from tidal stdout/stderr."""
    steps = None
    solver_time = None
    mcells = None
    mass_change = None

    for line in output.splitlines():
        m = re.search(r"Total steps:\s+(\d+)", line)
        if m:
            steps = int(m.group(1))
        m = re.search(r"Solver time:\s+([\d.]+)", line)
        if m:
            solver_time = float(m.group(1))
        m = re.search(r"Mcells/s:\s+([\d.]+)", line)
        if m:
            mcells = float(m.group(1))
        m = re.search(r"Mass change:\s+([\d.eE+-]+)", line)
        if m:
            mass_change = float(m.group(1))

    return steps, solver_time, mcells, mass_change


def run_case(binary, n, profile=False, profile_line=None):
    """Run a single NxN grid resolution and return parsed results."""
    nml_path = make_namelist(n)
    dx = DOMAIN_LENGTH / n

    cmd = [binary, nml_path]
    if profile_line:
        # User-provided profiler command, e.g. 'ncu -k "my_kernel" -c 1 --set=full'
        import shlex
        cmd = shlex.split(profile_line) + cmd
    elif profile:
        report_name = f"castro_{n}x{n}"
        cmd = [
            "nsys", "profile",
            "--stats=true",
            "-o", report_name,
            "--force-overwrite=true",
        ] + cmd

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,
        )
        output = result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT for {n}x{n}")
        return None
    except FileNotFoundError as e:
        print(f"  Command not found: {e}")
        return None
    finally:
        os.unlink(nml_path)

    if profile:
        print(output)

    steps, solver_time, mcells, mass_change = parse_output(output)

    if steps is None or solver_time is None:
        print(f"  Failed to parse output for {n}x{n}")
        print(output[-500:] if len(output) > 500 else output)
        return None

    return {
        "nx": n,
        "ny": n,
        "cells": n * n,
        "dx": dx,
        "steps": steps,
        "wall_s": solver_time,
        "mcells_s": mcells,
        "mass_change": mass_change,
        "cell_updates_s": n * n * steps / solver_time if solver_time > 0 else 0,
    }


def run_single(binary, n, profile=False, profile_line=None):
    """Run a single grid size and print results."""
    now = datetime.now().strftime("%H:%M:%S")
    ncells = n * n
    dx = DOMAIN_LENGTH / n
    tag = ""
    if profile_line:
        tag = f" [{profile_line.split()[0]}]"
    elif profile:
        tag = " [nsys profile]"
    print(f"[{now}] Running {n}x{n} ({ncells/1e6:.1f}M cells, dx={dx:.6f}m)"
          f"{tag}", flush=True)

    r = run_case(binary, n, profile=profile, profile_line=profile_line)
    if r:
        print(
            f"  {r['steps']} steps, {r['wall_s']:.2f}s, "
            f"{r['mcells_s']:.1f} Mcells/s, "
            f"mass err = {r['mass_change']:.2e}"
        )
    return r


def run_sweep(binary):
    """Run the full resolution sweep and print a summary table."""
    print(f"Binary: {binary}")
    print(f"Domain: [{-DOMAIN_LENGTH/2}, {DOMAIN_LENGTH/2}]^2")
    print(f"t_end = 1.0 s, CFL = 0.9, wall BCs, no friction")
    print()

    results = []
    for n in SWEEP_SIZES:
        r = run_single(binary, n)
        if r:
            results.append(r)
        print()

    if not results:
        print("No results collected.")
        return

    print("=" * 80)
    print("Castro Dam Break Benchmark — Summary")
    print("=" * 80)
    print(
        f"{'Grid':>12s}  {'Cells':>12s}  {'dx (m)':>10s}  "
        f"{'Steps':>7s}  {'Wall (s)':>10s}  {'Mcells/s':>10s}  "
        f"{'Mass err':>12s}"
    )
    print("-" * 80)
    for r in results:
        print(
            f"{r['nx']}x{r['ny']:>5d}  {r['cells']:>12,d}  {r['dx']:>10.6f}  "
            f"{r['steps']:>7d}  {r['wall_s']:>10.2f}  {r['mcells_s']:>10.1f}  "
            f"{r['mass_change']:>12.2e}"
        )
    print("=" * 80)

    # CSV output
    csv_path = Path(__file__).parent / "castro_results.csv"
    with open(csv_path, "w") as f:
        f.write("nx,ny,cells,dx,steps,wall_s,mcells_s,mass_change\n")
        for r in results:
            f.write(
                f"{r['nx']},{r['ny']},{r['cells']},{r['dx']},"
                f"{r['steps']},{r['wall_s']},{r['mcells_s']},{r['mass_change']}\n"
            )
    print(f"\nResults written to {csv_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Castro dam break benchmark runner"
    )
    parser.add_argument(
        "--binary", default="../../build/rakali",
        help="Path to tidal binary (default: ../../build/rakali)"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--grid", type=int, metavar="N",
        help="Run a single NxN grid (e.g. --grid 4096)"
    )
    group.add_argument(
        "--sweep", action="store_true",
        help="Run full sweep over all resolutions (100-8192)"
    )
    parser.add_argument(
        "--profile", action="store_true",
        help="Run under nsys profile (only with --grid)"
    )
    parser.add_argument(
        "--profile-line", type=str, default=None, metavar="CMD",
        help='Custom profiler command prefix, e.g. \'ncu -k "my_kernel" -c 1 --set=full\' (only with --grid)'
    )
    args = parser.parse_args()

    if (args.profile or args.profile_line) and args.sweep:
        parser.error("--profile/--profile-line is only compatible with --grid, not --sweep")

    binary = str(Path(args.binary).resolve())
    if not Path(binary).exists():
        print(f"Error: binary not found: {binary}")
        sys.exit(1)

    if args.grid:
        run_single(binary, args.grid, profile=args.profile,
                   profile_line=args.profile_line)
    else:
        run_sweep(binary)


if __name__ == "__main__":
    main()
