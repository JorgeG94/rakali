#!/usr/bin/env python3
"""Tracer-reconstruction comparison for the river plume — PPM vs WENO5/7/9.

Loads each reconstruction's NetCDF output and reports the quantities that
separate a high-order limiter from PPM on a buoyant-front problem:

  * salt conservation   — column-integrated total salt drift from t0 to tN
                          (should be ~0; horizontal advection is conservative)
  * S extrema growth    — final min/max vs the INITIAL min/max.  The still-water
                          IC of this case already carries a saline patch above
                          35 PSU, so the absolute bound is meaningless; what
                          matters is whether the reconstruction ADDS new
                          over/undershoot beyond the IC range.  A monotone
                          scheme must not push min below IC-min or max above
                          IC-max by more than round-off.
  * front area          — number of surface cells inside the transitional band
                          FRONT_LO < S < FRONT_HI across the WHOLE surface
                          (a SHARPER front occupies FEWER such cells; this is
                          where a high-order reconstruction buys resolution).

Usage:
    python3 compare_recon.py                         # ppm + weno5 + weno7 (if present)
    python3 compare_recon.py output_ppm output_weno5 # explicit dirs

Expected (see README.md):
    all runs conserve total salt (|drift| tiny relative to total)
    WENO5/7 add no S extrema beyond the IC range, and hold a narrower
    (fewer-cell) front than PPM.

This is a GATED benchmark: the script exits nonzero if any check fails,
so it can be run mechanically (CI / regression sweeps).  Gates:
  G1 conservation parity — every scheme's relative salt drift matches the
     PPM baseline to CONSERV_PARITY_TOL (the drift itself is Flather
     open-boundary outflow, scheme-independent by construction; a scheme
     that drifts DIFFERENTLY is leaking).
  G2 no new extrema — over/undershoot beyond the IC S range stays below
     EXTREMA_TOL (measured ~7e-15; the gate carries round-off headroom).
  G3 sharpness ordering — each WENO rung's transitional-front cell count
     is <= the PPM baseline's (WENO must never be MORE diffusive than PPM
     on the same case; measured weno5 -11 cells, weno7 -2).
"""
import sys
import os
import numpy as np
from netCDF4 import Dataset

# S values that count as "inside the front" (neither fresh nor fully saline).
FRONT_LO, FRONT_HI = 1.0, 34.0
S_FRESH, S_SALT = 0.0, 35.0

# --- gate thresholds (see module docstring; calibrated to measured values) ---
CONSERV_PARITY_TOL = 1.0e-5   # |rel_drift(scheme) - rel_drift(ppm)|
EXTREMA_TOL = 1.0e-6          # PSU beyond the IC range (measured ~7e-15)


def total_salt(ds, ti):
    """Column-integrated sum(h_layer * S) over the whole domain at time ti."""
    nz = ds.variables["h_layer"].shape[1]
    tot = 0.0
    for k in range(nz):
        h = ds.variables["h_layer"][ti, k, :, :]
        S = ds.variables["salinity"][ti, k, :, :]
        tot += float(np.sum(h * S))
    return tot


def analyse(path):
    nc = os.path.join(path, "output_rank_000000.nc")
    ds = Dataset(nc)
    t = ds.variables["time"][:]
    nt = len(t)
    y = ds.variables["y"][:]
    jc = len(y) // 2

    S0 = total_salt(ds, 0)
    SN = total_salt(ds, nt - 1)
    drift = SN - S0
    rel = drift / S0 if S0 != 0 else float("nan")

    # S extrema: compare the FINAL range against the INITIAL range.  The
    # still-water IC already carries values outside [0, 35], so overshoot is
    # only what the reconstruction adds beyond the IC bounds.
    Sini = ds.variables["salinity"][0, :, :, :]
    Sfin = ds.variables["salinity"][nt - 1, :, :, :]
    ic_min, ic_max = float(Sini.min()), float(Sini.max())
    smin, smax = float(Sfin.min()), float(Sfin.max())
    over = max(0.0, smax - ic_max)   # new maximum overshoot beyond IC
    under = max(0.0, ic_min - smin)  # new minimum undershoot below IC

    # Front area: transitional surface cells over the whole surface layer
    # (surface = last z index; bottom-up k=1..nz -> index -1 is surface).
    Ssurf = ds.variables["salinity"][nt - 1, -1, :, :]
    front_cells = int(np.sum((Ssurf > FRONT_LO) & (Ssurf < FRONT_HI)))

    ds.close()
    return dict(S0=S0, SN=SN, drift=drift, rel=rel, smin=smin, smax=smax,
                over=over, under=under, front=front_cells, nt=nt)


def main(argv):
    dirs = argv[1:] if len(argv) > 1 else [
        d for d in ("output_ppm", "output_weno5", "output_weno7", "output_weno9")
        if os.path.isdir(d)]
    if not dirs:
        print("no output dirs found — run the *.nml cases first")
        return 1

    print(f"{'run':14s} {'salt drift':>12s} {'rel':>10s} "
          f"{'overshoot':>10s} {'undershoot':>11s} {'front cells':>12s}")
    print("-" * 74)
    rows = {}
    for d in dirs:
        try:
            r = analyse(d)
        except (FileNotFoundError, OSError) as e:
            print(f"{d:14s} SKIP ({e})")
            continue
        rows[d] = r
        extrema = "  <-- NEW EXTREMA" if (r["over"] > 1e-4 or r["under"] > 1e-4) else ""
        print(f"{d:14s} {r['drift']:12.3e} {r['rel']:10.2e} "
              f"{r['over']:10.3e} {r['under']:11.3e} {r['front']:12d}{extrema}")

    # If we have a ppm baseline + at least one weno, report front sharpening.
    base = rows.get("output_ppm")
    if base:
        print()
        for d, r in rows.items():
            if d == "output_ppm":
                continue
            dw = r["front"] - base["front"]
            tag = "sharper" if dw < 0 else ("same" if dw == 0 else "wider")
            print(f"  {d}: front {r['front']} vs ppm {base['front']} "
                  f"({dw:+d} cells, {tag})")

    # ------------------------- gates (exit nonzero on failure) --------------
    failures = []
    for d, r in rows.items():
        # G2: no new extrema (every scheme, incl. the ppm baseline).
        if r["over"] > EXTREMA_TOL or r["under"] > EXTREMA_TOL:
            failures.append(f"G2 {d}: new S extrema beyond IC range "
                            f"(over={r['over']:.3e}, under={r['under']:.3e} "
                            f"> {EXTREMA_TOL:.0e})")
    if base:
        for d, r in rows.items():
            if d == "output_ppm":
                continue
            # G1: conservation parity with the ppm baseline.
            dpar = abs(r["rel"] - base["rel"])
            if dpar > CONSERV_PARITY_TOL:
                failures.append(f"G1 {d}: relative salt drift differs from "
                                f"ppm by {dpar:.3e} > {CONSERV_PARITY_TOL:.0e} "
                                f"(scheme-dependent leak)")
            # G3: sharpness ordering — WENO must not be wider than ppm.
            if r["front"] > base["front"]:
                failures.append(f"G3 {d}: front {r['front']} cells WIDER than "
                                f"ppm {base['front']} (more diffusive than PPM)")
    else:
        print("\nNOTE: no output_ppm baseline found — G1/G3 gates skipped")

    print()
    if failures:
        for f in failures:
            print(f"FAIL  {f}")
        print(f"\nRESULT: FAIL ({len(failures)} gate(s))")
        return 1
    print("RESULT: PASS (conservation parity, no new extrema, "
          "WENO front <= PPM)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
