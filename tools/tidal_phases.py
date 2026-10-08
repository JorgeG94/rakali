#!/usr/bin/env python3
"""Compute tidal constituent amplitudes and phases for a given epoch.

Uses EOT20 harmonic constants (amplitude, Greenwich phase lag) and computes
the astronomical arguments V0+u and nodal factors f at the requested epoch.

Model convention:  eta(t) = sum_n  A_n * cos(omega_n * t  -  phi_n)
where t is seconds since the epoch and phi is a phase LAG.

    A_n   = f_n * H_n
    phi_n = G_n - V0_n(t0) - u_n

References:
    Schureman (1958), Manual of Harmonic Analysis and Prediction of Tides
    Foreman (1977), Manual of Tidal Heights Analysis and Prediction
"""

import numpy as np
from datetime import datetime

# --- Astronomical arguments ---

def julian_century(dt):
    """Julian centuries since J2000.0 (2000-01-01 12:00 UTC)."""
    jd = (dt - datetime(2000, 1, 1, 12, 0, 0)).total_seconds() / 86400.0
    return jd / 36525.0

def gmst_degrees(dt):
    """Greenwich Mean Sidereal Time in degrees at datetime dt.

    Uses the IAU formula for GMST from Meeus (1991).
    """
    jd = (dt - datetime(2000, 1, 1, 12, 0, 0)).total_seconds() / 86400.0
    # GMST in hours
    gmst_h = 18.697374558 + 24.06570982441908 * jd
    return (gmst_h % 24) / 24.0 * 360.0

def astro_angles(T):
    """Mean astronomical longitudes (degrees) at Julian century T.

    Returns s, h, p, N (all in degrees).
      s = mean longitude of the Moon
      h = mean longitude of the Sun
      p = longitude of lunar perigee
      N = longitude of ascending lunar node
    """
    s = 218.3164477 + 481267.88123421 * T - 0.0015786 * T**2
    h = 280.46646 + 36000.76983 * T + 0.0003032 * T**2
    p = 83.3532465 + 4069.0137287 * T - 0.0103200 * T**2
    N = 125.0445479 - 1934.1362891 * T + 0.0020754 * T**2
    return s % 360, h % 360, p % 360, N % 360

def equilibrium_args(s, h, p, N, gmst):
    """Equilibrium argument V0 (degrees) for each constituent.

    Uses Schureman (1958) / Foreman (1977) formulas expressed in terms of
    GMST.  Since GMST = T + h  (T = mean solar hour angle), every formula
    below can be verified by substituting T = GMST - h and checking rates
    against the known constituent speeds.

    Phase constants (±90°) arise from the equilibrium-tide expansion
    (diurnal constituents pick up ±90° from the sin 2φ latitude factor).

    Constituent     Schureman V0            In terms of GMST
    ───────────     ────────────            ────────────────
    M2              2T + 2h − 2s            2·GMST − 2s
    S2              2T                      2·GMST − 2h
    K1              T + h − 90              GMST − 90
    O1              T + h − 2s − 90         GMST − 2s − 90
    N2              2T + 2h − 3s + p        2·GMST − 3s + p
    K2              2T + 2h                 2·GMST
    P1              T − h − 90              GMST − 2h − 90
    Q1              T + h − 3s + p − 90     GMST − 3s + p − 90
    """
    V0 = {
        'M2': 2*gmst - 2*s,
        'S2': 2*gmst - 2*h,
        'K1': gmst - 90.0,
        'O1': gmst - 2*s - 90.0,
        'N2': 2*gmst - 3*s + p,
        'K2': 2*gmst,
        'P1': gmst - 2*h - 90.0,
        'Q1': gmst - 3*s + p - 90.0,
    }
    return {k: v % 360 for k, v in V0.items()}

def nodal_corrections(N_deg):
    """Nodal factor f and nodal angle u (degrees) for each constituent.

    Simplified Schureman formulas using longitude of lunar node N.
    """
    N = np.radians(N_deg)
    cosN = np.cos(N)
    sinN = np.sin(N)
    cos2N = np.cos(2*N)
    sin2N = np.sin(2*N)

    f = {}
    u = {}

    # M2:  f = 1 - 0.0373*cosN;  u = -2.14*sinN  (degrees)
    # More precise: use Schureman Table 2
    f['M2'] = 1.0 - 0.03731*cosN + 0.00052*cos2N
    u['M2'] = np.degrees(-0.03731*sinN + 0.00052*sin2N)  # atan approx

    # Actually, standard Schureman:
    # f(M2) = (1 - 0.04742 cosN)^0.5  ≈ 1 - 0.02371 cosN  (linearised)
    # Let me use the more standard ones:

    # Foreman (1977) / Pawlowicz et al. nodal corrections:
    # M2: f = 1.0004 - 0.0373*cos(N) + 0.0002*cos(2N)
    #     u = -2.14*sin(N) (degrees)
    f['M2'] = 1.0004 - 0.0373*cosN + 0.0002*cos2N
    u['M2'] = -2.14*sinN  # degrees (approx from Schureman)

    # S2: purely solar, no nodal correction
    f['S2'] = 1.0
    u['S2'] = 0.0

    # K1: f = 1.006 + 0.115*cosN - 0.009*cos(2N)
    #     u = -8.86*sinN + 0.68*sin(2N)  (degrees)
    f['K1'] = 1.006 + 0.1150*cosN - 0.0088*cos2N
    u['K1'] = np.degrees(-0.1554*sinN + 0.0029*sin2N)

    # O1: f = 1.009 + 0.187*cosN - 0.015*cos(2N)
    #     u = 10.80*sinN - 1.34*sin(2N)  (degrees)
    f['O1'] = 1.0089 + 0.1871*cosN - 0.0147*cos2N
    u['O1'] = np.degrees(0.1885*sinN - 0.0234*sin2N)

    # N2: same nodal corrections as M2
    f['N2'] = f['M2']
    u['N2'] = u['M2']

    # K2: f = 1.024 + 0.286*cosN + 0.008*cos(2N)
    #     u = -17.74*sinN + 0.68*sin(2N)  (degrees)
    f['K2'] = 1.0241 + 0.2863*cosN + 0.0083*cos2N
    u['K2'] = np.degrees(-0.3093*sinN + 0.0029*sin2N)

    # P1: purely solar, no nodal correction
    f['P1'] = 1.0
    u['P1'] = 0.0

    # Q1: same nodal corrections as O1
    f['Q1'] = f['O1']
    u['Q1'] = u['O1']

    return f, u


def compute_model_params(epoch_str, eot20_amp, eot20_phase_deg):
    """Compute model amplitudes and phases for a given epoch.

    Parameters
    ----------
    epoch_str : str
        Epoch date string, e.g. "2025-03-26 00:00:00"
    eot20_amp : dict
        EOT20 amplitudes (m) keyed by constituent name
    eot20_phase_deg : dict
        EOT20 Greenwich phase lags (degrees) keyed by constituent name

    Returns
    -------
    dict with keys: amp, phase_rad, omega_rad_s for each constituent
    """
    dt = datetime.strptime(epoch_str, "%Y-%m-%d %H:%M:%S")
    T = julian_century(dt)
    s, h, p, N = astro_angles(T)
    gmst = gmst_degrees(dt)

    print(f"Epoch: {epoch_str}")
    print(f"Julian century T = {T:.10f}")
    print(f"s = {s:.4f}°  h = {h:.4f}°  p = {p:.4f}°  N = {N:.4f}°")
    print(f"GMST = {gmst:.4f}°")
    print()

    V0 = equilibrium_args(s, h, p, N, gmst)
    f, u_deg = nodal_corrections(N)

    # Standard angular frequencies (rad/s)
    # Derived from constituent periods in solar hours
    omega = {
        'M2': 1.405189e-4,   # 12.4206 h
        'S2': 1.454441e-4,   # 12.0000 h
        'K1': 7.292117e-5,   # 23.9345 h
        'O1': 6.759774e-5,   # 25.8193 h
        'N2': 1.378797e-4,   # 12.6583 h
        'K2': 1.458423e-4,   # 11.9672 h
        'P1': 7.252295e-5,   # 24.0659 h
        'Q1': 6.495854e-5,   # 26.8684 h
    }

    constituents = list(eot20_amp.keys())

    print(f"{'Const':>5}  {'H_eot(m)':>9}  {'G_eot(°)':>9}  {'V0(°)':>8}  "
          f"{'u(°)':>7}  {'f':>6}  {'A_model(m)':>10}  {'φ_model(rad)':>12}  "
          f"{'ω(rad/s)':>12}")
    print("-" * 105)

    results = {}
    for c in constituents:
        H = eot20_amp[c]
        G = eot20_phase_deg[c]
        A = f[c] * H
        # Model uses cos(omega*t - phi), so phi is a phase LAG:
        #   phi = G - V0 - u  (all in degrees), then convert to radians
        # This ensures cos(omega*t - phi) = cos(omega*t + V0 + u - G)
        # which is the standard tidal convention.
        phi_deg = (G - V0[c] - u_deg[c]) % 360
        phi_rad = np.radians(phi_deg)
        # Normalise to [-pi, pi]
        if phi_rad > np.pi:
            phi_rad -= 2*np.pi

        print(f"{c:>5}  {H:9.6f}  {G:9.4f}  {V0[c]:8.4f}  "
              f"{u_deg[c]:7.4f}  {f[c]:6.4f}  {A:10.6f}  {phi_rad:12.6f}  "
              f"{omega[c]:12.6e}")

        results[c] = {'amp': A, 'phase_rad': phi_rad, 'omega': omega[c]}

    return results


if __name__ == "__main__":
    import sys
    # TICON/GESLA harmonic constants for Fort Denison, Sydney Harbour
    # Station: (-33.8500, 151.2333), UHSLC record 1965-2012 (48 years)
    # Source: Piccioni et al. (2018), PANGAEA doi:10.1594/PANGAEA.896587
    # Amplitudes in cm in TICON -> converted to metres here
    # Phases in degrees (TICON convention: [-180, 180] -> converted to [0, 360])
    ticon_amp_cm = {
        'M2': 50.616,
        'S2': 12.528,
        'K1': 14.887,
        'O1':  9.686,
        'N2': 11.366,
        'K2':  3.742,
        'P1':  4.425,
        'Q1':  2.302,
    }
    ticon_phase_deg = {
        'M2': -52.547,
        'S2': -38.786,
        'K1': -30.749,
        'O1': -59.662,
        'N2': -60.654,
        'K2': -50.287,
        'P1': -33.615,
        'Q1': -81.869,
    }

    # Convert to metres and positive Greenwich phase lag [0, 360]
    amp = {k: v / 100.0 for k, v in ticon_amp_cm.items()}
    phase_deg = {k: v % 360 for k, v in ticon_phase_deg.items()}

    print("TICON/GESLA harmonic constants for Fort Denison (-33.85, 151.23)")
    print("=" * 50)
    for c in amp:
        print(f"  {c:>3}:  H = {amp[c]:.4f} m,  G = {phase_deg[c]:.3f}°")
    print()

    epoch = sys.argv[1] if len(sys.argv) > 1 else "2025-03-26 00:00:00"
    results = compute_model_params(epoch, amp, phase_deg)

    print()
    print("=== Namelist snippet (5 major constituents) ===")
    print()
    constituents = ['M2', 'S2', 'K1', 'O1', 'N2']
    amps = ", ".join(f"{results[c]['amp']:.6f}" for c in constituents)
    phases = ", ".join(f"{results[c]['phase_rad']:.6f}" for c in constituents)
    omegas = ", ".join(f"{results[c]['omega']:.6e}" for c in constituents)
    print(f"  ! Tidal forcing: TICON/GESLA Fort Denison (-33.85, 151.23)")
    print(f"  ! Constituents: {', '.join(constituents)}")
    print(f"  ! Epoch: {epoch}")
    print(f"  n_tidal_constituents = {len(constituents)}")
    print(f"  tidal_amp   = {amps}")
    print(f"  tidal_phase = {phases}")
    print(f"  tidal_omega = {omegas}")

    print()
    print("=== Namelist snippet (8 constituents) ===")
    print()
    constituents8 = ['M2', 'S2', 'K1', 'O1', 'N2', 'K2', 'P1', 'Q1']
    amps8 = ", ".join(f"{results[c]['amp']:.6f}" for c in constituents8)
    phases8 = ", ".join(f"{results[c]['phase_rad']:.6f}" for c in constituents8)
    omegas8 = ", ".join(f"{results[c]['omega']:.6e}" for c in constituents8)
    print(f"  ! Tidal forcing: TICON/GESLA Fort Denison (-33.85, 151.23)")
    print(f"  ! Constituents: {', '.join(constituents8)}")
    print(f"  ! Epoch: {epoch}")
    print(f"  n_tidal_constituents = {len(constituents8)}")
    print(f"  tidal_amp   = {amps8}")
    print(f"  tidal_phase = {phases8}")
    print(f"  tidal_omega = {omegas8}")
