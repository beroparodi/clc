"""Predictions: C/2014 UN271 scenarios (incl. a K2-like case), C/2025 Y1 and C/2024 T5."""
import os, sys
import numpy as np
from datetime import datetime, timedelta
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from un271_model import (eta, epoch_values, EP1, EP2, slope, m_nucleus, mag_from_afrho, add_mags)
from un271_orbit import ephem
KG = 0.01720209895

def r_parabolic(q, tp, t):
    """Heliocentric distance on a parabola (Barker's equation), q in au, datetimes."""
    dt = (t - tp).total_seconds() / 86400.0
    W = 3 * KG * dt / np.sqrt(2 * q**3)          # D + D^3/3 = W/3 ... solve D^3 + 3D - W = 0
    Y = np.cbrt(W / 2 + np.sqrt(W**2 / 4 + 1)); D = Y - 1 / Y
    return q * (1 + D**2)

def elems_ephem(t, q, e, inc, Om, w, tp):
    """Near-parabolic orbit (e ~ 1) via Barker's equation, plus the Earth from un271_orbit.
    Adequate for |e - 1| < 0.01 within a few years of perihelion."""
    from un271_orbit import earth_helio, jd
    D2R = np.pi / 180
    dt = (t - tp).total_seconds() / 86400.0
    W = 3 * KG * dt / np.sqrt(2 * q**3)
    Y = np.cbrt(W / 2 + np.sqrt(W**2 / 4 + 1)); Dn = Y - 1 / Y
    nu = 2 * np.arctan(Dn); r = q * (1 + Dn**2)
    i, O, om = inc * D2R, Om * D2R, w * D2R; u = nu + om
    C = r * np.array([np.cos(O)*np.cos(u) - np.sin(O)*np.sin(u)*np.cos(i),
                      np.sin(O)*np.cos(u) + np.cos(O)*np.sin(u)*np.cos(i),
                      np.sin(u)*np.sin(i)])
    E = earth_helio(jd(t)); G = C - E; Dl = np.linalg.norm(G); R = np.linalg.norm(E)
    eps = 23.4392911 * D2R
    ze = G[1]*np.sin(eps) + G[2]*np.cos(eps)
    return dict(r=r, D=Dl, dec=np.degrees(np.arcsin(ze / Dl)),
                phase=np.degrees(np.arccos((r*r + Dl*Dl - R*R) / (2*r*Dl))),
                elong=np.degrees(np.arccos((R*R + Dl*Dl - r*r) / (2*R*Dl))))

if __name__ == "__main__":
    e1, e2 = epoch_values(EP1), epoch_values(EP2); k = slope(e1, e2)
    A0, r0 = e2["afrho"], e2["r"]; et0 = eta(r0, "CO2")
    scen = {"K2-like (k=0.48)": lambda r: A0 * (r / r0)**-0.48,
            "CO, constant area": lambda r: A0 * (r / r0)**-1.5,
            "CO2 dayside, constant area": lambda r: A0 * eta(r, "CO2") / et0 * (r / r0)**-1.5,
            f"power law (k={k:.1f})": lambda r: A0 * (r / r0)**-k}
    print("== UN271 predictions (r' mag, 6\" aperture) ==")
    dates = [datetime(2026, 10, 1), datetime(2027, 7, 1), datetime(2028, 7, 1), datetime(2029, 7, 1),
             datetime(2031, 1, 25), datetime(2031, 4, 5)]
    for d in dates:
        g = ephem(d); r, D, a = g["r"], g["D"], g["phase"]; mN = m_nucleus(r, D, a)
        s = " & ".join(f"{add_mags(mN, mag_from_afrho(f(r), r, D, 6.0)):.2f}" for f in scen.values())
        print(d.date(), f"{r:.2f} {D:.2f} {mN:.2f} &", s)
    # quarterly rows for the K2-like line in the figure
    t = datetime(2021, 7, 1); out = []
    while t <= datetime(2032, 7, 1):
        g = ephem(t); r, D, a = g["r"], g["D"], g["phase"]; mN = m_nucleus(r, D, a)
        out.append(f"({t.year + (t.month - 1) / 12 - 2020:.3f},{add_mags(mN, mag_from_afrho(scen['K2-like (k=0.48)'](r), r, D, 6.0)):.3f})")
        t = (t.replace(day=1) + timedelta(days=95)).replace(day=1)
    print("%k2like\n" + " ".join(out))

    print("\n== C/2025 Y1 (ATLAS): q = 2.6 au, T = 2027 Oct 15 (parabolic) ==")
    tp = datetime(2027, 10, 15)
    for d in [datetime(2025, 12, 26), datetime(2026, 2, 15), datetime(2026, 11, 10), datetime(2027, 3, 1),
              datetime(2027, 6, 1), datetime(2027, 8, 1), datetime(2027, 10, 15)]:
        print(d.date(), f"r = {r_parabolic(2.6, tp, d):.2f} au")
    for rt in (6, 5, 4, 3):
        lo, hi = datetime(2024, 1, 1), tp
        for _ in range(60):
            mid = lo + (hi - lo) / 2
            if r_parabolic(2.6, tp, mid) > rt: lo = mid
            else: hi = mid
        print(f"  r = {rt} au on {lo.date()}")
    print("  Afrho growth factor r -> q=2.6 au (constant area, v ~ r^-1/2):")
    for rs in (6, 5, 4, 3):
        fac = lambda ice, chi: eta(2.6, ice, chi) / eta(rs, ice, chi) * (rs / 2.6)**1.5
        print(f"   from {rs} au: CO {fac('CO',1):.2f}  CO2 sub {fac('CO2',1):.2f}  H2O sub {fac('H2O',1):.1f}  H2O iso {fac('H2O',0.25):.3g}")
    for r in (6, 5, 4, 3, 2.6):
        d = 0.01
        kk = lambda ice, chi: 1.5 - (np.log(eta(r*(1+d), ice, chi)) - np.log(eta(r*(1-d), ice, chi))) / (np.log(1+d) - np.log(1-d))
        print(f"   local k at {r} au: CO {kk('CO',1):.2f} CO2 {kk('CO2',1):.2f} H2O sub {kk('H2O',1):.2f} H2O iso {kk('H2O',0.25):.1f}")

    # one-line formula, Eq. (14): m = H_N + 5 log(r D) + beta a - 2.5 log[1 + mu eta r^-1.5 D^j 10^(0.4 beta a)]
    from un271_model import H_V, VR, rR, BETA_N, M_SUN_R, AU_cm, ARCSEC
    HN = H_V - VR + rR; K = A0 / et0 * r0**1.5
    mu = K * 6.0 * ARCSEC / (4 * AU_cm) * 10**(0.4 * (HN - M_SUN_R))
    print(f"\n== Eq. (14), CO2 at the subsolar point: mu = {mu:.0f} (j = 1) ==")
    for d in [datetime(2026, 10, 1), datetime(2028, 7, 1), datetime(2031, 1, 25)]:
        g = ephem(d); r, D, a = g["r"], g["D"], g["phase"]
        m14 = HN + 5 * np.log10(r * D) + BETA_N * a - 2.5 * np.log10(1 + mu * eta(r, "CO2") * r**-1.5 * D * 10**(0.4 * BETA_N * a))
        print(d.date(), f"Eq. (14): {m14:.2f}   full model: {add_mags(m_nucleus(r, D, a), mag_from_afrho(scen['CO2 dayside, constant area'](r), r, D, 6.0)):.2f}")

    # pre-discovery test (Bernardinelli et al. 2021: Afrho ~1 m -> ~150 m between 28 and 20 au)
    g = eta(20, "CO2") / eta(28, "CO2") * (20 / 28)**-1.5
    print(f"\n== Pre-discovery test: eta_CO2(28) = {eta(28,'CO2'):.4f}, eta_CO2(20) = {eta(20,'CO2'):.3f}; "
          f"Afrho growth 28->20 au = {g:.0f} (k = {np.log(g)/np.log(28/20):.1f}); observed ~150 (k = {np.log(150)/np.log(1.4):.1f})")
    print(f"   model Afrho (normalized Feb 2024): {A0*eta(28,'CO2')/et0*(28/r0)**-1.5:.0f} cm at 28 au, "
          f"{A0*eta(20,'CO2')/et0*(20/r0)**-1.5:.0f} cm at 20 au")

    print("\n== C/2024 T5 (ATLAS) ==")
    el = (3.84251371, 1.00240368, 52.38681511, 100.67985066075, 352.41723883275,
          datetime(2027, 5, 6, 2, 38))
    for d in [datetime(2025, 11, 2), datetime(2026, 11, 10), datetime(2027, 2, 1), datetime(2027, 5, 6),
              datetime(2027, 8, 1), datetime(2027, 11, 1)]:
        g = elems_ephem(d, *el)
        print(d.date(), f"r={g['r']:.2f} D={g['D']:.2f} dec={g['dec']:.1f} elong={g['elong']:.0f} phase={g['phase']:.1f}")
    rn = elems_ephem(datetime(2026, 11, 10), *el)["r"]
    for ice, chi in (("CO", 1), ("CO2", 1), ("H2O", 1), ("H2O", 0.25)):
        f = eta(3.8425, ice, chi) / max(eta(rn, ice, chi), 1e-300) * (rn / 3.8425)**1.5
        print(f"  Afrho factor now ({rn:.2f} au) -> q: {ice} chi={chi}: {f:.3g}")
