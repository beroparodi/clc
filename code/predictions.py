"""Predictions: C/2014 UN271 scenarios (incl. a K2-like case with the effective k = 0.49 of C/2017 K2,
backtest.py B3), C/2025 Y1 and C/2024 T5 (JPL orbits, as in forecast_test.py)."""
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

if __name__ == "__main__":
    e1, e2 = epoch_values(EP1), epoch_values(EP2); k = slope(e1, e2)
    A0, r0 = e2["afrho"], e2["r"]; et0 = eta(r0, "CO2")
    scen = {"K2-like (k=0.49)": lambda r: A0 * (r / r0)**-0.49,
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
        out.append(f"({t.year + (t.month - 1) / 12 - 2020:.3f},{add_mags(mN, mag_from_afrho(scen['K2-like (k=0.49)'](r), r, D, 6.0)):.3f})")
        t = (t.replace(day=1) + timedelta(days=95)).replace(day=1)
    print("%k2like\n" + " ".join(out))

    # robustness of the CO2 inference (Sec. IV A): subsolar point vs sunlit hemisphere of a slow
    # rotator, constant v, and the thermal parameters A_B, eps of Eq. (8)
    import un271_model as UM
    from un271_model import eta_hemi
    r1, r2 = e1["r"], e2["r"]
    kof = lambda f: 1.5 + np.log(f(r2) / f(r1)) / np.log(r1 / r2)
    grow = lambda f: f(20) / f(28) * (20 / 28)**-1.5
    print("\n== UN271: robustness of the CO2 inference (k between the epochs; growth 28->20 au) ==")
    for lab, f in (("CO2 subsolar", lambda r: eta(r, "CO2")), ("CO2 sunlit hemisphere", lambda r: eta_hemi(r, "CO2"))):
        print(f"  {lab:22s}: k = {kof(f):.2f} (v const: {kof(f)+0.5:.2f}); growth 28->20 au = {grow(f):.0f}")
    for AB, EPS in ((0.04, 0.9), (0.10, 0.9), (0.01, 1.0)):
        UM.A_B, UM.EPS = AB, EPS
        f = lambda r: eta(r, "CO2")
        print(f"  CO2 subsolar, A_B={AB}, eps={EPS}: k = {kof(f):.2f}; growth 28->20 au = {grow(f):.0f}")
    UM.A_B, UM.EPS = 0.01, 0.9
    eh0 = eta_hemi(r0, "CO2")
    hemi = lambda r: A0 * eta_hemi(r, "CO2") / eh0 * (r / r0)**-1.5
    for d in (datetime(2026, 10, 1), datetime(2031, 1, 25)):
        g = ephem(d); r, D, a = g["r"], g["D"], g["phase"]; mN = m_nucleus(r, D, a)
        print(f"  {d.date()}: r' = {add_mags(mN, mag_from_afrho(scen['CO2 dayside, constant area'](r), r, D, 6.0)):.2f} (subsolar), "
              f"{add_mags(mN, mag_from_afrho(hemi(r), r, D, 6.0)):.2f} (hemisphere)")
    print(f"  absolute Afrho normalized to Feb 2024, subsolar: {scen['CO2 dayside, constant area'](28)/100:.2f} m at 28 au, "
          f"{scen['CO2 dayside, constant area'](20)/100:.0f} m at 20 au (observed ~1 m and ~150 m)")
    print(f"  absolute Afrho, hemisphere: {hemi(28)/100:.2f} m at 28 au, {hemi(20)/100:.0f} m at 20 au")
    # assumption (4): if gas drag sets the grain speed, v ~ Q^(1/2), Afrho ~ Q^(1/2) and k = -p/2 = (k + 1/2)/2
    for lab, f in (("CO2 subsolar", lambda r: eta(r, "CO2")), ("CO2 sunlit hemisphere", lambda r: eta_hemi(r, "CO2")),
                   ("CO2 isothermal", lambda r: eta(r, "CO2", 0.25)), ("CO subsolar", lambda r: eta(r, "CO"))):
        print(f"  gas-drag grain speed (Afrho ~ Q^1/2), {lab:22s}: k = {(kof(f) + 0.5) / 2:.2f}")
    # pre-discovery photometry: coma-to-nucleus flux ratio of the subsolar CO2 model on the way in, in the
    # 6" radius of the paper and in the ~1" of the DES photometry of Bernardinelli et al. (2021)
    print("  pre-discovery: coma/nucleus flux ratio, subsolar CO2 normalized to Feb 2024 (6\" | 1\" radius):")
    def date_at(rt):
        lo, hi = datetime(2008, 1, 1), datetime(2024, 1, 1)
        for _ in range(60):
            mid = lo + (hi - lo) / 2
            if ephem(mid)["r"] > rt: lo = mid
            else: hi = mid
        return lo
    for rt in (34.1, 30.0, 28.0, 27.4, 27.0, 26.0, 25.0, 24.0, 23.9, 23.7):
        g = ephem(date_at(rt)); r, D, a = g["r"], g["D"], g["phase"]; FN = 10**(-0.4 * m_nucleus(r, D, a))
        f6, f1 = [10**(-0.4 * mag_from_afrho(scen['CO2 dayside, constant area'](r), r, D, th)) / FN for th in (6.0, 1.0)]
        print(f"   r = {rt:5.1f} au ({date_at(rt).date()}): {f6:.3f} | {f1:.3f}")
    t = datetime(2021, 7, 1); out = []
    while t <= datetime(2032, 7, 1):
        g = ephem(t); r, D, a = g["r"], g["D"], g["phase"]; mN = m_nucleus(r, D, a)
        out.append(f"({t.year + (t.month - 1) / 12 - 2020:.3f},{add_mags(mN, mag_from_afrho(hemi(r), r, D, 6.0)):.3f})")
        t = (t.replace(day=1) + timedelta(days=95)).replace(day=1)
    print("%co2hemi\n" + " ".join(out))

    # C/2025 Y1 and C/2024 T5: JPL orbits (forecast_test.py), Kepler/Barker/hyperbolic as in the note on the orbit
    import forecast_test as FT
    def r_at(key, t):
        return FT.anomaly(FT.COMETS[key], FT.jd(t.isoformat()))[0]
    def when(key, rt):
        lo, hi = datetime(2024, 1, 1), FT.datetime(2000, 1, 1) + timedelta(days=FT.COMETS[key]["Tp"] - 2451544.5)
        for _ in range(60):
            mid = lo + (hi - lo) / 2
            if r_at(key, mid) > rt: lo = mid
            else: hi = mid
        return lo
    fac = lambda ice, chi, ra, rb: eta(rb, ice, chi) / max(eta(ra, ice, chi), 1e-300) * (ra / rb)**1.5
    for key in ("y1", "t5"):
        c = FT.COMETS[key]; tp = datetime(2000, 1, 1) + timedelta(days=c["Tp"] - 2451544.5)
        print(f"\n== {c['name']}: q = {c['q']:.2f} au, e = {c['e']:.5f}, perihelion {tp:%Y-%m-%d %H:%M} (TT) ==")
        for rt in (6, 5, 4.14, 4, 3):
            if rt > c["q"]: print(f"  r = {rt} au on {when(key, rt).date()}")
        t0 = datetime(2026, 11, 10); print(f"  r on {t0.date()}: {r_at(key, t0):.2f} au")
        r_start = 4.0 if key == "y1" else r_at(key, t0)
        print(f"  Afrho growth factor {r_start:.2f} au -> q (constant area, v ~ r^-1/2): " + ", ".join(
            f"{lab} {fac(x, ch, r_start, c['q']):.3g}" for lab, (x, ch) in
            (("CO", ("CO", 1)), ("CO2 sub", ("CO2", 1)), ("CO2 iso", ("CO2", 0.25)), ("H2O sub", ("H2O", 1)), ("H2O iso", ("H2O", 0.25)))))
        # Barker's parabolic equation (note on the orbit) as a check of the near-parabolic approximation
        dmax = max(abs(r_at(key, t) - r_parabolic(c["q"], tp, t)) for t in
                   (datetime(2025, 6, 1) + timedelta(days=10 * i) for i in range(120)) if r_at(key, t) < 6)
        print(f"  Barker's parabola with the same q and T: |r - r_exact| <= {dmax:.4f} au inside 6 au")

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
