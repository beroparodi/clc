"""Check of the C/2014 UN271 forecasts against COBS photometry (2024-2026).

COBS CCD/CMOS magnitudes are measured in large circular apertures (diameter D_ap, 8"-73").
Each is converted to the paper's r' magnitude in a 6" radius aperture:
  1. subtract the nucleus (V band, H_V = 8.62, beta = 0.04 mag/deg) from the total flux,
  2. scale the coma flux by 6"/rho_ap (steady coma, F_C ~ rho; Afrho constant),
  3. convert V -> r' with V - r' = V - R - (r' - R) = 0.40 - 0.18 = 0.22 (R: r' - R = 0.18),
  4. add the nucleus back in r'.
Unfiltered magnitudes calibrated to APASS V (method C, catalog AV) are treated as V.
Data: COBS Comet Observation Database, CC BY-NC-SA 4.0 (cobs/un271_cobs.csv)."""
import os, sys, csv
import numpy as np
from datetime import datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from un271_model import (eta, epoch_values, EP1, EP2, slope, m_nucleus, mag_from_afrho, add_mags, H_V, VR, rR)
from un271_orbit import ephem

V_MINUS_RP = VR - rR                 # V - r' = 0.22
HERE = os.path.dirname(os.path.abspath(__file__))

def load(path=os.path.join(HERE, "cobs", "un271_cobs.csv")):
    rows = []
    with open(path) as f:
        for row in csv.DictReader(line for line in f if not line.startswith("#")):
            rows.append(row)
    return rows

def to_rprime_6(row):
    """r' magnitude in a 6" radius aperture from one COBS row (None if not convertible)."""
    meth, cat = row["method"], row["refcat"]
    if meth == "V" or (meth == "C" and cat == "AV"):
        dcol, band = V_MINUS_RP, "V"
    elif meth == "k":
        dcol, band = -rR, "R"                       # r' = R + 0.18
    else:
        return None
    t = datetime.fromisoformat(row["date"]); g = ephem(t); r, D, a = g["r"], g["D"], g["phase"]
    m, rho = float(row["mag"]), float(row["ap_diam_arcsec"]) / 2
    mN_band = m_nucleus(r, D, a, vr=0.0 if band == "V" else VR, rr=0.0)   # nucleus in V or R
    fc = 10**(-0.4 * m) - 10**(-0.4 * mN_band)
    if fc <= 0: return None
    mc6 = -2.5 * np.log10(fc * 6.0 / rho) - dcol       # coma in r', 6" radius
    return t, add_mags(mc6, m_nucleus(r, D, a)), row["observer"], meth

def scenarios():
    e1, e2 = epoch_values(EP1), epoch_values(EP2); k = slope(e1, e2)
    A0, r0 = e2["afrho"], e2["r"]; et0 = eta(r0, "CO2")
    return {"K2-like": lambda r: A0 * (r / r0)**-0.48,
            "CO": lambda r: A0 * (r / r0)**-1.5,
            "CO2": lambda r: A0 * eta(r, "CO2") / et0 * (r / r0)**-1.5,
            "power law": lambda r: A0 * (r / r0)**-k}

def forecast(t, f):
    g = ephem(t); r, D, a = g["r"], g["D"], g["phase"]
    return add_mags(m_nucleus(r, D, a), mag_from_afrho(f(r), r, D, 6.0))

if __name__ == "__main__":
    pts = [p for p in (to_rprime_6(r) for r in load()) if p is not None]
    sc = scenarios()
    print(f"{len(pts)} convertible COBS points")
    # monthly medians of the converted magnitudes vs the forecasts
    from collections import defaultdict
    mon = defaultdict(list)
    for t, m, obs, meth in pts: mon[(t.year, t.month)].append(m)
    print("month     N  median r'(6\")   " + "  ".join(f"{k:>9s}" for k in sc))
    for (y, mo), v in sorted(mon.items()):
        if y < 2024: continue
        t = datetime(y, mo, 15)
        print(f"{y}-{mo:02d}  {len(v):3d}  {np.median(v):6.2f}         " +
              "  ".join(f"{forecast(t, f):9.2f}" for f in sc.values()))
    # residuals per scenario, 2024 Jun - 2026 Nov
    sel = [(t, m) for t, m, o, me in pts if t >= datetime(2024, 6, 1)]
    for name, f in sc.items():
        res = np.array([m - forecast(t, f) for t, m in sel])
        print(f"  {name:9s}: median O-C = {np.median(res):+.2f} mag, scatter {np.std(res):.2f} (N = {len(res)})")
    # filtered observers only (MAS01 V and R, Johnson-Cousins)
    print("  MAS01 (filtered V/R, 21-54\" apertures):")
    for t, m, o, me in pts:
        if o == "MAS01": print(f"   {t.date()} {me}  r'(6\") = {m:.2f}   CO2 {forecast(t, sc['CO2']):.2f}  CO {forecast(t, sc['CO']):.2f}")
