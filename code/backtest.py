"""Backtest of the sublimation-based light-curve model against archival data.
B1: SOHO/SWAN water-production exponents (Combi et al. 2019) vs constant-area envelope.
B2: active fractions recomputed from Q1 and r_N vs Combi's values (implementation check).
B3: C/2017 K2 (PANSTARRS): distant dust, CO and water ratios, long-range hindcast.
B4: C/2020 F3 (NEOWISE): exponents, perihelion hindcast, energy budget.
B5: 3I/ATLAS: exponents and energy budget.
B6: median brightening of long-period comets (Lacerda et al. 2025).
B7: local activity index vs the empirical model of Holt et al. (2026).  Needs numpy, scipy."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from un271_model import eta, eta_geom, latent, A_B, S0, AMU, afrho_from_mag, balance_mantle, M_SUN_V
from swan_data import LONG_PERIOD, SHORT_PERIOD, LATER

M_H2O = 18.015 * AMU; L_W = latent("H2O")
F1 = (1 - A_B) * S0                       # absorbed flux at 1 au [W m^-2]

def lnQ(r, chi, ice="H2O"):               # constant-area production, up to a constant
    return np.log(eta(r, ice, chi)) - 2 * np.log(r)

def envelope(q, rmax=None):
    """Mean power-law exponent of constant-area water production between q and rmax,
    for the subsolar (chi=1) and isothermal (chi=1/4) limits."""
    rmax = rmax or max(2.0, q + 1.0)
    p = lambda chi: (lnQ(rmax, chi) - lnQ(q, chi)) / (np.log(rmax) - np.log(q))
    return p(0.25), p(1.0)                # (steep, flat)

H_MANTLE = (np.inf, 100, 30, 10, 3, 1, 0.3, 0.1, 0.03)   # conductances h [W m^-2 K^-1]

def envelope_mantle(q, rmax=None):
    """Range of the constant-area exponent when the ice may lie under a dust mantle of any
    conductance h >= 0.03 W m^-2 K^-1, in both limits (note on the mantle)."""
    rmax = rmax or max(2.0, q + 1.0)
    lz = lambda r, c, h: np.log(balance_mantle(r, "H2O", c, h)[2])
    ps = [(lz(rmax, c, h) - lz(q, c, h)) / np.log(rmax / q) for h in H_MANTLE for c in (0.25, 1.0)]
    return min(ps), max(ps)

def judge(pobs, err, lo, hi, nsig=2.0):
    if pobs is None: return None
    return (pobs + nsig * err >= lo) and (pobs - nsig * err <= hi)

rows = []
for name, q, pp, ep, pq, eq in LONG_PERIOD:
    lo, hi = envelope(q)
    rows.append((name, "LP", q, pp, ep, pq, eq, lo, hi))
for name, q, Qa, pp, ep, Qb, pq, eq, rN, AFa, AFb in SHORT_PERIOD:
    lo, hi = envelope(q)
    rows.append((name, "SP", q, pp, ep, pq, eq, lo, hi))
for name, q, pp, ep, pq, eq in LATER:                # after the survey (3I: interstellar)
    lo, hi = envelope(q)
    rows.append((name, "LT", q, pp, ep, pq, eq, lo, hi))

def tally(sel):
    n = ok = 0; flat = steep = 0
    for name, cls, q, pp, ep, pq, eq, lo, hi in sel:
        for p, e in ((pp, ep), (pq, eq)):
            j = judge(p, e, lo, hi)
            if j is None: continue
            n += 1; ok += j
            if not j: flat += (p > hi); steep += (p < lo)
    return n, ok, flat, steep

BINS = ((0, .3), (.3, .6), (.6, .9), (.9, 1.1), (1.1, 1.4), (1.4, 1.9))

def binned_medians(seed=1, ndraw=10000, sel=None):
    """Median exponent in bins of q with a 68% bootstrap interval (squares in Fig. 3);
    sel: optional set of classes ("LP", "SP", "LT") to restrict the population."""
    P = np.array([(r[2], p) for r in rows if sel is None or r[1] in sel for p in (r[3], r[5]) if p is not None])
    rng = np.random.default_rng(seed); out = []
    for a, b in BINS:
        sel_ = P[(P[:, 0] >= a) & (P[:, 0] < b)]
        if len(sel_) == 0: out.append((a, b, 0, np.nan, np.nan, np.nan, np.nan)); continue
        sel = sel_
        boot = np.median(rng.choice(sel[:, 1], (ndraw, len(sel))), axis=1)
        out.append((a, b, len(sel), np.exp(np.log(sel[:, 0]).mean()), np.median(sel[:, 1]),
                    np.percentile(boot, 16), np.percentile(boot, 84)))
    return out

def fig2_coordinates(ymin=-10.0):
    """pgfplots coordinates of Fig. 3 (p-q diagram): band edges, mantle edge, exponents, binned medians."""
    qs = np.concatenate([np.geomspace(0.005, 0.1, 6), np.linspace(0.15, 1.85, 35)])
    lines = {"sub": [], "iso": [], "mantle": []}
    for q in qs:
        steep, flat = envelope(q)
        lines["iso"].append((q, steep)); lines["sub"].append((q, flat))
        lines["mantle"].append((q, envelope_mantle(q)[1]))
    for k, v in lines.items():
        print(f"%{k}\n" + " ".join(f"({q:.4f},{p:.3f})" for q, p in v))
    off = 0
    for pop, cls in (("lp", ("LP", "LT")), ("sp", ("SP",))):
        for k, (ip, ie) in (("pre", (3, 4)), ("post", (5, 6))):
            pts = []
            for r in rows:
                if r[1] not in cls or r[ip] is None: continue
                if r[ip] < ymin: off += 1; continue
                pts.append(f"({r[2]:.4f},{r[ip]:.2f}) +- (0,{r[ie]:.2f})")
            print(f"%{k}_{pop}\n" + " ".join(pts))
    print("%medians\n" + " ".join(f"({qg:.3f},{m:.2f}) += (0,{b84 - m:.2f}) -= (0,{m - b16:.2f})"
                                   for a, b, n, qg, m, b16, b84 in binned_medians()))
    print(f"% off scale (p < {ymin}): {off}")

def neff(r, x, c, h=1e-3):
    """Local activity index of a dominant coma, n = 7/2 - d ln(eta)/d ln(r) (note on Holt et al.);
    c = 1 (subsolar), 0.25 (isothermal) or "hemi" (sunlit hemisphere)."""
    return 3.5 - (np.log(eta_geom(r*(1+h), x, c)) - np.log(eta_geom(r*(1-h), x, c))) / (np.log(1+h) - np.log(1-h))

def fign_coordinates(nmax=18.0):
    """pgfplots coordinates of Fig. 2 (activity index): n(r) per ice and limit, and Holt et al.'s fit."""
    rs = np.geomspace(1.0, 35.0, 70)
    for lab, (x, c) in {"co": ("CO", 1), "co2sub": ("CO2", 1), "co2hemi": ("CO2", "hemi"), "co2iso": ("CO2", 0.25),
                        "h2osub": ("H2O", 1), "h2oiso": ("H2O", 0.25)}.items():
        pts = [(r, neff(r, x, c)) for r in rs]
        pts = [(r, n) for r, n in pts if n <= nmax]
        print(f"%{lab}\n" + " ".join(f"({r:.3f},{n:.2f})" for r, n in pts))
    print("%holt\n" + " ".join(f"({r:.2f},{0.40*r*(np.log(r)+1)-1.78:.2f})" for r in np.linspace(2.5, 10.5, 17)))

if __name__ == "__main__":
    if sys.argv[1:] == ["fig2"]:
        fig2_coordinates(); sys.exit()
    if sys.argv[1:] == ["fign"]:
        fign_coordinates(); sys.exit()
    print("== B1: SWAN exponents vs constant-area envelope (2 sigma) ==")
    for lab, sel in [("long-period", [r for r in rows if r[1] == "LP"]),
                     ("short-period", [r for r in rows if r[1] == "SP"]),
                     ("after 2016", [r for r in rows if r[1] == "LT"]),
                     ("survey only", [r for r in rows if r[1] != "LT"]), ("all", rows)]:
        n, ok, fl, st = tally(sel)
        print(f"{lab:13s}: {ok}/{n} pass ({100*ok/n:.0f}%); failures flatter {fl}, steeper {st}")
    for q in (0.1, 0.3, 0.5, 1.0, 1.3, 1.5):
        lo, hi = envelope(q); print(f"  envelope q={q}: [{lo:.2f}, {hi:.2f}]")
    # the same test when the ice may lie under a dust mantle (note on the mantle)
    mrows = [r[:7] + envelope_mantle(r[2]) for r in rows]
    n, ok, fl, st = tally(mrows)
    print(f"with a dust mantle (h >= 0.03 W m^-2 K^-1): {ok}/{n} pass ({100*ok/n:.0f}%); "
          f"failures flatter {fl}, steeper {st}")
    for q in (0.1, 0.5, 1.0):
        print(f"  mantle envelope q={q}: [{envelope_mantle(q)[0]:.2f}, {envelope_mantle(q)[1]:.2f}]",
              " thick mantle (h=0.03) Ts, Ti at 1 au, subsolar:",
              [round(v) for v in balance_mantle(1.0, "H2O", 1.0, 0.03)[:2]] if q == 1.0 else "")
    # robust statistics: medians
    allp = [p for r in rows for p in (r[3], r[5]) if p is not None]
    print("  median observed p:", np.median(allp), " N =", len(allp))
    for k, lab in ((3, "before"), (5, "after")):
        print(f"  median {lab} perihelion:", np.median([r[k] for r in rows if r[k] is not None]))
    # medians in bins of q (open squares in the figure), 68% bootstrap intervals
    print("  q bin      N  q_geo  median  boot 16-84%      band")
    for a, b, n, qg, med, b16, b84 in binned_medians():
        lo, hi = envelope(qg)
        print(f"  {a:.1f}-{b:.1f}  {n:3d}  {qg:.3f}  {med:6.2f}  [{b16:.2f}, {b84:.2f}]  [{lo:.2f}, {hi:.2f}]")
    # null model: a fixed exponent p = -2 (Q ~ r^-2) with the same 2-sigma criterion
    nfix = sum(judge(p, e, -2.0, -2.0) for r in rows for p, e in ((r[3], r[4]), (r[5], r[6]))
               if p is not None)
    print(f"  fixed p = -2: {nfix}/{len(allp)} pass")
    # null model 2: a fixed band of the model's median width, centred at c (no dependence on q)
    obs = [(r[2], p, e, r[7], r[8]) for r in rows for p, e in ((r[3], r[4]), (r[5], r[6])) if p is not None]
    w = round(float(np.median([r[8] - r[7] for r in rows])), 1)   # 0.8; edges rounded like the data
    for c in (-2.2, -2.3, -2.4, -2.5, -2.6):
        nb = sum(bool(judge(p, e, round(c - w/2, 2), round(c + w/2, 2))) for q, p, e, lo, hi in obs)
        print(f"  fixed band [{c - w/2:.2f}, {c + w/2:.2f}] (width {w:.2f}): {nb}/{len(obs)} pass ({100*nb/len(obs):.0f}%)")
    for lab, sel in (("q < 0.9", lambda q: q < 0.9), ("0.9-1.1", lambda q: 0.9 <= q < 1.1), ("q >= 1.1", lambda q: q >= 1.1)):
        s = [o for o in obs if sel(o[0])]
        nm = sum(bool(judge(p, e, lo, hi)) for q, p, e, lo, hi in s)
        nb = sum(bool(judge(p, e, -2.6, -1.8)) for q, p, e, lo, hi in s)
        print(f"  {lab:8s}: N = {len(s)}, model {nm}, fixed band [-2.60, -1.80] {nb}")
    # sensitivity to the outer distance r2 of the mean model exponent
    for lab, f in (("max(1.5, q+0.5)", lambda q: max(1.5, q + 0.5)), ("max(2, q+1)", lambda q: max(2.0, q + 1.0)),
                   ("max(3, q+1)", lambda q: max(3.0, q + 1.0))):
        nr = sum(bool(judge(p, e, *envelope(q, f(q)))) for q, p, e, lo, hi in obs)
        print(f"  r2 = {lab:15s}: {nr}/{len(obs)} pass ({100*nr/len(obs):.0f}%)")

    # the trend with q by population: short-period (SP) vs long-period and later comets (LP, LT)
    from scipy.stats import mannwhitneyu
    for lab, cl in (("long-period + later", ("LP", "LT")), ("long-period, survey", ("LP",)), ("short-period", ("SP",))):
        print(f"  medians by q bin, {lab}: " + "; ".join(f"{a:.1f}-{b:.1f}: {m:.2f} (N={n})"
              for a, b, n, qg, m, b16, b84 in binned_medians(sel=set(cl)) if n))
    for lab, cl in (("all", ("LP", "SP", "LT")), ("long-period + later", ("LP", "LT")), ("survey only", ("LP", "SP"))):
        inner = [p for r in rows if r[1] in cl and r[2] < 0.9 for p in (r[3], r[5]) if p is not None]
        outer = [p for r in rows if r[1] in cl and r[2] >= 1.1 for p in (r[3], r[5]) if p is not None]
        print(f"  Mann-Whitney q<0.9 vs q>=1.1 ({lab}): N = {len(inner)}, {len(outer)}; "
              f"medians {np.median(inner):.2f}, {np.median(outer):.2f}; p (one-sided) = "
              f"{mannwhitneyu(inner, outer, alternative='greater').pvalue:.3f}")
    # criterion (ii): active fraction of the projected area (subsolar, the most generous case) at q
    print("  criterion (ii), f = Q(q) / Q_max(q) with Q_max from Eq. (11), subsolar, f_X = 1:")
    nfail = set()
    for name, q, Qa, pp, ep, Qb, pq, eq, rN, AFa, AFb in SHORT_PERIOD:
        if rN is None: continue
        Qmax = np.pi * (rN * 1e3)**2 * F1 / q**2 * eta(q, "H2O") / L_W / M_H2O
        fs = [Q1 * q**p / Qmax for Q1, p in ((Qa, pp), (Qb, pq)) if Q1 is not None and p is not None]
        if fs and max(fs) > 1:
            nfail.add(name.split(" (")[0]); print(f"    {name:38s} rN = {rN} km: f = " + ", ".join(f"{f:.1f}" for f in fs))
    print(f"    comets with f > 1: {len(nfail)}: {', '.join(sorted(nfail))}")

    print("\n== B2: active fractions (isothermal, 4 pi rN^2) recomputed vs Combi ==")
    # Combi et al. normalize to a different reference sublimation rate, hence a constant ratio (5.13);
    # 46P differs by a further factor 4: Combi's table implies a radius of 1.2 km, we use 0.6 km.
    Ziso = eta(1.0, "H2O", 0.25) * F1 * 0.25 / L_W / M_H2O       # molecules m^-2 s^-1
    print(f"  Z_iso(1 au) = {Ziso/1e4:.3e} molecules cm^-2 s^-1")
    ratios = []
    for name, q, Qa, pp, ep, Qb, pq, eq, rN, AFa, AFb in SHORT_PERIOD:
        for Q, AF in ((Qa, AFa), (Qb, AFb)):
            if Q is None or AF is None or rN is None: continue
            ours = Q / (Ziso * 4 * np.pi * (rN * 1e3)**2)
            ratios.append(ours / AF)
    ratios = np.array(ratios)
    print(f"  N = {len(ratios)}, ours/Combi: median {np.median(ratios):.3f}, "
          f"range {ratios.min():.3f}-{ratios.max():.3f}, scatter {ratios.std()/ratios.mean()*100:.1f}%")

    print("\n== B3: C/2017 K2 (PANSTARRS) ==")
    # (a) distant dust: H in 40,000 km aperture, 2013 (CFHT) vs 2017 (HST), Jewitt et al. 2017
    r1, r2, dH, sH = 23.742, 15.874, 8.6 - 7.18, 0.5
    k = dH * 0.4 * np.log(10) / np.log(r1 / r2); sk = sH * 0.4 * np.log(10) / np.log(r1 / r2)
    print(f"  (a) Afrho ratio {10**(0.4*dH):.2f} -> k = {k:.2f} +- {sk:.2f}")
    for ice in ("CO", "CO2"):
        for chi in (1.0, 0.25):
            e1, e2 = eta(r1, ice, chi), eta(r2, ice, chi)
            kk = 1.5 + np.log(e2 / e1) / np.log(r1 / r2) if e1 > 1e-30 else np.inf
            print(f"      {ice} chi={chi}: k = {kk:.2f}")
    # (b) HST Afrho (V, 40,000 km) for the long-range hindcast
    D = 15.816; m = 19.33; rho_arcsec = 4e4 / (D * 1.495978707e8) * 206264.8
    from un271_model import afrho_from_mag as afm
    A_hst = afm(m, 15.874, D, rho_arcsec, msun=M_SUN_V)
    print(f"  (b) HST Afrho(15.87 au) = {A_hst:.0f} cm (V, rho = 40,000 km, alpha 3.7 deg)")
    for kk, lab in ((1.5, "CO constant area"), (k, "measured 2013-17 slope")):
        print(f"      hindcast Afrho(1.79 au) with k={kk:.2f} ({lab}): {A_hst*(15.874/1.79)**kk:.0f} cm;"
              f" observed A(0)frho ~15000 cm -> ratio {A_hst*(15.874/1.79)**kk/15000:.1f}")
    print(f"      effective k 15.87->1.79 au: {np.log(15000/A_hst)/np.log(15.874/1.79):.2f}")
    # (c) CO energy-limited test, CRIRES+ (Hmiddouch et al. 2025)
    ra, rb = 2.73, 2.12
    pred = np.exp(lnQ(rb, 1.0, "CO") - lnQ(ra, 1.0, "CO"))
    obs, sobs = 60 / 36, 60 / 36 * np.hypot(10 / 36, 10 / 60)
    print(f"  (c) CO ratio 2.73->2.12 au: predicted {pred:.2f}, observed {obs:.2f} +- {sobs:.2f}")
    # (d) water ratios
    for (rA, rB, QA, QB, sA, sB, lab) in [(3.23, 2.73, 3.96e27, 1.37e28, 0.79e27, 0.14e28, "OH (UVES)"),
                                          (2.73, 2.12, 1.37e28, 5.41e28, 0.14e28, 0.54e28, "OH (UVES)"),
                                          (2.73, 2.12, 335, 446, 96, 82, "H2O (CRIRES+)")]:
        ob = QB / QA; so = ob * np.hypot(sA / QA, sB / QB)
        ps = np.exp(lnQ(rB, 1.0) - lnQ(rA, 1.0)); pi = np.exp(lnQ(rB, 0.25) - lnQ(rA, 0.25))
        print(f"  (d) {lab} {rA}->{rB} au: observed x{ob:.2f}+-{so:.2f}; subsolar x{ps:.2f}, isothermal x{pi:.1f}")

    print("\n== B4: C/2020 F3 (NEOWISE) ==")
    q = 0.295
    lo, hi = envelope(q, 1.55); print(f"  envelope [{lo:.2f}, {hi:.2f}] vs pre -2.5+-0.2, post -3.5+-0.1")
    # hindcast: pre-perihelion law (r 0.53-1.30) extrapolated vs constant-area model, to 0.337 au
    Qpre = lambda r: 6.9e28 * r**-2.5
    for chi in (1.0, 0.25):
        Qm = Qpre(1.30) * np.exp(lnQ(0.337, chi) - lnQ(1.30, chi))
        print(f"  chi={chi}: model from 1.30 au -> Q(0.337) = {Qm:.2e}; observed peak 5.27e30 (ratio {5.27e30/Qm:.1f})")
    print(f"  pre-law at 0.337: {Qpre(0.337):.2e}")
    for RN in (2.5, 5.0):
        Qmax = np.pi * (RN * 1e3)**2 * F1 / 0.337**2 * eta(0.337, "H2O") / L_W / M_H2O
        print(f"  R_N = {RN} km: Qmax(0.337) = {Qmax:.2e} -> f = {5.27e30/Qmax:.2f}")
    print(f"  minimum R_N for f<=1: {2.5*np.sqrt(5.27e30/(np.pi*(2.5e3)**2*F1/0.337**2*eta(0.337,'H2O')/L_W/M_H2O)):.2f} km")

    print("\n== B5: 3I/ATLAS ==")
    lo, hi = envelope(1.36, 3.0); print(f"  pre envelope [1.36, 3.0]: [{lo:.2f}, {hi:.2f}] vs -5.9+-0.8")
    lo, hi = envelope(1.46, 2.24); print(f"  post envelope [1.46, 2.24]: [{lo:.2f}, {hi:.2f}] vs -3.3+-0.3")
    for RN in (2.8, 1.6, 1.5, 1.3, 1.1):   # HST: upper limit 2.8 km (2025); measured 1.3 +- 0.2 km (Hui et al. 2026)
        Qmax = np.pi * (RN * 1e3)**2 * F1 / 1.49**2 * eta(1.49, "H2O") / L_W / M_H2O
        print(f"  R_N = {RN} km: Qmax(1.49 au) = {Qmax:.2e} -> f for 4-5e28: {4e28/Qmax:.2f}-{5.15e28/Qmax:.2f}")

    print("\n== B6: median brightening of long-period comets from r to 1 au (Lacerda et al. 2025, Table 10) ==")
    obs = {"new": [1.68, 3.70, 5.32, 6.40, 6.73], "int": [3.25, 4.85, 5.80, 7.60, 8.20], "old": [3.90, 5.25, 6.22, 6.31, 12.70]}
    rs = [2, 3, 4, 5, 6]
    dm_j = lambda r, x, c: (np.log10(eta(1, x, c) / eta(r, x, c)) + 2 * np.log10(r)) / 0.2453       # Jorda: Q -> total magnitude
    dm_d = lambda r, x, c: 2.5 * np.log10(eta(1, x, c) / eta(r, x, c)) + 8.75 * np.log10(r)        # dust, fixed aperture size
    print("  r (au)            " + "".join(f"{r:>12}" for r in rs))
    for lab, (x, c) in {"CO subsolar": ("CO", 1), "CO2 subsolar": ("CO2", 1), "H2O subsolar": ("H2O", 1), "H2O isothermal": ("H2O", 0.25)}.items():
        print(f"  {lab:16s}  " + "".join(f"{min(dm_j(r,x,c),dm_d(r,x,c)):6.2f}-{max(dm_j(r,x,c),dm_d(r,x,c)):5.2f}" for r in rs))
    for g, v in obs.items():
        print(f"  observed, {g:5s}   " + "".join(f"{x:12.2f}" for x in v))

    print("\n== B7: local activity index n_eff = 7/2 - dln(eta)/dln(r) vs Holt et al. (2026) Hab fit ==")
    rs = [2.5, 3, 4, 5, 6, 8, 10, 12, 15]
    print("  r (au)          " + "".join(f"{r:>7}" for r in rs))
    print("  Holt (Fixed-ab) " + "".join(f"{0.40*r*(np.log(r)+1)-1.78:7.1f}" for r in rs))
    for lab, (x, c) in {"CO subsolar": ("CO", 1), "CO2 subsolar": ("CO2", 1), "CO2 hemisphere": ("CO2", "hemi"),
                        "CO2 isothermal": ("CO2", 0.25), "H2O subsolar": ("H2O", 1), "H2O isothermal": ("H2O", 0.25)}.items():
        print(f"  {lab:16s}" + "".join(f"{neff(r,x,c):7.2f}" for r in rs))
