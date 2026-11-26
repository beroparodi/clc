"""Forecasts of Eq. (14) against the empirical law, Eq. (1), for 25 COBS light curves (Sec. III C).

The sample was assembled in three steps, each rule fixed before the forecasts of its comets were computed:
  * the ten comets of the interactive tool (FIRST);
  * the ten comets not yet in the sample with the most COBS observations in one apparition among comets
    with perihelion in 2010-2025 (SECOND; census of the COBS Observation List API, 2026 November 17);
  * the five comets not yet in the sample with the most COBS observations in one such apparition whose
    record starts beyond 5 au before perihelion, with at least three records there (THIRD).
All light curves except C/2014 UN271 are binned by one rule (bin_cobs.py: bins of max(3 d, span/50)).

For each comet the points are sorted in time and split in half.  On the first half, Eq. (14) is fitted
(one constant K) for each of the five cases (CO, CO2 and H2O at the subsolar point, CO2 and H2O
isothermal) and the case with the lowest chi^2 is kept; the empirical law m = m0 + 5 lg D + 2.5 n lg r
is fitted (m0 and n) by weighted least squares.  Both then forecast the second half (rms in mag).
An ice "switches on" within the fitted half when the kept case gives a local activity index
n = 7/2 - d ln(eta)/d ln(r) above 4.5 at some point of that half (Sec. II C); this flag was fixed
before the forecasts of the 25-comet sample were computed.
Two ices (added after the one-ice test; design fixed before its forecasts were computed): a volatile ice
V and water W share one heat spread and act on separate patches, Afrho = (K_V eta_V + K_W eta_W) r^-3/2,
with K_V, K_W >= 0; the pairs are CO + H2O (subsolar and isothermal water), CO2 + H2O subsolar and
CO2 + H2O isothermal; both constants and the pair are fitted on the first half (two constants, like
Eq. (1) with m0 and n).  Whole apparitions: one ice with one mu, and the same ice with one mu before and
one after perihelion (the change of the active area).
C/2014 UN271 as in Sec. IV A: CO2 at the subsolar point with K from the 2024 LCO point, the empirical
law through the two professional points, tested on the COBS monthly medians (r' in 6").
Total magnitudes of the other comets are compared in a fixed aperture of 20 000 km (j = 0) with
beta_N = 0.04 and the nucleus magnitudes below (15 = placeholder for an unmeasured nucleus).
Data: COBS Comet Observation Database, CC BY-NC-SA 4.0 (cobs/cobs_lightcurves.txt).  Needs numpy, scipy."""
import os, sys
import numpy as np
from datetime import datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from un271_model import eta

HERE = os.path.dirname(os.path.abspath(__file__))
AU_CM, ARCSEC, D2R, KG = 1.495978707e13, np.pi / 180 / 3600, np.pi / 180, 0.01720209895

def jd(s):
    return (datetime.fromisoformat(s) - datetime(2000, 1, 1, 12)).total_seconds() / 86400 + 2451545.0

# orbital elements (JPL Small-Body Database; Hale-Bopp: NASA NSSDC), H_N, ice of the tool's preset
COMETS = {
    "un271": dict(name="C/2014 UN271", q=10.94525581, e=0.99869477, i=95.55335892, Om=190.0846361196, w=326.33161880945, Tp=2462892.9068, HN=8.40),
    "hb":    dict(name="C/1995 O1 (Hale-Bopp)", q=0.9141, e=0.9951, i=89.43, Om=282.47, w=130.59, Tp=jd("1997-04-01T03:19"), HN=9.77),
    "t5":    dict(name="C/2024 T5 (ATLAS)", q=3.84251371, e=1.00240368, i=52.38681511, Om=100.67985066075, w=352.41723883275, Tp=jd("2027-05-06T02:38"), HN=15),
    "y1":    dict(name="C/2025 Y1 (ATLAS)", q=2.601924069, e=0.999552254, i=14.55698337, Om=166.5518278289, w=124.2978493120, Tp=2461692.5495, HN=15),
    "neo":   dict(name="C/2020 F3 (NEOWISE)", q=0.295, e=0.99921, i=128.94, Om=61.010, w=37.279, Tp=jd("2020-07-03T12:00"), HN=15),
    "p12":   dict(name="12P/Pons-Brooks", q=0.78086475, e=0.95456023, i=74.19083929, Om=255.85505059528, w=198.98800614657, Tp=2460421.6316, HN=13),
    "i3":    dict(name="3I/ATLAS", q=1.35638454, e=6.13941774, i=175.11310480, Om=322.15684249105, w=128.01051366675, Tp=2460977.9814, HN=17.04),
    "a3":    dict(name="C/2023 A3 (Tsuchinshan-ATLAS)", q=0.39145906, e=1.00012475, i=139.11502467, Om=21.558609111275, w=308.4872848828, Tp=2460581.2353, HN=15),
    "sw1":   dict(name="29P/Schwassmann-Wachmann", q=5.77512802, e=0.04479814, i=9.36477141, Om=312.39091654503, w=50.662700071797, Tp=2458589.9259, HN=10.22),
    "hyak":  dict(name="C/1996 B2 (Hyakutake)", q=0.23022353, e=0.99989165, i=124.92204939, Om=188.04513199216, w=130.1751209781, Tp=2450204.8941, HN=15.71),
    # the second ten (JPL Small-Body Database, osculating elements near perihelion; no measured nucleus: H_N = 15)
    "q2":  dict(name="C/2014 Q2 (Lovejoy)", q=1.290358238881392, e=0.9977738520356274, i=80.30130181438155, Om=94.97505909428631, w=12.39526184803142, Tp=2457052.569249606, HN=15),
    "k2":  dict(name="C/2017 K2 (PANSTARRS)", q=1.797281961376049, e=1.000616562373503, i=87.55487509681637, Om=88.24109127986824, w=236.1923489989567, Tp=2459933.227066954, HN=15),
    "t2":  dict(name="C/2017 T2 (PANSTARRS)", q=1.615095955888523, e=0.9996765980850288, i=57.23142385337199, Om=64.37877215685057, w=92.99331253710214, Tp=2458974.449996825, HN=15),
    "p1":  dict(name="C/2009 P1 (Garradd)", q=1.550535584223489, e=1.000997435576857, i=106.1772916768743, Om=325.9972285528648, w=90.74706541748093, Tp=2455919.177555371, HN=15),
    "e3":  dict(name="C/2022 E3 (ZTF)", q=1.11224437022534, e=1.000301905819192, i=109.169480756749, Om=302.5550197168474, w=145.81492879, Tp=2459957.285198830, HN=15),
    "l3":  dict(name="C/2019 L3 (ATLAS)", q=3.554492702308672, e=1.001405275505796, i=48.36374698279935, Om=290.7888867057587, w=171.6137051016775, Tp=2459589.132623246, HN=15),
    "l4":  dict(name="C/2011 L4 (PANSTARRS)", q=0.3015442238164768, e=1.000032767108972, i=84.20819702923397, Om=65.66588639150341, w=333.6516376633611, Tp=2456361.669960002, HN=15),
    "a6":  dict(name="C/2025 A6 (Lemmon)", q=0.5299041302912412, e=0.9955362293967411, i=143.663707214442, Om=108.0979873301086, w=132.9688404751593, Tp=2460988.039387102, HN=15),
    "v2j": dict(name="C/2015 V2 (Johnson)", q=1.636950053071025, e=1.001702614530096, i=49.87481167153636, Om=69.8521003063327, w=164.8980978916873, Tp=2457916.847458616, HN=15),
    "v2z": dict(name="C/2020 V2 (ZTF)", q=2.227959873700067, e=1.00089953037377, i=131.6110418896687, Om=212.3712351547342, w=162.4243172377264, Tp=2460073.056589952, HN=15),
    # the third group: first seen beyond 5 au (JPL Small-Body Database; H_N = 15)
    "us10": dict(name="C/2013 US10 (Catalina)", q=0.8229694514684749, e=1.000280139519764, i=148.8782967756639, Om=186.1445797203693, w=340.358754276709, Tp=2457342.221641829, HN=15),
    "k1":   dict(name="C/2012 K1 (PANSTARRS)", q=1.054597098294, e=1.000152915493971, i=142.4282581726916, Om=317.7376491689354, w=203.1066008825931, Tp=2456897.160363006, HN=15),
    "g2":   dict(name="C/2021 G2 (ATLAS)", q=4.982360288479666, e=1.000603552652205, i=48.46954508656729, Om=221.095294828915, w=343.2672327577748, Tp=2460562.634007078, HN=15),
    "e2":   dict(name="C/2022 E2 (ATLAS)", q=3.666295808424358, e=1.001060418168188, i=137.1316685887328, Om=125.3761151010601, w=41.72559864492467, Tp=2460567.633887190, HN=15),
    "n2":   dict(name="C/2022 N2 (PANSTARRS)", q=3.825529163602361, e=1.003236448782429, i=5.503022530288324, Om=319.7406957022052, w=75.38965400185295, Tp=2460888.270675301, HN=15),
}
FIRST = ["un271", "hb", "t5", "y1", "neo", "p12", "i3", "a3", "sw1", "hyak"]
SECOND = ["q2", "k2", "t2", "p1", "e3", "l3", "l4", "a6", "v2j", "v2z"]
THIRD = ["us10", "k1", "g2", "e2", "n2"]
SAMPLE = FIRST + SECOND + THIRD
# apparition windows for bin_cobs.py (None = all records); C/2014 UN271 is converted by cobs_un271.py
WINDOWS = {k: (None, None) for k in SAMPLE if k != "un271"}
WINDOWS["p12"] = (datetime(2020, 1, 1), None)
WINDOWS["sw1"] = (datetime(2010, 1, 1), datetime(2026, 11, 17))
COBS_DES = {"hb": "C/1995 O1", "hyak": "C/1996 B2", "neo": "C/2020 F3", "p12": "12P", "a3": "C/2023 A3", "i3": "3I",
            "t5": "C/2024 T5", "y1": "C/2025 Y1", "sw1": "29P",
            "q2": "C/2014 Q2", "k2": "C/2017 K2", "t2": "C/2017 T2", "p1": "C/2009 P1", "e3": "C/2022 E3",
            "l3": "C/2019 L3", "l4": "C/2011 L4", "a6": "C/2025 A6", "v2j": "C/2015 V2", "v2z": "C/2020 V2",
            "us10": "C/2013 US10", "k1": "C/2012 K1", "g2": "C/2021 G2", "e2": "C/2022 E2", "n2": "C/2022 N2"}
N_SWITCH = 4.5           # threshold of the switch-on flag (Sec. II C)

def n_local(r, x, chi, h=1e-3):
    """Local activity index of a dominant coma, n = 7/2 - d ln(eta)/d ln(r)."""
    return 3.5 - (np.log(eta(r * (1 + h), x, chi)) - np.log(eta(r * (1 - h), x, chi))) / (np.log(1 + h) - np.log(1 - h))
CASES = [("CO", 1.0), ("CO2", 1.0), ("CO2", 0.25), ("H2O", 1.0), ("H2O", 0.25)]

def anomaly(c, J):
    """Heliocentric distance and true anomaly for an ellipse, parabola or hyperbola (note on the orbit)."""
    q, e, dt = c["q"], c["e"], J - c["Tp"]
    if e < 1:
        a = q / (1 - e); M = KG / a**1.5 * dt; M -= 2 * np.pi * np.round(M / (2 * np.pi))
        E = M + e * np.sin(M) if e < 0.8 else np.pi * np.sign(M) * min(1, abs(M) / np.pi) ** (1 / 3)
        for _ in range(200):
            dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E)); E -= dE
            if abs(dE) < 1e-14: break
        x, y = a * (np.cos(E) - e), a * np.sqrt(1 - e * e) * np.sin(E)
    else:
        a = q / (e - 1); M = KG / a**1.5 * dt
        H = np.arcsinh(M / e)
        for _ in range(200):
            dH = (e * np.sinh(H) - H - M) / (e * np.cosh(H) - 1); H -= dH
            if abs(dH) < 1e-14 * max(1, abs(H)): break
        x, y = a * (e - np.cosh(H)), a * np.sqrt(e * e - 1) * np.sinh(H)
    return np.hypot(x, y), np.arctan2(y, x)

def earth(J):
    T = (J - 2451545.0) / 36525; L0 = (280.46646 + 36000.76983 * T) % 360; M = (357.52911 + 35999.05029 * T) % 360
    Mr = M * D2R; C = (1.914602 - 0.004817 * T) * np.sin(Mr) + 0.019993 * np.sin(2 * Mr) + 0.000289 * np.sin(3 * Mr)
    lam = (L0 + C) * D2R; ecc = 0.016708634 - 0.000042037 * T
    R = 1.000001018 * (1 - ecc**2) / (1 + ecc * np.cos(Mr + C * D2R))
    return -R * np.array([np.cos(lam), np.sin(lam), 0.0])

def ephem(c, J):
    r, nu = anomaly(c, J); i, O, u = c["i"] * D2R, c["Om"] * D2R, nu + c["w"] * D2R
    P = r * np.array([np.cos(O) * np.cos(u) - np.sin(O) * np.sin(u) * np.cos(i),
                      np.sin(O) * np.cos(u) + np.cos(O) * np.sin(u) * np.cos(i), np.sin(u) * np.sin(i)])
    E = earth(J); D = np.linalg.norm(P - E); R = np.linalg.norm(E)
    return r, D, np.degrees(np.arccos(np.clip((r * r + D * D - R * R) / (2 * r * D), -1, 1)))

def load(path=os.path.join(HERE, "cobs", "cobs_lightcurves.txt")):
    data, key = {}, None
    for line in open(path):
        if line.startswith("#"): key = line[1:].strip(); data[key] = []
        elif line.strip() and not line.startswith("%"):
            d, m, s = line.split()[:3]; data[key].append((jd(d + "T00:00"), float(m), float(s)))
    return data

def model(c, pts, x, chi, lk, setup):
    """Eq. (14) at the points (r, D, alpha) for ice x, <cos zeta> = chi and K = 10**lk cm."""
    HN, beta, j, ap, msun = c["HN"], 0.04, setup["j"], setup["ap"], setup["msun"]
    mu = 10**lk * (ap * ARCSEC / (4 * AU_CM) if j == 1 else ap * 1e5 / (4 * AU_CM**2)) * 10**(0.4 * (HN - msun))
    r, D, al = pts[:, 0], pts[:, 1], pts[:, 2]
    et = np.array([eta(rv, x, chi) for rv in r]); ph = al * beta
    return HN + 5 * np.log10(r * D) + ph - 2.5 * np.log10(1 + mu * et * r**-1.5 * D**j * 10**(0.4 * ph))

def fit_K(c, pts, m, s, x, chi, setup):
    f = lambda lk: np.sum(((m - model(c, pts, x, chi, lk, setup)) / s) ** 2)
    grid = np.arange(-1, 10, 0.05); lk = grid[np.argmin([f(v) for v in grid])]
    a, b, g = lk - 0.05, lk + 0.05, (np.sqrt(5) - 1) / 2
    for _ in range(60):
        c1, c2 = b - g * (b - a), a + g * (b - a)
        if f(c1) < f(c2): b = c2
        else: a = c1
    lk = (a + b) / 2; return lk, f(lk)

PAIRS = [("CO", 1.0, "H2O", 1.0), ("CO", 1.0, "H2O", 0.25), ("CO2", 1.0, "H2O", 1.0), ("CO2", 0.25, "H2O", 0.25)]

def model2(c, pts, pair, KV, KW, setup, etas=None):
    """Two ices on separate patches: Afrho = (K_V eta_V + K_W eta_W) r^-3/2 in Eq. (14)."""
    xv, cv, xw, cw = pair
    HN, beta, j, ap, msun = c["HN"], 0.04, setup["j"], setup["ap"], setup["msun"]
    conv = (ap * ARCSEC / (4 * AU_CM) if j == 1 else ap * 1e5 / (4 * AU_CM**2)) * 10**(0.4 * (HN - msun))
    r, D, al = pts[:, 0], pts[:, 1], pts[:, 2]
    ev, ew = etas if etas is not None else (np.array([eta(rv, xv, cv) for rv in r]), np.array([eta(rv, xw, cw) for rv in r]))
    ph = al * beta
    return HN + 5 * np.log10(r * D) + ph - 2.5 * np.log10(1 + conv * (KV * ev + KW * ew) * r**-1.5 * D**j * 10**(0.4 * ph))

def fit2(c, pts, m, s, setup):
    """Best pair and constants K_V, K_W >= 0 by weighted least squares in magnitude (on the given points)."""
    from scipy.optimize import least_squares
    best = None
    for pair in PAIRS:
        lv, _ = fit_K(c, pts, m, s, pair[0], pair[1], setup); lw, _ = fit_K(c, pts, m, s, pair[2], pair[3], setup)
        etas = (np.array([eta(rv, pair[0], pair[1]) for rv in pts[:, 0]]), np.array([eta(rv, pair[2], pair[3]) for rv in pts[:, 0]]))
        for u0, v0 in ((1.0, 1e-3), (1e-3, 1.0), (0.5, 0.5)):
            f = lambda q: (m - model2(c, pts, pair, q[0] * 10**lv, q[1] * 10**lw, setup, etas)) / s
            res = least_squares(f, [u0, v0], bounds=([0, 0], [np.inf, np.inf]), x_scale=[1, 1])
            c2 = float(np.sum(res.fun**2))
            if best is None or c2 < best[3]: best = (pair, res.x[0] * 10**lv, res.x[1] * 10**lw, c2)
    return best

def fit_empirical(pts, m, s):
    X = np.column_stack([np.ones(len(m)), 2.5 * np.log10(pts[:, 0])]); y = m - 5 * np.log10(pts[:, 1]); w = 1 / s
    (m0, n), *_ = np.linalg.lstsq(X * w[:, None], y * w, rcond=None); return m0, n

def fit_m0(pts, m, s, n):
    """Empirical law with a fixed activity index n: only m0 is fitted (one constant, like Eq. (14))."""
    w = 1 / s**2; return float(np.sum(w * (m - 5 * np.log10(pts[:, 1]) - 2.5 * n * np.log10(pts[:, 0]))) / np.sum(w))

def empirical(pts, m0, n):
    return m0 + 5 * np.log10(pts[:, 1]) + 2.5 * n * np.log10(pts[:, 0])

rms = lambda d: float(np.sqrt(np.mean(np.square(d))))

def best_case(c, pts, m, s, setup):
    out = None
    for x, chi in CASES:
        lk, c2 = fit_K(c, pts, m, s, x, chi, setup)
        if out is None or c2 < out[3]: out = (x, chi, lk, c2)
    return out

def run(verbose=True):
    data, rows = load(), []
    pro = [(jd("2021-12-15T00:00"), 18.99, 0.05), (jd("2024-02-15T00:00"), 17.38, 0.03)]   # LOOK, LCO (r' in 6")
    for key in SAMPLE:
        c = COMETS[key]
        setup = dict(j=1, ap=6.0, msun=-26.93) if key == "un271" else dict(j=0, ap=20000.0, msun=-26.762)
        obs = sorted(data[key]) if key != "un271" else pro + sorted(data[key])
        J = np.array([o[0] for o in obs]); m = np.array([o[1] for o in obs]); s = np.array([o[2] for o in obs])
        pts = np.array([ephem(c, v) for v in J])
        if key == "un271":
            tr, te = np.arange(2), np.arange(2, len(J))
            x, chi = "CO2", 1.0; lk, _ = fit_K(c, pts[1:2], m[1:2], s[1:2], x, chi, setup)    # K from the 2024 point
        else:
            h = len(J) // 2; tr, te = np.arange(h), np.arange(h, len(J))
            x, chi, lk, _ = best_case(c, pts[tr], m[tr], s[tr], setup)
        m0, n = fit_empirical(pts[tr], m[tr], s[tr])
        f_mod = rms(m[te] - model(c, pts[te], x, chi, lk, setup)); f_emp = rms(m[te] - empirical(pts[te], m0, n))
        f_fix = {nf: rms(m[te] - empirical(pts[te], fit_m0(pts[tr], m[tr], s[tr], nf), nf)) for nf in (3.5, 4.0)}
        # hindsight: both fitted to all points (case chosen on all points)
        if key == "un271": xa, ca, la = x, chi, fit_K(c, pts, m, s, x, chi, setup)[0]
        else: xa, ca, la, _ = best_case(c, pts, m, s, setup)
        a0, an = fit_empirical(pts, m, s)
        h_mod = rms(m - model(c, pts, xa, ca, la, setup)); h_emp = rms(m - empirical(pts, a0, an))
        # two ices, fitted on the first half (C/2014 UN271: as for one ice, from its two professional points)
        if key == "un271": f2, pair2, kv2, kw2 = f_mod, (x, chi, "H2O", chi), 10**lk, 0.0
        else:
            pair2, kv2, kw2, _ = fit2(c, pts[tr], m[tr], s[tr], setup)
            f2 = rms(m[te] - model2(c, pts[te], pair2, kv2, kw2, setup))
        # whole apparition with one mu before and one after perihelion (same ice as the one-mu fit)
        jj = np.array([o[0] for o in obs]); bpre, bpost = jj < c["Tp"], jj >= c["Tp"]
        if bpre.sum() >= 3 and bpost.sum() >= 3:
            lpre, _ = fit_K(c, pts[bpre], m[bpre], s[bpre], xa, ca, setup); lpost, _ = fit_K(c, pts[bpost], m[bpost], s[bpost], xa, ca, setup)
            mm = np.where(bpre, model(c, pts, xa, ca, lpre, setup), model(c, pts, xa, ca, lpost, setup))
            h_two, ratio = rms(m - mm), 10**(lpost - lpre)
        else: h_two, ratio = np.nan, np.nan
        n_max = max(n_local(rv, x, chi) for rv in pts[tr, 0])
        rows.append((c["name"], len(tr), len(te), x + ("" if chi == 1 else " iso"), f_mod, f_emp, n, h_mod, h_emp,
                     f_fix[3.5], f_fix[4.0], n_max, n_max > N_SWITCH, f2, pair2, kv2, kw2, h_two, ratio,
                     xa + ("" if ca == 1 else " iso")))
    if verbose:
        print(f"{'comet':30s} fit test  case     n_max switch  forecast rms: model  n fitted (n)  n = 3.5  n = 4   hindsight: model  emp.")
        for r in rows:
            print(f"{r[0]:30s} {r[1]:3d} {r[2]:4d}  {r[3]:8s} {r[11]:5.1f}  {'yes' if r[12] else 'no ':3s}             {r[4]:5.2f}   {r[5]:5.2f} ({r[6]:4.1f})  {r[9]:5.2f}  {r[10]:5.2f}            {r[7]:5.2f}  {r[8]:5.2f}")
        from scipy.stats import binomtest, mannwhitneyu
        def score(sel, j):
            w = sum(r[4] < r[j] - 0.05 for r in sel); l = sum(r[4] > r[j] + 0.05 for r in sel)
            p = binomtest(w, w + l, 0.5, alternative="greater").pvalue if w + l else float("nan")
            return w, len(sel) - w - l, l, p, float(np.median([r[j] - r[4] for r in sel]))
        on = [r for r in rows if r[12]]; off = [r for r in rows if not r[12]]
        for grp, sel in (("all 25", rows), (f"ice switches on (n > {N_SWITCH})", on), ("no switch-on", off)):
            print(f"-- {grp}: {len(sel)} comets")
            for lab, j in (("n fitted (m0, n)", 5), ("n = 3.5 fixed (m0)", 9), ("n = 4 fixed (m0)", 10)):
                w, e, l, p, g = score(sel, j)
                print(f"model vs empirical law, {lab:18s}: better {w}, equal within 0.05 mag {e}, worse {l}; "
                      f"sign test p = {p:.3f}; median rms gain {g:+.2f} mag")
            print(f"criterion (iii), rms <= 0.75 mag: model {sum(r[4] <= 0.75 for r in sel)}, n fitted "
                  f"{sum(r[5] <= 0.75 for r in sel)}, n = 3.5 {sum(r[9] <= 0.75 for r in sel)}, n = 4 {sum(r[10] <= 0.75 for r in sel)} of {len(sel)}")
        for lab, j in (("n fitted", 5), ("n = 3.5", 9)):
            u = mannwhitneyu([r[j] - r[4] for r in on], [r[j] - r[4] for r in off], alternative="greater")
            print(f"gain over {lab}: switch-on comets vs the rest, one-sided Mann-Whitney p = {u.pvalue:.3f}")
        print(f"\ntwo ices (K_V, K_W >= 0, pair chosen on the first half)")
        print(f"{'comet':30s} pair                  K_W/K_V   forecast rms: two ices  one ice  n fitted  n = 3.5")
        for r in rows:
            pr = r[14]; lab = f"{pr[0]}{'' if pr[1] == 1 else ' iso'} + {pr[2]}{'' if pr[3] == 1 else ' iso'}"
            print(f"{r[0]:30s} {lab:20s} {r[16] / r[15] if r[15] > 0 else np.inf:9.3g}   {r[13]:5.2f}   {r[4]:5.2f}   {r[5]:5.2f}   {r[9]:5.2f}")
        for grp, sel in (("all 25", rows), ("switch-on", on), ("no switch-on", off)):
            for lab, j in (("one ice", 4), ("n fitted", 5), ("n = 3.5", 9)):
                w = sum(r[13] < r[j] - 0.05 for r in sel); l = sum(r[13] > r[j] + 0.05 for r in sel)
                p = binomtest(w, w + l, 0.5, alternative="greater").pvalue if w + l else float("nan")
                print(f"  {grp:12s} two ices vs {lab:8s}: better {w}, equal {len(sel) - w - l}, worse {l}; sign test p = {p:.3f}; "
                      f"median gain {float(np.median([r[j] - r[13] for r in sel])):+.2f} mag")
        print(f"\nwhole apparitions (hindsight rms): one mu, mu before/after perihelion, empirical law; mu_post/mu_pre")
        for r in rows:
            print(f"{r[0]:30s} {r[19]:8s} {r[7]:5.2f}  {r[17]:5.2f}  {r[8]:5.2f}   {r[18]:7.3g}")
        rt = [r[18] for r in rows if np.isfinite(r[18])]
        print(f"mu_post/mu_pre: median {np.median(rt):.2f}, range {min(rt):.2f}-{max(rt):.2f}, below 1 for {sum(v < 1 for v in rt)} of {len(rt)}")
    return rows

def water_forced(keys=("hb", "a3")):
    """Forecast rms when the expected ice, subsolar water, is imposed instead of the fitted case (note)."""
    data = load()
    for key in keys:
        c = COMETS[key]; setup = dict(j=0, ap=20000.0, msun=-26.762)
        obs = sorted(data[key]); J = np.array([o[0] for o in obs])
        m = np.array([o[1] for o in obs]); s = np.array([o[2] for o in obs])
        pts = np.array([ephem(c, v) for v in J]); h = len(J) // 2; tr, te = np.arange(h), np.arange(h, len(J))
        lk, _ = fit_K(c, pts[tr], m[tr], s[tr], "H2O", 1.0, setup)
        print(f"{c['name']:30s} subsolar water imposed: forecast rms {rms(m[te] - model(c, pts[te], 'H2O', 1.0, lk, setup)):.2f} mag")

if __name__ == "__main__":
    run()
    water_forced()
