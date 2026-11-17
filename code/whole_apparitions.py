"""Whole-apparition fits of the COBS light curves (Sec. III C, Fig. 4 and the supplementary figure).

For each comet: Eq. (14) with the ice chosen on all points and one coma constant mu; the same ice with
one mu before and one after perihelion (their ratio measures the change of the active area); and the
empirical law, Eq. (1), with m0 and n fitted.  C/2014 UN271 uses its r' magnitudes in 6" (j = 1).
Usage: python whole_apparitions.py          (table and supplementary figure whole_apparitions.pdf)
       python whole_apparitions.py fig4     (pgfplots coordinates of Fig. 4)
Needs numpy, scipy and, for the figure, matplotlib."""
import os, sys
import numpy as np
from datetime import datetime, timedelta
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forecast_test as F

FIG4 = ["hb", "sw1", "hyak", "q2", "a3", "p12"]       # the six light curves with the most COBS records

def year(J):
    t = datetime(2000, 1, 1, 12) + timedelta(days=float(J) - 2451545.0)
    a, b = datetime(t.year, 1, 1), datetime(t.year + 1, 1, 1)
    return t.year + (t - a).total_seconds() / (b - a).total_seconds()

def fits(key, data):
    c = F.COMETS[key]
    setup = dict(j=1, ap=6.0, msun=-26.93) if key == "un271" else dict(j=0, ap=20000.0, msun=-26.762)
    pro = [(F.jd("2021-12-15T00:00"), 18.99, 0.05), (F.jd("2024-02-15T00:00"), 17.38, 0.03)]
    obs = sorted(data[key]) if key != "un271" else pro + sorted(data[key])
    J = np.array([o[0] for o in obs]); m = np.array([o[1] for o in obs]); s = np.array([o[2] for o in obs])
    pts = np.array([F.ephem(c, v) for v in J])
    if key == "un271": x, chi = "CO2", 1.0; lk = F.fit_K(c, pts, m, s, x, chi, setup)[0]
    else: x, chi, lk, _ = F.best_case(c, pts, m, s, setup)
    m0, n = F.fit_empirical(pts, m, s)
    pre, post = J < c["Tp"], J >= c["Tp"]
    two = pre.sum() >= 3 and post.sum() >= 3
    if two:
        lpre = F.fit_K(c, pts[pre], m[pre], s[pre], x, chi, setup)[0]
        lpost = F.fit_K(c, pts[post], m[post], s[post], x, chi, setup)[0]
    # dense grid over the observed span
    Jg = np.linspace(J.min(), J.max(), 160); pg = np.array([F.ephem(c, v) for v in Jg])
    out = dict(name=c["name"], case=x + ("" if chi == 1 else " iso"), J=J, m=m, s=s, Jg=Jg, Tp=c["Tp"],
               one=F.model(c, pg, x, chi, lk, setup), emp=F.empirical(pg, m0, n), n=n,
               rms_one=F.rms(m - F.model(c, pts, x, chi, lk, setup)), rms_emp=F.rms(m - F.empirical(pts, m0, n)))
    if two:
        out["two"] = np.where(Jg < c["Tp"], F.model(c, pg, x, chi, lpre, setup), F.model(c, pg, x, chi, lpost, setup))
        out["rms_two"] = F.rms(m - np.where(pre, F.model(c, pts, x, chi, lpre, setup), F.model(c, pts, x, chi, lpost, setup)))
        out["ratio"] = 10**(lpost - lpre)
    return out

def table(res):
    print(f"{'comet':30s} {'ice':8s} one mu  two mu  Eq. (1)   n    mu_post/mu_pre")
    for r in res:
        print(f"{r['name']:30s} {r['case']:8s} {r['rms_one']:5.2f}  {r.get('rms_two', np.nan):5.2f}   {r['rms_emp']:5.2f}  {r['n']:5.1f}   {r.get('ratio', np.nan):6.2f}")
    two = [r for r in res if "rms_two" in r]
    rat = np.array([r["ratio"] for r in two])
    b1 = sum(r["rms_two"] < r["rms_one"] - 0.05 for r in two)
    be = sum(r["rms_two"] < r["rms_emp"] - 0.05 for r in two); we = sum(r["rms_two"] > r["rms_emp"] + 0.05 for r in two)
    o1 = sum(r["rms_one"] < r["rms_emp"] - 0.05 for r in res); w1 = sum(r["rms_one"] > r["rms_emp"] + 0.05 for r in res)
    print(f"one mu vs Eq. (1): better {o1}, equal {len(res) - o1 - w1}, worse {w1} of {len(res)}")
    print(f"two mu: {len(two)} comets with points on both sides of perihelion; better than one mu (by > 0.05 mag) for {b1};"
          f" vs Eq. (1): better {be}, equal {len(two) - be - we}, worse {we}")
    print(f"mu_post/mu_pre: median {np.median(rat):.2f}, range {rat.min():.2f}-{rat.max():.2f}; "
          f"below 0.5 for {sum(rat < 0.5)}, above 2 for {sum(rat > 2)}")

def figure(res, path):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.size": 7})
    nc = 5; nr = int(np.ceil(len(res) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(7.2, 1.7 * nr), constrained_layout=True)
    for ax, r in zip(axs.flat, res):
        yg = [year(v) for v in r["Jg"]]
        ax.errorbar([year(v) for v in r["J"]], r["m"], yerr=r["s"], fmt="o", ms=1.8, lw=0.5, color="#c05a28")
        ax.plot(yg, r["one"], color="#1f5fbf", lw=1.0)
        if "two" in r:
            pre = r["Jg"] < r["Tp"]
            ax.plot(np.array(yg)[pre], r["two"][pre], color="#1f5fbf", lw=0.9, ls="--")
            ax.plot(np.array(yg)[~pre], r["two"][~pre], color="#1f5fbf", lw=0.9, ls="--")
        ax.plot(yg, r["emp"], color="k", lw=0.8, ls="-.")
        if r["Jg"][0] < r["Tp"] < r["Jg"][-1]: ax.axvline(year(r["Tp"]), color="0.6", lw=0.5, ls=":")
        ax.invert_yaxis(); ax.set_title(f"{r['name'].split(' (')[0]}: {r['case']}", fontsize=6.5)
        ax.tick_params(labelsize=5.5, length=2)
        from matplotlib.ticker import MaxNLocator, FormatStrFormatter
        ax.xaxis.set_major_locator(MaxNLocator(3)); ax.ticklabel_format(axis="x", useOffset=False)
        span = yg[-1] - yg[0]
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f" if span > 3 else "%.1f"))
    for ax in axs.flat[len(res):]: ax.axis("off")
    fig.suptitle("Whole-apparition fits (magnitude against year). Solid: Eq. (14) with one $\\mu$; dashed: one $\\mu$ before "
                 "and one after perihelion;\ndash-dot: Eq. (1) with $m_0$ and $n$; dots: COBS bins; dotted: perihelion. "
                 "The ice of Eq. (14) is chosen on all points.", fontsize=7)
    fig.savefig(path)

def fig4(res):
    for r in res:
        print(f"%{r['name']}  case {r['case']}  rms one {r['rms_one']:.2f} two {r.get('rms_two', np.nan):.2f} emp {r['rms_emp']:.2f}")
        print("% data"); print(" ".join(f"({year(a):.3f},{b:.2f}) +- (0,{e:.2f})" for a, b, e in zip(r["J"], r["m"], r["s"])))
        for lab in ("one", "emp"):
            print(f"% {lab}"); print(" ".join(f"({year(a):.3f},{b:.2f})" for a, b in zip(r["Jg"], r[lab])))
        if "two" in r:
            pre = r["Jg"] < r["Tp"]
            for part, sel in (("two pre", pre), ("two post", ~pre)):
                print(f"% {part}"); print(" ".join(f"({year(a):.3f},{b:.2f})" for a, b in zip(r["Jg"][sel], r["two"][sel])))
        print(f"% perihelion {year(r['Tp']):.3f}")

if __name__ == "__main__":
    data = F.load()
    if sys.argv[1:2] == ["fig4"]:
        fig4([fits(k, data) for k in FIG4])
    else:
        res = [fits(k, data) for k in F.SAMPLE]
        table(res)
        figure(res, os.path.join(F.HERE, "whole_apparitions.pdf"))
