"""Light-curve model for comet C/2014 UN271 (Bernardinelli-Bernstein).
Two-component (nucleus + coma) magnitude model with Afrho and an energy-balance
sublimation share eta_X(r) for CO, CO2 and H2O ice (subsolar, sunlit hemisphere, isothermal).
Reproduces the numbers of Secs. II and IV A; backtest.py, predictions.py, cobs_un271.py and
forecast_test.py reproduce the rest.  Requires numpy and scipy."""
import os, sys
import numpy as np
from scipy.optimize import brentq
from datetime import datetime, timedelta
# un271_orbit.py must sit in the same folder as this file
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from un271_orbit import ephem

lg = np.log10
AU_cm = 1.495978707e13; AU_m = 1.495978707e11
ARCSEC = np.pi / 180 / 3600
S0, SIG, KB, AMU = 1361.0, 5.670374e-8, 1.380649e-23, 1.66054e-27
A_B, EPS = 0.01, 0.9
M_SUN_R = -26.93                     # Sun, SDSS r' (AB), Willmer 2018 (ApJS 236, 47)
M_SUN_V = -26.762                    # Sun, Johnson V (Willmer 2018: -26.76); gives the 1329 km of Eq. (4)
ICES = {"H2O": (3.56e12, 6141.667, 18.015), "CO2": (1.23e12, 3167.8, 44.01),
        "CO": (1.26e9, 764.16, 28.01)}  # A_X [Pa], B_X [K], molar mass

def latent(x):
    A, B, mu = ICES[x]; return B * KB / (mu * AMU)

def Z(T, x):
    A, B, mu = ICES[x]; m = mu * AMU
    return A * np.exp(-B / T) * np.sqrt(m / (2 * np.pi * KB * T))

def balance(r, x, chi=1.0):
    """Surface temperature, sublimation flux [kg m^-2 s^-1], energy share eta.
    chi is <cos zeta> of the paper: 1 subsolar, 1/4 isothermal."""
    fin = (1 - A_B) * S0 * chi / r**2
    T = brentq(lambda T: fin - EPS * SIG * T**4 - Z(T, x) * latent(x), 5, 400)
    return T, Z(T, x), Z(T, x) * latent(x) / fin

def eta(r, x, chi=1.0):
    return balance(r, x, chi)[2]

def eta_hemi(r, x, n=48):
    """Sunlit hemisphere of a slowly rotating nucleus: every surface element balances at its own
    cos(zeta) = chi, so the share of the absorbed power spent on sublimation is
    eta_hemi = 2 * int_0^1 chi * eta(r, x, chi) dchi (Gauss-Legendre in chi); Eq. (11) keeps its form."""
    t, w = np.polynomial.legendre.leggauss(n)
    chi, w = 0.5 * (t + 1), 0.5 * w
    return 2 * float(sum(wi * ci * eta(r, x, ci) for ci, wi in zip(chi, w)))

def eta_geom(r, x, g):
    """eta for the three geometries of Sec. II B: g = 1 (subsolar), 0.25 (isothermal) or "hemi"."""
    return eta_hemi(r, x) if g == "hemi" else eta(r, x, g)

def balance_mantle(r, x, chi=1.0, h=np.inf):
    """Ice under a dust mantle of thermal conductance h = kappa/d [W m^-2 K^-1] (note on the
    mantle): surface (1-A_B) S0 chi/r^2 = eps sigma Ts^4 + h (Ts - Ti), ice front
    h (Ts - Ti) = Z(Ti) L.  Returns Ts, Ti, flux Z(Ti), eta; h = inf is Eq. (8)."""
    if np.isinf(h):
        T, z, e = balance(r, x, chi); return T, T, z, e
    fin, L = (1 - A_B) * S0 * chi / r**2, latent(x)
    def res(Ti):
        zl = Z(Ti, x) * L; return fin - EPS * SIG * (Ti + zl / h)**4 - zl
    Ti = brentq(res, 5, 400); zl = Z(Ti, x) * L
    return Ti + zl / h, Ti, Z(Ti, x), zl / fin

# --- nucleus (Hui et al. 2022) --------------------------------------------
H_V, VR, rR, BETA_N = 8.62, 0.40, 0.18, 0.04
R_N_KM = 68.5

def m_nucleus(r, D, alpha, HV=H_V, vr=VR, rr=rR, beta=BETA_N):
    return HV - vr + rr + 5 * lg(r * D) + beta * alpha

def afrho_from_mag(m_coma, r, D, theta_arcsec, msun=M_SUN_R):
    rho = theta_arcsec * ARCSEC * D * AU_cm
    return (2 * D * AU_cm * r)**2 / rho * 10**(0.4 * (msun - m_coma))

def mag_from_afrho(afrho, r, D, theta_arcsec, msun=M_SUN_R):
    rho = theta_arcsec * ARCSEC * D * AU_cm
    return msun - 2.5 * lg(afrho * rho / (4 * r**2 * (D * AU_cm)**2))

def add_mags(m1, m2):
    return -2.5 * lg(10**(-0.4 * m1) + 10**(-0.4 * m2))

# --- data ------------------------------------------------------------------
EP1 = dict(label="LOOK mean Oct 2021-Feb 2022", date=datetime(2021, 12, 15),
           H=5.98, sH=0.05, phase_corr=True)        # Kelley et al. 2022
EP2 = dict(label="LCO mean Feb 2024", date=datetime(2024, 2, 15),
           m=17.38, sm=0.03)                        # Roth et al. 2025
BETA_C, S_BETA_C = 0.035, 0.010                     # coma phase coefficient (approx.)

def epoch_values(ep, rng=None):
    g = ephem(ep["date"]); r, D, a = g["r"], g["D"], g["phase"]
    if rng is None:
        HV, vr, rr, bn, bc = H_V, VR, rR, BETA_N, BETA_C; dH = dm = 0
    else:
        HV = rng.normal(H_V, 0.11); vr = rng.normal(VR, 0.02); rr = rng.normal(rR, 0.05)
        bn = rng.normal(BETA_N, 0.01); bc = rng.normal(BETA_C, S_BETA_C)
        dH = rng.normal(0, ep.get("sH", 0)); dm = rng.normal(0, ep.get("sm", 0))
    if "H" in ep:
        m = ep["H"] + dH + 5 * lg(r * D) + bc * a   # undo the dust phase correction
    else:
        m = ep["m"] + dm
    mN = m_nucleus(r, D, a, HV, vr, rr, bn)
    Fc = 10**(-0.4 * m) - 10**(-0.4 * mN)
    mc = -2.5 * lg(Fc)
    return dict(r=r, D=D, alpha=a, m=m, mN=mN, share=10**(-0.4 * (mN - m)), mc=mc,
                rho_km=6 * ARCSEC * D * AU_m / 1e3, afrho=afrho_from_mag(mc, r, D, 6.0))

def slope(e1, e2):
    return -np.log(e2["afrho"] / e1["afrho"]) / np.log(e2["r"] / e1["r"])

if __name__ == "__main__":
    print("Latent heats [J/kg]:", {x: f"{latent(x):.3e}" for x in ICES})
    print("\nEnergy balance (subsolar chi=1 | isothermal chi=1/4):")
    print(" r     " + "  ".join(f"{x:>4s}: T     eta      " for x in ICES))
    for r in [1, 3, 5, 7.5, 10.95, 13, 16.7, 19.5, 25, 34.1]:
        s = f"{r:5.2f} "
        for x in ICES:
            T, z, e = balance(r, x); T4, z4, e4 = balance(r, x, 0.25)
            s += f" {x:>4s}: {T:5.1f} {e:6.3f}|{e4:6.3f}"
        print(s)
    # note on the root of Eq. (8): dropping one term of the balance (subsolar CO2 at the two epochs)
    for r in (16.70, 19.53):
        T, z, e = balance(r, "CO2")
        T_rad = ((1 - A_B) * S0 / (EPS * SIG * r**2))**0.25          # radiation only
        z_lim = (1 - A_B) * S0 / (r**2 * latent("CO2"))               # sublimation only (eta = 1)
        print(f"CO2 at {r} au: eta={e:.3f}; radiation only: Z x {Z(T_rad, 'CO2')/z:.1f}, "
              f"sublimation only: Z x {z_lim/z:.1f}")
    e1, e2 = epoch_values(EP1), epoch_values(EP2)
    for lab, e in [("EP1", e1), ("EP2", e2)]:
        print(f"\n{lab}: r={e['r']:.2f} D={e['D']:.2f} alpha={e['alpha']:.1f} m={e['m']:.2f} "
              f"m_N={e['mN']:.2f} share={100*e['share']:.1f}% m_c={e['mc']:.2f} "
              f"rho={e['rho_km']:.0f} km Afrho={e['afrho']:.0f} cm")
    k = slope(e1, e2); print(f"\nk (nominal) = {k:.2f}")
    rng = np.random.default_rng(2031)
    ks, A1, A2 = [], [], []
    for _ in range(20000):
        a, b = epoch_values(EP1, rng), epoch_values(EP2, rng)
        ks.append(slope(a, b)); A1.append(a["afrho"]); A2.append(b["afrho"])
    q = lambda v: np.percentile(v, [16, 50, 84])
    print("k MC 16/50/84:", q(ks)); print("Afrho1:", q(A1)); print("Afrho2:", q(A2))
    # model slopes between the two epochs
    r1, r2 = e1["r"], e2["r"]
    for x in ["CO", "CO2"]:
        for chi, lab in [(1.0, "subsolar"), (0.25, "isothermal")]:
            et1, et2 = eta(r1, x, chi), eta(r2, x, chi)
            kx = 1.5 + np.log(et2 / et1) / np.log(r1 / r2) if et1 > 0 else np.nan
            print(f"{x} {lab}: eta1={et1:.3e} eta2={et2:.3e} -> k(v~r^-0.5)={kx:.2f}, k(v const)={kx+0.5:.2f}")
    # CO active fraction from ALMA (Roth et al. 2025)
    R = R_N_KM * 1e3
    mmax = (1 - A_B) * S0 / 16.6**2 * np.pi * R**2 / latent("CO") * eta(16.6, "CO")
    for Q in (3.44e27, 5.33e27):
        print(f"Q_CO={Q:.2e} -> {Q*28.01*AMU:.0f} kg/s, fraction of max {Q*28.01*AMU/mmax*100:.3f}%")
    mq = (1 - A_B) * S0 / 10.95**2 * np.pi * R**2 / latent("CO") * eta(10.95, "CO")
    print(f"max CO at q: {mq:.3e} kg/s; at f=0.07%: {7e-4*mq:.0f} kg/s = {7e-4*mq/(28.01*AMU):.2e} /s")
    print(f"max CO at 16.6: {mmax:.3e} kg/s = {mmax/(28.01*AMU):.2e} /s")
    # Scenario light curves (6" aperture, r'), anchored on EP2
    A0, r0 = e2["afrho"], e2["r"]; eta0 = eta(r0, "CO2")
    def afrho_s(r, s):
        if s == "co": return A0 * (r / r0)**-1.5
        if s == "co2": return A0 * (eta(r, "CO2") / eta0) * (r / r0)**-1.5
        if s == "pl": return A0 * (r / r0)**-k
    rows = []; t = datetime(2021, 7, 1)
    while t <= datetime(2032, 7, 1):
        g = ephem(t); r, D, a = g["r"], g["D"], g["phase"]
        mN = m_nucleus(r, D, a)
        row = dict(t=t.year + (t.month - 1) / 12, date=t.strftime("%Y-%m"), r=r, D=D, nuc=mN)
        for s in ("co", "co2", "pl"):
            row[s] = add_mags(mN, mag_from_afrho(afrho_s(r, s), r, D, 6.0))
        rows.append(row)
        t = (t.replace(day=1) + timedelta(days=95)).replace(day=1)
    print("\n date      r     D    nuc    co    co2    pl")
    for w in rows:
        print(f"{w['date']} {w['r']:5.2f} {w['D']:5.2f} {w['nuc']:6.2f} {w['co']:6.2f} {w['co2']:6.2f} {w['pl']:6.2f}")
    np.save("scen.npy", rows, allow_pickle=True)
    for s in ("co", "co2", "pl"):
        print(s, "Afrho(q) =", round(afrho_s(10.95, s)))
