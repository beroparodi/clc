"""The six comets of Parodi (2013) in the energy-balance framework of the paper.
(1) Energy budget of Eq. (4): implied active fraction f.
(2) Independent water rates from the Jorda et al. (2008) relation.
(3) Effective activity index n_eff(r) of the 2013 fits vs. energy-balance prediction.
Standalone: needs only numpy and scipy (Appendix C of the paper)."""
import numpy as np
from scipy.optimize import brentq

# --- surface energy balance (same as un271_model.py, Eqs. 6-8 of the paper) ---
S0, SIG, KB, AMU = 1361.0, 5.670374e-8, 1.380649e-23, 1.66054e-27
A_B, EPS = 0.01, 0.9
ICES = {"H2O": (3.56e12, 6141.667, 18.015), "CO2": (1.23e12, 3167.8, 44.01),
        "CO": (1.26e9, 764.16, 28.01)}  # A_X [Pa], B_X [K], molar mass

def latent(x):
    A, B, mu = ICES[x]; return B * KB / (mu * AMU)

def Z(T, x):
    A, B, mu = ICES[x]; m = mu * AMU
    return A * np.exp(-B / T) * np.sqrt(m / (2 * np.pi * KB * T))

def eta(r, x, chi=1.0):
    """Share of absorbed power spent on sublimating ice x at distance r [au]."""
    fin = (1 - A_B) * S0 * chi / r**2
    T = brentq(lambda T: fin - EPS * SIG * T**4 - Z(T, x) * latent(x), 5, 400)
    return Z(T, x) * latent(x) / fin

lg = np.log10
L = latent("H2O"); m_h2o = 18.015 * AMU
# comet, phase, m00, mu, nu, q_tab, q_cat, dq, mq, RN
C = [("Hale-Bopp","",0.85,3.70,4.19,0.914,0.914,1.351,-0.7,30),
     ("Pan-STARRS","",6.23,0.76,2.43,0.300,0.300,1.109,0.9,1.0),
     ("Hyakutake","",11.91,389,1.78,0.277,0.2302,1.082,2.5,1.3),
     ("Halley","pre",13.72,2639,3.36,0.575,0.575,1.550,2.9,4.2),
     ("Halley","post",12.52,2639,1.31,0.575,0.575,1.550,2.9,4.2),
     ("Swift-Tuttle","pre",16.19,97305,5.31,1.157,0.9595,1.204,4.8,12),
     ("Swift-Tuttle","post",12.72,9258,5.00,1.157,0.9595,1.204,4.8,12),
     ("Tempel 2","pre",14.58,5.4e9,49.1,1.450,1.450,0.770,9.2,4.7),
     ("Tempel 2","post",13.55,3919,9.72,1.450,1.450,0.770,9.2,4.7)]
def mmax(R_km, r, chi=1.0):
    """Energy-limited water production of the projected nucleus (f=1), kg/s."""
    return np.pi*(R_km*1e3)**2*(1-A_B)*S0*eta(r,"H2O",chi)/(r**2*L)
def eq4(R_km, r, nu): return 600*R_km**2/r**(2+nu)
def jorda(mq, dq):
    mh = mq - 5*lg(dq); Q = 10**(30.675 - 0.2453*mh); return mh, Q, Q*m_h2o
def n_eff(r, mu, nu):
    x = mu/r**nu; return 2 + nu*x/(1+x)
def dln_eta(r, chi=1.0, d=0.01):
    return (np.log(eta(r*(1+d),"H2O",chi))-np.log(eta(r*(1-d),"H2O",chi)))/(np.log(1+d)-np.log(1-d))
def n_pred(r, chi=1.0, slope=1/(2.5*0.2453)):
    """Index n predicted from Q ~ eta/r^2 and flux ~ Q^slope (slope=1: Parodi 2013, Eq. 5)."""
    return slope*(2 - dln_eta(r, chi))

print("Jorda flux-production exponent:", 1/(2.5*0.2453))
print("\nTable C1: energy budget and water rates at perihelion")
for c,ph,m00,mu,nu,qt,qc,dq,mq,R in C:
    if ph=="post": continue
    mx = mmax(R,qc); e4 = eq4(R,qc,nu); mh,Q,mj = jorda(mq,dq)
    print(f"{c:13s} q={qc:.3f} R={R:5.1f} eta={eta(qc,'H2O'):.3f} mmax={mx/1e3:7.1f} t/s  Eq4={e4/1e3:7.1f} f4={e4/mx:6.3f}"
          f"  m_h={mh:5.2f} Q_J={Q:.2e} mJ={mj/1e3:6.1f} t/s fJ={mj/mx:6.3f}  Eq4/J={e4/mj:5.2f}")
# post rows for Halley, Swift-Tuttle, Tempel 2 (Eq. 4 with post nu)
for c,ph,m00,mu,nu,qt,qc,dq,mq,R in C:
    if ph!="post": continue
    mx=mmax(R,qc); e4=eq4(R,qc,nu); print(f"{c} post: Eq4={e4/1e3:.2f} t/s f4={e4/mx:.3f}")
print("\nTable C2: effective activity index")
for c,ph,m00,mu,nu,qt,qc,dq,mq,R in C:
    print(f"{c:13s}{ph:5s} n_eff(q)={n_eff(qc,mu,nu):5.2f} n_eff(1)={n_eff(1,mu,nu):5.2f} n_eff(2)={n_eff(2,mu,nu):5.2f} n_eff(3)={n_eff(3,mu,nu):5.2f}"
          f" | pred(q) Jorda={n_pred(qc):.2f} lin={n_pred(qc,slope=1):.2f} iso={n_pred(qc,0.25):.2f}")
print("\nPredicted n(r):")
for r in [0.3,0.5,0.75,1,1.5,2,2.5,3,3.5,4]:
    print(r, "subsolar Jorda", round(n_pred(r),2), " lin", round(n_pred(r,slope=1),2), " iso Jorda", round(n_pred(r,0.25),2))
# nucleus magnitudes for p=0.04 vs m00
print("\nNucleus H_V (p=0.04) vs m00, m0:")
for c,ph,m00,mu,nu,qt,qc,dq,mq,R in C:
    H=5*lg(1329/(2*R*np.sqrt(0.04))); print(f"{c:13s}{ph:5s} H_N={H:6.2f} m00={m00:6.2f} m0={m00-2.5*lg(1+mu):6.2f} coma/nucleus at 1 AU per m00: {10**(0.4*(H-m00)):.3g}")

# ---- LaTeX rows ---------------------------------------------------------------
from scipy.optimize import brentq
printed = {"Hale-Bopp":"940","Pan-STARRS":"120","Hyakutake":"130","Halley":"200","Swift-Tuttle":"30","Tempel 2":"$<0.2$"}
print("\n% Table C1 rows")
for c,ph,m00,mu,nu,qt,qc,dq,mq,R in C:
    if c=="Tempel 2" and ph=="pre": continue
    if ph=="pre" and c in("Halley","Swift-Tuttle"): pass
    if ph=="post" and c!="Tempel 2": continue
    mx=mmax(R,qc); e4=eq4(R,qc,nu); mh,Q,mj=jorda(mq,dq)
    qs=f"{qc:.3f}" + ("$^a$" if qc!=qt else "")
    name=c+(" (post)" if c=="Tempel 2" else "")
    print(f"{name} & {qs} & {R:g} & {eta(qc,'H2O'):.2f} & {mx/1e3:.1f} & {printed[c]} & {e4/1e3:.{0 if e4>1e4 else 1}f} & {mj/1e3:.1f} & {e4/mx:.2f} & {mj/mx:.2f} \\\\")
print("\n% Table C2 rows")
def r_floor(m00,mu,nu,H):
    if m00<=H: return None
    target=10**(0.4*(m00-H))-1
    return (mu/target)**(1/nu)
for c,ph,m00,mu,nu,qt,qc,dq,mq,R in C:
    H=5*lg(1329/(2*R*np.sqrt(0.04)))
    rf=r_floor(m00,mu,nu,H)
    lab=c+(f" ({ph})" if ph else "")
    ne=n_eff(qc,mu,nu); ne2=n_eff(2,mu,nu)
    print(f"{lab} & {ne:.1f} & {n_pred(qc):.1f}--{n_pred(qc,0.25):.1f} & {ne2:.1f} & {H:.1f} & {m00:.2f} & {'--' if rf is None else f'{rf:.1f}'} \\\\")
print("n_pred(2) subsolar/iso:", round(n_pred(2),2), round(n_pred(2,0.25),2))
# figure data
grid=np.geomspace(0.2,4,40)
print("%fig sub\n"+" ".join(f"({r:.4f},{n_pred(r):.3f})" for r in grid))
print("%fig iso\n"+" ".join(f"({r:.4f},{min(n_pred(r,0.25),40):.3f})" for r in grid))
print("%fig lin\n"+" ".join(f"({r:.4f},{n_pred(r,slope=1):.3f})" for r in grid))
