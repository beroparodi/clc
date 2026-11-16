import numpy as np
from datetime import datetime, timedelta
D2R=np.pi/180
k=0.01720209895
# Elements (JPL SBDB via theskylive, epoch 2018-01-13)
q=10.94525581; e=0.99869477; inc=95.55335892*D2R; Om=190.0846361196*D2R; w=326.33161880945*D2R
Tp=2462892.9068
a=q/(1-e); n=k/a**1.5
def jd(dt):
    return (dt - datetime(2000,1,1,12)).total_seconds()/86400+2451545.0
def comet_helio(J):
    M=n*(J-Tp)
    # solve Kepler for near-parabolic ellipse: E - e sinE = M
    E=np.cbrt(6*M) if abs(M)>0 else 0.0
    for _ in range(100):
        f=E-e*np.sin(E)-M; fp=1-e*np.cos(E)
        dE=-f/fp; E+=dE
        if abs(dE)<1e-15: break
    x=a*(np.cos(E)-e); y=a*np.sqrt(1-e*e)*np.sin(E)
    r=np.hypot(x,y); nu=np.arctan2(y,x)
    u=nu+w
    X=r*(np.cos(Om)*np.cos(u)-np.sin(Om)*np.sin(u)*np.cos(inc))
    Y=r*(np.sin(Om)*np.cos(u)+np.cos(Om)*np.sin(u)*np.cos(inc))
    Z=r*(np.sin(u)*np.sin(inc))
    return np.array([X,Y,Z]), r
def earth_helio(J):
    T=(J-2451545.0)/36525
    L0=(280.46646+36000.76983*T)%360
    M=(357.52911+35999.05029*T)%360
    Mr=M*D2R
    C=(1.914602-0.004817*T)*np.sin(Mr)+0.019993*np.sin(2*Mr)+0.000289*np.sin(3*Mr)
    lam=(L0+C)*D2R
    ecc=0.016708634-0.000042037*T
    R=1.000001018*(1-ecc**2)/(1+ecc*np.cos(Mr+C*D2R))
    # sun geocentric ecliptic -> earth heliocentric = -sun
    return -R*np.array([np.cos(lam),np.sin(lam),0.0])
eps=23.4392911*D2R
def ephem(dt):
    J=jd(dt)
    C,r=comet_helio(J); Ea=earth_helio(J)
    G=C-Ea; Dl=np.linalg.norm(G)
    # light-time ignored (~1.5 h per 10 AU; negligible here)
    xe=G[0]; ye=G[1]*np.cos(eps)-G[2]*np.sin(eps); ze=G[1]*np.sin(eps)+G[2]*np.cos(eps)
    ra=np.degrees(np.arctan2(ye,xe))%360; dec=np.degrees(np.arcsin(ze/Dl))
    R=np.linalg.norm(Ea)
    phase=np.degrees(np.arccos((r*r+Dl*Dl-R*R)/(2*r*Dl)))
    elong=np.degrees(np.arccos((R*R+Dl*Dl-r*r)/(2*R*Dl)))
    return dict(r=r,D=Dl,ra=ra,dec=dec,phase=phase,elong=elong)
if __name__=="__main__":
    checks=[("ALMA 2024-03-08",datetime(2024,3,8)),("today 2026-11-10",datetime(2026,11,10)),("2027-01-15",datetime(2027,1,15)),
            ("HST 2022-01-08",datetime(2022,1,8)),("2021-09-09 outburst",datetime(2021,9,9)),("Discovery 2014-10-20",datetime(2014,10,20)),
            ("Precovery 2010-11-15",datetime(2010,11,15)),("TESS 2018-09-15",datetime(2018,9,15)),("perihelion 2031-01-23",datetime(2031,1,23)),("2031-04-05",datetime(2031,4,5))]
    for lab,dt in checks:
        x=ephem(dt); print(f"{lab:24s} r={x['r']:6.2f} D={x['D']:6.2f} RA={x['ra']/15:5.2f}h Dec={x['dec']:6.1f} phase={x['phase']:4.1f} elong={x['elong']:5.1f}")
    # min Delta and min dec over 2026-2032
    best=None; mind=None
    t=datetime(2026,1,1)
    while t<datetime(2033,1,1):
        x=ephem(t)
        if best is None or x['D']<best[1]: best=(t,x['D'],x['r'])
        if mind is None or x['dec']<mind[1]: mind=(t,x['dec'])
        t+=timedelta(days=1)
    print("closest approach:",best); print("most southern:",mind)
    # max phase angle
    print("max phase 2026-2032:", max(ephem(datetime(2026,1,1)+timedelta(days=d))['phase'] for d in range(0,2557,5)))
