#!/usr/bin/env python3
"""ESTIMATES: Lubow & Shu 1975 r_min, r_circ (q = Mdon/Macc) and accretor-L1 distance (Eggleton-free,
exact root). Input rows: name P_d Macc Mdon Racc [a_Rsun or 0 -> Kepler]."""
import sys, numpy as np
from scipy.optimize import brentq
G, Msun, Rsun = 6.674e-8, 1.989e33, 6.957e10
rows = [l.split() for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
print(f"{'name':10s} {'P':>6s} {'Ma':>5s} {'Md':>5s} {'q':>5s} {'a':>6s} {'Ra':>5s} {'Ra/a':>6s} "
      f"{'rmin':>6s} {'rcirc':>6s} {'dL1':>6s} {'rmin/R':>6s} {'rcirc/R':>7s} {'aOm':>5s}")
for r in rows:
    n, P, Ma, Md, Ra, a = r[0], *map(float, r[1:6])
    q = Md / Ma
    if a <= 0:
        a = (G * (Ma + Md) * Msun * (P * 86400)**2 / (4 * np.pi**2))**(1/3) / Rsun
    mu = Md / (Ma + Md); xa, xd = -mu, 1 - mu
    gx = lambda x: -(1-mu)*(x-xa)/abs(x-xa)**3 - mu*(x-xd)/abs(x-xd)**3 + x
    dL1 = brentq(gx, xa+1e-4, xd-1e-4) - xa
    rmin = 0.0488 * q**-0.464 * a; rc = 0.0859 * q**-0.426 * a
    aOm = 2*np.pi*a*Rsun/(P*86400)/1e5
    print(f"{n:10s} {P:6.2f} {Ma:5.2f} {Md:5.2f} {q:5.3f} {a:6.1f} {Ra:5.2f} {Ra/a:6.3f} "
          f"{rmin:6.2f} {rc:6.2f} {dL1*a:6.1f} {rmin/Ra:6.2f} {rc/Ra:7.2f} {aOm:5.0f}")
