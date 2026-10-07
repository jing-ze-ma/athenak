#!/usr/bin/env python3
"""ESTIMATE: phi Per progenitor during Case-B transfer, CONSERVATIVE (beta=1) track from the
Schootemeijer+2018 best fit M1,0=7.2, M2,0=3.8, P0=16 d (arXiv:1803.02379 s4.2):
P = P0 (M1,0 M2,0 / (Md Ma))^3. Gainer radius bracket: R_TE = M^0.6 (ZAMS-like B star, ESTIMATE)
and 2 R_TE (swollen by thermal-timescale accretion, ESTIMATE). r_min from LS75 fit."""
import numpy as np
from scipy.optimize import brentq
G, Msun, Rsun = 6.674e-8, 1.989e33, 6.957e10
M10, M20, P0 = 7.2, 3.8, 16.0
print("Md   Ma    q     P[d]   a[Rs]  dL1  rmin  rcirc  R_TE rmin/R_TE rcirc/R_TE rmin/(2R_TE)")
for Md in (6.5, 5.5, 4.5, 3.5, 2.5, 1.8, 1.2):
    Ma = M10 + M20 - Md; q = Md/Ma
    P = P0*(M10*M20/(Md*Ma))**3
    a = (G*(Md+Ma)*Msun*(P*86400)**2/(4*np.pi**2))**(1/3)/Rsun
    mu = Md/(Md+Ma); xa, xd = -mu, 1-mu
    gx = lambda x: -(1-mu)*(x-xa)/abs(x-xa)**3 - mu*(x-xd)/abs(x-xd)**3 + x
    dL1 = (brentq(gx, xa+1e-4, xd-1e-4)-xa)*a
    rmin = 0.0488*q**-0.464*a; rc = 0.0859*q**-0.426*a; R = Ma**0.6
    print(f"{Md:4.1f} {Ma:5.1f} {q:5.2f} {P:7.1f} {a:6.1f} {dL1:5.1f} {rmin:5.2f} {rc:5.2f} {R:5.2f} "
          f"{rmin/R:8.2f} {rc/R:9.2f} {rmin/(2*R):9.2f}")
