#!/usr/bin/env python3
"""ESTIMATE: Plaskett's star progenitor on the CONSERVATIVE track of Wade+2026 (arXiv:2609.25526)
Fig. 3: M1,i=18.2, M2,i=16.8, P_i=3.69 d; their Eq. 2: P/Pi = (q/qi)^-3 ((1+q)/(1+qi))^6, q=M2/M1.
Gainer radius: 9 Rsun at M=16 (quoted, start of accretion) scaled as R = 9 (M/16)^0.6 (ESTIMATE,
TE main-sequence slope; accretion swelling would make it larger). Last row: today (Table 1).
Eddington factor (electron scattering, X=0.7): Gamma_e = 10^-4.813 (1+X) L/M (ESTIMATE)."""
import numpy as np
from scipy.optimize import brentq
G, Msun, Rsun = 6.674e-8, 1.989e33, 6.957e10
def geo(P, Ma, Md, a=None):
    if a is None:
        a = (G*(Ma+Md)*Msun*(P*86400)**2/(4*np.pi**2))**(1/3)/Rsun
    mu = Md/(Ma+Md); xa, xd = -mu, 1-mu
    gx = lambda x: -(1-mu)*(x-xa)/abs(x-xa)**3 - mu*(x-xd)/abs(x-xd)**3 + x
    dL1 = (brentq(gx, xa+1e-4, xd-1e-4)-xa)*a
    q = Md/Ma
    return a, dL1, 0.0488*q**-0.464*a, 0.0859*q**-0.426*a
M1i, M2i, Pi = 18.2, 16.8, 3.69
qi = M2i/M1i
print("label        Md    Ma   qdon   P[d]   a     dL1  rmin  rcirc  R_acc rmin/R rcirc/R")
for Md in (17.0, 15.0, 12.0, 10.0, 8.0, 6.62):
    Ma = M1i+M2i-Md; q = Ma/Md
    P = Pi*(q/qi)**-3*((1+q)/(1+qi))**6
    a, d, rm, rc = geo(P, Ma, Md); R = 9*(Ma/16)**0.6
    print(f"track       {Md:5.2f} {Ma:5.1f} {Md/Ma:5.3f} {P:6.2f} {a:5.1f} {d:5.1f} {rm:5.2f} {rc:5.2f} {R:5.1f} {rm/R:6.2f} {rc/R:6.2f}")
for lab, Md, Ma, a, R in (("today", 5.9, 40.8, 78.4+11.4, 10.4), ("today_q4", 6.62, 26.4, None, 10.4)):
    P = 14.39626
    a, d, rm, rc = geo(P, Ma, Md, a)
    print(f"{lab:11s} {Md:5.2f} {Ma:5.1f} {Md/Ma:5.3f} {P:6.2f} {a:5.1f} {d:5.1f} {rm:5.2f} {rc:5.2f} {R:5.1f} {rm/R:6.2f} {rc/R:6.2f}")
for lab, L, M, X in (("today", 10**5.12, 40.8, 0.73), ("today Xs", 10**5.12, 40.8, 0.42)):
    print(lab, "Gamma_e =", round(10**-4.813*(1+X)*L/M, 3))
