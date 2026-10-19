#!/usr/bin/env python3
"""Ballistic L1 stream in the restricted 3-body problem (co-rotating frame).
Units: a=1, Omega=1, G(M1+M2)=1. Accretor (mass fraction mu_a) at x=0 after shift.
Prints ESTIMATES for DESIGN.md. Usage: roche_stream.py Macc Mdon a_Rsun P_day Racc_Rsun Tdon"""
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

G, Msun, Rsun, kB, mH = 6.674e-8, 1.989e33, 6.957e10, 1.381e-16, 1.6726e-24
Ma, Md, a_rs, Pd, Ra_rs, Tdon = map(float, sys.argv[1:7])
q = Md / Ma                      # donor / accretor
mu = Md / (Ma + Md)              # donor mass fraction
xa, xd = -mu, 1 - mu             # accretor, donor positions (COM origin)
Om = 2 * np.pi / (Pd * 86400.0)
a = a_rs * Rsun

def gx(x):  # x-acceleration on the line of centres (y=0), corotating
    return (-(1 - mu) * (x - xa) / abs(x - xa)**3 - mu * (x - xd) / abs(x - xd)**3 + x)

xL1 = brentq(gx, xa + 1e-3, xd - 1e-3)
dL1 = xL1 - xa                   # distance accretor -> L1 in a

def rhs(t, s):
    x, y, vx, vy = s
    r1 = np.hypot(x - xa, y); r2 = np.hypot(x - xd, y)
    ax = -(1 - mu) * (x - xa) / r1**3 - mu * (x - xd) / r2**3 + x + 2 * vy
    ay = -(1 - mu) * y / r1**3 - mu * y / r2**3 + y - 2 * vx
    return [vx, vy, ax, ay]

cs_d = np.sqrt(kB * Tdon / (1.27 * mH))      # neutral-ish donor photosphere
eps = cs_d / (a * Om)
# launch from L1 toward the accretor with speed eps (Flannery-type ballistic start)
s0 = [xL1 - 1e-4, 0.0, -eps, 0.0]
ev_r = lambda t, s: np.hypot(s[0] - xa, s[1]) - 1e-3
ev_r.terminal = True
sol = solve_ivp(rhs, [0, 20], s0, rtol=1e-10, atol=1e-12, dense_output=True, events=ev_r,
                max_step=1e-3)
x, y, vx, vy = sol.y
r = np.hypot(x - xa, y)
# first periastron = r_min (Lubow-Shu varpi_min)
i_min = np.argmax(np.diff(np.sign(np.diff(r))) > 0) + 1
rmin = r[i_min]
Ra = Ra_rs / a_rs
print(f"q=Md/Ma={q:.4f}  mu_don={mu:.4f}  Omega={Om:.4e} s^-1  a={a:.4e} cm  a*Om={a*Om/1e5:.1f} km/s")
print(f"L1 distance from accretor = {dL1:.4f} a = {dL1*a_rs:.2f} Rsun")
qa = Ma / Md
RL = 0.49 * qa**(2/3) / (0.6 * qa**(2/3) + np.log(1 + qa**(1/3)))
print(f"Eggleton RL(accretor) = {RL:.4f} a = {RL*a_rs:.2f} Rsun; donor RL = "
      f"{0.49*q**(2/3)/(0.6*q**(2/3)+np.log(1+q**(1/3)))*a_rs:.2f} Rsun")
print(f"cs_donor(T={Tdon:.0f},mu=1.27)={cs_d/1e5:.2f} km/s eps=cs/(a Om)={eps:.4f} "
      f"cs/Om={cs_d/Om/Rsun:.2f} Rsun")
print(f"ballistic r_min = {rmin:.4f} a = {rmin*a_rs:.2f} Rsun ; LS fit 0.0488 q^-0.464 = "
      f"{0.0488*q**-0.464:.4f} a; r_circ fit 0.0859 q^-0.426 = {0.0859*q**-0.426:.4f} a = "
      f"{0.0859*q**-0.426*a_rs:.2f} Rsun")
print(f"R_acc/a = {Ra:.4f}; direct impact: {Ra > rmin}")
GMa = G * Ma * Msun
for lab, rr in [("impact R_acc", Ra), ("r_out=0.9 dL1", 0.9 * dL1), ("r_out=0.8 dL1", 0.8 * dL1),
                ("r_out=0.7 dL1", 0.7 * dL1)]:
    idx = np.where(r[:i_min + 1] <= rr)[0]
    if len(idx) == 0:
        print(f"{lab}: stream does not reach r={rr:.3f} a"); continue
    k = idx[0]
    dx, dy = x[k] - xa, y[k] - ya if (ya := 0.0) == 0.0 else 0
    # inertial-frame velocity relative to accretor: v_rot + Om x r (accretor-centred)
    vxi = vx[k] - y[k]; vyi = vy[k] + (x[k])
    # remove accretor's own orbital velocity (accretor at (xa,0) moves with vy = xa)
    vyi -= xa
    rx, ry = dx / r[k], dy / r[k]
    vr_rot = vx[k] * rx + vy[k] * ry; vt_rot = -vx[k] * ry + vy[k] * rx
    vr_in = vxi * rx + vyi * ry; vt_in = -vxi * ry + vyi * rx
    phi = np.degrees(np.arctan2(dy, dx))  # 0 = towards donor (+x)
    psi = np.degrees(np.arctan2(abs(vt_rot), -vr_rot))
    jkep = np.sqrt(GMa * r[k] * a)
    print(f"{lab}: r={r[k]:.4f} a={r[k]*a_rs:.2f} Rsun phi={phi:.1f} deg (0=donor, + = "
          f"leading y) | rot frame v_r={vr_rot*a*Om/1e5:.1f} v_t={vt_rot*a*Om/1e5:.1f} km/s "
          f"angle from inward radial={psi:.1f} deg | inertial(acc-centred) v_r={vr_in*a*Om/1e5:.1f}"
          f" v_t={vt_in*a*Om/1e5:.1f} km/s | j_in/j_kep={vt_in*a*Om*r[k]*a/jkep:.3f}"
          f" | v_ff={np.sqrt(2*GMa/(r[k]*a))/1e5:.0f} km/s")
vcrit = np.sqrt(GMa / (Ra * a))
print(f"v_crit(R_acc)=sqrt(GM/R)={vcrit/1e5:.0f} km/s ; v_sync(R_acc)={Om*Ra*a/1e5:.1f} km/s")
