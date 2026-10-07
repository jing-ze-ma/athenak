#!/usr/bin/env python3
"""ESTIMATE of cost per orbit, scaled from RY Per stage 1 (IMPL.md: 448x4x2048, dlnr=dphi=3.07e-3,
r_in=4.06, r_out=0.85 dL1, 1.29e5 cycles/orbit, 23 min/orbit/apu node, dt set by the stream
speed v_ff at r_in over r_in*dphi). Same cell shape; cycles/orbit ~ P*v_ff(r_in)/r_in; cells ~ nr.
Also: stream width cs/Omega in cells at r_out (resolution check) and inner Kepler/orbit ratio."""
import numpy as np
from scipy.optimize import brentq
G, Msun, Rsun = 6.674e-8, 1.989e33, 6.957e10
D = 2*np.pi/2048
def geo(P, Ma, Md, a):
    mu = Md/(Ma+Md); xa, xd = -mu, 1-mu
    gx = lambda x: -(1-mu)*(x-xa)/abs(x-xa)**3 - mu*(x-xd)/abs(x-xd)**3 + x
    return (brentq(gx, xa+1e-4, xd-1e-4)-xa)*a
def metric(P, Ma, Md, a, rin, Tdon):
    dL1 = geo(P, Ma, Md, a); rout = 0.85*dL1
    vff = np.sqrt(2*G*Ma*Msun*(1/rin - 1/dL1)/Rsun)/1e5
    nr = np.log(rout/rin)/D
    X = P*vff/rin
    Om = 2*np.pi/(P*86400); cs = np.sqrt(1.381e-16*Tdon/(1.27*1.6726e-24))
    w = cs/Om/Rsun
    return nr, X, w/(rout*D), dL1
ref = metric(6.8636, 6.24, 1.69, 30.3, 4.06, 6250)
sysl = [  # name P Macc Mdon a Racc rmin(ballistic, stream_runs.txt) Tdon
 ("RY_Per", 6.8636, 6.24, 1.69, 30.3, 4.06, 2.69, 6250),
 ("RX_Cas", 32.33, 5.8, 1.8, 84.0, 2.5, 7.07, 4400),
 ("SX_Cas", 36.6, 5.1, 1.5, 87.0, 3.0, 7.53, 4000),
 ("W_Cru", 198.5, 7.8, 1.2, 297.9, 4.0, 36.52, 5500),
 ("bet_Lyr", 12.94, 13.2, 3.0, 58.7, 6.0, 5.75, 13200),
 ("RS_Cep", 12.42, 2.8, 0.4, 33.3, 2.65, 4.09, 4610),
 ("HD170582", 16.871, 9.0, 1.9, 61.4, 5.5, 6.23, 8000),
 ("AU_Mon", 11.113, 7.0, 1.2, 42.3, 5.1, 4.72, 5750)]
print(f"{'name':9s} {'r_in':>6s} {'choice':>9s} {'nr':>5s} {'cyc/orb':>8s} {'min/orb':>8s} {'strm_cells':>10s}")
for n, P, Ma, Md, a, R, rmin, T in sysl:
    for lab, rin in (("R_acc", R), ("0.6rmin", 0.6*rmin)):
        if lab == "0.6rmin" and 0.6*rmin <= R: continue
        nr, X, wc, dL1 = metric(P, Ma, Md, a, rin, T)
        cyc = 1.29e5*X/ref[1]
        cost = 23*(nr/448)*(cyc/1.29e5)
        print(f"{n:9s} {rin:6.2f} {lab:>9s} {nr:5.0f} {cyc:8.2e} {cost:8.0f} {wc:10.1f}")
