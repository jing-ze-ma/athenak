#!/usr/bin/env python3
"""tests_r9: the ENERGY AND MASS BUDGET of the shells around the pile-up, from the run's
OWN t=0 face dump (mltfaces_*.txt: F_req, F_2s, F_cond = w_blend F_raddiff, F_used = the
MLT closure's applied flux, F_conv_res) plus the profile record's measured d/dt.
Usage: budget9.py <mltfaces.txt> <rt_profile.bin> [--raw]"""
import sys
import numpy as np
sys.path.insert(0, '.')
from prof9 import read_prof, TURN, TCGS, GM, RIN, RSTAR, LSTAR

mf = np.loadtxt(sys.argv[1], comments='#')
recs = read_prof(sys.argv[2])
raw = '--raw' in sys.argv
traw = raw or ('--traw' in sys.argv)
rf_c, freq, f2s, fcond, fused, fres = (mf[:, 1], mf[:, 10], mf[:, 17],
                                       mf[:, 8]*mf[:, 12], mf[:, 9], mf[:, 14])
pf = mf[:, 3]

print("t = 0 FACE BUDGET (the run's own dump), fluxes / F_req = L/(4 pi r^2):")
print(" r/R      F_req[cgs]   F_2s/F_req  F_cond/F_req  F_MLTused/F_req  F_res/F_req"
      "   SUM   L_face/L")
for s in [0.90, 0.94, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01]:
    k = int(np.argmin(np.abs(rf_c/RSTAR - s)))
    tot = (f2s[k]+fcond[k]+fused[k]+fres[k])/freq[k]
    print(" %6.4f %.4e %11.4f %13.4f %16.4f %12.3e %7.4f %8.4f"
          % (rf_c[k]/RSTAR, freq[k], f2s[k]/freq[k], fcond[k]/freq[k],
             fused[k]/freq[k], fres[k]/freq[k], tot,
             4*np.pi*rf_c[k]**2*(f2s[k]+fcond[k]+fused[k])/LSTAR))

# the divergence of the TOTAL carried flux across each shell at t = 0, per unit volume
print("\nt = 0 SHELL DIVERGENCE of the carried flux (negative = the shell is HEATED):")
print(" cell r/R   rf_in/R    rf_out/R   -div F [erg/cm3/s]  (L_in - L_out)[erg/s]"
      "   /L")
fc = f2s + fcond + fused
for s in [0.94, 0.96, 0.97, 0.98, 0.99]:
    k = int(np.argmin(np.abs(rf_c/RSTAR - s)))
    a = k; b = k+1
    La = 4*np.pi*rf_c[a]**2*fc[a]
    Lb = 4*np.pi*rf_c[b]**2*fc[b]
    V = (4*np.pi/3.0)*(rf_c[b]**3 - rf_c[a]**3)
    print(" %8.4f  %8.4f  %8.4f   %+.4e        %+.4e   %+.4f"
          % (0.5*(rf_c[a]+rf_c[b])/RSTAR, rf_c[a]/RSTAR, rf_c[b]/RSTAR,
             (La-Lb)/V, La-Lb, (La-Lb)/LSTAR))

# the MEASURED evolution of the same shells, off the profile record
print("\nMEASURED, per 0.25 turnover, for the cells at 0.96/0.97/0.98/0.99/1.00 R:")
want = np.arange(0.0, 3.001, 0.25)*TURN
sel = []
for tw in want:
    k = int(np.argmin([abs(rr[0]-tw) for rr in recs]))
    if k not in sel and abs(recs[k][0]-tw) < 0.13*TURN:
        sel.append(k)
prev = {}
for k in sel:
    t, x1v, q = recs[k]
    rho = q[0]; v1 = q[1]; rv1 = q[2]
    T = q[5]*TCGS if traw else q[5]
    eint = q[6] - (rho*GM*(1.0/RIN - 1.0/x1v) if raw else 0.0)
    out = ["t/turn %5.3f" % (t/TURN)]
    for s in [0.96, 0.97, 0.98, 0.99, 1.00]:
        i = int(np.argmin(np.abs(x1v/RSTAR - s)))
        de = dr_ = float('nan')
        if s in prev:
            t0, e0, r0 = prev[s]
            de = (eint[i]-e0)/(t-t0); dr_ = (rho[i]-r0)/(t-t0)
        prev[s] = (t, eint[i], rho[i])
        out.append("| %4.2fR rho %.3e (%+.2e/s) eint %.3e (%+.3e/s) v1 %+.2e"
                   % (s, rho[i], dr_, eint[i], de, v1[i]))
    print("  ".join(out))
