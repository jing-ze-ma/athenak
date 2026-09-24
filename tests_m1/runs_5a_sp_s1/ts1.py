#!/usr/bin/env python3
"""T-S1 steady spherical diffusion: E(r) vs E_a = e_out + b (1/r - 1/r_ie),
b = 3 k F_in r_in^2 / c, and the face luminosity r_f^2 F_f (rebuilt from the cell
means F1_i = (F_l + F_r)/2 starting from the imposed F_in) against r_in^2 F_in.

usage: ts1.py k F_in e_out c r_in r_out [poly|none] dump1 [dump2 ...]
"""
import sys

import numpy as np

import s1lib

k, fin, eout, c, rin, rout = map(float, sys.argv[1:7])
poly = None if sys.argv[7] == 'none' else s1lib.HE4_POLY
res = []
for fn in sys.argv[8:]:
    d = s1lib.load_tab(fn)
    e = d['m1_e']
    n = e.shape[0]
    rf, x1v = s1lib.rgrid(rin, rout, n, poly)
    b = 3.0*k*fin*rin*rin/c
    ea = eout + b*(1.0/x1v - 1.0/x1v[-1])
    em = e
    l1 = np.mean(np.abs(em - ea)/ea)
    linf = np.max(np.abs(em - ea)/ea)
    f1 = d['m1_f1']
    ff = np.zeros(n + 1)
    ff[0] = fin
    for i in range(n):
        ff[i+1] = 2.0*f1[i] - ff[i]
    lum = rf*rf*ff/(rin*rin*fin)
    ldev = np.max(np.abs(lum[:-1] - 1.0))   # the outer (efix) face carries no flux eq.
    t2 = np.max(np.abs(d['m1_f2']))/np.max(np.abs(f1))
    t3 = np.max(np.abs(d['m1_f3']))/np.max(np.abs(f1))
    spread = 0.0
    res.append((n, l1, linf))
    print(f"{fn}: n={n} L1(E)={l1:.4e} Linf={linf:.4e} "
          f"max|r^2F/(rin^2Fin)-1|={ldev:.3e} shell spread={spread:.2e} "
          f"|F2|/|F1|={t2:.1e} |F3|/|F1|={t3:.1e}")
for a, bb in zip(res[:-1], res[1:]):
    print(f"order n={a[0]}->{bb[0]}: L1 {np.log(a[1]/bb[1])/np.log(bb[0]/a[0]):.3f} "
          f"Linf {np.log(a[2]/bb[2])/np.log(bb[0]/a[0]):.3f}")
