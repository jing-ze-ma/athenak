#!/usr/bin/env python3
"""T-S2: free-streaming point source (kappa = 0, M1 closure) on the sp wedge.

Imposed F_in at r_in (implicit_bc_x1min = flux), vacuum Marshak outflow F = c E_ie
(marshak_q = 1) at r_out.  At the steady state the integrating-factor face gives
r^2 E = const to round-off and the E row gives r_f^2 F_f = r_in^2 F_in (telescoping).
F = c E then holds with the face value E_f = (r^2 E)/r_f^2 up to the Marshak end-cell
offset, F_f/(c E_f) = r_in^2 F_in/(c r^2 E), which converges at first order in dr/r.
(The stored cell F is limited to |F| <= c E, so F/(cE) is not read from the dump.)

usage: ts2.py c F_in r_in r_out poly|none dump [dump ...]
"""
import sys

import numpy as np

import s2lib

c, fin, rin, rout = map(float, sys.argv[1:5])
poly = None if sys.argv[5] == 'none' else s2lib.HE4_POLY
res = []
for fn in sys.argv[6:]:
    d = s2lib.load_tab(fn)
    e = d['m1_e']
    n = e.shape[0]
    rf, x1v = s2lib.rgrid(rin, rout, n, poly)
    q = x1v*x1v*e
    dq = np.max(np.abs(q/q.mean() - 1.0))
    ratio = rin*rin*fin/(c*q.mean())
    t2 = np.max(np.abs(d['m1_f2']))/np.max(np.abs(d['m1_f1']))
    t3 = np.max(np.abs(d['m1_f3']))/np.max(np.abs(d['m1_f1']))
    res.append((n, abs(ratio - 1.0)))
    print(f"{fn}: n={n} t={d['time']:.4g} max|r^2E/<r^2E>-1|={dq:.3e} "
          f"F/(cE)-1={ratio - 1.0:.4e} (r_out/r_ie)^2-1={(rout/x1v[-1])**2 - 1:.4e} "
          f"|F2|/|F1|={t2:.1e} |F3|/|F1|={t3:.1e}")
for a, b in zip(res[:-1], res[1:]):
    print(f"order n={a[0]}->{b[0]}: F/cE-1 {np.log(a[1]/b[1])/np.log(b[0]/a[0]):.3f}")
