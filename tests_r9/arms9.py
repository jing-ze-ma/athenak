#!/usr/bin/env python3
"""tests_r9 task 2: the bisection table.  arms9.py <dir> [<dir> ...]"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prof9 import read_prof, TCGS, GM, RIN, RSTAR, TURN, taper_w

print("%-12s %8s %8s %9s %9s %9s %9s %9s %9s %8s %8s" % (
    "arm", "t/turn", "dt_end", "rho(1t)", "rho(2t)", "rho(3t)", "T(1t)", "T(3t)",
    "v1(last)", "dM/M", "collapse"))
for d in sys.argv[1:]:
    a = d.rstrip('/')
    hst = os.path.join(a, 'he4.hydro.hst')
    if not os.path.exists(hst):
        print("%-12s  (no hst)" % a); continue
    h = np.loadtxt(hst)
    if h.ndim == 1:
        h = h[None, :]
    tend, dtend = h[-1, 0], h[-1, 1]
    dm = h[-1, 2]/h[0, 2] - 1.0
    log = os.path.join(a, '..', os.path.basename(a).replace('_', '_') + '.log')
    coll = "?"
    for cand in [a + '.log', os.path.join(os.path.dirname(a) or '.', a + '.log')]:
        pass
    pf = os.path.join(a, 'rt_profile.bin')
    vals = {}
    v1last = float('nan')
    if os.path.exists(pf):
        recs = read_prof(pf)
        for want in (1.0, 2.0, 3.0):
            k = int(np.argmin([abs(r[0]-want*TURN) for r in recs]))
            if abs(recs[k][0]-want*TURN) > 0.12*TURN:
                continue
            t, x, q = recs[k]
            i = int(np.argmin(np.abs(x/RSTAR - 0.97)))
            vals[want] = (q[0][i], q[5][i]*TCGS, q[1][i])
        t, x, q = recs[-1]
        i = int(np.argmin(np.abs(x/RSTAR - 0.97)))
        v1last = q[1][i]
    def g(w, j, f="%9.3e"):
        return (f % vals[w][j]) if w in vals else "        -"
    print("%-12s %8.3f %8.3g %s %s %s %s %s %9.2e %8.1e" % (
        a, tend/TURN, dtend, g(1.0, 0), g(2.0, 0), g(3.0, 0),
        g(1.0, 1), g(3.0, 1), v1last, dm))
