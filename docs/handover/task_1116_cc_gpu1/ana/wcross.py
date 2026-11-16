#!/usr/bin/env python3
"""w at beam crossings vs single-beam cells from an M1_AKDBG dump. usage: wcross.py ANG RUNDIR..."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cb as C
ang = float(sys.argv[1])
for d in sys.argv[2:]:
    p = C.readin(d + '/in.athinput')
    p['rad_m1/implicit_beam0_ang'] = str(-ang); p['rad_m1/implicit_beam1_ang'] = str(ang)
    w = np.loadtxt(d + '/w_r0.txt')
    x1 = np.unique(w[:, 0]); x2 = np.unique(w[:, 1])
    e0, e1 = C.exact(p, x1, x2, each=True, ns=8, nsy=8)
    W1 = w[:, 2].reshape(len(x2), len(x1)); W2 = w[:, 3].reshape(len(x2), len(x1))
    both = (e0 > 0.5) & (e1 > 0.5); one = ((e0 > 0.5) ^ (e1 > 0.5)); none = (e0 + e1) < 1e-6
    s = '  '.join('%s n=%d w1 %.3f w2 %.3f' % (nm, m.sum(), W1[m].mean() if m.any() else np.nan,
                  W2[m].mean() if m.any() else np.nan) for nm, m in (('cross', both), ('single', one), ('vac', none)))
    print(d.split('/p2/')[-1], s)
