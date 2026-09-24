#!/usr/bin/env python3
"""T-S4: Richardson extrapolation 2 E_2n - E_n (E_n interpolated to the 2n centres) of
each closure vs the exact-transfer reference of ts4v.py (ts4_exact.npz); removes the
first-order Marshak end-cell error common to all closures.
usage: ts4rich.py <root> <prefix-closure> <n> [...]   e.g. q4v 128"""
import glob
import os
import sys

import numpy as np
from scipy.interpolate import CubicSpline

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'runs_5b_sp_s2'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'vis', 'python'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import s2lib  # noqa: E402

arad, rin, rout = 1e8, 1.0, 5.0
here = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(here, 'ts4_exact.npz'))
ref = CubicSpline(d['r'], d['J'])
root = sys.argv[1]
for pre, n in zip(sys.argv[2::2], sys.argv[3::2]):
    n = int(n)
    ee = []
    for m in (n, 2*n):
        fm = sorted(glob.glob(f'{root}/{pre}{m}/tab/*.m1.*.tab'))[-1]
        e = s2lib.load_tab(fm)['m1_e']
        ee.append((s2lib.rgrid(rin, rout, m)[1], e))
    x2, e2 = ee[1]
    e1 = CubicSpline(ee[0][0], np.log(ee[0][1]))(x2)
    er = 2.0*e2 - np.exp(e1)
    sel = (x2 > ee[0][0][0]) & (x2 < ee[0][0][-1])
    tr = (ref(x2[sel])/arad)**0.25
    dv = np.abs((er[sel]/arad)**0.25/tr - 1.0)
    d2 = np.abs((e2[sel]/arad)**0.25/tr - 1.0)
    print(f"{pre} n={n},{2*n}: Richardson L1 {np.mean(dv):.3e} Linf {np.max(dv):.3e}"
          f"   (n={2*n} alone: L1 {np.mean(d2):.3e} Linf {np.max(d2):.3e})")
