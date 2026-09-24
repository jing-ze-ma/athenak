#!/usr/bin/env python3
"""T-sym: radial slices at several (theta, phi) cells of a spherically symmetric run
must agree.  usage: tsym.py dir basename num   (reads <dir>/tab/<base>.j*k*.<num>.tab
and the matching .00000 start for the amount of evolution)"""
import glob
import sys

import numpy as np

import s1lib

d, base, num = sys.argv[1:4]
fs = sorted(glob.glob(f"{d}/tab/{base}.j*k*.{num}.tab"))
es = np.array([s1lib.load_tab(f)['m1_e'] for f in fs])
f1 = np.array([s1lib.load_tab(f)['m1_f1'] for f in fs])
f2 = np.array([s1lib.load_tab(f)['m1_f2'] for f in fs])
f3 = np.array([s1lib.load_tab(f)['m1_f3'] for f in fs])
e0 = s1lib.load_tab(fs[0].replace(f'.{num}.', '.00000.'))['m1_e']
t = s1lib.load_tab(fs[0])['time']
spread = np.max((es.max(axis=0) - es.min(axis=0))/es.mean(axis=0))
fsp = np.max(np.abs(f1.max(axis=0) - f1.min(axis=0)))/np.max(np.abs(f1))
print(f"{d}: t={t:g} slices={len(fs)} max_i spread(E)/E={spread:.3e} "
      f"spread(F1)/max|F1|={fsp:.3e} "
      f"max|F2|/max|F1|={np.max(np.abs(f2))/np.max(np.abs(f1)):.3e} "
      f"max|F3|/max|F1|={np.max(np.abs(f3))/np.max(np.abs(f1)):.3e} "
      f"evolution max|E-E0|/max E0={np.max(np.abs(es[0]-e0))/np.max(e0):.3e}")
hs = sorted(glob.glob(f"{d}/tab/{base}.hj*k*.{num}.tab"))
if hs:
    h = [s1lib.load_tab(f) for f in hs]
    for v in ('dens', 'eint', 'velx'):
        a = np.array([x[v] for x in h])
        sc = np.max(np.abs(a)) if v == 'velx' else a.mean(axis=0)
        print(f"  gas {v}: max_i spread/{'max|v1|' if v == 'velx' else v}="
              f"{np.max((a.max(axis=0) - a.min(axis=0))/sc):.3e}")
    for v in ('vely', 'velz'):
        a = np.array([x[v] for x in h])
        vx = np.max(np.abs([x['velx'] for x in h]))
        print(f"  gas max|{v}|/max|v1|={np.max(np.abs(a))/vx:.3e}")
