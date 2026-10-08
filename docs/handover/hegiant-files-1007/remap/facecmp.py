#!/usr/bin/env python3
"""Face fluxes vs cell state: for restarts on the same grid, |f1_B - f1_A| / (c E) at x1 faces, and the
consistency of each file's stored f1 with its own cell F1 (|f1 - mean of the two cells' F1| / (c E)), by zone.
usage: facecmp.py A.rst B.rst"""
import re
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/he_giant_1006/remap')
import he_remap_rst_w256 as H
RS = 6.957e10
R = [H.read_rst(f) for f in sys.argv[1:3]]
G = [H.gather(r) for r in R]
r = R[0]
ng, n1 = r.ng, r.bind[1]
p = [float(H.get_param(r.text, 'mesh', k) or 0.0) for k in H.SKEYS]
fa = H.faces_of(p, n1, ng, r.msize[0], r.msize[3]) / RS
cl = float(re.search(r'<rad_m1>\n(?:(?!<).*\n)*?c_light\s*=\s*(\S+)', r.text).group(1))
zones = [(3, 30), (30, 60), (60, 65), (65, 70), (70, 200.1)]
for nm, g in zip(['A', 'B'], G):
    E, F1, f1 = g['m1'][0], g['m1'][1], g['f1']
    Ef = 0.5 * (E[..., 1:] + E[..., :-1])
    Fm = 0.5 * (F1[..., 1:] + F1[..., :-1])
    d = np.abs(f1[..., 1:-1] - Fm) / (cl * Ef)
    x = fa[1:-1]
    print(nm, 'stored f1 vs mean cell F1 (/cE):', ' | '.join('%g-%g med %.1e max %.1e' % (
        lo, hi, np.median(d[..., (x >= lo) & (x < hi)]), d[..., (x >= lo) & (x < hi)].max()) for lo, hi in zones))
E = G[0]['m1'][0]
Ef = 0.5 * (E[..., 1:] + E[..., :-1])
d = np.abs(G[1]['f1'][..., 1:-1] - G[0]['f1'][..., 1:-1]) / (cl * Ef)
x = fa[1:-1]
print('f1 B - A (/cE):', ' | '.join('%g-%g med %.1e max %.1e' % (
    lo, hi, np.median(d[..., (x >= lo) & (x < hi)]), d[..., (x >= lo) & (x < hi)].max()) for lo, hi in zones))
for k in ('f2', 'f3'):
    a, b = G[0][k], G[1][k]
    sc = cl * np.median(E)
    print(k, 'max |B - A| / (c median E) %.2e, max |A| / (c median E) %.2e' % (np.abs(b - a).max() / sc,
                                                                            np.abs(a).max() / sc))
