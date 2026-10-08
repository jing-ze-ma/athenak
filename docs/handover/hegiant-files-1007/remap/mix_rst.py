#!/usr/bin/env python3
"""Diagnostic: write a restart on grid A whose gathered arrays come from file A except the listed keys (and
optionally only for r > RMIN), which come from file B (same grid).  Caches of A are kept (identity writer, shown
bitwise for A itself).  usage: mix_rst.py A.rst B.rst OUT.rst key[,key..] [RMIN_Rsun]
keys: hyd m1 f1 f2 f3 ctr (or hyd0..hyd4 / m10..m13 for single components)"""
import re
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/he_giant_1006/remap')
import he_remap_giant as H
RS = 6.957e10
ra, rb = H.read_rst(sys.argv[1]), H.read_rst(sys.argv[2])
Ga, Gb = H.gather(ra), H.gather(rb)
keys = sys.argv[4].split(',')
rmin = float(sys.argv[5]) * RS if len(sys.argv) > 5 else 0.0
ng, n1 = ra.ng, ra.bind[1]
p = [float(H.get_param(ra.text, 'mesh', k) or 0.0) for k in H.SKEYS]
f = H.faces_of(p, n1, ng, ra.msize[0], ra.msize[3])
xc = 0.5 * (f[1:] + f[:-1])
for k in keys:
    m = re.match(r'(hyd|m1)(\d)$', k)
    base, comp = (m.group(1), int(m.group(2))) if m else (k, None)
    A, B = np.array(Ga[base]), np.asarray(Gb[base])
    x = xc if A.shape[-1] == xc.size else f
    sel = x >= rmin
    if comp is None:
        A[..., sel] = B[..., sel]
    else:
        A[comp][..., sel] = B[comp][..., sel]
    Ga[base] = A
cl = float(re.search(r'<rad_m1>\n(?:(?!<).*\n)*?c_light\s*=\s*(\S+)', ra.text).group(1))
src_mb = {(int(ra.lloc[m, 1]), int(ra.lloc[m, 2])): m for m in range(ra.nmb)}
out, st = H.remap(ra, Ga, 1, cl)
H.write_rst(ra, out, sys.argv[3], 1, ra.bind[2], ra.bind[3], True, 'source', True, src_mb)
print('wrote', sys.argv[3], 'with', keys, 'from B for r >= %.1f Rsun' % (rmin / RS))
