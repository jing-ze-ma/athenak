#!/usr/bin/env python3
"""Flux saturation |F|/(cE) by radial zone for restarts (same or different grids): fraction of cells above
0.99 / 0.999 / 0.99999 and max.  usage: fsat.py A.rst [B.rst ...]"""
import re
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/he_giant_1006/remap')
import he_remap_rst_w256 as H
RS = 6.957e10
for fn in sys.argv[1:]:
    r = H.read_rst(fn)
    G = H.gather(r)
    ng, n1 = r.ng, r.bind[1]
    p = [float(H.get_param(r.text, 'mesh', k) or 0.0) for k in H.SKEYS]
    f = H.faces_of(p, n1, ng, r.msize[0], r.msize[3])
    xc = (0.5 * (f[1:] + f[:-1]) / RS)[ng:ng + n1]
    cl = float(re.search(r'<rad_m1>\n(?:(?!<).*\n)*?c_light\s*=\s*(\S+)', r.text).group(1))
    w = G['m1'][..., ng:ng + n1]
    fr = np.sqrt(w[1]**2 + w[2]**2 + w[3]**2) / (cl * w[0])
    print(fn.split('/')[-2], 'nx1', n1)
    for lo, hi in [(3, 60), (60, 65), (65, 70), (70, 100), (100, 150), (150, 200.1)]:
        m = (xc >= lo) & (xc < hi)
        q = fr[..., m]
        print('  r %5.0f-%5.0f: |F|/cE > 0.99 %.4f  > 0.999 %.4f  > 0.99999 %.4f  max %.6f  median %.4f' % (
            lo, hi, (q > 0.99).mean(), (q > 0.999).mean(), (q > 0.99999).mean(), q.max(), np.median(q)))
