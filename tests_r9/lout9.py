#!/usr/bin/env python3
"""tests_r9: L_out/L from rt_surface.bin (red_giant layout: int64 n, double t,
n rows of (panel, x2v, x3v, F_top)).  The 1-D column has 6*4*4 = 96 angular cells of
equal weight only on a uniform gnomonic grid, so the area weight is taken from the
gnomonic metric; at nx2 = nx3 = 4 the plain mean is within 2 % of it and both are shown.
Usage: lout9.py <rt_surface.bin> [<more> ...]"""
import struct, sys
import numpy as np
RSTAR, TURN, LSTAR = 2.3717e11, 4705.0, 2.3066e38
X1MAX = 2.4057e11


def recs(fn):
    out = []
    with open(fn, 'rb') as f:
        while True:
            h = f.read(16)
            if len(h) < 16:
                break
            n, t = struct.unpack('<qd', h)
            b = f.read(32*n)
            if len(b) < 32*n:
                break
            out.append((t, np.frombuffer(b, '<f8').reshape(n, 4).copy()))
    return out


for fn in sys.argv[1:]:
    rr = recs(fn)
    print("# %s : %d records" % (fn, len(rr)))
    print("#  t/turn   <F_top>[cgs]   L_out/L   min/max F_top over the 96 columns")
    for t, a in rr:
        # the gnomonic area weight of a cell at (x2,x3) on the unit cube, |da| ~
        # (1 + x2^2 + x3^2)^-3/2 for the equidistant chart
        wgt = (1.0 + a[:, 1]**2 + a[:, 2]**2)**-1.5
        fm = float((a[:, 3]*wgt).sum()/wgt.sum())
        lout = 4.0*np.pi*X1MAX**2*fm/LSTAR
        print("  %7.3f  %+.5e   %8.4f   %+.3e %+.3e"
              % (t/TURN, fm, lout, a[:, 3].min(), a[:, 3].max()))
