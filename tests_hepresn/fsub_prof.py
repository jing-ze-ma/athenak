#!/usr/bin/env python3
"""Adaptive MLT closure profiles (problem/mlt_closure_dout): <run>/hepresn.fsub.txt.
usage: fsub_prof.py RUNDIR [--tt 4700] [--shells]
Per output: FeCZ (F_MLT/F_IC > 0.02 faces) median and range of F_sub/F_MLT, median L_sub,
L_rad, L_conv (units of L_in) and max |L_rad + L_sub + L_conv - 1|; --shells prints the
per-face profile of the last output every 8th face."""
import sys

import numpy as np


def main():
    run = sys.argv[1]
    tt = float(sys.argv[sys.argv.index('--tt') + 1]) if '--tt' in sys.argv else 4700.0
    d = {}
    for ln in open(run + '/hepresn.fsub.txt'):
        r = ln.split()
        d.setdefault(r[1], {})[float(r[0])] = np.array(r[2:], float)
    fm = list(d['fmlt'].values())[0]
    fz = fm > 0.02 * fm.max() / 0.159
    print('  t/tt   Fsub/FMLT med  [min  max]   Lsub   Lrad   Lconv  max|sum-1|')
    for t in sorted(d['fsub']):
        fs, lr, lc = d['fsub'][t], d['lrad'][t], d['lconv'][t]
        q = fs[fz] / fm[fz]
        print('%6.2f  %8.3f  [%6.3f %6.3f] %6.3f %6.3f %+7.4f  %.2e'
              % (t / tt, np.median(q), q.min(), q.max(), np.median(fs[fz]),
                 np.median(lr[fz]), np.median(lc[fz]),
                 np.abs(fs + lr + lc - 1)[1:-1].max()))
    if '--shells' in sys.argv:
        t = max(d['fsub'])
        fs, lr, lc = d['fsub'][t], d['lrad'][t], d['lconv'][t]
        print('face  FMLT   Fsub   Lrad   Lconv  (t/tt %.2f)' % (t / tt))
        for i in range(0, len(fs), 8):
            print('%4d %6.3f %6.3f %6.3f %+7.4f' % (i, fm[i], fs[i], lr[i], lc[i]))


if __name__ == '__main__':
    main()
