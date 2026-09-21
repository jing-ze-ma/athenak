#!/usr/bin/env python3
"""Milestone 3a2, LIMIT 3: read the user history of the 1-D He column run and report
the bottom-boundary diagnostics -- max |v1| per 200 s bin in units of v_MLT, the
emergent/imposed flux ratios, the residual force and the column kinetic energy.

Usage: python3 he_v1.py <run dir>/feczrt.user.hst [--vmlt 1.86e4] [--fin 2.475202e15]
"""
import argparse
import json
import sys

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('hst')
    ap.add_argument('--vmlt', type=float, default=1.86e4)
    ap.add_argument('--fin', type=float, default=2.475202e15)
    ap.add_argument('--json', default=None)
    ap.add_argument('--label', default='run')
    args = ap.parse_args()
    d = np.loadtxt(args.hst)
    # columns: time, dt, F1top, F1mid, F1bot, V1max, Etot, Fres, V1mid, KEcol
    t = d[:, 0]
    f1top, f1mid, f1bot = d[:, 2], d[:, 3], d[:, 4]
    v1max, fres, v1mid, kecol = d[:, 5], d[:, 7], d[:, 8], d[:, 9]
    print('%-10s %-12s %-12s' % ('bin [s]', 'max|v1|', '/v_MLT'))
    for lo in range(0, int(t.max()), 200):
        m = (t >= lo) & (t < lo + 200)
        if not m.any():
            continue
        print('%-10s %-12.3e %-12.2f' % ('%d-%d' % (lo, lo + 200),
                                         v1max[m].max(), v1max[m].max()/args.vmlt))
    late = t >= 0.5*t.max()
    print('late half: max|v1| %.3e = %.2f v_MLT   median %.3e = %.2f v_MLT'
          % (v1max[late].max(), v1max[late].max()/args.vmlt,
             np.median(v1max[late]), np.median(v1max[late])/args.vmlt))
    print('F1/F_in   top %.6f  mid %.6f  bot %.6f   (last dump)'
          % (f1top[-1]/args.fin, f1mid[-1]/args.fin, f1bot[-1]/args.fin))
    print('F1bot/F_in over the late half: min %.6f max %.6f  |dev| max %.2e'
          % ((f1bot[late]/args.fin).min(), (f1bot[late]/args.fin).max(),
             np.abs(f1bot[late]/args.fin - 1.0).max()))
    print('residual force (last) %.3e    KEcol first/last %.3e %.3e'
          % (fres[-1], kecol[0], kecol[-1]))
    if args.json:
        n = max(1, len(t)//600)
        json.dump({'label': args.label, 't': t[::n].tolist(),
                   'v1max': v1max[::n].tolist(), 'v1mid': v1mid[::n].tolist(),
                   'f1bot_over_fin': (f1bot[::n]/args.fin).tolist(),
                   'kecol': kecol[::n].tolist(), 'vmlt': args.vmlt},
                  open(args.json, 'w'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
