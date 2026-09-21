"""Summarise a runs_3d_edd arm: KE_1/KE_2, max|v1|, F1top/Fin, Picard, fallbacks."""
import sys
import re
import numpy as np

FIN = 2.475202e15
VMLT = 1.86e4
TT = [1.0, 5.0, 10.0, 30.0, 100.0, 200.0]


def at(t, x, y):
    return float(np.interp(t, x, y)) if t <= x[-1] else float('nan')


for d in sys.argv[1:]:
    h = np.loadtxt(d + '/m1slab.hydro.hst')
    u = np.loadtxt(d + '/m1slab.user.hst')
    print('==== %s   t_end = %.1f' % (d, h[-1, 0]))
    print('  %6s %11s %11s %11s %11s' % ('t', 'KE_1', 'KE_2', 'V1max', 'F1top/Fin'))
    for t in TT:
        if t > h[-1, 0] + 1e-9:
            continue
        print('  %6.0f %11.3e %11.3e %11.3e %11.4f'
              % (t, at(t, h[:, 0], h[:, 7]), at(t, h[:, 0], h[:, 8]),
                 at(t, u[:, 0], u[:, 5]), at(t, u[:, 0], u[:, 2])/FIN))
    print('  max over run: |v1| %.3e (= %.3f v_MLT), KE_1 %.3e, KE_2 %.3e'
          % (u[:, 5].max(), u[:, 5].max()/VMLT, h[:, 7].max(), h[:, 8].max()))
    print('  F1top/Fin range %.4f .. %.4f' % ((u[:, 2]/FIN).min(), (u[:, 2]/FIN).max()))
    for ln in open(d + '/log.txt', errors='replace'):
        if re.search(r'implicit transport:|positivity fallbacks|zone-cycles', ln):
            print('  | ' + ln.strip())
