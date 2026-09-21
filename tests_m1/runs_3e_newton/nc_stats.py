"""Milestone 3e gate N-c: the seeded-slab numbers of one arm."""
import sys
import glob
import numpy as np

FIN = 2.475202e15


def one(d):
    h = np.loadtxt(d + '/m1slab.hydro.hst')
    u = np.loadtxt(d + '/m1slab.user.hst')
    print('==== %s  (%d hst rows, t_end=%.3f, dt_end=%.4g)'
          % (d, len(h), h[-1, 0], h[-1, 1]))
    print('  t      dt        KE_1        KE_2        F1top/Fin')
    for t in (1.0, 5.0, 10.0, 20.0, 30.0):
        if h[-1, 0] < t:
            continue
        r = int(np.argmin(np.abs(h[:, 0] - t)))
        print('  %5.1f  %8.4f  %10.3e  %10.3e  %8.4f'
              % (h[r, 0], h[r, 1], h[r, 7], h[r, 8], u[r, 2]/FIN))
    r = len(h) - 1
    print('  END   %8.4f  %10.3e  %10.3e  %8.4f'
          % (h[r, 1], h[r, 7], h[r, 8], u[r, 2]/FIN))
    nb = len(glob.glob(d + '/bin/*.m1.*.bin'))
    print('  m1 dumps: %d' % nb)


for a in sys.argv[1:]:
    one(a)
