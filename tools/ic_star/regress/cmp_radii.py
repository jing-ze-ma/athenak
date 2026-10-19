#!/usr/bin/env python3
"""cmp_radii.py NEW REF r1 r2 ... [Rsun]: ratios NEW/REF of rho, F_r, T (log-interpolated in r)
and both F_MLT/F at the given radii (two IC columns with different rows, e.g. the tool's mode
'mesa' 14p5M column vs hg145 ic3/ic_hg145_mesaP.txt)."""
import sys

import numpy as np

RSUN = 6.957e10


def main():
    a = np.loadtxt(sys.argv[1], comments='#')
    b = np.loadtxt(sys.argv[2], comments='#')
    print('# %s (%d rows) / %s (%d rows)' % (sys.argv[1], a.shape[0], sys.argv[2], b.shape[0]))
    print('# r[Rsun]  rho     F_r     T       F_MLT/F new  ref')
    for x in sys.argv[3:]:
        rr = float(x)*RSUN
        i = int(np.argmin(abs(a[:, 0] - rr)))
        ri = a[i, 0]

        def ib(j):
            return np.exp(np.interp(np.log(ri), np.log(b[:, 0]), np.log(np.maximum(b[:, j], 1e-300))))
        print('%8.3f  %.5f %.5f %.5f  %.4f %.4f' % (ri/RSUN, a[i, 1]/ib(1), a[i, 3]/ib(3), a[i, 5]/ib(5),
                                                   a[i, 6], np.interp(ri, b[:, 0], b[:, 6])))


if __name__ == '__main__':
    main()
