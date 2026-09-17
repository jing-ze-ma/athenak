#!/usr/bin/env python3
"""Read red_giant's rt_profile.bin (layout: src/pgen/red_giant.cpp:494) and print T(r,t).

Usage: prof.py <rt_profile.bin> [i0 i1]
"""
import sys
import numpy as np

SLOT = ['rho', 'v1', 'rhov1', 'v1sq', 'vhsq', 'T', 'eint', 'v1h']


def read(fn):
    recs = []
    with open(fn, 'rb') as f:
        while True:
            hb = f.read(8)
            if len(hb) < 8:
                break
            t = np.frombuffer(hb, '<f8')[0]
            nx1, nvar = np.frombuffer(f.read(8), '<i4')
            x1v = np.frombuffer(f.read(8*nx1), '<f8')
            q = np.frombuffer(f.read(8*nvar*nx1), '<f8').reshape(nvar, nx1)
            recs.append((t, x1v, q))
    return recs


def main():
    recs = read(sys.argv[1])
    t = np.array([r[0] for r in recs])
    x1v = recs[0][1]
    Q = np.array([r[2] for r in recs])          # (nt, nvar, nx1)
    print('# %d records, t = %.1f .. %.1f, nx1 = %d' % (len(recs), t[0], t[-1],
                                                        len(x1v)))
    R = 2.3717e11
    T = Q[:, 5, :]
    # the hottest cell of each record, and its radius
    print('# %-9s %-5s %-10s %-10s %-10s %-10s' %
          ('t', 'imax', 'r/R', 'Tmax', 'T(0.79R)', 'T(top)'))
    i79 = int(np.argmin(abs(x1v/R - 0.791)))
    for n in range(0, len(recs)):
        im = int(np.argmax(T[n]))
        print('%-11.1f %-5d %-10.4f %-10.4g %-10.4g %-10.4g' %
              (t[n], im+3, x1v[im]/R, T[n, im], T[n, i79], T[n, -1]))


if __name__ == '__main__':
    main()
