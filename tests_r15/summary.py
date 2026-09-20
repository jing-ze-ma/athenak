#!/usr/bin/env python3
"""tests_r15: one table per quantity across the discriminator arms. Read-only.

usage: summary.py <arm:lo:hi> ...     e.g. summary.py G0:1e-11:1e-10 G2:3e-11:3e-10
"""
import os
import sys

import numpy as np

import cliff

TURN = 4705.0
RS = 2.3717e11
TT = (0.25, 0.5, 1.0, 1.5)


def load(spec):
    a = spec.split(':')
    arm, lo, hi = a[0], float(a[1]), float(a[2])
    P = cliff.rd(os.path.join(arm, 'rt_profile.bin'))
    H = np.loadtxt(os.path.join(arm, 'he4.hydro.hst'), comments='#')
    return arm, lo, hi, P, H


def at(P, tt):
    k = int(np.argmin([abs(x[0]-tt*TURN) for x in P]))
    return P[k]


def main():
    A = [load(s) for s in sys.argv[1:]]
    # per-time cell tables
    for tt in TT:
        print('\n########## t = %.2f turnover (%.0f s)' % (tt, tt*TURN))
        for nm, sl, fmt in (('rho', 0, '%10.3e'), ('T', 5, '%10.3e'),
                            ('v1', 1, '%10.3e')):
            print(' %-4s  i=' % nm + ''.join('%10d' % i for i in range(96, 103))
                  + '   max|v1|i>=96  max|v1|i>=98')
            for arm, lo, hi, P, H in A:
                t, r, q = at(P, tt)
                if abs(t - tt*TURN) > 0.4*TURN:
                    print('  %-5s  (no profile)' % arm)
                    continue
                v = q[sl]
                ex = ''
                if nm == 'v1':
                    ex = '  %11.3e %13.3e' % (np.abs(q[1][96:]).max(),
                                              np.abs(q[1][98:]).max())
                print('  %-5s   ' % arm + ''.join(fmt % v[i] for i in range(96, 103))
                      + ex)
        print(' w     i=' + ''.join('%10d' % i for i in range(96, 103)))
        for arm, lo, hi, P, H in A:
            t, r, q = at(P, tt)
            w = cliff.weight(q[0], lo, hi)
            print('  %-5s   ' % arm + ''.join('%10.4f' % w[i] for i in range(96, 103)))
    # hst table
    print('\n########## hst (nearest time)')
    print(' arm    t/turn      dt          mass          1-mom        1-KE'
          '        M(rho<rho_hi)   M(rho<rho_lo)')
    for arm, lo, hi, P, H in A:
        for tt in TT:
            t, r, q = at(P, tt)
            if abs(t - tt*TURN) > 0.4*TURN:
                continue
            j = int(np.argmin(np.abs(H[:, 0] - t)))
            rf = cliff.np.empty(len(r)+1)
            rf[1:-1] = 0.5*(r[1:] + r[:-1])
            rf[0] = r[0] - (rf[1]-r[0])
            rf[-1] = r[-1] + (r[-1]-rf[-2])
            vol = 4*np.pi/3*(rf[1:]**3 - rf[:-1]**3)
            rho = q[0]
            mhi = (rho*vol)[rho < hi].sum()
            mlo = (rho*vol)[rho < lo].sum()
            print(' %-5s  %6.3f  %10.4g %14.7e %12.4e %12.4e  %12.4e  %12.4e'
                  % (arm, t/TURN, H[j, 1], H[j, 2], H[j, 3], H[j, 7], mhi, mlo))
        print('')


if __name__ == '__main__':
    main()
