#!/usr/bin/env python3
"""tests_r23: the two numbers that judge an ic_rad_* arm.

  (1) hst drift over 1.5 turnovers vs the tests_r15 references
      (Z_m0 = the OLD MLT IC with the closure switched off, S0.3 = the same with it on)
  (2) the face luminosity profile from problem/e_ledger (e_ledger_shell.txt), normalised
      to the innermost face: a radiative-equilibrium IC must be flat at 1 inside 1 R.

usage: an23.py [arm ...]        (default: the arms built on 2026-09-20)
"""
import struct
import sys

import numpy as np

B = '/viper/u2/jinma/ATHENAK/bench/wt_he4'
RS, LSTAR = 2.3717e11, 2.3066e38
X1MIN, X1MAX = 9.486800e10, 2.964625e11


def rd(fn):
    out, f = [], open(fn, 'rb')
    while True:
        h = f.read(16)
        if len(h) < 16:
            break
        t, n1, nv = struct.unpack('<dii', h)
        b = f.read(8*(n1 + nv*n1))
        if len(b) < 8*(n1 + nv*n1):
            break
        a = np.frombuffer(b, '<f8')
        out.append((t, a[:n1].copy(), a[n1:].reshape(nv, n1).copy()))
    return out


def hst(runs):
    print('%-7s %8s %8s %12s %12s %10s %9s'
          % ('arm', 't_end', 'dt_end', 'mom1', 'KEr', 'dE/E0', 'M/M0'))
    for n, d in runs:
        a = np.loadtxt(d + '/he4.hydro.hst', comments='#')
        t, dt, m, p1, E, ke = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 6], a[:, 7]
        k = len(t) - 1
        print('%-7s %8.1f %8.3f %12.4e %12.4e %10.3e %9.6f'
              % (n, t[k], dt[k], p1[k], ke[k], E[k]/E[0] - 1.0, m[k]/m[0]))


def ledger(arms):
    for arm in arms:
        rc = rd(arm + '/rt_profile.bin')[0][1]
        rf = np.concatenate(([X1MIN], 0.5*(rc[1:] + rc[:-1]), [X1MAX]))
        fn = arm + '/e_ledger_shell.txt'
        t = float(open(fn).readline().split('=')[1].split()[0])
        lf = np.loadtxt(fn)[:, 2]/(LSTAR*t)
        lf = lf/lf[1]
        idx = [int(np.argmin(np.abs(rf/RS - s)))
               for s in (0.55, 0.65, 0.75, 0.85, 0.90, 0.95, 0.97, 1.00, 1.05)]
        m = (rf[:len(lf)]/RS > 0.51) & (rf[:len(lf)]/RS < 1.02)
        print('%-7s t=%7.1f s  L_face/L_in: %s | 0.51-1.02 R min %.4f max %.4f'
              % (arm, t, ' '.join('%.2f:%.4f' % (rf[i]/RS, lf[i]) for i in idx),
                 lf[m].min(), lf[m].max()))


if __name__ == '__main__':
    arms = sys.argv[1:] or ['D100', 'Dmlt', 'D100s']
    hst([('V100ref', B + '/tests_r23/V100'), ('Z_m0', B + '/tests_r15/Z_m0')]
        + [(a, a) for a in arms])
    print()
    ledger(arms)
