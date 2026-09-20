#!/usr/bin/env python3
"""tests_r15: cliff diagnostic for the G0..G3 discriminator arms. Read-only.

usage: cliff.py <arm> [lo hi]      (window defaults 1e-11 / 1e-10)
prints, at t = 0.25/0.5/1.0/1.5 turnovers: v1, rho, T, w at i = 96..102,
max |v1| over i >= 96, and hst dt/mass/1-mom/1-KE at the nearest time.
"""
import os
import struct
import sys

import numpy as np

TURN = 4705.0
RS = 2.3717e11


def rd(fn):
    out = []
    f = open(fn, 'rb')
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
    f.close()
    return out


def weight(rho, lo, hi):
    s = np.clip((np.log10(np.maximum(rho, 1e-300)) - np.log10(lo))
                / (np.log10(hi) - np.log10(lo)), 0.0, 1.0)
    return s*s*(3.0 - 2.0*s)


def main():
    arm = sys.argv[1]
    lo = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-11
    hi = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-10
    P = rd(os.path.join(arm, 'rt_profile.bin'))
    print('== %s : %d profiles, t = %.1f .. %.1f s (%.3f turn)'
          % (arm, len(P), P[0][0], P[-1][0], P[-1][0]/TURN))
    H = np.loadtxt(os.path.join(arm, 'he4.hydro.hst'), comments='#')
    for tt in (0.25, 0.5, 1.0, 1.5):
        k = int(np.argmin([abs(x[0]-tt*TURN) for x in P]))
        t, r, q = P[k]
        if abs(t - tt*TURN) > 0.6*TURN:
            continue
        rho, v1, T = q[0], q[1], q[5]
        w = weight(rho, lo, hi)
        j = int(np.argmin(np.abs(H[:, 0] - t)))
        print('\n-- t = %.1f s = %.3f turn | hst t=%.1f dt=%.4g mass=%.8e '
              '1-mom=%.5e 1-KE=%.5e'
              % (t, t/TURN, H[j, 0], H[j, 1], H[j, 2], H[j, 3], H[j, 7]))
        print('     i   r/R      rho         T          v1          w')
        for i in range(96, 103):
            if i < len(r):
                print('   %3d %7.4f %11.3e %10.3e %11.3e %8.4f'
                      % (i, r[i]/RS, rho[i], T[i], v1[i], w[i]))
        m = np.abs(v1[96:])
        print('   max|v1| i>=96 : %.3e at i=%d ; i>=98 : %.3e at i=%d'
              % (m.max(), 96+int(np.argmax(m)), np.abs(v1[98:]).max(),
                 98+int(np.argmax(np.abs(v1[98:])))))


if __name__ == '__main__':
    main()
