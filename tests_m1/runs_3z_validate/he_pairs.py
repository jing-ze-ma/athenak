#!/usr/bin/env python3
"""runs_3z item 4: cfl 0.3 vs 0.15 per scheme (and h2 vs be at cfl 0.15) from a common
restart.  KE_h, KE_v ratios from the hydro hst and max |T(z)/T_ref - 1| (rt_profile row 5,
all layers / layers 0..n-21 / top 20) at common profile times."""
import numpy as np
import os
import struct


def prof(d):
    recs = []
    with open(os.path.join(d, 'rt_profile.bin'), 'rb') as fh:
        while True:
            h = fh.read(16)
            if len(h) < 16:
                break
            t, = struct.unpack('<d', h[:8])
            n1, nv = struct.unpack('<ii', h[8:])
            fh.read(8 * n1)
            b = fh.read(8 * n1 * nv)
            if len(b) < 8 * n1 * nv:
                break
            recs.append((t, None, np.frombuffer(b, '<f8').reshape(nv, n1).copy()))
    return recs


def at(t, y, tt):
    return float(np.interp(tt, t, y))


H = '/viper/ptmp2/jinma/validate_0923/he3d/'
TT = (21., 25., 40., 70., 100., 150., 200., 300., 400., 500.)


def load(d):
    a = np.loadtxt(H + d + '/m1slab.hydro.hst')
    _, i = np.unique(a[:, 0], return_index=True)
    a = a[i]
    return a[:, 0], a[:, 7], a[:, 8] + a[:, 9], prof(H + d)


for a, b in (('h2_c0.3', 'h2_c0.15'), ('be_c0.3', 'be_c0.15'), ('h2_c0.15', 'be_c0.15'),
             ('h2_c0.3', 'be_c0.3')):
    A, B = load(a), load(b)
    print('== %s / %s' % (a, b))
    for tt in TT:
        print('  t=%5.0f  KE_h ratio %.4f  KE_v ratio %.4f' % (
            tt, at(A[0], A[2], tt) / at(B[0], B[2], tt),
            at(A[0], A[1], tt) / at(B[0], B[1], tt)))
    for tt in (45., 95., 145., 245., 345., 495.):
        ra = min(A[3], key=lambda r: abs(r[0] - tt))
        rb = min(B[3], key=lambda r: abs(r[0] - tt))
        rel = ra[2][5] / rb[2][5] - 1.0
        n = rel.size
        print('  T(z) t=%6.1f/%6.1f max|dT/T| all %.3e  interior %.3e  top20 %.3e' % (
            ra[0], rb[0], abs(rel).max(), abs(rel[:n - 20]).max(),
            abs(rel[n - 20:]).max()))
