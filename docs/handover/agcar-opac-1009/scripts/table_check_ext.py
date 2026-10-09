#!/usr/bin/env python3
"""Checks of the extended tables against the old ones: finiteness, bitwise region, es limit,
Fe bump, jumps between neighbouring nodes (continuity) in the new parts."""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/data/stellar_opac/planck_tools')
from compare import read_repo
D = sys.argv[1]
TAG = sys.argv[2] if len(sys.argv) > 2 else 'ext'
O = '/viper/ptmp2/jinma/lbv_1008/agcar/tables/'
for nm, on in (('rosseland', 'rosseland_tops'), ('planck', 'planck_tops')):
    T, Dn, K = read_repo('%s/%s_%s_gs98_x0.36_z0.02.txt' % (D, nm, TAG))
    To, Do, Ko = read_repo(O + '%s_gs98_x0.36_z0.02.txt' % on)
    print('==', nm, K.shape, 'finite', np.isfinite(K).all(), 'range %.2f..%.2f' % (K.min(), K.max()))
    j0 = 140
    i42 = int(round((4.2 - 2.6)/0.025))
    same = K[i42:, j0:] == Ko[i42:, :]
    print('  bitwise log T>=4.2, log rho>=-14: %d/%d equal' % (same.sum(), same.size))
    d = K[:, j0:] - Ko
    for a, b in ((2.6, 3.764), (3.764, 4.0), (4.0, 4.2)):
        m = (T >= a - 1e-9) & (T < b - 1e-9)
        dd = d[m]
        print('  new-old (rho>=-14) log T %.3f-%.3f: median %+.3f, |d| 90%% %.3f, max %.3f' % (
            a, b, np.median(dd), np.percentile(abs(dd), 90), abs(dd).max()))
    # continuity: largest node-to-node jumps in the new rows/cols vs the old table's own
    dT = abs(np.diff(K, axis=0)); dD = abs(np.diff(K, axis=1))
    dTo = abs(np.diff(Ko, axis=0)); dDo = abs(np.diff(Ko, axis=1))
    sel = (T[:-1] >= 3.4) & (T[:-1] < 4.5)
    print('  max |dlogk| per T step (log T 3.4-4.5, all rho): new %.3f (old, rho>=-14: %.3f); '
          'per rho step: new %.3f (old %.3f)' % (dT[sel].max(), dTo[sel].max(),
                                                 dD[(T >= 3.4) & (T < 4.5)].max(),
                                                 dDo[(T >= 3.4) & (T < 4.5)].max()))
    for lt in (3.95, 4.0, 4.1, 4.2, 4.25):
        i = int(round((lt - 2.6)/0.025))
        print('  across log T %.3f->%.3f: max |dlogk| %.3f (rho<=-8)' % (
            T[i], T[i+1], dT[i, Dn <= -8].max()))
    if nm == 'rosseland':
        for lt, ld in ((7.0, -10), (7.0, -20), (6.0, -20), (5.0, -20), (4.5, -20), (4.0, -20),
                       (3.8, -20), (3.6, -20), (3.5, -16), (3.477, -16), (3.7, -16), (3.85, -16)):
            i = int(round((lt - 2.6)/0.025)); j = int(round((ld + 21)/0.05))
            print('  kR(log T %.3f, log rho %d) = %.4g' % (T[i], ld, 10**K[i, j]))
        for ld in (-10, -9, -8, -7.2, -6):
            j = int(round((ld + 21)/0.05)); m = (T > 4.9) & (T < 5.6)
            ii = np.argmax(np.where(m, K[:, j], -99))
            print('  Fe bump log rho %.1f: max kR %.3f at log T %.3f' % (ld, 10**K[ii, j], T[ii]))
