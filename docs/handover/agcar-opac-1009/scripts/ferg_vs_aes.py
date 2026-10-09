#!/usr/bin/env python3
"""(a) Ferguson 2005 GS98 vs AESOPUS 2.1 gas Rosseland (both X->0.36, Z 0.02), where both are DATA,
per log T band, all log R in [-8,1] and the AG Car-like log R <= -3; (b) old / ext / ext2 kR, kP and
kP/kR along the A and B columns (ICs and evolved geos/col runs)."""
import sys
import numpy as np
sys.argv += [] ; sys.path.insert(0, '/viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/scripts')
_argv = sys.argv; sys.argv = ['x']
import build_ext as b  # noqa: E402
import column_compare as cc  # noqa: E402
from convert_tops import bilin  # noqa: E402
from compare import read_repo  # noqa: E402
sys.argv = _argv
tF, rF, KF, _ = b.aes_field('fergR')
tA, rA, KA, _ = b.aes_field('xint')
print('== (a) Rosseland Ferguson - AESOPUS [dex], at Ferguson nodes, both data')
for a, c in ((3.4, 3.6), (3.6, 3.7), (3.7, 3.8), (3.8, 3.9), (3.9, 4.0), (4.0, 4.1), (4.1, 4.2), (4.2, 4.5)):
    m = (tF >= a) & (tF < c)
    TT, RR = np.meshgrid(tF[m], rF, indexing='ij')
    d = KF[m] - bilin(tA, rA, KA, TT, RR)
    lo = RR <= -3
    print('  log T %.1f-%.1f: all log R: median %+.3f |d| 90%% %.3f max %.3f | log R -8..-3: median %+.3f '
          '|d| 90%% %.3f max %.3f' % (a, c, np.median(d), np.percentile(abs(d), 90), abs(d).max(),
                                       np.median(d[lo]), np.percentile(abs(d[lo]), 90), abs(d[lo]).max()))
B = '/viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/'
T3 = dict(old=cc.T_['old'], ext=cc.T_['new'],
          ext2=[cc.code_extend(read_repo(B + f)) for f in ('rosseland_ext2_gs98_x0.36_z0.02.txt',
                                                         'planck_ext2_gs98_x0.36_z0.02.txt')])
print('== (b) along the columns: median (min-max) per band; log R max where log T < 4.2')
for X, rph in (('A', 388.3007), ('B', 101.2962)):
    ic = np.loadtxt(cc.B + '/geos/ic%s/ic_agcar_%s_ge.txt' % (X, X))
    cols = {'IC %s' % X: (ic[:, 0]/cc.RS/rph, ic[:, 1], ic[:, 5])}
    x, rho, T, t = cc.evolved(X, rph)
    cols['col%s t=%.2e' % (X, t)] = (x, rho, T)
    for nm, (x, rho, T) in cols.items():
        cool = np.log10(T) < 4.2
        lRmax = (np.log10(rho) - 3*np.log10(T) + 18)[cool].max() if cool.any() else np.nan
        print('-- %s (max log R at log T<4.2: %.2f)' % (nm, lRmax))
        for a, c in ((0.9, 1.0), (1.0, 1.1), (1.1, 1.5), (1.5, 3.2)):
            m = (x >= a) & (x < c)
            s = []
            for k in ('old', 'ext', 'ext2'):
                kr = cc.look(T3[k][0], T[m], rho[m]); kp = cc.look(T3[k][1], T[m], rho[m])
                q = kp/kr
                s.append('%s kR %.3g kP %.3g kP/kR %.3g (%.3g-%.3g)' % (k, np.median(kr), np.median(kp),
                                                                     np.median(q), q.min(), q.max()))
            kr1 = cc.look(T3['ext'][0], T[m], rho[m]); kr2 = cc.look(T3['ext2'][0], T[m], rho[m])
            print('  %.1f-%.1f R_ph T %.0f-%.0f K: %s | kR ext2/ext dlog %+.2f..%+.2f' % (
                a, c, T[m].min(), T[m].max(), ' | '.join(s), np.log10(kr2/kr1).min(),
                np.log10(kr2/kr1).max()))
