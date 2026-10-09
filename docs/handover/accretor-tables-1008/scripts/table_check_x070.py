#!/usr/bin/env python3
"""Checks of the X 0.70 Z 0.02 ext2 tables: finiteness, bitwise TOPS region, electron-scattering
limit, Fe bump vs the repo X 0.7 table (OPLIB Z 0.014) and the dev pair (TOPS X 0.7 Z 0.008),
continuity across the blend, kP/kR, and new/dev ratios in representative regions.
usage: table_check_x070.py <dir with *_ext2_gs98_x0.70_z0.02.txt> <tables dir with *_tops_*>"""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/data/stellar_opac/planck_tools')
from compare import read_repo   # noqa: E402

D, TB = sys.argv[1], sys.argv[2]
DEV = '/viper/ptmp2/jinma/bsg_1001/opac/'
REPO = '/viper/u2/jinma/ATHENAK/athenak/data/stellar_opac/rosseland_gs98_x0.7_z0.014.txt'
T, Dn, KR = read_repo(D + '/rosseland_ext2_gs98_x0.70_z0.02.txt')
_, _, KP = read_repo(D + '/planck_ext2_gs98_x0.70_z0.02.txt')
To, Do, KRo = read_repo(TB + '/rosseland_tops_gs98_x0.70_z0.02.txt')
_, _, KPo = read_repo(TB + '/planck_tops_gs98_x0.70_z0.02.txt')
Tv, Dv, KRv = read_repo(DEV + 'rosseland_tops_x0.7_z0.008.txt')
_, _, KPv = read_repo(DEV + 'planck_tops_x0.7_z0.008.txt')
Tr, Dr, KRr = read_repo(REPO)


def at(lT, lD, K, t, d):
    i = int(round((t - lT[0])/(lT[1] - lT[0])))
    j = int(round((d - lD[0])/(lD[1] - lD[0])))
    j = min(max(j, 0), len(lD) - 1)
    return K[i, j]


print('grid', KR.shape, 'log T %.3f..%.3f, log rho %.2f..%.2f' % (T[0], T[-1], Dn[0], Dn[-1]))
for nm, K in (('R', KR), ('P', KP)):
    print('finite %s %s, range %.2f..%.2f' % (nm, np.isfinite(K).all(), K.min(), K.max()))
j0 = int(round((-14 - Dn[0])/0.05))
i42 = int(round((4.2 - T[0])/0.025))
for nm, K, Ko in (('R', KR, KRo), ('P', KP, KPo)):
    s = K[i42:, j0:] == Ko[i42:, :]
    print('bitwise %s (log T>=4.2, log rho>=-14) vs TOPS-only: %d/%d' % (nm, s.sum(), s.size))

print('\n== electron-scattering limit, 0.2(1+X) = 0.34 (Thomson, full ionisation)')
for t, d in ((7.0, -10), (7.5, -10), (7.6, -15), (8.0, -10), (7.0, -18), (6.83, -15.3),
             (6.0, -20), (5.0, -20), (4.5, -20), (4.0, -20)):
    print('  kR(log T %.2f, log rho %d) = %.4f   dev Z0.008: %s   repo z0.014: %s' % (
        t, d, 10**at(T, Dn, KR, t, d),
        '%.4f' % 10**at(Tv, Dv, KRv, t, d) if d >= -14 else '(below its -14 edge)',
        '%.4f' % 10**at(Tr, Dr, KRr, t, d) if d >= Dr[0] else '-'))

print('\n== Fe bump (max kR at log T 4.9-5.6)')
for ld in (-12, -10, -9, -8, -7.2, -6):
    out = []
    for lab, lT, lD, K in (('new Z0.02', T, Dn, KR), ('dev TOPS Z0.008', Tv, Dv, KRv),
                           ('repo OPLIB Z0.014', Tr, Dr, KRr)):
        j = int(round((ld - lD[0])/(lD[1] - lD[0])))
        if j < 0:
            out.append('%s -' % lab)
            continue
        m = (lT > 4.9) & (lT < 5.6)
        ii = np.argmax(np.where(m, K[:, j], -99))
        out.append('%s %.3f @ %.3f' % (lab, 10**K[ii, j], lT[ii]))
    print('  log rho %5.1f: ' % ld + ' | '.join(out))

print('\n== continuity across the blend (max |dlog k| per 0.025 log T step, log rho <= -8 / all)')
dTR, dTP = abs(np.diff(KR, axis=0)), abs(np.diff(KP, axis=0))
for lt in (3.9, 3.95, 3.975, 4.0, 4.1, 4.175, 4.2, 4.225, 4.25):
    i = int(round((lt - T[0])/0.025))
    m = Dn <= -8
    print('  log T %.3f->%.3f: R %.3f / %.3f, P %.3f / %.3f' % (
        T[i], T[i+1], dTR[i, m].max(), dTR[i].max(), dTP[i, m].max(), dTP[i].max()))
sel = (T[:-1] >= 3.5) & (T[:-1] < 4.6)
mr = Dn <= -8
print('  context, log T 3.5-4.6 (H/He ionisation), rho<=-8: max R step %.3f, P step %.3f' % (
    dTR[sel][:, mr].max(), dTP[sel][:, mr].max()))
dDR, dDP = abs(np.diff(KR, axis=1)), abs(np.diff(KP, axis=1))
print('  max |dlog k| per 0.05 log rho step, log T 3.5-8: R %.3f, P %.3f' % (
    dDR[T >= 3.5].max(), dDP[T >= 3.5].max()))

print('\n== kP/kR')
for t, d in ((3.7, -10), (3.9, -10), (4.0, -10), (4.2, -10), (4.5, -10), (5.2, -10), (5.2, -8),
             (6.0, -6), (6.83, -15.3), (7.0, -12), (4.0, -16), (4.5, -18)):
    print('  log T %.2f log rho %6.1f: kP %.4g kR %.4g kP/kR %.3g' % (
        t, d, 10**at(T, Dn, KP, t, d), 10**at(T, Dn, KR, t, d),
        10**(at(T, Dn, KP, t, d) - at(T, Dn, KR, t, d))))
for a, b in ((3.5, 4.0), (4.0, 4.2), (4.2, 5.0), (5.0, 6.0), (6.0, 7.0), (7.0, 8.0)):
    m = (T >= a - 1e-9) & (T < b - 1e-9)
    r = (KP - KR)[m][:, Dn >= -18]
    print('  log T %.1f-%.1f, log rho -18..0: log kP/kR min %+.2f median %+.2f; frac kP<kR %.2f' % (
        a, b, r.min(), np.median(r), (r < 0).mean()))

print('\n== new (Z 0.02, ext2) - dev (TOPS X0.7 Z0.008) [dex], log rho -14..0')
for a, b in ((3.764, 4.0), (4.0, 4.2), (4.2, 5.0), (5.0, 5.5), (5.5, 6.5), (6.5, 7.065)):
    m = (T >= a - 1e-9) & (T < b - 1e-9)
    dd = KR[m][:, j0:] - KRv[m]
    pp = KP[m][:, j0:] - KPv[m]
    print('  log T %.3f-%.3f: R median %+.3f range %+.3f..%+.3f | P median %+.3f range %+.3f..%+.3f' % (
        a, b, np.median(dd), dd.min(), dd.max(), np.median(pp), pp.min(), pp.max()))
m = T > 7.065
print('  log T > 7.065 (dev edge-filled): new kR(8.0,-10) %.4f vs dev %.4f' % (
    10**at(T, Dn, KR, 8.0, -10), 10**at(Tv, Dv, KRv, 8.0, -10)))
print('\n== ext2 - TOPS-only (below log T 4.2, log rho >= -14) [dex]')
for a, b in ((3.5, 3.764), (3.764, 4.0), (4.0, 4.2)):
    m = (T >= a - 1e-9) & (T < b - 1e-9)
    for nm, K, Ko in (('R', KR, KRo), ('P', KP, KPo)):
        dd = (K[m][:, j0:] - Ko[m])[:, Do <= -8]
        print('  %s log T %.3f-%.3f, rho<=1e-8: median %+.3f, max |d| %.3f' % (
            nm, a, b, np.median(dd), abs(dd).max()))
