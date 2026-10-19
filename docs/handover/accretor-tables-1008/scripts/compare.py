#!/usr/bin/env python3
"""(1) TOPS Rosseland vs the repo's OPLIB-based Rosseland table;
(2) kappa_P/kappa_R from TOPS at the bumps; (3) Ferguson 2005 low-T Planck
vs Ferguson Rosseland."""
import numpy as np


def read_repo(fn):
    nT = nD = None
    vals = []
    for ln in open(fn):
        if ln.startswith('#'):
            if nT is None:
                p = ln[1:].split()
                try:
                    nT, nD = int(p[0]), int(p[1])
                    lt0, dlt, ld0, dld = (float(x) for x in p[2:6])
                except (ValueError, IndexError):
                    nT = None
            continue
        if ln.strip():
            vals.append(float(ln))
    K = np.array(vals).reshape(nT, nD)
    return lt0 + dlt*np.arange(nT), ld0 + dld*np.arange(nD), K


def read_ferg(fn):
    ln = open(fn).read().splitlines()
    i0 = [i for i, s in enumerate(ln) if s.lstrip().startswith('log T')][0]
    lR = np.array(ln[i0].split()[2:], float)
    rows = [s.split() for s in ln[i0+1:] if s.strip()]
    lT = np.array([r[0] for r in rows], float)
    K = np.array([r[1:] for r in rows], float)
    o = np.argsort(lT)
    return lT[o], lR, K[o]


def at(lT, lD, K, t, d):
    i = np.argmin(abs(lT-t))
    j = np.argmin(abs(lD-d))
    return K[i, j]


def main():
    print('=== (1) TOPS Rosseland vs repo OPLIB Rosseland, GS98 X=0.7 Z=0.014 ===')
    rT, rD, rK = read_repo('/viper/u2/jinma/ATHENAK/athenak/data/stellar_opac/'
                           'rosseland_gs98_x0.7_z0.014.txt')
    tT, tD, tK = read_repo('rosseland_tops_gs98_x0.7_z0.014.txt')
    # overlap where BOTH are data: logT 3.80..7.05 and logR = lrho-3lT+18 in -8..1
    TT, DD = np.meshgrid(rT, rD, indexing='ij')
    RR = DD - 3*TT + 18
    m = (TT >= 3.80) & (TT <= 7.05) & (RR >= -8) & (RR <= 1)
    d = tK[m] - rK[m]
    print('overlap nodes %d  median |dlog| %.3f  90%% %.3f  max %.3f  '
          'median signed %+.3f (TOPS - repo)'
          % (m.sum(), np.median(abs(d)), np.percentile(abs(d), 90),
             abs(d).max(), np.median(d)))
    mb = m & (TT >= 5.20) & (TT <= 5.40) & (DD >= -10) & (DD <= -6)
    print('Fe bump (logT 5.2-5.4, log rho -10..-6): median |dlog| %.3f max %.3f'
          % (np.median(abs(tK[mb]-rK[mb])), abs(tK[mb]-rK[mb]).max()))
    mh = m & (TT >= 4.60) & (TT <= 4.80) & (DD >= -11) & (DD <= -8)
    print('He II bump (logT 4.6-4.8, log rho -11..-8): median |dlog| %.3f max %.3f'
          % (np.median(abs(tK[mh]-rK[mh])), abs(tK[mh]-rK[mh]).max()))

    print()
    print('=== (2) TOPS kappa_P / kappa_R (dex, and ratio) ===')
    pT, pD, pK = read_repo('planck_gs98_x0.7_z0.014.txt')
    p2T, p2D, p2K = read_repo('planck_he_x0.0_z0.02.txt')
    r2T, r2D, r2K = read_repo('rosseland_tops_he_x0.0_z0.02.txt')
    hdr = '%-34s %10s %10s' % ('regime (logT, log rho)', 'A: X=.7', 'B: He')
    print(hdr)
    for lab, t, dd in [('photosphere-ish 4.0, -9', 4.0, -9),
                       ('photosphere-ish 4.0, -11', 4.0, -11),
                       ('He II bump 4.7, -10', 4.7, -10),
                       ('He II bump 4.7, -8', 4.7, -8),
                       ('Fe bump 5.3, -9', 5.3, -9),
                       ('Fe bump 5.3, -7', 5.3, -7),
                       ('deep 6.0, -6', 6.0, -6),
                       ('deep 6.5, -5', 6.5, -5)]:
        a = at(pT, pD, pK, t, dd) - at(tT, tD, tK, t, dd)
        b = at(p2T, p2D, p2K, t, dd) - at(r2T, r2D, r2K, t, dd)
        print('%-34s %10s %10s' % (lab, '%.2f dex (x%.0f)' % (a, 10**a),
                                   '%.2f dex (x%.0f)' % (b, 10**b)))

    print()
    print('=== (3) Ferguson 2005 GS98 low-T: kappa_P/kappa_R ===')
    for tag, pf, rf in [('X=0.7 Z=0.02', 'ferguson05/g98.pl.7.02.tpon',
                         'ferguson05/ross/g98.7.02.tron'),
                        ('X=0.0 Z=0.02', 'ferguson05/g98.pl.0.02.tpon',
                         'ferguson05/ross/g98.0.02.tron')]:
        lT, lR, KP = read_ferg(pf)
        lT2, lR2, KR = read_ferg(rf)
        assert np.allclose(lR, lR2)
        print(' %s   Planck log T %.2f..%.2f (%d rows); Rosseland log T %.2f..%.2f '
              '(%d rows); log R %.1f..%.1f'
              % (tag, lT[0], lT[-1], len(lT), lT2[0], lT2[-1], len(lT2), lR[0],
                 lR[-1]))
        print('   %-8s' % 'logT' + ''.join('%12s' % ('logR=%g' % r)
                                           for r in (-7, -5, -3, -1)))
        for t in (3.5, 3.8, 4.0, 4.3, 4.5):
            row = ''
            for r in (-7, -5, -3, -1):
                i, i2 = np.argmin(abs(lT-t)), np.argmin(abs(lT2-t))
                j = np.argmin(abs(lR-r))
                row += '%12s' % ('x%.3g' % 10**(KP[i, j]-KR[i2, j]))
            print('   %-8.2f' % t + row)


if __name__ == '__main__':
    main()
