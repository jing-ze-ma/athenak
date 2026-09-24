"""T-S6: f_K(r) of vet_col (column dumps of the first build) vs sphref.py.
usage: ts6.py <run-root>   (runs t6_<case>_<tag>/dump.00000.txt)"""
import sys
import numpy as np
from sphref import moments, eddington_E

root = sys.argv[1]
c, fin, rin, rout, q = 100.0, 1.0, 1.0, 5.0, 0.5
KT = {'FS': 1.0e-12, 'EX': 10.0, 'TK': 1.0e4}


def load(tag):
    d = np.loadtxt('%s/t6_%s/dump.00000.txt' % (root, tag))
    return d[:, 0], d[:, 6], d[:, 8], d[:, 9]


def reference(case, rr, Eb, Fb):
    kt = KT[case]
    C = kt*rin*rin
    if case == 'FS':
        C = 0.0
    Sf = lambda r: eddington_E(r, 1.0, kt, fin, rin, rout, c, q)
    out = []
    for r in rr:
        J, H, K = moments(r, rin, rout, C, Sf, Eb, Fb, c)
        out.append(K/J)
    return np.array(out)


for case in ['FS', 'EX', 'TK']:
    # the fine reference grid, interpolated to every run's radii
    r0, _, E0, F0 = load(case + '_r256')
    if case == 'FS':
        Eb, Fb = E0[0], F0[0]   # the code's own bottom-cell values (angular test only)
    else:
        Eb = eddington_E(rin, 1.0, KT[case], fin, rin, rout, c, q)
        Fb = fin
    rg = np.concatenate([np.linspace(1.0, 4.9, 157), np.linspace(4.905, 5.0, 20)])
    fg = reference(case, rg, Eb, Fb)
    print('== %s (reference f_K: r=1.0 %.6f, 2.0 %.6f, 3.0 %.6f, 4.0 %.6f, 5.0 %.6f)' % (
        case, fg[0], np.interp(2.0, rg, fg), np.interp(3.0, rg, fg),
        np.interp(4.0, rg, fg), fg[-1]))
    for tag in ['r32', 'r64', 'r128', 'r256', 'c2', 'c4', 'c16', 'c32', 'p2', 'p4']:
        try:
            r, fk, E, F = load(case + '_' + tag)
        except OSError:
            continue
        if case == 'FS':
            fr = reference(case, r, E[0], F[0])
        else:
            fr = np.interp(r, rg, fg)
        e = fk - fr
        print('%-5s n=%3d  L1 %.3e  Linf %.3e (r=%.3f)  top fK %.5f ref %.5f  '
              'max|fK-1/3| %.2e' % (tag, len(r), np.mean(np.abs(e)), np.max(np.abs(e)),
                                    r[np.argmax(np.abs(e))], fk[-1], fr[-1],
                                    np.max(np.abs(fk - 1/3.))))
