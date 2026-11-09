#!/usr/bin/env python3
"""acc.py DIR REF ARM [ARM ...]: accuracy of each arm's final m1 state vs the reference arm REF (all in DIR, newest
<basename>.m1.*.bin of each), by region of the reference reduced flux f = |F|/(c E):
  interior f < 0.03 (tau >~ 10), transition 0.03 <= f < 0.3, thin f >= 0.3.
Per region: E rel error (median, p99, max), F_1 error |dF1|/(c E_ref) (median, p99, max); and the radial luminosity
L_i = sum_{j,k} w F1 (w = solid-angle weight of the cell; cs equiangular or sp), error |dL_i|/|L_ref,i| by region of
the shell-mean f (median, max)."""
import sys
import os
import glob
import numpy as np
sys.path.insert(0, os.environ.get('ATHENAK_VIS', os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../vis/python')))
import bin_convert as bc  # noqa: E402

C = 2.99792458e10


def load(d):
    f = sorted(glob.glob(os.path.join(d, 'bin', '*.m1.*.bin')))[-1]
    r = bc.read_binary(f)
    G = np.asarray(r['mb_geometry'])
    D = r['mb_data']
    return f, G, D


def weights(G, shp):
    nk, nj, ni = shp
    W = []
    cs = np.all(np.abs(G[:, 2:6]) <= 1.0001)
    for m in range(len(G)):
        x2 = G[m, 2] + (G[m, 3] - G[m, 2]) * (np.arange(nj) + 0.5) / nj
        x3 = G[m, 4] + (G[m, 5] - G[m, 4]) * (np.arange(nk) + 0.5) / nk
        d2 = (G[m, 3] - G[m, 2]) / nj
        d3 = (G[m, 5] - G[m, 4]) / nk
        if cs:
            X = np.tan(np.pi / 4 * x2)[None, :]
            Y = np.tan(np.pi / 4 * x3)[:, None]
            w = (1 + X * X) * (1 + Y * Y) / (1 + X * X + Y * Y) ** 1.5 * d2 * d3
        else:
            w = np.sin(x2)[None, :] * d2 * d3 * np.ones((nk, 1))
        W.append(np.broadcast_to(w[:, :, None], shp))
    return np.array(W)


def arr(D, v):
    return np.array([np.asarray(x, float) for x in D[v]])


def stats(x):
    if x.size == 0:
        return 'n=0'
    return 'med %.2e p99 %.2e max %.2e' % (np.median(x), np.percentile(x, 99), x.max())


def main():
    d, ref = sys.argv[1], sys.argv[2]
    fr, G, Dr = load(os.path.join(d, ref))
    Er = arr(Dr, 'm1_e')
    F1r, F2r, F3r = arr(Dr, 'm1_f1'), arr(Dr, 'm1_f2'), arr(Dr, 'm1_f3')
    f = np.sqrt(F1r**2 + F2r**2 + F3r**2) / (C * np.maximum(Er, 1e-300))
    if np.percentile(f, 50) > 1.5:   # flux stored as F/c
        f = f / C
        scl = 1.0
    else:
        scl = C
    W = weights(G, Er.shape[1:])
    # shells: blocks grouped by their x1 range, then the x1 index
    keys = sorted(set((round(g[0], 6), round(g[1], 6)) for g in np.asarray(G)[:, :2].tolist()))
    grp = [keys.index((round(g[0], 6), round(g[1], 6))) for g in np.asarray(G)[:, :2].tolist()]

    def shell(x):
        return np.concatenate([x[[m for m in range(len(grp)) if grp[m] == q]].sum(axis=(0, 1, 2))
                               for q in range(len(keys))])
    Lr = shell(W * F1r)
    fsh = shell(W * f) / shell(W)
    reg = {'interior': f < 0.03, 'transition': (f >= 0.03) & (f < 0.3), 'thin': f >= 0.3}
    regs = {'interior': fsh < 0.03, 'transition': (fsh >= 0.03) & (fsh < 0.3), 'thin': fsh >= 0.3}
    print('ref %s: cells %d interior %d transition %d thin %d; shells %d/%d/%d' % (
        ref, f.size, reg['interior'].sum(), reg['transition'].sum(), reg['thin'].sum(),
        regs['interior'].sum(), regs['transition'].sum(), regs['thin'].sum()))
    for a in sys.argv[3:]:
        try:
            fa, Ga, Da = load(os.path.join(d, a))
        except IndexError:
            print('%-5s MISSING' % a)
            continue
        Ea, F1a = arr(Da, 'm1_e'), arr(Da, 'm1_f1')
        eE = np.abs(Ea - Er) / np.maximum(np.abs(Er), 1e-300)
        eF = np.abs(F1a - F1r) / (scl * np.maximum(Er, 1e-300))
        La = shell(W * F1a)
        eL = np.abs(La - Lr) / np.maximum(np.abs(Lr), 1e-300)
        print('%-5s' % a)
        for k in ['interior', 'transition', 'thin']:
            m = reg[k]
            ms = regs[k]
            print('   %-10s E: %s | F1/(cE): %s | L: %s' % (
                k, stats(eE[m]), stats(eF[m]),
                ('med %.2e max %.2e' % (np.median(eL[ms]), eL[ms].max())) if ms.any() else 'n=0'))


if __name__ == '__main__':
    main()
