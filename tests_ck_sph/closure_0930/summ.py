"""VARIANTS.md table rows.  usage: python3 summ.py SUF [SUF ...]  (reference = f0)"""
import sys
import os
import numpy as np
R = '/viper/ptmp2/jinma/cksph_test_0930/runs'
NB, NG = 11, 8
CASES = ['tr_x1.2', 'tr_x1.7', 'tr_x3', 'st_x1.2', 'st_x1.7',
         'w_T2000', 'w_T3000', 'o_T2000', 'o_T3000']


def exact(case, suf):
    f = f'{R}/{case}_{suf}/exact.txt'
    if not os.path.exists(f):
        return None
    L = open(f).read().splitlines()
    band = {}
    tot = np.nan
    dep = []
    for ln in L:
        w = ln.split()
        if ln.startswith('band'):
            band[int(w[1])] = float(w[3])
        elif ln.startswith('TOTAL'):
            tot = float(w[2])
        elif len(w) == 7 and w[0].isdigit():
            dep.append([float(x) for x in w])
    return tot, band, np.array(dep)


def dump(case, suf):
    d = np.loadtxt(f'{R}/{case}_{suf}/dump.txt')
    return d


def grad(suf, ref='f0'):
    d, e = dump('grad', suf), dump('grad', ref)
    rf = d[:, 1]
    nc = len(rf) - 1
    rho = d[:nc, 13]
    kap = d[:nc, 41:41+NB*NG].reshape(nc, NB, NG)
    gw = d[0, 41+NB*NG:41+NB*NG+NG]
    Bb = d[:nc, 25:36]
    kr = kap*rho[:, None, None]
    dtc = kr*np.diff(rf)[:, None, None]
    ttop = np.cumsum(dtc[::-1], axis=0)[::-1]        # tau from face i to the top
    tmin = ttop.min(axis=(1, 2))
    # deep faces: min-g tau from the face to the top > 50 (every g thick)
    idx = [i for i in range(1, nc) if tmin[i] > 50.0]
    F, F0 = d[:nc, 4], e[:nc, 4]
    # non-grey diffusion flux at face i: -(4 pi/3) sum w_g dB/dtau over the two half cells
    fd = np.zeros(nc)
    for i in idx:
        dt = 0.5*(dtc[i-1] + dtc[i])
        fd[i] = -(4*np.pi/3)*np.sum(gw[None, :]*(Bb[i] - Bb[i-1])[:, None]/dt)
    idx = np.array(idx)
    r = F[idx]/F0[idx] - 1
    rd = F0[idx]/fd[idx] - 1
    rv = F[idx]/fd[idx] - 1
    return len(idx), np.max(np.abs(r)), np.mean(r), np.max(np.abs(rd)), np.max(np.abs(rv))


def newton(suf):
    f = f'{R}/nw10_{suf}/rt.txt'
    S = [ln.split() for ln in open(f) if ln.startswith('S ')]
    p = [int(s[7]) for s in S]
    nc = [int(s[8]) for s in S]
    bud = [float(s[16]) for s in S]
    return len(S), max(p), sum(nc), max(bud), bud[-1]


if __name__ != "__main__":
    pass
for suf in (sys.argv[1:] if __name__ == "__main__" else []):
    print(f'## {suf}')
    for c in CASES:
        x = exact(c, suf)
        if x is None:
            print(c, 'missing')
            continue
        tot, band, dep = x
        b9 = {k: v for k, v in band.items() if k < 10}
        wb = max(b9, key=lambda k: abs(b9[k] - 1))
        s = f'{c:8s} L {tot:.4f}  worst b{wb} {b9[wb]:.4f}  b10 {band[10]:.4f}'
        if c[0] in 'wo' and len(dep):
            thin = dep[dep[:, 3] < 40]
            q = thin[:, 6]
            s += (f'  Q {q[0]:.3f}->{q[-1]:.3f} [min {q[:-3].min():.3f} max '
                  f'{q[:-3].max():.3f} excl top3]  r/rin '
                  f'{thin[0, 1]:.2f}-{thin[-1, 1]:.2f}')
        print(s)
    try:
        n, mx, mn, rd, rv = grad(suf)
        print(f'grad deep faces {n}: |F/F_edd-1| max {mx:.2e} mean {mn:+.2e};'
              f' F_edd/F_diff-1 max {rd:.2e}; F/F_diff-1 max {rv:.2e}')
    except Exception as ex:
        print('grad', ex)
    try:
        n, p, nc, bm, bl = newton(suf)
        print(f'nw10 cycles {n} passes max {p} nonconv {nc} budmax max {bm:.2e} '
              f'last {bl:.2e}')
    except Exception as ex:
        print('nw10', ex)
