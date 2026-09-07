"""Merge a high-temperature atomic Rosseland table (OPLIB/OPAL type 1) with a
low-temperature molecular one (AESOPUS 2.0) into ONE table of log10 kappa_R on a
regular (log10 T, log10 rho) grid, in the ASCII layout AthenaK's radiative conduction
reads through rad_kappa_src = table_rho (see src/pgen/red_giant.cpp).

Both sources are tabulated against logR = log10 rho - 3 log10 T + 18.  Below --tlo the
low-T table is used, above --thi the atomic one, and between them log10 kappa is
blended linearly in log10 T.  Nodes outside a source's logR range are clamped to its
edge (both edges are in the optically thin / thick regimes a giant envelope never
visits; the join region is checked and printed).

usage: merge_rosseland.py OPLIB_FILE AESOPUS_FILE OUT [--tlo 4.0 --thi 4.2]
"""
import re
import argparse
import numpy as np


def read_oplib(fn):
    lines = open(fn).read().splitlines()
    i0 = [i for i, ln in enumerate(lines) if ln.startswith('logT')][0]
    lR = np.array(lines[i0+1].split(), float)
    rows = [ln.split() for ln in lines[i0+2:] if ln.strip()]
    lT = np.array([r[0] for r in rows], float)
    K = np.array([r[1:] for r in rows], float)
    return lT, lR, K


def read_aesopus(fn):
    txt = open(fn).read()
    m = re.search(r'log10\(R\) range: nre=\s*(\d+)\s+values from\s*([-\d.]+)'
                  r'\s+to\s*([-\d.]+)', txt)
    nR, r0, r1 = int(m.group(1)), float(m.group(2)), float(m.group(3))
    lR = np.linspace(r0, r1, nR)
    rows = [ln.split() for ln in txt.splitlines()
            if ln.strip() and not ln.startswith('#')]
    lT = np.array([r[0] for r in rows], float)
    K = np.array([r[1:] for r in rows], float)
    o = np.argsort(lT)
    return lT[o], lR, K[o]


def interp2(lT, lR, K, T, R):
    """bilinear in (logT, logR), clamped to the grid"""
    T = np.clip(T, lT[0], lT[-1])
    R = np.clip(R, lR[0], lR[-1])
    i = np.clip(np.searchsorted(lT, T) - 1, 0, len(lT)-2)
    j = np.clip(np.searchsorted(lR, R) - 1, 0, len(lR)-2)
    fx = (T - lT[i])/(lT[i+1]-lT[i])
    fy = (R - lR[j])/(lR[j+1]-lR[j])
    return ((1-fx)*(1-fy)*K[i, j] + fx*(1-fy)*K[i+1, j]
            + (1-fx)*fy*K[i, j+1] + fx*fy*K[i+1, j+1])


ap = argparse.ArgumentParser()
ap.add_argument('oplib')
ap.add_argument('aesopus')
ap.add_argument('out')
ap.add_argument('--tlo', type=float, default=4.0)
ap.add_argument('--thi', type=float, default=4.2)
ap.add_argument('--lt', type=float, nargs=3, default=[2.6, 8.0, 0.025],
                help='log10 T min max step')
ap.add_argument('--ld', type=float, nargs=3, default=[-14.0, 0.0, 0.05],
                help='log10 rho min max step')
a = ap.parse_args()
Th, Rh, Kh = read_oplib(a.oplib)
Tl, Rl, Kl = read_aesopus(a.aesopus)
lT = np.arange(a.lt[0], a.lt[1] + 1e-9, a.lt[2])
lD = np.arange(a.ld[0], a.ld[1] + 1e-9, a.ld[2])
TT, DD = np.meshgrid(lT, lD, indexing='ij')
RR = DD - 3.0*TT + 18.0
kh = interp2(Th, Rh, Kh, TT, RR)
kl = interp2(Tl, Rl, Kl, TT, RR)
w = np.clip((TT - a.tlo)/(a.thi - a.tlo), 0.0, 1.0)     # 0 = low-T table, 1 = atomic
K = (1.0 - w)*kl + w*kh
# report the join
sel = ((TT >= a.tlo) & (TT <= a.thi)
       & (RR >= max(Rh[0], Rl[0])) & (RR <= min(Rh[-1], Rl[-1])))
d = np.abs(kh - kl)[sel]
print('join %.2f-%.2f: |log kappa(atomic) - log kappa(low-T)| median %.3f, '
      '90%% %.3f, max %.3f over %d nodes'
      % (a.tlo, a.thi, np.median(d), np.percentile(d, 90), d.max(), sel.sum()))
with open(a.out, 'w') as f:
    f.write('# AthenaK stellar Rosseland opacity table, log10 kappa_R [cm^2/g] on '
            '(log10 T[K], log10 rho[g/cm3])\n')
    f.write('# sources: %s (atomic, T >= %.2f) + %s (molecular+grains, T <= %.2f), '
            'linear blend in log T between\n'
            % (a.oplib.split('/')[-1], a.thi, a.aesopus.split('/')[-1], a.tlo))
    f.write('# grid: nT nD lTmin dlT lDmin dlD\n')
    f.write('# %d %d %g %g %g %g\n' % (len(lT), len(lD), lT[0], a.lt[2], lD[0], a.ld[2]))
    f.write('# then nT*nD rows, T slowest: log10 kappa_R\n')
    for i in range(len(lT)):
        for j in range(len(lD)):
            f.write('%.5f\n' % K[i, j])
print('wrote', a.out, K.shape)
