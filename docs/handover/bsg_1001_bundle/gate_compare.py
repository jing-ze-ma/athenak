#!/usr/bin/env python3
"""BSG arm-2 port gate: compare two runs of bsg_col_arm2.athinput (e.g. CPU vs GPU).

usage: python3 gate_compare.py <run dir A> <run dir B> [path to athenak/vis/python]
Prints, per run: final dt, mean Picard passes, L_top/L_in (last row and range), FATAL/nan
counts in the logs; then the max relative difference A vs B of every hst column and of
every variable in the last hydro_w and m1 dumps (all cells, and only cells with
rho > 1e-17, i.e. without the floor gas in 54-61 Rsun, whose velocities are noise).
Transverse v / F of the 4 x 4 column are round-off-seeded noise: gate on the others.
"""
import glob
import os
import sys

import numpy as np


def hst(d):
    f = os.path.join(d, 'bsg.user.hst')
    names = []
    with open(f) as fh:
        for line in fh:
            if line.startswith('#') and '[1]=' in line:
                names = [t.split('=')[1] for t in line[1:].split() if '=' in t]
    return names, np.loadtxt(f)


def logcounts(d):
    n_fatal = n_nan = 0
    for f in glob.glob(os.path.join(d, '*log*')) + glob.glob(os.path.join(d, '*.out')):
        if os.path.isdir(f):
            continue
        txt = open(f, errors='replace').read()
        n_fatal += txt.count('FATAL')
        n_nan += txt.lower().count(' nan') + txt.count('=nan')
    return n_fatal, n_nan


def summary(d):
    names, a = hst(d)
    c = {n: i for i, n in enumerate(names)}
    ratio = a[:, c['L_top']] / a[:, c['L_in']]
    pic = a[1:, c['Picard']]
    fat, nan = logcounts(d)
    print(f'{d}: rows {len(a)}  t_end {a[-1, 0]:.6e}  dt_end {a[-1, c["dt"]]:.6e}  '
          f'Picard mean {pic.mean():.3f}  L_top/L_in end {ratio[-1]:.6f} '
          f'range [{ratio[1:].min():.6f}, {ratio[1:].max():.6f}]  FATAL {fat}  nan {nan}')
    return names, a


def main():
    da, db = sys.argv[1], sys.argv[2]
    if len(sys.argv) > 3:
        sys.path.insert(0, sys.argv[3])
    names, a = summary(da)
    _, b = summary(db)
    n = min(len(a), len(b))
    print('hst max |A-B|/max(|A|,tiny) per column:')
    for i, nm in enumerate(names):
        x, y = a[:n, i], b[:n, i]
        s = np.maximum(np.abs(x), 1e-300)
        print(f'  {nm:12s} {np.max(np.abs(x - y) / s):.3e}')
    try:
        import bin_convert as bc
    except ImportError:
        print('bin_convert not found: pass athenak/vis/python as the 3rd argument')
        return
    mask = None
    for var in ('hydro_w', 'm1'):
        fa = sorted(glob.glob(os.path.join(da, 'bin', f'*.{var}.*.bin')))
        fb = sorted(glob.glob(os.path.join(db, 'bin', f'*.{var}.*.bin')))
        if not fa or not fb:
            continue
        A, B = bc.read_binary(fa[-1]), bc.read_binary(fb[-1])
        if var == 'hydro_w':
            mask = np.concatenate([np.ravel(m) for m in A['mb_data']['dens']]) > 1e-17
        print(f'last {var} dump {os.path.basename(fa[-1])}: per variable, all cells: '
              'max|A-B|/max|A|, max|A-B|, max|A|;  rho > 1e-17 only: max|A-B|/max|A|')
        for v in A['var_names']:
            x = np.concatenate([np.ravel(m) for m in A['mb_data'][v]])
            y = np.concatenate([np.ravel(m) for m in B['mb_data'][v]])
            s = np.abs(x).max() + 1e-300
            d = np.max(np.abs(x - y))
            dm = np.max(np.abs(x - y)[mask]) / (np.abs(x[mask]).max() + 1e-300)
            print(f'  {v:8s} {d / s:.3e}  {d:.3e}  {s:.3e}   {dm:.3e}')


if __name__ == '__main__':
    main()
