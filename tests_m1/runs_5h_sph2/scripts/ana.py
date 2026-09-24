#!/usr/bin/env python3
"""m1-sph2 analysis.
  ana.py order <glob> <col> <sub> <ref> <d1> <d2> ... [--pert X]
      L1 of column col of the LAST file matching <glob> (relative to each run dir) vs
      the reference run, divided by the L1 of (ref - X) when --pert X is given (else by
      L1 of ref); successive orders log2(e_k/e_{k+1}) (the dirs halve dt).
  ana.py steady <glob> <dirA> <dirB>
      max |a/b - 1| of every column of the LAST matching file (m1_e, m1_f1, ...)."""
import glob
import os
import sys

import numpy as np


def tab(fn):
    with open(fn) as f:
        h1 = f.readline()
        cols = f.readline().lstrip('#').split()
    d = np.loadtxt(fn, ndmin=2)
    out = {c: d[:, k] for k, c in enumerate(cols)}
    out['time'] = float(h1.split('time=')[1].split()[0])
    return out


def last(d, g):
    fs = sorted(glob.glob(os.path.join(d, g)))
    if not fs:
        raise SystemExit(f'no {g} in {d}')
    return fs[-1]


if __name__ == "__main__" and sys.argv[1] == "order":
    g, col, ref = sys.argv[2], sys.argv[3], sys.argv[4]
    rest = sys.argv[5:]
    pert = None
    if '--pert' in rest:
        k = rest.index('--pert')
        pert = float(rest[k+1])
        rest = rest[:k] + rest[k+2:]
    r = tab(last(ref, g))
    nrm = np.mean(np.abs(r[col] - pert)) if pert is not None else np.mean(np.abs(r[col]))
    es = []
    for d in rest:
        a = tab(last(d, g))
        if abs(a['time'] - r['time']) > 1e-9*max(1.0, abs(r['time'])):
            print(f'  WARNING time {a["time"]} vs ref {r["time"]}')
        es.append(np.mean(np.abs(a[col] - r[col]))/nrm)
    s = ' '.join(f'{e:.3e}' for e in es)
    o = ' '.join(f'{np.log2(es[k]/es[k+1]):.2f}' for k in range(len(es)-1)
                 if es[k+1] > 0)
    print(f'{col}: rel L1 err {s} | orders {o} (t={r["time"]:.6g})')
elif __name__ == '__main__' and sys.argv[1] == 'steady':
    g, a, b = sys.argv[2:5]
    ta, tb = tab(last(a, g)), tab(last(b, g))
    out = []
    for c in ta:
        if c in ('time', 'gid', 'i', 'x1v'):
            continue
        sc = np.max(np.abs(tb[c]))
        if sc == 0:
            continue
        out.append(f'{c} {np.max(np.abs(ta[c]-tb[c]))/sc:.2e}')
    print(f'max|A-B|/max|B| at t={ta["time"]:.6g}/{tb["time"]:.6g}: ' + ', '.join(out))
elif __name__ == '__main__' and sys.argv[1] == 'steadyf':
    ta, tb = tab(sys.argv[2]), tab(sys.argv[3])
    print(f"t={ta['time']:.6g}->{tb['time']:.6g}: " + ', '.join(
        f"{c} {np.max(np.abs(ta[c]-tb[c]))/np.max(np.abs(tb[c])):.2e}"
        for c in ('m1_e', 'm1_f1')))
