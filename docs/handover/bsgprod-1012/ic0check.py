#!/usr/bin/env python3
"""ic0check.py make|check RUN_DIR REF_JSON  (bsgrec-1011)
Cycle-0 IC fingerprint of a BSG run: bsg3d.{hydro_w,m1}.00000.bin in RUN_DIR/bin.  Per variable of the assembled
(nx3, nx2, nx1) float32 array: md5 of the bytes (layout-independent: blocks placed by their geometry), float64
sum/min/max, and the shell mean over (k, j) per radial index.  'make' writes REF_JSON; 'check' compares and prints,
per variable: BITWISE (md5 equal) or the max relative difference of the shell means and of sum/min/max.
Needs numpy and bin_convert_repro.py (docs/handover/bsg-hrdet-1009/page/, put that dir on PYTHONPATH)."""
import sys, os, json, hashlib
import numpy as np
import bin_convert_repro as bc


def assemble(d):
    geo = np.asarray(d['mb_geometry'])
    nx1, nx2, nx3 = d['Nx1'], d['Nx2'], d['Nx3']
    b1, b2, b3 = d['nx1_mb'], d['nx2_mb'], d['nx3_mb']
    out = {}
    for v in d['var_names']:
        a = np.empty((nx3, nx2, nx1), np.float32)
        for m, blk in enumerate(d['mb_data'][v]):
            i0 = int(round((geo[m, 0] - d['x1min']) / (d['x1max'] - d['x1min']) * nx1))
            j0 = int(round((geo[m, 2] - d['x2min']) / (d['x2max'] - d['x2min']) * nx2))
            k0 = int(round((geo[m, 4] - d['x3min']) / (d['x3max'] - d['x3min']) * nx3))
            a[k0:k0 + b3, j0:j0 + b2, i0:i0 + b1] = blk
        out[v] = a
    return out


def fp(run):
    res = {}
    for kind in ('hydro_w', 'm1'):
        d = bc.read_binary(os.path.join(run, 'bin/bsg3d.%s.00000.bin' % kind))
        w = assemble(d)
        for v, a in w.items():
            a = np.ascontiguousarray(a, np.float32)
            a64 = a.astype(np.float64)
            res[v] = dict(md5=hashlib.md5(a.tobytes()).hexdigest(), sum=float(a64.sum()),
                          min=float(a64.min()), max=float(a64.max()),
                          shell=[float(x) for x in a64.mean(axis=(0, 1))])
    return res


def rel(x, y):
    x, y = np.asarray(x), np.asarray(y)
    s = np.maximum(np.abs(x), np.abs(y))
    return float(np.max(np.where(s > 0, np.abs(x - y) / np.where(s > 0, s, 1), 0.0)))


mode, run, ref = sys.argv[1:4]
if mode == 'make':
    json.dump(fp(run), open(ref, 'w'))
    print('wrote', ref)
else:
    R = json.load(open(ref)); N = fp(run)
    for v in R:
        if v not in N:
            print('%-8s MISSING' % v); continue
        if R[v]['md5'] == N[v]['md5']:
            print('%-8s BITWISE' % v); continue
        print('%-8s differs: max rel shell-mean %.3e sum %.3e min %.3e max %.3e' % (
            v, rel(R[v]['shell'], N[v]['shell']), rel(R[v]['sum'], N[v]['sum']),
            rel(R[v]['min'], N[v]['min']), rel(R[v]['max'], N[v]['max'])))
