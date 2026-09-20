#!/usr/bin/env python3
"""Compare two red_giant history files line by line.

The restarted leg rewrites the history line at the restart time, so the comparison is
made on the (time -> line) map and the duplicate is dropped.
"""
import sys


def load(fn):
    d = {}
    for ln in open(fn):
        if ln.startswith('#'):
            continue
        f = ln.split()
        d[f[0]] = [float(x) for x in f]
    return d


a, b = load(sys.argv[1]), load(sys.argv[2])
common = sorted(set(a) & set(b), key=float)
print('# %d lines in A, %d in B, %d shared times' % (len(a), len(b), len(common)))
nbit, worst, wt, wc = 0, 0.0, '', 0
for t in common:
    ra, rb = a[t], b[t]
    if ra == rb:
        nbit += 1
        continue
    for c, (x, y) in enumerate(zip(ra, rb)):
        e = abs(x-y)/max(abs(x), abs(y), 1e-300)
        if e > worst:
            worst, wt, wc = e, t, c+1
print('# bitwise-identical shared lines: %d / %d' % (nbit, len(common)))
print('# worst relative difference: %.3e at t = %s, column %d' % (worst, wt, wc))
