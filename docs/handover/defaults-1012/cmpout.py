#!/usr/bin/env python3
"""cmpout.py <dirA> <dirB> [maxdump]: bitwise comparison of two run directories: every bin dump
(arrays via bin_convert; the header holds the parameters, so only the data are compared), every
tab file (text after the header line) and every hst file (data rows).  Files present in both.
maxdump: only dumps with number <= maxdump (e.g. a longer reference run)."""
import os, sys, glob
import numpy as np
sys.path.insert(0, os.environ.get('ATHENA_VIS', '/raven/u/jinma/ATHENAK/wt-prodfinal/vis/python'))
import bin_convert
A, B = sys.argv[1], sys.argv[2]
mx = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9
mind = int(sys.argv[4]) if len(sys.argv) > 4 else -1


def rel(d, pat):
    return {os.path.relpath(f, d): f for f in glob.glob(os.path.join(d, pat), recursive=True)}


def num(f):
    try:
        return int(os.path.basename(f).split('.')[-2])
    except ValueError:
        return -1


n = bad = 0
for pat in ('**/*.bin', '**/*.tab'):
    fa, fb = rel(A, pat), rel(B, pat)
    for k in sorted(set(fa) & set(fb)):
        if num(k) > mx or num(k) < mind:
            continue
        n += 1
        if k.endswith('.bin'):
            da, db = bin_convert.read_binary(fa[k]), bin_convert.read_binary(fb[k])
            same = abs(da['time'] - db['time']) == 0
            for v in da['var_names']:
                for m in range(len(da['mb_data'][v])):
                    if not np.array_equal(da['mb_data'][v][m], db['mb_data'][v][m]):
                        same = False
        else:
            la = open(fa[k]).read().split('\n')[1:]
            lb = open(fb[k]).read().split('\n')[1:]
            same = la == lb
        if not same:
            bad += 1
            print('DIFF', k)
for k in sorted(set(rel(A, '**/*.hst')) & set(rel(B, '**/*.hst'))):
    # rows matched by their time column (a restarted run's history starts later)
    ra = {l.split()[0]: l for l in open(os.path.join(A, k)) if l.strip() and not l.startswith('#')}
    rb = {l.split()[0]: l for l in open(os.path.join(B, k)) if l.strip() and not l.startswith('#')}
    com = sorted(set(ra) & set(rb), key=float)
    n += 1
    nd = sum(ra[t] != rb[t] for t in com)
    if nd or not com:
        bad += 1
        print('DIFF', k, 'rows', len(com), 'differ', nd)
print('%s: %d files compared, %d differ -> %s' % (os.path.basename(A.rstrip('/')) + ' vs ' + B, n, bad,
                                                  'BITWISE' if bad == 0 and n > 0 else 'NOT BITWISE'))
