#!/usr/bin/env python3
"""cmpout.py REF NEW: compare the outputs of two run dirs of the same deck (fin4-1014).
bin/*.bin: bytes after the <par_end> line (the parameter header may differ: FIN4 adds its default keys);
tab/*.tab: bytes after the leading '#' lines; *.hst: the non-'#' lines.  Prints one line:
  BITWISE nbin ntab nhst | DIFF nbin ntab nhst ndiff <first differing file> maxrel=<max over the differing numeric data>
  | MISSING <files present in REF but not in NEW>"""
import sys, os, glob
import numpy as np
R, N = sys.argv[1], sys.argv[2]
def after_par(b):
    i = b.find(b'<par_end>')
    return b if i < 0 else b[b.index(b'\n', i) + 1:]
def tabdata(b):
    return b'\n'.join(l for l in b.split(b'\n') if not l.startswith(b'#'))
def nums(b):
    out = []
    for l in b.decode('latin-1').split('\n'):
        if l.startswith('#'): continue
        for w in l.split():
            try: out.append(float(w))
            except ValueError: pass
    return np.array(out)
files = [os.path.relpath(f, R) for pat in ('bin/*.bin', 'tab/*.tab', '*.hst') for f in glob.glob(os.path.join(R, pat))]
miss = [f for f in files if not os.path.exists(os.path.join(N, f))]
cnt = {'bin': 0, 'tab': 0, 'hst': 0}; diff = []; mr = 0.0
for f in sorted(set(files) - set(miss)):
    a = open(os.path.join(R, f), 'rb').read(); b = open(os.path.join(N, f), 'rb').read()
    k = f.split('.')[-1]; cnt[k] += 1
    if k == 'bin':
        da, db = after_par(a), after_par(b)
        if da != db:
            diff.append(f)
            try:
                sys.path.insert(0, os.environ.get('BINCONV', ''))
                import bin_convert
                xa = bin_convert.read_binary(os.path.join(R, f)); xb = bin_convert.read_binary(os.path.join(N, f))
                for v in xa['var_names']:
                    x = np.concatenate([np.ravel(q) for q in xa['mb_data'][v]]); y = np.concatenate([np.ravel(q) for q in xb['mb_data'][v]])
                    s = np.maximum(np.abs(x), np.abs(y)); mr = max(mr, float(np.max(np.where(s > 0, np.abs(x - y) / np.where(s > 0, s, 1), 0))))
            except Exception:
                mr = float('inf')   # bin_convert not found (set BINCONV=<athenak>/vis/python): maxrel unknown
    else:
        da, db = tabdata(a), tabdata(b)
        if da != db:
            diff.append(f); x, y = nums(da), nums(db)
            if x.size == y.size and x.size:
                s = np.maximum(np.abs(x), np.abs(y)); mr = max(mr, float(np.max(np.where(s > 0, np.abs(x - y) / np.where(s > 0, s, 1), 0))))
            else: mr = float('inf')
n = '%d %d %d' % (cnt['bin'], cnt['tab'], cnt['hst'])
if miss: print('MISSING %s %s' % (n, ' '.join(miss[:5])))
elif diff: print('DIFF %s ndiff=%d first=%s maxrel=%.3e' % (n, len(diff), diff[0], mr))
else: print('BITWISE %s' % n)
