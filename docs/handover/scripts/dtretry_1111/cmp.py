#!/usr/bin/env python3
# cmp.py A B : compare restart files of two run dirs (body after the input dump), per file
import sys, os, glob
import numpy as np
def body(p):
    b = open(p, 'rb').read(); i = b.find(b'<par_end>'); return b[i:] if i >= 0 else b
for fa in sorted(glob.glob(sys.argv[1] + '/rst/*.rst')):
    fb = sys.argv[2] + '/rst/' + os.path.basename(fa)
    if not os.path.exists(fb): print(os.path.basename(fa), 'MISSING'); continue
    a, b = body(fa), body(fb)
    if a == b: print(os.path.basename(fa), 'IDENTICAL'); continue
    n = min(len(a), len(b)) // 8 * 8
    x = np.frombuffer(a[len(a) % 8:], dtype=np.float64); y = np.frombuffer(b[len(b) % 8:], dtype=np.float64)
    m = min(len(x), len(y)); x, y = x[:m], y[:m]
    ok = np.isfinite(x) & np.isfinite(y)
    d = np.abs(x[ok] - y[ok]) / (np.abs(x[ok]) + np.abs(y[ok]) + 1e-300)
    print(os.path.basename(fa), 'DIFFER len %d/%d ndiff %d maxrel %.3e (8-byte words, header-aligned approx)'
          % (len(a), len(b), int((x != y).sum()), d.max() if d.size else 0))
