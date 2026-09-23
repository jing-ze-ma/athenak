# usage: python3 hstcmp.py dirA dirB : max relative difference per hst column (same rows)
import sys, glob, numpy as np
a = sorted(glob.glob(sys.argv[1] + '/*.hst'))
for fa in a:
    fb = fa.replace(sys.argv[1], sys.argv[2], 1)
    A = np.loadtxt(fa); B = np.loadtxt(fb)
    n = min(len(A), len(B)); A = A[:n]; B = B[:n]
    hdr = [l for l in open(fa) if l.startswith('#')][-1].split()
    d = np.abs(A - B) / np.maximum(np.abs(A).max(axis=0), 1e-300)
    m = d.max(axis=0)
    idx = np.argsort(m)[::-1][:6]
    print(fa.split('/')[-1], 'rows', n, ' '.join(f'{hdr[i+1] if i+1 < len(hdr) else i}:{m[i]:.2e}' for i in idx))
