"""max relative difference of the history files of runs against a reference.
usage: python3 hstdiff.py REFDIR RUNDIR [RUNDIR ...]"""
import glob
import os
import sys
import numpy as np

ref = sys.argv[1]
for d in sys.argv[2:]:
    out = []
    for f in sorted(glob.glob(ref + '/*.hst')):
        a = np.loadtxt(f)
        b = np.loadtxt(os.path.join(d, os.path.basename(f)))
        n = min(len(a), len(b))
        a, b = a[:n], b[:n]
        sc = np.maximum(np.abs(a).max(axis=0), 1e-300)
        rel = (np.abs(a - b) / sc).max(axis=0)
        out.append('%s rows %d: max rel diff (col-scaled) %.2e, last row %.2e'
                   % (os.path.basename(f), n, rel[1:].max(),
                      (np.abs(a[-1] - b[-1]) / sc)[1:].max()))
    print(d + ': ' + ' | '.join(out))
