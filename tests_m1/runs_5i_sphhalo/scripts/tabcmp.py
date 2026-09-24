"""tabcmp.py A B: m1 tab slices: max rel diff of m1_e (vs max e) and of the flux
vector (vs max |F|) over every m1 tab file; also the transverse flux level
max|F2,F3|/max|F1|."""
import glob
import os
import sys

import numpy as np

a, b = sys.argv[1], sys.argv[2]
de = dfl = tr = 0.0
n = 0
for fa in sorted(glob.glob(a + "/tab/*.tab")):
    fb = b + fa[len(a):]
    da = np.loadtxt(fa, ndmin=2)
    db = np.loadtxt(fb, ndmin=2)
    e, f = da[:, 3], da[:, 4:7]
    de = max(de, np.max(np.abs(e - db[:, 3]))/np.max(np.abs(e)))
    fm = np.max(np.linalg.norm(f, axis=1))
    if fm > 0:
        dfl = max(dfl, np.max(np.linalg.norm(f - db[:, 4:7], axis=1))/fm)
        tr = max(tr, np.max(np.abs(db[:, 5:7]))/np.max(np.abs(f[:, 0])))
    n += 1
print("%s vs %s: %d m1 tabs, dE/E %.1e, dF/|F| %.1e, B transverse |F2,F3|/|F1| %.1e"
      % (os.path.basename(a), os.path.basename(b), n, de, dfl, tr))
