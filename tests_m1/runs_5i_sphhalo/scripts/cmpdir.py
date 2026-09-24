"""cmpdir.py A B: every output file of run A vs run B.  bitwise (rst/bin: payload after
<par_end> resp. the whole file; tab/hst: the numbers) and, for tab/hst, the max relative
difference per column (against the column's max |a|)."""
import glob
import os
import sys

import numpy as np

a, b = sys.argv[1], sys.argv[2]
fs = sorted(glob.glob(a + "/**/*.*", recursive=True))
fs = [f for f in fs if f.endswith((".tab", ".hst", ".rst", ".bin"))]
nbit, worst, wat = 0, 0.0, ""
for fa in fs:
    fb = b + fa[len(a):]
    if not os.path.exists(fb):
        print("missing", fb)
        continue
    ra, rb = open(fa, "rb").read(), open(fb, "rb").read()
    if fa.endswith((".rst", ".bin")):
        same = ra.split(b"<par_end>", 1)[-1] == rb.split(b"<par_end>", 1)[-1]
    else:
        da = np.loadtxt(fa, comments="#", ndmin=2)
        db = np.loadtxt(fb, comments="#", ndmin=2)
        same = da.shape == db.shape and np.array_equal(da, db)
        if da.shape == db.shape:
            s = np.max(np.abs(da), axis=0)
            s[s == 0] = 1.0
            r = np.max(np.abs(da - db), axis=0)/s
            if r.max() > worst:
                worst, wat = r.max(), os.path.basename(fa) + " col %d" % r.argmax()
    nbit += same
print("%s vs %s: %d/%d files bitwise; worst tab/hst rel diff %.2e (%s)"
      % (os.path.basename(a), os.path.basename(b), nbit, len(fs), worst, wat))
