"""sp_pairs.py A B: final restart state of two runs: max rel diff of dens, gas E, E_rad,
F_rad (vector, against max |F|), and whether the payload is bitwise."""
import glob
import sys

import numpy as np

sys.path.insert(0, "/viper/ptmp2/jinma/wt_sphhalo/tests_m1/gates")
import cmp  # noqa: E402

a, b = sys.argv[1], sys.argv[2]
fa = sorted(glob.glob(a + "/rst/*.rst"))[-1]
fb = b + fa[len(a):]
sa, sb = cmp._rst_state(fa), cmp._rst_state(fb)
out = []
for n, g in (("dens", [0]), ("Egas", [4]), ("Erad", [5]), ("Frad", [6, 7, 8])):
    s = np.max(np.sqrt(sum(sa[:, i]**2 for i in g)))
    d = np.max(np.sqrt(sum((sa[:, i] - sb[:, i])**2 for i in g)))
    out.append("%s %.1e" % (n, d/s))
pa = open(fa, "rb").read().split(b"<par_end>", 1)[-1]
pb = open(fb, "rb").read().split(b"<par_end>", 1)[-1]
print("  ".join(out) + "  bitwise %s" % (pa == pb))
