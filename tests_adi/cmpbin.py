"""Compare two directories of AthenaK binary dumps by DATA.

The dump embeds the effective input, and this branch adds parameters to it, so `cmp` on
the files always differs; what has to be identical is every variable of every dump.
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "vis", "python"))
import bin_convert  # noqa: E402

a, b = sys.argv[1], sys.argv[2]
fa = sorted(glob.glob(os.path.join(a, "*.bin")))
fb = sorted(glob.glob(os.path.join(b, "*.bin")))
assert len(fa) == len(fb) and fa, (len(fa), len(fb))
bad = 0
for x, y in zip(fa, fb):
    ra, rb = bin_convert.read_binary(x), bin_convert.read_binary(y)
    for k in ra["mb_data"]:
        u, v = np.asarray(ra["mb_data"][k]), np.asarray(rb["mb_data"][k])
        if not np.array_equal(u, v):
            bad += 1
            print("DIFF", os.path.basename(x), k, "max|d| =", np.abs(u - v).max())
    if ra["time"] != rb["time"]:
        bad += 1
        print("DIFF time", os.path.basename(x), ra["time"], rb["time"])
print(("BITWISE IDENTICAL over %d dumps" % len(fa)) if bad == 0
      else "%d DIFFERENCES" % bad)
