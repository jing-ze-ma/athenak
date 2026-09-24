#!/usr/bin/env python3
"""max over i of the transverse (j,k) spread of E / mean, per m1 bin dump (every 4th)"""
import glob
import sys
import numpy as np
import binlib
fs = sorted(glob.glob(sys.argv[1] + '/bin/*m1*.bin'))
out = []
for f in fs:
    h = binlib.load(f)
    e = h['m1_e']
    s = np.max((e.max(axis=(0, 1)) - e.min(axis=(0, 1)))/e.mean(axis=(0, 1)))
    out.append(f"{h['cycle']}:{s:.1e}")
print('spread(E)/E by cycle ' + ' '.join(out[::4] + ([out[-1]] if (len(out)-1) % 4 else [])))
