#!/usr/bin/env python3
"""Angular spread of E per radial index over the bin dumps of a wedge run:
max_i (max - min)/mean of E over (theta, phi) (c is unused, kept for the call form).
usage: spread.py c rundir [rundir ...]"""
import glob
import sys

import s2lib

c = float(sys.argv[1])
for rd in sys.argv[2:]:
    out = []
    for f in sorted(glob.glob(rd + '/bin/*.m1.*.bin')):
        d = s2lib.load(f)
        e = d['m1_e']
        sp = (e.max(axis=(0, 1)) - e.min(axis=(0, 1)))/e.mean(axis=(0, 1))
        out.append(f"{d['cycle']}:{sp.max():.1e}")
    print(rd, ' '.join(out))
