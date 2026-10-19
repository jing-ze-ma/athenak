#!/usr/bin/env python3
"""Join two TOPS fetches in T (0.0005-1 keV and 1-10 keV, same rho grid); the 1 keV rows of
the second file are dropped (checked identical to the first file's).  This is the step done
for the X 0.70 tables (10-09).
usage: merge_tops_T.py LOW_T.dat HIGH_T.dat OUT.dat"""
import sys
import numpy as np

lo, hi, out = sys.argv[1:4]
L = open(lo).read().splitlines()
H = open(hi).read().splitlines()
a = np.loadtxt(lo)
b = np.loadtxt(hi)
assert np.array_equal(a[a[:, 0] == 1.0], b[b[:, 0] == 1.0]), '1 keV rows differ'
body = ([s for s in L if not s.startswith('#')]
        + [s for s in H if not s.startswith('#') and not s.startswith('1.000000e+00 ')])
hdr = [s for s in L if s.startswith('#')]
hdr.append('# merged: %s (T 0.0005-1 keV) + %s (T 1-10 keV, its 1 keV rows dropped: identical)'
           % (lo, hi))
open(out, 'w').write('\n'.join(hdr + body) + '\n')
d = np.loadtxt(out)
print(out, 'nT', len(np.unique(d[:, 0])), 'nrho', len(np.unique(d[:, 1])))
