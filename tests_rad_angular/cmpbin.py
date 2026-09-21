#!/usr/bin/env python3
"""Bitwise comparison of the DATA payload of two AthenaK binary dumps.

The header of a .bin dump embeds the full parameter input, so two runs that differ
only in an input key have different headers while the physics payload may be
identical.  This compares the decoded variable arrays instead.
"""
import sys
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import numpy as np
import bin_convert


def payload(fn):
    d = bin_convert.read_binary(fn)
    out = []
    for mb in d['mb_data'].values():
        out.append(np.asarray(mb).ravel())
    return np.concatenate(out)


a, b = sys.argv[1], sys.argv[2]
pa, pb = payload(a), payload(b)
same = (pa.tobytes() == pb.tobytes())
if same:
    print('%s  %s  BITWISE IDENTICAL payload (%d values)' % (a, b, pa.size))
else:
    d = np.abs(pa - pb)
    r = d / np.maximum(np.abs(pa), 1e-300)
    print('%s  %s  DIFFER: max|d| = %.6e  max rel = %.6e' % (a, b, d.max(), r.max()))
