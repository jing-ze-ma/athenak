"""Compare two directories of AthenaK binary dumps cycle by cycle.

Usage: python3 cmpdumps.py <dir A> <dir B>

Dumps are matched by the cycle number in their own header, not by file number:
the restarted run keeps counting output files from where the first leg stopped.
For every cycle present in both directories the DATA of the two files are
compared bit for bit (the text header is not compared: it embeds the parameter
dump, which legitimately differs after a restart -- the command line and every
output block's last_time are not the same strings).  Where the data differ, the
largest absolute and relative difference over all variables is reported, so that
a failing gate carries a number and not just a verdict.
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'vis', 'python'))
import bin_convert as bc              # noqa: E402


def index(d):
    out = {}
    for f in sorted(glob.glob(os.path.join(d, '*.bin'))):
        out[int(bc.read_binary(f)['cycle'])] = f
    return out


def main():
    ia, ib = index(sys.argv[1]), index(sys.argv[2])
    cycles = sorted(set(ia) & set(ib))
    if not cycles:
        print('NO COMMON CYCLES')
        return 2
    first_bad = None
    for c in cycles:
        da = bc.read_binary(ia[c])['mb_data']
        db = bc.read_binary(ib[c])['mb_data']
        same = all(np.array_equal(np.asarray(da[k]), np.asarray(db[k])) for k in da)
        if same:
            print('cycle %3d  bitwise identical' % c)
            continue
        worst_a, worst_r, wname = 0.0, 0.0, ''
        for k in da:
            x = np.asarray(da[k], dtype=np.float64)
            y = np.asarray(db[k], dtype=np.float64)
            ad = np.abs(x - y)
            sc = np.maximum(np.abs(x), np.abs(y))
            sc[sc == 0.0] = 1.0
            if ad.max() > worst_a:
                worst_a, worst_r, wname = ad.max(), (ad / sc).max(), k
        print('cycle %3d  DIFFERS: worst %s  abs %.3e  rel %.3e'
              % (c, wname, worst_a, worst_r))
        if first_bad is None:
            first_bad = c
    if first_bad is None:
        print('PASS: all %d common cycles bitwise identical' % len(cycles))
        return 0
    print('FAIL: first differing cycle %d of %d compared' % (first_bad, len(cycles)))
    return 1


if __name__ == '__main__':
    sys.exit(main())
