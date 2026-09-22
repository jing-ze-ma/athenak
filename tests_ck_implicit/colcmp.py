#!/usr/bin/env python3
"""Gate (d)/(f): max |T_a - T_b| by pressure decade between two ck_dump_file columns."""
import sys
import numpy as np


def load(path):
    p, t = [], []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.split()
        try:
            v = [float(x) for x in f]
        except ValueError:
            continue
        if len(v) < 4:
            continue
        p.append(v[2])
        t.append(v[3])
    return np.array(p), np.array(t)


def main():
    pa, ta = load(sys.argv[1])
    pb, tb = load(sys.argv[2])
    n = min(len(ta), len(tb))
    pa, ta, tb = pa[:n], ta[:n], tb[:n]
    bins = [(10., 1e9), (1., 10.), (1e-2, 1.), (1e-4, 1e-2), (0., 1e-4)]
    out = []
    for lo, hi in bins:
        m = (pa > lo) & (pa <= hi)
        out.append('%.2f' % np.max(np.abs(ta[m] - tb[m])) if m.any() else '-')
    print('\t'.join(out))


if __name__ == '__main__':
    main()
