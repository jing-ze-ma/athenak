"""Compare two Athena history (.hst) files column by column.

Usage: python3 cmp.py A.hst B.hst
Prints, per column, the maximum absolute and maximum relative difference over the
rows both files share.  Used by the milestone 3b2 phase-A decomposition gates.
"""
import sys

import numpy as np


def rd(fname):
    rows = []
    for line in open(fname):
        if line.startswith('#'):
            continue
        rows.append([float(x) for x in line.split()])
    return np.array(rows)


def main():
    a = rd(sys.argv[1])
    b = rd(sys.argv[2])
    n = min(len(a), len(b))
    a = a[:n]
    b = b[:n]
    d = np.abs(a - b)
    s = np.maximum(np.abs(a), np.abs(b))
    rel = np.where(s > 0, d/np.maximum(s, 1e-300), 0.0)
    col = int(np.unravel_index(rel.argmax(), rel.shape)[1])
    print("rows %d  max abs diff %g  max rel %g  (col %d)"
          % (n, d.max(), rel.max(), col))
    for c in range(a.shape[1]):
        print("  col %d: maxrel %.3e  maxabs %.3e" % (c, rel[:, c].max(), d[:, c].max()))


if __name__ == '__main__':
    main()
