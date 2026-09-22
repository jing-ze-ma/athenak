#!/usr/bin/env python
"""Max relative difference, column dump against column dump.

    python coldiff.py A/col.txt B/col.txt

Compares the correlated-k column dump's physical columns -- T, the NET longwave face
flux, the longwave deposit and the stellar heating -- and prints, for each, the largest
relative difference and where it sits.  Fluxes and deposits change sign down a column,
so they are normalised by the column maximum of |quantity| rather than pointwise.
"""
import sys
import numpy as np

COLS = {'T': 3, 'F_lw': 4, 'Q_sw': 5, 'Src_lw': 10}


def load(path):
    return np.loadtxt(path)


def main():
    a = load(sys.argv[1])
    b = load(sys.argv[2])
    if a.shape != b.shape:
        print('shape mismatch', a.shape, b.shape)
        return
    print('%-40s %-40s' % (sys.argv[1], sys.argv[2]))
    for name, c in COLS.items():
        x, y = a[:, c], b[:, c]
        scale = max(np.max(np.abs(x)), 1e-300)
        d = np.abs(x - y)
        i = int(np.argmax(d))
        rel = d[i]/scale
        ptw = np.max(d/np.maximum(np.abs(x), 1e-300))
        print('  %-7s max|dA|/max|A| = %10.3e  at i=%3d   max pointwise = %10.3e'
              % (name, rel, int(a[i, 0]), ptw))


if __name__ == '__main__':
    main()
