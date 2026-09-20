#!/usr/bin/env python3
"""The FIRST resumed history line of a chained run against the continuous run's line at
the nearest time: the radial kinetic energy is what the closure re-seed kicks."""
import sys


def load(fn):
    return [[float(x) for x in ln.split()] for ln in open(fn) if not ln.startswith('#')]


a, b = load(sys.argv[1]), load(sys.argv[2])
NAME = ['time', 'dt', 'mass', '1-mom', '2-mom', '3-mom', 'tot-E', '1-KE', '2-KE', '3-KE']
for n in range(min(3, len(b))):
    rb = b[n]
    ra = min(a, key=lambda r: abs(r[0]-rb[0]))
    print('# resumed line %d, t = %.10g (continuous line at t = %.10g)'
          % (n, rb[0], ra[0]))
    for c in range(1, min(len(ra), len(rb))):
        d = abs(ra[c]-rb[c])/max(abs(ra[c]), abs(rb[c]), 1e-300)
        print('    %-7s cont %+.10e  rst %+.10e  rel %.3e'
              % (NAME[c], ra[c], rb[c], d))
