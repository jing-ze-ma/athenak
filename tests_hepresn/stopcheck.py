#!/usr/bin/env python3
"""exit 1 (stop) if RUNDIR's user.hst shows NaN in any row, a dt jump, dt < 0.1 s,
M_tot < 0.5 M0, or time >= TLIM.  usage: stopcheck.py RUNDIR TLIM"""
import sys

import numpy as np

run, tlim = sys.argv[1], float(sys.argv[2])
try:
    d = np.loadtxt(run + '/hepresn.user.hst', comments='#', ndmin=2)
except OSError:
    sys.exit(0)
t, dt, m = d[-1, 0], d[-1, 1], d[-1, 15]
# any non-finite row (not only the last), or a dt jump > 20x the median of the
# preceding 50 rows (a blown-up state runs away in dt before the next check)
jump = len(d) > 50 and dt > 20 * np.median(d[-51:-1, 1])
bad = (not np.all(np.isfinite(d)) or jump or dt < 0.1 or m < 0.5 * d[0, 15]
       or t >= tlim - 1.0)
print('stopcheck %s: t %.0f dt %.3g M/M0 %.4f -> %s'
      % (run, t, dt, m / d[0, 15], 'STOP' if bad else 'go'))
sys.exit(1 if bad else 0)
