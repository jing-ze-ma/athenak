#!/usr/bin/env python3
"""Risk 6 (T-S4 with gas momentum on): a static pure-scattering atmosphere in radiative
equilibrium, restarted with the gas free to move.  Per run: rms and max |v| over all
cells for each hydro_w bin dump (every 10 cycles), and the time.  With the m1 tab at the
restart state: the radiation force the coupling deposits per unit mass, (kappa_t/c)
<F_f>_cell (uniform rho kappa_t: the mean of the two face fluxes = the cell F1), against
the well-balanced reference kappa_t F_in r_in^2/(c r^2) of force_reference = wb_arad.

usage: tr6.py F_in r_in r_out m1tab rundir [rundir ...]
"""
import glob
import sys

import numpy as np

import s2lib

fin, rin, rout = map(float, sys.argv[1:4])
d = s2lib.load_tab(sys.argv[4])
f1 = d['m1_f1']
n = f1.shape[0]
rf, x1v = s2lib.rgrid(rin, rout, n, None)
res = f1*x1v*x1v/(fin*rin*rin) - 1.0
print(f"force residual (deposit/wb_arad reference - 1): max {np.max(np.abs(res)):.3e} "
      f"rms {np.sqrt(np.mean(res**2)):.3e} (n={n}); cell-centred r^2 F/(r_in^2 F_in)-1 "
      f"inner {res[0]:+.3e} outer {res[-1]:+.3e}")
for rd in sys.argv[5:]:
    fs = sorted(glob.glob(rd + '/bin/*.bin'))
    out = []
    for f in fs:
        h = s2lib.load(f)
        v2 = h['velx']**2 + h['vely']**2 + h['velz']**2
        out.append((h['time'], h['cycle'], np.sqrt(v2.mean()), np.sqrt(v2.max())))
    t0 = out[0][0]
    s = ' '.join(f"{c}:{r:.2e}" for _, c, r, _ in out[::2])
    print(f"{rd}: cycles {out[0][1]}..{out[-1][1]} dt_total={out[-1][0] - t0:.4g} "
          f"final rms|v|={out[-1][2]:.3e} max|v|={out[-1][3]:.3e}  rms by cycle: {s}")
