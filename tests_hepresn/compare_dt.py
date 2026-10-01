#!/usr/bin/env python3
"""compare_dt.py RUN_A RUN_B R_INT : force f/(rho g) = [rho v_r(t) - rho v_r(0)]/t/(rho g)
by zone at (nearly) the same times in two column runs that differ only in the time step
(cfl 0.3 vs 0.15); a dt-linear force (split-lag bug) shows as a ratio of 2 between the
two at equal t.  Also the same for L(r)/L(0) - 1."""
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_hepresn/tests_hepresn')
sys.path.insert(0, '/viper/ptmp2/jinma/wt_hepresn/vis/python')
import analyze_gate1a as A  # noqa: E402

a, b, rint = sys.argv[1], sys.argv[2], float(sys.argv[3])
ha, hb = A.load(a, 'hydro_u'), A.load(b, 'hydro_u')
ma, mb = A.load(a, 'm1'), A.load(b, 'm1')
rc = ha[0][2]
g = A.GM / rc**2
rho0 = ha[0][3]['dens']
zones = [(0.0, 0.8), (0.8, 0.9), (0.9, 0.97), (0.97, 9.9)]
below = rc < rint
below[0] = False
for tt in (13.4, 40.3, 100.0, 200.0, 299.9):
    ia = int(np.argmin([abs(h[0] - tt) for h in ha]))
    ib = int(np.argmin([abs(h[0] - tt) for h in hb]))
    print('t_A = %.2f (dt %.2f)  t_B = %.2f (dt %.2f)' %
          (ha[ia][0], ha[ia][0] / max(ia, 1), hb[ib][0], hb[ib][0] / max(ib, 1)))
    for lo, hi in zones:
        z = below & (rc / A.RSTAR >= lo) & (rc / A.RSTAR < hi)
        fa = (ha[ia][3]['mom1'] - ha[0][3]['mom1']) / ha[ia][0] / (rho0 * g)
        fb = (hb[ib][3]['mom1'] - hb[0][3]['mom1']) / hb[ib][0] / (rho0 * g)
        La = 4 * np.pi * rc**2 * ma[ia][3]['m1_f1'] / \
            (4 * np.pi * rc**2 * ma[0][3]['m1_f1']) - 1
        Lb = 4 * np.pi * rc**2 * mb[ib][3]['m1_f1'] / \
            (4 * np.pi * rc**2 * mb[0][3]['m1_f1']) - 1
        print('   %.2f-%.2f R: max|f| A %.2e B %.2e ratio %.2f | max|dL/L| A %.2e B %.2e'
              % (lo, min(hi, 1.0), np.abs(fa[z]).max(), np.abs(fb[z]).max(),
                 np.abs(fa[z]).max() / max(np.abs(fb[z]).max(), 1e-30),
                 np.abs(La[z]).max(), np.abs(Lb[z]).max()))
