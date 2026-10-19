#!/usr/bin/env python3
"""usage: budget.py <user.hst> [t_start] : mean rates over t >= t_start (last run block)."""
import sys
import numpy as np
fn = sys.argv[1]
t0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
rows = []
for line in open(fn):
    if line.startswith('# Athena'):
        rows = []          # keep the last run only (hst files are appended)
    elif not line.startswith('#'):
        rows.append([float(x) for x in line.split()])
a = np.array(rows)
t = a[:, 0]
s = t >= t0
# columns: 0 t 1 dt 2 Mdom 3 Jdom 4 Min 5 Mout 6 Macc 7 Jacc 8 Jstr 9 Jout, 10.. rates
dt = t[s][-1] - t[s][0]
d = a[s][-1] - a[s][0]
GMa = 6.674e-8*1.989e33/(6.957e10*1e10)*6.24
jk = np.sqrt(GMa*4.06)
print(f"{fn}: t {t[s][0]:.4f}..{t[s][-1]:.4f} (P_orb 0.8519)")
print(f"  Mdot_in {d[4]/dt:.4g} Mdot_out {d[5]/dt:.4g} Mdot_acc {d[6]/dt:.4g} "
      f"Mdom {a[-1,2]:.4g}")
print(f"  Jdot_acc {d[7]/dt:.4g}  j_acc = {d[7]/d[6]:.1f} Rsun km/s = {d[7]/d[6]/jk:.4f} "
      f"j_Kep(R_acc)  ; wall-stress part Jstr/Jacc = {d[8]/d[7]:.4f}; Jout/dt {d[9]/dt:.4g}")
print(f"  last-interval rates: dMacc {a[-1,12]:.4g} dJacc {a[-1,13]:.4g} dJstr {a[-1,14]:.4g}")
