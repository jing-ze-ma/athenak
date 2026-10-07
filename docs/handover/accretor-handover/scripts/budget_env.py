#!/usr/bin/env python3
"""Stage-2 torque budget. usage: budget_env.py <user.hst> [t_start]
Means over t >= t_start of the cumulative boundary integrals (ry_per_accretor, envelope)."""
import sys
import numpy as np
fn = sys.argv[1]
t0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
rows, lab = [], None
for line in open(fn):
    if line.startswith('# Athena'):
        rows = []
    elif line.startswith('#'):
        lab = [w.split('=')[1] for w in line[1:].split() if '=' in w]
    else:
        rows.append([float(x) for x in line.split()])
a = np.array(rows)
c = {n: i for i, n in enumerate(lab)}
t = a[:, c['time']]
s = t >= t0
dt = t[s][-1] - t[s][0]
d = {n: (a[s][-1, c[n]] - a[s][0, c[n]])/dt for n in lab}
GMa = 6.674e-8*1.989e33/(6.957e10*1e10)*6.24
jk = np.sqrt(GMa*4.06)
print(f"{fn}: t {t[s][0]:.4f}..{t[s][-1]:.4f} (P_orb 0.8519), j_Kep(R_acc) {jk:.1f}")
print(f"  Mdot_in {d['Min']:.4g}  Mdot_out {d['Mout']:.4g}  Mdot through R_acc {d['MR']:.4g}"
      f"  dMenv/dt {d['Menv']:.4g}  Mwal {d['Mwal']:.3g}")
print(f"  Jdot_in(stream) {d['Jin']:.4g} (j_in {d['Jin']/d['Min']/jk:.4f} j_K)  "
      f"Jdot through R_acc {d['JR']:.4g} (j {d['JR']/d['MR']/jk:.4f} j_K per accreted mass)")
print(f"  dJenv/dt {d['Jenv']:.4g}  Jwal-rate {d['Jwal']:.3g}  Jout-rate {d['Jout']:.4g}"
      f"  dJdom/dt {d['Jdom']:.4g}")
print(f"  JR / (Mdot_in j_K) = {d['JR']/d['Min']/jk:.4f}")
