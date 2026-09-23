#!/usr/bin/env python3
"""Shock displacement implied by each conserved total: python3 cons.py <run dir> [...]
For a steady shock held in its own frame the end fluxes balance, so the domain totals of
mass and total energy (eint + rho v^2/2 + E_rad) can only change by the shock moving:
dU = (U_up - U_down) * s.  Printed: s/dx from each total at the last dump, the x1
momentum budget (rho v is continuous across the shock, so its total changes only by an
end-flux imbalance: dP/(t * rho v^2 upstream)), and the boundary-cell E_rad relative to
the t = 0 end states.  If s(mass) and s(energy) disagree, energy is not conserved."""
import glob
import os
import sys
import numpy as np


def load(f):
    return np.loadtxt(f)


for d in sys.argv[1:]:
    tm = sorted(glob.glob(os.path.join(d, 'tab', '*.m1.*.tab')))
    th = sorted(glob.glob(os.path.join(d, 'tab', '*.hydro_w.*.tab')))
    if len(tm) < 2:
        continue
    out = []
    for f1, f2 in [(tm[0], th[0]), (tm[-1], th[-1])]:
        r = load(f1)
        h = load(f2)
        x = h[:, 2]
        rho, v, ei = h[:, 3], h[:, 4], h[:, 7]
        er = r[:, 3]
        ut = np.array([rho, rho*v, ei + 0.5*rho*v*v + er])
        out.append((x, ut, er))
    tt = float(open(tm[-1]).readline().split('time=')[1].split()[0])
    x, u0, e0 = out[0]
    _, u1, e1 = out[1]
    dx = x[1] - x[0]
    s = [(u1[q].sum() - u0[q].sum())*dx/(u0[q][0] - u0[q][-1])/dx for q in (0, 2)]
    pm = (u1[1].sum() - u0[1].sum())*dx/(tt*u0[1][0]**2/u0[0][0])
    print('%-26s s/dx mass %+.3f energy %+.3f  mom-flux imbalance %+.2e | Erad end-cell'
          ' rel: lo %+.2e hi %+.2e' % (os.path.basename(d.rstrip('/')), s[0], s[1], pm,
                                       e1[0]/e0[0]-1, e1[-1]/e0[-1]-1))
