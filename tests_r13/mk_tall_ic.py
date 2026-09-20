#!/usr/bin/env python3
"""tests_r13: extend ic_he4_presn_sph.txt above its last radius (1.0256 R) with a thin
FLOOR ATMOSPHERE so x1max can sit at ~1.25 R: rho decays from the file's last value to
RHO_ATM over 0.005 R and stays there, T = the file's last T (4e4 K class, below the EOS
temperature gate).  Not hydrostatic: it rains on the star (4e-5 of the envelope mass)."""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r12')
import relax_ic as R
RHO_ATM, RTOP = 3.0e-13, 1.32*R.RSTAR
d = np.loadtxt(R.IC0, comments='#')
r, rho, e = d[:, 0], d[:, 1], d[:, 2]
Tend = float(R.temp_of(rho[-1:], e[-1:])[0])
rn = np.arange(r[-1] + 2.0e7, RTOP, 2.0e7)
x = (rn - r[-1])/(0.005*R.RSTAR)
rhon = RHO_ATM + (rho[-1] - RHO_ATM)*np.exp(-x)
Tn = np.full_like(rn, Tend)
en = R.egas(rhon, Tn) + R.taper_w(rhon, Tn)*R.A_RAD*Tn**4
# continuity of eint at the join: scale the offline value onto the file's last entry
e0 = (R.egas(rho[-1:], np.array([Tend])) + R.taper_w(rho[-1:], np.array([Tend]))*R.A_RAD*Tend**4)[0]
en *= e[-1]/e0
out = 'ic_he4_tall.txt'
with open(out, 'w') as fh:
    fh.write('# tests_r13/mk_tall_ic.py: ic_he4_presn_sph.txt + floor atmosphere rho=%g above 1.0256 R, T=%.0f K\n' % (RHO_ATM, Tend))
    for a, b, c in zip(np.concatenate((r, rn)), np.concatenate((rho, rhon)), np.concatenate((e, en))):
        fh.write('%.10e %.10e %.10e\n' % (a, b, c))
print('T_end = %.0f K, e scale %.4f, rows %d -> %d, top r = %.4f R, atmosphere mass = %.2e g'
      % (Tend, e[-1]/e0, len(r), len(r) + len(rn), rn[-1]/R.RSTAR, np.trapz(4*np.pi*rn**2*rhon, rn)))
