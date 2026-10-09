"""hrup-1009: the grey atmosphere of t9_atmosphere.py against the EXACT TRANSPORT (Hopf)
solution, the right reference for a VET (ray) closure: E(tau) = 3 (F/c) (tau + q(tau)),
q(tau) the Hopf function (q(0) = 1/sqrt 3, q(inf) = 0.710446; fit 0.710446 - 0.133054
exp(-3.4488 tau), |err| < 2e-3).  tau from the top face as in t9: kappa rho_top H
(exp[(ztop - z)/H] - 1).  Prints max|dE/E| over the column, over the top 5 cells, L1.
usage: t9h.py <tab> [--kappa 1 --rho-top 0.128 --scale-h 0.113 --flux 1 --c 1 --ztop 1]"""
import argparse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('tab')
ap.add_argument('--kappa', type=float, default=1.0)
ap.add_argument('--rho-top', type=float, default=0.128)
ap.add_argument('--scale-h', type=float, default=0.113)
ap.add_argument('--flux', type=float, default=1.0)
ap.add_argument('--c', type=float, default=1.0)
ap.add_argument('--ztop', type=float, default=1.0)
ap.add_argument('--label', default='')
a = ap.parse_args()
d = np.loadtxt(a.tab)
z, E = d[:, 2], d[:, 3]
tau = a.kappa*a.rho_top*a.scale_h*(np.exp((a.ztop - z)/a.scale_h) - 1.0)
q = 0.710446 - 0.133054*np.exp(-3.4488*tau)
Ex = 3.0*a.flux/a.c*(tau + q)
r = np.abs(E/Ex - 1.0)
print('%s HOPF: max|dE/E| %.3e  top5 %.3e  L1 %.3e  (E_top %.4f exact %.4f, tau_top-cell %.2e)'
      % (a.label, r.max(), r[-5:].max(), np.abs(E - Ex).sum()/np.abs(Ex).sum(), E[-1], Ex[-1],
         tau[-1]))
