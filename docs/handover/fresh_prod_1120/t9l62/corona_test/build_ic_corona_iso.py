# build_ic_corona_iso.py (CORONA9 11-21): T9L6.2 corona IC variant: the corona of ic_t9l62_corona.txt (same join
# radius r_j, same rho and P at r_j) made ISOTHERMAL at T_j = T(r_j) (1.94 MK) instead of T ~ 1/r. Reason: the
# T ~ 1/r corona cools to 2.5e5 K at x1max, i.e. through log T 5.4-5.9 where the Planck table has kappa_P ~ 0.3-1
# (vs 1e-7 above 1.6 MK): the frozen-gas clean_init steady start then fills the corona with radiation (E_rad 8e46
# -> 3.3e49 erg on the GPU). Isothermal HSE: d ln rho/dr = -GM (1 - kappa_R L_IC/(4 pi c GM))/(r^2 cs_iso^2).
# Rows below r_j are unchanged; F_r, E, fmlt columns unchanged.
import numpy as np, sys
from scipy.interpolate import RegularGridInterpolator as RGI
F = '/orion/ptmp/jinma/clean_1115/stag2/freshgo/fresh_prod_1120/'
SRC = F + 't9l62/corona_test/ic_t9l62_corona.txt'
OUT = sys.argv[1] if len(sys.argv) > 1 else 'ic_t9l62_corona_iso.txt'
TAB = F + 'tables/rosseland_mesakap_gs98_x0.70_z0.02.txt'
mu = 0.674; mu_u = 1.66053906660e-24; kB = 1.380649e-16; c = 2.99792458e10; gam = 5/3
GM = 7.4033817458e27; LIC = 2.47875e39; RJ = 6.641422e12; R_TOP = 5.632159e13
A = open(TAB).read().splitlines(); gi = [i for i, l in enumerate(A) if l.startswith('# grid')][0]
h = A[gi+1].lstrip('#').split(); nT, nD = int(h[0]), int(h[1]); lT0, dlT, lD0, dlD = map(float, h[2:6])
v = np.array([float(x) for x in [l for l in A if l.strip() and not l.startswith('#')][:nT*nD]]).reshape(nT, nD)
fk = RGI((lT0+dlT*np.arange(nT), lD0+dlD*np.arange(nD)), v)
kap = lambda T, rho: 10**fk([[np.clip(np.log10(T), lT0, lT0+dlT*(nT-1)), np.clip(np.log10(rho), lD0, lD0+dlD*(nD-1))]])[0]
lines = open(SRC).read().splitlines(); hdr = [l for l in lines if l.startswith('#')]
d = np.loadtxt(SRC); r, rho, eint, Fr, E, T, fm = d.T
j = np.nonzero(r >= RJ*(1 - 1e-9))[0][0]
Tj = T[j]; cs2 = kB*Tj/(mu*mu_u)
nr = rho.copy(); lr = np.log(rho[j])
for n in range(j+1, r.size):          # midpoint (RK2) in r, isothermal
    x0, x1 = r[n-1], r[n]; hh = x1 - x0
    f = lambda x, l: -GM*(1 - kap(Tj, np.exp(l))*LIC/(4*np.pi*c*GM))/(x*x*cs2)
    k1 = f(x0, lr); k2 = f(x0 + hh/2, lr + hh/2*k1); lr = lr + hh*k2; nr[n] = np.exp(lr)
nT_ = np.where(np.arange(r.size) >= j, Tj, T)
ne = nr*kB*nT_/(mu*mu_u)/(gam - 1)
print('r_j %.4e T_j %.3e rho_j %.3e rho(x1max) %.3e rho(rmax) %.3e' % (r[j], Tj, nr[j], np.interp(R_TOP, r, nr), nr[-1]))
hdr = hdr[:-1] + ['# CORONA9 11-21 build_ic_corona_iso.py: the corona above r_j made ISOTHERMAL at T_j = %.4e K (HSE with kappa_R Gamma), rho_j unchanged' % Tj, hdr[-1]]
np.savetxt(OUT, np.c_[r, nr, ne, Fr, E, nT_, fm], fmt='%.12e', header='\n'.join(h_.lstrip('# ') if False else h_[2:] for h_ in hdr), comments='# ')
print('wrote', OUT)
