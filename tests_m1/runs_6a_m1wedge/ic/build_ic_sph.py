#!/usr/bin/env python3
"""
Spherical initial condition for the m1-wedge He-envelope test (m1_test = sph_wedge,
wg_ic = file): the He-box FeCZ star of bench/m1_stage2/ic (V3edd box: g0 = 3.98107e5,
F = 2.475202e15 = sigma (81283 K)^4) put on a sphere of radius R0 = sqrt(GM/g0) with
M = <--mstar> (3.15 Msun default), and re-solved as a SPHERICAL grey radiative- and
hydrostatic-equilibrium column (Eddington closure, gas-only EOS, Rosseland table):

    F(r)    = F0 (R0/r)^2
    dE/dr   = -3 rho kappa_R(rho,T) F/c,       E_top = F_top/(c q), q = 1/2
    dpg/dr  = -rho (GM/r^2 - kappa_R F/c),     T = (E/a)^(1/4),  rho = rho(pg, T)

integrated inward (and a few cells outward, for the ghosts) from the box top face
z = 4.1290910e7 cm with the V3edd top density.  EOS, tables and constants are imported
from bench/m1_stage2/ic/build_ic.py (read-only).  Output: r rho eint E (cgs = code units
with <units> 1 1 1), ascending r, on a fine grid of nsub nodes per box cell.

usage: build_ic_sph.py [--mstar 3.15] [--zbot -1.8149339e8] [--out ic_wedge_he.txt]
"""
import argparse, sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/m1_stage2/ic')
import build_ic as B   # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--mstar', type=float, default=3.15)
ap.add_argument('--zbot', type=float, default=-1.8149339e8)
ap.add_argument('--ztop', type=float, default=4.1290910e7)
ap.add_argument('--nx1', type=int, default=84)
ap.add_argument('--nsub', type=int, default=16)
ap.add_argument('--pad', type=int, default=8, help='cells beyond each face')
ap.add_argument('--out', default='ic_wedge_he.txt')
ap.add_argument('--v3edd', default='/viper/ptmp2/jinma/caltech_handover_0926/athenak_data/'
                'he_box/ic_m1_V3edd_pgen.txt')
a = ap.parse_args()

G = 6.67430e-8
MSUN = 1.98847e33
GM = G*a.mstar*MSUN
R0 = np.sqrt(GM/B.G0)
F0 = B.FLUX
c, ar = B.C_LIGHT, B.A_RAD
tabR = B.read_opac(B.OPAC_R)
dz = (a.ztop - a.zbot)/a.nx1
rtop = R0 + a.ztop
zc, rc, ec = np.loadtxt(a.v3edd, unpack=True)
rho_top = float(np.exp(np.interp(a.ztop, zc, np.log(rc))))


def kap(rho, t):
    return float(B.kappa_of(tabR, np.array([t]), np.array([rho]))[0])


def state(r, E, pg, rg):
    t = (E/ar)**0.25
    rho = float(B.rho_from_p_t(np.array([pg]), np.array([t]), np.array([rg]))[0])
    return t, rho


def rhs(r, E, pg, rg):
    t, rho = state(r, E, pg, rg)
    F = F0*(R0/r)**2
    k = kap(rho, t)
    return -3.0*rho*k*F/c, -rho*(GM/r**2 - k*F/c), rho, t


def march(r0, E, pg, rg, h, n):
    out = []
    r = r0
    for _ in range(n):
        k1 = rhs(r, E, pg, rg)
        rg = k1[2]
        out.append((r, k1[2], k1[3], E, pg))
        E1, p1 = E + h*k1[0], pg + h*k1[1]
        k2 = rhs(r + h, E1, p1, rg)
        E, pg = E + 0.5*h*(k1[0] + k2[0]), pg + 0.5*h*(k1[1] + k2[1])
        r += h
    k1 = rhs(r, E, pg, rg)
    out.append((r, k1[2], k1[3], E, pg))
    return out


Etop = F0*(R0/rtop)**2/(c*0.5)
ttop = (Etop/ar)**0.25
pgtop = float(B.p_gas_of(np.array([rho_top]), np.array([ttop]))[0])
h = dz/a.nsub
nin = (a.nx1 + a.pad)*a.nsub
nout = a.pad*a.nsub
dn = march(rtop, Etop, pgtop, rho_top, -h, nin)
up = march(rtop, Etop, pgtop, rho_top, +h, nout)
rows = sorted(dn + up[1:], key=lambda x: x[0])
r = np.array([x[0] for x in rows]); rho = np.array([x[1] for x in rows])
t = np.array([x[2] for x in rows]); E = np.array([x[3] for x in rows])
eg = B.e_gas_of(rho, t)
kR = B.kappa_of(tabR, t, rho)
F = F0*(R0/r)**2
gam = kR*F/(c*GM/r**2)
tau = np.zeros_like(r)
m = r <= rtop
x = r[m][::-1]; y = (rho*kR)[m][::-1]
tau[m] = np.concatenate([[0.0], np.cumsum(0.5*(y[1:] + y[:-1])*(x[:-1] - x[1:]))])[::-1]
rbot = R0 + a.zbot
with open(a.out, 'w') as f:
    f.write('# m1-wedge He IC (build_ic_sph.py): M = %.4f Msun, GM = %.6e, R0 = %.6e cm, '
            'F0 = %.6e, rho_top = %.6e\n' % (a.mstar, GM, R0, F0, rho_top))
    f.write('# mesh: x1min = %.9e  x1max = %.9e  nx1 = %d  (box z %.6e .. %.6e)\n'
            % (rbot, rtop, a.nx1, a.zbot, a.ztop))
    f.write('# implicit_flux_x1min = %.9e   wg_gm = %.9e\n' % (F0*(R0/rbot)**2, GM))
    f.write('# columns: r[cm] rho[g/cm^3] eint_gas[erg/cm^3] E_rad[erg/cm^3]\n')
    for i in range(len(r)):
        f.write('%.12e %.12e %.12e %.12e\n' % (r[i], rho[i], eg[i], E[i]))
print('R0 = %.6e cm, GM = %.6e, x1min = %.9e, x1max = %.9e, F_in = %.9e'
      % (R0, GM, rbot, rtop, F0*(R0/rbot)**2))
for zz in [a.ztop, 0.0, -5e7, -1e8, -1.5e8, a.zbot]:
    i = np.argmin(abs(r - (R0 + zz)))
    print('z %+.3e rho %.3e T %.4e tau %.3e Gamma %.3f Prad/Pg %.3f' %
          (zz, rho[i], t[i], tau[i], gam[i], E[i]/3/B.p_gas_of(rho[i:i+1], t[i:i+1])[0]))
print('max Gamma %.3f at z = %.4e' % (gam[m].max(), r[m][np.argmax(gam[m])] - R0))
