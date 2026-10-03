# flake8: noqa
"""dt of the hydro CFL with the rad signal speed (hydro_newdt.cpp formula), from the restart's
own caches (Gamma_1, p, opac, tau_ten chi/n): old grid at 64^2, and with replaced dr / angular
resolution.  usage: python3 dtest.py [newgrid_faces.npy]"""
import he_remap_rst as H
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/he_mltpp_1002/remap128')
r = H.read_rst('/viper/ptmp2/jinma/he_mltpp_1002/remap128/src/hepresn.00033.rst')
G = H.gather(r)
ng = r.ng
sl = slice(ng, ng + r.bind[1])
rho = G['hyd'][0][..., sl]
v = [G['hyd'][1 + n][..., sl] / rho for n in range(3)]
p = G['wd'][0][..., sl]
g1 = G['wd'][1][..., sl]
E = G['m1'][0][..., sl]
opac = G['rss'][0][..., sl]
chi = G['rss'][1][..., sl]
nv = [G['rss'][2 + n][..., sl] for n in range(3)]
cg2 = g1 * p / rho
rows = [l.split() for l in open('/viper/ptmp2/jinma/he_mltpp_1002/race/pp2/hepresn.x1grid.txt')
        if l.strip() and l[0] != '#' and len(l.split()) >= 4]
xf = np.array([float(q[1]) for q in rows])
dx = np.array([float(q[3]) for q in rows])
dr = dx[3:473]
rv = np.array([float(q[2]) for q in rows])[3:473]
R = 9.4868e10 / 0.4
dth = r.msize[7]
dph = r.msize[8]
N2 = r.mind[2]
th = r.msize[1] + (np.arange(N2) + 0.5) * dth


def ddt(dxd, d, vd):
    pdd = E * ((1 - chi) / 2 + (3 * chi - 1) / 2 * nv[d]**2)
    ceq2 = cg2 + 4.0 / 3.0 * pdd / rho
    tau = opac * dxd
    c2 = cg2 + (1 - np.exp(-tau)) * (ceq2 - cg2)
    return dxd / (np.abs(vd) + np.sqrt(c2))


def report(dr_, f, tag):
    d1 = ddt(dr_[None, None, :], 0, v[0])
    d2 = ddt(rv[None, None, :] * dth / f, 1, v[1])
    d3 = ddt(rv[None, None, :] * np.sin(th)[None, :, None] * dph / f, 2, v[2])
    out = []
    for nm, a in (('r', d1), ('th', d2), ('ph', d3)):
        i = np.unravel_index(np.argmin(a), a.shape)
        out.append('%s %.4g s at r/R %.4f' % (nm, 0.3 * a.min(), rv[i[2]] / R))
    print('%-28s cfl*min dx/(|v|+c): ' % tag + '; '.join(out))


print('restart dt %.5g s (cfl 0.3)' % r.dt)
report(dr, 1, 'old 470, 64^2')
report(dr, 2, 'old 470, 128^2')
if len(sys.argv) > 1:
    fn = np.load(sys.argv[1])
    drn = np.interp(rv, 0.5 * (fn[1:] + fn[:-1]), np.diff(fn))
    report(drn, 1, 'new grid, 64^2')
    report(drn, 2, 'new grid, 128^2')
