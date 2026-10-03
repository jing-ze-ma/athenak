# flake8: noqa
"""Sector radial grid: the 522-grid faces up to 1.1 R*, then cubic cells dr = r dtheta
(dtheta = pi/512, geometric r_{k+1} = r_k exp(dtheta)) to x1max; fit StretchRPoly (18 p)."""
import sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '/viper/ptmp2/jinma/he_mltpp_1002/remap128')
import he_remap_rst as H
RS = 2.3717e11
X0, X1 = 9.4868e10, 7.1151e11
p522 = np.load('/viper/ptmp2/jinma/he_mltpp_1002/remap128/out/grid522.npy')
p522 = np.array([float('%.12e' % v) for v in p522])
f522 = H.faces_of(p522, 522, 0, X0, X1)
DTH = np.pi/512
ks = np.searchsorted(f522, 1.1*RS)            # first 522 face >= 1.1 R*
rin = f522[:ks+1]
nout = int(round(np.log(X1/rin[-1])/DTH))
q = np.log(X1/rin[-1])/nout
rout = rin[-1]*np.exp(q*np.arange(1, nout+1))
ft = np.concatenate([rin, rout])
N = ft.size - 1
print('inner cells %d (to %.4f R*), outer %d (dr/r %.5f vs dtheta %.5f), N %d' % (
    ks, rin[-1]/RS, nout, np.expm1(q), DTH, N))
drt = np.diff(ft)
import os
PW = float(os.environ.get('PW', '20'))
if os.environ.get('P0'):
    p522 = np.load(os.environ['P0'])


def res(p):
    f = H.faces_of(p, N, 0, X0, X1)
    dr = np.diff(f)
    if np.any(dr <= 0):
        return np.full(3*N+1, 1e3)
    lr = np.log(dr/drt)
    pen = PW*np.maximum(0.0, lr - float(os.environ.get("THR", "0.02")))*(0.5*(ft[1:] + ft[:-1]) < 1.1*RS)
    return np.concatenate([lr, 0.3*(f - ft)/(0.003*RS), pen])


best = None
p0 = p522.copy()
for it in range(8):
    lb = np.full(18, -np.inf)
    ub = np.full(18, np.inf)
    for b in range(2):
        lb[8+3*b], ub[8+3*b] = -3.0, 3.0
        lb[9+3*b], ub[9+3*b] = 0.0, 1.0
        lb[10+3*b], ub[10+3*b] = 0.005, 0.5
    lb[14], ub[14] = -0.99, 0.0
    lb[15], ub[15] = -0.2, 1.0
    lb[16], ub[16] = 0.0, 1.0
    lb[17], ub[17] = 0.003, 0.2
    p0 = np.clip(p0, lb + 1e-9, ub - 1e-9)
    s = least_squares(res, p0, method='trf', max_nfev=6000, x_scale='jac', bounds=(lb, ub))
    if best is None or s.cost < best.cost:
        best = s
    p0 = best.x*(1 + 0.02*np.random.default_rng(it).standard_normal(18))
p = np.array([float('%.12e' % v) for v in best.x])
f = H.faces_of(p, N, 0, X0, X1)
dr = np.diff(f)
rc = 0.5*(f[1:]+f[:-1])
e = dr/np.interp(rc, 0.5*(ft[1:]+ft[:-1]), drt) - 1
print('cost %.3g; dr/dr_target - 1: max %.3f min %.3f; max adj ratio %.3f' % (
    best.cost, e.max(), e.min(), np.max(np.maximum(dr[1:]/dr[:-1], dr[:-1]/dr[1:]))))
for x in (0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0, 1.05, 1.1, 1.15, 1.2,
          1.5, 2.0, 2.5, 2.99):
    i = np.argmin(abs(rc - x*RS))
    j = np.argmin(abs(0.5*(f522[1:]+f522[:-1]) - x*RS))
    print('  r/R* %.2f dr/R* %.5f (522: %.5f, target %.5f) dr/(r dth) %.2f' % (
        x, dr[i]/RS, np.diff(f522)[j]/RS, np.interp(rc[i], 0.5*(ft[1:]+ft[:-1]), drt)/RS,
        dr[i]/(rc[i]*DTH)))
np.save('/viper/ptmp2/jinma/he_mltpp_1002/sector/grid/p_%d.npy' % N, p)
print('p =', ' '.join('%.12e' % v for v in p))
