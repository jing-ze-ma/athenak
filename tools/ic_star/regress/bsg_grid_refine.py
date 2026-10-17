#!/usr/bin/env python3
"""bsgres-1017 grid_refine.py FAC [NX]: grid_2xatm.py generalized: the g480 fine plateau dr multiplied by FAC
(FAC 0.5 = the a776 2x-atmosphere grid; FAC 0.0320/0.03363 = the refined deck grid r496-class, item 1 of gate (a')).
Original docstring: '2x atmosphere' radial grid = the g480 grid (grid_nx480.txt) with dr halved in its fine
photosphere/atmosphere plateau (dr 0.0336 -> 0.0168 Rsun over ~50-58.3 Rsun) and unchanged elsewhere.
Target dr_new(r) = dr480(r) * f(r), f = 1 - 0.5 clip((0.2 - dr480)/(0.2 - 0.0336), 0, 1) (smooth in the transitions);
nx from int dr/dr_new; then the 8 grid_stretch.hpp two-plateau parameters are fitted (least squares on the face radii).
usage: grid_2xatm.py   (prints the mesh keys + checks)"""
import sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '/raven/ptmp/jinma/bsg_newsetup_1015/grid')
R0, R1 = 20.0, 80.0


def lch(x):
    ax = np.abs(x)
    return ax + np.log1p(np.exp(-2*ax)) - np.log(2.0)


def G(xa, xb, w, xi):
    return 0.5*w*(lch((xi - xa)/w) - lch(-xa/w) - lch((xi - xb)/w) + lch(-xb/w))


def grid(nx, a1, xa1, xb1, w1, a2, xa2, xb2, w2):
    xi = np.linspace(0, 1, nx + 1)
    u = xi + a1*(G(xa1, xb1, w1, xi) - xi*G(xa1, xb1, w1, 1.0)) + a2*(G(xa2, xb2, w2, xi) - xi*G(xa2, xb2, w2, 1.0))
    return R0 + (R1 - R0)*u


P480 = (-1.270035907858e+00, 9.258713492978e-02, 8.768396646040e-01, 1.5e-02,
        -1.233181510246e+00, 2.478286835644e-01, 8.768396646040e-01, 1.5e-02)
rf480 = grid(480, *P480)
dr480 = np.diff(rf480); rc480 = 0.5*(rf480[1:] + rf480[:-1])
DF = dr480.min()
FAC = float(sys.argv[1])
f = 1 - (1 - FAC)*np.clip((0.2 - dr480)/(0.2 - DF), 0, 1)
drt = dr480*f
# target face radii: integrate cells of size drt(r) on a fine r grid
rr = np.linspace(R0, R1, 600001)
dri = np.interp(rr, rc480, drt)
N = np.concatenate([[0], np.cumsum(np.diff(rr)/(0.5*(dri[1:] + dri[:-1])))])
nx = int(sys.argv[2]) if len(sys.argv) > 2 else int(8*np.ceil(N[-1]/8))
rt = np.interp(np.linspace(0, N[-1], nx + 1), N, rr)
s = nx/480.0


def res(q):
    a1, xa1, a2, xa2, xb, w = q
    try:
        r = grid(nx, a1, xa1, xb, w, a2, xa2, xb, w)
    except FloatingPointError:
        return 1e3*np.ones(nx + 1)
    return (r - rt)/np.interp(r, rr, dri)


# initial: same edge radii as g480 -> index fractions on the new target
x0 = []
for xa in (P480[1], P480[5], P480[2]):
    redge = np.interp(xa, np.linspace(0, 1, 481), rf480)
    x0.append(np.interp(redge, rt, np.linspace(0, 1, nx + 1)))
lev0 = drt.max()*nx/60; lev1 = 0.1875*nx/60; lev2 = DF*FAC*nx/60
q0 = [lev1 - lev0, x0[0], lev2 - lev1, x0[1], x0[2], 0.015*480/nx]
L1, L2 = x0[2] - x0[0], x0[2] - x0[1]
sc = 1.0/(1.0 + (lev1/lev0 - 1)*L1 + (lev2/lev0 - lev1/lev0)*L2)
q0[0] = (lev1/lev0 - 1)*sc; q0[2] = (lev2/lev0 - lev1/lev0)*sc
sol = least_squares(res, q0, x_scale='jac', max_nfev=4000)
a1, xa1, a2, xa2, xb, w = sol.x
p = (a1, xa1, xb, w, a2, xa2, xb, w)
rf = grid(nx, *p); dr = np.diff(rf); rc = 0.5*(rf[1:] + rf[:-1])
assert np.all(dr > 0), 'fold'
print('nx1 %d (target cells %.1f); fit rms %.3f cells, max %.3f cells' % (nx, N[-1], np.sqrt(np.mean(sol.fun**2)),
                                                                           np.abs(sol.fun).max()))
print('mesh/nx1=%d meshblock/nx1=%d mesh/f_stretch_r_p_amp=%.12e mesh/f_stretch_r_p_xa=%.12e mesh/f_stretch_r_p_xb=%.12e '
      'mesh/f_stretch_r_p_w=%.12e mesh/f_stretch_r_p2_amp=%.12e mesh/f_stretch_r_p2_xa=%.12e mesh/f_stretch_r_p2_xb=%.12e '
      'mesh/f_stretch_r_p2_w=%.12e' % ((nx, nx) + p))
for x in [20, 25, 30, 34, 36, 40, 45, 48, 49, 50, 50.4, 51, 52, 53, 54, 55, 56, 57, 58, 58.3, 59, 60, 62, 65, 70, 80]:
    i = np.argmin(abs(rc - x)); j = np.argmin(abs(rc480 - x))
    print('  r %6.2f dr %.4f   (g480 %.4f, ratio %.3f)' % (rc[i], dr[i], dr480[j], dr[i]/dr480[j]))
band = (rc >= 50.38) & (rc <= 58.33)
print('0.95-1.1 R_ph (50.38-58.33 Rsun): %d cells (g480 %d), max dr %.4f -> cells per 0.1424 Rsun min %.2f; '
      'max neighbour ratio %.3f; min dr %.4f' % (band.sum(), ((rc480 >= 50.38) & (rc480 <= 58.33)).sum(), dr[band].max(),
                                                 0.1424/dr[band].max(),
                                                 np.max(np.maximum(dr[1:]/dr[:-1], dr[:-1]/dr[1:])), dr.min()))
