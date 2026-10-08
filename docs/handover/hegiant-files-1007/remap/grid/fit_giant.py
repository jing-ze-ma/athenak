#!/usr/bin/env python3
"""Radial grid for the He giant: fit the code's StretchRPoly (8 poly coefs + 2 sech^2 bumps +
tanh plateau, mesh/f_stretch_r_*) to a target dr(r) on x1min = 3 Rsun .. x1max = RTOP.
Target (env DLN, DLNX, RTOP, NPH):
  * r < 57 Rsun: dr = max(DLN r, DMIN) (log grid; DLN 0.015, DMIN 0.10 Rsun);
  * plateau 60.3 .. 64.2 Rsun: dr = H/NPH, H = 0.117 Rsun = the IC atmosphere's density
    scale height above R_ph (61.78 Rsun) (NPH default 3);
  * a smooth geometric (5 %/cell) transition into and out of the plateau;
  * r > 64.2 Rsun: grows 5 %/cell up to dr = DLNX r (default 0.01) to 130 Rsun, then
    dr/r rises linearly to 0.025 at RTOP.
usage: fit_giant.py [N]   -> prints the keys and writes p_N.npy, faces_N.txt"""
import os
import sys
import numpy as np
from scipy.optimize import least_squares

RS = 6.957e10
DLN = float(os.environ.get('DLN', '0.015'))
DLNX = float(os.environ.get('DLNX', '0.01'))
RTOP = float(os.environ.get('RTOP', '200'))*RS
NPH = float(os.environ.get('NPH', '3'))
R0, R1 = 3.0*RS, RTOP
H_ATM = 0.117*RS
rr = np.linspace(R0, R1, 400001)
x = rr/RS
drp = H_ATM/NPH
GR = 1.05
d_in = np.where(x < 57, DLN*rr, np.inf)
# into the plateau from below: geometric approach (dr shrinks 5 %/cell towards 60.3)
dist_lo = np.maximum(60.3*RS - rr, 0.0)
d_lo = drp + (GR - 1.0)*dist_lo
d_hi_side = drp + (GR - 1.0)*np.maximum(rr - 64.2*RS, 0.0)
outer = np.where(x < 130, DLNX*rr, rr*(DLNX + (0.025 - DLNX)*(x - 130)/(RTOP/RS - 130)))
tgt = np.where(x < 60.3, np.minimum(np.maximum(d_in, 0) if False else np.minimum(d_in, d_lo),
                                    d_lo), drp)
DMIN = float(os.environ.get('DMIN', '0.10'))*RS     # deep floor on dr (r_in H_p 1.0 Rsun)
tgt = np.where(x < 57, np.minimum(np.maximum(DLN*rr, DMIN), d_lo), tgt)
tgt = np.where((x >= 57) & (x < 60.3), d_lo, tgt)
tgt = np.where(x > 64.2, np.minimum(d_hi_side, outer), tgt)
N_t = np.trapz(1.0/tgt, rr)


def lch(z):
    a = np.abs(z)
    return a + np.log1p(np.exp(-2*a)) - np.log(2.0)


def G(xa, xb, w, xi):
    return 0.5*w*(lch((xi-xa)/w) - lch(-xa/w) - lch((xi-xb)/w) + lch(-xb/w))


def ufun(p, xi):
    u = xi.copy()
    xik = xi.copy()
    for k in range(8):
        u = u + p[k]*xik*(1-xi)
        xik = xik*xi
    for b in range(2):
        a, xb, w = p[8+3*b:11+3*b]
        if a != 0:
            u = u + a*w*(np.tanh((xi-xb)/w) - (1-xi)*np.tanh(-xb/w) - xi*np.tanh((1-xb)/w))
    a, xa, xb, w = p[14:18]
    if a != 0:
        u = u + a*(G(xa, xb, w, xi) - xi*G(xa, xb, w, 1.0))
    return u


def faces(p, N):
    return R0 + (R1-R0)*ufun(p, np.linspace(0, 1, N+1))


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else int(round(N_t))
    print('target cells %.1f, plateau dr %.4f Rsun; fitting N = %d' % (N_t, drp/RS, N))
    cum = np.concatenate([[0], np.cumsum(0.5*(1/tgt[1:]+1/tgt[:-1])*np.diff(rr))])
    xi_t = np.linspace(0, 1, N+1)
    rf_t = np.interp(xi_t*cum[-1], cum, rr)

    def res(p):
        rf = faces(p, N)
        dr = np.diff(rf)
        if np.any(dr <= 0):
            return np.full(3*N + 1, 1e3)
        rcn = 0.5*(rf[1:]+rf[:-1])
        tg = np.interp(rcn, rr, tgt)
        r1 = np.log(dr/tg)
        r2 = 10*np.maximum(0, np.log(dr/(1.15*tg)))*(rcn < 70*RS)
        r3 = 0.3*(rf - rf_t)/(0.003*(R1 - R0))
        return np.concatenate([r1, r2, r3])
    # initial guess: least-squares polynomial to the target u(xi), then the plateau
    p0 = np.zeros(18)
    xa0 = np.interp(60.3*RS, rf_t, xi_t)
    xb0 = np.interp(64.2*RS, rf_t, xi_t)
    # smooth part: the target WITHOUT the photospheric valley (log + outer law), as a
    # polynomial in xi; the valley then comes from the plateau (flat part) and bump 1
    tgs = np.where(x < 64.2, np.maximum(DLN*rr, DMIN), outer)
    cs = np.concatenate([[0], np.cumsum(0.5*(1/tgs[1:]+1/tgs[:-1])*np.diff(rr))])
    xs = np.interp(rr, rr, cs)/cs[-1]
    A = np.array([xs**(k+1)*(1-xs) for k in range(8)]).T
    p0[:8] = np.linalg.lstsq(A, (rr - R0)/(R1 - R0) - xs, rcond=None)[0]
    p0[14:18] = [float(os.environ.get('PA', '-0.9')), xa0, xb0, 0.004]
    p0[8:11] = [float(os.environ.get('BA', '-0.2')), 0.5*(xa0 + xb0), 0.5*(xb0 - xa0) + 0.05]
    p0[11:14] = [0.0, 0.9, 0.05]
    lb = np.full(18, -np.inf)
    ub = np.full(18, np.inf)
    for b in range(2):
        lb[8+3*b], ub[8+3*b] = -3.0, 3.0
        lb[9+3*b], ub[9+3*b] = 0.0, 1.0
        lb[10+3*b], ub[10+3*b] = 0.003, 0.5
    lb[14], ub[14] = -0.999, 0.0
    lb[15:17], ub[15:17] = 0.0, 1.0
    lb[17], ub[17] = 0.001, 0.2
    best = None
    for it in range(8):
        p0 = np.clip(p0, lb + 1e-9, ub - 1e-9)
        sol = least_squares(res, p0, method='trf', max_nfev=6000, x_scale='jac',
                            bounds=(lb, ub))
        if best is None or sol.cost < best.cost:
            best = sol
        p0 = best.x*(1 + 0.03*np.random.default_rng(it).standard_normal(18))
    p = best.x
    rf = faces(p, N)
    dr = np.diff(rf)
    rcn = 0.5*(rf[1:]+rf[:-1])
    tg = np.interp(rcn, rr, tgt)
    print('N %d cost %.4g; dr/target min %.3f max %.3f; max neighbour ratio %.3f; dr min '
          '%.4f Rsun' % (N, best.cost, (dr/tg).min(), (dr/tg).max(),
                         np.max(np.maximum(dr[1:]/dr[:-1], dr[:-1]/dr[1:])), dr.min()/RS))
    for xx in (3, 5, 10, 20, 30, 40, 50, 55, 58, 60, 61, 61.5, 61.8, 62.5, 63.5, 64.5, 66,
               70, 80, 100, 130, 160, 199):
        i = np.argmin(abs(rcn - xx*RS))
        print('  r %6.2f Rsun: dr %.4f Rsun dr/r %.5f target %.4f' % (
            rcn[i]/RS, dr[i]/RS, dr[i]/rcn[i], tg[i]/RS))
    np.save('p_%d.npy' % N, p)
    np.savetxt('faces_%d.txt' % N, rf)
    names = ['f_stretch_r_c%d' % (k+1) for k in range(8)]
    for b in range(2):
        names += ['f_stretch_r_b%d_amp' % (b+1), 'f_stretch_r_b%d_x' % (b+1),
                  'f_stretch_r_b%d_w' % (b+1)]
    names += ['f_stretch_r_p_amp', 'f_stretch_r_p_xa', 'f_stretch_r_p_xb',
              'f_stretch_r_p_w']
    for n_, v in zip(names, p):
        print('%s = %.12e' % (n_, v))
