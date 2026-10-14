#!/usr/bin/env python3
"""fit_giant3.py (hegiant_large_1014, worker HGL): radial grid for the LARGE-domain He giant.

Same target and same StretchRPoly family as fit_giant2.py (Viper he_giant_1006/grid, -> fit2_897.out,
the current 200-Rsun grid), extended outward:
  * r < 57 Rsun: dr = max(DLN r, DMIN) (DLN 0.015, DMIN 0.10 Rsun), geometric (5 %/cell) approach to
    the fine band;                                                              [= fit_giant2]
  * fine band PA (60) .. RB: dr = H_atm (r/61.78)^2 / NPH (H_atm 0.117 Rsun, NPH 3) = the IC photospheric
    scale height scaled as r^2 at fixed T                                       [= fit_giant2 on 60-95]
    -- continued from 95 Rsun up to RB (default 0.75 RTOP = the rmax-guard radius), where the
    photosphere can move during the run;
  * r > RB: dr grows 5 %/cell, capped at DLNO r (0.025), to RTOP.
Fit: 8 polynomial coefs + 2 sech^2 bumps + NPLAT (1|2) tanh plateaus (mesh/f_stretch_r_*,
f_stretch_r_p2_*), least squares in log(dr/target) + a 15 % overshoot penalty inside RB + face
positions; staged (positions first, then widths) from a monotone start.
usage: RTOP=700 [RB=525 NPLAT=1] fit_giant3.py [N]  -> keys on stdout, p3_N.npy, faces3_N.txt"""
import os
import sys
import numpy as np
from scipy.optimize import least_squares

RS = 6.957e10
PA = float(os.environ.get('PA', '60.0'))
DLN = float(os.environ.get('DLN', '0.015'))
RTOP = float(os.environ.get('RTOP', '700'))*RS
RB = float(os.environ.get('RB', '%g' % (0.75*RTOP/RS)))
DLNO = float(os.environ.get('DLNO', '0.025'))
NPH = float(os.environ.get('NPH', '3'))
NPLAT = int(os.environ.get('NPLAT', '1'))
NP = 14 + 4*NPLAT
R0, R1 = 3.0*RS, RTOP
H_ATM = 0.117*RS
rr = np.linspace(R0, R1, 800001)
x = rr/RS
drp = H_ATM/NPH
GR = 1.05
DMIN = float(os.environ.get('DMIN', '0.10'))*RS
d_lo = drp + (GR - 1.0)*np.maximum(PA*RS - rr, 0.0)
band = drp*(np.maximum(rr, PA*RS)/(61.78*RS))**2
d_hi = drp*(RB/61.78)**2 + (GR - 1.0)*np.maximum(rr - RB*RS, 0.0)
tgt = np.where(x < 57, np.minimum(np.maximum(DLN*rr, DMIN), d_lo), band)
tgt = np.where((x >= 57) & (x < PA), d_lo, tgt)
tgt = np.where(x > RB, np.minimum(d_hi, DLNO*rr), tgt)
N_t = np.trapezoid(1.0/tgt, rr) if hasattr(np, 'trapezoid') else np.trapz(1.0/tgt, rr)


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
    for q in range(NPLAT):
        a, xa, xb, w = p[14+4*q:18+4*q]
        if a != 0:
            u = u + a*(G(xa, xb, w, xi) - xi*G(xa, xb, w, 1.0))
    return u


def faces(p, N):
    return R0 + (R1-R0)*ufun(p, np.linspace(0, 1, N+1))


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else int(round(N_t))
    print('RTOP %.1f RB %.1f: target cells %.1f, plateau dr %.4f Rsun; fitting N = %d'
          % (RTOP/RS, RB, N_t, drp/RS, N))
    cum = np.concatenate([[0], np.cumsum(0.5*(1/tgt[1:]+1/tgt[:-1])*np.diff(rr))])
    xi_t = np.linspace(0, 1, N+1)
    rf_t = np.interp(xi_t*cum[-1], cum, rr)
    scale = 0.003*(R1 - R0)

    def res(p, wdr=1.0):
        rf = faces(p, N)
        dr = np.diff(rf)
        rcn = 0.5*(rf[1:]+rf[:-1])
        tg = np.interp(rcn, rr, tgt)
        neg = np.maximum(0.0, -dr)/tg
        drs = np.maximum(dr, 1e-6*tg)
        r1 = wdr*np.log(drs/tg) + 1e3*neg
        r2 = wdr*10*np.maximum(0, np.log(drs/(1.15*tg)))*(rcn < (RB + 5)*RS)
        r3 = 0.3*(rf - rf_t)/scale
        return np.concatenate([r1, r2, r3])
    xa0 = np.interp(57.0*RS, rf_t, xi_t)
    # monotone start: polynomial lstsq to the target u(xi) with a plateau over [xa0, 1)
    p0 = np.zeros(NP)
    p0[14:18] = [-0.8, xa0, 0.999, 0.004]
    if NPLAT == 2:
        p0[18:22] = [-0.3, np.interp(RB*RS, rf_t, xi_t), 0.999, 0.004]
    base = ufun(p0, xi_t)
    A = np.array([xi_t**(k+1)*(1-xi_t) for k in range(8)]).T
    p0[:8] = np.linalg.lstsq(A, (rf_t - R0)/(R1 - R0) - base, rcond=None)[0]
    lb = np.full(NP, -np.inf)
    ub = np.full(NP, np.inf)
    for b in range(2):
        lb[8+3*b], ub[8+3*b] = -3.0, 3.0
        lb[9+3*b], ub[9+3*b] = 0.0, 1.0
        lb[10+3*b], ub[10+3*b] = 0.003, 0.5
    for q in range(NPLAT):
        lb[14+4*q], ub[14+4*q] = -0.999, 0.0
        lb[15+4*q:17+4*q], ub[15+4*q:17+4*q] = 0.0, 1.0
        lb[17+4*q], ub[17+4*q] = 0.001, 0.2
    p0[8:14] = [-0.05, xa0, 0.05, -0.05, 0.98, 0.02]
    p0 = np.clip(p0, lb + 1e-9, ub - 1e-9)
    best = None
    rng = np.random.default_rng(0)
    for it in range(int(os.environ.get('NIT', '10'))):
        s1 = least_squares(lambda p: res(p, 0.05), p0, method='trf', max_nfev=3000,
                           x_scale='jac', bounds=(lb, ub))
        sol = least_squares(res, s1.x, method='trf', max_nfev=8000, x_scale='jac',
                            bounds=(lb, ub))
        print('  it %d cost %.4g' % (it, sol.cost), flush=True)
        if best is None or sol.cost < best.cost:
            best = sol
        p0 = np.clip(best.x*(1 + 0.03*rng.standard_normal(NP)), lb + 1e-9, ub - 1e-9)
    p = best.x
    rf = faces(p, N)
    dr = np.diff(rf)
    rcn = 0.5*(rf[1:]+rf[:-1])
    tg = np.interp(rcn, rr, tgt)
    print('N %d cost %.4g; dr/target min %.3f max %.3f; max neighbour ratio %.3f; dr min '
          '%.4f Rsun; monotone %s' % (N, best.cost, (dr/tg).min(), (dr/tg).max(),
                                      np.max(np.maximum(dr[1:]/dr[:-1], dr[:-1]/dr[1:])),
                                      dr.min()/RS, bool(np.all(dr > 0))))
    for xx in (3, 5, 10, 20, 30, 40, 50, 55, 58, 60, 62, 66, 70, 75, 80, 85, 90, 95, 100,
               110, 130, 160, 200, 250, 300, 350, 400, 450, 500, 525, 550, 600, 650, 700,
               800, 900, 1000):
        if xx > RTOP/RS:
            continue
        i = np.argmin(abs(rcn - xx*RS))
        print('  r %7.2f Rsun: dr %.4f Rsun dr/r %.5f target %.4f  cells below %d' % (
            rcn[i]/RS, dr[i]/RS, dr[i]/rcn[i], tg[i]/RS, i))
    np.save('p3_%d.npy' % N, p)
    np.savetxt('faces3_%d.txt' % N, rf)
    names = ['f_stretch_r_c%d' % (k+1) for k in range(8)]
    for b in range(2):
        names += ['f_stretch_r_b%d_amp' % (b+1), 'f_stretch_r_b%d_x' % (b+1),
                  'f_stretch_r_b%d_w' % (b+1)]
    for q in range(NPLAT):
        s = '' if q == 0 else '2'
        names += ['f_stretch_r_p%s_amp' % s, 'f_stretch_r_p%s_xa' % s,
                  'f_stretch_r_p%s_xb' % s, 'f_stretch_r_p%s_w' % s]
    for n_, v in zip(names, p):
        print('%s = %.12e' % (n_, v))
