#!/usr/bin/env python3
"""tests_r27: NONLINEAR refit of the 4-coefficient radial stretch for a deeper wall.

tests_r13/fitgrid.py fits u(xi) by linear least squares; with two fine regions (the
new inner edge and the photosphere) that fit goes non-monotonic.  Here the same
4 coefficients are fitted to log dr directly, with a monotonicity penalty.

usage: fit2.py <rin/R> <N> [hpfac] [dr_phot] [dr_top]
"""
import sys

import numpy as np
from scipy.optimize import least_squares

import deep
import grid

RS = deep.RS
ROUT = grid.ROUT


def run(rinf, N, hpfac=0.25, drph=3.9e8, drtop=6.0e8):
    r0 = rinf*RS
    cm, cr = deep.col(deep.ICM), deep.col(deep.ICR)
    rr = cm['r'][1:-1]
    hm = -cm['p'][1:-1]/np.gradient(cm['p'], cm['r'])[1:-1]
    hr = -cr['p'][1:-1]/np.gradient(cr['p'], cr['r'])[1:-1]
    hp = np.minimum(hm, np.interp(rr, cr['r'][1:-1], hr))

    def target(r):
        x = r/RS
        deepv = np.maximum(hpfac*np.interp(r, rr, hp), 2.5e8)
        old = np.interp(r, grid.RC_OLD, grid.DR_OLD)
        out = np.where(x < 0.50, deepv, old)
        out = np.where(x < 0.97, out, drph)
        out = np.where(x < 1.00, out,
                       drph + (drtop-drph)*np.clip((x-1.00)/0.08, 0, 1))
        out = np.where(x < 1.10, out, drtop*np.exp((x-1.10)/0.06))
        return out

    xic = (np.arange(N)+0.5)/N

    def resid(c):
        f = grid.faces(r0, ROUT, c, N)
        dr = np.diff(f)
        rc = 0.5*(f[1:]+f[:-1])
        bad = np.minimum(dr, 0.0)
        tg = target(np.maximum(rc, r0))
        w = np.where(rc/RS < 0.50, 6.0, 1.0)
        r1 = w*(np.log(np.maximum(dr, 1.0))-np.log(tg))
        return np.concatenate((r1, 1e3*bad/1e8))

    best = None
    for s in range(12):
        c0 = np.array([0.99, -0.40, -0.37, -1.50]) + 0.5*np.random.RandomState(s).randn(4)
        try:
            sol = least_squares(resid, c0, max_nfev=4000)
        except Exception:
            continue
        f = grid.faces(r0, ROUT, sol.x, N)
        if (np.diff(f) <= 0).any():
            continue
        if best is None or sol.cost < best.cost:
            best = sol
    if best is None:
        print('NO MONOTONIC FIT')
        return None
    c = best.x
    print('rin=%.4f R = %.6e cm   N=%d   c = %s'
          % (rinf, r0, N, ' '.join('%+.6f' % v for v in c)))
    grid.show(grid.faces(r0, ROUT, c, N), '  ')
    grid.dtcalc(rinf, N, c)
    print('  rad_flux_inner = %.6e' % (2.3066e38/(4*np.pi*r0**2)))
    return c


if __name__ == '__main__':
    run(float(sys.argv[1]), int(sys.argv[2]),
        *[float(v) for v in sys.argv[3:]])
