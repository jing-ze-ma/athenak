#!/usr/bin/env python3
"""tests_r27: radial grid refit for a DEEPER inner wall, plus the dt estimate.

  grid.py hp                       -> H_p and the present grid's dr, for reference
  grid.py fit <rin/R> <N> [dra]    -> refit f_stretch_r_c1..4 on [rin, 1.25 R]
  grid.py dt  <rin/R> <N> c1 c2 c3 c4  -> per-cell CFL on that grid, both ICs
"""
import sys

import numpy as np

import deep

RS = deep.RS
ROUT = 2.964625e11              # 1.25 R, unchanged
R0_OLD = 1.18585e11
C_OLD = [0.992525, -0.404821, -0.372255, -1.499759]
N_OLD = 144


def u_of(c, xi):
    return xi + (1-xi)*sum(c[k]*xi**(k+1) for k in range(len(c)))


def faces(r0, r1, c, n):
    return r0 + (r1-r0)*u_of(c, np.linspace(0, 1, n+1))


F_OLD = faces(R0_OLD, ROUT, C_OLD, N_OLD)
DR_OLD = np.diff(F_OLD)
RC_OLD = 0.5*(F_OLD[1:] + F_OLD[:-1])


def show(f, tag):
    dr = np.diff(f)
    rc = 0.5*(f[1:]+f[:-1])
    print('%s: N=%d  dr[1e8 cm] at r/R: %s' % (tag, len(dr), ' '.join(
        '%.3f:%.1f' % (rc[i]/RS, dr[i]/1e8) for i in range(0, len(dr), 8))))
    m = (rc/RS > 0.97) & (rc/RS < 1.10)
    print('   0.97-1.10 R: %d cells, dr %.2f-%.2f e8 | at 1.02 R dr=%.2f e8'
          % (m.sum(), dr[m].min()/1e8, dr[m].max()/1e8,
             np.interp(1.02*RS, rc, dr)/1e8))
    m2 = (rc/RS > 0.635) & (rc/RS < 0.97)
    print('   FeCZ 0.635-0.97 R: %d cells, dr %.2f-%.2f e8'
          % (m2.sum(), dr[m2].min()/1e8, dr[m2].max()/1e8))
    m3 = rc/RS < 0.50
    if m3.sum():
        print('   deep  < 0.50 R : %d cells, dr %.2f-%.2f e8'
              % (m3.sum(), dr[m3].min()/1e8, dr[m3].max()/1e8))


def hp():
    show(F_OLD, 'PRESENT 144')
    for nm, fn in (('MLT', deep.ICM), ('RAD', deep.ICR)):
        c = deep.col(fn)
        r, p = c['r'], c['p']
        hpv = -p[1:-1]/np.gradient(p, r)[1:-1]
        rr = r[1:-1]
        print('%s  H_p [1e8 cm] : %s' % (nm, ' '.join(
            '%.3f:%.1f' % (x, np.interp(x*RS, rr, hpv)/1e8)
            for x in (0.355, 0.37, 0.39, 0.40, 0.42, 0.45, 0.50, 0.635, 0.90))))


def fit(rinf, N, dra=6.0e8):
    r0 = rinf*RS

    def dr_target(r):
        x = r/RS
        if x < 0.50:
            # deep stable interior: resolve H_p with >= 4 cells, capped
            hpd = min(np.interp(r, RRHP, HPV_MIN), 3.0e9)
            return max(0.25*hpd, 2.5e8)
        old = np.interp(r, RC_OLD, DR_OLD)
        if x < 0.97:
            return old
        if x < 1.00:
            return min(old, 3.9e8)
        if x < 1.10:
            return 3.9e8 + (dra-3.9e8)*min(1.0, (x-1.00)/0.08)
        return dra*np.exp((x-1.10)/0.06)
    f = [r0]
    while f[-1] < ROUT:
        f.append(f[-1] + dr_target(f[-1]))
    f = np.array(f)
    print('target spacing needs %d cells (asked %d)' % (len(f)-1, N))
    s = np.linspace(0, 1, len(f))
    f[-1] = ROUT
    xi = np.linspace(0, 1, 4001)
    ut = (np.interp(xi, s, f) - r0)/(ROUT - r0)
    A = np.stack([(1-xi)*xi**(k+1) for k in range(4)], 1)
    c, _, _, _ = np.linalg.lstsq(A, ut-xi, rcond=None)
    fn = faces(r0, ROUT, c, N)
    drn = np.diff(fn)
    print('c =', ['%+.6f' % v for v in c], 'monotonic:', bool((drn > 0).all()))
    show(fn, 'NEW rin=%.4fR' % rinf)
    return c


def dtcalc(rinf, N, c):
    r0 = rinf*RS
    fn = faces(r0, ROUT, c, N)
    dr = np.diff(fn)
    rc = 0.5*(fn[1:]+fn[:-1])
    for nm, fnm in (('MLT', deep.ICM), ('RAD', deep.ICR)):
        col = deep.col(fnm)
        cs = np.interp(rc, col['r'], col['cs'])
        for nang in (96, 192):
            dth = (np.pi/2)/nang
            inv = cs/dr + 2.0*cs/(rc*dth)
            dt = 0.3/inv
            k = int(np.argmin(dt))
            print('  %s n_ang=%3d : dt_min = %6.3f s at r/R %.3f (dr %.2fe8, '
                  'cs %.2e) | radial-only %.3f s'
                  % (nm, nang, dt[k], rc[k]/RS, dr[k]/1e8, cs[k],
                     (0.3*dr/cs).min()))


if __name__ == '__main__':
    if sys.argv[1] == 'hp':
        hp()
    else:
        cm = deep.col(deep.ICM)
        cr = deep.col(deep.ICR)
        RRHP = cm['r'][1:-1]
        hm = -cm['p'][1:-1]/np.gradient(cm['p'], cm['r'])[1:-1]
        hr = -cr['p'][1:-1]/np.gradient(cr['p'], cr['r'])[1:-1]
        HPV_MIN = np.minimum(hm, np.interp(RRHP, cr['r'][1:-1], hr))
        if sys.argv[1] == 'fit':
            c = fit(float(sys.argv[2]), int(sys.argv[3]),
                    float(sys.argv[4]) if len(sys.argv) > 4 else 6.0e8)
            dtcalc(float(sys.argv[2]), int(sys.argv[3]), c)
        else:
            dtcalc(float(sys.argv[2]), int(sys.argv[3]),
                   [float(v) for v in sys.argv[4:8]])
