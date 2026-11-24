#!/usr/bin/env python3
"""ohmdrag-1029: reference solution of the kinematic Ohmic-drag problem.

Shell r in [1,2], dipole B_p = (2 m cos/r^3, m sin/r^3) held fixed (v_p = 0, B_p
potential), shellular rotation Omega(r) = Om0 sin^4(pi (r-1)), eta = eta(r).  The induced
field is exactly B_phi = b(r,t) sin(th) cos(th) (a single l = 2 mode) with u = r b obeying

    u_t = (eta u')' - 6 eta u / r^2 + 2 m Omega'(r) / r,   u' = 0 at r = 1, 2 (conducting
    walls at rest: E_theta = eta J_theta = 0),  u(r,0) = 0.

Solved EXACTLY in time on an N-cell 2nd-order finite-volume grid by eigendecomposition
of the symmetric tridiagonal operator, numpy only
(u(t) = V diag((e^{lam t}-1)/lam) V^T s); time
integrals of u^2 and u'^2 (Joule heating) in closed form.  Grid convergence is checked by
running two N.

usage: ref1d.py --eta_bot 0.1 --ndec 4 --rt 1.6 --w 0.06 --om0 0.05 --m 1 --T 2 --N 2000
          --out ref.dat [--series ref_series.dat --nt 101]
"""
import argparse
import numpy as np


def eta_of_r(r, a):
    return a.eta_bot*10.0**(-a.ndec*0.5*(1.0 + np.tanh((r - a.rt)/a.w)))


def domega(r, a):
    x = np.pi*(r - 1.0)
    return 4.0*np.pi*a.om0*np.sin(x)**3*np.cos(x)


def build(a, N):
    h = 1.0/N
    rc = 1.0 + h*(np.arange(N) + 0.5)
    rf = 1.0 + h*np.arange(N + 1)
    ef = eta_of_r(rf, a)
    ef[0] = 0.0
    ef[-1] = 0.0                      # zero diffusive flux through the walls (u' = 0)
    ec = eta_of_r(rc, a)
    diag = -(ef[:-1] + ef[1:])/h**2 - 6.0*ec/rc**2
    off = ef[1:-1]/h**2
    s = 2.0*a.m*domega(rc, a)/rc
    return rc, h, diag, off, s


def E(x, T):
    # (e^{xT}-1)/x, x <= 0, stable
    out = np.empty_like(x)
    small = np.abs(x*T) < 1e-8
    out[small] = T
    out[~small] = np.expm1(x[~small]*T)/x[~small]
    return out


def solve(a, N, T, series_t=None):
    rc, h, diag, off, s = build(a, N)
    A = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    lam, V = np.linalg.eigh(A)
    c = V.T @ s
    g = np.expm1(lam*T)/lam
    u = V @ (g*c)
    # derivative operator (cell centres), u' = 0 at the walls
    D = np.zeros((N, N))
    i = np.arange(1, N - 1)
    D[i, i + 1] = 0.5/h
    D[i, i - 1] = -0.5/h
    # one-sided with the wall value u'(wall) = 0: quadratic through wall slope 0
    D[0, 0] = -0.5/h
    D[0, 1] = 0.5/h
    D[-1, -1] = 0.5/h
    D[-1, -2] = -0.5/h
    du = D @ u
    # time integrals: int g_k g_l dt
    L = lam[:, None]
    M = lam[None, :]
    G = (E(L + M, T) - E(L, T) - E(M, T) + T)/(L*M)
    VC = V*c[None, :]
    Ib = np.einsum('ik,kl,il->i', VC, G, VC, optimize=True)/rc**2
    DVC = D @ VC
    Id = np.einsum('ik,kl,il->i', DVC, G, DVC, optimize=True)
    res = dict(r=rc, b=u/rc, up=du, Ib=Ib, Id=Id, h=h)
    if series_t is not None:
        ser = []
        for t in series_t:
            gt = np.expm1(lam*t)/lam if t > 0 else np.zeros_like(lam)
            ut = V @ (gt*c)
            uw1 = (9*ut[0] - ut[1])/8.0              # wall values (u' = 0)
            uw2 = (9*ut[-1] - ut[-2])/8.0
            b1, b2 = uw1/1.0, uw2/2.0
            dut = D @ ut
            eta = eta_of_r(rc, a)
            jdot = 2*np.pi*np.sum(eta*(1.6*(ut/rc)**2 + (4.0/15.0)*dut**2))*h
            emag = np.pi*(4.0/15.0)*np.sum((ut/rc)**2*rc**2)*h   # int B_phi^2/2 dV
            torque = 16*np.pi*a.m*(b2 - b1)/15.0
            ser.append((t, torque, jdot, emag, b1, b2))
        res['series'] = np.array(ser)
    return res


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--eta_bot', type=float, default=0.1)
    p.add_argument('--ndec', type=float, default=4.0)
    p.add_argument('--rt', type=float, default=1.6)
    p.add_argument('--w', type=float, default=0.06)
    p.add_argument('--om0', type=float, default=0.05)
    p.add_argument('--m', type=float, default=1.0)
    p.add_argument('--T', type=float, default=2.0)
    p.add_argument('--N', type=int, default=2000)
    p.add_argument('--out', default='ref.dat')
    p.add_argument('--series', default=None)
    p.add_argument('--nt', type=int, default=101)
    a = p.parse_args()
    res = solve(a, a.N, a.T)
    lo = solve(a, a.N//2, a.T)
    # grid self-check: N/2 vs N interpolated at the N/2 points
    bi = np.interp(lo['r'], res['r'], res['b'])
    ib = np.interp(lo['r'], res['r'], res['Ib'])
    idd = np.interp(lo['r'], res['r'], res['Id'])
    print('# ref self-check N/2 vs N: max|db|/max|b| = %.2e  Ib %.2e  Id %.2e' % (
        np.max(np.abs(bi - lo['b']))/np.max(np.abs(res['b'])),
        np.max(np.abs(ib - lo['Ib']))/np.max(np.abs(res['Ib'])),
        np.max(np.abs(idd - lo['Id']))/np.max(np.abs(res['Id']))))
    eta = eta_of_r(res['r'], a)
    h = res['h']
    Q = 2*np.pi*np.sum(eta*(1.6*res['Ib'] + (4.0/15.0)*res['Id']))*h
    b1 = (9*res['b'][0]*res['r'][0] - res['b'][1]*res['r'][1])/8.0
    b2 = (9*res['b'][-1]*res['r'][-1] - res['b'][-2]*res['r'][-2])/8.0/2.0
    tq = 16*np.pi*a.m*(b2 - b1)/15.0
    hdr = ('ohmdrag ref: eta_bot %g ndec %g rt %g w %g om0 %g m %g T %g N %d\n'
           'Q_total %.10e torque %.10e b_wall1 %.10e b_wall2 %.10e\n'
           'columns: r b u_prime Ib Id' % (a.eta_bot, a.ndec, a.rt, a.w, a.om0, a.m,
                                           a.T, a.N, Q, tq, b1, b2))
    np.savetxt(a.out, np.column_stack([res['r'], res['b'], res['up'], res['Ib'],
                                       res['Id']]), header=hdr, fmt='%.12e')
    print(hdr)
    if a.series:
        ts = np.linspace(0.0, a.T, a.nt)
        sres = solve(a, a.N, a.T, series_t=ts)
        np.savetxt(a.series, sres['series'], fmt='%.12e',
                   header='t torque jouledot emag_phi b_wall1 b_wall2')


if __name__ == '__main__':
    main()
