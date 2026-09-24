#!/usr/bin/env python3
"""Plane-parallel Milne gate for vet_col_surface_q (runs_5e): the Cartesian grey
radiative-equilibrium slab (inp/milne_vc.athinput: rho kappa_t = 100 on x = 0..0.8,
tau = 80, flux F_in = 1 at the bottom, c = 100, Marshak top) vs the EXACT Milne solution
J = 3 H (tau + q(tau)), H = F_in/c, the Hopf function q(tau) computed here by a VET
iteration with the exact E_n kernels (piecewise-linear S): q(0) = 1/sqrt(3),
q(inf) = 0.710446.

usage: milne_v.py <root> <run> [<run> ...]          L1 / Linf of T_rad/T_ref - 1
       milne_v.py <root> -r <prefix> <n>            Richardson (n, 2n)
"""
import glob
import os
import sys

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import expn

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'runs_5b_sp_s2'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'vis', 'python'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import s2lib  # noqa: E402

c, arad, fin, chi, x1max = 100.0, 1e8, 1.0, 100.0, 0.8
H = fin/c


def seg(n, xa, xb, sa, sb):
    """int_xa^xb S(x) E_n(x) dx, S linear from sa (xa) to sb (xb), elementwise"""
    dx = xb - xa
    sl = np.where(dx > 0, (sb - sa)/np.where(dx > 0, dx, 1.0), 0.0)
    i0 = expn(n+1, xa) - expn(n+1, xb)
    i1 = -dx*expn(n+1, xb) - expn(n+2, xb) + expn(n+2, xa)
    return sa*i0 + sl*i1


def moments(t, s):
    """J, H (outward), K at the nodes t of a piecewise-linear S (semi-infinite, S linear
    beyond t[-1])"""
    nt = len(t)
    J = np.zeros(nt)
    Hh = np.zeros(nt)
    K = np.zeros(nt)
    bt = (s[-1] - s[-2])/(t[-1] - t[-2])
    for k in range(nt):
        tau = t[k]
        # deeper side: x = t - tau
        xa, xb = t[k:-1] - tau, t[k+1:] - tau
        sa, sb = s[k:-1], s[k+1:]
        xn = t[-1] - tau
        dj = seg(1, xa, xb, sa, sb).sum() + s[-1]*expn(2, xn) + bt*expn(3, xn)
        dh = seg(2, xa, xb, sa, sb).sum() + s[-1]*expn(3, xn) + bt*expn(4, xn)
        dk = seg(3, xa, xb, sa, sb).sum() + s[-1]*expn(4, xn) + bt*expn(5, xn)
        # upper side: x = tau - t
        if k > 0:
            xa, xb = tau - t[1:k+1], tau - t[:k]
            sa, sb = s[1:k+1], s[:k]
            uj = seg(1, xa, xb, sa, sb).sum()
            uh = seg(2, xa, xb, sa, sb).sum()
            uk = seg(3, xa, xb, sa, sb).sum()
        else:
            uj = uh = uk = 0.0
        J[k], Hh[k], K[k] = 0.5*(dj + uj), 0.5*(dh - uh), 0.5*(dk + uk)
    return J, Hh, K


def hopf(cache):
    if cache and os.path.exists(cache):
        d = np.load(cache)
        return d['t'], d['q']
    t = np.concatenate(([0.0], np.logspace(-6, np.log10(120.0), 1200)))
    jj = 3.0*(t + 0.71)
    for it in range(30):
        J, Hh, K = moments(t, jj)
        f = K/J
        g = Hh[0]/J[0]
        jn = (f[0]/g + t)/f
        dq = np.max(np.abs(jn - jj)/jn)
        jj = jn
        if dq < 1e-10:
            break
    q = jj/3.0 - t
    if cache:
        np.savez(cache, t=t, q=q)
    return t, q


here = os.path.dirname(os.path.abspath(__file__))
th, qh = hopf(os.path.join(here, 'milne_hopf.npz'))
qspl = CubicSpline(th, qh)


def eref(x):
    tau = chi*(x1max - x)
    return 3.0*H*(tau + qspl(tau))


def load(root, run):
    fm = sorted(glob.glob(f'{root}/{run}/tab/*.m1.*.tab'))[-1]
    d = s2lib.load_tab(fm)
    e = d['m1_e']
    n = e.shape[0]
    x = (np.arange(n) + 0.5)*x1max/n
    return x, e, d['time']


def main():
    root = sys.argv[1]
    print(f"# Hopf q(0) {qh[0]:.6f} (1/sqrt3 0.577350), q(0.01) {qspl(0.01):.6f} "
          f"(0.588236), q(0.1) {qspl(0.1):.6f} (0.627919), q(1) {qspl(1.0):.6f} "
          f"(0.697729), q(inf) {qh[-1]:.6f} (0.710446)")
    if sys.argv[2] == '-r':
        pre, n = sys.argv[3], int(sys.argv[4])
        x1, e1, _ = load(root, f'{pre}{n}')
        x2, e2, _ = load(root, f'{pre}{2*n}')
        er = 2.0*e2 - np.exp(CubicSpline(x1, np.log(e1))(x2))
        sel = (x2 > x1[0]) & (x2 < x1[-1])
        tr = (eref(x2[sel])/arad)**0.25
        dv = np.abs((er[sel]/arad)**0.25/tr - 1.0)
        print(f"{pre} n={n},{2*n}: Richardson L1 {np.mean(dv):.3e} Linf {np.max(dv):.3e}")
        return
    for run in sys.argv[2:]:
        x, e, t = load(root, run)
        tr = (eref(x)/arad)**0.25
        dv = np.abs((e/arad)**0.25/tr - 1.0)
        print(f"{run:8s} n={len(x):3d} t={t:.4g}  L1 {np.mean(dv):.3e} Linf {np.max(dv):.3e}"
              f"  top E/E_ref-1 {e[-1]/eref(x[-1])-1:+.3e}")


main()
