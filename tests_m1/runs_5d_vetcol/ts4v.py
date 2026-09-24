#!/usr/bin/env python3
"""T-S4 with vet_col: the grey spherical atmosphere (runs_5b_sp_s2 inp/sph_atm) vs the
EXACT-TRANSFER 1-D reference and the moment references of the M1 and Eddington closures.

Exact transfer (radiative equilibrium S = J, H = F_in r_in^2/(c r^2)): VET iteration
    d(f J)/dr + (3 f - 1) J / r = -chi H,   J(r_out) = H(r_out)/g,
    f = K/J, g = H/J(r_out) from the long-characteristics formal solution of sphref.py
    with S = J (inner sphere: I = J + 3 H mu, top: vacuum), until f changes < 1e-8.
The Marshak variant of the same reference (J(r_out) = H(r_out)/q, the code's outer BC,
f from the formal solution) is also given: it isolates the closure from the outer BC.

usage: ts4v.py <root> <run> [<run> ...]   (root/run/tab: last m1 tab slice)
"""
import glob
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'runs_5b_sp_s2'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'vis', 'python'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import s2lib  # noqa: E402
from sphref import moments  # noqa: E402

c, arad, fin, q, kap, rho0, en, rin, rout = 100.0, 1e8, 1.0, 0.5, 100.0, 1.0, 2.0, 1.0, 5.0
C = rho0*kap*rin**2


def hflux(r):
    return fin*rin*rin/(c*r*r)


def solve_j(fs, gtop):
    def rhs(r, y):
        f = fs(r)
        return [-(3.0*f - 1.0)/(f*r)*y[0] - rho0*kap*(r/rin)**(-en)*hflux(r)]
    y0 = fs(rout)*hflux(rout)/gtop
    sol = solve_ivp(rhs, (rout, rin), [y0], rtol=1e-12, atol=1e-16, dense_output=True,
                    method='DOP853')
    return lambda r: sol.sol(r)[0]/fs(r)


def exact(marshak=False, cache=None):
    if cache and os.path.exists(cache):
        d = np.load(cache)
        return CubicSpline(d['r'], d['J'])
    rg = np.concatenate([np.linspace(rin, 4.8, 97), np.linspace(4.81, rout, 20)])
    f = np.full(rg.size, 1.0/3.0)
    g = q
    for it in range(40):
        fs = CubicSpline(rg, f)
        J = solve_j(fs, q if marshak else g)
        Sf = lambda r: J(np.clip(r, rin, rout))
        fn = np.zeros(rg.size)
        for n, r in enumerate(rg):
            jj, hh, kk = moments(r, rin, rout, C, Sf, J(rin), c*hflux(rin), c, ng=32,
                                 nt=1201)
            fn[n] = kk/jj
            if n == rg.size - 1:
                g = hh/jj
        dmax = np.max(np.abs(fn - f))
        f = fn
        print(f"# VET iteration {it}: max|df| {dmax:.2e}  g_top {g:.5f}  f_top {f[-1]:.5f}",
              flush=True)
        if dmax < 1e-8:
            break
    Jg = solve_j(CubicSpline(rg, f), q if marshak else g)(rg)
    if cache:
        np.savez(cache, r=rg, J=Jg, f=f, g=g)
    return CubicSpline(rg, Jg)


def main():
    root = sys.argv[1]
    here = os.path.dirname(os.path.abspath(__file__))
    ref = exact(False, os.path.join(here, 'ts4_exact.npz'))
    refm = exact(True, os.path.join(here, 'ts4_exact_marshak.npz'))
    for run in sys.argv[2:]:
        fm = sorted(glob.glob(f'{root}/{run}/tab/*.m1.*.tab'))[-1]
        d = s2lib.load_tab(fm)
        e = d['m1_e']
        nn = e.shape[0]
        poly = s2lib.HE4_POLY if 'vs' in run else None
        rf, x1v = s2lib.rgrid(rin, rout, nn, poly)
        trad = (e/arad)**0.25
        out = [f"{run:8s} n={nn:3d} t={d['time']:.4g}"]
        for nm, rr in (('exact', ref), ('exactMarshak', refm)):
            tr = (rr(x1v)/arad)**0.25
            dv = np.abs(trad/tr - 1.0)
            out.append(f"{nm}: L1 {np.mean(dv):.3e} Linf {np.max(dv):.3e}")
        f1 = d['m1_f1']
        ff = np.zeros(nn + 1)
        ff[0] = fin
        for i in range(nn):
            ff[i+1] = 2.0*f1[i] - ff[i]
        lum = rf*rf*ff/(rin*rin*fin)
        out.append(f"max|L/L_in-1| {np.max(np.abs(lum - 1.0)):.2e}")
        print('  '.join(out))


main()
