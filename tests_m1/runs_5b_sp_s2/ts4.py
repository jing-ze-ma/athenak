#!/usr/bin/env python3
"""T-S4: steady grey extended atmosphere on the sp wedge vs a 1-D spherical reference.

The reference integrates the steady moment equations of the SAME closure inward from
the Marshak outer face,
    d(chi E)/dr = -(3 chi - 1) E/r - rho kappa_t F/c,   F = F_in r_in^2/r^2,
    F(r_out) = c q E(r_out),   chi = chi(f), f = F/(c E),
with d(chi E)/dr = (chi - f chi') dE/dr + chi' dF/dr / c.  Radiative equilibrium: the
gas is at T_rad = (E/a)^(1/4).  The face luminosity is rebuilt from the cell means
(F_{i+1/2} = 2 F1_i - F_{i-1/2}, F_{1/2} = F_in).

usage: ts4.py closure(m1|eddington|kershaw) c a F_in q kap_t rho0 n rin rout
              poly|none m1dump hydrodump [m1dump hydrodump ...]
"""
import sys

import numpy as np
from scipy.integrate import solve_ivp

import s2lib


def chi_of(f, kind):
    f2 = min(f*f, 1.0)
    if kind == 'eddington':
        return 1.0/3.0
    if kind == 'kershaw':
        return (1.0 + 2.0*f2)/3.0
    s = np.sqrt(max(4.0 - 3.0*f2, 0.0))
    return (3.0 + 4.0*f2)/(5.0 + 2.0*s)


def reference(kind, c, fin, q, kap, rho0, n, rin, rout):
    def flux(r):
        return fin*rin*rin/(r*r)

    def rhs(r, y):
        e = y[0]
        ff = flux(r)
        f = ff/(c*e)
        h = 1.0e-6
        chi = chi_of(f, kind)
        dchi = (chi_of(f + h, kind) - chi_of(f - h, kind))/(2.0*h)
        dfdr = -2.0*ff/r
        dg = -(3.0*chi - 1.0)*e/r - rho0*(r/rin)**(-n)*kap*ff/c
        return [(dg - dchi*dfdr/c)/(chi - f*dchi)]
    e0 = flux(rout)/(c*q)
    sol = solve_ivp(rhs, (rout, rin), [e0], rtol=1e-12, atol=1e-14,
                    dense_output=True, method='DOP853')
    return lambda r: sol.sol(r)[0]


def main():
    a = sys.argv
    kind = a[1]
    c, arad, fin, q, kap, rho0, n, rin, rout = map(float, a[2:11])
    poly = None if a[11] == 'none' else s2lib.HE4_POLY
    ref = reference(kind, c, fin, q, kap, rho0, n, rin, rout)
    res = []
    for fm, fh in zip(a[12::2], a[13::2]):
        d = s2lib.load_tab(fm)
        h = s2lib.load_tab(fh)
        e = d['m1_e']
        nn = e.shape[0]
        rf, x1v = s2lib.rgrid(rin, rout, nn, poly)
        er = ref(x1v)
        tr = (er/arad)**0.25
        trad = (e/arad)**0.25
        tgas = (2.0/3.0)*h['eint']/h['dens']   # gamma = 5/3, T = p/rho
        l1 = np.mean(np.abs(trad - tr)/tr)
        linf = np.max(np.abs(trad - tr)/tr)
        f1 = d['m1_f1']
        ff = np.zeros(nn + 1)
        ff[0] = fin
        for i in range(nn):
            ff[i+1] = 2.0*f1[i] - ff[i]
        lum = rf*rf*ff/(rin*rin*fin)
        ldev = np.max(np.abs(lum - 1.0))
        gdev = np.max(np.abs(tgas - trad)/trad)
        res.append((nn, l1, linf))
        print(f"{fm}: n={nn} t={d['time']:.4g} L1(T/Tref-1)={l1:.4e} Linf={linf:.4e} "
              f"max|L_f/L_in-1|={ldev:.3e} max|Tgas/Trad-1|={gdev:.2e} "
              f"f_out={f1[-1]/(c*e[-1]):.4f}")
    for x, y in zip(res[:-1], res[1:]):
        print(f"order n={x[0]}->{y[0]}: L1 {np.log(x[1]/y[1])/np.log(y[0]/x[0]):.3f} "
              f"Linf {np.log(x[2]/y[2])/np.log(y[0]/x[0]):.3f}")


main()
