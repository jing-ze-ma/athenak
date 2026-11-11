#!/usr/bin/env python3
"""Exact LINEAR dispersion relations of the radiation-modified acoustic wave.

Background: rho0, T0 (P_gas = rho0 T0, ideal gas, gamma), E0 = a T0^4, F0 = 0,
v = 0, one grey opacity kappa (kappa_P = kappa_E = kappa_F), no scattering, the
true c.  Perturbations ~ exp(i(k x - omega t)); state q = (rho1, v1, T1, E1, F1).
The gas + moment equations are the <rad_m1> mixed-frame system
(docs/dev/rad_m1_design.md sect. 1) linearised about the static equilibrium:

  d rho1/dt = -i k rho0 v1
  d v1/dt   = -i k (T0 rho1 + rho0 T1)/rho0 + kappa (F1 - (4/3) E0 v1)/c
  d T1/dt   = -i k (gamma-1) T0 v1 + (gamma-1) c kappa (E1 - 4 a T0^3 T1)
  d E1/dt   = -i k F1 - c rho0 kappa (E1 - 4 a T0^3 T1)
  d F1/dt   = -i k c^2 P1 - c rho0 kappa (F1 - (4/3) E0 v1)

with three closures for P1:

  edd     P1 = E1/3                                   (closure = eddington, and M1
                                                       at linear order)
  vetqs   P1 = E1/3 + E0 dD, dD = (K1 - J1/3)/E0 from a STATIC formal solution
          mu.grad I = chi (S - I), S = a T^4 (the vet_sc model of rad_m1_vet.cpp:
          no 1/c dI/dt, no Doppler terms in S), with continuous angles ('vetqs')
          or the code's double-Gauss(nmu) x nphi quadrature ('vetqs_d')
  transp  the EXACT time-dependent grey transfer equation to O(v/c),
          (1/c) dI/dt + mu dI/dx = kappa rho [ S (1 + 3 mu beta) - I (1 - mu beta) ]
          (lab frame, comoving isotropic grey emission/absorption), whose
          first two moments are the system above with the exact closure.

The acoustic root is the right-going root whose phase speed lies nearest the
interval [isothermal gas speed, adiabatic mixture speed].
"""

import argparse
import json
import math

import numpy as np


class Bg(object):
    """Background and parameters, all in code units."""

    def __init__(self, prat, tau_lam, lam=1.0, c=1.0e3, gamma=5.0 / 3.0,
                 rho0=1.0, t0=1.0):
        self.rho0, self.t0, self.gamma, self.c = rho0, t0, gamma, c
        self.pgas = rho0 * t0
        self.arad = 3.0 * prat * self.pgas / t0 ** 4
        self.e0 = self.arad * t0 ** 4
        self.lam = lam
        self.k = 2.0 * math.pi / lam
        self.kappa = tau_lam / (rho0 * lam)
        self.prat = prat
        self.tau_lam = tau_lam

    def mixture_cs(self):
        ptot = self.pgas + self.e0 / 3.0
        bet = self.pgas / ptot
        gm1 = self.gamma - 1.0
        g1 = bet + (4.0 - 3.0 * bet) ** 2 * gm1 / (bet + 12.0 * gm1
                                                   * (1.0 - bet))
        return math.sqrt(g1 * ptot / self.rho0)

    def speeds(self):
        """(isothermal gas, adiabatic gas, equilibrium mixture) sound speeds."""
        return (math.sqrt(self.t0), math.sqrt(self.gamma * self.t0),
                self.mixture_cs())


def matrix(bg, q_vet=0.0, form="lab"):
    """d q/dt = L q for the moment system; q_vet = coefficient of T1 in P1.

    form = 'lab': the flux variable is the lab F with d F/dt = -c^2 dP/dx - c rho kappa
    (F - (4/3) E0 v) (explicit transport).  form = 'f0': the variable is the comoving
    F0 = F - (4/3) E0 v evolved as d F0/dt = -c^2 dP/dx - c rho kappa F0, i.e. WITHOUT
    the -(4/3) E0 dv/dt term: the face-flux equation of the implicit transport
    (rad_m1_implicit.cpp, 'fn = th*(f0n - c^2 dt grad(D E) ...)'), with the energy flux
    F0 + (4/3) E0 v.  The returned eigenvector always carries the LAB F in slot 4."""
    if form == "f0":
        return matrix_f0(bg, q_vet)
    r0, t0, g, c, kap, k = bg.rho0, bg.t0, bg.gamma, bg.c, bg.kappa, bg.k
    e0, a = bg.e0, bg.arad
    b = 4.0 * a * t0 ** 3
    ik = 1j * k
    L = np.zeros((5, 5), dtype=complex)
    L[0, 1] = -ik * r0
    L[1, 0] = -ik * t0 / r0
    L[1, 2] = -ik
    L[1, 1] = -kap * (4.0 / 3.0) * e0 / c
    L[1, 4] = kap / c
    L[2, 1] = -ik * (g - 1.0) * t0
    L[2, 3] = (g - 1.0) * c * kap
    L[2, 2] = -(g - 1.0) * c * kap * b
    L[3, 4] = -ik
    L[3, 3] = -c * r0 * kap
    L[3, 2] = c * r0 * kap * b
    L[4, 3] = -ik * c * c / 3.0
    L[4, 2] = -ik * c * c * q_vet
    L[4, 4] = -c * r0 * kap
    L[4, 1] = c * r0 * kap * (4.0 / 3.0) * e0
    return L


def matrix_f0(bg, q_vet=0.0):
    r0, t0, g, c, kap, k = bg.rho0, bg.t0, bg.gamma, bg.c, bg.kappa, bg.k
    e0, a = bg.e0, bg.arad
    b = 4.0 * a * t0 ** 3
    ik = 1j * k
    L = np.zeros((5, 5), dtype=complex)     # q = (rho, v, T, E, F0)
    L[0, 1] = -ik * r0
    L[1, 0] = -ik * t0 / r0
    L[1, 2] = -ik
    L[1, 4] = kap / c
    L[2, 1] = -ik * (g - 1.0) * t0
    L[2, 3] = (g - 1.0) * c * kap
    L[2, 2] = -(g - 1.0) * c * kap * b
    L[3, 4] = -ik
    L[3, 1] = -ik * (4.0 / 3.0) * e0
    L[3, 3] = -c * r0 * kap
    L[3, 2] = c * r0 * kap * b
    L[4, 3] = -ik * c * c / 3.0
    L[4, 2] = -ik * c * c * q_vet
    L[4, 4] = -c * r0 * kap
    return L


def pick_acoustic(bg, omegas):
    lo, _, hi = bg.speeds()
    lo, hi = min(lo, hi) * 0.5, max(lo, hi) * 1.5
    best, bd = None, None
    for w in omegas:
        if w.real <= 0.0:
            continue
        cph = w.real / bg.k
        d = 0.0 if lo <= cph <= hi else min(abs(cph - lo), abs(cph - hi))
        d += 1e-12 * abs(w.imag)
        if bd is None or d < bd:
            best, bd = w, d
    return best


def moment_root(bg, q_vet=0.0, form="lab"):
    """(omega, eigenvector normalised to rho1/rho0 = 1) of the moment system."""
    L = matrix(bg, q_vet, form)
    lam, vec = np.linalg.eig(L)
    omegas = 1j * lam                 # -i omega q = L q
    w = pick_acoustic(bg, omegas)
    idx = int(np.argmin(np.abs(omegas - w)))
    v = vec[:, idx] / vec[0, idx] * bg.rho0
    if form == "f0":
        v = v.copy()
        v[4] = v[4] + (4.0 / 3.0) * bg.e0 * v[1]     # F = F0 + (4/3) E0 v
    return omegas[idx], v


def jfun(a, k):
    """J_n = (1/2) int_{-1}^{1} mu^n/(a + i k mu) dmu, n = 0, 1, 2, 3."""
    j0 = np.arctan(k / a) / k
    j1 = (1.0 - a * j0) / (1j * k)
    j2 = -a * j1 / (1j * k)
    j3 = (1.0 / 3.0 - a * j2) / (1j * k)
    return j0, j1, j2, j3


def jfun_quad(a, k, mus, wts):
    """The same moments with a discrete angular quadrature (weights sum to 1)."""
    den = a + 1j * k * mus
    return tuple(np.sum(wts * mus ** n / den) for n in range(4))


def quad_code(nmu, nphi, kdir):
    """(mu_k, w) of the code's double-Gauss(nmu) x nphi rays for direction kdir."""
    xg, wg = np.polynomial.legendre.leggauss(nmu)
    xg = 0.5 * (xg + 1.0)
    wg = 0.5 * wg
    mus, wts = [], []
    kd = np.asarray(kdir, dtype=float)
    kd = kd / np.linalg.norm(kd)
    for h in (1.0, -1.0):
        for a in range(nmu):
            m1 = h * xg[a]
            st = math.sqrt(1.0 - m1 * m1)
            for ll in range(nphi):
                ph = 2.0 * math.pi * (ll + 0.5) / nphi
                om = np.array([m1, st * math.cos(ph), st * math.sin(ph)])
                mus.append(float(om @ kd))
                wts.append(0.5 * wg[a] / nphi)
    return np.array(mus), np.array(wts)


def vetqs_coef(bg, quad=None):
    """q_vet: P1 = E1/3 + q_vet T1 from the static thermal-source formal solution."""
    chi = bg.rho0 * bg.kappa
    b = 4.0 * bg.arad * bg.t0 ** 3
    if quad is None:
        j0, _, j2, _ = jfun(chi + 0j, bg.k)
    else:
        j0, _, j2, _ = jfun_quad(chi + 0j, bg.k, *quad)
        # the quadrature's own background: J = S, K = S sum(w mu^2) (= 1/3 exactly)
    return chi * b * (j2 - j0 / 3.0)


def transport_det(bg, w, quad=None):
    """Gas 3x3 determinant with E1, F1 from the exact transfer equation."""
    r0, t0, g, c, kap, k = bg.rho0, bg.t0, bg.gamma, bg.c, bg.kappa, bg.k
    e0 = bg.e0
    b = 4.0 * bg.arad * t0 ** 3
    chi = r0 * kap
    a = chi - 1j * w / c
    if quad is None:
        j0, j1, j2, _ = jfun(a, k)
    else:
        j0, j1, j2, _ = jfun_quad(a, k, *quad)
    # E1 = chi [b J0 T1 + 4 E0 J1 v1/c],  F1 = c chi [b J1 T1 + 4 E0 J2 v1/c]
    eT, ev = chi * b * j0, chi * 4.0 * e0 * j1 / c
    fT, fv = c * chi * b * j1, c * chi * 4.0 * e0 * j2 / c
    ik = 1j * k
    M = np.zeros((3, 3), dtype=complex)   # rows: cont, mom, energy; cols rho, v, T
    M[0, 0] = -1j * w
    M[0, 1] = ik * r0
    M[1, 0] = ik * t0 / r0
    M[1, 1] = -1j * w - kap * (fv - (4.0 / 3.0) * e0) / c
    M[1, 2] = ik - kap * fT / c
    M[2, 1] = ik * (g - 1.0) * t0 - (g - 1.0) * c * kap * ev
    M[2, 2] = -1j * w - (g - 1.0) * c * kap * (eT - b)
    return M, (eT, ev, fT, fv)


def transport_root(bg, w0, quad=None):
    """Newton (secant) on det M(omega) = 0 from w0; (omega, eigvec rho,v,T,E,F)."""
    def f(w):
        return np.linalg.det(transport_det(bg, w, quad)[0])
    wa, wb = w0, w0 * (1.0 + 1e-6) + 1e-9
    fa, fb = f(wa), f(wb)
    for _ in range(200):
        if fb == fa:
            break
        wn = wb - fb * (wb - wa) / (fb - fa)
        wa, fa = wb, fb
        wb, fb = wn, f(wn)
        if abs(wb - wa) < 1e-14 * abs(wb):
            break
    M, (eT, ev, fT, fv) = transport_det(bg, wb, quad)
    # null vector with rho1 = rho0
    u, s, vh = np.linalg.svd(M)
    nv = vh.conj()[-1]
    nv = nv / nv[0] * bg.rho0
    rho1, v1, t1 = nv
    e1 = eT * t1 + ev * v1
    f1 = fT * t1 + fv * v1
    return wb, np.array([rho1, v1, t1, e1, f1]), s[-1] / s[0]


def all_roots(bg, kdir=(1, 0, 0), nmu=4, nphi=8):
    """Every reference root for one (P_rad/P_gas, tau_lambda)."""
    out = {}
    w, v = moment_root(bg)
    out["edd"] = (w, v)
    out["eddc"] = moment_root(bg, 0.0, "f0")
    qv = vetqs_coef(bg)
    out["vetqs"] = moment_root(bg, qv)
    out["vetqsc"] = moment_root(bg, qv, "f0")
    qd = vetqs_coef(bg, quad_code(nmu, nphi, kdir))
    out["vetqs_d"] = moment_root(bg, qd)
    out["vetqsc_d"] = moment_root(bg, qd, "f0")
    wt, vt, res = transport_root(bg, w)
    out["transp"] = (wt, vt)
    return out


def eig_params(bg, v):
    """<problem>/radwave_eig_* for an eigenvector (rho1, v1, T1, E1, F1)."""
    rho1, v1, t1, e1, f1 = v / v[0] * bg.rho0
    vals = {"v": v1 / bg.rho0, "t": t1 / bg.t0, "e": e1 / bg.e0,
            "f": f1 / (bg.c * bg.e0)}
    # rho1 = rho0 * 1 per unit amplitude; v^ etc. per unit drho/rho
    vals = {k: complex(x) for k, x in vals.items()}
    out = {}
    for kk, x in vals.items():
        out["radwave_eig_%s_re" % kk] = x.real
        out["radwave_eig_%s_im" % kk] = x.imag
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--c", type=float, default=1.0e3)
    p.add_argument("--json", default=None)
    args = p.parse_args()
    rows = {}
    keys = ("edd", "eddc", "vetqsc", "vetqsc_d", "transp")
    print("%6s %7s |" % ("Pr/Pg", "tau") + "".join(" %-25s" % kk for kk in keys))
    for prat in (0.1, 1.0, 10.0, 100.0):
        for tau in (0.1, 10.0, 1.0e3):
            bg = Bg(prat, tau, c=args.c)
            r = all_roots(bg)
            s = "%6g %7g |" % (prat, tau)
            for key in keys:
                w = r[key][0]
                s += " %12.8f %12.5e" % (w.real / bg.k, w.imag)
            print(s)
            rows["%g_%g" % (prat, tau)] = {
                kk: [r[kk][0].real, r[kk][0].imag] for kk in r}
            rows["%g_%g" % (prat, tau)]["speeds"] = bg.speeds()
    if args.json:
        with open(args.json, "w") as fp:
            json.dump(rows, fp, indent=1)


if __name__ == "__main__":
    main()
