#!/usr/bin/env python3
"""le08_shock.py -- Lowrie & Edwards (2008) grey non-equilibrium-diffusion radiative
shock, semi-analytic (ODE) solution.  Written for the IDORT radiative-shock family
(idort-rs-1017); independent of tests_m1/t7_radshock.py (different state variables,
different shock-connection criterion), which it is cross-checked against.

Model (LE08 non-dimensionalisation: rho0 = T0 = 1, v0 = M0, p = rho T / gamma,
E_r = P0 Theta, Theta = T_r^4, F_r = -P0 kappa dTheta/dx):
  integrals     rho v = M0
                rho v^2 + rho T/gamma + P0 Theta/3 = K_m/(3 gamma),
                K_m = 3(gamma M0^2 + 1) + gamma P0
                M0 (v^2/2 + T/(gamma-1)) + (4/3) P0 v Theta + F_r = const
  LE08 ODEs in Mach space (M = v / sqrt(T)), their closed form:
    dx/dM = -6 M0 rho T (M^2 - 1) / ((gamma+1) P0 M Z_d)
    dT/dM = -2 (gamma-1)/(gamma+1) (T/M) Z_n / Z_d
    Z_d = M0 dTheta/dx + (gamma-1)/(gamma+1) (gamma M^2 + 1) r
    Z_n = M0 dTheta/dx + (gamma M^2 - 1) r,     r = 3 rho sigma_a (Theta - T^4)
State here = (M, T): rho = M0/(M sqrt T), Theta from the momentum integral, dTheta/dx
from the energy integral.  Near the equilibria the x-space form (dM/dx, dT/dx), regular
there, is used for the eigen-directions and the first leg; the Mach-space form near the
adiabatic sonic point M = 1, where dM/dx is singular.

Shock connection (LE08 sect. 5): the embedded hydrodynamic shock conserves Theta and
F_r, so the GAS Rankine-Hugoniot jump maps a precursor state (M_p > 1, T_p, Theta_p)
onto (M_s, T_s, Theta_p); the solution is the precursor point whose RH image lies on the
relaxation trajectory, i.e. the root of g(M_p) = T_relax(Theta_p) - T_s(M_p).
(t7_radshock.py instead intersects the two branches in the (Theta, F_r) plane.)  If no
root with M_p > 1 exists and both branches end at the same state at M = 1, the solution
is continuous (no embedded shock; the small-M0 cases).

usage:
  python le08_shock.py --m0 27 [--p0 1e-4 --sigma-a 1e6 --kappa 1] [--write-ref FILE --nref N]
        [--athinput] [--json FILE]          # metrics + reference table (t7 cgs scaling)
  python le08_shock.py --validate           # published / cross-code checks (see VALIDATE)
"""
import argparse
import json
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve

KB = 1.380649e-16
MH = 1.67353e-24
AR = 7.5657332502e-15
CL = 2.99792458e10
# cgs scaling of the IDORT / t7 decks (Skinner & Ostriker 2013 / QUOKKA RadhydroShockCGS)
RHO0, T0, SIG_CGS = 5.69, 2.18e6, 577.16   # deck: kappa_p = kappa_e = kappa_f = 577.16 /cm


def implied_params(gamma=5.0 / 3.0):
    """(P0, sigma_a, kappa=1) of the PHYSICAL problem the cgs decks solve: a0 from
    gamma k T0 / mu, C = c / a0, L chosen so that kappa = C/(3 sigma L) = 1, then
    sigma_a = sigma L C = C^2/3 and P0 = a T0^4 / (rho0 a0^2).  The nominal LE08
    values are P0 = 1e-4, sigma_a = 1e6; the decks are off by +0.19 % / -0.05 %."""
    a0 = np.sqrt(gamma * KB * T0 / MH)
    cc = CL / a0
    return AR * T0 ** 4 / (RHO0 * a0 * a0), cc * cc / 3.0, 1.0


class LE08(object):
    def __init__(self, m0, p0=1e-4, gamma=5.0 / 3.0, sigma_a=1e6, kappa=1.0):
        self.m0, self.p0, self.g, self.sa, self.kap = float(m0), p0, gamma, sigma_a, kappa
        g = gamma
        self.km = 3.0 * (g * m0 ** 2 + 1.0) + g * p0
        self.e0 = m0 * (0.5 * m0 ** 2 + 1.0 / (g - 1.0)) + (4.0 / 3.0) * p0 * m0

    # ---- algebra --------------------------------------------------------------
    def state(self, M, T):
        rho = self.m0 / (M * np.sqrt(T))
        v = self.m0 / rho
        th = (self.km - 3.0 * self.g * self.m0 ** 2 / rho - 3.0 * T * rho) / (self.g * self.p0)
        dth = (self.m0 * (0.5 * v * v + T / (self.g - 1.0)) + (4.0 / 3.0) * self.p0 * v * th
               - self.e0) / (self.p0 * self.kap)
        return rho, v, th, dth

    def zz(self, M, T):
        g = self.g
        rho, v, th, dth = self.state(M, T)
        r = 3.0 * rho * self.sa * (th - T ** 4)
        zd = self.m0 * dth + (g - 1.0) / (g + 1.0) * (g * M * M + 1.0) * r
        zn = self.m0 * dth + (g * M * M - 1.0) * r
        return rho, zd, zn

    def rhs_x(self, M, T):
        """(dM/dx, dT/dx): regular at the equilibria, singular at M = 1."""
        g = self.g
        rho, zd, zn = self.zz(M, T)
        dm = -(g + 1.0) * self.p0 * M * zd / (6.0 * self.m0 * rho * T * (M * M - 1.0))
        dt = (g - 1.0) * self.p0 * zn / (3.0 * self.m0 * rho * (M * M - 1.0))
        return np.array([dm, dt])

    def rhs_m(self, M, T):
        """(dx/dM, dT/dM): LE08 closed form; regular at M = 1 (if Z_d != 0)."""
        g = self.g
        rho, zd, zn = self.zz(M, T)
        dx = -6.0 * self.m0 * rho * T * (M * M - 1.0) / ((g + 1.0) * self.p0 * M * zd)
        dt = -2.0 * (g - 1.0) / (g + 1.0) * (T / M) * zn / zd
        return np.array([dx, dt])

    def downstream(self):
        """Equilibrium jump (rho1, T1, v1): Theta = T^4, F_r = 0 in the 3 integrals.
        1-D: for a given rho the momentum integral fixes T (monotone quartic), then the
        root of the energy residual in rho > 1 (scan + Brent)."""
        g, m0, p0 = self.g, self.m0, self.p0
        pm = self.km / (3.0 * g)

        def t_of(rho):
            rhs = pm - m0 * m0 / rho
            return brentq(lambda T: rho * T / g + p0 * T ** 4 / 3.0 - rhs, 0.0,
                          max(1.0, 10.0 * rhs * g / rho + (3.0 * rhs / p0) ** 0.25),
                          xtol=1e-16, rtol=1e-15)

        def fe(rho):
            T = t_of(rho)
            v = m0 / rho
            return (m0 * (0.5 * v * v + T / (g - 1.0)) + (4.0 / 3.0) * p0 * v * T ** 4) / self.e0 - 1.0
        rr = np.linspace(1.0 + 1e-6, 7.0 - 1e-9, 20001)
        fv = np.array([fe(r) for r in rr])
        k = np.nonzero(np.diff(np.sign(fv)))[0]
        if k.size == 0:
            raise RuntimeError("no downstream equilibrium")
        r1 = brentq(fe, rr[k[-1]], rr[k[-1] + 1], xtol=1e-15, rtol=1e-15)
        return float(r1), float(t_of(r1)), m0 / float(r1)

    def downstream_published(self):
        """LE08 closed-form T1 relation (their eqs. for rho1(T1) and the T1 root),
        as transcribed in QUOKKA extern/LowrieEdwards/radshock.py (shock_jump)."""
        g, p0, m0 = self.g, self.p0, self.m0

        def f1(T1):
            return 3.0 * (g + 1.0) * (T1 - 1.0) - p0 * g * (g - 1.0) * (7.0 + T1 ** 4)

        def f2(T1):
            return 12.0 * (g - 1.0) ** 2 * T1 * (3.0 + g * p0 * (1.0 + 7.0 * T1 ** 4))

        def rho(T1):
            return (f1(T1) + np.sqrt(f1(T1) ** 2 + f2(T1))) / (6.0 * (g - 1.0) * T1)

        def fn(T1):
            r1 = rho(T1)
            return (3.0 * r1 * (r1 * T1 - 1.0) + g * p0 * r1 * (T1 ** 4 - 1.0)
                    - 3.0 * g * (r1 - 1.0) * m0 ** 2)
        tg = np.geomspace(1.0 + 1e-9, 4.0 * m0 ** 2, 20000)
        fv = np.array([fn(t) for t in tg])
        k = np.nonzero(np.diff(np.sign(fv)))[0]
        T1 = brentq(fn, tg[k[-1]], tg[k[-1] + 1], xtol=1e-15, rtol=1e-15)
        return rho(T1), T1, m0 / rho(T1)


def _eig(p, M, T, growing):
    h = 1e-7 * max(1.0, abs(M))
    J = np.zeros((2, 2))
    for k, d in enumerate(([h, 0.0], [0.0, h * T])):
        J[:, k] = (p.rhs_x(M + d[0], T + d[1]) - p.rhs_x(M - d[0], T - d[1])) / (2.0 * np.linalg.norm(d))
    w, V = np.linalg.eig(J)
    idx = [i for i in range(2) if (w[i].real > 0) == growing]
    i = idx[0]
    vec = V[:, i].real
    return vec / np.linalg.norm(vec), w[i].real


def _branch(p, M_eq, T_eq, vec, sgn, eps, toward):
    """Integrate from the equilibrium + eps*vec to M = 1. x-space leg until |M - 1| < 0.05
    or T stalls, then Mach-space leg to M = 1.  Returns arrays (x, M, T) ordered from the
    equilibrium end (x measured from the start point)."""
    M = M_eq + sgn * eps * vec[0]
    T = T_eq + sgn * eps * vec[1]
    xs, Ms, Ts = [0.0], [M], [T]
    direction = toward  # +1: integrate forward in x (precursor), -1: backward (relaxation)
    sw = 0.05 if abs(M_eq - 1.0) > 0.1 else 0.3 * abs(M_eq - 1.0)

    def f(x, y):
        return direction * p.rhs_x(y[0], y[1])

    def near(x, y):
        return abs(y[0] - 1.0) - sw
    near.terminal = True
    span = 1e3 * abs(1.0 / _eig(p, M_eq, T_eq, toward > 0)[1]) + 10.0
    s1 = solve_ivp(f, (0.0, span), [M, T], method="LSODA", rtol=1e-11, atol=1e-14,
                   events=near, dense_output=False, max_step=np.inf)
    xs = list(direction * s1.t)
    Ms = list(s1.y[0])
    Ts = list(s1.y[1])
    # Mach-space leg
    Mstart, Tstart, xstart = Ms[-1], Ts[-1], xs[-1]

    def g(Mv, y):
        d = p.rhs_m(Mv, y[1])
        return [d[0], d[1]]
    mend = 1.0
    s2 = solve_ivp(g, (Mstart, mend), [xstart, Tstart], method="LSODA", rtol=1e-12,
                   atol=1e-14, dense_output=True)
    mm = np.linspace(Mstart, mend, 4001)[1:]
    yy = s2.sol(mm)
    xs += list(yy[0])
    Ms += list(mm)
    Ts += list(yy[1])
    return np.array(xs), np.array(Ms), np.array(Ts), s1.status, s2.status


def rh_gas(g, Mp, Tp):
    Ms2 = (2.0 + (g - 1.0) * Mp ** 2) / (2.0 * g * Mp ** 2 - (g - 1.0))
    Ts = Tp * (2.0 * g * Mp ** 2 - (g - 1.0)) * ((g - 1.0) * Mp ** 2 + 2.0) / ((g + 1.0) ** 2 * Mp ** 2)
    return np.sqrt(Ms2), Ts


class Solution(object):
    def __init__(self, p, eps=1e-7):
        self.p = p
        g = p.g
        self.rho1, self.T1, self.v1 = p.downstream()
        self.M1 = self.v1 / np.sqrt(self.T1)
        vu, self.lam_up = _eig(p, p.m0, 1.0, True)
        if vu[0] > 0:
            vu = -vu            # precursor: M decreases downstream
        vd, self.lam_dn = _eig(p, self.M1, self.T1, False)
        if vd[0] < 0:
            vd = -vd            # relaxation, read upstream from the far field: M increases to 1
        xa, Ma, Ta, *_ = _branch(p, p.m0, 1.0, vu, 1.0, eps, +1)
        xb, Mb, Tb, *_ = _branch(p, self.M1, self.T1, vd, 1.0, eps, -1)
        tha = p.state(Ma, Ta)[2]
        thb = p.state(Mb, Tb)[2]
        self.raw = (xa, Ma, Ta, tha, xb, Mb, Tb, thb)
        # RH image of the precursor in the (Theta, T) plane: (Theta_p, T_s(M_p)); it must
        # meet the relaxation curve (Theta_b, T_b) at a point with the same Mach number.
        sel = (Ma > 1.0 + 1e-9)
        Mp, Tp, thp, xpa = Ma[sel], Ta[sel], tha[sel], xa[sel]
        # Theta from the momentum integral carries round-off ~1e-16 K_m/P0 next to the
        # upstream equilibrium; keep the strictly increasing envelope (drops a few points
        # whose dip is < 1e-7 Theta_max; anything larger is an error)
        keep = thp > np.maximum.accumulate(np.concatenate([[-np.inf], thp[:-1]]))
        dip = np.maximum.accumulate(thp) - thp
        if dip.max() > 1e-7 * thp.max():
            raise RuntimeError("Theta not monotone on the precursor")
        Mp, Tp, thp, xpa = Mp[keep], Tp[keep], thp[keep], xpa[keep]
        Msh, Tsh = rh_gas(g, Mp, Tp)
        cands = []
        for j in range(thb.size - 1):
            t0_, t1_ = thb[j], thb[j + 1]
            lo, hi = min(t0_, t1_), max(t0_, t1_)
            if hi < thp[0] or lo > thp[-1] or hi == lo:
                continue
            h0 = Tb[j] - np.interp(t0_, thp, Tsh)
            h1 = Tb[j + 1] - np.interp(t1_, thp, Tsh)
            if h0 == 0.0 or h0 * h1 < 0.0:
                w = h0 / (h0 - h1)
                th = t0_ + w * (t1_ - t0_)
                if th < thp[0] or th > thp[-1]:
                    continue
                mb = Mb[j] + w * (Mb[j + 1] - Mb[j])
                xb_ = xb[j] + w * (xb[j + 1] - xb[j])
                mpi = np.interp(th, thp, Mp)
                ms_ = np.interp(th, thp, Msh)
                cands.append((abs(mb - ms_), th, mpi, mb, xb_))
        good = [c for c in cands if c[0] < 1e-3 and c[2] > 1.0 + 1e-6]
        self.ncand = len(cands)
        self.nroots = len(good)
        self.shock = bool(good)
        if self.shock:
            c = min(good, key=lambda c: c[0])
            self.Ms_check = c[0]
            self.thp = c[1]
            self.Mp = c[2]
            self.Tp = float(np.interp(c[1], thp, Tp))
            self.xp = float(np.interp(c[1], thp, xpa))
            self.Ms, self.Ts = rh_gas(g, self.Mp, self.Tp)
            self.xs_rel = c[4]
            self.Ms_rel = c[3]
        else:
            self.Ms_check = 0.0
            self.Mp = self.Ms = 1.0
            self.Tp = float(Ta[-1])
            self.Ts = float(Tb[-1])
            self.thp = float(tha[-1])
            self.xp = float(xa[-1])
            self.xs_rel = float(xb[-1])
            self.Ms_rel = 1.0
        self.cont_gap = (abs(Ta[-1] - Tb[-1]), abs(tha[-1] - thb[-1]))
        # assemble profile, x = 0 at the shock
        ka = xa <= self.xp
        kb = xb >= self.xs_rel
        self.xa = xa[ka] - self.xp
        self.Ma, self.Ta = Ma[ka], Ta[ka]
        self.xb = xb[kb] - self.xs_rel
        self.Mb, self.Tb = Mb[kb], Tb[kb]
        oa, ob2 = np.argsort(self.xa), np.argsort(self.xb)
        self.xa, self.Ma, self.Ta = self.xa[oa], self.Ma[oa], self.Ta[oa]
        self.xb, self.Mb, self.Tb = self.xb[ob2], self.Mb[ob2], self.Tb[ob2]
        if self.shock:   # add the exact end points at the shock
            self.xa = np.append(self.xa, 0.0)
            self.Ma = np.append(self.Ma, self.Mp)
            self.Ta = np.append(self.Ta, self.Tp)
            self.xb = np.insert(self.xb, 0, 0.0)
            self.Mb = np.insert(self.Mb, 0, self.Ms)
            self.Tb = np.insert(self.Tb, 0, self.Ts)

    # ---- sampled profile -------------------------------------------------------
    def prim(self, M, T):
        rho, v, th, dth = self.p.state(M, T)
        return rho, v, T, np.maximum(th, 0.0) ** 0.25, -self.p.p0 * self.p.kap * dth

    def sample(self, x):
        x = np.asarray(x, float)
        out = np.empty((5, x.size))
        L = x < 0.0
        for msk, xx, MM, TT in ((L, self.xa, self.Ma, self.Ta), (~L, self.xb, self.Mb, self.Tb)):
            if msk.any():
                M = np.interp(x[msk], xx, MM)
                T = np.interp(x[msk], xx, TT)
                out[:, msk] = np.array(self.prim(M, T))
        return out   # rho, v, T, Tr, F

    def conservation(self):
        p = self.p
        M = np.concatenate([self.Ma, self.Mb])
        T = np.concatenate([self.Ta, self.Tb])
        rho, v, th, dth = p.state(M, T)
        F = -p.p0 * p.kap * dth
        c2 = (rho * v * v + rho * T / p.g + p.p0 * th / 3.0) / (p.km / (3 * p.g)) - 1.0
        c3 = (p.m0 * (0.5 * v * v + T / (p.g - 1.0)) + (4.0 / 3.0) * p.p0 * v * th + F) / p.e0 - 1.0
        return float(np.abs(c2).max()), float(np.abs(c3).max())

    def domain(self, frac=1e-4):
        """x range where |T - T_eq| >= frac*(T1 - 1) (both sides)."""
        d1 = frac * (self.T1 - 1.0)
        a = np.nonzero(self.Ta - 1.0 >= d1)[0]
        xmin = float(self.xa[a[0]]) if a.size else float(self.xa[0])
        b = np.nonzero(np.abs(self.Tb - self.T1) >= d1)[0]
        xmax = float(self.xb[b[-1]]) if b.size else float(self.xb[-1])
        return xmin, xmax

    def metrics(self, frac_pre=(0.1, 0.01)):
        out = {"M0": self.p.m0, "rho1": self.rho1, "T1": self.T1, "v1": self.v1,
               "shock": self.shock, "Mp": float(self.Mp), "Ms": float(self.Ms),
               "Tp": float(self.Tp), "Ts": float(self.Ts), "Trs": float(self.thp ** 0.25),
               "Tmax": float(max(self.Ta.max(), self.Tb.max())),
               "Trmax": float(np.max(self.sample(np.concatenate([self.xa, self.xb]))[3]))}
        out["subcritical"] = bool(self.shock and self.Tp < 0.999 * self.T1)
        out["kind"] = ("continuous" if not self.shock else
                       ("subcritical" if out["subcritical"] else "supercritical"))
        for f in frac_pre:
            lev = 1.0 + f * (self.Tp - 1.0)
            k = np.nonzero(self.Ta >= lev)[0]
            out["Lpre_%g" % f] = float(-self.xa[k[0]]) if k.size else 0.0
        # Zel'dovich spike: excess above T1 behind the shock
        exc = self.Tb - self.T1
        if exc.max() > 1e-6 * self.T1:
            half = 0.5 * exc.max()
            k = np.nonzero(exc >= half)[0]
            out["spike_fwhm"] = float(self.xb[k[-1]] - self.xb[k[0]])
            pos = np.clip(exc, 0.0, None)
            out["spike_int"] = float(np.trapz(pos, self.xb))
        else:
            out["spike_fwhm"] = 0.0
            out["spike_int"] = 0.0
        out["cons"] = self.conservation()
        out["cont_gap"] = [float(v) for v in self.cont_gap]
        out["nroots"] = self.nroots
        out["Ms_check"] = float(self.Ms_check)
        out["ncand"] = self.ncand
        return out


class Units(object):
    def __init__(self, p):
        self.v = np.sqrt(p.g * KB * T0 / MH)
        self.L = (CL / self.v) / (3.0 * p.kap * SIG_CGS)
        self.p0 = AR * T0 ** 4 / (RHO0 * self.v ** 2)


def write_ref(path, sol, nref, frac=1e-4, dpad=0.25):
    """Uniform table over [xmin, xmax]; the downstream side is extended with the
    (constant) far-field state to at least dpad x the precursor length, so that the
    Dirichlet x1max face sits well behind the relaxation region (for M0 >= 5 the
    relaxation region alone is < 1e-3 of the domain, i.e. sub-cell)."""
    p = sol.p
    u = Units(p)
    xmin, xmax = sol.domain(frac)
    xmax = max(1.5 * xmax, -dpad * xmin)
    xg = np.linspace(xmin, xmax, nref)
    # put x = 0 (the shock) exactly on a grid point so the jump is not smeared
    xg = xg - xg[np.argmin(np.abs(xg))]
    s = sol.sample(xg)
    cols = np.column_stack([xg * u.L, s[0] * RHO0, s[1] * u.v, s[2] * T0, s[3] * T0])
    head = ("Lowrie-Edwards steady radiative shock (le08_shock.py, idort-rs-1017)\n"
            "M0=%.10g P0=%.10g gamma=%.10g sigma_a=%.10g kappa=%.10g\n"
            "scales: rho0=%.10e T0=%.10e v0=%.10e L=%.10e\n"
            "shock at x=%.10e ; columns: x rho v T_gas T_rad"
            % (p.m0, p.p0, p.g, p.sa, p.kap, RHO0, T0, u.v, u.L, -xg[0] * u.L))
    np.savetxt(path, cols, header=head, fmt="%.12e")
    return xg[0] * u.L, xg[-1] * u.L


def athinput(sol, xlo, xhi):
    p = sol.p
    u = Units(p)
    return "\n".join([
        "m1_shock_rho_l = %.8e" % RHO0, "m1_shock_v_l = %.8e" % (p.m0 * u.v),
        "m1_shock_t_l = %.8e" % T0, "m1_shock_rho_r = %.8e" % (sol.rho1 * RHO0),
        "m1_shock_v_r = %.8e" % (sol.v1 * u.v), "m1_shock_t_r = %.8e" % (sol.T1 * T0),
        "m1_shock_xs = %.8e" % (-xlo), "x1max = %.8e" % (xhi - xlo),
        "ebath_l = %.8e" % (AR * T0 ** 4), "ebath_r = %.8e" % (AR * (sol.T1 * T0) ** 4)])


# ---------------------------------------------------------------------------------
def _cross(xx, q, lev):
    k = np.nonzero(q >= lev)[0][0]
    return xx[k - 1] + (lev - q[k - 1]) / (q[k] - q[k - 1]) * (xx[k] - xx[k - 1])


def validate(args):
    """VALIDATE (nominal LE08 parameters P0 = 1e-4, sigma_a = 1e6, kappa = 1):
    (a) downstream state vs the LE08 closed-form T1 relation (as transcribed in QUOKKA
        extern/LowrieEdwards/radshock.py shock_jump);
    (b) profile features vs t7_radshock.py (independent derivation; it connects the
        branches in the (Theta, F_r) plane) -- run it separately (it fails for 1.05 and
        hangs for M0 >= 27), numbers in the idort-rs-1017 NOTE;
    (c) Mach 30 vs QUOKKA's LE08 reference table (extern/LowrieEdwards/
        shock_Mach_30.0.txt, CVODE, independent code): x of the level crossings
        T, T_r = 1.5, 10, 30, 50, 60 (--quokka DIR);
    (d) Mach 3 vs the published-figure digitisation shipped with QUOKKA
        (Tmat_fig30.csv / Trad_fig30.csv): peak T_gas, peak T_rad;
    (e) LE08 regime: continuous for 1.05, sub-critical 1.2/2/3, super-critical 5/27,
        continuous (no embedded hydro shock) for 30 and 50; the embedded shock
        vanishes between M0 = 27 and 28 (Mp -> 1)."""
    rows = []
    for m0 in (1.05, 1.2, 2.0, 3.0, 5.0, 27.0, 30.0, 50.0):
        p = LE08(m0)
        s = Solution(p)
        mt = s.metrics()
        r1, t1, v1 = p.downstream_published()
        rows.append((m0, mt, abs(s.rho1 / r1 - 1), abs(s.T1 / t1 - 1)))
        print("M0=%5.2f rho1=%.6f T1=%.6f | LE08 closed form rho1=%.6f T1=%.6f (rel %.1e %.1e) | "
              "%-13s Mp=%.4f Tp=%.4f Ts=%.4f Tmax=%.4f Trs=%.4f Lpre0.1=%.4e spike_fwhm=%.3e "
              "int=%.3e cons=%.1e/%.1e Mscheck=%.1e"
              % (m0, s.rho1, s.T1, r1, t1, rows[-1][2], rows[-1][3], mt["kind"], mt["Mp"],
                 mt["Tp"], mt["Ts"], mt["Tmax"], mt["Trs"], mt["Lpre_0.1"], mt["spike_fwhm"],
                 mt["spike_int"], mt["cons"][0], mt["cons"][1], mt["Ms_check"]))
        sys.stdout.flush()
        if args.quokka and m0 == 30.0:
            q = np.loadtxt(args.quokka + "/shock_Mach_30.0.txt")
            xg = np.linspace(q[0, 0], q[-1, 0], 400001)
            ref = s.sample(xg)
            for lev in (1.5, 10.0, 30.0, 50.0, 60.0):
                print("   QUOKKA M30 level %4g: x(T_r) %.6f vs %.6f  x(T) %.6f vs %.6f"
                      % (lev, _cross(q[:, 0], q[:, 4], lev), _cross(xg, ref[3], lev),
                         _cross(q[:, 0], q[:, 3], lev), _cross(xg, ref[2], lev)))
            print("   QUOKKA M30 max T %.5f vs %.5f" % (q[:, 3].max(), mt["Tmax"]))
        if args.quokka and m0 == 3.0:
            f1 = np.loadtxt(args.quokka + "/Tmat_fig30.csv", delimiter=",", skiprows=1)
            f2 = np.loadtxt(args.quokka + "/Trad_fig30.csv", delimiter=",", skiprows=1)
            print("   digitised M3 figure: max T_gas %.4f (le08 %.4f), max T_rad %.4f (le08 "
                  "T_rad at the shock %.4f, T1 %.4f)"
                  % (f1[:, 1].max(), mt["Tmax"], f2[:, 1].max(), mt["Trs"], s.T1))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--m0", type=float, default=2.0)
    ap.add_argument("--p0", type=float, default=1e-4)
    ap.add_argument("--gamma", type=float, default=5.0 / 3.0)
    ap.add_argument("--sigma-a", type=float, default=1e6)
    ap.add_argument("--kappa", type=float, default=1.0)
    ap.add_argument("--eps", type=float, default=1e-7)
    ap.add_argument("--write-ref", default=None)
    ap.add_argument("--nref", type=int, default=65537)
    ap.add_argument("--frac", type=float, default=1e-4, help="domain edge |T-Teq|/(T1-1)")
    ap.add_argument("--athinput", action="store_true")
    ap.add_argument("--dpad", type=float, default=0.25,
                    help="downstream extent >= dpad x precursor length")
    ap.add_argument("--json", default=None)
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--quokka", default=None,
                    help="dir with QUOKKA extern/LowrieEdwards files (validate c, d)")
    ap.add_argument("--nominal", action="store_true",
                    help="LE08 nominal P0/sigma_a/kappa (default: the deck-implied values)")
    a = ap.parse_args()
    if a.validate:
        validate(a)
        return 0
    if not a.nominal:
        a.p0, a.sigma_a, a.kappa = implied_params(a.gamma)
    p = LE08(a.m0, a.p0, a.gamma, a.sigma_a, a.kappa)
    s = Solution(p, a.eps)
    mt = s.metrics()
    u = Units(p)
    mt["units"] = {"L_cm": u.L, "v0": u.v, "P0_implied": u.p0}
    mt["domain"] = list(s.domain(a.frac))
    print(json.dumps(mt, indent=1))
    if a.write_ref:
        lo, hi = write_ref(a.write_ref, s, a.nref, a.frac, a.dpad)
        if a.athinput:
            print(athinput(s, lo, hi))
    if a.json:
        json.dump(mt, open(a.json, "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
