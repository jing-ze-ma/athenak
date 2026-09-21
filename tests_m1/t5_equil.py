#!/usr/bin/env python3
"""T5 (equilibration): single-zone matter-radiation coupling.

Design note sections 1 and 4(a).  With v = 0 and no transport the cell obeys

    dE/dt       = -chat * G0,      G0 = rho (kappa_E E - kappa_P a T^4)
    d(rho e)/dt =  c    * G0

so that ``E_gas + (c/chat) E`` is conserved exactly.  The reference is the
stiff integration of that pair (reduced to one ODE using the invariant) plus
the exact fixed point from a root find, for a pluggable EOS e(rho, T).

Pass (design note section 8, T5): the final state matches the exact
equilibrium to --eq-tol (1e-10 relative), the invariant is conserved to
--cons-tol (1e-13), and (T_gas - T_rad) never changes sign.
"""

import numpy as np

import common

K_B = 1.380649e-16
M_U = 1.66053906660e-24

try:
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    HAVE_SCIPY = True
except ImportError:                                    # pragma: no cover
    HAVE_SCIPY = False


class IdealEOS(object):
    """e(rho, T) = c_v T.  Replace this object to plug in a table."""

    def __init__(self, gamma=5.0 / 3.0, mu=0.6, cv=None):
        self.cv = cv if cv is not None else K_B / ((gamma - 1.0) * mu * M_U)
        self.name = "ideal(c_v=%.6e)" % self.cv

    def e_of_T(self, rho, temp):
        return self.cv * temp

    def T_of_e(self, rho, eint):
        return eint / self.cv


def bisect(fun, lo, hi, tol=1e-14, itmax=300):
    """Plain bisection fallback when scipy is unavailable."""
    flo, fhi = fun(lo), fun(hi)
    if flo * fhi > 0.0:
        raise ValueError("root not bracketed")
    for _ in range(itmax):
        mid = 0.5 * (lo + hi)
        fm = fun(mid)
        if flo * fm <= 0.0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
        if hi - lo <= tol * max(abs(hi), 1.0):
            break
    return 0.5 * (lo + hi)


def solve_root(fun, lo, hi):
    if HAVE_SCIPY:
        return brentq(fun, lo, hi, xtol=1e-14, rtol=1e-15, maxiter=300)
    return bisect(fun, lo, hi)


class Zone(object):
    def __init__(self, args, eos):
        self.rho = args.rho
        self.c = args.c
        self.chat = args.chat if args.chat is not None else args.c
        self.a_r = args.a_rad
        self.kP = args.kappa_p
        self.kE = args.kappa_e if args.kappa_e is not None else args.kappa_p
        self.eos = eos

    def temp(self, egas):
        return self.eos.T_of_e(self.rho, egas / self.rho)

    def t_rad(self, erad):
        return (max(erad, 0.0) / self.a_r) ** 0.25

    def invariant(self, egas, erad):
        return egas + (self.c / self.chat) * erad

    def rhs(self, erad, total):
        egas = total - (self.c / self.chat) * erad
        temp = self.temp(max(egas, 1e-300))
        g0 = self.rho * (self.kE * erad - self.kP * self.a_r * temp ** 4)
        return -self.chat * g0

    def equilibrium(self, total):
        """Exact fixed point: kappa_E E = kappa_P a T^4 at fixed invariant."""
        def resid(temp):
            erad = (self.kP / self.kE) * self.a_r * temp ** 4
            egas = self.rho * self.eos.e_of_T(self.rho, temp)
            return self.invariant(egas, erad) - total
        lo, hi = 1e-6, 1.0
        while resid(hi) < 0.0 and hi < 1e30:
            hi *= 10.0
        while resid(lo) > 0.0 and lo > 1e-30:
            lo *= 0.1
        temp = solve_root(resid, lo, hi)
        erad = (self.kP / self.kE) * self.a_r * temp ** 4
        egas = total - (self.c / self.chat) * erad
        return erad, egas, temp

    def integrate(self, erad0, total, times):
        """Stiff integration of dE/dt at fixed invariant."""
        times = np.asarray(times, dtype=float)
        if HAVE_SCIPY:
            sol = solve_ivp(lambda t, y: [self.rhs(y[0], total)],
                            (times[0], times[-1]), [erad0], method="BDF",
                            t_eval=times, rtol=1e-12, atol=1e-30 * abs(total))
            return np.asarray(sol.y[0])
        out = [erad0]
        erad = erad0
        for t0, t1 in zip(times[:-1], times[1:]):
            nsub = 200
            dt = (t1 - t0) / nsub
            for _ in range(nsub):
                def res(en, e_old=erad, step=dt):
                    return en - e_old - step * self.rhs(en, total)
                lo = min(erad, 0.0)
                hi = max(erad, total * self.chat / self.c) + abs(erad) + 1.0
                erad = solve_root(res, lo, hi)
            out.append(erad)
        return np.array(out)


def load_code_series(args):
    if args.hist:
        names, arr = common.read_history(args.hist)
        t = common.hist_column(names, arr, args.hist_time)
        erad = common.hist_column(names, arr, args.hist_erad)
        egas = common.hist_column(names, arr, args.hist_egas)
        return t, erad, egas
    dumps = common.load_series(args.dumps, args.bin_convert_dir,
                               args.all_ranks)
    if args.hydro:
        hyd = common.load_series(args.hydro, args.bin_convert_dir,
                                 args.all_ranks)
        for d in dumps:
            best = min(hyd, key=lambda h: abs(h.time - d.time))
            for k, v in best.data.items():
                d.data.setdefault(k, v)
    t, erad, egas = [], [], []
    for d in dumps:
        t.append(d.time)
        erad.append(float(d.var("m1_e").mean()))
        if d.has("eint"):
            egas.append(float(d.var("eint").mean()))
        elif d.has("ener"):
            # hydro_u output: the TOTAL energy density; these zones are at rest
            egas.append(float(d.var("ener").mean()))
        else:
            egas.append(float(d.var("press").mean()) / (args.gamma - 1.0))
    return np.array(t), np.array(erad), np.array(egas)


def selftest(args, zone):
    """Synthetic 'code' output = the reference trajectory itself."""
    tcoup = 1.0 / (zone.chat * zone.rho * zone.kE)
    times = np.linspace(0.0, 200.0 * tcoup, 12)
    total = zone.invariant(args.egas0, args.erad0)
    erad = zone.integrate(args.erad0, total, times)
    egas = total - (zone.c / zone.chat) * erad
    if args.selftest_fail:
        erad = erad * (1.0 + 1e-6)
    return times, erad, egas


def main():
    p = common.base_parser(__doc__.splitlines()[0])
    p.add_argument("dumps", nargs="*", help="single-zone dumps (time series)")
    p.add_argument("--hydro", nargs="*", default=[],
                   help="parallel series of hydro dumps (file_type = tab "
                        "writes one file per <output> block)")
    p.add_argument("--hist", default=None,
                   help="history file instead of dumps")
    p.add_argument("--hist-time", default="time")
    p.add_argument("--hist-erad", default="m1_e",
                   help="history column (name or 1-based index) for E")
    p.add_argument("--hist-egas", default="eint",
                   help="history column for the gas energy density")
    p.add_argument("--rho", type=float, default=1.0)
    p.add_argument("--c", type=float, default=common.C_CGS,
                   help="speed of light in code units")
    p.add_argument("--chat", type=float, default=None,
                   help="reduced speed of light (default: = c)")
    p.add_argument("--a-rad", type=float, default=common.AR_CGS,
                   help="radiation constant a_r in code units")
    p.add_argument("--kappa-p", type=float, default=4e-8,
                   help="Planck (emission) opacity per unit mass")
    p.add_argument("--kappa-e", type=float, default=None,
                   help="energy (absorption) mean; default: = kappa_P")
    p.add_argument("--gamma", type=float, default=5.0 / 3.0)
    p.add_argument("--mu", type=float, default=0.6,
                   help="mean molecular weight of the ideal EOS")
    p.add_argument("--cv", type=float, default=None,
                   help="specific heat at constant volume; overrides --mu")
    p.add_argument("--erad0", type=float, default=1e12,
                   help="initial E (HERACLES/T5 value)")
    p.add_argument("--egas0", type=float, default=1e2,
                   help="initial gas energy density (1e2 or 1e10 in T5)")
    p.add_argument("--eq-tol", type=float, default=1e-10,
                   help="relative tolerance on the final equilibrium")
    p.add_argument("--cons-tol", type=float, default=1e-13,
                   help="relative tolerance on E_gas + (c/chat) E")
    p.add_argument("--ode-tol", type=float, default=None,
                   help="if given, also require the trajectory to match the "
                        "stiff reference to this relative tolerance")
    p.add_argument("--sign-eps", type=float, default=1e-6,
                   help="|T_gas - T_rad| below this fraction of T counts as "
                        "converged, not as a sign flip")
    p.add_argument("--no-equilibrium", action="store_true",
                   help="skip the final-equilibrium check (short runs)")
    p.add_argument("--selftest-fail", action="store_true",
                   help="selftest variant that must FAIL")
    args = p.parse_args()

    eos = IdealEOS(args.gamma, args.mu, args.cv)
    zone = Zone(args, eos)

    if args.selftest:
        t, erad, egas = selftest(args, zone)
    else:
        if not args.dumps and not args.hist:
            p.error("dumps or --hist are required unless --selftest")
        t, erad, egas = load_code_series(args)

    inv = zone.invariant(egas, erad)
    cons = float(np.abs(inv / inv[0] - 1.0).max())
    total = float(inv[0])

    eq_e, eq_g, eq_t = zone.equilibrium(total)
    err_e = abs(erad[-1] / eq_e - 1.0)
    err_g = abs(egas[-1] / eq_g - 1.0)

    tgas = np.array([zone.temp(g) for g in egas])
    trad = np.array([zone.t_rad(e) for e in erad])
    d = tgas - trad
    nz = d[np.abs(d) > args.sign_eps * np.maximum(tgas, trad)]
    nsign = int((np.diff(np.sign(nz)) != 0).sum()) if nz.size > 1 else 0

    ode_err = float("nan")
    if args.ode_tol is not None or not args.quiet:
        ref = zone.integrate(erad[0], total, t)
        ode_err = float(np.abs(erad / ref - 1.0).max())

    common.report(args, "EOS=%s  scipy=%s  T_eq=%.10e  E_eq=%.10e  "
                        "Egas_eq=%.10e" % (eos.name, HAVE_SCIPY, eq_t, eq_e,
                                           eq_g))
    common.report(args, "final: E=%.10e  Egas=%.10e  t=%.6e"
                  % (erad[-1], egas[-1], t[-1]))
    ok = cons <= args.cons_tol and nsign == 0
    if not args.no_equilibrium:
        ok = ok and max(err_e, err_g) <= args.eq_tol
    if args.ode_tol is not None:
        ok = ok and ode_err <= args.ode_tol
    common.verdict(ok, "T5 equil: dE_eq=%.3e dEgas_eq=%.3e (<= %.3g)  "
                       "conservation=%.3e (<= %.3g)  sign flips=%d  "
                       "ODE dev=%.3e"
                   % (err_e, err_g, args.eq_tol, cons, args.cons_tol,
                      nsign, ode_err))


if __name__ == "__main__":
    main()
