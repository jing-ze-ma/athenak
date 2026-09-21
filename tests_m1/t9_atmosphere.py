#!/usr/bin/env python3
"""t9_atmosphere.py -- gate I6b of milestone 3a2: the STEADY grey plane-parallel column.

The run is `<problem>/m1_test = atmosphere` (src/pgen/tests/rad_m1_tests2.cpp,
inputs/tests/rad_m1_atmosphere.athinput): constant imposed flux F at the bottom face,
a free surface at the top, a prescribed exponential density and a pure-scattering grey
opacity, so the converged state is the closed-form M1 moment solution

    chi(f)/f = q0 + tau,   E(tau) = F/(c f(tau)),   P_rad = (F/c)(q0 + tau),

with tau measured DOWNWARD from the top mesh face and

    tau(z) = kappa rho_top H (exp[(ztop - z)/H] - 1).

The CONSISTENT PAIR.  The Marshak condition at the surface is F = c q E_top, i.e.
f(0) = q = <rad_m1>/marshak_q exactly; hence q0 = chi(q)/q, and for q = 1/2,
chi(1/2) = 0.46481586 and q0 = 0.92963172.  Nothing else in the problem fixes q0.

chi(f)/f falls from +inf at f -> 0 to a minimum 0.89806 at f = 0.82504 and rises back
to 1 at f = 1.  The physical branch is the decreasing one (f falls with depth); q0 for
any q <= 0.8 is above the minimum, so the solution never leaves it.

Usage:
    python3 t9_atmosphere.py <dir>/tab/*.m1.*.tab --kappa 1 --rho-top 0.128 \\
        --scale-h 0.113 --flux 1 --c 1 --ztop 1 [--q 0.5] [--json out.json] [--quiet]
    python3 t9_atmosphere.py --selftest [--selftest-fail]

Verdict: max relative error of E(tau) and of f(tau) against the closed form, over the
whole column and over the top 5 cells separately (--tol, default 0.01).
"""

import argparse
import json
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import common  # noqa: E402


def m1_chi(f):
    """Levermore closure factor, the same expression as rad_m1_closure.hpp M1Chi."""
    f2 = np.minimum(np.asarray(f, dtype=float)**2, 1.0)
    return (3.0 + 4.0*f2)/(5.0 + 2.0*np.sqrt(np.maximum(4.0 - 3.0*f2, 0.0)))


def f_from_chiover(q):
    """Invert chi(f)/f = q on the DECREASING branch f in (0, f*], by bisection."""
    q = np.asarray(q, dtype=float)
    lo = np.full(q.shape, 1.0e-14)
    hi = np.full(q.shape, 0.82504)
    for _ in range(200):
        mid = 0.5*(lo + hi)
        big = (m1_chi(mid)/mid) > q
        lo = np.where(big, mid, lo)
        hi = np.where(big, hi, mid)
    return 0.5*(lo + hi)


def exact(tau, q0, flux, clight):
    f = f_from_chiover(q0 + tau)
    return flux/(clight*f), f


def analyse(x1, e, f1, args):
    tau = args.kappa*args.rho_top*args.scale_h*np.expm1((args.ztop - x1)/args.scale_h)
    q0 = float(m1_chi(args.q)/args.q)
    eex, fex = exact(tau, q0, args.flux, args.c)
    fnum = f1/(args.c*e)
    ee = np.abs(e/eex - 1.0)
    ef = np.abs(fnum/fex - 1.0)
    # the top 5 cells are the ones the gate exists for: tau there is ~1e-3
    top = np.argsort(tau)[:5]
    return dict(tau=tau, q0=q0, e=e, e_exact=eex, f=fnum, f_exact=fex,
                max_e=float(ee.max()), max_f=float(ef.max()),
                top_e=float(ee[top].max()), top_f=float(ef[top].max()),
                l1_e=float(np.mean(ee)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='*')
    ap.add_argument('--c', type=float, default=1.0)
    ap.add_argument('--kappa', type=float, default=1.0)
    ap.add_argument('--rho-top', type=float, default=0.128)
    ap.add_argument('--scale-h', type=float, default=0.113)
    ap.add_argument('--flux', type=float, default=1.0)
    ap.add_argument('--ztop', type=float, default=1.0)
    ap.add_argument('--q', type=float, default=0.5)
    ap.add_argument('--tol', type=float, default=0.01)
    ap.add_argument('--json', default=None)
    ap.add_argument('--label', default='run')
    ap.add_argument('--quiet', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--selftest-fail', action='store_true')
    args = ap.parse_args()

    if args.selftest:
        # build the exact solution on the default grid and feed it back in
        n = 128
        x1 = (np.arange(n) + 0.5)/n*args.ztop
        tau = args.kappa*args.rho_top*args.scale_h*np.expm1(
            (args.ztop - x1)/args.scale_h)
        q0 = float(m1_chi(args.q)/args.q)
        e, f = exact(tau, q0, args.flux, args.c)
        if args.selftest_fail:
            e = e*(1.0 + 0.05*np.sin(np.arange(n)))
        r = analyse(x1, e, args.c*e*f, args)
        ok = (r['max_e'] <= args.tol) and (r['max_f'] <= args.tol)
        print('%s T9 atmosphere selftest: max|dE/E|=%.3e max|df/f|=%.3e (<= %.3g)'
              % ('PASS' if ok else 'FAIL', r['max_e'], r['max_f'], args.tol))
        return 0 if ok else 1

    if not args.files:
        print('FAIL T9 atmosphere: no dump given')
        return 1
    d = common.load_dump(sorted(args.files)[-1])
    x1, e = common.extract_1d(d, 'm1_e', axis=1)
    _, f1 = common.extract_1d(d, 'm1_f1', axis=1)
    r = analyse(np.asarray(x1), np.asarray(e), np.asarray(f1), args)
    ok = (r['max_e'] <= args.tol) and (r['max_f'] <= args.tol)
    if args.json:
        with open(args.json, 'w') as fp:
            json.dump({k: (v.tolist() if isinstance(v, np.ndarray) else v)
                       for k, v in r.items()} | {'label': args.label}, fp)
    print('%s T9 atmosphere: q0=%.8f max|dE/E|=%.3e max|df/f|=%.3e  '
          'top5: dE=%.3e df=%.3e  L1(E)=%.3e (tol %.3g)'
          % ('PASS' if ok else 'FAIL', r['q0'], r['max_e'], r['max_f'],
             r['top_e'], r['top_f'], r['l1_e'], args.tol))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
