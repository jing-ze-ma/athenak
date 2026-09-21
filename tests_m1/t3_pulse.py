#!/usr/bin/env python3
"""T3 (static thick pulse): the asymptotic-preserving gate.

Design note section 8, (T3): 1-D Gaussian in E, scattering dominated, no
hydro.  The variance of E must grow at the physical rate

    d(sigma^2)/dt = 2 D,   D = c / (3 rho kappa)      (1-D diffusion)

Pass: within --tol (default 2 %) for ``thick_flux = ap_hll``.  For
``thick_flux = scaled`` the predicted ratio is 1 + X with
X = 3<|mu|>/(2 * scaled_prefactor) (4 % at the default prefactor 20): pass
``--expect-excess 0.04``.  ``thick_flux = none`` is the failing control at
tau_cell >= 10.

With ``--nyquist`` the script instead measures the decay rate of a seeded
odd-even mode and compares it with the rate of the compact 3-point diffusion
operator at the Nyquist wavenumber, 4 D / dx^2.
"""

import numpy as np

import common


def series_variance(dumps, args):
    t, s2 = [], []
    for d in dumps:
        x, e = common.extract_1d(d, "m1_e", axis=args.axis,
                                 reduce=args.reduce)
        bg = e.min() if args.subtract_min else args.background
        _, var = common.moments(x, e, background=bg)
        t.append(d.time)
        s2.append(var)
    return np.array(t), np.array(s2)


def series_nyquist(dumps, args):
    t, amp = [], []
    for d in dumps:
        _, e = common.extract_1d(d, "m1_e", axis=args.axis,
                                 reduce=args.reduce)
        t.append(d.time)
        amp.append(abs(common.nyquist_amplitude(e)))
    return np.array(t), np.array(amp)


def selftest(args):
    """Synthetic diffusing Gaussian (or decaying odd-even mode)."""
    n = 256
    x = common.centres(0.0, 1.0, n)
    dx = x[1] - x[0]
    diff = args.c / (3.0 * args.rho * args.kappa)
    fac = 1.10 if args.selftest_fail else 1.0
    sig0 = 0.05
    gam = 4.0 * diff / dx ** 2
    tend = 3.0 / gam if args.nyquist else 0.25 * sig0 ** 2 / diff
    dumps = []
    for k in range(6):
        tt = tend * k / 5.0
        if args.nyquist:
            sgn = np.where(np.arange(n) % 2 == 0, 1.0, -1.0)
            e = 1.0 + 1e-3 * np.exp(-gam * fac * tt) * sgn
        else:
            var = sig0 ** 2 + 2.0 * diff * tt * fac
            e = np.exp(-0.5 * (x - 0.5) ** 2 / var) / np.sqrt(var)
        dumps.append(common.Dump(tt, x, np.array([0.5]), np.array([0.5]),
                                 {"m1_e": e[None, None, :]}, "<selftest>"))
    return dumps


def main():
    p = common.base_parser(__doc__.splitlines()[0])
    p.add_argument("dumps", nargs="*", help="time series of bin dumps")
    p.add_argument("--c", type=float, default=1.0,
                   help="speed of light in code units")
    p.add_argument("--rho", type=float, default=1.0,
                   help="background density (uniform)")
    p.add_argument("--kappa", type=float, default=1.0,
                   help="scattering opacity kappa_s (per unit mass)")
    p.add_argument("--tol", type=float, default=0.02,
                   help="allowed relative deviation of the measured ratio")
    p.add_argument("--expect-excess", type=float, default=0.0,
                   help="predicted excess X; target ratio is 1 + X")
    p.add_argument("--skip", type=int, default=0,
                   help="drop this many early dumps from the fit")
    p.add_argument("--background", type=float, default=0.0,
                   help="constant background subtracted from E")
    p.add_argument("--subtract-min", action="store_true",
                   help="subtract min(E) of each dump instead")
    p.add_argument("--nyquist", action="store_true",
                   help="measure the odd-even decay rate instead")
    p.add_argument("--axis", type=int, default=1, choices=(1, 2, 3))
    p.add_argument("--reduce", default="mid", choices=("mid", "mean"))
    p.add_argument("--selftest-fail", action="store_true",
                   help="selftest variant that must FAIL")
    args = p.parse_args()

    if args.selftest:
        dumps = selftest(args)
    else:
        if len(args.dumps) < 3:
            p.error("at least three dumps are needed for the fit")
        dumps = common.load_series(args.dumps, args.bin_convert_dir,
                                   args.all_ranks)
    dumps = dumps[args.skip:]
    diff = args.c / (3.0 * args.rho * args.kappa)
    target = 1.0 + args.expect_excess

    if args.nyquist:
        t, amp = series_nyquist(dumps, args)
        if (amp <= 0).any():
            common.verdict(False, "T3 nyquist: non-positive mode amplitude")
        slope, _ = common.linfit(t, np.log(amp))
        rate = -slope
        expect = 4.0 * diff / dumps[0].dx1 ** 2
        ratio = rate / expect
        for tt, aa in zip(t, amp):
            common.report(args, "t=%.6e  |E_nyq|=%.6e" % (tt, aa))
        ok = abs(ratio / target - 1.0) <= args.tol
        common.verdict(ok, "T3 nyquist: rate=%.6e  4D/dx^2=%.6e  "
                           "ratio=%.6f (target %.6f +- %.3g)"
                       % (rate, expect, ratio, target, args.tol))

    t, s2 = series_variance(dumps, args)
    slope, _ = common.linfit(t, s2)
    ratio = slope / (2.0 * diff)
    for tt, vv in zip(t, s2):
        common.report(args, "t=%.6e  sigma^2=%.6e" % (tt, vv))
    ok = abs(ratio / target - 1.0) <= args.tol
    common.verdict(ok, "T3 pulse: d(sigma^2)/dt=%.6e  2D=%.6e  "
                       "ratio=%.6f (target %.6f +- %.3g)"
                   % (slope, 2.0 * diff, ratio, target, args.tol))


if __name__ == "__main__":
    main()
