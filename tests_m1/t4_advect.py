#!/usr/bin/env python3
"""T4 (advected thick pulse): dynamic diffusion.

Design note section 8, (T4): periodic box, prescribed v.  Pass: the pulse
centre sits at x0 + v t within one cell, the width equals the static solution
sqrt(sigma0^2 + 2 D t) within 2 %, and there is no net gas heating beyond
1e-10 relative (the test that the beta^2 terms of section 1 cancel).  D is
c/(3 rho kappa).  An odd-even diagnostic (power in the Nyquist mode of E
relative to the total fluctuation power) is reported, and enters the verdict
only if --nyq-tol is given: QUOKKA needed 1024 cells to suppress it, we
require none at 512 with ap_hll.

Moments are taken with periodic (circular) statistics, so a pulse that has
wrapped around the box is handled correctly.
"""

import numpy as np

import common


def circular_moments(x, w, length, background=0.0):
    """Circular mean and variance of a periodic profile."""
    wp = np.maximum(np.asarray(w, dtype=float) - background, 0.0)
    tot = wp.sum()
    if tot <= 0.0:
        raise ValueError("no positive weight left after background removal")
    th = 2.0 * np.pi * (x - x[0]) / length
    mean_th = np.arctan2((wp * np.sin(th)).sum(), (wp * np.cos(th)).sum())
    mean = x[0] + length * (mean_th % (2.0 * np.pi)) / (2.0 * np.pi)
    d = (x - mean + 0.5 * length) % length - 0.5 * length
    var = float((wp * d ** 2).sum() / tot)
    return float(mean), var


def gas_energy(dump, gamma):
    if dump.has("eint"):
        return float(dump.var("eint").sum())
    if dump.has("press"):
        return float(dump.var("press").sum()) / (gamma - 1.0)
    return None


def gas_temperature(dump, gamma):
    """Volume mean of the code temperature T = P/rho of an ideal gas."""
    if not dump.has("dens"):
        return None
    rho = dump.var("dens")
    if dump.has("eint"):
        tt = dump.var("eint") * (gamma - 1.0) / rho
    elif dump.has("press"):
        tt = dump.var("press") / rho
    else:
        return None
    return float(tt.mean())


def merge_hydro(dumps, hydro):
    """Fold the variables of a parallel hydro series into the m1 dumps.

    ``file_type = tab`` writes one file per output block, so the radiation
    and the hydro variables of the same snapshot are in two files; they are
    matched by time.
    """
    for d in dumps:
        best = min(hydro, key=lambda h: abs(h.time - d.time))
        if abs(best.time - d.time) > 1e-8 * max(abs(d.time), 1.0):
            raise ValueError("no hydro dump at t = %.17g" % d.time)
        for k, v in best.data.items():
            d.data.setdefault(k, v)
    return dumps


def measure(dump, args):
    x, e = common.extract_1d(dump, "m1_e", axis=args.axis,
                             reduce=args.reduce)
    dx = x[1] - x[0]
    length = x.size * dx
    bg = e.min() if args.subtract_min else args.background
    mean, var = circular_moments(x, e, length, background=bg)
    return dict(t=dump.time, x=x, e=e, dx=dx, length=length, centre=mean,
                sigma=np.sqrt(var), nyq=common.nyquist_power_fraction(e),
                egas=gas_energy(dump, args.gamma),
                tgas=gas_temperature(dump, args.gamma))


def selftest(args):
    n = 512
    x = common.centres(0.0, 1.0, n)
    diff = args.c / (3.0 * args.rho * args.kappa)
    sig0 = 0.05
    tend = 0.2 * sig0 ** 2 / diff
    dumps = []
    for k in range(4):
        tt = tend * k / 3.0
        var = sig0 ** 2 + 2.0 * diff * tt
        cen = 0.3 + args.v * tt
        egas = 1.0
        if args.selftest_fail:
            cen += 3.0 * (x[1] - x[0]) * k / 3.0
            var *= 1.0 + 0.2 * k / 3.0
            egas = 1.0 + 1e-6 * k
        d = (x - cen + 0.5) % 1.0 - 0.5
        e = np.exp(-0.5 * d ** 2 / var) / np.sqrt(var)
        data = {"m1_e": e[None, None, :],
                "eint": np.full((1, 1, n), egas / n)}
        dumps.append(common.Dump(tt, x, np.array([0.5]), np.array([0.5]),
                                 data, "<selftest>"))
    return dumps


def main():
    p = common.base_parser(__doc__.splitlines()[0])
    p.add_argument("dumps", nargs="*", help="time series of bin dumps")
    p.add_argument("--hydro", nargs="*", default=[],
                   help="parallel series of hydro dumps (file_type = tab "
                        "writes one file per <output> block), matched by time")
    p.add_argument("--static", default=None,
                   help="final dump of the STATIC run at the same time; the "
                        "advected profile is compared with it shifted by v t "
                        "(QUOKKA report < 0.03 %% static, < 0.06 %% dynamic)")
    p.add_argument("--static-tol", type=float, default=None,
                   help="pass threshold on that relative L1 difference")
    p.add_argument("--temp-tol", type=float, default=None,
                   help="if given, also require the relative drift of the "
                        "mean gas temperature to stay below this (T4b)")
    p.add_argument("--c", type=float, default=1.0,
                   help="speed of light in code units")
    p.add_argument("--v", type=float, default=0.0,
                   help="prescribed advection velocity")
    p.add_argument("--rho", type=float, default=1.0)
    p.add_argument("--kappa", type=float, default=1.0,
                   help="opacity per unit mass (scattering + absorption)")
    p.add_argument("--x0", type=float, default=None,
                   help="initial centre (default: measured from dump 0)")
    p.add_argument("--sigma0", type=float, default=None,
                   help="initial width (default: measured from dump 0)")
    p.add_argument("--gamma", type=float, default=5.0 / 3.0,
                   help="gas gamma, used only if the dump carries 'press'")
    p.add_argument("--centre-tol", type=float, default=1.0,
                   help="allowed centre error in cells")
    p.add_argument("--width-tol", type=float, default=0.02,
                   help="allowed relative width error")
    p.add_argument("--drift-tol", type=float, default=1e-10,
                   help="allowed relative gas-energy drift")
    p.add_argument("--nyq-tol", type=float, default=None,
                   help="if given, also require the Nyquist power fraction "
                        "of E to stay below this")
    p.add_argument("--background", type=float, default=0.0)
    p.add_argument("--subtract-min", action="store_true")
    p.add_argument("--axis", type=int, default=1, choices=(1, 2, 3))
    p.add_argument("--reduce", default="mid", choices=("mid", "mean"))
    p.add_argument("--selftest-fail", action="store_true",
                   help="selftest variant that must FAIL")
    args = p.parse_args()

    if args.selftest:
        dumps = selftest(args)
    else:
        if len(args.dumps) < 2:
            p.error("at least two dumps are needed (initial and final)")
        dumps = common.load_series(args.dumps, args.bin_convert_dir,
                                   args.all_ranks)
        if args.hydro:
            dumps = merge_hydro(dumps, common.load_series(
                args.hydro, args.bin_convert_dir, args.all_ranks))

    m = [measure(d, args) for d in dumps]
    first, last = m[0], m[-1]
    x0 = args.x0 if args.x0 is not None else first["centre"]
    sig0 = args.sigma0 if args.sigma0 is not None else first["sigma"]
    diff = args.c / (3.0 * args.rho * args.kappa)
    dt = last["t"] - first["t"]

    expect_c = x0 + args.v * dt
    length = last["length"]
    dcen = (last["centre"] - expect_c + 0.5 * length) % length - 0.5 * length
    dcen_cells = dcen / last["dx"]
    expect_s = np.sqrt(sig0 ** 2 + 2.0 * diff * dt)
    dwidth = last["sigma"] / expect_s - 1.0

    drift = 0.0
    if first["egas"] is not None and first["egas"] != 0.0:
        drift = abs(last["egas"] / first["egas"] - 1.0)

    for mm in m:
        common.report(args, "t=%.6e  centre=%.6e  sigma=%.6e  nyq=%.3e"
                      % (mm["t"], mm["centre"], mm["sigma"], mm["nyq"]))
    sdiff = float("nan")
    if args.static:
        sd = common.load_dump(args.static, args.bin_convert_dir,
                              args.all_ranks)
        xs, es = common.extract_1d(sd, "m1_e", axis=args.axis,
                                   reduce=args.reduce)
        # the static profile evaluated at x - v t, periodically
        shift = args.v * (last["t"] - first["t"])
        xq = (last["x"] - shift - xs[0]) % length + xs[0]
        ref = np.interp(xq, np.append(xs, xs[0] + length),
                        np.append(es, es[0]))
        sdiff = common.l1_rel(last["e"], ref)

    tdrift = float("nan")
    if first["tgas"] is not None and first["tgas"] != 0.0:
        tdrift = abs(last["tgas"] / first["tgas"] - 1.0)

    ok = (abs(dcen_cells) <= args.centre_tol
          and abs(dwidth) <= args.width_tol
          and drift <= args.drift_tol)
    if args.nyq_tol is not None:
        ok = ok and last["nyq"] <= args.nyq_tol
    if args.temp_tol is not None:
        ok = ok and (tdrift == tdrift) and tdrift <= args.temp_tol
    if args.static_tol is not None:
        ok = ok and (sdiff == sdiff) and sdiff <= args.static_tol
    common.verdict(ok, "T4 advect: dcentre=%.3f cells (<= %.3g)  "
                       "dwidth=%.3e (<= %.3g)  gas drift=%.3e (<= %.3g)  "
                       "Tgas drift=%.3e  vs-static=%.3e  nyquist=%.3e"
                   % (dcen_cells, args.centre_tol, dwidth, args.width_tol,
                      drift, args.drift_tol, tdrift, sdiff, last["nyq"]))


if __name__ == "__main__":
    main()
