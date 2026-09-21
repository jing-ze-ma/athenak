#!/usr/bin/env python3
"""T10: the LINEAR RADIATION-MODIFIED ACOUSTIC WAVE (milestone 3b phase E).

One wavelength of a linear sound wave in a uniform, optically thick,
radiation-pressure-significant medium, run along x1, along x2, along x3 or on
the x1-x2 diagonal (``<problem>/radwave_dir``).  In the equilibrium-diffusion
limit the mixture has Chandrasekhar's generalised exponents, with
``beta = P_gas/P_tot``,

    Gamma_1   = beta + (4 - 3 beta)^2 (gamma-1)/[beta + 12 (gamma-1)(1-beta)]
    Gamma_3-1 = (Gamma_1 - beta)/(4 - 3 beta)
    c_s^2     = Gamma_1 (P_gas + P_rad)/rho

and the wave travels at ``c_s``, weakly damped by radiative diffusion.

The measurement is the complex Fourier amplitude of the density perturbation
at the fundamental mode of the box,

    Z(t) = < (rho - <rho>) exp(-i k.x) >   ->   (A/2) exp(-i omega t - Gamma t)

so that the unwrapped phase gives ``omega`` (hence the phase speed
``omega/|k|``) and ``ln|Z|`` gives the growth/decay rate.

PASS: the phase speed is within ``--tol`` (1 % by default) of ``c_s`` and the
mode does not grow (``Gamma >= -tol_growth/period``).

The same analysis path runs on synthetic analytic data with ``--selftest``.
"""

import glob
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402


def mixture(gamma, beta):
    """Chandrasekhar's Gamma_1 and Gamma_3 - 1 for a gas + radiation mix."""
    gm1 = gamma - 1.0
    g1 = beta + (4.0 - 3.0 * beta) ** 2 * gm1 / (beta + 12.0 * gm1
                                                 * (1.0 - beta))
    return g1, (g1 - beta) / (4.0 - 3.0 * beta)


def analytic(rho0, t0, arad, gamma):
    """(beta, Gamma_1, Gamma_3-1, c_s) of the uniform background."""
    pgas = rho0 * t0
    prad = arad * t0 ** 4 / 3.0
    ptot = pgas + prad
    beta = pgas / ptot
    g1, g3m1 = mixture(gamma, beta)
    return beta, g1, g3m1, math.sqrt(g1 * ptot / rho0)


def wave_vector(direction, lx, ly, lz):
    """The box fundamental wave vector for a named direction."""
    modes = {"x1": (1, 0, 0), "x2": (0, 1, 0), "x3": (0, 0, 1),
             "xy": (1, 1, 0)}
    if direction not in modes:
        raise ValueError("radwave_dir '%s' is not x1|x2|x3|xy" % direction)
    mx, my, mz = modes[direction]
    return (2.0 * math.pi * mx / lx, 2.0 * math.pi * my / ly,
            2.0 * math.pi * mz / lz)


def fourier_amplitude(dump, kvec, name="dens"):
    """Complex Fourier coefficient of the fundamental mode of a dump."""
    arr = np.asarray(dump.var(name), dtype=float)
    x3, x2, x1 = np.meshgrid(dump.x3, dump.x2, dump.x1, indexing="ij")
    phase = kvec[0] * x1 + kvec[1] * x2 + kvec[2] * x3
    fluc = arr - arr.mean()
    return complex((fluc * np.exp(-1j * phase)).mean())


def measure(dumps, kvec, name="dens"):
    """(times, omega, gamma, |Z|) from a series of dumps."""
    times = np.array([d.time for d in dumps], dtype=float)
    zz = np.array([fourier_amplitude(d, kvec, name) for d in dumps])
    mag = np.abs(zz)
    if np.any(mag <= 0.0):
        raise ValueError("the fundamental mode is empty in some dump")
    phase = np.unwrap(np.angle(zz))
    slope, _ = common.linfit(times, phase)
    omega = -slope
    gslope, _ = common.linfit(times, np.log(mag))
    return times, omega, gslope, mag


def synth_dumps(kvec, cs, gamma_rate, amp, nx, nper, ndump, lx, ly, lz):
    """Analytic dumps of the travelling mode, for --selftest."""
    kmag = math.sqrt(sum(k * k for k in kvec))
    omega = cs * kmag
    period = 2.0 * math.pi / omega
    n1 = nx if kvec[0] != 0.0 else 4
    n2 = nx if kvec[1] != 0.0 else 4
    n3 = nx if kvec[2] != 0.0 else 1
    x1 = common.centres(0.0, lx, n1)
    x2 = common.centres(0.0, ly, n2)
    x3 = common.centres(0.0, lz, n3)
    g3, g2, g1 = np.meshgrid(x3, x2, x1, indexing="ij")
    phase = kvec[0] * g1 + kvec[1] * g2 + kvec[2] * g3
    out = []
    for n in range(ndump):
        tt = nper * period * n / (ndump - 1)
        dens = 1.0 + amp * math.exp(gamma_rate * tt) * np.cos(phase
                                                              - omega * tt)
        out.append(common.Dump(tt, x1, x2, x3, {"dens": dens}, "synthetic"))
    return out


def main():
    p = common.base_parser(__doc__)
    p.add_argument("--arm", action="append", default=None,
                   metavar="PATH[,DIR[,LABEL[,PEER]]]",
                   help="one arm: the directory holding its bin/ dumps, its "
                        "radwave_dir (default x1), a label (default the "
                        "directory name) and what it is gated against -- "
                        "empty = the analytic c_s, a LABEL = that arm (to "
                        "--match-tol), 'none' = reported but not gated.  "
                        "Repeat for several arms.")
    p.add_argument("--match-tol", type=float, default=1.0e-4,
                   help="allowed relative difference against a peer arm")
    p.add_argument("--glob", default="bin/*.hydro_w.*.bin",
                   help="dump glob inside each rundir")
    p.add_argument("--rho0", type=float, default=1.0)
    p.add_argument("--t0", type=float, default=1.0)
    p.add_argument("--arad", type=float, default=3.0)
    p.add_argument("--gamma", type=float, default=5.0 / 3.0)
    p.add_argument("--lx", type=float, default=1.0)
    p.add_argument("--ly", type=float, default=1.0)
    p.add_argument("--lz", type=float, default=1.0)
    p.add_argument("--tol", type=float, default=1.0e-2,
                   help="allowed |c_phase/c_s - 1|")
    p.add_argument("--growth-tol", type=float, default=1.0e-3,
                   help="allowed growth per period")
    p.add_argument("--json", default=None, help="write the numbers here")
    args = p.parse_args()

    beta, g1, g3m1, cs = analytic(args.rho0, args.t0, args.arad, args.gamma)
    common.report(args, "background: beta=%.6f Gamma_1=%.9f Gamma_3-1=%.9f "
                        "c_s=%.12g" % (beta, g1, g3m1, cs))

    rows = []
    peers = {}
    if args.selftest:
        for name, rate in (("x1", 0.0), ("x2", -0.05), ("xy", 0.0)):
            kvec = wave_vector(name, args.lx, args.ly, args.lz)
            dumps = synth_dumps(kvec, cs, rate, 1.0e-4, 64, 2.0, 33,
                                args.lx, args.ly, args.lz)
            rows.append((name + "/synth", name, kvec, dumps))
    else:
        if not args.arm:
            print("FAIL no --arm given (and no --selftest)")
            return 1
        for spec in args.arm:
            parts = spec.split(",")
            rd = parts[0]
            dd = parts[1] if len(parts) > 1 and parts[1] else "x1"
            lb = (parts[2] if len(parts) > 2 and parts[2]
                  else os.path.basename(os.path.normpath(rd)))
            peer = parts[3] if len(parts) > 3 else ""
            paths = sorted(glob.glob(os.path.join(rd, args.glob)))
            if len(paths) < 4:
                print("FAIL %s: only %d dumps under %s"
                      % (rd, len(paths), args.glob))
                return 1
            kvec = wave_vector(dd, args.lx, args.ly, args.lz)
            dumps = common.load_series(paths, args.bin_convert_dir,
                                       args.all_ranks)
            peers[lb] = peer
            rows.append((lb, dd, kvec, dumps))

    out = {"c_s": cs, "beta": beta, "Gamma_1": g1, "Gamma_3m1": g3m1,
           "arms": {}}
    ok = True
    common.report(args, "%-22s %-4s %12s %10s %12s %10s"
                  % ("arm", "dir", "c_phase", "c/c_s-1", "gamma[1/t]",
                     "decay/per"))
    for lb, dd, kvec, dumps in rows:
        kmag = math.sqrt(sum(k * k for k in kvec))
        times, omega, grate, mag = measure(dumps, kvec)
        cph = omega / kmag
        rel = cph / cs - 1.0
        period = 2.0 * math.pi / (cs * kmag)
        per_period = grate * period
        peer = peers.get(lb, "")
        if peer == "none":
            good = True
        elif peer:
            good = (per_period <= args.growth_tol)
        else:
            good = (abs(rel) <= args.tol) and (per_period
                                               <= args.growth_tol)
        ok = ok and good
        common.report(args, "%-22s %-4s %12.8g %10.2e %12.4e %10.2e%s"
                      % (lb, dd, cph, rel, grate, per_period,
                         "" if good else "   <-- FAIL"))
        out["arms"][lb] = {
            "dir": dd, "kmag": kmag, "c_phase": cph, "rel_err": rel,
            "gamma": grate, "decay_per_period": per_period,
            "period": period, "pass": bool(good), "peer": peer,
            "t": [float(t) for t in times],
            "amp": [float(a) for a in mag],
        }

    # the cross-arm comparison: every arm that names a PEER against that peer
    for lb, _, _, _ in rows:
        peer = peers.get(lb, "")
        if not peer or peer == "none" or peer not in out["arms"]:
            continue
        dc = out["arms"][lb]["c_phase"] / out["arms"][peer]["c_phase"] - 1.0
        g0 = out["arms"][peer]["gamma"]
        dg = (out["arms"][lb]["gamma"] - g0) / abs(g0) if g0 != 0.0 else 0.0
        hit = (abs(dc) <= args.match_tol) and (abs(dg) <= args.match_tol)
        ok = ok and hit
        out["arms"][lb]["dc_vs_peer"] = dc
        out["arms"][lb]["dgamma_vs_peer"] = dg
        out["arms"][lb]["match"] = bool(hit)
        common.report(args, "%-22s vs %-18s dc/c=%9.2e dgamma/gamma=%9.2e%s"
                      % (lb, peer, dc, dg, "" if hit else "   <-- FAIL"))

    if args.json:
        import json
        with open(args.json, "w") as fp:
            json.dump(out, fp)
        common.report(args, "wrote %s" % args.json)

    common.verdict(ok, "T10 radwave: %d arm(s), |c/c_s-1| <= %.1e, growth "
                       "<= %.1e per period, peers within %.1e"
                   % (len(rows), args.tol, args.growth_tol, args.match_tol))


if __name__ == "__main__":
    sys.exit(main())
