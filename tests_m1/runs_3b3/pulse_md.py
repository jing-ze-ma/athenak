#!/usr/bin/env python3
"""Milestone 3b phase B, gate G1: the ISOTROPY and the RATE of the thick pulse in
2-D / 3-D with <rad_m1>/transport = implicit.

The pulse is a Gaussian in a static pure-scattering medium, so the variance of the
marginal distribution of E must grow at 2 D in EVERY direction, with
D = c/(3 rho kappa).  The script reports, per axis, the fitted d(sigma^2)/dt divided by
2 D and the spread between the axes.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import common  # noqa: E402


def marginal_variance(dump, bg):
    e = dump.var("m1_e") - bg
    e = np.maximum(e, 0.0)
    out = []
    axes = [(dump.x1, (0, 1)), (dump.x2, (0, 2)), (dump.x3, (1, 2))]
    for x, red in axes:
        if len(x) < 2:
            out.append(None)
            continue
        w = e.sum(axis=red)
        s = w.sum()
        mu = (w * x).sum() / s
        out.append(((w * (x - mu) ** 2).sum()) / s)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("files", nargs="+")
    p.add_argument("--c", type=float, default=1.0)
    p.add_argument("--rho", type=float, default=1.0)
    p.add_argument("--kappa", type=float, required=True)
    p.add_argument("--tol", type=float, default=0.02)
    p.add_argument("--iso-tol", type=float, default=1.0e-3)
    p.add_argument("--label", default="")
    a = p.parse_args()

    dumps = common.load_series(sorted(a.files))
    bg = min(d.var("m1_e").min() for d in dumps)
    t = np.array([d.time for d in dumps])
    var = [marginal_variance(d, bg) for d in dumps]
    diff = a.c / (3.0 * a.rho * a.kappa)
    ratios = []
    for ax in range(3):
        if var[0][ax] is None:
            ratios.append(None)
            continue
        s2 = np.array([v[ax] for v in var])
        slope, _ = common.linfit(t, s2)
        ratios.append(slope / (2.0 * diff))
    got = [r for r in ratios if r is not None]
    iso = (max(got) - min(got)) / max(abs(max(got)), 1e-300) if len(got) > 1 else 0.0
    ok = all(abs(r - 1.0) <= a.tol for r in got) and iso <= a.iso_tol
    txt = " ".join("%s=%.6f" % (n, r) for n, r in
                   zip("xyz", ratios) if r is not None)
    print("%-28s %s iso=%.3e %s" % (a.label, txt, iso, "PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
