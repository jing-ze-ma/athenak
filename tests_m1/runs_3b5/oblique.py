#!/usr/bin/env python3
"""Milestone 3b phase D, gate G-oblique: how much the OFF-DIAGONAL Eddington terms
matter when the radiation flux is not along a grid axis.

The test problem is the 2-D thick-pulse pgen run at a SMALL cell optical depth, where
the reduced flux of the expanding pulse reaches f ~ 1 and the flux direction n sweeps
through every angle: the Eddington tensor P = D(chi, n) E is then strongly anisotropic
and its off-diagonal components D_12 = (3 chi - 1)/2 n_1 n_2 are largest at 45 degrees
to the grid.  The exact solution of a circular pulse in a uniform medium stays
CIRCULARLY SYMMETRIC, so the deviation from circular symmetry of the computed E is a
direct measure of the error a scheme makes on an oblique flux.

Reported per run, from the LAST dump:
  aniso = rms_r,phi( E - <E>_phi(r) ) / max E     (0 for the exact solution)
  ax45  = max_r |<E>_diag(r) - <E>_axis(r)| / max E, the 4-fold part of it alone
and, with --ref, the max and L1 relative difference of E against a reference run.

Usage: oblique.py --label NAME --dir RUNDIR [--ref REFDIR]
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import common  # noqa: E402


def last_dump(d):
    bd = os.path.join(d, "bin")
    fs = sorted(os.path.join(bd, f) for f in os.listdir(bd))
    return common.load_dump(fs[-1])


def radial(dump, rmax):
    """E(k=0) on the 2-D plane, its radius, and the azimuthally averaged profile."""
    e = dump.var("m1_e")[0]
    x = dump.x1 - 0.5 * (dump.x1[0] + dump.x1[-1])
    y = dump.x2 - 0.5 * (dump.x2[0] + dump.x2[-1])
    xx, yy = np.meshgrid(x, y)
    rr = np.sqrt(xx * xx + yy * yy)
    dr = float(dump.x1[1] - dump.x1[0])
    nb = int(rmax / dr)
    ib = np.minimum((rr / dr).astype(int), nb)
    prof = np.zeros(nb + 1)
    cnt = np.zeros(nb + 1)
    np.add.at(prof, ib, e)
    np.add.at(cnt, ib, 1.0)
    prof /= np.maximum(cnt, 1.0)
    return e, rr, ib, prof, dr, nb


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--label", required=True)
    p.add_argument("--dir", required=True)
    p.add_argument("--ref", default=None)
    p.add_argument("--rmax", type=float, default=0.4)
    a = p.parse_args()

    d = last_dump(a.dir)
    e, rr, ib, prof, dr, nb = radial(d, a.rmax)
    emax = float(e.max())
    m = rr < a.rmax
    aniso = float(np.sqrt(np.mean((e[m] - prof[ib][m]) ** 2)) / emax)

    # the 4-fold part: a wedge around the axes against a wedge around the diagonals
    x = d.x1 - 0.5 * (d.x1[0] + d.x1[-1])
    y = d.x2 - 0.5 * (d.x2[0] + d.x2[-1])
    xx, yy = np.meshgrid(x, y)
    ph = np.arctan2(yy, xx)
    c4 = np.cos(4.0 * ph)
    ax45 = 0.0
    for ibin in range(2, nb):
        sel = (ib == ibin)
        if sel.sum() < 16:
            continue
        de = e[sel] - prof[ibin]
        # amplitude of the cos(4 phi) component, normalised by the peak
        amp = 2.0 * np.mean(de * c4[sel])
        ax45 = max(ax45, abs(amp) / emax)

    txt = "%-26s t=%.4g aniso=%.4e cos4phi=%.4e" % (a.label, d.time, aniso, ax45)
    if a.ref:
        r = last_dump(a.ref)
        er = r.var("m1_e")[0]
        sc = float(np.abs(er).max())
        txt += " vs ref: max=%.4e L1=%.4e" % (float(np.abs(e - er).max() / sc),
                                              float(np.abs(e - er).mean() / sc))
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
