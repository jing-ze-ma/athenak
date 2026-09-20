#!/usr/bin/env python3
"""tests_r22: build a problem/ic_profile file from a (finished or running) 3-D wedge.

The low-resolution relaxation arm is stopped once its MEAN state is settled; this
script turns that mean state into a new 1-D initial stratification for the
high-resolution run:

    shell-mean rho and eint  (rt_profile.bin slots 0 and 6, already solid-angle
    weighted by the dump itself -- see red_giant.cpp:505-550)

on the run's own radial grid (the stretched cell centres stored in every record),
optionally averaged over a window of records and optionally boxcar-smoothed in
radius, written back onto the ORIGINAL IC file's rows: rows outside the run's
radial range keep their original values, and the two joins are blended over a few
cells so the file stays smooth.  Format and row density are the original's.

Units: the dump's slot 6 is the internal energy DENSITY in erg/cm^3 (the
EintFromCons value, rho*Phi already removed), i.e. exactly the third column of the
IC file; slot 0 is g/cm^3, the first-column companion.  --check verifies this by
round-tripping the t = 0 record against the IC file (see README.md).

Examples
    python3 mk_ic_from_run.py --check tests_r11/wedge10
    python3 mk_ic_from_run.py tests_r22/lr1 --tmin 5.0 --tmax 6.0 --smooth 3 \\
        --ic tests_r14/ic_he4_tall_neww.txt --out tests_r22/ic_he4_relaxed.txt
"""
import argparse
import os
import struct
import sys

import numpy as np

TURN = 4705.0            # one convective turnover [s]
RSTAR = 2.3717e11        # tau = 2/3 photosphere [cm]


def rd(fn):
    """tests_r14/dtdiag.py rd(): [(t, x1v[nx1], q[nvar][nx1]), ...]."""
    out = []
    with open(fn, 'rb') as f:
        while True:
            h = f.read(16)
            if len(h) < 16:
                break
            t, n1, nv = struct.unpack('<dii', h)
            b = f.read(8*(n1 + nv*n1))
            if len(b) < 8*(n1 + nv*n1):
                break
            a = np.frombuffer(b, '<f8')
            out.append((t, a[:n1].copy(), a[n1:].reshape(nv, n1).copy()))
    return out


def loglin(r, rf, yf):
    """The code's own IC interpolant: linear in r, logarithmic in the value
    (red_giant.cpp:2325-2331)."""
    ly = np.interp(r, rf, np.log(np.maximum(yf, 1e-300)))
    return np.exp(ly)


def boxcar(y, n):
    """Boxcar of n cells in LOG space, edge-clamped.  n <= 1 is a no-op."""
    if n is None or n <= 1:
        return y
    k = int(n) | 1                      # force odd
    p = k//2
    ly = np.log(np.maximum(y, 1e-300))
    lp = np.concatenate([np.full(p, ly[0]), ly, np.full(p, ly[-1])])
    ker = np.ones(k)/k
    return np.exp(np.convolve(lp, ker, mode='valid'))


def mean_profile(recs, tmin, tmax):
    """Record-average of slots 0 (rho) and 6 (eint) over t in [tmin, tmax] turnovers.

    The average is taken in LOG space, which is what the stratification is smooth in
    and what the file's interpolant assumes.  Returns (r, rho, eint, tlist)."""
    sel = [x for x in recs
           if (tmin is None or x[0] >= tmin*TURN) and
              (tmax is None or x[0] <= tmax*TURN)]
    if not sel:
        sys.exit("no rt_profile records in the requested window")
    r = sel[0][1]
    for t, rr, _ in sel:
        if rr.shape != r.shape or not np.allclose(rr, r):
            sys.exit("the radial grid changes between records: refusing to average")
    ld = np.mean([np.log(np.maximum(q[0], 1e-300)) for _, _, q in sel], axis=0)
    le = np.mean([np.log(np.maximum(q[6], 1e-300)) for _, _, q in sel], axis=0)
    return r, np.exp(ld), np.exp(le), [x[0] for x in sel]


def ramp(x):
    """Smoothstep on [0, 1]."""
    s = np.clip(x, 0.0, 1.0)
    return s*s*(3.0 - 2.0*s)


def build(rf_ic, y_ic, r, y, nblend):
    """Original rows outside [r[0], r[-1]] keep y_ic; inside they take y, with a
    smoothstep blend over nblend profile cells at each join.  Log space throughout."""
    y_new = np.array(y_ic, dtype=float)
    inside = (rf_ic >= r[0]) & (rf_ic <= r[-1])
    y_new[inside] = loglin(rf_ic[inside], r, y)
    if nblend > 0 and inside.any():
        dlo = np.mean(np.diff(r[:max(2, nblend+1)]))
        dhi = np.mean(np.diff(r[-max(2, nblend+1):]))
        wlo = ramp((rf_ic - r[0])/(nblend*dlo))
        whi = ramp((r[-1] - rf_ic)/(nblend*dhi))
        w = np.minimum(wlo, whi)
        w[~inside] = 0.0
        lo = np.log(np.maximum(y_ic, 1e-300))
        ln = np.log(np.maximum(y_new, 1e-300))
        y_new = np.exp(w*ln + (1.0 - w)*lo)
    return y_new


def check(run, icfile):
    """Round-trip: the t = 0 record of a run against the IC file it started from."""
    recs = rd(os.path.join(run, 'rt_profile.bin'))
    t, r, q = recs[0]
    d = np.loadtxt(icfile)
    d0 = loglin(r, d[:, 0], d[:, 1])
    e0 = loglin(r, d[:, 0], d[:, 2])
    print("round-trip %s  t = %.3f s  nx1 = %d  r = %.4e .. %.4e cm"
          % (run, t, len(r), r[0], r[-1]))
    for nm, a, b in (("rho ", q[0], d0), ("eint", q[6], e0)):
        rel = np.abs(a/b - 1.0)
        j = int(np.argmax(rel))
        print("  %s  median %.3e   max %.3e at i = %d, r/R = %.4f"
              % (nm, np.median(rel), rel[j], j, r[j]/RSTAR))
        for lo, hi in ((0.50, 0.90), (0.90, 0.97), (0.97, 1.01), (1.01, 1.30)):
            m = (r/RSTAR >= lo) & (r/RSTAR < hi)
            if m.any():
                print("      r/R %.2f-%.2f : median %.2e  max %.2e"
                      % (lo, hi, np.median(rel[m]), rel[m].max()))


def bands(r, rho, eint, d_ic, e_ic):
    print("  band        n   <rho/rho_IC - 1>   max|.|      <eint/eint_IC - 1>  max|.|")
    for lo, hi in ((0.50, 0.70), (0.70, 0.90), (0.90, 0.97),
                   (0.97, 1.01), (1.01, 1.10), (1.10, 1.30)):
        m = (r/RSTAR >= lo) & (r/RSTAR < hi)
        if not m.any():
            continue
        a = rho[m]/d_ic[m] - 1.0
        b = eint[m]/e_ic[m] - 1.0
        print("  %.2f-%.2f %4d   %+10.3e      %9.3e   %+10.3e      %9.3e"
              % (lo, hi, m.sum(), np.mean(a), np.abs(a).max(),
                 np.mean(b), np.abs(b).max()))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('run', help='run directory holding rt_profile.bin')
    p.add_argument('--ic', default=os.path.join(os.path.dirname(__file__), '..',
                                                'tests_r14', 'ic_he4_tall_neww.txt'),
                   help='the ORIGINAL ic_profile file (row grid + the rows outside '
                        'the mesh are taken from it)')
    p.add_argument('--out', default=None, help='output ic_profile file')
    p.add_argument('--tmin', type=float, default=None, help='window start [turnovers]')
    p.add_argument('--tmax', type=float, default=None, help='window end [turnovers]')
    p.add_argument('--smooth', type=int, default=0,
                   help='radial boxcar in CELLS of the run grid (default 0 = none)')
    p.add_argument('--blend', type=int, default=4,
                   help='blend width in run cells at each join (default 4)')
    p.add_argument('--check', action='store_true',
                   help='only round-trip the t = 0 record against --ic and exit')
    a = p.parse_args()

    icfile = os.path.abspath(a.ic)
    if a.check:
        check(a.run, icfile)
        return

    recs = rd(os.path.join(a.run, 'rt_profile.bin'))
    if not recs:
        sys.exit("no records in %s/rt_profile.bin" % a.run)
    tmin, tmax = a.tmin, a.tmax
    if tmin is None and tmax is None:            # default: the last record alone
        tmin = tmax = recs[-1][0]/TURN
    r, rho, eint, tl = mean_profile(recs, tmin, tmax)
    rho = boxcar(rho, a.smooth)
    eint = boxcar(eint, a.smooth)

    d = np.loadtxt(icfile)
    rf_ic = d[:, 0]
    d_new = build(rf_ic, d[:, 1], r, rho, a.blend)
    e_new = build(rf_ic, d[:, 2], r, eint, a.blend)

    print("mk_ic_from_run: %s, %d record(s), t = %.1f .. %.1f s = %.3f .. %.3f turnovers"
          % (a.run, len(tl), tl[0], tl[-1], tl[0]/TURN, tl[-1]/TURN))
    print("  run grid %d cells, %.4e .. %.4e cm (r/R %.4f .. %.4f); smooth %d, blend %d"
          % (len(r), r[0], r[-1], r[0]/RSTAR, r[-1]/RSTAR, a.smooth, a.blend))
    bands(r, rho, eint, loglin(r, rf_ic, d[:, 1]), loglin(r, rf_ic, d[:, 2]))

    if a.out is None:
        print("  (no --out given: nothing written)")
        return
    hdr = ("# tests_r22/mk_ic_from_run.py: shell-mean (rho, eint) of %s\n"
           "# records t = %.1f .. %.1f s (%.3f .. %.3f turnovers), %d record(s),\n"
           "# radial boxcar %d cells, join blend %d cells, rows outside\n"
           "# [%.6e, %.6e] cm taken unchanged from %s\n"
           "# r[cm]  rho[g/cm^3]  eint[erg/cm^3]"
           % (os.path.abspath(a.run), tl[0], tl[-1], tl[0]/TURN, tl[-1]/TURN, len(tl),
              a.smooth, a.blend, r[0], r[-1], icfile))
    np.savetxt(a.out, np.column_stack([rf_ic, d_new, e_new]),
               fmt='%.10e', header=hdr, comments='')
    print("  wrote %s (%d rows)" % (a.out, len(rf_ic)))


if __name__ == '__main__':
    main()
