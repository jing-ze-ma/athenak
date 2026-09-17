#!/usr/bin/env python3
"""Cycle-0 read-back check: the run's own shell means at t = 0 against the ic file.

Usage: icchk.py <rt_profile.bin> <ic_file>

rt_profile.bin slots (red_giant.cpp:494): 0 rho 1 v1 2 rho*v1 3 v1^2 4 v2^2+v3^2
5 wtemp 6 eint 7 v1(eint+p).  This does NOT re-derive the EOS; it asks whether the
profile the code built at cycle 0 is the profile the ic file specifies, which is the
half of the ic check that the offline solve in make_ic_sph.py cannot do for itself.
"""
import sys
import numpy as np

RSTAR = 2.3717e11


def read_prof(fn):
    ts, qs, x1v = [], [], None
    with open(fn, "rb") as f:
        while True:
            b = f.read(8)
            if len(b) < 8:
                break
            t = np.frombuffer(b, "<f8")[0]
            nx1, nvar = np.frombuffer(f.read(8), "<i4")
            x = np.frombuffer(f.read(8 * nx1), "<f8")
            q = np.frombuffer(f.read(8 * nvar * nx1), "<f8").reshape(nvar, nx1)
            if x1v is None:
                x1v = x.copy()
            ts.append(t)
            qs.append(q.copy())
    return np.array(ts), x1v, np.array(qs)


def main():
    t, x1v, q = read_prof(sys.argv[1])
    ic = np.loadtxt(sys.argv[2], comments="#")
    print("%d records, %d radial nodes, t[0] = %.4g s" % (len(t), len(x1v), t[0]))
    rho = np.interp(x1v, ic[:, 0], ic[:, 1])
    eint = np.interp(x1v, ic[:, 0], ic[:, 2])
    dr = np.abs(q[0, 0] / rho - 1.0)
    de = np.abs(q[0, 6] / eint - 1.0)
    print("  |d rho/rho| : max %.3e  rms %.3e" % (dr.max(), np.sqrt((dr ** 2).mean())))
    print("  |d eint/eint|: max %.3e  rms %.3e" % (de.max(), np.sqrt((de ** 2).mean())))
    i = de.argmax()
    print("  worst eint node: r/R %.4f  code %.6e  ic %.6e"
          % (x1v[i] / RSTAR, q[0, 6][i], eint[i]))
    # the taper window, which is the only place the gate changes the ic at all
    w = (x1v > 2.3386e11) & (x1v < 2.3826e11)
    if w.any():
        print("  inside the taper window (%d nodes): |d eint/eint| max %.3e"
              % (w.sum(), de[w].max()))


if __name__ == "__main__":
    main()
