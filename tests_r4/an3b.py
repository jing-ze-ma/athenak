#!/usr/bin/env python3
"""G3: F_2s/F_req of the optically thick spherical two-stream test.

Restricted to the interior faces with dtau_cell > 0.3 -- below that the diffusion
limit itself stops holding -- with the wall face and the top ghost mirror dropped,
exactly as tests_r2/thick/README.md section 3 does.
"""
import os
import sys
import numpy as np

KB = 1.380649e-16
MH = 1.6726219e-24
HERE = os.path.dirname(os.path.abspath(__file__))
CASES = os.path.join(HERE, "..", "tests_r2", "thick", "cases.txt")


def load_cases():
    out = {}
    for row in open(CASES):
        if row.startswith("#") or not row.strip():
            continue
        f = row.split()
        out[f[0]] = (float(f[4]), float(f[5]))       # x1max, kappa_const
    return out


def main(base):
    print("%-16s %5s %9s %9s %9s %9s"
          % ("case", "N", "min", "max", "median", "max|1-x|"))
    for tag, (x1max, kap) in load_cases().items():
        d = np.loadtxt(os.path.join(base, tag, "mltfaces.txt"))
        t, p = d[:, 2], d[:, 3]
        freq, f2s = d[:, 10], d[:, 17]
        dr = (x1max - 1.0e12)/128.0
        dtau = kap*(p*MH/(KB*t))*dr                  # mu = 1, ideal gas
        ok = (dtau > 0.3) & (freq != 0)
        ok[:6] = False
        ok[-6:] = False
        r = (f2s/freq)[ok]
        print("%-16s %5d %9.5f %9.5f %9.5f %9.5f"
              % (tag, ok.sum(), r.min(), r.max(), np.median(r),
                 np.abs(r - 1).max()))


if __name__ == "__main__":
    for a in sys.argv[1:]:
        main(a)
