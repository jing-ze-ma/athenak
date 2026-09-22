#!/usr/bin/env python3
"""Tables for tests_cs_loop_growth: per-region L1(B) vs t, its growth rate, KE and max|v|.

Reads the `*.user.hst` files that CSTestLoopHistory writes (columns are documented in
their own header line) and prints (a) the time series per region, (b) a LINEAR fit
dL1/dt and an EXPONENTIAL fit d(ln L1)/dt over a chosen window, so that a secular
truncation drift (linear, with an h^2 rate) can be told from an instability
(exponential, with a rate that does not converge).
"""
import sys
import numpy as np

COLS = {"t": 0, "dt": 1, "l1in": 2, "nin": 3, "l1sm": 4, "nsm": 5, "l1vx": 6,
        "nvx": 7, "kein": 8, "kesm": 9, "kevx": 10, "me": 11, "mvin": 12,
        "mvsm": 13, "mvvx": 14, "divb": 15, "bmax": 16}


def load(path):
    d = np.loadtxt(path, comments="#")
    out = {k: d[:, i] for k, i in COLS.items()}
    b = out["bmax"]
    for r, n in (("in", "nin"), ("sm", "nsm"), ("vx", "nvx")):
        out["L1" + r] = out["l1" + r] / out[n] / b
    out["KE"] = out["kein"] + out["kesm"] + out["kevx"]
    return out


def fits(t, y, t0, t1):
    m = (t >= t0) & (t <= t1) & (y > 0)
    if m.sum() < 3:
        return float("nan"), float("nan")
    lin = np.polyfit(t[m], y[m], 1)[0]
    exp = np.polyfit(t[m], np.log(y[m]), 1)[0]
    return lin, exp


def series(tag, path, times=(0.0, 0.125, 0.25, 0.5, 1.0, 1.5, 2.0)):
    d = load(path)
    t = d["t"]
    row = [tag]
    for r in ("in", "sm", "vx"):
        vals = [d["L1" + r][np.argmin(abs(t - tt))] for tt in times]
        row.append((r, vals))
    return d, row


def main(args):
    times = (0.0, 0.125, 0.25, 0.5, 1.0, 1.5, 2.0)
    print("t =", "  ".join("%.3f" % x for x in times))
    print()
    for path in args:
        tag = path.split("/")[-1].replace(".user.hst", "")
        d, row = series(tag, path, times)
        t = d["t"]
        print("== %s" % tag)
        for r, vals in row[1:]:
            lin, exp = fits(t, d["L1" + r], 0.25, 2.0)
            print("  L1(B) %-3s " % r + " ".join("%.4e" % v for v in vals) +
                  "   dL1/dt=%.3e  dlnL1/dt=%.3f" % (lin, exp))
        klin, kexp = fits(t, d["KE"], 0.25, 2.0)
        print("  KE        " +
              " ".join("%.4e" % d["KE"][np.argmin(abs(t - tt))] for tt in times) +
              "   dKE/dt=%.3e  dlnKE/dt=%.3f" % (klin, kexp))
        print("  max|v| in/sm/vx at t=2: %.4e %.4e %.4e  (t=0.25: %.4e %.4e %.4e)"
              % (d["mvin"][-1], d["mvsm"][-1], d["mvvx"][-1],
                 d["mvin"][np.argmin(abs(t - 0.25))],
                 d["mvsm"][np.argmin(abs(t - 0.25))],
                 d["mvvx"][np.argmin(abs(t - 0.25))]))
        print("  ME(0)=%.8e  ME(end)=%.8e  dME/ME=%.2e   max|divB|dx/|B| end=%.2e"
              % (d["me"][0], d["me"][-1], d["me"][-1] / d["me"][0] - 1.0,
                 d["divb"][-1]))
        print()


if __name__ == "__main__":
    main(sys.argv[1:])
