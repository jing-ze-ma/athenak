#!/usr/bin/env python3
"""runs_3s_space2: the SPACE table of tables.py (Richardson in time, nt/2nt) for every
implicit_enthalpy arm, a headline (orders 32->64 and 64->128, errors at N = 64) and the
Picard statistics of every run (mean/max passes, NON-CONVERGED) from its log.

  python3 space2_table.py [--runs DIR] [--txt OUT]
"""
import argparse
import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tables as tb  # noqa: E402


def picard(path):
    try:
        s = open(os.path.join(path, "log.txt")).read()
    except OSError:
        return None
    m = re.search(r"solves=(\S+) Picard iterations mean=(\S+) max=(\S+) "
                  r"NON-CONVERGED=(\S+)", s)
    if not m:
        return None
    q = re.search(r"inner iterations mean=\S+ max=\S+ total=(\S+)", s)
    inner = float(q.group(1)) / float(m.group(1)) if q else float("nan")
    return float(m.group(2)), float(m.group(3)), float(m.group(4)), inner


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--runs", default="/viper/ptmp2/jinma/space2_3s/runs")
    p.add_argument("--txt", default=os.path.join(HERE, "RESULTS_space.txt"))
    a = p.parse_args()
    rs = tb.load(a.runs)
    with open(a.txt, "w") as fp:
        sp = tb.space_table(rs, fp)
        fp.write("\n== HEADLINE  p = log2 ratio of the Richardson errors; e_c64, e_g64, "
                 "L1_64 at N = 64\n")
        fp.write("%-28s %6s %7s %6s %6s %6s %6s %10s %10s %10s\n"
                 % ("arm", "Pr/Pg", "tau", "pw1", "pw2", "pL1", "pL2", "e_c64",
                    "e_g64", "L1_64"))
        for arm in sorted(sp):
            rows = sp[arm]
            for pr, tau in tb.PT:
                g = {n: rows.get((pr, tau, n)) for n in (32, 64, 128)}
                if g[64] is None:
                    continue

                def od(k, n1, n2):
                    if g[n1] is None or g[n2] is None:
                        return "-"
                    return "%.2f" % math.log2(g[n1][k] / g[n2][k])
                fp.write("%-28s %6g %7g %6s %6s %6s %6s %10.2e %10.2e %10.2e\n"
                         % (arm, pr, tau, od("err_w", 32, 64), od("err_w", 64, 128),
                            od("l1", 32, 64), od("l1", 64, 128), g[64]["err_c"],
                            g[64]["err_g"], g[64]["l1"]))
            fp.write("\n")
        fp.write("\n== PICARD  per arm: passes mean over runs, max, NON-CONVERGED sum, BiCGStab inner iterations per step (multi-D only)\n")
        agg = {}
        for r in rs:
            q = picard(r["path"])
            if q is None:
                continue
            ar = tb.arm_of(r)
            x = agg.setdefault(ar, [0.0, 0, 0.0, 0.0, 0.0])
            x[0] += q[0]
            x[1] += 1
            x[2] = max(x[2], q[1])
            x[3] += q[2]
            x[4] += q[3]
        for ar in sorted(agg):
            x = agg[ar]
            fp.write("%-40s runs=%4d mean=%.4f max=%g NONCONV=%g inner/step=%.2f\n"
                     % (ar, x[1], x[0] / x[1], x[2], x[3], x[4] / x[1]))
    print(open(a.txt).read())


if __name__ == "__main__":
    main()
