#!/usr/bin/env python3
"""Compact G1 table: per (P, tau) the median successive-difference order and e128/e1024
for two run sets side by side.   [EXCL=vetf|edd_] python3 med_table.py RUNS_A RUNS_B [p_R=2]"""
import subprocess
import os
import sys
import statistics

pr = sys.argv[3] if len(sys.argv) > 3 else "2"


def rows(d):
    out = subprocess.run([sys.executable, __file__.replace("med_table", "time_table"), d,
                          "tshesdirk2", pr, os.environ.get("EXCL", "")],
                         capture_output=True, text=True).stdout
    r = {}
    for ln in out.splitlines():
        a, b = ln.split("|")[0], ln.split("|")[1]
        ent = b.split()
        ps = [float(e.split("(")[1].rstrip(")")) for e in ent if "(" in e and "()" not in e]
        es = {int(e.split(":")[0]): e.split(":")[1].split("(")[0].rstrip("h") for e in ent}
        r[tuple(a.split())] = (statistics.median(ps), es.get(128), es.get(1024), min(ps))
    return r


A, B = rows(sys.argv[1]), rows(sys.argv[2])
print("%6s %7s | %-30s | %-30s" % ("P", "tau", "A: med p (min) e128 e1024",
                                    "B: med p (min) e128 e1024"))
for k in sorted(A, key=lambda x: (float(x[0]), float(x[1]))):
    a, b = A[k], B.get(k, (float("nan"),) * 4)
    print("%6s %7s | %.2f (%.2f) %s %s | %.2f (%.2f) %s %s" % (k[0], k[1], a[0], a[3], a[1],
                                                               a[2], b[0], b[3], b[1], b[2]))
