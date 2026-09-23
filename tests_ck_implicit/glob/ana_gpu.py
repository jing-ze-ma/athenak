"""GPU cost table of submit_gpu.sh: cpu time used per arm and repeat, x semi, and the
ck_implicit report lines (both ranks, calls after the first 20 cycles): mean / max
passes, non-converged calls, lsrej / subrst / subfail totals.
usage: python3 ana_gpu.py [dir]   (default /viper/ptmp2/jinma/ckglob_0924/gpu/hyd)"""
import os
import re
import sys

D = sys.argv[1] if len(sys.argv) > 1 else "/viper/ptmp2/jinma/ckglob_0924/gpu/hyd"


def cpu(d):
    for ln in open(os.path.join(D, d, "run.log")):
        if "cpu time used" in ln:
            return float(ln.split("=")[1])
    return float("nan")


def conv(d, n0=873600):
    ps, nc, tot = [], 0, [0, 0, 0]
    for ln in open(os.path.join(D, d, "run.log")):
        if not ln.startswith("### ck_implicit"):
            continue
        if int(re.search(r"ncycle=(\d+)", ln).group(1)) < n0:
            continue
        ps.append(int(re.search(r"passes=(\d+)", ln).group(1)))
        nc += "NOT-CONVERGED" in ln
        for n, k in enumerate(("lsrej", "subrst", "subfail")):
            m = re.search(k + r"=(\d+)", ln)
            if m:
                tot[n] += int(m.group(1))
    return ps, nc, tot


s = [cpu("s_r%d" % r) for r in (1, 2, 3)]
sm = sum(s) / 3
print("%-6s %24s %7s %8s %5s %7s %6s %6s %6s" % ("arm", "cpu r1/r2/r3 [s]", "x semi",
                                                 "passes", "maxp", "nonconv", "lsrej",
                                                 "subrst", "sfail"))
for a in ("s", "t4", "ls", "sub"):
    t = [cpu("%s_r%d" % (a, r)) for r in (1, 2, 3)]
    if a == "s":
        print("%-6s %24s %7.2f" % (a, "/".join("%.2f" % x for x in t), sum(t) / 3 / sm))
        continue
    ps, nc, tot = conv("%s_r1" % a)
    print("%-6s %24s %7.2f %8.2f %5d %7d %6d %6d %6d" % (
        a, "/".join("%.2f" % x for x in t), sum(t) / 3 / sm, sum(ps) / max(len(ps), 1),
        max(ps + [0]), nc, *tot))
