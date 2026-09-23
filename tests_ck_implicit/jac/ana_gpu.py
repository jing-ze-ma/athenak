"""GPU cost table of jac/submit_gpu.sh: cpu time used per arm and repeat, x semi, and
the ck_implicit report lines (both ranks, calls after the first 20 cycles): mean / max
passes, non-converged calls, subrst / subfail / escx / aa totals.
usage: python3 ana_gpu.py [dir] [arms]
(default dir /viper/ptmp2/jinma/ckjac_0923/gpu/hyd)"""
import os
import re
import sys

D = sys.argv[1] if len(sys.argv) > 1 else "/viper/ptmp2/jinma/ckjac_0923/gpu/hyd"
ARMS = (sys.argv[2].split(",") if len(sys.argv) > 2
        else ("s", "t4", "sub", "esc", "aa", "aae"))
N0 = 986896 + 20


def cpu(d):
    f = os.path.join(D, d, "run.log")
    if not os.path.exists(f):
        return float("nan")
    for ln in open(f):
        if "cpu time used" in ln:
            return float(ln.split("=")[1])
    return float("nan")


def conv(d):
    ps, nc, tot = [], 0, [0, 0, 0, 0]
    for ln in open(os.path.join(D, d, "run.log")):
        if not ln.startswith("### ck_implicit"):
            continue
        if int(re.search(r"ncycle=(\d+)", ln).group(1)) < N0:
            continue
        ps.append(int(re.search(r"passes=(\d+)", ln).group(1)))
        nc += "NOT-CONVERGED" in ln
        for n, k in enumerate(("subrst", "subfail", "escx", "aa")):
            m = re.search(" " + k + r"=(\d+)", ln)
            if m:
                tot[n] += int(m.group(1))
    return ps, nc, tot


s = [cpu("s_r%d" % r) for r in (1, 2, 3)]
sm = sum(s) / 3
print("%-6s %24s %7s %8s %5s %7s %6s %6s %6s %6s" % (
    "arm", "cpu r1/r2/r3 [s]", "x semi", "passes", "maxp", "nonconv", "subrst", "sfail",
    "escx", "aa"))
for a in ARMS:
    t = [cpu("%s_r%d" % (a, r)) for r in (1, 2, 3)]
    if a == "s":
        print("%-6s %24s %7.2f" % (a, "/".join("%.2f" % x for x in t), sum(t) / 3 / sm))
        continue
    if not os.path.exists(os.path.join(D, "%s_r1" % a, "run.log")):
        continue
    ps, nc, tot = conv("%s_r1" % a)
    print("%-6s %24s %7.2f %8.2f %5d %7d %6d %6d %6d %6d" % (
        a, "/".join("%.2f" % x for x in t), sum(t) / 3 / sm, sum(ps) / max(len(ps), 1),
        max(ps + [0]), nc, *tot))
