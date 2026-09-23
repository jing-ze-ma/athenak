"""Tables for the ck_impl_esc / ck_impl_aa gates.  usage: python3 ana_jac.py A1|A2|D|E
[root] [arms].  root defaults to /viper/ptmp2/jinma/ckjac_0923; reads rec.txt
(wellposed/rec.py) and the ck_implicit report lines of run.log (subrst / subfail / escx /
aa totals, and the per-pass residual history of ck_impl_debug = -2: "rate" is the
median over the calls with >= 6 passes and no sub-step of (r_last/r_3)^(1/(n-4)), the
mean contraction of the max residual per pass after pass 3).
python3 ana_jac.py H <rundir> <call> [<call> ..]: the per-pass max residual of those
calls (0-based report-line index)."""
import os
import re
import sys

import numpy as np

sys.path.insert(0, "/viper/u2/jinma/ATHENAK/athenak/tests_ck_implicit/wellposed")
import rec  # noqa: E402

ROOT = (sys.argv[2] if len(sys.argv) > 2 and sys.argv[1] != "H"
        else "/viper/ptmp2/jinma/ckjac_0923")
WP = "/viper/ptmp2/jinma/wellposed_0923"
ARMS = (sys.argv[3].split(",") if len(sys.argv) > 3
        else ("t4", "sub", "esc", "aa", "aas", "aae"))


def rd(d, root=ROOT):
    f = os.path.join(root, d, "rec.txt")
    return rec.read(f) if os.path.exists(f) else None


def glob_counts(d):
    """sum of subrst, subfail, escx, aa over the report lines, the max pass count and
    the median per-pass contraction of the max residual from pass 3 on, over the calls
    that took >= 6 passes without a sub-step"""
    f = os.path.join(ROOT, d, "run.log")
    tot = [0, 0, 0, 0]
    mp = 0
    rates = []
    if not os.path.exists(f):
        return tot, mp
    for ln in open(f):
        if not ln.startswith("### ck_implicit"):
            continue
        m = re.search(r"passes=(\d+)", ln)
        mp = max(mp, int(m.group(1)))
        for n, k in enumerate(("subrst", "subfail", "escx", "aa")):
            m = re.search(" " + k + r"=(\d+)", ln)
            if m:
                tot[n] += int(m.group(1))
        m = re.search(r"hist=(\S+)", ln)
        sr = re.search(r"subrst=(\d+)", ln)
        if m and (sr is None or int(sr.group(1)) == 0):
            h = [float(x.split("/")[0]) for x in m.group(1).split(",")]
            h = [x for x in h if x > 0.0]
            if len(h) >= 6:
                rates.append((h[-1] / h[3]) ** (1.0 / (len(h) - 4)))
    glob_counts.rate = float(np.median(rates)) if rates else float("nan")
    return tot, mp


def rel(a, b, ic):
    e = np.abs(a[ic:] - b[ic:]) / b[ic:]
    return e.max(), e.mean()


HDR = "%-4s %5s %10s %10s %8s %5s %7s %9s %9s %6s %6s %6s %6s %5s"


def A1():
    t2 = rd("A/a1_t4x_20", WP)
    ic = t2["icut"]
    ts = rec.last(t2)
    print("A1: from T** (wellposed_0923 a1_t4x_20 rst 3), 200 calls; errors vs T**")
    print(HDR % ("arm", "dt", "max@200", "mean@200", "passes", "maxp", "nonconv",
                 "resfinal", "|bud0|", "subrst", "sfail", "escx", "aa", "rate"))
    for dt in (20, 200, 2000):
        for a in ARMS:
            d = "A/a1_%s_%d" % (a, dt)
            r = rd(d)
            if r is None:
                continue
            q = r["s"]
            e = rel(rec.last(r), ts, ic)
            c, mp = glob_counts(d)
            print(("%-4s %5d %10.3e %10.3e %8.2f %5d %7d %9.2e %9.2e %6d %6d %6d %6d"
                   " %5.2f") % (
                a, dt, e[0], e[1], q["passes"].mean(), mp, int(q["nonconv"].sum()),
                q["resfinal"].max(), np.abs(q["bud0"]).max(), *c,
                glob_counts.rate))


def A2():
    ref = rd("A/a2_ref", WP)
    tr = rec.last(ref)
    ic = ref["icut"]
    print("A2: vs t4x dt=4 (wellposed_0923 A/a2_ref), t = 2000 s")
    print("%-4s %5s %10s %10s %8s %5s %7s %9s %6s %6s %6s %6s %5s" % (
        "arm", "dt", "max", "mean", "passes", "maxp", "nonconv", "resfinal", "subrst",
        "sfail", "escx", "aa", "rate"))
    for dt in (20, 200, 2000):
        for a in ARMS:
            d = "A/a2_%s_%d" % (a, dt)
            r = rd(d)
            if r is None:
                continue
            q = r["s"]
            e = rel(rec.last(r), tr, ic)
            c, mp = glob_counts(d)
            print("%-4s %5d %10.3e %10.3e %8.2f %5d %7d %9.2e %6d %6d %6d %6d %5.2f" % (
                a, dt, e[0], e[1], q["passes"].mean(), mp, int(q["nonconv"].sum()),
                q["resfinal"].max(), *c, glob_counts.rate))


def D():
    ref = rd("D/d_true_ref20")
    print("D (eos_h2 true, mu0 1): 100 calls; err vs t4 at dt 20 at the same t "
          "(ref calls 1000 / 10000)" if ref else "D: no dt=20 reference yet")
    print("%-4s %5s %8s %5s %7s %9s %10s %10s %11s %6s %6s %6s %6s %5s" % (
        "arm", "dt", "passes", "maxp", "nonconv", "resfinal", "err@t1", "err@t2",
        "Trange@100", "subrst", "sfail", "escx", "aa", "rate"))
    for dt in (20, 200, 2000):
        for a in ARMS:
            d = "D/d_true_%d_%s" % (dt, a)
            r = rd(d)
            if r is None:
                continue
            ic = r["icut"]
            q = r["s"]
            e1 = e2 = float("nan")
            if ref is not None:
                # t1 = 10 calls, t2 = 100 calls of this arm
                for n, (c, e) in enumerate(((10, None), (100, None))):
                    rc = c * dt // 20
                    if rc in ref["T"] and c in r["T"]:
                        v = rel(r["T"][c], ref["T"][rc], ic)[0]
                        if n == 0:
                            e1 = v
                        else:
                            e2 = v
            tl = r["T"][max(r["T"])][ic:]
            c, mp = glob_counts(d)
            print(("%-4s %5d %8.2f %5d %7d %9.2e %10.3e %10.3e %5.0f-%5.0f %6d %6d %6d"
                   " %6d %5.2f") % (
                a, dt, q["passes"].mean(), mp, int(q["nonconv"].sum()),
                q["resfinal"].max(), e1, e2, tl.min(), tl.max(), *c,
                glob_counts.rate))


def E():
    print("%-16s %6s %8s %11s %8s %7s %6s %6s %6s %6s" % (
        "run", "calls", "ndiffT", "maxrelT", "passes", "nonconv", "subrst", "sfail",
        "escx", "aa"))
    for g in ("tm", "tmb"):
        for dt in (20, 2000):
            for a in ARMS:
                d = "E/e_%s_%d_%s" % (g, dt, a)
                r = rd(d)
                if r is None:
                    continue
                s = r["s"]
                c, mp = glob_counts(d)
                print("%-16s %6d %8d %11.3e %8.2f %7d %6d %6d %6d %6d" % (
                    d[2:], len(s["ncall"]), s["ndiffT"].max(), s["maxrelT"].max(),
                    s["passes"].mean(), int(s["nonconv"].sum()), *c))


def H():
    lines = [ln for ln in open(os.path.join(ROOT, sys.argv[2], "run.log"))
             if ln.startswith("### ck_implicit")]
    for c in sys.argv[3:]:
        ln = lines[int(c)]
        h = re.search(r"hist=(\S+)", ln).group(1).split(",")
        print("%s call %s: %s" % (sys.argv[2], c, " ".join(x.split("/")[0] for x in h)))


if __name__ == "__main__":
    globals()[sys.argv[1]]()
