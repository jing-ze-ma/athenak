"""Tables for the ck_impl_glob gates.  usage: python3 ana_glob.py A1|A2|D|E [root]
root defaults to /viper/ptmp2/jinma/ckglob_0924; reads rec.txt (wellposed/rec.py) and
the ck_implicit report lines of run.log (lsrej / subrst / subfail totals)."""
import os
import re
import sys

import numpy as np

sys.path.insert(0, "/viper/u2/jinma/ATHENAK/athenak/tests_ck_implicit/wellposed")
import rec  # noqa: E402

ROOT = sys.argv[2] if len(sys.argv) > 2 else "/viper/ptmp2/jinma/ckglob_0924"
WP = "/viper/ptmp2/jinma/wellposed_0923"
ARMS = ("t4", "ls", "sub", "s32")


def rd(d, root=ROOT):
    f = os.path.join(root, d, "rec.txt")
    return rec.read(f) if os.path.exists(f) else None


def glob_counts(d):
    """sum of lsrej, subrst, subfail over the report lines, and the max pass count"""
    f = os.path.join(ROOT, d, "run.log")
    tot = [0, 0, 0]
    mp = 0
    if not os.path.exists(f):
        return tot, mp
    for ln in open(f):
        if not ln.startswith("### ck_implicit"):
            continue
        m = re.search(r"passes=(\d+)", ln)
        mp = max(mp, int(m.group(1)))
        for n, k in enumerate(("lsrej", "subrst", "subfail")):
            m = re.search(k + r"=(\d+)", ln)
            if m:
                tot[n] += int(m.group(1))
    return tot, mp


def rel(a, b, ic):
    e = np.abs(a[ic:] - b[ic:]) / b[ic:]
    return e.max(), e.mean()


HDR = "%-4s %5s %10s %10s %8s %5s %7s %9s %9s %6s %6s %6s"


def A1():
    t2 = rd("A/a1_t4x_20", WP)
    ic = t2["icut"]
    ts = rec.last(t2)
    print("A1: from T** (wellposed_0923 a1_t4x_20 rst 3), 200 calls; errors vs T**")
    print(HDR % ("arm", "dt", "max@200", "mean@200", "passes", "maxp", "nonconv",
                 "resfinal", "|bud0|", "lsrej", "subrst", "sfail"))
    for dt in (20, 200, 2000):
        for a in ARMS:
            d = "A/a1_%s_%d" % (a, dt)
            r = rd(d)
            if r is None:
                continue
            q = r["s"]
            e = rel(rec.last(r), ts, ic)
            c, mp = glob_counts(d)
            print("%-4s %5d %10.3e %10.3e %8.2f %5d %7d %9.2e %9.2e %6d %6d %6d" % (
                a, dt, e[0], e[1], q["passes"].mean(), mp, int(q["nonconv"].sum()),
                q["resfinal"].max(), np.abs(q["bud0"]).max(), *c))


def A2():
    ref = rd("A/a2_ref", WP)
    tr = rec.last(ref)
    ic = ref["icut"]
    print("A2: vs t4x dt=4 (wellposed_0923 A/a2_ref), t = 2000 s")
    print("%-4s %5s %10s %10s %8s %5s %7s %9s %6s %6s %6s" % (
        "arm", "dt", "max", "mean", "passes", "maxp", "nonconv", "resfinal", "lsrej",
        "subrst", "sfail"))
    for dt in (20, 200, 2000):
        for a in ARMS:
            d = "A/a2_%s_%d" % (a, dt)
            r = rd(d)
            if r is None:
                continue
            q = r["s"]
            e = rel(rec.last(r), tr, ic)
            c, mp = glob_counts(d)
            print("%-4s %5d %10.3e %10.3e %8.2f %5d %7d %9.2e %6d %6d %6d" % (
                a, dt, e[0], e[1], q["passes"].mean(), mp, int(q["nonconv"].sum()),
                q["resfinal"].max(), *c))


def D():
    ref = rd("D/d_true_ref20")
    print("D (eos_h2 true, mu0 1): 100 calls; err vs t4 at dt 20 at the same t "
          "(ref calls 1000 / 10000)" if ref else "D: no dt=20 reference yet")
    print("%-4s %5s %8s %5s %7s %9s %10s %10s %11s %6s %6s %6s" % (
        "arm", "dt", "passes", "maxp", "nonconv", "resfinal", "err@t1", "err@t2",
        "Trange@100", "lsrej", "subrst", "sfail"))
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
            print("%-4s %5d %8.2f %5d %7d %9.2e %10.3e %10.3e %5.0f-%5.0f %6d %6d %6d" % (
                a, dt, q["passes"].mean(), mp, int(q["nonconv"].sum()),
                q["resfinal"].max(), e1, e2, tl.min(), tl.max(), *c))


def E():
    print("%-16s %6s %8s %11s %8s %7s %6s %6s %6s" % (
        "run", "calls", "ndiffT", "maxrelT", "passes", "nonconv", "lsrej", "subrst",
        "sfail"))
    for g in ("tm", "tmb"):
        for dt in (20, 2000):
            for a in ARMS:
                d = "E/e_%s_%d_%s" % (g, dt, a)
                r = rd(d)
                if r is None:
                    continue
                s = r["s"]
                c, mp = glob_counts(d)
                print("%-16s %6d %8d %11.3e %8.2f %7d %6d %6d %6d" % (
                    d[2:], len(s["ncall"]), s["ndiffT"].max(), s["maxrelT"].max(),
                    s["passes"].mean(), int(s["nonconv"].sum()), *c))


if __name__ == "__main__":
    globals()[sys.argv[1]]()
