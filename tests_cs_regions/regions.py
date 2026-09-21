#!/usr/bin/env python3
"""Build the per-region cubed-sphere error tables from the run logs.

Usage:  regions.py <logdir>

Parses the finaliser output of cs_test (CSTestConvErrors for iprob 9/13,
CSTestResistCheck for iprob 11) and prints, per test and per configuration:

  * L1 and Linf by REGION (panel INTERIOR / panel SEAM / CUBE VERTEX),
  * the measured convergence ORDER per region across the resolution sweep,
  * the vertex/interior and seam/interior error RATIOS,
  * max |v| / c_s for the at-rest atmosphere,
  * the final dt and the wall clock, so a "fix" that costs the time step shows.
"""
import math
import os
import re
import sys

REGIONS = ["panel INTERIOR", "panel SEAM", "CUBE VERTEX"]
RSHORT = {"panel INTERIOR": "INTERIOR", "panel SEAM": "SEAM", "CUBE VERTEX": "VERTEX"}
# c_s of the iprob-13 atmosphere: gamma = 5/3, p0 = 0.1, d0 = 1
CS_STRAT = math.sqrt(1.666667 * 0.1 / 1.0)

RE_CONV = re.compile(r"CS CONVERGENCE \(iprob=(\d+)\): nx2=(\d+) t=([0-9.eE+-]+)")
RE_MAX = re.compile(r"MAX norms over active cells:\s+Linf\(v\)=(\S+)\s+Linf\(p\)=(\S+)")
RE_REG = re.compile(
    r"REGION (panel INTERIOR|panel SEAM   |CUBE VERTEX  )\s*\(\s*(\d+) cells\)\s+"
    r"L1\(v\)=(\S+) L1\(p\)=(\S+) L1\(B\)=(\S+)\s+"
    r"Linf\(v\)=(\S+) Linf\(p\)=(\S+) Linf\(B\)=(\S+)")
RE_TAN = re.compile(
    r"REGION (panel INTERIOR|panel SEAM   |CUBE VERTEX  )\s+TANGENTIAL ONLY\s+"
    r"L1\(v_t\)=(\S+)\s+Linf\(v_t\)=(\S+)")
RE_FREG = re.compile(
    r"REGION (panel INTERIOR|panel SEAM   |CUBE VERTEX  )\s*\(\s*(\d+) faces\)\s+"
    r"L1=(\S+)\s+Linf=(\S+)")
RE_FTOT = re.compile(r"EVOLVED FIELD error vs the exact static B.*?L1=(\S+)\s+Linf=(\S+)")
# the per-cycle diagnostic line, and the termination block.  ndiag is large in the gate
# inputs, so the diagnostic line appears only at cycle 0 -- the FINAL cycle count and time
# come from the termination block instead, and the final dt from the last diagnostic.
RE_CYC = re.compile(r"cycle=(\d+)\s+time=(\S+)\s+dt=(\S+)")
RE_END = re.compile(r"^time=(\S+) cycle=(\d+)\s*$", re.M)
RE_WALL = re.compile(r"WALLCLOCK ([0-9.]+) s")
RE_TAGN = re.compile(r"^(\w+?)_(\w+?)_wb(\d)_n(\d+)(?:_t([0-9.]+))?(.*)\.log$")


def parse(path):
    d = {"regions": {}, "tan": {}, "freg": {}}
    with open(path, errors="replace") as f:
        txt = f.read()
    if "FATAL" in txt or "nan" in txt.lower():
        d["bad"] = True
    m = RE_CONV.search(txt)
    if m:
        d["iprob"], d["n"] = int(m.group(1)), int(m.group(2))
        d["time"] = float(m.group(3))
    m = RE_MAX.search(txt)
    if m:
        d["linf_v"], d["linf_p"] = float(m.group(1)), float(m.group(2))
    for m in RE_REG.finditer(txt):
        r = m.group(1).strip()
        d["regions"][r] = dict(ncell=int(m.group(2)), l1v=float(m.group(3)),
                               l1p=float(m.group(4)), l1b=float(m.group(5)),
                               liv=float(m.group(6)), lip=float(m.group(7)),
                               lib=float(m.group(8)))
    for m in RE_TAN.finditer(txt):
        d["tan"][m.group(1).strip()] = dict(l1vt=float(m.group(2)),
                                            livt=float(m.group(3)))
    for m in RE_FREG.finditer(txt):
        d["freg"][m.group(1).strip()] = dict(nface=int(m.group(2)),
                                             l1=float(m.group(3)),
                                             linf=float(m.group(4)))
    m = RE_FTOT.search(txt)
    if m:
        d["fl1"], d["flinf"] = float(m.group(1)), float(m.group(2))
    cyc = RE_CYC.findall(txt)
    if cyc:
        d["dt"] = float(cyc[-1][2])
    m = RE_END.search(txt)
    if m:
        d["tend"], d["ncycle"] = float(m.group(1)), int(m.group(2))
        # mean dt actually taken, which is what a "does this cost the time step?" question
        # is about; the cycle-0 dt above is only the start-up value
        if d["ncycle"] > 0:
            d["dt"] = d["tend"] / d["ncycle"]
    m = RE_WALL.search(txt)
    if m:
        d["wall"] = float(m.group(1))
    return d


def order(e1, e2):
    """order between successive resolutions (factor 2 refinement)"""
    if e1 is None or e2 is None or e1 <= 0.0 or e2 <= 0.0:
        return float("nan")
    return math.log(e1 / e2) / math.log(2.0)


def fmt(x, w=10):
    if x is None:
        return " " * w
    if isinstance(x, float) and (x != x):
        return "%*s" % (w, "-")
    return "%*.3e" % (w, x)


def main():
    logdir = sys.argv[1]
    runs = {}
    for fn in sorted(os.listdir(logdir)):
        m = RE_TAGN.match(fn)
        if not m:
            continue
        test, recon, wb, n, tlim, extra = m.groups()
        key = (test, recon, int(wb), tlim or "short", extra)
        runs.setdefault(key, {})[int(n)] = parse(os.path.join(logdir, fn))

    print("=" * 90)
    print("CUBED-SPHERE PER-REGION ERROR GATE")
    print("  tests: strat = iprob 13 hydrostatic atmosphere at rest (|v| IS the error)")
    print("         rot   = iprob 9  rigid rotation + uniform field (exact at all t)")
    print("         loop  = iprob 11 toroidal field from its vector potential (static)")
    print("  regions by angle; band half-width = 2/4/8 cells at n = 16/32/64")
    print("=" * 90)

    # ---- per-configuration resolution sweeps -----------------------------------------
    for key in sorted(runs, key=lambda k: (k[0], k[3], k[1], k[2])):
        test, recon, wb, tlim, extra = key
        res = runs[key]
        ns = sorted(res)
        head = "%s  recon=%-5s  cs_wellbalanced_src=%s  tlim=%s%s" % (
            test, recon, "ON " if wb else "OFF", tlim, extra)
        print("\n" + "-" * 90)
        print(head)
        d0 = res[ns[0]]
        if d0.get("bad"):
            print("  *** run reported FATAL / NaN")
        info = []
        for n in ns:
            d = res[n]
            info.append("n=%d: cyc=%s dtbar=%s t=%s wall=%ss" % (
                n, d.get("ncycle", "?"), ("%.3e" % d["dt"]) if "dt" in d else "?",
                ("%.4f" % d["tend"]) if "tend" in d else "?",
                ("%.1f" % d["wall"]) if "wall" in d else "?"))
        print("  " + " | ".join(info))

        if test == "loop":
            metrics = [("L1(B_face)", "l1", "freg"), ("Linf(B_face)", "linf", "freg")]
        else:
            metrics = [("L1(v)", "l1v", "regions"), ("Linf(v)", "liv", "regions"),
                       ("L1(v_tang)", "l1vt", "tan"), ("Linf(v_tang)", "livt", "tan"),
                       ("L1(p)", "l1p", "regions")]
            if any(res[n]["regions"].get("CUBE VERTEX", {}).get("l1b", 0.0) > 0
                   for n in ns):
                metrics += [("L1(B)", "l1b", "regions"), ("Linf(B)", "lib", "regions")]

        for label, field, where in metrics:
            print("  %-12s %s" % (label, "".join("%12s%7s" % ("n=%d" % n, "order")
                                                 for n in ns)))
            for r in REGIONS:
                vals = [res[n].get(where, {}).get(r, {}).get(field) for n in ns]
                row = "    %-9s" % RSHORT[r]
                for i, n in enumerate(ns):
                    o = order(vals[i - 1], vals[i]) if i > 0 else None
                    row += "%12s" % (("%.3e" % vals[i]) if vals[i] is not None else "-")
                    row += "%7s" % (("%.2f" % o) if (o is not None and o == o) else "")
                print(row)
            # vertex / interior ratio
            row = "    VTX/INT  "
            for i, n in enumerate(ns):
                vi = res[n].get(where, {}).get("panel INTERIOR", {}).get(field)
                vv = res[n].get(where, {}).get("CUBE VERTEX", {}).get(field)
                row += "%12s%7s" % (("%.2f" % (vv / vi)) if (vi and vv) else "-", "")
            print(row)
            row = "    SEAM/INT "
            for i, n in enumerate(ns):
                vi = res[n].get(where, {}).get("panel INTERIOR", {}).get(field)
                vs = res[n].get(where, {}).get("panel SEAM", {}).get(field)
                row += "%12s%7s" % (("%.2f" % (vs / vi)) if (vi and vs) else "-", "")
            print(row)

        if test == "strat":
            row = "  max|v|/c_s   "
            for n in ns:
                lv = res[n].get("linf_v")
                row += "%12s%7s" % (("%.3e" % (lv / CS_STRAT)) if lv else "-", "")
            print(row)

    # ---- cross-configuration comparison at each resolution ---------------------------
    print("\n" + "=" * 90)
    print("CONFIGURATION COMPARISON  (error relative to the plm / wb OFF baseline")
    print("at the same test and resolution; < 1 is an improvement)")
    print("=" * 90)
    for test in ("strat", "rot", "loop"):
        for n in (16, 32, 64):
            base = runs.get((test, "plm", 0, "short", ""), {}).get(n)
            if base is None:
                continue
            if test == "loop":
                mets = [("L1", "l1", "freg"), ("Linf", "linf", "freg")]
            else:
                mets = [("L1(v_t)", "l1vt", "tan"), ("Linf(v_t)", "livt", "tan"),
                        ("L1(v)", "l1v", "regions")]
            print("\n  %s, n=%d" % (test, n))
            print("    %-20s %s" % ("config", "".join(
                "%13s" % ("%s.%s" % (m[0], RSHORT[r][:3]))
                for m in mets for r in REGIONS)))
            for key in sorted(runs, key=lambda k: (k[1], k[2])):
                if key[0] != test or key[3] != "short" or key[4] != "":
                    continue
                d = runs[key].get(n)
                if d is None:
                    continue
                row = "    %-20s" % ("%s/wb%d" % (key[1], key[2]))
                for label, field, where in mets:
                    for r in REGIONS:
                        v = d.get(where, {}).get(r, {}).get(field)
                        b = base.get(where, {}).get(r, {}).get(field)
                        row += "%13s" % (("%.3f" % (v / b)) if (v and b) else "-")
                print(row)
    print()


if __name__ == "__main__":
    main()
