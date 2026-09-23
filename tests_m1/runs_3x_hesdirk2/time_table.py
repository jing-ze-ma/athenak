#!/usr/bin/env python3
"""Time-order table of radwave runs (runs_3r tables.py time_table): for each (P, tau) of
an arm, omega at nt = 128 ... 2048; err = |w(nt) - w_inf|/|w_ref| with w_inf the
Richardson value of the two finest nt for order p_R (2 for hesdirk2, 1 for be);
p = log2 of successive differences.   python3 time_table.py RUNS_DIR ARM_SUBSTRING p_R"""
import glob
import json
import math
import os
import sys
from collections import defaultdict

runs, sub, pr = sys.argv[1], sys.argv[2], float(sys.argv[3])
excl = [x for x in sys.argv[4].split(",") if x] if len(sys.argv) > 4 else []
g = defaultdict(list)
for f in glob.glob(os.path.join(runs, "*", "results.json")):
    tag = os.path.basename(os.path.dirname(f))
    if sub not in tag or any(x in tag for x in excl):
        continue
    if "tshesdirk2" in tag and "tshesdirk2" not in sub:
        continue
    if "nper" in tag or "amp" in tag or "tol" in tag or "_v0" in tag:
        continue
    r = json.load(open(f))
    c = r["case"]
    g[(c["prat"], c["tau"])].append(r)
for key in sorted(g):
    lst = sorted(g[key], key=lambda x: x["case"]["nt"])
    if len(lst) < 3:
        continue
    ws = [complex(x["w_re"], x["w_im"]) for x in lst]
    nts = [x["case"]["nt"] for x in lst]
    f = 2.0 ** pr
    winf = (f * ws[-1] - ws[-2]) / (f - 1.0)
    wr = complex(*lst[0]["refs"][lst[0]["ref_key"]])
    out = "%6g %7g |" % key
    ps = []
    for i in range(len(ws)):
        e = abs(ws[i] - winf) / abs(wr)
        p = ""
        if 0 < i < len(ws) - 1:
            d0, d1 = abs(ws[i - 1] - ws[i]), abs(ws[i] - ws[i + 1])
            if d0 > 0 and d1 > 0:
                p = "%.2f" % math.log2(d0 / d1)
                ps.append(float(p))
        h = "h" if (lst[i].get("ncycle") or 0) > 1.01 * nts[i] else ""
        out += " %d:%.2e%s(%s)" % (nts[i], e, h, p)
    out += "  | min p %.2f" % (min(ps) if ps else float("nan"))
    print(out)
