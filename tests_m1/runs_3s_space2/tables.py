#!/usr/bin/env python3
"""Tables of runs_3r_radwave from the results.json of every run.

SPACE: for each (arm, P_rad/P_gas, tau_lambda, N) the two runs at nt and 2 nt are
Richardson-extrapolated to dt -> 0 (omega_0 = 2 omega(2nt) - omega(nt); the same on the
error fields for L1), which removes the first-order time error (measured separately in
TIME).  err_w = |omega_0 - omega_ref|/|omega_ref| (complex), err_c = relative error of
the phase speed, err_g = error of the damping rate over |omega_ref|; the order is the
log2 ratio between successive N.  Explicit arms have no Richardson (dt = 0.4 dx/c).

TIME: at N = 64, omega(nt) for nt = 16 .. 8192; order = log2 of the ratio of successive
differences |omega(nt) - omega(2nt)|; err(nt) = |omega(nt) - omega_inf|/|omega_ref| with
omega_inf the Richardson value of the two finest nt.  Runs whose actual cycle count
exceeds nt were limited by the HYDRO CFL and are flagged 'h'.

  python3 tables.py [--runs DIR] [--json OUT]
"""

import argparse
import glob
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

PT = [(0.1, 0.1), (0.1, 10.0), (0.1, 1000.0), (1.0, 0.1), (1.0, 10.0), (1.0, 1000.0),
      (10.0, 0.1), (10.0, 10.0), (10.0, 1000.0), (100.0, 0.1), (100.0, 10.0),
      (100.0, 1000.0)]
AMP = 1.0e-5


def load(runs):
    out = []
    for p in glob.glob(os.path.join(runs, "*", "results.json")):
        with open(p) as fp:
            r = json.load(fp)
        r["path"] = os.path.dirname(p)
        out.append(r)
    return out


def arm_of(r):
    c = r["case"]
    ex = dict(c["extra"])
    a = "%s_%s_%s" % (c["clos"], c["trans"], c["dir"])
    for k in sorted(ex):
        a += "_%s%s" % (k, ex[k])
    return a


def wref(r):
    return complex(*r["refs"][r["ref_key"]])


def errfield(r, name):
    """exact-minus-numerical rho or E along the dumped line, per unit amplitude."""
    c = r["case"]
    kd = {"x1": (1, 0, 0), "x2": (0, 1, 0), "xy": (1, 1, 0), "xyz": (1, 1, 1)}[c["dir"]]
    kd = np.array(kd, dtype=float)
    lam = 1.0 / np.linalg.norm(kd)
    k = 2.0 * math.pi / lam
    kvec = k * kd / np.linalg.norm(kd)
    x1 = np.array(r["x1"])
    n = c["N"]
    x2s = 0.5 / n if c["dir"] in ("xy", "xyz") else (0.5 / 4 if kd[1] else 0.0)
    x3s = 0.5 / n if c["dir"] == "xyz" else 0.0
    ph = kvec[0] * x1 + kvec[1] * x2s + kvec[2] * x3s
    w = wref(r)
    ex = np.real(np.exp(1j * ph - 1j * w * r["t_end"]))
    if name == "rho":
        return (np.array(r["rho_end"]) - 1.0) / AMP - ex
    return None


def space_table(rs, fp):
    groups = defaultdict(dict)
    for r in rs:
        c = r["case"]
        if c["dir"] == "x1" and c["N"] == 64 and "ny" in c["extra"]:
            continue
        groups[arm_of(r)].setdefault((c["prat"], c["tau"], c["N"]), []).append(r)
    res = {}
    for arm in sorted(groups):
        g = groups[arm]
        rows = {}
        for key, lst in g.items():
            lst = sorted(lst, key=lambda x: x["case"]["nt"])
            if lst[0]["case"]["trans"] == "expl":
                r = lst[0]
                w0 = complex(r["w_re"], r["w_im"])
                e1 = errfield(r, "rho")
                l1 = float(np.mean(np.abs(e1)))
            else:
                # the pair (nt, 2 nt) with the largest nt
                pair = None
                for a in lst:
                    for b in lst:
                        if b["case"]["nt"] == 2 * a["case"]["nt"]:
                            pair = (a, b)
                if pair is None:
                    continue
                a, b = pair
                wa = complex(a["w_re"], a["w_im"])
                wb = complex(b["w_re"], b["w_im"])
                w0 = 2.0 * wb - wa
                e1 = 2.0 * errfield(b, "rho") - errfield(a, "rho")
                l1 = float(np.mean(np.abs(e1)))
                r = b
            wr = wref(r)
            rows[key] = {"w0": [w0.real, w0.imag], "wref": [wr.real, wr.imag],
                         "err_w": abs(w0 - wr) / abs(wr),
                         "err_c": (w0.real - wr.real) / wr.real,
                         "err_g": (w0.imag - wr.imag) / abs(wr),
                         "l1": l1, "wraw": [r["w_re"], r["w_im"]]}
        res[arm] = rows
        fp.write("\n== SPACE  arm %s  (dt -> 0 by Richardson; errors vs %s)\n"
                 % (arm, next(iter(g.values()))[0]["ref_key"]))
        fp.write("%6s %7s %5s %11s %11s %11s %6s %11s %6s\n"
                 % ("Pr/Pg", "tau", "N", "err_c", "err_g", "err_w", "p_w", "L1(rho)",
                    "p_L1"))
        for pr, tau in PT:
            ns = sorted(n for (a, b, n) in rows if (a, b) == (pr, tau))
            prev = None
            for n in ns:
                x = rows[(pr, tau, n)]
                pw = pl = ""
                if prev is not None and prev[0] * 2 == n:
                    pw = "%6.2f" % math.log2(prev[1]["err_w"] / x["err_w"])
                    pl = "%6.2f" % math.log2(prev[1]["l1"] / x["l1"])
                fp.write("%6g %7g %5d %11.3e %11.3e %11.3e %6s %11.3e %6s\n"
                         % (pr, tau, n, x["err_c"], x["err_g"], x["err_w"], pw,
                            x["l1"], pl))
                prev = (n, x)
    return res


def time_table(rs, fp):
    groups = defaultdict(dict)
    for r in rs:
        c = r["case"]
        if c["N"] != 64 or c["dir"] != "x1" or c["trans"] != "impl":
            continue
        if len([1 for x in rs if x["case"]["N"] == 64 and arm_of(x) == arm_of(r)
                and x["case"]["prat"] == c["prat"] and x["case"]["tau"] == c["tau"]]) < 5:
            continue
        groups[arm_of(r)].setdefault((c["prat"], c["tau"]), []).append(r)
    res = {}
    for arm in sorted(groups):
        fp.write("\n== TIME  arm %s  N=64  (err vs omega_inf = Richardson of the two "
                 "finest nt; p from successive differences; h = hydro-CFL limited)\n"
                 % arm)
        res[arm] = {}
        for (pr, tau), lst in sorted(groups[arm].items()):
            lst = sorted(lst, key=lambda x: x["case"]["nt"])
            ws = [complex(x["w_re"], x["w_im"]) for x in lst]
            nts = [x["case"]["nt"] for x in lst]
            winf = 2.0 * ws[-1] - ws[-2]
            wr = wref(lst[0])
            fp.write("%6g %7g |" % (pr, tau))
            row = []
            for i, x in enumerate(lst):
                e = abs(ws[i] - winf) / abs(wr)
                p = ""
                if 0 < i < len(ws) - 1:
                    d0 = abs(ws[i - 1] - ws[i])
                    d1 = abs(ws[i] - ws[i + 1])
                    if d1 > 0 and d0 > 0:
                        p = "%.2f" % math.log2(d0 / d1)
                h = "h" if (x.get("ncycle") or 0) > 1.01 * nts[i] else ""
                row.append({"nt": nts[i], "err": e, "p": p, "hydro": bool(h),
                            "w": [ws[i].real, ws[i].imag],
                            "dre": (ws[i].real - winf.real) / abs(wr),
                            "dim": (ws[i].imag - winf.imag) / abs(wr)})
                fp.write(" %d:%.2e%s(%s)" % (nts[i], e, h, p))
            fp.write("\n")
            res[arm]["%g_%g" % (pr, tau)] = {"rows": row, "w_inf": [winf.real, winf.imag],
                                             "wref": [wr.real, wr.imag]}
    return res


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runs", default="/viper/ptmp2/jinma/radwave_3r/runs")
    p.add_argument("--json", default=os.path.join(HERE, "tables.json"))
    p.add_argument("--txt", default=os.path.join(HERE, "RESULTS_tables.txt"))
    args = p.parse_args()
    rs = load(args.runs)
    with open(args.txt, "w") as fp:
        sp = space_table(rs, fp)
        tm = time_table(rs, fp)
    with open(args.json, "w") as fp:
        json.dump({"space": {a: {"%g_%g_%d" % k: v for k, v in d.items()}
                             for a, d in sp.items()}, "time": tm}, fp)
    print(open(args.txt).read())


if __name__ == "__main__":
    main()
