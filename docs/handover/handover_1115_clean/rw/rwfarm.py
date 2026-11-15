#!/usr/bin/env python3
"""Farm front end of the radiation-wave ladder driver rw/run_radwave.py.

The decks, references and analysis are those of run_radwave.py; this file only splits
"write deck + run + analyse" into steps that can run on different machines.

  rwfarm.py emit  CASELIST ROOT [--rel]  > runlist   write ROOT/<tag>/rw.athinput, print one
                                                     run command per case ($ATHENA = binary)
  rwfarm.py score CASELIST ROOT                      analyse finished runs -> <tag>/farm.json
                                                     (+ conservation of mass, energy, momentum)
  rwfarm.py table CASELIST ROOT [--self]             space / time ladder table with orders

--rel prints commands relative to the directory holding ROOT (packages for Orion).
"""

import argparse
import glob
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, "/viper/ptmp2/jinma/clean_1114/rw")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "rw"))
import run_radwave as rr  # noqa: E402
import common  # noqa: E402,F401  (run_radwave.line imports it)


def cases(path):
    with open(path) as fp:
        return [s.strip() for s in fp if s.strip() and not s.startswith("#")]


def setup(line):
    c = rr.parse_case(line)
    bg, refs = rr.references(c)
    rk = rr.ref_key(c)
    w, vec = refs[rk]
    return c, bg, refs, rk, w, vec


def emit(args):
    root = os.path.abspath(args.root)
    for line in cases(args.caselist):
        c, bg, refs, rk, w, vec = setup(line)
        tag = rr.case_tag(c)
        rdir = os.path.join(root, tag)
        os.makedirs(rdir, exist_ok=True)
        with open(os.path.join(rdir, "rw.athinput"), "w") as fp:
            fp.write(rr.make_input(c, bg, w, vec)[0])
        d = os.path.relpath(rdir, os.path.dirname(root)) if args.rel else rdir
        print('cd %s && nice -n 5 "$ATHENA" -i rw.athinput -d . > log.txt 2>&1; '
              'echo "rc=$?" >> log.txt' % d)


def conservation(rdir, nx2, nx3):
    """Relative change of the box totals between the first and the last dump:
    mass, total energy (gas + radiation), momentum rho v + F/c^2 (per AMP rho0)."""
    def tot(i):
        h = common.load_tab(sorted(glob.glob(os.path.join(rdir, "tab", "rw.hydro_w.*.tab")))[i])
        m = common.load_tab(sorted(glob.glob(os.path.join(rdir, "tab", "rw.m1.*.tab")))[i])
        rho, v, e = h.var("dens"), h.var("velx"), h.var("eint")
        return (rho.mean(), (e + 0.5 * rho * v * v + m.var("m1_e")).mean(),
                (rho * v + m.var("m1_f1") / rr.C_LIGHT ** 2).mean())
    a, b = tot(0), tot(-1)
    return {"dmass": (b[0] - a[0]) / a[0], "denergy": (b[1] - a[1]) / a[1],
            "dmom": b[2] - a[2]}


def score(args):
    root = os.path.abspath(args.root)
    for line in cases(args.caselist):
        c, bg, refs, rk, w, vec = setup(line)
        tag = rr.case_tag(c)
        rdir = os.path.join(root, tag)
        log = os.path.join(rdir, "log.txt")
        done = os.path.exists(os.path.join(rdir, "results.json")) or (
            os.path.exists(log) and "rc=0" in open(log).read())
        if not done:
            print("%-40s MISSING/FAILED" % tag[:40])
            continue
        _, period, dt, icfl, nx2, nx3 = rr.make_input(c, bg, w, vec)
        out = rr.analyse(rdir, c, bg, w, vec, nx2, nx3)
        amp = float(c["extra"].get("amp", rr.AMP))
        cons = conservation(rdir, nx2, nx3)
        cons["dmom"] /= amp * bg.rho0
        tail = open(log).read().splitlines()
        out.update(case=c, ref_key=rk, period=period, dt_target=dt, implicit_cfl=icfl,
                   refs={k: [refs[k][0].real, refs[k][0].imag] for k in refs},
                   cons=cons,
                   clean_last=next((s for s in reversed(tail) if s.startswith("CLEAN")), ""))
        with open(os.path.join(rdir, "farm.json"), "w") as fp:
            json.dump(out, fp)
        print("%-40s l1_rho %.3e l1_E %.3e dM %.1e dE %.1e dP %.1e" % (
            tag[:40], out["l1_rho"], out["l1_E"], cons["dmass"], cons["denergy"],
            cons["dmom"]))


def order(a, b, r=2.0):
    return math.log(a / b) / math.log(r) if a > 0 and b > 0 else float("nan")


def table(args):
    """Groups = case lines differing only in N (space) or only in nt (time)."""
    root = os.path.abspath(args.root)
    groups = {}
    for line in cases(args.caselist):
        c = rr.parse_case(line)
        p = os.path.join(root, rr.case_tag(c), "farm.json")
        if not os.path.exists(p):
            continue
        r = json.load(open(p))
        lay = c["extra"].get("rm_clean_layout", "staggered")
        amp = c["extra"].get("amp", "1e-5")
        space = c["nt"] == 2 * c["N"]
        key = ("space" if space else "time", lay, amp, c["prat"], c["tau"])
        groups.setdefault(key, []).append((c["N"] if space else c["nt"], r))
    for key in sorted(groups):
        lv = sorted(groups[key], key=lambda t: t[0])
        kind, lay, amp, prat, tau = key
        print("\n### %s ladder, layout %s, amp %s, P %.3g, tau %g\n" % (kind, lay, amp, prat, tau))
        hdr = "| %s | L1 rho | ord | L1 E | ord | self dRho | ord | dE_tot | dP |" % (
            "N" if kind == "space" else "nt")
        print(hdr + "\n" + "|---" * 9 + "|")
        prev, dself = None, []
        for q, (n, r) in enumerate(lv):
            ds = float("nan")
            if q + 1 < len(lv) and kind == "time":    # same mesh: profile differences
                ds = float(np.mean(np.abs(np.array(r["rho_end"]) - np.array(lv[q + 1][1]["rho_end"]))))
                ds /= float(r["case"]["extra"].get("amp", rr.AMP))
            dself.append(ds)
            o_r = order(prev["l1_rho"], r["l1_rho"]) if prev else float("nan")
            o_e = order(prev["l1_E"], r["l1_E"]) if prev else float("nan")
            o_s = order(dself[-2], ds) if q > 0 else float("nan")
            print("| %d | %.3e | %.2f | %.3e | %.2f | %.3e | %.2f | %.1e | %.1e |" % (
                n, r["l1_rho"], o_r, r["l1_E"], o_e, ds, o_s, r["cons"]["denergy"],
                r["cons"]["dmom"]))
            prev = r


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=("emit", "score", "table"))
    p.add_argument("caselist")
    p.add_argument("root")
    p.add_argument("--rel", action="store_true")
    a = p.parse_args()
    {"emit": emit, "score": score, "table": table}[a.mode](a)


if __name__ == "__main__":
    main()
