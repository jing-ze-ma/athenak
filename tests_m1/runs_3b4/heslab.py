#!/usr/bin/env python3
"""Analyse a 2-D He slab run of milestone 3b phase C and emit its JSON record.

Usage:  heslab.py <label> <rundir> [<label> <rundir> ...] --out <file.json>

Per run it collects, from the history file and the <rad_m1> `bin` dumps,
  * max |v1|, max |v2| in units of v_MLT, max |v|/c_s and min E
  * F1top/Fin, F1bot/Fin
  * the kinetic energy in x1 and x2 against time (the hydro history)
  * the horizontal-mean F_rad/F_in against depth at the last dump
  * a thinned map of v2 (or of E'/<E> when no velocity dump exists)
and keeps the file well under 200 kB by thinning the time series and the map.
"""

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import common  # noqa: E402

# vis/python/bin_convert.py splits a header line on "=" and unpacks two fields, which
# trips over the He slab's input file (its comments carry "=" as well).  Split on the
# FIRST "=" only; nothing else about the reader changes.
sys.path.insert(0, os.path.join(HERE, "..", "..", "vis", "python"))
import bin_convert  # noqa: E402


def _get_from_header(header, blockname, keyname):
    blockname = blockname.strip()
    keyname = keyname.strip()
    if not blockname.startswith("<"):
        blockname = "<" + blockname
    if blockname[-1] != ">":
        blockname += ">"
    block = "<none>"
    for line in [entry for entry in header]:
        if line.startswith("<"):
            block = line.split("#")[0].strip()
        if block != blockname:
            continue
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key.strip() == keyname:
            return value.split("#")[0].strip()
    raise KeyError("no parameter called %s/%s" % (blockname, keyname))


bin_convert.get_from_header = _get_from_header

VMLT = 1.86e4        # cm/s, the MLT velocity of the He FeCZ box
FIN = 2.475202e15    # code = cgs, the imposed bottom flux


def thin(a, n):
    a = np.asarray(a, dtype=float)
    if a.size <= n:
        return a.tolist()
    idx = np.unique(np.linspace(0, a.size - 1, n).astype(int))
    return a[idx].tolist()


def binfiles(d, tag):
    p = os.path.join(d, "bin")
    if not os.path.isdir(p):
        return []
    return sorted(os.path.join(p, f) for f in os.listdir(p)
                  if f.endswith(".bin") and tag in f)


def one(label, d):
    rec = {"label": label, "dir": os.path.basename(d)}
    # ---- the history: KE_1, KE_2 and the three <rad_m1> columns
    hh = [f for f in os.listdir(d) if f.endswith(".hydro.hst")]
    if hh:
        names, arr = common.read_history(os.path.join(d, hh[0]))
        t = common.hist_column(names, arr, "time")
        rec["t"] = thin(t, 400)
        for k, key in (("ke1", "1-KE"), ("ke2", "2-KE"), ("ke3", "3-KE")):
            try:
                rec[k] = thin(common.hist_column(names, arr, key), 400)
            except Exception:
                pass
        rec["dt"] = thin(common.hist_column(names, arr, "dt"), 400)
    hu = [f for f in os.listdir(d) if f.endswith(".user.hst")]
    if hu:
        names, arr = common.read_history(os.path.join(d, hu[0]))
        for k, key in (("F1top", "F1top"), ("F1bot", "F1bot"), ("V1max", "V1max")):
            rec[k] = thin(common.hist_column(names, arr, key), 400)
        rec["F1top_over_Fin_end"] = float(
            common.hist_column(names, arr, "F1top")[-1]) / FIN
        rec["F1bot_over_Fin_end"] = float(
            common.hist_column(names, arr, "F1bot")[-1]) / FIN
        rec["V1max_over_vMLT_end"] = float(
            common.hist_column(names, arr, "V1max")[-1]) / VMLT
        rec["t_end"] = float(common.hist_column(names, arr, "time")[-1])
    # ---- the last hydro dump: v1, v2, c_s proxy
    hw = binfiles(d, "hydro_w")
    if hw:
        du = common.load_dump(hw[-1])
        v1 = du.var("velx")
        v2 = du.var("vely")
        rec["maxv1_over_vMLT"] = float(np.max(np.abs(v1))) / VMLT
        rec["maxv2_over_vMLT"] = float(np.max(np.abs(v2))) / VMLT
        # a gamma = 5/3 sound-speed proxy from the dumped primitive p and rho
        if du.has("dens") and du.has("eint"):
            cs = np.sqrt(np.maximum((5.0 / 3.0) * (2.0 / 3.0)
                                    * du.var("eint") / du.var("dens"), 1e-300))
            rec["maxv_over_cs"] = float(np.max(np.sqrt(v1**2 + v2**2) / cs))
        # a thinned v2 map (x1 x x2), at most 84 x 32
        m = np.asarray(v2[0]) / VMLT
        rec["v2_map"] = [[round(float(x), 6) for x in row] for row in m]
        rec["map_x1"] = thin(du.x1, 200)
        rec["map_x2"] = thin(du.x2, 200)
        rec["map_time"] = du.time
    # ---- the last <rad_m1> dump: the horizontal-mean F_rad/F_in profile
    mb = binfiles(d, "m1")
    if mb:
        du = common.load_dump(mb[-1])
        f1 = np.asarray(du.var("m1_f1"))
        rec["z"] = [float(x) for x in du.x1]
        rec["Frad_over_Fin"] = [round(float(x), 8)
                                for x in (f1[0].mean(axis=0) / FIN)]
        e = np.asarray(du.var("m1_e"))
        rec["minE"] = float(np.min(e))
        if "v2_map" not in rec:
            eb = e[0].mean(axis=0)
            rec["Eprime_map"] = [[round(float(x), 6) for x in row]
                                 for row in (e[0] / eb - 1.0)]
        rec["m1_time"] = du.time
    # ---- the solver statistics, straight from the run log
    lg = os.path.join(d, "run.log")
    if os.path.exists(lg):
        txt = open(lg).read()
        for tag, key in (("implicit transport: solves=", "picard_line"),
                         ("bicgstab: outer passes=", "bicg_line"),
                         ("bicgstab: breakdowns=", "bicg_break_line"),
                         ("implicit transverse (", "lres_line")):
            i = txt.find(tag)
            if i >= 0:
                rec[key] = txt[i:txt.find("\n", i)].strip()
        rec["dt_collapses"] = txt.count("dt COLLAPSE")
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pairs", nargs="+")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = {}
    for i in range(0, len(a.pairs), 2):
        lab, d = a.pairs[i], a.pairs[i + 1]
        out[lab] = one(lab, d)
        print("%s: %s" % (lab, {k: v for k, v in out[lab].items()
                                if not isinstance(v, list)}))
    with open(a.out, "w") as f:
        json.dump(out, f)
    print("wrote %s (%.1f kB)" % (a.out, os.path.getsize(a.out) / 1024.0))


if __name__ == "__main__":
    main()
