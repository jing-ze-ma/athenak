#!/usr/bin/env python3
"""Write the profile data behind the milestone-1c-B results page as JSON.

No plotting library is involved: each file is a dict

    {"x": [...], "series": {"label": [...], ...}, "meta": {...}}

with at most ``--max-points`` samples per series (strided, peaks kept) and
every number in double precision, read from ``file_type = tab`` dumps.  The
beam file carries an extra ``"map"``, a coarsened 2-D E field.

Run it after tests_m1/runs_1cB/run_plots.sh, run_recon.sh, run_shear.sh and
the runs_1cB/step0 batch; it writes into tests_m1/plots/.
"""

import argparse
import json
import math
import os

import numpy as np

import common

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs_1cB")
RUNS_1C = os.path.join(HERE, "runs_1c")
OUT = os.path.join(HERE, "plots")


# ---------------------------------------------------------------------------
def prof(path, name="m1_e"):
    """(x, q) of a 1-D tab (or bin) dump."""
    d = common.load_dump(path)
    x, q = common.extract_1d(d, name, axis=1)
    return np.asarray(x, float), np.asarray(q, float)


def dump_time(path):
    return common.load_dump(path).time


def thin(x, ys, maxn):
    """Stride x and every series in ys down to <= maxn points, keeping the
    global extrema of the FIRST series so that peaks survive."""
    n = x.size
    if n <= maxn:
        return x, ys
    step = int(math.ceil(n / float(maxn)))
    keep = set(range(0, n, step))
    keep.add(n - 1)
    ref = ys[next(iter(ys))]
    keep.add(int(np.argmax(ref)))
    keep.add(int(np.argmin(ref)))
    idx = np.array(sorted(keep))
    return x[idx], {k: v[idx] for k, v in ys.items()}


def write(name, x, series, meta, extra=None, maxn=260):
    x, series = thin(np.asarray(x, float),
                     {k: np.asarray(v, float) for k, v in series.items()},
                     maxn)
    doc = {"x": [float(v) for v in x],
           "series": {k: [float(v) for v in a] for k, a in series.items()},
           "meta": meta}
    if extra:
        doc.update(extra)
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name + ".json")
    with open(path, "w") as fp:
        json.dump(doc, fp, separators=(",", ":"))
    print("%-16s %7.1f kB  %d points" % (name, os.path.getsize(path) / 1024.0,
                                         len(doc["x"])))


# ---------------------------------------------------------------------------
def beam2d(maxn):
    """(i) A cut of E across the 45-degree beam at mid-domain, plus a
    coarsened map.  The exact solution is a top-hat of the beam's own width
    measured perpendicular to it: the inflow patch spans beam_y1 - beam_y0
    in y, so 0.125/sqrt(2) across a 45-degree beam."""
    d = common.load_dump(os.path.join(RUNS, "beam_plm", "bin",
                                      "m1_beam.m1.00006.bin"))
    x1 = np.asarray(d.x1, float)
    x2 = np.asarray(d.x2, float)
    e = np.asarray(d.var("m1_e"), float)
    e = e.reshape(e.shape[-2], e.shape[-1])          # (x2, x1)
    # cut perpendicular to the beam through the mid-point of the domain, i.e.
    # the anti-diagonal s -> (0.5 + s/sqrt2, 0.5 - s/sqrt2)
    dx = x1[1] - x1[0]
    ns = 121
    s = (np.arange(ns) - (ns - 1) / 2.0) * dx
    px = 0.5 + s / math.sqrt(2.0)
    py = 0.5 - s / math.sqrt(2.0)
    cut = common.sample_2d(x1, x2, e, px, py)
    # Along this cut y - x = -sqrt(2) s, so the signed perpendicular distance
    # from the beam axis y = x is u = (y - x)/sqrt(2) = -s.  The inflow patch
    # y in [y0, y1] at x = 0 fills u in [y0/sqrt2, y1/sqrt2], i.e. the exact
    # free-streaming beam is the top-hat |s + (y0+y1)/(2 sqrt2)| <=
    # (y1-y0)/(2 sqrt2).
    y0, y1 = 0.0, 0.125
    s0 = -0.5 * (y0 + y1) / math.sqrt(2.0)
    half = 0.5 * (y1 - y0) / math.sqrt(2.0)
    exact = np.where(np.abs(s - s0) <= half, 1.0, 0.0)
    # the coarse map: 128^2 -> 64^2 block average
    m = e.reshape(64, 2, 64, 2).mean(axis=(1, 3))
    meta = {"test": "T1 beam", "t": d.time, "nx": int(e.shape[1]),
            "cut": "perpendicular to the beam through the domain centre",
            "x_units": "cells along the cut", "dx": float(dx),
            "map_extent": [float(x1[0]), float(x1[-1]),
                           float(x2[0]), float(x2[-1])],
            "map_shape": [64, 64], "scheme": "plm, thick_flux=none"}
    write("beam2d", s / dx, {"E": cut, "exact_tophat": exact}, meta,
          extra={"map": [[float(v) for v in row] for row in m]}, maxn=maxn)


def pulse1d(maxn):
    """(ii) Free-streaming pulse at t=0 and after one crossing."""
    series, meta = {}, {"test": "T1 pulse1d", "exact": "the initial profile"}
    x = None
    for tag, rec, n in [("plm_64", "plm", 64), ("plm_256", "plm", 256),
                        ("wenoz_64", "wenoz", 64), ("ppm4_64", "ppm4", 64)]:
        d = os.path.join(RUNS, "p1d_%s_%d" % (rec, n), "tab")
        xi, e0 = prof(os.path.join(d, "m1_pulse1d.m1.00000.tab"))
        _, e1 = prof(os.path.join(d, "m1_pulse1d.m1.00001.tab"))
        if n == 256:
            # put the 256-cell run on the 64-cell abscissa by block averaging
            e0 = e0.reshape(64, 4).mean(axis=1)
            e1 = e1.reshape(64, 4).mean(axis=1)
        if x is None:
            x = xi
        series["init_" + tag] = e0
        series["final_" + tag] = e1
    meta["note"] = ("the 256-cell run is block averaged onto the 64-cell "
                    "abscissa so that all five curves share one x")
    write("pulse1d", x, series, meta, maxn=maxn)


def thick_pulse(maxn):
    """(iii) T3 at tau_cell = 1e3 for the three thick_flux choices."""
    series = {}
    x = None
    t = None
    for tf in ["ap_hll", "scaled", "none"]:
        p = os.path.join(RUNS, "pl_t3_%s" % tf, "tab")
        f = sorted(os.listdir(p))[-1]
        xi, e = prof(os.path.join(p, f))
        t = dump_time(os.path.join(p, f))
        x = xi if x is None else x
        series[tf] = e
    # analytic: a Gaussian of variance w^2/2 spreading at 2D, D = c/(3 rho k)
    amp, w, x0, bg = 1.0, 0.05, 0.5, 1.0e-3
    dcoef = 1.0 / (3.0 * 1.0 * 128000.0)
    s0sq = 0.5 * w * w
    ssq = s0sq + 2.0 * dcoef * t
    series["analytic"] = bg + amp * math.sqrt(s0sq / ssq) * \
        np.exp(-(x - x0) ** 2 / (2.0 * ssq))
    meta = {"test": "T3 thick pulse", "tau_cell": 1000.0, "t": t,
            "D": dcoef, "scheme": "plm, ap_form=alpha2",
            "analytic": "Gaussian spreading at 2D"}
    write("thick_pulse", x, series, meta, maxn=maxn)


def tophat(maxn):
    """(iv) T3c clipped top-hat at tau_cell = 1e4."""
    series = {}
    x = None
    t = None
    for tf in ["ap_hll", "none"]:
        p = os.path.join(RUNS, "pl_tophat_%s" % tf, "tab")
        f = sorted(os.listdir(p))[-1]
        xi, e = prof(os.path.join(p, f))
        t = dump_time(os.path.join(p, f))
        x = xi if x is None else x
        series[tf] = e
    amp, w, x0, bg = 1.0, 0.05, 0.5, 1.0e-3
    dcoef = 1.0 / (3.0 * 1.0 * 1280000.0)
    sq = math.sqrt(4.0 * dcoef * t)
    from math import erf
    ev = np.array([0.5 * amp * (erf((x0 + w - xx) / sq) +
                                erf((xx - x0 + w) / sq)) for xx in x])
    series["analytic"] = bg + ev
    meta = {"test": "T3c clipped top-hat", "tau_cell": 10000.0, "t": t,
            "D": dcoef, "scheme": "plm",
            "analytic": "diffusing top-hat (error functions)"}
    write("tophat", x, series, meta, maxn=maxn)


def t4_advect(maxn):
    """(v) T4 dynamic at 512 cells: split on, split off, static shifted."""
    base = os.path.join(RUNS, "step0")
    xs, es = prof(os.path.join(base, "t4_s512", "tab",
                               "m1_advect_pulse.m1.00004.tab"))
    # split_on comes from the 1c-B default scheme (split_vel = recon)
    x, e_on = prof(os.path.join(RUNS, "recheck", "t4_recon", "tab",
                                "m1_advect_pulse.m1.00004.tab"))
    _, e_off = prof(os.path.join(base, "t4_d512_nosplit", "tab",
                                 "m1_advect_pulse.m1.00004.tab"))
    t = dump_time(os.path.join(RUNS, "recheck", "t4_recon", "tab",
                               "m1_advect_pulse.m1.00004.tab"))
    v = 3.0e7
    length = xs[-1] - xs[0] + (xs[1] - xs[0])
    xq = (xs - v * t - xs[0]) % length + xs[0]
    order = np.argsort(xq)
    shifted = np.interp(x, xq[order], es[order], period=length)
    meta = {"test": "T4 advected thick pulse (dynamic)", "nx": int(x.size),
            "t": t, "v": v, "beta_tau": 14.0,
            "scheme": "plm, ap_hll, ap_form=alpha2, split_vel=recon",
            "static_shifted": "the v = 0 run translated by v t"}
    write("t4_advect", x, {"split_on": e_on, "split_off": e_off,
                           "static_shifted": shifted}, meta, maxn=maxn)


def marshak(maxn):
    """(vi) T6 at 64 cells: E and the material energy vs the S_N reference."""
    series = {}
    x = None
    t = None
    for tf in ["ap_hll", "none"]:
        p = os.path.join(RUNS, "pl_marshak_%s" % tf, "tab")
        fm = [f for f in sorted(os.listdir(p)) if ".m1." in f][-1]
        fh = [f for f in sorted(os.listdir(p)) if ".hydro_w." in f][-1]
        xi, e = prof(os.path.join(p, fm))
        _, ei = prof(os.path.join(p, fh), "eint")
        t = dump_time(os.path.join(p, fm))
        x = xi if x is None else x
        series["E_" + tf] = e
        # the reference's material column is a_r T_mat^4, so convert the code's
        # internal energy density the same way: T = e_int/(rho c_v), c_v = 1.5
        series["material_" + tf] = 1.0e30 * (ei / 1.5) ** 4
    ref = np.loadtxt(os.path.join(RUNS_1C, "t6_ref_sn.txt"))
    series["E_reference"] = np.interp(x, ref[:, 0], ref[:, 1])
    series["material_reference"] = np.interp(x, ref[:, 0], ref[:, 2])
    meta = {"test": "T6 constant-c_v Marshak wave", "nx": int(x.size),
            "t": t, "scheme": "plm, 64 cells",
            "reference": "tests_m1/runs_1c/t6_ref_sn.txt (S_N, same c_v)",
            "material": "a_r T_mat^4, with T = e_int/(rho c_v), c_v = 1.5, a_r = 1e30"}
    write("marshak", x, series, meta, maxn=maxn)


def t3b(maxn):
    """(vii) T3b case A: E near the jump, and the face flux."""
    p = os.path.join(RUNS, "pl_t3b_ap_hll", "tab")
    f = sorted(os.listdir(p))[-1]
    x, e = prof(os.path.join(p, f))
    t = dump_time(os.path.join(p, f))
    c, kf, flux, xj, ratio, e_l = 1.0, 64.0, 1.0e-6, 0.5, 1.0e3, 1.0
    xl = 0.0
    s_l = -3.0 * 1.0 * kf * flux / c
    s_r = ratio * s_l
    ej = e_l + s_l * (xj - xl)
    exact = np.where(x < xj, e_l + s_l * (x - xl), ej + s_r * (x - xj))
    # face flux -c dE/dx/(3 rho kappa_F) with the arithmetic-mean face opacity
    dx = x[1] - x[0]
    rk = np.where(x < xj, kf, ratio * kf)
    rkf = 0.5 * (rk[:-1] + rk[1:])
    fface = -c * (e[1:] - e[:-1]) / (3.0 * dx * rkf)
    xface = 0.5 * (x[1:] + x[:-1])
    ij = int(np.argmin(np.abs(x - xj)))
    lo, hi = max(0, ij - 20), min(x.size, ij + 20)
    write("t3b", x[lo:hi], {"E": e[lo:hi], "analytic": exact[lo:hi]},
          {"test": "T3b opacity jump, case A", "tau_cell": [1.0, 1000.0],
           "t": t, "x_jump": xj, "imposed_flux": flux,
           "scheme": "plm, ap_hll, ap_form=alpha2",
           "analytic": "two-slope constant-flux steady state"}, maxn=maxn)
    lo, hi = max(0, ij - 20), min(xface.size, ij + 20)
    write("t3b_flux", xface[lo:hi],
          {"face_flux": fface[lo:hi],
           "imposed": np.full(hi - lo, flux)},
          {"test": "T3b opacity jump, case A: the face flux",
           "t": t, "x_jump": xj, "imposed_flux": flux,
           "definition": "-c (E_{i+1}-E_i)/(3 dx <rho kappa_F>_face)"},
          maxn=maxn)


def convergence(maxn):
    """(viii) The two convergence tables."""
    import subprocess
    tbl = {}
    for rec in ["plm", "ppm4", "wenoz"]:
        l1 = []
        for n in [32, 64, 128, 256, 512]:
            d = os.path.join(RUNS, "p1d_%s_%d" % (rec, n), "tab")
            x, e0 = prof(os.path.join(d, "m1_pulse1d.m1.00000.tab"))
            _, e1 = prof(os.path.join(d, "m1_pulse1d.m1.00001.tab"))
            l1.append(float(np.abs(e1 - e0).sum() / np.abs(e0).sum()))
        tbl[rec] = l1
    meta = {"test": "pulse1d L1 vs N", "N": [32, 64, 128, 256, 512],
            "quantity": "relative L1 error of E after one box crossing"}
    write("conv_pulse1d", [32, 64, 128, 256, 512], tbl, meta, maxn=maxn)

    tbl2 = {}
    for rec in ["plm", "wenoz"]:
        errs = []
        for n in [64, 128, 256, 512]:
            k = 1280.0 * n / 128.0
            d = os.path.join(RUNS, "t3n_%s_%d" % (rec, n), "bin")
            f = sorted(os.listdir(d))
            out = subprocess.run(
                ["python3", os.path.join(HERE, "t3_pulse.py")] +
                [os.path.join(d, ff) for ff in f] +
                ["--c", "1", "--rho", "1", "--kappa", repr(k),
                 "--subtract-min", "--quiet"],
                capture_output=True, text=True).stdout
            r = float(out.split("ratio=")[1].split()[0])
            errs.append(abs(r - 1.0))
        tbl2[rec] = errs
    meta2 = {"test": "T3 diffusion-rate error vs N at tau_cell = 10",
             "N": [64, 128, 256, 512],
             "quantity": "|measured/exact - 1| of d(sigma^2)/dt"}
    write("conv_t3_rate", [64, 128, 256, 512], tbl2, meta2, maxn=maxn)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--max-points", type=int, default=260)
    p.add_argument("--only", default=None, help="run one section only")
    args = p.parse_args()
    todo = [("beam2d", beam2d), ("pulse1d", pulse1d),
            ("thick_pulse", thick_pulse), ("tophat", tophat),
            ("t4_advect", t4_advect), ("marshak", marshak), ("t3b", t3b),
            ("convergence", convergence)]
    for name, fn in todo:
        if args.only and args.only != name:
            continue
        fn(args.max_points)


if __name__ == "__main__":
    main()
