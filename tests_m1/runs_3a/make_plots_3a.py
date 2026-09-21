#!/usr/bin/env python3
"""Milestone 3a: write tests_m1/plots/impl_*.json from the run directories of
runs_3a (see RESULTS.txt).  Same format as make_plots.py: one
{"x": [...], "series": {...}, "meta": {...}} per figure, double precision,
no plotting library, <= 260 points per series."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(HERE)
RUNS = os.path.join(TESTS, "tests_out_pl")
OUT = os.path.join(TESTS, "plots")
MAXN = 260


def prof(path, var="m1_e"):
    d = common.load_dump(path)
    return common.extract_1d(d, var)


def thin(x, n=MAXN):
    x = np.asarray(x, dtype=float)
    if x.size <= n:
        return x
    idx = np.linspace(0, x.size - 1, n).round().astype(int)
    return x[idx]


def write(name, x, series, meta):
    os.makedirs(OUT, exist_ok=True)
    obj = {"x": [float(v) for v in thin(x)],
           "series": {k: [float(v) for v in thin(np.asarray(s))]
                      for k, s in series.items()},
           "meta": meta}
    p = os.path.join(OUT, "impl_%s.json" % name)
    with open(p, "w") as fp:
        json.dump(obj, fp)
    print("%-28s %6.1f kB" % (os.path.basename(p), os.path.getsize(p) / 1024.0))


def thick_pulse():
    """(i) the thick pulse at tau_cell = 1e3 after t = 12, vs the analytic Gaussian."""
    series = {}
    x = None
    for tag, d in [("implicit_cfl1", "tp_c1"), ("implicit_cfl100", "tp_c100"),
                   ("implicit_cfl1e4", "tp_c10000"), ("explicit", "tp_exp")]:
        p = os.path.join(RUNS, d, "tab")
        f = [q for q in sorted(os.listdir(p)) if ".m1." in q][-1]
        xi, e = prof(os.path.join(p, f))
        x = xi if x is None else x
        series[tag] = e
    # analytic: a Gaussian of variance sigma0^2 + 2 D t on a background
    dd = 1.0 / (3.0 * 128000.0)
    s2 = 0.05 ** 2 + 2.0 * dd * 12.0
    amp = 1.0 * 0.05 / np.sqrt(s2)
    series["analytic"] = 1.0e-3 + amp * np.exp(-(x - 0.5) ** 2 / (2.0 * s2))
    write("thick_pulse", x, series,
          {"test": "I1 thick pulse, tau_cell = 1e3, t = 12, 128 cells",
           "D": dd, "analytic": "Gaussian, variance sigma0^2 + 2 D t, bg 1e-3",
           "note": "implicit curves at radiation CFL c dt/dx = 1, 100, 1e4"})


def t3b():
    """(ii) T3b opacity jump: E(x) explicit vs implicit vs the analytic steady state."""
    series = {}
    x = None
    for tag, d in [("explicit", "jump_exp"), ("implicit_cfl10", "jump_imp")]:
        p = os.path.join(RUNS, d, "tab")
        f = [q for q in sorted(os.listdir(p)) if ".m1." in q][-1]
        xi, e = prof(os.path.join(p, f))
        x = xi if x is None else x
        series[tag] = e
    flux, kap, ratio, xj = 1.0e-3, 0.64, 1000.0, 0.5
    sl = -3.0 * kap * flux
    ex = np.where(x < xj, 1.0 + sl * x, 1.0 + sl * xj + ratio * sl * (x - xj))
    series["analytic"] = ex
    write("t3b_jump", x, series,
          {"test": "I2 T3b opacity jump (x1000 at x = 0.5), steady state, 64 cells",
           "analytic": "dE/dx = -3 rho kappa F/c, E(0) = 1, F = 1e-3"})


def marshak():
    """(iii) T6 Marshak at 128 cells, implicit at CFL 1/10/100, vs the S_N table."""
    series = {}
    x = None
    for tag, d in [("explicit", "ma_exp"), ("implicit_cfl1", "ma_c1"),
                   ("implicit_cfl10", "ma_c10"), ("implicit_cfl100", "ma_c100")]:
        p = os.path.join(RUNS, d, "tab")
        f = [q for q in sorted(os.listdir(p)) if ".m1." in q][-1]
        xi, e = prof(os.path.join(p, f))
        x = xi if x is None else x
        series[tag] = e
    ref = np.loadtxt(os.path.join(HERE, "t6_ref_sn.txt"))
    series["reference_SN"] = np.interp(x, ref[:, 0], ref[:, 1])
    write("marshak", x, series,
          {"test": "I4 constant-c_v Marshak wave, 128 cells, t = 202",
           "reference": "tests_m1/runs_3a/t6_ref_sn.txt (S_N, same c_v)",
           "bc": "implicit: F_f = -c q (E - E_bath), q = 0.5, E_bath = 1e-10"})


if __name__ == "__main__":
    thick_pulse()
    t3b()
    marshak()
