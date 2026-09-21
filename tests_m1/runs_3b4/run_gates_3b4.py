#!/usr/bin/env python3
"""Milestone 3b phase C gates: <rad_m1>/implicit_solver = bicgstab against line_jacobi.

Every case runs both solvers on the SAME command line, is analysed and has its dumps
deleted immediately (the viper inode quota is at its limit).  Run from tests_m1/runs_3b4
with no arguments, or with the gate names as arguments (g1 g2 g4 g5 g6).
"""

import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import common  # noqa: E402,F401
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
XPL = os.path.join(ROOT, "build_cpu_m1", "src", "athena")
INP = os.path.join(HERE, "pulse_md.athinput")
# the 3-D copy: <output3> is a formatted table, which needs BOTH slice planes in 3-D
# and refuses the extra one in 2-D, so the two meshes need two input files
INP3 = os.path.join(HERE, "pulse_md3.athinput")
PULSE = os.path.join(ROOT, "tests_m1", "runs_3b3", "pulse_md.py")
OUT = os.environ.get("M1_3B4_OUT", os.path.join(HERE, "out"))
SOLVERS = ("line_jacobi", "bicgstab")
TABLE = []


def run(name, args, tlim, odt, extra=(), inp=None):
    d = os.path.join(OUT, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    cmd = ([XPL, "-i", inp or INP, "time/cfl_number=1e6", "rad_m1/transport=implicit",
            "time/tlim=%r" % tlim, "output1/dt=%r" % odt] + list(args) + list(extra))
    t0 = time.time()
    with open(os.path.join(d, "run.log"), "w") as f:
        rc = subprocess.call(cmd, cwd=d, stdout=f, stderr=subprocess.STDOUT)
    return d, rc, time.time() - t0


def stats(d):
    txt = open(os.path.join(d, "run.log")).read()
    out = {"inner": "-", "imax": "-", "brk": "-", "red": "-", "ncyc": "0"}
    i = txt.find("implicit transport: solves=")
    if i >= 0:
        seg = txt[i:i + 220]
        out["picard"] = seg.split("iterations mean=")[1].split()[0]
        out["pmax"] = seg.split("max=")[1].split()[0]
    i = txt.find("bicgstab: outer passes=")
    if i >= 0:
        seg = txt[i:i + 400]
        out["inner"] = seg.split("inner iterations mean=")[1].split()[0]
        out["imax"] = seg.split("max=")[1].split()[0]
    i = txt.find("bicgstab: breakdowns=")
    if i >= 0:
        seg = txt[i:i + 300]
        out["brk"] = seg.split("breakdowns=")[1].split()[0]
        out["red"] = seg.split("per inner iteration")[0].split("(")[-1].strip()
    i = txt.find("cycles = ")
    if i >= 0:
        out["ncyc"] = txt[i + 9:].split()[0]
    return out


def analyse(d, kappa, label, tol=0.02):
    files = sorted(os.path.join(d, "bin", f)
                   for f in os.listdir(os.path.join(d, "bin")))
    p = subprocess.run([sys.executable, PULSE] + files +
                       ["--kappa", repr(kappa), "--label", label, "--tol", repr(tol)],
                       capture_output=True, text=True)
    return p.stdout.strip() + p.stderr.strip()


def pair(log, name, args, tlim, odt, kappa, tol=0.02, inp=None):
    """run the same case with both solvers and print the two lines side by side"""
    row = {"case": name}
    for s in SOLVERS:
        d, rc, wall = run("%s_%s" % (name, s), args + ["rad_m1/implicit_solver=%s" % s],
                          tlim, odt, inp=inp)
        st = stats(d)
        txt = analyse(d, kappa, "%s/%s" % (name, s), tol=tol) if rc == 0 else "RUN FAILED"
        ncyc = max(int(float(st["ncyc"])), 1)
        log("%-58s %-11s outer=%s/%s inner=%s/%s red/it=%s brk=%s  %.3f s/step rc=%d" %
            (txt, s, st.get("picard"), st.get("pmax"), st["inner"], st["imax"],
             st["red"], st["brk"], wall / ncyc, rc))
        row[s] = {"outer": st.get("picard"), "outer_max": st.get("pmax"),
                  "inner": st["inner"], "inner_max": st["imax"],
                  "red_per_it": st["red"], "breakdowns": st["brk"],
                  "wall_per_step": wall / ncyc, "nstep": ncyc, "rate": txt}
        shutil.rmtree(os.path.join(d, "bin"), ignore_errors=True)
    TABLE.append(row)


def g1(log):
    log("\n=== G1  isotropy and rate, bicgstab vs line_jacobi ===")
    for nd, nx in ((2, 64), (3, 32)):
        for tau in (10.0, 1.0e3):
            kappa = tau * nx
            diff = 1.0 / (3.0 * kappa)
            for cfl in (1.0, 1.0e2, 1.0e4):
                dt = cfl / nx
                tlim = max(3.1e-4 * kappa, 12.0 * dt)
                name = "g1_%dd_tau%g_cfl%g" % (nd, tau, cfl)
                if 2.0 * diff * tlim > 0.02:
                    log("%-58s NOT MEASURABLE (the diffusion length of 12 steps "
                        "exceeds the box)" % name)
                    continue
                args = ["rad_m1/kappa_s=%r" % kappa, "rad_m1/implicit_cfl=%r" % cfl,
                        "mesh/nx1=%d" % nx, "meshblock/nx1=%d" % nx,
                        "mesh/nx2=%d" % nx, "meshblock/nx2=%d" % nx,
                        "mesh/x2min=0.0", "mesh/x2max=1.0", "problem/pulse_y0=0.5"]
                if nd == 3:
                    args += ["mesh/nx3=%d" % nx, "meshblock/nx3=%d" % nx,
                             "mesh/x3min=0.0", "mesh/x3max=1.0",
                             "problem/pulse_z0=0.5"]
                pair(log, name, args, tlim, tlim / 12.0, kappa,
                     inp=(INP3 if nd == 3 else None))


def g2(log):
    log("\n=== G2  anisotropic cells, dx2 = 4 dx1 and dx2 = dx1/4 ===")
    nx, tau = 64, 1.0e3
    kappa = tau * nx
    tlim = 3.1e-4 * kappa
    for tag, ny in (("dx2_4dx1", 16), ("dx2_dx1_4", 256)):
        args = ["rad_m1/kappa_s=%r" % kappa, "rad_m1/implicit_cfl=1",
                "mesh/nx1=%d" % nx, "meshblock/nx1=%d" % nx,
                "mesh/nx2=%d" % ny, "meshblock/nx2=%d" % ny,
                "mesh/x2min=0.0", "mesh/x2max=1.0", "problem/pulse_y0=0.5"]
        pair(log, "g2_%s" % tag, args, tlim, tlim / 12.0, kappa, tol=0.05)


def tabs(d):
    import numpy as np
    fs = sorted(f for f in os.listdir(os.path.join(d, "tab")) if f.endswith(".tab"))
    return np.loadtxt(os.path.join(d, "tab", fs[-1]))


def maxdiff(a, b):
    import numpy as np
    sc = np.maximum(np.abs(a).max(axis=0), np.abs(b).max(axis=0))
    sc[sc == 0.0] = 1.0
    return float(np.max(np.abs(a - b).max(axis=0) / sc))


def etot(d):
    import numpy as np
    files = sorted(os.path.join(d, "bin", f)
                   for f in os.listdir(os.path.join(d, "bin")))
    du = common.load_series([files[0], files[-1]])
    return [float(np.sum(x.var("m1_e"))) for x in du]


def g3(log):
    log("\n=== G3  bicgstab == line_jacobi on the SAME 2-D problem ===")
    nx, tau = 64, 1.0e3
    kappa = tau * nx
    tlim = 3.1e-4 * kappa
    for cfl in (1.0, 1.0e4):
        base = ["rad_m1/kappa_s=%r" % kappa, "rad_m1/implicit_cfl=%r" % cfl,
                "mesh/nx1=%d" % nx, "meshblock/nx1=%d" % nx,
                "mesh/nx2=%d" % nx, "meshblock/nx2=%d" % nx,
                "mesh/x2min=0.0", "mesh/x2max=1.0", "problem/pulse_y0=0.5",
                "output1/dt=1.0e30", "output3/dt=%r" % tlim]
        ds = []
        for s in SOLVERS:
            d, rc, _ = run("g3_cfl%g_%s" % (cfl, s),
                           base + ["rad_m1/implicit_solver=%s" % s], tlim, 1.0e30)
            ds.append(d)
        log("g3 cfl=%g  max relative difference of the final (E,F1,F2) slice = %.3e" %
            (cfl, maxdiff(tabs(ds[0])[:, 3:6], tabs(ds[1])[:, 3:6])))


def g4(log):
    log("\n=== G4  MeshBlock decomposition in 2-D, bicgstab ===")
    nx, tau = 64, 1.0e3
    kappa = tau * nx
    tlim = 3.1e-4 * kappa
    base = ["rad_m1/kappa_s=%r" % kappa, "rad_m1/implicit_cfl=1e4",
            "rad_m1/implicit_solver=bicgstab",
            "mesh/nx1=%d" % nx, "mesh/nx2=%d" % nx,
            "mesh/x2min=0.0", "mesh/x2max=1.0", "problem/pulse_y0=0.5",
            "output3/dt=%r" % tlim, "output1/dt=%r" % tlim]
    refl = ["mesh/ix1_bc=reflect", "mesh/ox1_bc=reflect",
            "rad_m1/implicit_partition=gather"]
    cases = (("a", [], ["meshblock/nx1=%d" % nx, "meshblock/nx2=%d" % nx],
              ["meshblock/nx1=%d" % nx, "meshblock/nx2=%d" % (nx // 2)]),
             ("b", refl, ["meshblock/nx1=%d" % nx, "meshblock/nx2=%d" % nx],
              ["meshblock/nx1=%d" % (nx // 2), "meshblock/nx2=%d" % (nx // 2)]))
    for tag, com, one, many in cases:
        d1, r1, _ = run("g4%s_b1" % tag, base + com + one, tlim, tlim)
        d2, r2, _ = run("g4%s_bn" % tag, base + com + many, tlim, tlim)
        if r1 or r2:
            log("g4%s  RUN FAILED rc=%d/%d" % (tag, r1, r2))
            continue
        e1, e2 = etot(d1), etot(d2)
        log("g4%s  maxdiff(single vs split) = %.3e   dEtot/Etot: single %.3e, split "
            "%.3e   outer %s / %s" %
            (tag, maxdiff(tabs(d1)[:, 3:7], tabs(d2)[:, 3:7]),
             (e1[1] - e1[0]) / e1[0], (e2[1] - e2[0]) / e2[0],
             stats(d1).get("picard"), stats(d2).get("picard")))
        for d in (d1, d2):
            shutil.rmtree(os.path.join(d, "bin"), ignore_errors=True)


def g5(log):
    log("\n=== G5  restart on a SINGLE MeshBlock, bicgstab ===")
    nx, tau = 64, 1.0e3
    kappa = tau * nx
    tlim = 3.1e-4 * kappa
    base = ["rad_m1/kappa_s=%r" % kappa, "rad_m1/implicit_cfl=1e4",
            "rad_m1/implicit_solver=bicgstab",
            "mesh/nx1=%d" % nx, "mesh/nx2=%d" % nx,
            "meshblock/nx1=%d" % nx, "meshblock/nx2=%d" % nx,
            "mesh/x2min=0.0", "mesh/x2max=1.0", "problem/pulse_y0=0.5",
            "output1/dt=1.0e30", "output3/dt=%r" % tlim,
            "output4/dt=%r" % (0.5 * tlim)]
    d1, _, _ = run("g5_full", base, tlim, 1.0e30)
    rst = sorted(f for f in os.listdir(os.path.join(d1, "rst")) if f.endswith(".rst"))
    mid = os.path.join(d1, "rst", rst[len(rst) // 2])
    d2 = os.path.join(OUT, "g5_restart")
    shutil.rmtree(d2, ignore_errors=True)
    os.makedirs(d2)
    with open(os.path.join(d2, "run.log"), "w") as f:
        rc = subprocess.call([XPL, "-r", mid], cwd=d2, stdout=f,
                             stderr=subprocess.STDOUT)
    f1 = sorted(os.listdir(os.path.join(d1, "tab")))[-1]
    f2 = sorted(os.listdir(os.path.join(d2, "tab")))[-1]
    same = subprocess.call(["cmp", "-s", os.path.join(d1, "tab", f1),
                            os.path.join(d2, "tab", f2)]) == 0
    log("g5  restart from %s (rc=%d): final slice %s" %
        (os.path.basename(mid), rc, "BYTE-IDENTICAL" if same else "DIFFERS"))
    if not same:
        log("g5  maxdiff = %.3e" % maxdiff(tabs(d1)[:, 3:7], tabs(d2)[:, 3:7]))


def g6(log):
    log("\n=== G6  regression: tests_gate_merge/postmerge.sh ===")
    p = subprocess.run(["bash", os.path.join(ROOT, "tests_gate_merge",
                                             "postmerge.sh")],
                       capture_output=True, text=True)
    log(p.stdout.strip())
    log("\n--- implicit_x1 thick pulse, tau_cell 1e3, implicit_cfl 100, central ---")
    kappa = 1.28e5
    d, _, _ = run("g6_i1", ["rad_m1/kappa_s=%r" % kappa, "rad_m1/implicit_cfl=100",
                            "mesh/nx1=128", "meshblock/nx1=128"], 12.0, 1.2,
                  extra=["rad_m1/transport=implicit_x1"])
    log(analyse(d, kappa, "g6_i1_implicit_x1"))
    shutil.rmtree(os.path.join(d, "bin"), ignore_errors=True)
    log("\n--- line_jacobi 2-D pulse, tau_cell 1e3, implicit_cfl 1e4 (3b3 G1 line) ---")
    nx, kappa = 64, 6.4e4
    tlim = 3.1e-4 * kappa
    d, _, _ = run("g6_lj", ["rad_m1/kappa_s=%r" % kappa, "rad_m1/implicit_cfl=1e4",
                            "mesh/nx1=%d" % nx, "meshblock/nx1=%d" % nx,
                            "mesh/nx2=%d" % nx, "meshblock/nx2=%d" % nx,
                            "mesh/x2min=0.0", "mesh/x2max=1.0",
                            "problem/pulse_y0=0.5",
                            "rad_m1/implicit_solver=line_jacobi"], tlim, tlim / 12.0)
    log("%s  outer=%s/%s" % (analyse(d, kappa, "g6_lj"), stats(d).get("picard"),
                             stats(d).get("pmax")))
    shutil.rmtree(os.path.join(d, "bin"), ignore_errors=True)


def main():
    which = sys.argv[1:] or ["g1", "g2", "g3", "g4", "g5", "g6"]
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)
    for g in which:
        globals()[g](log)
    with open(os.path.join(HERE, "gates_%s.txt" % "_".join(which)), "w") as f:
        f.write("\n".join(lines) + "\n")
    if TABLE:
        pj = os.path.join(ROOT, "tests_m1", "plots", "impl3b4_iters.json")
        old = {}
        if os.path.exists(pj):
            old = json.load(open(pj))
        old.update({r["case"]: {k: v for k, v in r.items() if k != "case"}
                    for r in TABLE})
        with open(pj, "w") as f:
            json.dump(old, f, indent=1)


if __name__ == "__main__":
    main()
