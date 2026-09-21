#!/usr/bin/env python3
"""Milestone 3b phase B gates: <rad_m1>/transport = implicit, the TRANSVERSE (x2/x3)
implicit solve with implicit_solver = line_jacobi.

Every case runs, is analysed and has its dumps deleted immediately (the viper inode
quota is at its limit).  Run from tests_m1/runs_3b3 with no arguments, or with the gate
names to run as arguments (g1 g2 g3 g4 g5 g6).
"""

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import common  # noqa: E402
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
XPL = os.path.join(ROOT, "build_cpu_m1", "src", "athena")
INP = os.path.join(HERE, "pulse_md.athinput")
OUT = os.environ.get("M1_3B3_OUT", os.path.join(HERE, "out"))


def run(name, args, tlim, odt, extra=()):
    d = os.path.join(OUT, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    cmd = ([XPL, "-i", INP, "time/cfl_number=1e6", "rad_m1/transport=implicit",
            "time/tlim=%r" % tlim, "output1/dt=%r" % odt] + list(args) + list(extra))
    with open(os.path.join(d, "run.log"), "w") as f:
        rc = subprocess.call(cmd, cwd=d, stdout=f, stderr=subprocess.STDOUT)
    return d, rc


def stats(d):
    txt = open(os.path.join(d, "run.log")).read()
    out = {}
    for key, tag in (("picard", "iterations mean="), ("pmax", "max="),
                     ("lres", "linear residual mean=")):
        i = txt.find(tag)
        if i >= 0:
            out[key] = txt[i + len(tag):].split()[0]
    i = txt.find("implicit transport: solves=")
    if i >= 0:
        seg = txt[i:i + 200]
        out["picard"] = seg.split("iterations mean=")[1].split()[0]
        out["pmax"] = seg.split("max=")[1].split()[0]
    return out


def analyse(d, kappa, label, tol=0.02):
    files = sorted(os.path.join(d, "bin", f)
                   for f in os.listdir(os.path.join(d, "bin")))
    p = subprocess.run([sys.executable, os.path.join(HERE, "pulse_md.py")] + files +
                       ["--kappa", repr(kappa), "--label", label, "--tol", repr(tol)],
                       capture_output=True, text=True)
    return p.stdout.strip() + p.stderr.strip()


def g1(log):
    log("\n=== G1  isotropy and rate, transport = implicit, line_jacobi ===")
    log("%-30s %-34s %-24s" % ("case", "d(sigma^2)/dt / 2D per axis", "Picard mean/max"))
    for nd, nx in ((2, 64), (3, 32)):
        for tau in (10.0, 1.0e3):
            kappa = tau * nx
            diff = 1.0 / (3.0 * kappa)
            for cfl in (1.0, 1.0e2, 1.0e4):
                dt = cfl / nx
                tsmall = 3.1e-4 * kappa
                tlim = max(tsmall, 12.0 * dt)
                if 2.0 * diff * tlim > 0.02:
                    log("%-30s NOT MEASURABLE (the diffusion length of 12 steps "
                        "exceeds the box)" % ("g1_%dd_tau%g_cfl%g" % (nd, tau, cfl)))
                    continue
                name = "g1_%dd_tau%g_cfl%g" % (nd, tau, cfl)
                args = ["rad_m1/kappa_s=%r" % kappa, "rad_m1/implicit_cfl=%r" % cfl,
                        "mesh/nx1=%d" % nx, "meshblock/nx1=%d" % nx,
                        "mesh/nx2=%d" % nx, "meshblock/nx2=%d" % nx,
                        "mesh/x2min=0.0", "mesh/x2max=1.0", "problem/pulse_y0=0.5"]
                if nd == 3:
                    args += ["mesh/nx3=%d" % nx, "meshblock/nx3=%d" % nx,
                             "mesh/x3min=0.0", "mesh/x3max=1.0",
                             "problem/pulse_z0=0.5"]
                d, rc = run(name, args, tlim, tlim / 12.0)
                s = stats(d)
                log("%s  picard=%s/%s lres=%s rc=%d" %
                    (analyse(d, kappa, name), s.get("picard"), s.get("pmax"),
                     s.get("lres"), rc))
                shutil.rmtree(os.path.join(d, "bin"), ignore_errors=True)


def g2(log):
    log("\n=== G2  anisotropic cells, dx2 = 4 dx1 and dx2 = dx1/4 ===")
    nx = 64
    tau = 1.0e3
    kappa = tau * nx
    tlim = 3.1e-4 * kappa
    # the SAME physical box, a different x2 cell count: dx2 = 4 dx1 (the He box ratio)
    # and dx2 = dx1/4.  The variance-growth rate of the 3-point diffusion operator is
    # exactly 2 D at any resolution, so a coarse transverse direction is a legitimate
    # accuracy test and not a discretisation artefact.
    for tag, ny in (("dx2_4dx1", 16), ("dx2_dx1_4", 256)):
        name = "g2_%s" % tag
        args = ["rad_m1/kappa_s=%r" % kappa, "rad_m1/implicit_cfl=1",
                "mesh/nx1=%d" % nx, "meshblock/nx1=%d" % nx,
                "mesh/nx2=%d" % ny, "meshblock/nx2=%d" % ny,
                "mesh/x2min=0.0", "mesh/x2max=1.0", "problem/pulse_y0=0.5"]
        d, rc = run(name, args, tlim, tlim / 12.0)
        s = stats(d)
        log("%s  picard=%s/%s lres=%s rc=%d" %
            (analyse(d, kappa, name, tol=0.05), s.get("picard"), s.get("pmax"),
             s.get("lres"), rc))
        shutil.rmtree(os.path.join(d, "bin"), ignore_errors=True)


def tabs(d):
    """the final double-precision x2 slice of a run, as a flat array"""
    import numpy as np
    fs = sorted(f for f in os.listdir(os.path.join(d, "tab")) if f.endswith(".tab"))
    return np.loadtxt(os.path.join(d, "tab", fs[-1]))


def maxdiff(a, b):
    """max |a-b| per column, divided by the COLUMN's own max: the far tail of F1 sits
    20 orders below its peak and a per-element relative norm is meaningless there."""
    import numpy as np
    sc = np.maximum(np.abs(a).max(axis=0), np.abs(b).max(axis=0))
    sc[sc == 0.0] = 1.0
    return float(np.max(np.abs(a - b).max(axis=0) / sc))


def etot(d):
    """the volume sum of E of the first and last bin dump"""
    import numpy as np
    files = sorted(os.path.join(d, "bin", f)
                   for f in os.listdir(os.path.join(d, "bin")))
    du = common.load_series([files[0], files[-1]])
    return [float(np.sum(x.var("m1_e"))) for x in du]


def g3(log):
    log("\n=== G3  a problem UNIFORM in x2 (nx2 = 4): implicit == implicit_x1 ===")
    nx, ny, tau = 64, 4, 1.0e3
    kappa = tau * nx
    tlim = 3.1e-4 * kappa
    base = ["rad_m1/kappa_s=%r" % kappa, "rad_m1/implicit_cfl=1",
            "mesh/nx1=%d" % nx, "meshblock/nx1=%d" % nx,
            "mesh/nx2=%d" % ny, "meshblock/nx2=%d" % ny,
            "mesh/x2min=0.0", "mesh/x2max=1.0", "problem/pulse_1d=true",
            "output1/dt=1.0e30", "output3/dt=%r" % tlim]
    d1, _ = run("g3_implicit", base, tlim, 1.0e30)
    d2, _ = run("g3_implicit_x1",
                base + ["rad_m1/implicit_allow_multid=true"], tlim, 1.0e30,
                extra=["rad_m1/transport=implicit_x1"])
    import numpy as np
    a, b = tabs(d1), tabs(d2)
    log("g3  max relative difference of the final (E,F1) slice = %.3e   "
        "max |F2| = %.3e (implicit), %.3e (implicit_x1)" %
        (maxdiff(a[:, 3:5], b[:, 3:5]), np.max(np.abs(a[:, 5])),
         np.max(np.abs(b[:, 5]))))
    log("g3  picard implicit=%s implicit_x1=%s" %
        (stats(d1).get("picard"), stats(d2).get("picard")))


def g4(log):
    log("\n=== G4  MeshBlock decomposition in 2-D ===")
    log("    (a) 1 x 2 blocks, periodic x1;  (b) 2 x 2 blocks, REFLECTING x1 with")
    log("    implicit_partition = gather -- periodic x1 across several blocks is the")
    log("    clean fatal of 3b sect. 8 and cannot be decomposed.")
    nx, tau = 64, 1.0e3
    kappa = tau * nx
    tlim = 3.1e-4 * kappa
    for lt in ("1.0e-10", "1.0e-13"):
        base = ["rad_m1/kappa_s=%r" % kappa, "rad_m1/implicit_cfl=1",
                "rad_m1/implicit_lin_tol=%s" % lt,
                "mesh/nx1=%d" % nx, "mesh/nx2=%d" % nx,
                "mesh/x2min=0.0", "mesh/x2max=1.0", "problem/pulse_y0=0.5",
                "output3/dt=%r" % tlim, "output1/dt=%r" % tlim]
        refl = ["mesh/ix1_bc=reflect", "mesh/ox1_bc=reflect",
                "rad_m1/implicit_partition=gather"]
        cases = (("a", [], ["meshblock/nx1=%d" % nx, "meshblock/nx2=%d" % nx],
                  ["meshblock/nx1=%d" % nx, "meshblock/nx2=%d" % (nx // 2)]),
                 ("b", refl, ["meshblock/nx1=%d" % nx, "meshblock/nx2=%d" % nx],
                  ["meshblock/nx1=%d" % (nx // 2),
                   "meshblock/nx2=%d" % (nx // 2)]))
        for tag, com, one, many in cases:
            d1, r1 = run("g4%s_b1_%s" % (tag, lt), base + com + one, tlim, tlim)
            d2, r2 = run("g4%s_bn_%s" % (tag, lt), base + com + many, tlim, tlim)
            if r1 or r2:
                log("g4%s lin_tol=%s  RUN FAILED rc=%d/%d" % (tag, lt, r1, r2))
                continue
            e1, e2 = etot(d1), etot(d2)
            log("g4%s lin_tol=%s  maxdiff = %.3e   dEtot/Etot: single %.3e, split "
                "%.3e   picard %s / %s" %
                (tag, lt, maxdiff(tabs(d1)[:, 3:7], tabs(d2)[:, 3:7]),
                 (e1[1] - e1[0]) / e1[0], (e2[1] - e2[0]) / e2[0],
                 stats(d1).get("picard"), stats(d2).get("picard")))
            for d in (d1, d2):
                shutil.rmtree(os.path.join(d, "bin"), ignore_errors=True)


def g5(log):
    log("\n=== G5  restart, 2-D, 2x2 MeshBlocks ===")
    nx, tau = 64, 1.0e3
    kappa = tau * nx
    tlim = 3.1e-4 * kappa
    base = ["rad_m1/kappa_s=%r" % kappa, "rad_m1/implicit_cfl=1",
            "mesh/nx1=%d" % nx, "mesh/nx2=%d" % nx,
            "meshblock/nx1=%d" % (nx // 2), "meshblock/nx2=%d" % (nx // 2),
            "mesh/ix1_bc=reflect", "mesh/ox1_bc=reflect",
            "rad_m1/implicit_partition=gather",
            "mesh/x2min=0.0", "mesh/x2max=1.0", "problem/pulse_y0=0.5",
            "output1/dt=1.0e30", "output3/dt=%r" % tlim,
            "output4/dt=%r" % (0.5 * tlim)]
    d1, _ = run("g5_full", base, tlim, 1.0e30)
    rst = sorted(f for f in os.listdir(os.path.join(d1, "rst"))
                 if f.endswith(".rst"))
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
    d, _ = run("g6_i1", ["rad_m1/kappa_s=%r" % kappa, "rad_m1/implicit_cfl=100",
                         "mesh/nx1=128", "meshblock/nx1=128"], 12.0, 1.2,
               extra=["rad_m1/transport=implicit_x1"])
    log(analyse(d, kappa, "g6_i1_implicit_x1"))
    shutil.rmtree(os.path.join(d, "bin"), ignore_errors=True)


def main():
    which = sys.argv[1:] or ["g1", "g2"]
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)
    for g in which:
        globals()[g](log)
    with open(os.path.join(HERE, "gates_%s.txt" % "_".join(which)), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
