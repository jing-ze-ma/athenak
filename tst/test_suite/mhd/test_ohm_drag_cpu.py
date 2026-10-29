"""
Ohmic (Lorentz) drag with a steep eta(r) on the cubed sphere (src/pgen/ohm_drag.cpp,
inputs/tests/ohm_drag_cs.athinput), against the exact 1-D reference of
inputs/tests/ohm_drag_ref1d.py.

Shell r in [1,2], axial dipole (m = 1) threading two conducting walls at rest, shellular
rotation Omega = 0.25 sin^4(pi (r-1)) HELD by a user source (kinematic), and
eta(r) = 0.1 * 10^(-4 (1 + tanh((r - 1.6)/0.08))/2): four decades across ~0.3 in r, low
Rm below the transition, ~1e4 above.  The induced field is exactly
B_phi = b(r,t) sin cos, u = r b obeying u_t = (eta u')' - 6 eta u/r^2 + 2 m Omega'/r with
u' = 0 on both walls; the reference solves it exactly in time on a fine FV grid.

Gated at t = 0.5, nx1 = nx2 = nx3 = N:
  * B_phi L1 (volume-weighted, all tangential faces, relative to |ref|) at N = 16, 32;
  * the total Lorentz drag torque (the L_z the kinematic hold removes per stage / beta dt)
    over 16 pi m (b(2) - b(1))/15;
  * the L1 order of B_phi between N = 16 and 32;
  * RKG super-time-stepping (mhd/use_rkg_sts) against the explicit path at N = 16.
The thresholds below carry the values measured on stswall-1029 c306af525 +
ohmdrag-1029 (viper login node, MPI build) in the comments, with a modest margin.

The binary is built here with -D PROBLEM=ohm_drag (and MPI when mpirun exists, as
rad_m1/m1_common.py does); the runs use up to 12 ranks.
"""

import os
import re
import shutil
import subprocess
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
INPUT = os.path.join(REPO, "inputs", "tests", "ohm_drag_cs.athinput")
REF = os.path.join(REPO, "inputs", "tests", "ohm_drag_ref1d.py")
BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_ohm_drag")
HAVE_MPI = shutil.which("mpirun") is not None
TLIM = 0.5

# Thresholds (measured value in the comment)
BPHI_L1_MAX = {16: 0.065,     # measured 0.05673
               32: 0.025}     # measured 0.02159
TORQUE_TOL = {16: 0.045,      # |ratio - 1|, measured 0.0332 (ratio 0.96683)
              32: 0.018}      # measured 0.0125 (ratio 0.98755)
ORDER_MIN = 1.25              # B_phi L1 order 16 -> 32, measured 1.394
RKG_REL_MAX = 0.10            # |Bphi_L1(rkg)/Bphi_L1(ex) - 1| at N = 16, measured 0.055
#                               (RKG temporal error at t = 0.5; 0.005 at t = 2)
RKG_TORQUE_MAX = 0.006        # |torque ratio rkg - explicit| at N = 16, measured 0.0031


def ncpu():
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count()


def ranks():
    """24 MeshBlocks per run: the largest divisor of 24 the cores allow (max 12)."""
    if not HAVE_MPI:
        return 0
    n = ncpu()
    for r in (12, 8, 6, 4, 3, 2):
        if r <= n:
            return r
    return 0


def build():
    exe = os.path.join(BUILD, "src", "athena")
    if os.path.exists(exe):
        return exe
    flags = ["-D", "Athena_ENABLE_MPI=ON"] if HAVE_MPI else []
    subprocess.run(["cmake", "-S", REPO, "-B", BUILD, "-D", "PROBLEM=ohm_drag",
                    "-D", "CMAKE_BUILD_TYPE=Release"] + flags,
                   check=True, capture_output=True, text=True)
    subprocess.run(["make", "-C", BUILD, "-j", str(min(16, ncpu()))],
                   check=True, capture_output=True, text=True)
    return exe


def reference(rundir):
    table = os.path.join(rundir, "ref.dat")
    subprocess.run([sys.executable, REF, "--eta_bot", "0.1", "--ndec", "4", "--rt", "1.6",
                    "--w", "0.08", "--om0", "0.25", "--m", "1", "--T", repr(TLIM),
                    "--N", "1000", "--out", table],
                   check=True, capture_output=True, text=True)
    return table


def run(exe, rundir, table, n, rkg):
    tag = "n%d_%s" % (n, "rkg" if rkg else "ex")
    d = os.path.join(rundir, tag)
    os.makedirs(d, exist_ok=True)
    errs = os.path.join(d, "errs.dat")
    nr = ranks()
    launch = ["mpirun", "-np", repr(nr), "--oversubscribe"] if nr else []
    cmd = launch + [exe, "-i", INPUT, "-d", d,
                    "mesh/nx1=%d" % n, "mesh/nx2=%d" % n, "mesh/nx3=%d" % n,
                    "meshblock/nx1=%d" % n, "meshblock/nx2=%d" % (n // 2),
                    "meshblock/nx3=%d" % (n // 2), "time/tlim=%r" % TLIM,
                    "mhd/use_rkg_sts=%s" % ("true" if rkg else "false"),
                    "problem/ref_file=%s" % table, "problem/err_file=%s" % errs]
    env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run(cmd, cwd=d, capture_output=True, text=True, env=env)
    if proc.returncode != 0 or not os.path.exists(errs):
        pytest.fail(f"{tag} failed:\n{proc.stdout[-3000:]}\n{proc.stderr[-2000:]}")
    line = open(errs).read().strip().splitlines()[-1]
    bphi = float(re.search(r"Bphi_L1 (\S+)", line).group(1))
    torque = float(re.findall(r"ratio (\S+)", line)[1])
    print(f"ohm_drag {tag}: Bphi_L1 {bphi:.4e} torque ratio {torque:.5f}")
    return bphi, torque


def test_ohm_drag_cpu():
    rundir = os.path.join(BUILD, "runs")
    try:
        exe = build()
        os.makedirs(rundir, exist_ok=True)
        table = reference(rundir)
        res = {n: run(exe, rundir, table, n, False) for n in (16, 32)}
        rkg = run(exe, rundir, table, 16, True)
        for n, (bphi, torque) in res.items():
            assert bphi < BPHI_L1_MAX[n], \
                f"N={n}: B_phi L1 {bphi:.4e} > {BPHI_L1_MAX[n]}"
            assert abs(torque - 1.0) < TORQUE_TOL[n], \
                f"N={n}: torque ratio {torque:.5f}, tolerance {TORQUE_TOL[n]}"
        order = np.log2(res[16][0] / res[32][0])
        assert order > ORDER_MIN, f"B_phi L1 order 16->32 {order:.3f} < {ORDER_MIN}"
        assert abs(rkg[0] / res[16][0] - 1.0) < RKG_REL_MAX, \
            f"RKG B_phi L1 {rkg[0]:.4e} vs explicit {res[16][0]:.4e}"
        assert abs(rkg[1] - res[16][1]) < RKG_TORQUE_MAX, \
            f"RKG torque ratio {rkg[1]:.5f} vs explicit {res[16][1]:.5f}"
    finally:
        shutil.rmtree(BUILD, ignore_errors=True)
