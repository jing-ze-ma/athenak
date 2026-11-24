"""
The UNIFIED super-time-stepped radiative conduction operator (<hydro>/rad_sts_all,
Conduction::StsConductionUpdate in src/diffusion/conduction_transverse.cpp): ONE RKL1
loop over x1, x2 AND x3, against the analytic solution of a 3-D ANISOTROPIC spreading
Gaussian, and against the SPLIT operator (rad_implicit_x1 + rad_implicit_ang) it
replaces.

THE PROBLEM (pgen rad_diff2d with problem/sig1 > 0, inputs/tests/rad_sts_all_gauss.
athinput).  A Cartesian box at rest with a 1e-4 Gaussian temperature perturbation spread
in BOTH x1 and x2 with DIFFERENT widths (sig1 = 1.5 sig0) and uniform in x3, on a uniform
density.  At that amplitude the radiative conductivity is constant to a few parts in 1e4,
so the perturbation obeys dT/dt = D grad^2 T with a constant D and the exact solution is
the product of two spreading 1-D Gaussians.  x2/x3 are periodic; x1 is closed by
zero-flux (reflecting) walls, and a blob centred on the midpoint of a symmetric interval
has its Neumann images at exactly the periodic positions, so summing images makes the
reference exact in x1 as well.  D is formed in the problem generator from the same
Freedman Rosseland mean and the same unit conversion the conduction module uses.

WHAT THIS GUARDS, and what the transverse test cannot.  test_rad_transverse_gauss_cpu.py
runs a state that is UNIFORM IN x1: the radial operator moves nothing there, so that test
is blind to everything about the x1 direction and to the splitting error between the two
operators.  Here both x1 and x2 carry flux, and a different amount of it, so:
  * the x1 faces of the 7-point stencil are exercised (a wrong face conductance, a
    missing 1/dx1^2 or a closed interior face all show up as a first-order error),
  * the two arms -- unified and split -- are run at exactly the same timestep from
    exactly the same initial data, so their MUTUAL difference is the splitting error of
    the split arm and nothing else, and it has to converge away at 2nd order too.

THE TIMESTEP.  rad_sts_all drops the conduction constraint in all three directions, so dt
is the hydro CFL, which is proportional to dx and would make a first-order-in-time RKL1
super-step converge at 1st order overall.  The harness therefore scales cfl_number with
1/n, which restores dt proportional to dx^2 -- the same thing the transverse test does by
scaling x1max -- and sets the step at about 100x the explicit diffusion limit at n = 32.
Both arms use the same cfl, so both take the same steps.
"""

# Modules
import os
import shutil
import subprocess

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
INPUT = os.path.join(REPO, "inputs", "tests", "rad_sts_all_gauss.athinput")
BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_sts_all")

# cfl_number at nx = 32; halved with every doubling of the resolution so that dt scales
# as dx^2.  At this value the super-step is ~100 explicit diffusion steps.
CFL32 = 0.00465

# L1 error in the total energy against the analytic solution.  MEASURED on CPU:
#   unified  1.299634e-07 (32)  3.281713e-08 (64)  8.545311e-09 (128), orders 1.99, 1.94
#   split    5.594047e-07 (32)  1.493199e-07 (64)  3.781496e-08 (128), orders 1.91, 1.98
# Both arms are clean 2nd order, and the unified operator is 4.3-4.6x MORE ACCURATE at
# every resolution: that factor is the splitting error the split arm carries and the
# unified one does not.  The gates below carry 1.15x of margin.
L1_MAX = {32: 1.50e-7, 64: 3.78e-8}
L1_MAX_SPLIT = {32: 6.44e-7, 64: 1.72e-7}
ORDER_MIN = 1.85
# ... and the unified operator must stay clearly the more accurate of the two (measured
# 4.30x at 32 and 4.55x at 64): if this sags to 1 the x1 faces have stopped being part of
# one linearisation and the unified path has silently become the split one.
SPLIT_GAIN_MIN = 3.0
# |sum V de| / sum V |de| of one call of the operator.  The stencil is in flux form and
# every physical face is closed, so this is round-off and nothing else (measured 4e-16).
CONS_MAX = 1.0e-12
# the super-step the solver is asked to cover, in RKL1 substages: measured 10 for the
# unified operator against 8 for the transverse-only one, i.e. 25 % more substages to
# cover all three directions instead of two (z_max 95.97 against 63.98, exactly 3/2)
SUBSTAGE_MAX = 12


def run(binary, n, unified):
    """run at nx1 = nx2 = n; returns (L1 energy error, stdout)."""
    tag = "sts" if unified else "split"
    rundir = os.path.join(BUILD, "run_%s_%d" % (tag, n))
    os.makedirs(rundir, exist_ok=True)
    args = ["mesh/nx1=%d" % n, "mesh/nx2=%d" % n,
            "meshblock/nx1=%d" % n, "meshblock/nx2=%d" % n,
            "time/cfl_number=%.10g" % (CFL32 * 32.0 / n)]
    if not unified:
        args += ["hydro/rad_sts_all=false", "hydro/rad_implicit_x1=true",
                 "hydro/rad_implicit_ang=true"]
    proc = subprocess.run([binary, "-i", INPUT] + args, cwd=rundir,
                          capture_output=True, text=True)
    if proc.returncode != 0:
        pytest.fail("run failed:\n" + proc.stdout[-4000:])
    errs = os.path.join(rundir, "sts3d-errs.dat")
    assert os.path.exists(errs), "no error file was written"
    row = [ln for ln in open(errs).read().splitlines() if not ln.startswith("#")][-1]
    return float(row.split()[-1]), proc.stdout


def parse_report(stdout, tag):
    """(max substage count, max conservation residual) over the reported calls."""
    subs, viol = [], []
    for line in stdout.splitlines():
        if tag not in line:
            continue
        assert "CLAMPED" not in line, "the substage count was clamped: " + line
        subs.append(int(line.split("substages = ")[1].split(",")[0]))
        viol.append(float(line.rsplit("=", 1)[1]))
    assert subs, "the operator never reported: is the switch on?"
    return max(subs), max(viol)


def test_run():
    try:
        subprocess.run(["cmake", "-S", REPO, "-B", BUILD,
                        "-D", "CMAKE_BUILD_TYPE=Release"],
                       check=True, capture_output=True, text=True)
        subprocess.run(["make", "-C", BUILD, "-j", str(os.cpu_count()), "athena"],
                       check=True, capture_output=True, text=True)
        binary = os.path.join(BUILD, "src", "athena")
        l1, l1s = {}, {}
        for n in (32, 64):
            l1[n], out = run(binary, n, True)
            nsub, viol = parse_report(out, "### rad_sts_all")
            assert nsub <= SUBSTAGE_MAX, \
                f"nx={n}: {nsub} RKL1 substages, expected <= {SUBSTAGE_MAX}"
            assert viol < CONS_MAX, \
                f"nx={n}: the unified operator is not conservative, " \
                f"|sum V de|/sum V|de| = {viol:.3g}"
            assert l1[n] < L1_MAX[n], \
                f"nx={n}: unified L1 energy error {l1[n]:.4g}, gate {L1_MAX[n]:.4g}"
            # the split arm it replaces, at the same timestep
            l1s[n], out = run(binary, n, False)
            _, viol = parse_report(out, "### rad_implicit_ang")
            assert viol < CONS_MAX, \
                f"nx={n}: the transverse operator is not conservative, " \
                f"|sum V de|/sum V|de| = {viol:.3g}"
            assert l1s[n] < L1_MAX_SPLIT[n], \
                f"nx={n}: split L1 energy error {l1s[n]:.4g}, " \
                f"gate {L1_MAX_SPLIT[n]:.4g}"
        order = np.log2(l1[32] / l1[64])
        assert order > ORDER_MIN, \
            f"unified convergence order {order:.3f} " \
            f"(L1 {l1[32]:.4g} -> {l1[64]:.4g}), gate {ORDER_MIN}"
        orders = np.log2(l1s[32] / l1s[64])
        assert orders > ORDER_MIN, \
            f"split convergence order {orders:.3f} " \
            f"(L1 {l1s[32]:.4g} -> {l1s[64]:.4g}), gate {ORDER_MIN}"
        for n in (32, 64):
            gain = l1s[n] / l1[n]
            assert gain > SPLIT_GAIN_MIN, \
                f"nx={n}: the unified operator is only {gain:.2f}x more accurate than " \
                f"the split one, gate {SPLIT_GAIN_MIN}"
    finally:
        shutil.rmtree(BUILD, ignore_errors=True)
