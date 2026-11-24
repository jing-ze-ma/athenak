"""
The two CHEAPENING switches on the super-time-stepped radiative conduction operator,
<hydro>/rad_sts_split and <hydro>/rad_sts_once (src/diffusion/conduction.cpp,
src/diffusion/conduction_transverse.cpp), on the same 3-D anisotropic spreading Gaussian
test_rad_sts_all_cpu.py runs the unified operator on.

WHAT THE TWO SWITCHES DO.

rad_sts_split is the STIFFNESS SPLIT.  Every face conductance is written as
C_f = C_exp,f + C_sts,f with C_exp,f = min(C_f, C_max,f), C_max,f the largest conductance
the current step can carry EXPLICITLY on that face:

    C_max,f = rad_sts_split_x min(V_i/alpha_i, V_j/alpha_j)/(nf beta_dt),

i.e. the explicit row-sum budget x_i <= rad_sts_split_x (0.5) of NewTimeStep's own
stability limit, shared evenly over the nf = 2 ndim faces of the row and held to the
smaller of the two cells' budgets.  (With rad_sts_once the budget is formed from the
FULL dt instead of the stage's beta_dt, so that every stage removes the same C_exp and
the two parts still add up to the whole operator.)
C_exp is added to the ordinary face fluxes inside the
RK stages -- so the non-stiff part of the operator is coupled to the hydro at the
integrator's order instead of being operator-split from it at first order -- and only
C_sts goes into the RKL1 loop.

rad_sts_once applies that RKL1 loop ONCE per cycle, after the last stage and over the
full dt, instead of once per stage over beta_dt.  An RKL1 super-step of tau costs
sqrt(tau) substages, so one call over dt is cheaper than two over dt/2 (MEASURED: 1920 ->
1152 substages over this test at 64^2, -40 %).  The splitting against the hydro is first
order either way; this only doubles the size of a term that is already there, which is
why the two switches belong in one test: `once` raises the first-order error and `split`
lowers it, and together they still converge.

WHAT THIS GUARDS.
  * the explicit part is in FLUX form and the RKL1 part still is, so the operator stays
    conservative to round-off -- the same gate the other two rad tests apply, on a path
    where the conductance is now shared between two operators,
  * the split cannot open a face the RKL1 loop treats as closed, and cannot double-count
    one: either would show here as a first-order error, since the x1 walls are closed and
    x2/x3 are periodic,
  * rad_sts_once really does run once per CYCLE and not once per stage (the operator's
    own call counter is checked against the cycle count),
  * and the accuracy of the combination is held to what was measured, with 15 % of margin.

test_rad_sts_all_cpu.py covers the same input with both switches OFF, so the two tests
together also pin that the switches are inert when they are not set.
"""

# Modules
import os
import shutil
import subprocess

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
INPUT = os.path.join(REPO, "inputs", "tests", "rad_sts_all_gauss.athinput")
BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_sts_split")

# cfl_number at nx = 32, halved with every doubling of the resolution so that dt scales
# as dx^2 -- the same scaling test_rad_sts_all_cpu.py uses, and what makes a
# first-order-in-time super-step converge at the 2nd order of the spatial stencil.
CFL32 = 0.00465

# L1 error in the total energy against the analytic solution, with split+once.
# MEASURED on CPU: 5.959181e-07 (32), 1.347863e-07 (64), 3.365457e-08 (128), orders
# 2.14 and 2.00.  (Both switches off: 1.299634e-07, 3.281713e-08, 8.545311e-09.  `once`
# alone costs a factor ~4 in L1 at these resolutions -- it is a bigger first-order
# splitting term -- and `split` alone BUYS 14 %, 1.119904e-07 against 1.299634e-07 at 32:
# the non-stiff part of the operator is no longer split off from the hydro at all.)
# The gates carry 15 % of margin.
L1_MAX = {32: 6.85e-7, 64: 1.55e-7}
ORDER_MIN = 1.85
# |sum V de| / sum V |de| of one call of the RKL1 part.  Flux form on both sides of the
# split, every physical face closed: round-off and nothing else (measured 1e-15).
CONS_MAX = 1.0e-12
# RKL1 substages of one call.  With rad_sts_once the call covers the FULL dt, so it takes
# MORE substages than a per-stage call (measured 10 at 32, 9 at 64) while taking fewer of
# them in total.
SUBSTAGE_MAX = 12


def run(binary, n):
    """run at nx1 = nx2 = n with split+once; returns (L1 energy error, stdout)."""
    rundir = os.path.join(BUILD, "run_%d" % n)
    os.makedirs(rundir, exist_ok=True)
    args = ["mesh/nx1=%d" % n, "mesh/nx2=%d" % n,
            "meshblock/nx1=%d" % n, "meshblock/nx2=%d" % n,
            "time/cfl_number=%.10g" % (CFL32 * 32.0 / n),
            "hydro/rad_sts_split=true", "hydro/rad_sts_once=true"]
    proc = subprocess.run([binary, "-i", INPUT] + args, cwd=rundir,
                          capture_output=True, text=True)
    if proc.returncode != 0:
        pytest.fail("run failed:\n" + proc.stdout[-4000:])
    errs = os.path.join(rundir, "sts3d-errs.dat")
    assert os.path.exists(errs), "no error file was written"
    row = [ln for ln in open(errs).read().splitlines() if not ln.startswith("#")][-1]
    return float(row.split()[-1]), proc.stdout


def parse_report(stdout):
    """(max substages of a call, max conservation residual, calls, cycles)."""
    subs, viol = [], []
    for line in stdout.splitlines():
        if "### rad_sts_all" not in line:
            continue
        assert "CLAMPED" not in line, "the substage count was clamped: " + line
        subs.append(int(line.split("substages = ")[1].split(",")[0]))
        viol.append(float(line.rsplit("=", 1)[1]))
    assert subs, "the operator never reported: is the switch on?"
    calls = [ln for ln in stdout.splitlines() if "### rad_sts totals:" in ln]
    assert calls, "the operator never reported its totals"
    ncall = int(calls[-1].split("totals: ")[1].split(" calls")[0])
    ncycle = int([ln for ln in stdout.splitlines()
                  if ln.startswith("time=")][-1].split("cycle=")[1])
    return max(subs), max(viol), ncall, ncycle


def test_run():
    try:
        subprocess.run(["cmake", "-S", REPO, "-B", BUILD,
                        "-D", "CMAKE_BUILD_TYPE=Release"],
                       check=True, capture_output=True, text=True)
        subprocess.run(["make", "-C", BUILD, "-j", str(os.cpu_count()), "athena"],
                       check=True, capture_output=True, text=True)
        binary = os.path.join(BUILD, "src", "athena")
        l1 = {}
        for n in (32, 64):
            l1[n], out = run(binary, n)
            nsub, viol, ncall, ncycle = parse_report(out)
            assert nsub <= SUBSTAGE_MAX, \
                f"nx={n}: {nsub} RKL1 substages, expected <= {SUBSTAGE_MAX}"
            assert viol < CONS_MAX, \
                f"nx={n}: the split operator is not conservative, " \
                f"|sum V de|/sum V|de| = {viol:.3g}"
            # rad_sts_once: ONE call per cycle, not one per RK stage
            assert ncall == ncycle, \
                f"nx={n}: the operator ran {ncall} times in {ncycle} cycles; " \
                f"rad_sts_once should make that one call per cycle"
            assert l1[n] < L1_MAX[n], \
                f"nx={n}: split+once L1 energy error {l1[n]:.4g}, " \
                f"gate {L1_MAX[n]:.4g}"
        order = np.log2(l1[32] / l1[64])
        assert order > ORDER_MIN, \
            f"convergence order {order:.3f} " \
            f"(L1 {l1[32]:.4g} -> {l1[64]:.4g}), gate {ORDER_MIN}"
    finally:
        shutil.rmtree(BUILD, ignore_errors=True)
