"""
The IMPLICIT TRANSVERSE radiative diffusion operator (<hydro>/rad_implicit_ang,
src/diffusion/conduction_transverse.cpp) against the analytic spreading Gaussian.

A Cartesian box at rest, uniform in x1, with a 1e-4 Gaussian temperature perturbation in
(x2,x3) on a uniform density (pgen rad_diff2d, with
inputs/tests/rad_transverse_gauss.athinput).
At that amplitude the radiative conductivity is constant to a few parts in 1e4, so the
perturbation obeys dT/dt = D grad^2 T and the exact solution is

    dT(r,t) = amp T0 (s0^2/s^2) exp(-r^2/(2 s^2)),   s^2 = s0^2 + 2 D t

(with the periodic images summed).  D is formed in the problem generator from the same
Freedman Rosseland mean and the same unit conversion the conduction module uses, so the
reference carries no separately tuned constant, and the conductivity is made large enough
that the blob diffuses ~100x faster than sound crosses it: the gas cannot respond and
conduction is the only operator that moves anything.

WHAT THIS GUARDS.  dx1 = 10 dx2, and x1 -- where the state is uniform and no flux ever
flows -- is the direction that sets the timestep, so the transverse solver is asked to
cover exactly 100 explicit transverse steps per step, in 7 RKL1 substages.  Halving dx2
and dx1 together holds that ratio and scales dt as dx^2, so the scheme converges at 2nd
order.  The gates are the L1 error in the energy at both resolutions, the observed order,
and the conservation residual sum V de / sum V |de| of every call, which must be at
round-off.

THE ERROR FLOOR, and why the gates are where they are.  Measured over four resolutions
n = 32, 64, 128, 256 (L1 = 1.0373e-7, 2.7855e-8, 8.8556e-9, 4.1603e-9), the error is not
a pure power law; it fits

    L1(n) = 6.32e-9 (128/n)^2 + 2.56e-9

to within a percent.  The operator is therefore EXACTLY 2nd order, sitting on a
resolution-INDEPENDENT floor of 2.56e-9.  That floor belongs to the REFERENCE, not to the
scheme: the analytic solution assumes a constant D, while the real conductivity varies by
~1e-4 across the blob, and the gas back-reacts a little.  Refining cannot remove it, so at
high n the apparent order sags -- 64 -> 128 gives 1.65 and 128 -> 256 only 1.09 -- and any
gate calibrated up there is measuring the reference's error, not the operator's.  Hence
this test runs the two COARSEST resolutions, 32 and 64, where the floor is 2.5 % and 9 %
of the error and the measured order is 1.897.  Do not "improve" it by adding finer runs.

This test is deliberately BLIND to one thing, and the blindness is the point of the
design note: the state is uniform in x1, so the radial operator moves nothing and the
lag between the two split operators -- which is what destabilised the He-star FeCZ box --
shows up here as max |dT*|/T* = 3e-6 instead of the 5-7 % it reaches in that box.  Do not
treat a pass here as evidence that the operator splitting is consistent.
"""

# Modules
import os
import shutil
import subprocess

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
INPUT = os.path.join(REPO, "inputs", "tests", "rad_transverse_gauss.athinput")
BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_tr_gauss")

# L1 error in the total energy against the analytic solution, and the order.
# MEASURED 1.0373e-7 at 32^2 and 2.7855e-8 at 64^2, order 1.897, identical on CPU and on
# one and two MI300A.  The gates carry ~1.15x of margin.  These numbers moved once, when
# the operator stopped linearising about the stage-START state and started re-evaluating
# T* and alpha from the conserved energy it actually writes on (conduction_transverse.cpp)
# -- that refresh removed a splitting error and cut the 64^2 error 3.1x -- so a run that
# comes back ~3x high is the lag returning, and will trip the L1 gates.
L1_MAX = {32: 1.2e-7, 64: 3.3e-8}
ORDER_MIN = 1.85
# |sum V de| / sum V |de| of one call of the operator
CONS_MAX = 1.0e-12
# the super-step the solver is asked to cover, in explicit transverse steps
SUBSTAGE_MAX = 12


def run(binary, n):
    """run at nx2 = nx3 = n; returns (L1 energy error, stdout)."""
    rundir = os.path.join(BUILD, "run_%d" % n)
    os.makedirs(rundir, exist_ok=True)
    # x1max scales with the transverse cell size, so dx1/dx2 -- and with it the ratio of
    # the step to the explicit transverse limit -- is the same at both resolutions
    args = ["mesh/nx2=%d" % n, "mesh/nx3=%d" % n,
            "meshblock/nx2=%d" % n, "meshblock/nx3=%d" % n,
            "mesh/x1max=%.10g" % (0.625 * 64.0 / n)]
    proc = subprocess.run([binary, "-i", INPUT] + args, cwd=rundir,
                          capture_output=True, text=True)
    if proc.returncode != 0:
        pytest.fail("run failed:\n" + proc.stdout[-4000:])
    errs = os.path.join(rundir, "lap2d-errs.dat")
    assert os.path.exists(errs), "no error file was written"
    row = [ln for ln in open(errs).read().splitlines() if not ln.startswith("#")][-1]
    return float(row.split()[-1]), proc.stdout


def parse_report(stdout):
    """(max substage count, max conservation residual) over the reported calls."""
    subs, viol = [], []
    for line in stdout.splitlines():
        if "### rad_implicit_ang" not in line:
            continue
        assert "CLAMPED" not in line, "the substage count was clamped: " + line
        subs.append(int(line.split("substages = ")[1].split(",")[0]))
        viol.append(float(line.rsplit("=", 1)[1]))
    assert subs, "the operator never reported: is rad_implicit_ang on?"
    return max(subs), max(viol)


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
            nsub, viol = parse_report(out)
            assert nsub <= SUBSTAGE_MAX, \
                f"nx={n}: {nsub} RKL1 substages for a 100x step, expected <= " \
                f"{SUBSTAGE_MAX}"
            assert viol < CONS_MAX, \
                f"nx={n}: the operator is not conservative, |sum V de|/sum V|de| = " \
                f"{viol:.3g}"
            assert l1[n] < L1_MAX[n], \
                f"nx={n}: L1 energy error {l1[n]:.4g}, gate {L1_MAX[n]:.4g}"
        order = np.log2(l1[32] / l1[64])
        assert order > ORDER_MIN, \
            f"convergence order {order:.3f} (L1 {l1[32]:.4g} -> {l1[64]:.4g}), " \
            f"gate {ORDER_MIN}"
    finally:
        shutil.rmtree(BUILD, ignore_errors=True)
