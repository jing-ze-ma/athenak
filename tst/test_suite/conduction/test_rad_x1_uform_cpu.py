"""
<hydro>/rad_x1_uform and <hydro>/rad_x1_kiter, the two cures for the frozen-coefficient
error of the implicit radial radiative conduction (Conduction::ImplicitRadialUpdate,
src/diffusion/conduction.cpp), on the SPLIT operator (rad_implicit_x1 + rad_implicit_ang)
running the spreading-Gaussian problem of test_rad_sts_all_cpu.py.

WHAT THE TWO SWITCHES DO.

The solve linearises the column in T and freezes the face conductance
C_f = A_f K_f/dl_f, K ~ T^3/kappa_R(T,p), at the old state.  Where a cell carries a
through-flow that is a steady, one-signed error (MEASURED on the stellar columns this
operator was built for: 8 % of the through-flux on average and 38-50 % on the worst face
of a He-star column, 0.6 %/9-13 % at a B-star tau ramp).

rad_x1_uform solves the same column for du = d(T^4) instead of dT.  At frozen opacity
the flux is EXACTLY linear in u, F = -(a c/(3 kappa_R rho)) du/dl = -(K/(4 T^3)) du/dl,
so the T^3 half of the nonlinearity disappears outright.  The face recipe is untouched:
the off-diagonal is the SAME C_f divided by 4 T_f^3 at the frozen face temperature, and
the diagonal (rho c_v V)_i/(4 T_i^3 beta_dt) carries the same factor at the cell, so the
two forms are the same operator wherever T is uniform across the face.

rad_x1_kiter > 1 re-runs the whole assembly -- T, c_v, kappa_R(T,p), the limiter and the
face T^3 -- at the state the previous pass wrote and re-solves, a Picard iteration on
the backward-Euler balance, with the heat-capacity term held on T^n so that the passes
CORRECT one step instead of taking k of them.

WHAT THIS GUARDS.
  * both forms converge at the 2nd order of the stencil, and to the same solution: this
    problem is a weak (amp = 1e-4) perturbation on a uniform background, which is exactly
    where the two forms are supposed to coincide, and they do to 2e-6 in L1.  A face
    recipe that drifted between the forms -- a different T_f, a factor 4 out of place,
    an unscaled diagonal -- would show up here first and as a large effect;
  * the u-form rows are still a flux divergence, so pass 0 still telescopes and the
    column is conserved to round-off;
  * a further Picard pass leaves a converged solution alone, and cannot break the
    positivity clip or drive a row to the non-conservative fallback;
  * and with both switches at their defaults the solve is untouched, which
    test_rad_sts_all_cpu.py pins from the other side.

THE CONSERVATION NUMBER OF A LATER PASS.  Passes 2..k carry the T^n anchor, which is not
a flux divergence, so the routine compares the pass's column sum against the energy the
earlier passes put in rather than against zero.  On a converged pass both sides are
round-off and their RATIO is a 0/0; that is why the kiter arm is held to a loose gate and
the round-off gate is applied to the single-pass arms.
"""

# Modules
import math
import os
import shutil
import subprocess

import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
INPUT = os.path.join(REPO, "inputs", "tests", "rad_x1_uform_gauss.athinput")
BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_x1_uform")

# cfl_number at nx = 32, halved with every doubling so that dt scales as dx^2 -- the
# scaling test_rad_sts_all_cpu.py uses, and what makes the split operator converge at the
# 2nd order of its spatial stencil.
CFL32 = 0.00465

# L1 error of the total energy against the analytic solution.  MEASURED on CPU:
#   T-form (both switches off)  5.594047e-07 (32)  1.493199e-07 (64)  order 1.905
#   u-form, kiter 1             5.594037e-07       1.493192e-07      order 1.905
#   u-form, kiter 2             5.593995e-07       1.493184e-07      order 1.905
# The gates carry 15 % of margin.
L1_MAX = {32: 6.45e-7, 64: 1.72e-7}
ORDER_MIN = 1.85
# how far the u-form may sit from the T-form on THIS problem, where the perturbation is
# 1e-4 of the background and the two forms are the same operator to that order.
# Measured: 1.8e-6 (kiter 1) and 9.3e-6 (kiter 2) in relative L1.
AGREE_MAX = 1.0e-4
# |sum V de - expected|/sum V|de| of the solve.  Interior faces only, flux form: round-off
# (measured 6e-16 for the T-form and 9e-16 for the u-form).
CONS_MAX = 1.0e-12
# ...and for a second Picard pass, whose own increments are ~1e-8 of the first pass's, so
# that the same ratio is round-off over round-off (measured 1.0e-8).
CONS_MAX_KITER = 1.0e-5

ARMS = {
    "Tform": ["hydro/rad_x1_uform=false", "hydro/rad_x1_kiter=1"],
    "uform": ["hydro/rad_x1_uform=true", "hydro/rad_x1_kiter=1"],
    "uform_k2": ["hydro/rad_x1_uform=true", "hydro/rad_x1_kiter=2"],
}


def run(binary, arm, n):
    """run one arm at nx1 = nx2 = n; returns (L1 energy error, stdout)."""
    rundir = os.path.join(BUILD, "run_%s_%d" % (arm, n))
    os.makedirs(rundir, exist_ok=True)
    args = ["mesh/nx1=%d" % n, "mesh/nx2=%d" % n,
            "meshblock/nx1=%d" % n, "meshblock/nx2=%d" % n,
            "time/cfl_number=%.10g" % (CFL32 * 32.0 / n)] + ARMS[arm]
    proc = subprocess.run([binary, "-i", INPUT] + args, cwd=rundir,
                          capture_output=True, text=True)
    if proc.returncode != 0:
        pytest.fail("%s at nx=%d failed:\n%s" % (arm, n, proc.stdout[-4000:]))
    errs = os.path.join(rundir, "sts3d-errs.dat")
    assert os.path.exists(errs), "no error file was written"
    row = [ln for ln in open(errs).read().splitlines() if not ln.startswith("#")][-1]
    return float(row.split()[-1]), proc.stdout


def parse_report(stdout):
    """(worst conservation residual, worst fallback count) of the radial solve."""
    viol, bad = [], []
    for line in stdout.splitlines():
        if "### rad_implicit_x1 rank" not in line:
            continue
        viol.append(float(line.split("sum V|de| = ")[1].split(",")[0]))
        bad.append(int(line.split("fallback cells = ")[1].split(",")[0]))
    assert viol, "the radial solve never reported: is rad_implicit_x1 on?"
    return max(viol), max(bad)


def test_run():
    try:
        subprocess.run(["cmake", "-S", REPO, "-B", BUILD,
                        "-D", "CMAKE_BUILD_TYPE=Release"],
                       check=True, capture_output=True, text=True)
        subprocess.run(["make", "-C", BUILD, "-j", str(os.cpu_count()), "athena"],
                       check=True, capture_output=True, text=True)
        binary = os.path.join(BUILD, "src", "athena")
        l1 = {}
        for arm in ARMS:
            cmax = CONS_MAX_KITER if arm == "uform_k2" else CONS_MAX
            for n in (32, 64):
                l1[(arm, n)], out = run(binary, arm, n)
                viol, nbad = parse_report(out)
                assert nbad == 0, \
                    f"{arm} nx={n}: {nbad} cells fell back to a non-conservative update"
                assert viol < cmax, \
                    f"{arm} nx={n}: the radial solve is not conservative, " \
                    f"|sum V de - expected|/sum V|de| = {viol:.3g}"
                assert l1[(arm, n)] < L1_MAX[n], \
                    f"{arm} nx={n}: L1 = {l1[(arm, n)]:.6e}, expected < {L1_MAX[n]:.3g}"
            rate = math.log2(l1[(arm, 32)] / l1[(arm, 64)])
            assert rate > ORDER_MIN, \
                f"{arm}: convergence order {rate:.3f}, expected > {ORDER_MIN}"
        # the two forms are the same operator on a weak perturbation
        for arm in ("uform", "uform_k2"):
            for n in (32, 64):
                d = abs(l1[(arm, n)] - l1[("Tform", n)]) / l1[("Tform", n)]
                assert d < AGREE_MAX, \
                    f"{arm} nx={n}: L1 differs from the T-form by {d:.3g}, " \
                    f"expected < {AGREE_MAX:.3g}: the two forms should coincide where " \
                    f"T is nearly uniform across a face"
    finally:
        shutil.rmtree(BUILD, ignore_errors=True)
