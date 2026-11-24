"""
The IMPLICIT transverse radiative operator on the CUBED SPHERE (hydro/rad_implicit_ang,
src/diffusion/conduction_transverse.cpp; design note docs/dev/cs_implicit_transverse.md).

Same problem and same measurement as test_rad_cs_raddiff_cpu.py -- cs_test iprob 15,
inputs/tests/cubed_sphere_raddiff.athinput: a uniform density at rest with
T = T0 (1 + amp f), f a spherical harmonic of degree l with no radial dependence, so the
radial fluxes vanish identically and the ANGULAR operator is the only thing that moves
energy.  The source read off two dumps must equal -l(l+1) kappa (T - T0) r_fac, i.e. the
harmonic must decay as exp(-l(l+1) D t/r^2); the fitted amplitude is the rate error.

What this guards that the explicit test does not: the cell VOLUME, the metric CROSS term
carried explicitly through cap_g2/cap_g3, and the PANEL SEAM -- which is an open face for
rad_ang_solver = sts and, for = adi, a face no tridiagonal line can cross and which is
applied by the pair-implicit seam sub-step instead.  Dropping any of the three shows up
here as a rate error of tens of per cent, or as a seam-localised residual.

The time step is scaled as dx^2 between the two resolutions, so the operator-splitting
error (this operator runs after the RK update) is the same at both and what is compared
is the space discretisation.
"""

# Modules
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "vis", "python"))
import bin_convert  # noqa: E402

INPUT = os.path.join(REPO, "inputs", "tests", "cubed_sphere_raddiff.athinput")
BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_cs_implicit")

L_CGS, M_CGS, T_CGS, MU = 1.0e8, 1.0e18, 527.0, 2.3
AMU, KB, SIGMA = 1.66053906660e-24, 1.380649e-16, 5.670374419e-5
KFAC, GAMMA = 0.01, 1.666667

# (nx per panel edge, cfl, nlim): dt ~ dx^2 at a fixed final time
LADDER = ((16, "0.03", "20"), (32, "0.0075", "80"))


def freedman(tk, p):
    """Rosseland mean of Freedman et al. 2014 (the fit in src/utils/rosseland.hpp)."""
    t1 = min(max(tk, 75.0), 4000.0)
    p1 = min(max(p, 1.0), 3.0e8)
    lt, lp = np.log10(t1), np.log10(p1)
    c1, c2, c3, c4, c5, c7 = 10.602, 2.882, 6.09e-15, 2.954, -2.526, -5.490
    if t1 < 800.0:
        c8, c9, c10, c11, c12 = -14.051, 3.055, 0.024, 1.877, -0.445
    else:
        c8, c9, c10, c11, c12 = 82.241, -55.456, 8.754, 0.7048, -0.0414
    lkl = c1 * np.arctan(lt - c2) - c3 / (lp + c4) * np.exp((lt - c5) ** 2) + c7
    lkh = c8 + c9 * lt + c10 * lt ** 2 + lp * (c11 + c12 * lt)
    return 10.0 ** lkl + 10.0 ** lkh


def kappa_code(kfac=KFAC):
    """kappa_rad at (T0, p0, rho0) in code units."""
    v = L_CGS / T_CGS
    dens = M_CGS / L_CGS ** 3
    pres = dens * v ** 2
    temp = MU * AMU * v ** 2 / KB
    kap = 16.0 * SIGMA * temp ** 3 / (3.0 * kfac * freedman(temp, pres) * dens)
    return kap / (pres * v * L_CGS / temp)


def run(binary, rundir, args):
    os.makedirs(rundir, exist_ok=True)
    proc = subprocess.run([binary, "-i", INPUT] + args, cwd=rundir,
                          capture_output=True, text=True)
    if proc.returncode != 0:
        pytest.fail("run failed:\n" + proc.stdout[-4000:] + proc.stderr[-2000:])
    return proc.stdout


def measure(rundir, lh):
    """(rate amplitude, L1 residual, panel-edge residual, panel-interior residual)."""
    files = sorted(glob.glob(os.path.join(rundir, "bin", "lap.hydro_w.*.bin")))
    assert len(files) >= 2, "no dumps"
    r0, r1 = bin_convert.read_binary(files[0]), bin_convert.read_binary(files[-1])
    dt = r1["time"] - r0["time"]
    e0, e1 = np.asarray(r0["mb_data"]["eint"]), np.asarray(r1["mb_data"]["eint"])
    rho = np.asarray(r0["mb_data"]["dens"])
    src = (e1 - e0) / dt
    temp = e0 * (GAMMA - 1.0) / rho
    g = np.asarray(r0["mb_geometry"])
    ni = e0.shape[-1]
    rf = g[0, 0] + (g[0, 1] - g[0, 0]) * np.arange(ni + 1) / ni
    rl, rr = rf[:-1], rf[1:]
    rfac = 1.5 * (rr ** 2 - rl ** 2) / ((rr ** 3 - rl ** 3) * 0.5 * (rl + rr))
    lam0 = lh * (lh + 1.0)
    kap = kappa_code()
    expect = -lam0 * kap * (temp - 1.0) * rfac[None, None, None, :]
    slope = np.sum(src * expect) / np.sum(expect * expect)
    resid = src - slope * expect
    lam = lam0 * kap * np.mean(rfac) / (1.0 / (GAMMA - 1.0))
    decay = (1.0 - np.exp(-lam * dt)) / (lam * dt)
    rn = np.abs(resid) / np.abs(expect).max()
    edge = np.zeros(rn.shape, dtype=bool)
    edge[:, 0, :, :] = edge[:, -1, :, :] = True
    edge[:, :, 0, :] = edge[:, :, -1, :] = True
    return (slope / decay, np.abs(resid).sum() / np.abs(expect).sum(),
            rn[edge].max(), rn[~edge].max())


def violation(stdout):
    """The largest |sum V de|/sum V|de| the operator itself reported."""
    v = [float(x) for x in
         re.findall(r"\|sum V de\|/sum V\|de\| = ([0-9.eE+-]+)", stdout)]
    return max(v) if v else 0.0


def arm(binary, tag, n, cfl, nlim, solver, lharm=2, extra=()):
    args = ["mesh/nx2=%d" % n, "mesh/nx3=%d" % n,
            "meshblock/nx2=%d" % n, "meshblock/nx3=%d" % n,
            "time/cfl_number=" + cfl, "time/nlim=" + nlim,
            "problem/lharm=%d" % lharm,
            "hydro/rad_implicit_ang=true", "hydro/rad_ang_solver=" + solver,
            "hydro/rad_ang_verbose=true"] + list(extra)
    rundir = os.path.join(BUILD, "run_" + tag)
    out = run(binary, rundir, args)
    return measure(rundir, lharm) + (violation(out),)


def test_run():
    try:
        subprocess.run(["cmake", "-S", REPO, "-B", BUILD, "-D", "PROBLEM=cs_test",
                        "-D", "CMAKE_BUILD_TYPE=Release"],
                       check=True, capture_output=True, text=True)
        subprocess.run(["make", "-C", BUILD, "-j", str(os.cpu_count())],
                       check=True, capture_output=True, text=True)
        binary = os.path.join(BUILD, "src", "athena")

        # ---- l = 2, both solvers, two resolutions
        for solver in ("sts", "adi"):
            prev = None
            for n, cfl, nlim in LADDER:
                amp, l1, edge, intr, viol = arm(binary, "%s_%d" % (solver, n),
                                                n, cfl, nlim, solver)
                assert abs(amp - 1.0) < 0.02, \
                    f"{solver} nx={n}: decay rate off by {100*abs(amp-1):.2f} %"
                assert l1 < 0.05, f"{solver} nx={n}: residual L1 {l1:.3g}"
                # the seam must not be a special place: a dropped cross term or a
                # mishandled panel face shows up here first
                assert edge < 3.0 * intr, \
                    f"{solver} nx={n}: panel-edge residual {edge:.3g} vs " \
                    f"interior {intr:.3g}"
                # conservation: exact to round-off inside a panel, to the along-seam
                # resample at a seam (which converges with resolution)
                assert viol < 1.0e-4, \
                    f"{solver} nx={n}: |sum V de|/sum V|de| = {viol:.3g}"
                if prev is not None:
                    assert abs(amp - 1.0) < prev, \
                        f"{solver}: the rate error did not improve with resolution"
                prev = abs(amp - 1.0)

        # ---- l = 6: four times the eigenvalue and a much shorter wavelength, which is
        # where the metric cross term is worth tens of per cent
        for solver in ("sts", "adi"):
            amp, l1, edge, intr, viol = arm(binary, "l6_" + solver, 32, "0.0075",
                                            "80", solver, lharm=6)
            assert abs(amp - 1.0) < 0.03, \
                f"{solver} l=6: decay rate off by {100*abs(amp-1):.2f} %"
            assert l1 < 0.05, f"{solver} l=6: residual L1 {l1:.3g}"

        # ---- STABILITY at a transverse diffusion number ~ 2e3 (rad_implicit_x1 as
        # well, or the radial conduction step would shrink dt with the conductivity).
        # The ADI takes one step however stiff the row is and must stay bounded and
        # monotone; RKL1 is given the substage count its own bound asks for.
        for solver, maxit in (("adi", "200"), ("sts", "2000")):
            rundir = os.path.join(BUILD, "run_stab_" + solver)
            run(binary, rundir,
                ["mesh/nx2=32", "mesh/nx3=32", "meshblock/nx2=32", "meshblock/nx3=32",
                 "time/nlim=40", "hydro/rad_kappa_fac=3.0e-3",
                 "hydro/rad_implicit_x1=true", "hydro/rad_implicit_ang=true",
                 "hydro/rad_ang_solver=" + solver, "hydro/rad_ang_maxit=" + maxit,
                 "hydro/rad_ang_verbose=true"])
            hist = amplitude(rundir)
            assert np.all(np.isfinite(hist)), f"{solver}: stiff arm went non-finite"
            assert np.max(np.abs(hist)) <= 1.0 + 1.0e-8, \
                f"{solver}: stiff arm GREW, max |A/A0| = {np.max(np.abs(hist)):.3g}"
            assert abs(hist[-1]) < 0.05, \
                f"{solver}: stiff arm barely relaxed, A/A0 = {hist[-1]:.3g}"
    finally:
        shutil.rmtree(BUILD, ignore_errors=True)


def amplitude(rundir):
    """A(t)/A(0) of the initial harmonic, projected on the initial field."""
    files = sorted(glob.glob(os.path.join(rundir, "bin", "lap.hydro_w.*.bin")))
    r0 = bin_convert.read_binary(files[0])
    f = (np.asarray(r0["mb_data"]["eint"]) * (GAMMA - 1.0)
         / np.asarray(r0["mb_data"]["dens"]) - 1.0)
    den = np.sum(f * f)
    out = []
    for fn in files:
        r = bin_convert.read_binary(fn)
        t = (np.asarray(r["mb_data"]["eint"]) * (GAMMA - 1.0)
             / np.asarray(r["mb_data"]["dens"]) - 1.0)
        out.append(np.sum(t * f) / den)
    return np.asarray(out)
