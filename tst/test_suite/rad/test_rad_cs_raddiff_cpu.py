"""
Radiative diffusion on the cubed sphere: the angular operator against the exact
Laplacian of an l = 2 temperature field (cs_test iprob 15, inputs/tests/
cubed_sphere_raddiff.athinput).

Uniform density at rest, T = T0 (1 + amp f) with f an l = 2 combination in the direction
cosines and no radial dependence, so the radial fluxes vanish identically and the energy
source is the discrete ANGULAR operator alone, which must equal -6 kappa (T - T0) r_fac
with r_fac the discrete 1/r^2 the area-over-volume weighting produces.  The source is
read off two dumps 20 steps apart.  Gates: the fitted amplitude (corrected for the
field's own decay over the run) within 3 % of 1 at both resolutions and the residual
shape error small.  The orthogonal form of the face-normal derivative
(hydro/rad_cs_exact = false) fails these by 24-30 % in L1, which is what this guards.
"""

# Modules
import glob
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "vis", "python"))
import bin_convert  # noqa: E402

INPUT = os.path.join(REPO, "inputs", "tests", "cubed_sphere_raddiff.athinput")
BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_cs_raddiff")

# the input's units: T0 = 1 -> mu m_u v^2/k_B, p0 = 1 -> rho0 v^2 (the code uses the
# atomic mass unit, so does this)
L_CGS, M_CGS, T_CGS, MU = 1.0e8, 1.0e18, 527.0, 2.3
AMU, KB, SIGMA = 1.66053906660e-24, 1.380649e-16, 5.670374419e-5
KFAC, GAMMA = 0.01, 1.666667


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


def kappa_code():
    """kappa_rad at (T0, p0, rho0) in code units."""
    v = L_CGS / T_CGS
    dens = M_CGS / L_CGS ** 3
    pres = dens * v ** 2
    temp = MU * AMU * v ** 2 / KB
    kap = 16.0 * SIGMA * temp ** 3 / (3.0 * KFAC * freedman(temp, pres) * dens)
    return kap / (pres * v * L_CGS / temp)


def run(binary, rundir, n, exact):
    os.makedirs(rundir, exist_ok=True)
    args = ["mesh/nx2=%d" % n, "mesh/nx3=%d" % n, "meshblock/nx2=%d" % n,
            "meshblock/nx3=%d" % n,
            "hydro/rad_cs_exact=%s" % ("true" if exact else "false")]
    proc = subprocess.run([binary, "-i", INPUT] + args, cwd=rundir,
                          capture_output=True, text=True)
    if proc.returncode != 0:
        pytest.fail("run failed:\n" + proc.stdout[-4000:])
    return proc.stdout


def measure(rundir):
    """(decay-corrected amplitude ratio, L1 residual, Linf residual) of the source."""
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
    kap = kappa_code()
    expect = -6.0 * kap * (temp - 1.0) * rfac[None, None, None, :]
    slope = np.sum(src * expect) / np.sum(expect * expect)
    resid = src - slope * expect
    # the l = 2 field decays at 6 kappa r_fac/(rho c_v) during the run; the measured
    # source is its run average
    lam = 6.0 * kap * np.mean(rfac) / (1.0 / (GAMMA - 1.0))
    decay = (1.0 - np.exp(-lam * dt)) / (lam * dt)
    return (slope / decay, np.abs(resid).sum() / np.abs(expect).sum(),
            np.abs(resid).max() / np.abs(expect).max())


def test_run():
    try:
        subprocess.run(["cmake", "-S", REPO, "-B", BUILD, "-D", "PROBLEM=cs_test",
                        "-D", "CMAKE_BUILD_TYPE=Release"],
                       check=True, capture_output=True, text=True)
        subprocess.run(["make", "-C", BUILD, "-j", str(os.cpu_count())],
                       check=True, capture_output=True, text=True)
        binary = os.path.join(BUILD, "src", "athena")
        for n in (16, 32):
            amp, l1, linf = measure_run(binary, n, True)
            assert abs(amp - 1.0) < 0.03, \
                f"nx={n}: angular operator amplitude {amp:.4f}, expected 1 within 3 %"
            assert l1 < 0.05 and linf < 0.1, \
                f"nx={n}: residual shape error L1 {l1:.3g}, Linf {linf:.3g}"
        # the orthogonal form must be visibly wrong, or this test guards nothing
        amp, l1, linf = measure_run(binary, 16, False)
        assert l1 > 0.1, f"the orthogonal form has L1 residual {l1:.3g}; gate is vacuous"
    finally:
        shutil.rmtree(BUILD, ignore_errors=True)


def measure_run(binary, n, exact):
    rundir = os.path.join(BUILD, "run_%d_%s" % (n, "exact" if exact else "ortho"))
    run(binary, rundir, n, exact)
    return measure(rundir)
