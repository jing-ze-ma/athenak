"""
The IMPLICIT transverse radiative operator on the CUBED SPHERE (hydro/rad_implicit_ang,
src/diffusion/conduction_transverse.cpp; design note docs/dev/cs_implicit_transverse.md).

Same problem and same measurement as test_rad_cs_raddiff_cpu.py -- cs_test iprob 15,
tst/inputs/cubed_sphere_raddiff.athinput: a uniform density at rest with
T = T0 (1 + amp f), f a spherical harmonic of degree l with no radial dependence.  The
harmonic decays as R_i(t) on radial shell i, dR/dt = D (radial - l(l+1) r_fac) R: every
shell has its own rate and the shells conduct to each other, so R is the exact solution
of the code's own radial finite-volume system (no angular error in it).  The change of
e between two dumps is fitted to (R - 1) times the initial harmonic: rate error, residual.
kappa(T) makes the decay rate linear in amp: the runs come in pairs +-amp and the odd half
is measured.  rad_kappa_fac is small so that the sound waves do not respond.

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

INPUT = os.path.join(REPO, "tst", "inputs", "cubed_sphere_raddiff.athinput")
BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_cs_implicit")

L_CGS, M_CGS, T_CGS, MU = 1.0e8, 1.0e18, 527.0, 2.3
AMU, KB, SIGMA = 1.66053906660e-24, 1.380649e-16, 5.670374419e-5
KFAC, GAMMA, AMP, TLIM = 1.0e-3, 1.666667, 0.01, 5.0e-5

# (nx per panel edge, cfl): dt ~ dx^2 at a fixed final time; adi splits the seam faces
# from the lines and needs the smaller dt
LADDER = {"sts": ((16, "0.12"), (32, "0.03")), "adi": ((16, "0.03"), (32, "0.0075"))}


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


def dumps(rundir):
    """(time between the two dumps, e at the first, e at the last, x1 faces)."""
    files = sorted(glob.glob(os.path.join(rundir, "bin", "lap.hydro_w.*.bin")))
    assert len(files) >= 2, "no dumps"
    r0, r1 = bin_convert.read_binary(files[0]), bin_convert.read_binary(files[-1])
    e0 = np.asarray(r0["mb_data"]["eint"], dtype=float)
    e1 = np.asarray(r1["mb_data"]["eint"], dtype=float)
    g = np.asarray(r0["mb_geometry"])
    rf = g[0, 0] + (g[0, 1] - g[0, 0]) * np.arange(e0.shape[-1] + 1) / e0.shape[-1]
    return r1["time"] - r0["time"], e0, e1, rf


def decay(rf, lh, dt):
    """R(dt) - 1 per radial shell: the exact solution of the radial system of the code
    (centroid radii, insulated walls) for a harmonic of degree lh and constant kappa."""
    a, b, ni = rf[:-1], rf[1:], len(rf) - 1
    vol = (b ** 3 - a ** 3) / 3.0
    xc = 0.75 * (b ** 4 - a ** 4) / (b ** 3 - a ** 3)
    m = np.diag(-lh * (lh + 1.0) * 0.5 * (b ** 2 - a ** 2) / (vol * xc))
    for i in range(ni):
        for k in (i - 1, i + 1):
            if 0 <= k < ni:
                c = rf[max(i, k)] ** 2 / abs(xc[k] - xc[i]) / vol[i]
                m[i, k] += c
                m[i, i] -= c
    w, v = np.linalg.eig(m * kappa_code() * (GAMMA - 1.0) * dt)
    return np.real(v @ (np.exp(w) * np.linalg.solve(v, np.ones(ni)))) - 1.0


def measure(rundirs, lh):
    """(rate amplitude, L1 residual, residual max over the panel edges, over the 2 x 2
    cells at the cube vertices)."""
    (dt, ep0, ep1, rf), (dt2, em0, em1, _) = (dumps(d) for d in rundirs)
    assert abs(dt2 / dt - 1.0) < 1.0e-9, "the two dumps are at different times"
    expect = 0.5 * (ep0 - em0) * decay(rf, lh, dt)
    de = 0.5 * ((ep1 - ep0) - (em1 - em0))
    slope = np.sum(de * expect) / np.sum(expect * expect)
    resid = np.abs(de - expect)
    rn = resid / np.abs(expect).max(axis=(0, 1, 2))
    edge = np.zeros(rn.shape, dtype=bool)
    edge[:, [0, -1], :, :] = edge[:, :, [0, -1], :] = True
    d = np.minimum(np.arange(rn.shape[2]), rn.shape[2] - 1 - np.arange(rn.shape[2]))
    vert = (d[:, None] < 2) & (d[None, :] < 2)
    return (slope, resid.sum() / np.abs(expect).sum(), rn[edge].max(), rn[:, vert].max())


def violation(stdout):
    """The largest |sum V de|/sum V|de| the operator itself reported."""
    v = [float(x) for x in
         re.findall(r"\|sum V de\|/sum V\|de\| = ([0-9.eE+-]+)", stdout)]
    return max(v) if v else 0.0


def arm(binary, tag, n, cfl, solver, lharm=2):
    rundirs, viol = [], 0.0
    for sign in (1, -1):
        rundirs.append(os.path.join(BUILD, "run_%s_%d" % (tag, sign)))
        out = run(binary, rundirs[-1],
                  ["mesh/nx2=%d" % n, "mesh/nx3=%d" % n,
                   "meshblock/nx2=%d" % n, "meshblock/nx3=%d" % n,
                   "time/cfl_number=" + cfl, "time/nlim=-1", "time/tlim=%g" % TLIM,
                   "output1/dt=1.0e30", "problem/lharm=%d" % lharm,
                   "problem/amp=%g" % (sign * AMP), "hydro/rad_kappa_fac=%g" % KFAC,
                   "hydro/rad_implicit_ang=true", "hydro/rad_ang_solver=" + solver,
                   "hydro/rad_ang_verbose=true"])
        viol = max(viol, violation(out))
    return measure(rundirs, lharm) + (viol,)


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
            for n, cfl in LADDER[solver]:
                amp, l1, edge, vert, viol = arm(binary, "%s_%d" % (solver, n),
                                                n, cfl, solver)
                print("%s n %d rate-1 %.3e L1 %.3e edge %.3e vert %.3e viol %.1e"
                      % (solver, n, amp - 1.0, l1, edge, vert, viol))
                assert abs(amp - 1.0) < 0.005, \
                    f"{solver} nx={n}: decay rate off by {100*abs(amp-1):.2f} %"
                assert l1 < 0.005, f"{solver} nx={n}: residual L1 {l1:.3g}"
                # panel edges and cube vertices: a dropped cross term, a mishandled
                # panel face or a first-order corner value shows up here first
                assert edge < 0.02, f"{solver} nx={n}: panel-edge residual {edge:.3g}"
                assert vert < 0.02, f"{solver} nx={n}: cube-vertex residual {vert:.3g}"
                # conservation: the panel-seam faces carry one flux, shared by both sides
                assert viol < 1.0e-13, \
                    f"{solver} nx={n}: |sum V de|/sum V|de| = {viol:.3g}"
                if prev is not None:
                    assert abs(amp - 1.0) < prev[0] / 3.0, \
                        f"{solver}: the rate error did not fall at second order"
                    assert l1 < prev[1] / 2.0 and vert < prev[2] / 2.0, \
                        f"{solver}: the residual did not fall with resolution"
                prev = abs(amp - 1.0), l1, vert

        # ---- l = 6: seven times the eigenvalue and a much shorter wavelength, which is
        # where the metric cross term is worth tens of per cent
        for solver in ("sts", "adi"):
            amp, l1, edge, vert, viol = arm(binary, "l6_" + solver, 32, "0.03",
                                            solver, lharm=6)
            print("%s l 6 rate-1 %.3e L1 %.3e edge %.3e vert %.3e"
                  % (solver, amp - 1.0, l1, edge, vert))
            assert abs(amp - 1.0) < 0.01, \
                f"{solver} l=6: decay rate off by {100*abs(amp-1):.2f} %"
            assert l1 < 0.01, f"{solver} l=6: residual L1 {l1:.3g}"
            assert vert < 0.01, f"{solver} l=6: cube-vertex residual {vert:.3g}"
            assert viol < 1.0e-13, f"{solver} l=6: |sum V de|/sum V|de| = {viol:.3g}"

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
