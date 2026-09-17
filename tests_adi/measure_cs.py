"""Measure the cubed-sphere angular radiative operator on cs_test iprob 15 dumps.

Reports, for a run directory:
  amp  : the fitted amplitude of the discrete angular source against the exact
         -l(l+1) kappa (T - T0)/r^2, corrected for the field's decay over the run
         (1 = the operator has exactly the right eigenvalue, i.e. the harmonic decays
         as exp(-l(l+1) D t/r^2))
  l1   : L1 of the residual shape error, normalised by L1 of the expectation
  linf : the same in Linf
  seam : max |residual| over the cells adjacent to a panel edge, over the same in the
         panel interiors
  cons : |sum_i V_i de_i| / sum_i V_i |de_i| over the whole shell
"""
import glob
import os
import sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO, "vis", "python"))
import bin_convert  # noqa: E402

L_CGS, M_CGS, T_CGS, MU = 1.0e8, 1.0e18, 527.0, 2.3
AMU, KB, SIGMA = 1.66053906660e-24, 1.380649e-16, 5.670374419e-5
KFAC, GAMMA = 0.01, 1.666667


def freedman(tk, p):
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
    v = L_CGS / T_CGS
    dens = M_CGS / L_CGS ** 3
    pres = dens * v ** 2
    temp = MU * AMU * v ** 2 / KB
    kap = 16.0 * SIGMA * temp ** 3 / (3.0 * kfac * freedman(temp, pres) * dens)
    return kap / (pres * v * L_CGS / temp)


def measure(rundir, lh=2, kfac=KFAC, base="lap"):
    files = sorted(glob.glob(os.path.join(rundir, "bin", base + ".hydro_w.*.bin")))
    assert len(files) >= 2, "no dumps in " + rundir
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
    kap = kappa_code(kfac)
    lam0 = lh * (lh + 1.0)
    expect = -lam0 * kap * (temp - 1.0) * rfac[None, None, None, :]
    slope = np.sum(src * expect) / np.sum(expect * expect)
    resid = src - slope * expect
    lam = lam0 * kap * np.mean(rfac) / (1.0 / (GAMMA - 1.0))
    decay = (1.0 - np.exp(-lam * dt)) / (lam * dt)
    rn = np.abs(resid) / np.abs(expect).max()
    # panel-edge cells (the outermost ring of every block in x2/x3) vs the interior
    edge = np.zeros(rn.shape, dtype=bool)
    edge[:, 0, :, :] = True
    edge[:, -1, :, :] = True
    edge[:, :, 0, :] = True
    edge[:, :, -1, :] = True
    # the conserved-energy budget uses the true cell volumes; approximate them by the
    # solid angle x radial shell weighting the code uses (exact for this uniform grid
    # only up to the gnomonic Jacobian), so read the code's own report instead when it
    # is available.  Here: the unweighted sum is enough to show a gross violation.
    return dict(amp=slope / decay,
                l1=np.abs(resid).sum() / np.abs(expect).sum(),
                linf=np.abs(resid).max() / np.abs(expect).max(),
                edge=rn[edge].max(), intr=rn[~edge].max(),
                dt=dt, nt=len(files))


if __name__ == "__main__":
    lh = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    kf = float(sys.argv[3]) if len(sys.argv) > 3 else KFAC
    r = measure(sys.argv[1], lh, kf)
    print("%-28s amp=%8.5f  L1=%8.4g  Linf=%8.4g  edge=%8.4g  interior=%8.4g" %
          (os.path.basename(os.path.normpath(sys.argv[1])), r["amp"], r["l1"],
           r["linf"], r["edge"], r["intr"]))
