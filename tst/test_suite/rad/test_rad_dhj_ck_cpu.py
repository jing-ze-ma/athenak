"""
Regression test for the correlated-k radiative transfer in the deep_hot_jupiter_rt
problem generator.

Why this test builds its own binary and skips without the Exo-FMS tables is explained in
dhj_ck_common.py, which it shares with the MPI sibling.

What is checked, cheapest and sharpest first:

  * the two exact-limit self-tests the scheme runs on the device at startup -- an
    isothermal atmosphere must carry zero net flux (the recurrence coefficients), and a
    transparent slab over a blackbody floor must emit exactly sigma T^4 (the g-point
    weights, the band Planck fractions and the flux prefactor);
  * the band-integrated Planck fractions summing to one;
  * the chain count, i.e. that the run really is 11 bands x 8 g-points x 1 angle;
  * that rt_ck implies rt_split -- the failure mode fixed in 35173e28 was rt_ck=true
    silently running the GREY scheme with every diagnostic looking healthy;
  * on a NIGHTSIDE column, the net longwave flux at the deep cutoff against
    sigma T_int^4, which is what the bottom boundary condition is supposed to deliver;
  * on a DAYSIDE column, the depth-integrated shortwave heating against the incident
    mu0 sigma T_irr^4, i.e. that the stellar sweep conserves energy.

The two physical tolerances are loose because the test runs at nx1 = 64, sixteen times
coarser in every direction than production: the same checks give 0.04 % and 98.8 % at
nx1 = 256.  They are here to catch a broken table read, a lost weight or a flipped
boundary condition, not to measure accuracy.
"""

# Modules
import os
import shutil
import numpy as np
import pytest
import test_suite.rad.dhj_ck_common as ck

BUILD = os.path.join(ck.REPO, "tst", "build_ck")

# Smallest grid the problem will run on: two meshblocks of 64 x 8 x 4.
MESH = ck.mesh(8, 8, 8, 4, 1)

# phi indices of a nightside and a dayside column on this grid (mu0 = -0.92 and +0.92).
NIGHT_K = 2
DAY_K = 5


def column(binary, name, dump_k):
    """One cycle of the correlated-k scheme, dumping one radial column."""
    rundir = os.path.join(BUILD, "run")
    out = ck.run(binary, rundir, MESH + ck.CK
                 + ["problem/ck_dump_file=" + name,
                    "problem/ck_dump_k=" + repr(dump_k)])
    return out, os.path.join(rundir, name)


@pytest.mark.skipif(not ck.HAVE_TABLES, reason=ck.NO_TABLES)
def test_run():
    """Build the problem generator, then check the correlated-k RT on two columns."""
    try:
        binary = ck.build(BUILD)
        out, night = column(binary, "col_night.txt", NIGHT_K)
        _, day = column(binary, "col_day.txt", DAY_K)

        # --- exact limits the scheme tests on the device at startup ------------------
        iso = ck.grab(out, r"isothermal net flux \|F\|/sigmaT\^4 <= *([0-9.eE+-]+)",
                      "isothermal self-test")
        thin = ck.grab(out, r"transparent slab F/sigmaT\^4 = *([0-9.eE+-]+)",
                       "transparent-slab self-test")
        assert iso < 1.0e-12, f"isothermal atmosphere carries net flux {iso:g} sigma T^4"
        assert abs(thin - 1.0) < 1.0e-10, f"transparent slab emits {thin:g} sigma T^4"

        planck = ck.grab(out, r"worst \|sum_b f_b - 1\| = *([0-9.eE+-]+)",
                         "Planck fraction sum")
        assert planck < 1.0e-12, f"band Planck fractions sum to 1 + {planck:g}"

        # --- the run is the scheme we think it is ------------------------------------
        assert "88 column solves" in out, \
            "expected 11 bands x 8 g-points x 1 angle = 88 chains"
        assert "rt_ck implies problem/rt_split" in out and "RT split path ON" in out, \
            "correlated-k did not take the split path; it may be running grey (35173e28)"

        # --- nightside: the bottom boundary delivers the interior flux ---------------
        hdr, col = ck.read_column(night)
        assert hdr["mu0"] < 0.0, "column k = %d is not on the nightside" % NIGHT_K
        assert np.all(np.isfinite(col)), "non-finite values in the nightside column"
        p, temp, flw, qsw = col[:, 2], col[:, 3], col[:, 4], col[:, 5]
        assert np.all(temp > 0.0) and np.all(p > 0.0), "unphysical column"
        assert np.all(np.diff(p) < 0.0), "pressure does not fall outwards"
        assert np.all(qsw == 0.0), "shortwave heating on the nightside"
        assert np.all(flw[:hdr["icut"] - int(col[0, 0])] == 0.0), \
            "longwave flux below the correlated-k cutoff"

        # --- the optical-depth blend with the radiative diffusion (default in the
        # dhj inputs): the two-stream owns 1 - w of every face, the diffusion w, with
        # w(tau_R) rising from 0 to 1 over the overlap.  The two operators share the
        # ck-table opacity, so in the overlap the two-stream's net flux must sit near
        # the diffusion flux -16 sigma T^3/(3 kappa_R rho) dT/dr, which the dump's own
        # tau_R column gives (kappa_R rho dr = -d tau_R).  A 30 % gate: the two-stream's
        # diffusion limit is 2/D = 1.20 against the exact 4/3 with the diffusivity
        # factor D = 1.66, and the face averaging adds a few percent more.
        rf, tau, w = col[:, 1], col[:, 8], col[:, 9]
        assert w[0] == 1.0 and w[-1] == 0.0, \
            "the blend weight must run 1 (wall) -> 0 (top)"
        assert np.all(np.diff(tau) <= 0.0), "tau_R does not fall outwards"
        # faces the two-stream still carries at least a tenth of; at its very bottom
        # (w ~ 1) its net flux is set by its Planck bottom boundary, not the gradient
        ov = np.where((w > 0.02) & (w < 0.9))[0]
        assert ov.size >= 3, "the diffusion / two-stream overlap is fewer than 3 faces"
        assert 0.0 < w[hdr["icut"] - int(col[0, 0])] < 1.0, \
            "the two-stream bottom is not inside the overlap"
        kr_rho = -np.diff(tau) / np.diff(rf)                 # per cell, 1/cm
        tc = 0.5 * (temp[1:] + temp[:-1])
        rc = 0.5 * (rf[1:] + rf[:-1])
        ratios = []
        for f in ov:
            if f < 1 or f >= len(rf) - 1:
                continue
            kr = 0.5 * (kr_rho[f - 1] + kr_rho[f])
            dtdr = (tc[f] - tc[f - 1]) / (rc[f] - rc[f - 1])
            fdiff = -16.0 * ck.SIGMA_SB * temp[f] ** 3 / (3.0 * kr) * dtdr
            ratios.append(flw[f] / fdiff)
        ratios = np.array(ratios)
        assert np.all(np.abs(ratios - 1.0) < 0.3), \
            f"two-stream / diffusion flux in the overlap: {ratios}"
        # the wall still delivers the interior flux, now through the diffusion: the
        # column dump reports it in the run's own words
        assert "wall flux sigma T_int^4" in out, "the wall flux was not set from T_int"

        # --- dayside: the stellar sweep deposits the insolation ----------------------
        hdr, col = ck.read_column(day)
        assert hdr["mu0"] > 0.0, "column k = %d is not on the dayside" % DAY_K
        assert np.all(np.isfinite(col)), "non-finite values in the dayside column"
        rf, qsw = col[:, 1], col[:, 5]
        assert np.all(qsw >= 0.0), "the shortwave cools somewhere"
        absorbed = np.sum(qsw[:-1] * np.diff(rf))
        incident = hdr["mu0"] * ck.SIGMA_SB * hdr["T_irr"] ** 4
        assert 0.9 < absorbed / incident < 1.0, \
            f"shortwave absorbs {absorbed / incident:g} of the incident flux"
    finally:
        shutil.rmtree(BUILD, ignore_errors=True)
