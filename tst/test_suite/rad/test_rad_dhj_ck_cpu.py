"""
Regression test for the correlated-k radiative transfer in the deep_hot_jupiter_rt
problem generator.

This test is unlike the rest of the suite in two ways, both forced by what it covers:

  1. The scheme lives in a USER problem generator, so it needs its own binary built with
     -D PROBLEM=deep_hot_jupiter_rt.  The shared build that run_test_suite.py makes does
     not contain it, and a test cannot ask that build for a different PROBLEM, so this
     file configures and compiles a second binary into tst/build_ck and removes it again.
     That build is the expensive part of the test.

  2. It needs the Exo-FMS correlated-k tables, which are not in git (no upstream licence;
     data/exo_fms_ck/PROVENANCE.md records where to fetch them).  Without them the test
     SKIPS rather than fails, which is what happens in CI.

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
import re
import shutil
import subprocess
import numpy as np
import pytest

SIGMA_SB = 5.6704e-5                       # Stefan-Boltzmann, cgs

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
DATA = os.path.join(REPO, "data", "exo_fms_ck")
KTABLE = os.path.join(DATA, "ck", "Premixed_1x_g8_11.txt")
BUILD = os.path.join(REPO, "tst", "build_ck")
INPUT = os.path.join(REPO, "inputs", "mhd", "deep_hot_jupiter_rt_eos.athinput")

# Smallest grid the problem will run on: nx1 must equal the meshblock nx1 (the two-stream
# RT sweeps whole radial columns), and the polar boundary needs an even number of
# meshblocks in x3.  Two blocks of 64 x 8 x 4.
MESH = [
    "mesh/nx1=64", "mesh/nx2=8", "mesh/nx3=8",
    "meshblock/nx1=64", "meshblock/nx2=8", "meshblock/nx3=4",
    "time/nlim=1",
]

CK = [
    "problem/rt_ck=true",
    "problem/ck_table=" + KTABLE,
    "problem/ck_data_dir=" + DATA,
]

# phi indices of a nightside and a dayside column on this grid (mu0 = -0.92 and +0.92).
NIGHT_K = 2
DAY_K = 5


def build():
    """Configure and compile a deep_hot_jupiter_rt binary into tst/build_ck."""
    subprocess.run(
        ["cmake", "-S", REPO, "-B", BUILD, "-D", "PROBLEM=deep_hot_jupiter_rt",
         "-D", "CMAKE_BUILD_TYPE=Release"],
        check=True, capture_output=True, text=True,
    )
    subprocess.run(
        ["make", "-C", BUILD, "-j", str(os.cpu_count())],
        check=True, capture_output=True, text=True,
    )
    return os.path.join(BUILD, "src", "athena")


def run(binary, dumpfile, dump_k):
    """Run one cycle of the correlated-k scheme, dumping one radial column."""
    rundir = os.path.join(BUILD, "run")
    os.makedirs(rundir, exist_ok=True)
    proc = subprocess.run(
        [binary, "-i", INPUT] + MESH + CK
        + ["problem/ck_dump_file=" + dumpfile, "problem/ck_dump_k=" + repr(dump_k)],
        cwd=rundir, capture_output=True, text=True,
    )
    if proc.returncode != 0:
        pytest.fail(f"correlated-k run failed (k = {dump_k}):\n{proc.stdout[-4000:]}")
    return proc.stdout, os.path.join(rundir, dumpfile)


def read_column(path):
    """Return (header dict, columns) from a correlated-k column dump."""
    with open(path) as f:
        head = "".join([ln for ln in f if ln.startswith("#")])
    keys = ["mu0", "icut", "T_int", "T_irr"]
    hdr = {}
    for k in keys:
        m = re.search(k + r" = *([-0-9.eE+]+)", head)
        assert m is not None, f"'{k}' missing from the column dump header"
        hdr[k] = float(m.group(1))
    hdr["icut"] = int(hdr["icut"])
    data = np.loadtxt(path)
    return hdr, data


def grab(stdout, pattern, name):
    """Pull one float out of the run's diagnostic output."""
    m = re.search(pattern, stdout)
    if m is None:
        pytest.fail(f"{name} not reported by the run; the diagnostic may have moved")
    return float(m.group(1))


@pytest.mark.skipif(
    not os.path.exists(KTABLE),
    reason="Exo-FMS correlated-k tables absent; see data/exo_fms_ck/PROVENANCE.md",
)
def test_run():
    """Build the problem generator, then check the correlated-k RT on two columns."""
    try:
        binary = build()
        out, night = run(binary, "col_night.txt", NIGHT_K)
        _, day = run(binary, "col_day.txt", DAY_K)

        # --- exact limits the scheme tests on the device at startup ------------------
        iso = grab(out, r"isothermal net flux \|F\|/sigmaT\^4 <= *([0-9.eE+-]+)",
                   "isothermal self-test")
        thin = grab(out, r"transparent slab F/sigmaT\^4 = *([0-9.eE+-]+)",
                    "transparent-slab self-test")
        assert iso < 1.0e-12, f"isothermal atmosphere carries net flux {iso:g} sigma T^4"
        assert abs(thin - 1.0) < 1.0e-10, f"transparent slab emits {thin:g} sigma T^4"

        planck = grab(out, r"worst \|sum_b f_b - 1\| = *([0-9.eE+-]+)",
                      "Planck fraction sum")
        assert planck < 1.0e-12, f"band Planck fractions sum to 1 + {planck:g}"

        # --- the run is the scheme we think it is ------------------------------------
        assert "88 column solves" in out, \
            "expected 11 bands x 8 g-points x 1 angle = 88 chains"
        assert "rt_ck implies problem/rt_split" in out and "RT split path ON" in out, \
            "correlated-k did not take the split path; it may be running grey (35173e28)"

        # --- nightside: the bottom boundary delivers the interior flux ---------------
        hdr, col = read_column(night)
        assert hdr["mu0"] < 0.0, "column k = %d is not on the nightside" % NIGHT_K
        assert np.all(np.isfinite(col)), "non-finite values in the nightside column"
        p, temp, flw, qsw = col[:, 2], col[:, 3], col[:, 4], col[:, 5]
        assert np.all(temp > 0.0) and np.all(p > 0.0), "unphysical column"
        assert np.all(np.diff(p) < 0.0), "pressure does not fall outwards"
        assert np.all(qsw == 0.0), "shortwave heating on the nightside"
        assert np.all(flw[:hdr["icut"] - int(col[0, 0])] == 0.0), \
            "longwave flux below the correlated-k cutoff"
        fint = SIGMA_SB * hdr["T_int"] ** 4
        fbase = flw[hdr["icut"] - int(col[0, 0])]
        assert abs(fbase / fint - 1.0) < 0.2, \
            f"net flux at the cutoff is {fbase / fint:g} sigma T_int^4, not ~1"

        # --- dayside: the stellar sweep deposits the insolation ----------------------
        hdr, col = read_column(day)
        assert hdr["mu0"] > 0.0, "column k = %d is not on the dayside" % DAY_K
        assert np.all(np.isfinite(col)), "non-finite values in the dayside column"
        rf, qsw = col[:, 1], col[:, 5]
        assert np.all(qsw >= 0.0), "the shortwave cools somewhere"
        absorbed = np.sum(qsw[:-1] * np.diff(rf))
        incident = hdr["mu0"] * SIGMA_SB * hdr["T_irr"] ** 4
        assert 0.9 < absorbed / incident < 1.0, \
            f"shortwave absorbs {absorbed / incident:g} of the incident flux"
    finally:
        shutil.rmtree(BUILD, ignore_errors=True)
