"""
Regression test for the explicit-RT source limiter (problem/rt_de_max) in the
deep_hot_jupiter_rt problem generator.

The radiation is operator split and applied explicitly at the hydrodynamic timestep, so
it is stable only while the local radiative time e/|src| exceeds that step.  Nothing in
the scheme enforces it, and floors that pin the top of the atmosphere can drive the ratio
to 0.01 -- which is how an ideal-gas correlated-k run reached adjacent cells at 650, 2081
and 9806 K and then NaN'd.  LimitRTSource caps each update at rt_de_max * e_int.

What is checked:

  * THE IDENTITY PROPERTY, which is the whole reason the limiter is a hard clamp rather
    than a smooth function: on a healthy setup, running with the limiter on must give
    BITWISE the same history as running with it disabled.  A limiter that perturbed the
    valid regime would invalidate the Exo-FMS validation in docs/correlated_k_rt.md.
    Checked on both the correlated-k and the grey path, which reach the clamp through
    different kernels (rt_apply vs the monolithic 2stream_rt).
  * that a healthy run stays silent -- no warning is emitted when nothing is clipped;
  * the plumbing, by forcing every cell to clip with an absurdly small rt_de_max: the
    warning must appear, exactly once, and the run must still complete.

Not checked here: that the limiter rescues a real blow-up.  That takes ~40 min on a GPU
and half a rotation of simulated time; see docs/correlated_k_rt.md for the measurement.
"""

# Modules
import os
import shutil
import pytest
import test_suite.rad.dhj_ck_common as ck

BUILD = os.path.join(ck.REPO, "tst", "build_ck")

# Small grid, and enough cycles that the RT has acted many times.  The history is written
# every cycle so the two runs can be compared row by row.
MESH = ck.mesh(8, 8, 8, 4, 20)
EVERY_CYCLE = ["output1/dt=1.0e-30"]

WARNING = "the explicit radiative source was clipped"


def history(binary, tag, args):
    """Run once in its own directory and return the history file's contents."""
    rundir = os.path.join(BUILD, "run_" + tag)
    shutil.rmtree(rundir, ignore_errors=True)
    out = ck.run(binary, rundir, MESH + EVERY_CYCLE + args)
    hst = [f for f in os.listdir(rundir) if f.endswith(".hst")]
    assert len(hst) == 1, f"expected one history file in {rundir}, found {hst}"
    with open(os.path.join(rundir, hst[0])) as f:
        return out, f.read()


def identity(binary, tag, scheme):
    """The limiter must not move the answer on a setup that never trips it."""
    on, hon = history(binary, tag + "_on", scheme + ["problem/rt_de_max=0.5"])
    off, hoff = history(binary, tag + "_off", scheme + ["problem/rt_de_max=-1.0"])
    assert WARNING not in on, (
        f"the {tag} setup clipped with rt_de_max = 0.5; it is supposed to be healthy, "
        "so either the setup or the limiter has changed"
    )
    assert hon == hoff, (
        f"rt_de_max = 0.5 changed the {tag} history relative to the disabled limiter. "
        "The clamp must be the identity wherever |de| < rt_de_max * e_int."
    )
    return off


@pytest.mark.skipif(not ck.HAVE_TABLES, reason=ck.NO_TABLES)
def test_run():
    """Build the problem generator, then check the source limiter."""
    try:
        binary = ck.build(BUILD)

        # 1. identity, on both the correlated-k and the grey path
        identity(binary, "ck", ck.CK)
        identity(binary, "grey", ["problem/rt_ck=false", "problem/rt_split=false"])

        # 2. plumbing: an absurd cap clips every cell, and must warn exactly once
        out, _ = history(binary, "clip", ck.CK + ["problem/rt_de_max=1.0e-12"])
        assert out.count(WARNING) == 1, (
            "rt_de_max = 1e-12 must clip and warn exactly once; "
            f"the warning appeared {out.count(WARNING)} time(s)"
        )
    finally:
        shutil.rmtree(BUILD, ignore_errors=True)
