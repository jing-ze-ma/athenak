"""
Decomposition-independence test for the correlated-k radiative transfer, and the only
coverage in this suite of the POLAR boundary under MPI -- deep_hot_jupiter_rt is the only
problem in inputs/ that sets use_polar_boundary.

The correlated-k RT itself has no halo exchange: it is a per-column solve, and the radial
direction is never split (the meshblock nx1 must equal the mesh nx1), so a correct
implementation must give bit-identical answers however theta and phi are divided.  The
polar boundary around it does communicate, and that is where this bites.  Its azimuthal
EMF average reduces over all ranks through a pair of host mirrors that were left at
extent 0 for want of a realloc, so every multi-rank polar run aborted on the first cycle
between 006bae07 and the fix.  Serial never touched that branch.

The test therefore runs the same problem on one rank and on four and demands the history
output match byte for byte.  It builds its own MPI binary for the same reasons the CPU
test does -- see dhj_ck_common.py -- and skips without the Exo-FMS tables or mpirun.
"""

# Modules
import os
import shutil
import pytest
import test_suite.rad.dhj_ck_common as ck

BUILD = os.path.join(ck.REPO, "tst", "build_ck_mpi")

# 8 meshblocks of 64 x 4 x 4, so 1, 2, 4 and 8 ranks all divide the work differently.
MESH = ck.mesh(8, 16, 4, 4, 10)
RANKS = [1, 4]

HAVE_MPI = shutil.which("mpirun") is not None


@pytest.mark.skipif(not ck.HAVE_TABLES, reason=ck.NO_TABLES)
@pytest.mark.skipif(not HAVE_MPI, reason="mpirun not on PATH")
def test_run():
    """The answer must not depend on how theta and phi are divided between ranks."""
    try:
        binary = ck.build(BUILD, ["-D", "Athena_ENABLE_MPI=ON"])
        history = []
        for n in RANKS:
            rundir = os.path.join(BUILD, "run%d" % n)
            out = ck.run(binary, rundir, MESH + ck.CK + ["output1/dt=1.0"], ranks=n)
            assert "Terminating on cycle limit" in out, \
                f"the {n}-rank run did not reach the cycle limit"
            hst = os.path.join(rundir, "dhj.mhd.hst")
            assert os.path.exists(hst), f"the {n}-rank run wrote no history file"
            with open(hst, "rb") as f:
                history.append(f.read())
        for n, h in zip(RANKS[1:], history[1:]):
            if h != history[0]:
                pytest.fail(
                    f"{n} ranks and {RANKS[0]} rank disagree: the correlated-k RT is a "
                    "per-column solve and must be decomposition-independent"
                )
    finally:
        shutil.rmtree(BUILD, ignore_errors=True)
