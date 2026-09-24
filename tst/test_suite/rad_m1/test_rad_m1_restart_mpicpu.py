"""
Restart test of the implicit M1 He slab on 2 ranks: N cycles in one go must equal N/2
cycles + restart + N/2 cycles BITWISE: the restart data after the parameter dump, and
every history row the two runs share.  The dt column is left out of the row test: the
final row of a leg repeats the time of the row before it with the dt of the next step
(dt itself is in the restart data).
The configuration is the multi-rank be default since m1-accmerge (vimp_fold, one_pass,
predictor order 2, halo_mpi, overlap) with implicit_halo_ovl_faces; the one_pass
contraction state and the predictor history travel in the restart file, so a missing
field shows up as a difference.
"""

# Modules
import glob
import os
import pytest
import test_suite.rad_m1.m1_common as m1

N = 16
ARGS = ["meshblock/nx2=16", "rad_m1/implicit_halo_mpi=true",
        "rad_m1/implicit_halo_overlap=true",
        "rad_m1/implicit_halo_ovl_faces=true"] + m1.LEVERS


@pytest.mark.skipif(not m1.HAVE_DATA, reason=m1.NO_DATA)
@pytest.mark.skipif(not m1.HAVE_MPI, reason="mpirun not on PATH")
def test_run():
    """N cycles == N/2 + restart + N/2, bitwise."""
    exe = m1.build(True)
    a = os.path.join(m1.BUILD[True], "rst_a")
    b = os.path.join(m1.BUILD[True], "rst_b")
    m1.run(exe, a, ARGS + ["time/nlim=%d" % N], ranks=2)
    m1.run(exe, b, ARGS + ["time/nlim=%d" % (N // 2)], ranks=2)
    rst = sorted(glob.glob(os.path.join(b, "rst", "*.rst")))[-1]
    m1.run(exe, b, ["time/nlim=%d" % N], ranks=2, restart=rst)
    ra = sorted(glob.glob(os.path.join(a, "rst", "*.rst")))[-1]
    rb = sorted(glob.glob(os.path.join(b, "rst", "*.rst")))[-1]
    assert m1.payload(ra) == m1.payload(rb), "the final restart data differ"
    for name in ["m1slab.hydro.hst", "m1slab.user.hst"]:
        _, rows_a = m1.history(a, name)
        _, rows_b = m1.history(b, name)
        by_t = {r[0]: r for r in rows_a}
        common = [r for r in rows_b if r[0] in by_t]
        assert len(common) >= N - 1, f"{name}: only {len(common)} rows in common"
        by_t = {r[0]: r[:1] + r[2:] for r in rows_a}
        for r in common:
            r, ref = r[:1] + r[2:], by_t[r[0]]
            assert r == ref, f"{name}: the row at t = {r[0]} differs"
