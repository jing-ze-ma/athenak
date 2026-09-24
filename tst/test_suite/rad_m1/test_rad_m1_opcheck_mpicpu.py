"""
Operator-equivalence test of the implicit M1 Krylov operator on 2 and 4 ranks.

<rad_m1>/implicit_op_check = 2 (rad_m1_opcheck.cpp) applies, at the first two implicit
solves, every operator variant the configuration allows -- the stored stencil, the
overlapped interior + shell of implicit_halo_overlap (with and without
implicit_halo_ovl_faces), the legacy 7-point + M1VimpRow + off-diagonal kernels, each
halo path (ordinary exchange, implicit_halo_mpi), with implicit_vimp_fold on and off --
to the same pseudo-random vector, checks the ghost zones against the neighbours, and
exits with a fatal error if any variant differs from the reference by more than 1e-12
of its row.  The overlap + vimp_fold row counted the x2/x3 +-1 vimp terms twice before
m1-sync (89174004): that gives 1e-10 here and failed (tests_m1/gates/README.md), while
the solver-tolerance gate of the time let it through.
"""

# Modules
import os
import re
import pytest
import test_suite.rad_m1.m1_common as m1

CASES = [
    # (ranks, meshblock nx2, closure)
    (2, 16, m1.VET),     # vet_sc: non-zero edge coefficients
    (4, 8, []),          # Eddington: the configuration of the old bug
]


@pytest.mark.skipif(not m1.HAVE_DATA, reason=m1.NO_DATA)
@pytest.mark.skipif(not m1.HAVE_MPI, reason="mpirun not on PATH")
def test_run():
    """Every operator variant must give the reference y = A x to 1e-12 per row."""
    exe = m1.build(True)
    for ranks, nx2, clo in CASES:
        rundir = os.path.join(m1.BUILD[True], "opcheck_%d" % ranks)
        out = m1.run(exe, rundir, [
            "meshblock/nx2=%d" % nx2, "time/nlim=1", "rad_m1/implicit_op_check=2",
            "rad_m1/implicit_halo_mpi=true", "rad_m1/implicit_halo_overlap=true",
            "rad_m1/implicit_halo_ovl_faces=true"] + m1.LEVERS + clo, ranks=ranks)
        assert "M1OPCHK solve 2: PASS" in out, out[-3000:]
        assert not re.search(r"(  FAIL$|: FAIL)", out, re.M), out[-3000:]
        for v in ["mpi/overlap ", "mpi/overlap_faces", "mpi+nofold/overlap ",
                  "exch/legacy_7pt", "halo mpi"]:
            assert v in out, f"variant {v.strip()} was not checked ({ranks} ranks)"
