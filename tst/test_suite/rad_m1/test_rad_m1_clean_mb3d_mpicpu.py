"""
Clean implicit M1 (cell layout, kinetic consistency, vet_sc) on a 3-D PERIODIC mesh
split into several MeshBlocks: a hot ellipsoid radiating into a cooler medium with gas
feedback, 3 cycles.

The run is made on one MeshBlock (1 rank), on 8 MeshBlocks of one rank, on 2x1x1
MeshBlocks of 2 ranks and on 2x2x2 MeshBlocks of 8 ranks.  The multi-block sweep
(VetSweepMB: periodic x1 passes, x2/x3 band wrap) is exact, so the final state of each
split agrees with the single block to the accuracy of the linear solves (their dot
products are summed in another order).  Total energy (gas + radiation) and momentum
(gas + F/c^2) are conserved to round-off in every run.
"""

# Modules
import glob
import shutil
import numpy as np
import test_suite.testutils as testutils
import bin_convert

input_file = "inputs/rad_m1_clean_mb3d.athinput"
CASES = [   # (ranks, MeshBlock nx1, nx2, nx3)
    (1, 16, 16, 16), (1, 8, 8, 8), (2, 8, 16, 16), (8, 8, 8, 8)]


def read(path):
    """The global arrays of a binary dump (MeshBlocks placed by logical location)."""
    fd = bin_convert.read_binary(path)
    ll, out = fd["mb_logical"], {}
    for name, blocks in fd["mb_data"].items():
        nk, nj, ni = np.asarray(blocks[0]).shape
        a = np.zeros(((ll[:, 2].max() + 1)*nk, (ll[:, 1].max() + 1)*nj,
                      (ll[:, 0].max() + 1)*ni))
        for b, (i, j, k, _) in zip(blocks, ll):
            a[k*nk:(k + 1)*nk, j*nj:(j + 1)*nj, i*ni:(i + 1)*ni] = b
        out[name] = a
    return out


def totals(h, r, c):
    """Total energy and the three momentum components of one output time."""
    e = h["ener"].sum() + r["m1_e"].sum()
    p = [h["mom%d" % n].sum() + r["m1_f%d" % n].sum()/c**2 for n in (1, 2, 3)]
    return e, np.array(p), np.abs(h["mom1"]).sum()


def run_case(ranks, n1, n2, n3):
    """One split; returns its last state and its energy and momentum drifts."""
    shutil.rmtree("bin", ignore_errors=True)
    flags = ["meshblock/nx1=%d" % n1, "meshblock/nx2=%d" % n2, "meshblock/nx3=%d" % n3]
    testutils.mpi_run(input_file, flags, threads=ranks)
    hs = sorted(glob.glob("bin/clean_mb3d.hydro_u.*.bin"))
    rs = sorted(glob.glob("bin/clean_mb3d.m1.*.bin"))
    assert len(hs) == 2 and len(rs) == 2, "the first and last dumps were not written"
    first = [read(f) for f in (hs[0], rs[0])]
    last = [read(f) for f in (hs[-1], rs[-1])]
    e0, p0, _ = totals(*first, 100.0)
    e1, p1, pg = totals(*last, 100.0)
    shutil.rmtree("bin", ignore_errors=True)
    return last, abs(e1 - e0)/e0, np.abs(p1 - p0).max()/pg


def test_clean_mb3d_periodic():
    """Single vs multi-block agreement and conservation on a 3-D periodic mesh."""
    try:
        res = {case: run_case(*case) for case in CASES}
        ref = res[CASES[0]][0]
        for case, (state, de, dp) in res.items():
            assert de < 1.0e-14, f"{case}: total energy drifts by {de:g}"
            assert dp < 1.0e-11, f"{case}: total momentum drifts by {dp:g}"
            for a, b in zip(state, ref):
                for name in b:
                    dev = np.abs(a[name] - b[name]).max()/np.abs(b[name]).max()
                    assert dev < 1.0e-9, f"{case}: {name} differs by {dev:g}"
    finally:
        testutils.cleanup()
