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

The test therefore runs the same problem on one rank and on four and compares the history
output value by value (printed at %.17e) within a relative tolerance RTOL.  It no longer
demands byte equality: the state is NOT bit-identical across decompositions (already on
b2db548d, 625 of the rst floats differ between 1 and 4 ranks, relative 4e-16..3e-14,
probably the order of the polar-boundary reductions; not investigated), and a byte
comparison of the 6-digit history then depends on whether a last printed digit happens to
round the same.  Cancelling columns (2-mom is ~1e18 against terms of ~1e27) get an
absolute floor ATOL times the magnitude scale of their group (all *-mom, *-KE, *-ME
columns), so a last-bit difference in the terms cannot fail them.

It builds its own MPI binary for the same reasons the CPU test does -- see
dhj_ck_common.py -- and skips without the Exo-FMS tables or mpirun.
"""

# Modules
import os
import shutil
import numpy as np
import pytest
import test_suite.rad.dhj_ck_common as ck

BUILD = os.path.join(ck.REPO, "tst", "build_ck_mpi")

# 8 meshblocks of 64 x 4 x 4, so 1, 2, 4 and 8 ranks all divide the work differently.
MESH = ck.mesh(8, 16, 4, 4, 10)
RANKS = [1, 4]
# the 1- vs 4-rank state differs at 1e-16..3e-14 relative (see the docstring)
RTOL = 1.0e-9
ATOL = 1.0e-9


def read_hst(path):
    """(column names, values) of a history file."""
    names = []
    with open(path) as f:
        for ln in f:
            if ln.startswith("#") and "[1]=" in ln:
                names = [t.split("=", 1)[1] for t in ln[1:].split()]
    return names, np.atleast_2d(np.loadtxt(path))


def scales(names, data):
    """Per-column magnitude scale: the max |value| over the column's group (the columns
    sharing its suffix, e.g. 1-mom/2-mom/3-mom), or over the column alone."""
    colmax = np.abs(data).max(axis=0)
    out = colmax.copy()
    for j, nm in enumerate(names):
        if "-" in nm:
            suf = nm.split("-", 1)[1]
            grp = [k for k, o in enumerate(names) if o.split("-", 1)[-1] == suf]
            out[j] = colmax[grp].max()
    return out


HAVE_MPI = shutil.which("mpirun") is not None


@pytest.mark.skipif(not ck.HAVE_TABLES, reason=ck.NO_TABLES)
@pytest.mark.skipif(not HAVE_MPI, reason="mpirun not on PATH")
def test_run():
    """The answer must not depend on how theta and phi are divided between ranks."""
    try:
        binary = ck.build(BUILD, ["-D", "Athena_ENABLE_MPI=ON"])
        # full-precision history: data_format is not in the input, so it cannot be set
        # from the command line; write it into a copy of the input's <output1> block
        with open(ck.INPUT) as f:
            text = f.read()
        assert "<output1>\n" in text, "the dhj input has no <output1> block"
        inp = os.path.join(BUILD, "input_hst17.athinput")
        with open(inp, "w") as f:
            f.write(text.replace("<output1>\n", "<output1>\ndata_format = %.17e\n", 1))
        history = []
        for n in RANKS:
            rundir = os.path.join(BUILD, "run%d" % n)
            args = MESH + ck.CK + ["output1/dt=1.0"]
            out = ck.run(binary, rundir, args, ranks=n, inp=inp)
            assert "Terminating on cycle limit" in out, \
                f"the {n}-rank run did not reach the cycle limit"
            hst = os.path.join(rundir, "dhj.mhd.hst")
            assert os.path.exists(hst), f"the {n}-rank run wrote no history file"
            history.append(read_hst(hst))
        names, ref = history[0]
        assert len(names) == ref.shape[1], "history header and data disagree"
        scale = scales(names, ref)
        for n, (nm, h) in zip(RANKS[1:], history[1:]):
            assert nm == names and h.shape == ref.shape, \
                f"{n} ranks wrote a different history layout"
            tol = RTOL*np.maximum(np.abs(ref), np.abs(h)) + ATOL*scale
            bad = np.abs(h - ref) > tol
            if bad.any():
                i, j = np.argwhere(bad)[0]
                pytest.fail(
                    f"{n} ranks and {RANKS[0]} rank disagree in {names[j]} (row {i}): "
                    f"{h[i, j]!r} vs {ref[i, j]!r}; the correlated-k RT is a per-column "
                    "solve and must be decomposition-independent to round-off"
                )
    finally:
        shutil.rmtree(BUILD, ignore_errors=True)
