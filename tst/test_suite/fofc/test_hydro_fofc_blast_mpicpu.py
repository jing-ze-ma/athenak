"""
First-order flux correction (FOFC) across MPI boundaries and across an SMR level
boundary.

The strong 2D double-rarefaction blast of test_hydro_fofc_blast_cpu.py is cut into 16
root MeshBlocks (4 per rank on 4 ranks) and run

  - uniform    (<mesh_refinement>/refinement = none),   16 MeshBlocks,
  - refined    (<mesh_refinement>/refinement = static), 28 MeshBlocks, with the level
    boundary sitting inside the expanding vacuum so that the rarefaction heads and the
    FOFC-flagged region both cross it,

each with fofc=true on 1 and on 4 ranks.  Checks that

  - every run completes and nothing in it is a NaN,
  - the event log records a non-zero FOFC count (otherwise the test is vacuous),
  - the total mass agrees with the fofc=false control,
  - the 4-rank final field is BIT-IDENTICAL to the 1-rank one.  Nothing in the update
    depends on the rank decomposition -- the only global reduction that feeds back into
    the evolution is the dt minimum, which is exact -- so anything less than bitwise
    equality here means FOFC flagged, or corrected, a different set of faces depending
    on which MeshBlocks happened to share a rank.  (The history integrals are summed in
    rank order and so are only equal to round-off; they are compared separately.)
"""

# Modules
import os
import glob
import shutil
import pytest
import test_suite.testutils as testutils
import athena_read
import bin_convert
import numpy as np

input_file = "inputs/fofc_blast_mpi.athinput"
_refinement = ["none", "static"]
_logfile = "fofc_blast_mpi.log"
_hstfile = "fofc_blast_mpi.hydro.hst"
_binglob = "bin/fofc_blast_mpi.hydro_w.*.bin"
_nranks = 4


def cleanup_logs():
    """Remove the event log, history file and dumps left by a run."""
    for name in [_logfile, _hstfile]:
        if os.path.exists(name):
            os.remove(name)
    shutil.rmtree("bin", ignore_errors=True)


def max_fofc_count(path):
    """Largest value of the 'fofc' column of an event log (0 if never written)."""
    counts = [0]
    with open(path, "r") as fp:
        for line in fp:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) >= 8:
                counts.append(int(fields[7]))
    return max(counts)


def read_field(path):
    """Final field from a binary dump, MeshBlocks sorted by logical location so the
    arrays line up whatever the block-to-rank assignment was."""
    fd = bin_convert.read_binary(path)
    ll = fd["mb_logical"]
    order = np.lexsort((ll[:, 0], ll[:, 1], ll[:, 2], ll[:, 3]))
    return ll[order], {k: np.asarray(v)[order] for k, v in fd["mb_data"].items()}


def run_case(nranks, fofc, refinement):
    """Run one configuration on nranks ranks; return its counters and its final field."""
    cleanup_logs()
    flags = [
        "hydro/fofc=" + ("true" if fofc else "false"),
        "mesh_refinement/refinement=" + refinement,
    ]
    assert testutils.mpi_run(input_file, flags, threads=nranks), (
        f"FOFC blast failed on {nranks} ranks for fofc={fofc}, "
        f"refinement={refinement}"
    )
    dumps = sorted(glob.glob(_binglob))
    assert dumps, "no binary dump was written"
    logical, field = read_field(dumps[-1])
    out = {
        "nfofc": max_fofc_count(_logfile),
        # athena_read raises on a NaN (testutils sets check_nan_flag)
        "hst": athena_read.hst(_hstfile),
        "logical": logical,
        "field": field,
    }
    cleanup_logs()
    return out


@pytest.mark.parametrize("ref", _refinement)
def test_fofc_blast_mpicpu_matches_serial(ref):
    """FOFC must fire, stay finite, conserve mass and be decomposition independent."""
    try:
        serial = run_case(1, True, ref)
        parallel = run_case(_nranks, True, ref)
        control = run_case(_nranks, False, ref)

        for tag, case in (("1 rank", serial), (f"{_nranks} ranks", parallel)):
            if case["nfofc"] <= 0:
                pytest.fail(
                    f"FOFC never fired on {tag} for refinement={ref}: "
                    f"the test would be vacuous"
                )
            for key, col in case["hst"].items():
                assert np.all(np.isfinite(col)), (
                    f"non-finite {key} in the history of {tag}, refinement={ref}"
                )
            for key, arr in case["field"].items():
                assert np.all(np.isfinite(arr)), (
                    f"non-finite {key} in the field of {tag}, refinement={ref}"
                )

        assert serial["nfofc"] == parallel["nfofc"], (
            f"FOFC flagged {serial['nfofc']} cells on 1 rank but "
            f"{parallel['nfofc']} on {_nranks} ranks for refinement={ref}"
        )

        # mass created by the floors differs between fofc on and off, but only at the
        # level of the floor bookkeeping itself, not at the level of a lost flux
        m_on = parallel["hst"]["mass"][-1]
        m_off = control["hst"]["mass"][-1]
        dev = abs(m_on - m_off)/abs(m_off)
        assert dev < 1.0e-3, (
            f"mass differs by {dev:g} between fofc on and off for refinement={ref}"
        )

        # the two runs must describe the same mesh before their data can be compared
        assert np.array_equal(serial["logical"], parallel["logical"]), (
            f"the 1-rank and {_nranks}-rank meshes differ for refinement={ref}"
        )
        for key in ("dens", "velx", "vely", "eint"):
            ser = serial["field"][key]
            par = parallel["field"][key]
            if not np.array_equal(ser, par):
                dev = np.max(np.abs(par - ser)/(np.abs(ser) + 1.0e-30))
                pytest.fail(
                    f"{key} differs between 1 and {_nranks} ranks for "
                    f"refinement={ref}: max relative deviation {dev:g}"
                )

        # the history integrals are summed in rank order, so only to round-off
        for key in ("mass", "tot-E"):
            ser = serial["hst"][key]
            par = parallel["hst"][key]
            dev = np.max(np.abs(par - ser)/(np.abs(ser) + 1.0e-30))
            assert dev < 1.0e-12, (
                f"{key} differs by {dev:g} between 1 and {_nranks} ranks for "
                f"refinement={ref}"
            )
    finally:
        cleanup_logs()
        testutils.cleanup()
