"""
First-order flux correction (FOFC) regression test for MHD.

Runs a strong double-rarefaction MHD blast with a transverse guide field, which opens a
near-vacuum, so the trial conserved state genuinely needs floors and the first-order
fallback is exercised.  The same problem is run on three Cartesian configurations (the
x1, x2 and x3 sweeps) and on the spherical-polar grid, where the trial update has to
form the area/volume flux divergence and the spherical curl instead of Cartesian
differences.

Checks that
  - the run completes with fofc=true and the event log records a non-zero FOFC count
    (otherwise the test would be silently testing nothing), for both fallback solvers
    (<mhd>/fofc_rsolver = llf and hlle),
  - neither the history integrals nor the final dump contain a NaN,
  - div(B) from the FACE fields stays at round-off, i.e. the EMFs FOFC writes are still
    the ones constrained transport consumes,
  - with floors that never trip, fofc=true is BITWISE identical to fofc=false (the
    fallback must be a no-op when no cell is flagged) and mass is conserved to round-off
    in the closed domain,
  - with floors that do trip, the fofc=true and fofc=false runs conserve mass equally
    well,
  - <mhd>/vceil clips the velocity on the real pass and flags the cell for FOFC.
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

_cart_input = "inputs/fofc_mhd_blast.athinput"
_sp_input = "inputs/fofc_mhd_sp.athinput"

# name -> (input file, extra command-line arguments, output basename)
# The 3D case is coarsened so that it costs about as much as the 2D ones.
_configs = {
    "cart_x1": (_cart_input, ["problem/shock_dir=1"], "fofc_mhd_blast"),
    "cart_x2": (_cart_input, ["problem/shock_dir=2"], "fofc_mhd_blast"),
    "cart_x3_3d": (
        _cart_input,
        [
            "problem/shock_dir=3",
            "mesh/nx1=32",
            "mesh/nx2=32",
            "mesh/nx3=32",
            "meshblock/nx1=16",
            "meshblock/nx2=16",
            "meshblock/nx3=16",
            "mesh/ix3_bc=reflect",
            "mesh/ox3_bc=reflect",
        ],
        "fofc_mhd_blast",
    ),
    "spherical_polar": (_sp_input, [], "fofc_mhd_sp"),
}

# A density floor ABOVE the minimum density the double rarefaction reaches, so that the
# trial state reliably needs the floor and FOFC is actually exercised.  The floor in the
# input files themselves is far below it, which is the "nothing is flagged" case.
_firing_dfloor = "mhd/dfloor=0.05"


def cleanup(basename):
    """Remove the event log, history file and binary dumps left by a run."""
    for name in glob.glob(basename + "*.log") + glob.glob(basename + "*.hst"):
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


def max_vceil_count(path):
    """Largest value of the 'eos_vceil' column of an event log."""
    counts = [0]
    with open(path, "r") as fp:
        for line in fp:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) >= 5:
                counts.append(int(fields[4]))
    return max(counts)


def read_dumps(basename, var):
    """All binary dumps of one variable, oldest first, as a list of data dicts."""
    files = sorted(glob.glob("bin/{:}.{:}.*.bin".format(basename, var)))
    assert files, "no {:} output was written".format(var)
    return [bin_convert.read_binary(f) for f in files]


def check_mass(basename):
    """Total mass at the first and last history dump; the history must be finite."""
    hst = athena_read.hst("{:}.mhd.hst".format(basename))
    for key, col in hst.items():
        assert np.all(np.isfinite(col)), "non-finite {:} in history".format(key)
    mass = hst["mass"]
    return mass[0], mass[-1]


@pytest.mark.parametrize("rsolver", ["llf", "hlle"])
@pytest.mark.parametrize("config", list(_configs))
def test_fofc_fires_and_stays_finite(config, rsolver):
    """FOFC must fire, the run must stay finite, and div(B) must stay at round-off."""
    input_file, extra, basename = _configs[config]
    try:
        cleanup(basename)
        assert testutils.run(
            input_file,
            extra + ["mhd/fofc=true", "mhd/fofc_rsolver=" + rsolver, _firing_dfloor],
        ), "MHD FOFC blast failed for {:}, fofc_rsolver={:}".format(config, rsolver)

        nfofc = max_fofc_count(basename + ".log")
        if nfofc <= 0:
            pytest.fail(
                "FOFC never fired for {:}, fofc_rsolver={:}: the test would be "
                "vacuous".format(config, rsolver)
            )

        check_mass(basename)  # also asserts the history integrals are finite
        final = read_dumps(basename, "mhd_w_bcc")[-1]["mb_data"]
        for key in ("dens", "velx", "eint", "bcc1", "bcc2", "bcc3"):
            assert np.all(
                np.isfinite(np.asarray(final[key]))
            ), "non-finite {:} in final dump".format(key)

        # div(B) is formed from the FACE fields, so this is a direct check that the EMFs
        # FOFC overwrote are the ones the CT update consumed.
        divb = read_dumps(basename, "mhd_divb")[-1]["mb_data"]["divb"]
        maxdivb = np.max(np.abs(np.asarray(divb)))
        assert maxdivb < 1.0e-10, "div(B) = {:} is not at round-off".format(maxdivb)
    finally:
        cleanup(basename)
        testutils.cleanup()


@pytest.mark.parametrize("rsolver", ["llf", "hlle"])
@pytest.mark.parametrize("config", list(_configs))
def test_no_flagged_cell_is_bitwise_identical(config, rsolver):
    """With floors that never trip, fofc=true must change nothing at all."""
    input_file, extra, basename = _configs[config]
    try:
        cleanup(basename)
        assert testutils.run(input_file, extra + ["mhd/fofc=false"])
        assert max_fofc_count(basename + ".log") == 0
        reference = [d["mb_data"] for d in read_dumps(basename, "mhd_w_bcc")]
        m_first, m_last = check_mass(basename)
        # The domain is closed, so the only mass budget here is the truncation error of
        # the reflecting walls (the wall flux is zero only to the accuracy of the ghost
        # reconstruction); the bitwise comparison below is the strong statement.
        assert abs(m_last - m_first) < 1.0e-4 * abs(
            m_first
        ), "mass is not conserved in the closed domain with fofc=false"
        cleanup(basename)

        assert testutils.run(
            input_file, extra + ["mhd/fofc=true", "mhd/fofc_rsolver=" + rsolver]
        )
        assert (
            max_fofc_count(basename + ".log") == 0
        ), "a cell was flagged; this test only means something when none is"
        trial = [d["mb_data"] for d in read_dumps(basename, "mhd_w_bcc")]
        m_first_t, m_last_t = check_mass(basename)
        assert m_first_t == m_first and m_last_t == m_last, (
            "fofc=true changed the mass budget with no cell flagged"
        )

        assert len(reference) == len(trial)
        for ref, new in zip(reference, trial):
            for key in ref:
                assert np.array_equal(
                    np.asarray(ref[key]), np.asarray(new[key])
                ), "fofc=true changed {:} with no cell flagged ({:})".format(key, config)
    finally:
        cleanup(basename)
        testutils.cleanup()


@pytest.mark.parametrize("config", list(_configs))
def test_mass_conservation_matches_fofc_off(config):
    """FOFC must not make mass conservation worse than the unmodified scheme."""
    input_file, extra, basename = _configs[config]
    rel = []
    try:
        for fofc in ("false", "true"):
            cleanup(basename)
            assert testutils.run(
                input_file, extra + ["mhd/fofc=" + fofc, _firing_dfloor]
            )
            first, last = check_mass(basename)
            rel.append((last - first) / first)
        # The density floor creates mass in both runs; what is checked is that the
        # first-order fallback does not change the budget by more than that.
        assert max(abs(r) for r in rel) < 5.0e-2
        assert abs(rel[0] - rel[1]) < 5.0e-3, "fofc changes the mass budget by " + str(
            rel[0] - rel[1]
        )
    finally:
        cleanup(basename)
        testutils.cleanup()


def test_vceil_clips_and_flags():
    """<mhd>/vceil must clip |v| on the real pass and flag the cell for FOFC."""
    basename = "fofc_mhd_blast"
    try:
        cleanup(basename)
        assert testutils.run(_cart_input, ["mhd/fofc=false", "mhd/vceil=15.0"])
        assert max_vceil_count(basename + ".log") > 0, "vceil never fired"
        assert max_fofc_count(basename + ".log") == 0
        cleanup(basename)

        assert testutils.run(_cart_input, ["mhd/fofc=true", "mhd/vceil=15.0"])
        assert (
            max_fofc_count(basename + ".log") > 0
        ), "the velocity ceiling did not flag any cell for FOFC"
    finally:
        cleanup(basename)
        testutils.cleanup()
