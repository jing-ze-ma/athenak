"""
First-order flux correction (FOFC) for MHD on the CUBED SPHERE.

<mhd>/fofc used to be a startup FATAL on this grid.  This runs a strong radial
double-rarefaction MHD shock tube in a closed spherical shell (reflecting radial walls,
panel seams in the two angular directions) carrying a divergence-free AZIMUTHAL guide
field, with and without FOFC, and checks that

  - the run completes with fofc=true for both fallback solvers (<mhd>/fofc_rsolver =
    llf and hlle) and the event log records a non-zero FOFC count (otherwise the test
    would be silently testing nothing),
  - neither the history integrals nor the final dump contain a NaN,
  - div(B) from the FACE fields stays at round-off, i.e. the face EMFs FOFC writes --
    which on this grid have to be rotated with GnomonicEquiangleEmfX1 exactly as the
    high-order sweep's are -- are still the ones constrained transport consumes,
  - with floors that never trip, fofc=true is BITWISE identical to fofc=false: the
    fallback must be a no-op when no cell is flagged,
  - mass is conserved as well with FOFC as without,
  - FOFC REDUCES the energy-floor events the same run needs without it,
  - <mhd>/vceil, whose application on this grid is deferred to
    Coordinates::GnomonicEquiangleRaiseVelMHD, clips the velocity on the real pass and
    flags the cell for FOFC.
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

_input = "inputs/fofc_mhd_cs.athinput"
_basename = "fofc_mhd_cs"

# A pressure floor ABOVE the minimum internal energy the double rarefaction reaches, so
# that the trial state reliably needs the floor and FOFC is actually exercised.  The
# floor in the input file itself is far below it, which is the "nothing is flagged"
# case.  NOTE a raised DENSITY floor is not usable here: raising rho lowers the specific
# kinetic energy, which lifts the internal energy back above the pressure floor, and the
# two then never fire together.
_firing_floor = "mhd/pfloor=0.05"


def cleanup():
    """Remove the event log, history file and binary dumps left by a run."""
    for name in glob.glob(_basename + "*.log") + glob.glob(_basename + "*.hst"):
        if os.path.exists(name):
            os.remove(name)
    shutil.rmtree("bin", ignore_errors=True)


def read_log(path):
    """(max fofc count, summed eos_efloor, max eos_vceil) from an event log."""
    nfofc, nefloor, nvceil = 0, 0, 0
    with open(path, "r") as fp:
        for line in fp:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) >= 8:
                nefloor += int(fields[2])
                nvceil = max(nvceil, int(fields[4]))
                nfofc = max(nfofc, int(fields[7]))
    return nfofc, nefloor, nvceil


def read_dumps(var):
    """All binary dumps of one variable, oldest first, as a list of data dicts."""
    files = sorted(glob.glob("bin/{:}.{:}.*.bin".format(_basename, var)))
    assert files, "no {:} output was written".format(var)
    return [bin_convert.read_binary(f) for f in files]


def check_mass():
    """Total mass at the first and last history dump; the history must be finite."""
    hst = athena_read.hst("{:}.mhd.hst".format(_basename))
    for key, col in hst.items():
        assert np.all(np.isfinite(col)), "non-finite {:} in history".format(key)
    mass = hst["mass"]
    return mass[0], mass[-1]


@pytest.mark.parametrize("rsolver", ["llf", "hlle"])
def test_fofc_cs_fires_and_stays_finite_cpu(rsolver):
    """FOFC must fire on the cubed sphere and div(B) must stay at round-off."""
    try:
        cleanup()
        assert testutils.run(
            _input,
            ["mhd/fofc=true", "mhd/fofc_rsolver=" + rsolver, _firing_floor],
        ), "cubed-sphere MHD FOFC run failed for fofc_rsolver={:}".format(rsolver)

        nfofc, _, _ = read_log(_basename + ".log")
        assert nfofc > 0, (
            "FOFC never fired (fofc_rsolver={:}): the test would be "
            "vacuous".format(rsolver)
        )

        check_mass()  # also asserts the history integrals are finite
        final = read_dumps("mhd_w_bcc")[-1]["mb_data"]
        for key in ("dens", "velx", "eint", "bcc1", "bcc2", "bcc3"):
            assert np.all(
                np.isfinite(np.asarray(final[key]))
            ), "non-finite {:} in final dump".format(key)

        divb = read_dumps("mhd_divb")[-1]["mb_data"]["divb"]
        maxdivb = np.max(np.abs(np.asarray(divb)))
        assert maxdivb < 1.0e-10, "div(B) = {:} is not at round-off".format(maxdivb)
    finally:
        cleanup()
        testutils.cleanup()


@pytest.mark.parametrize("rsolver", ["llf", "hlle"])
def test_fofc_cs_no_flagged_cell_is_bitwise_identical_cpu(rsolver):
    """With floors that never trip, fofc=true must change nothing at all."""
    try:
        cleanup()
        assert testutils.run(_input, ["mhd/fofc=false"])
        assert read_log(_basename + ".log")[0] == 0
        reference = [d["mb_data"] for d in read_dumps("mhd_w_bcc")]
        m_first, m_last = check_mass()
        assert abs(m_last - m_first) < 1.0e-4 * abs(
            m_first
        ), "mass is not conserved in the closed shell with fofc=false"
        cleanup()

        assert testutils.run(
            _input, ["mhd/fofc=true", "mhd/fofc_rsolver=" + rsolver]
        )
        assert (
            read_log(_basename + ".log")[0] == 0
        ), "a cell was flagged; this test only means something when none is"
        trial = [d["mb_data"] for d in read_dumps("mhd_w_bcc")]
        m_first_t, m_last_t = check_mass()
        assert m_first_t == m_first and m_last_t == m_last, (
            "fofc=true changed the mass budget with no cell flagged"
        )

        assert len(reference) == len(trial)
        for ref, new in zip(reference, trial):
            for key in ref:
                assert np.array_equal(
                    np.asarray(ref[key]), np.asarray(new[key])
                ), "fofc=true changed {:} with no cell flagged".format(key)
    finally:
        cleanup()
        testutils.cleanup()


def test_fofc_cs_removes_floor_events_cpu():
    """FOFC must reduce the energy-floor events, and must not cost mass."""
    counts, mass = [], []
    try:
        for fofc in ("false", "true"):
            cleanup()
            assert testutils.run(_input, ["mhd/fofc=" + fofc, _firing_floor])
            nfofc, nefloor, _ = read_log(_basename + ".log")
            counts.append(nefloor)
            mass.append(check_mass())
            if fofc == "false":
                assert nfofc == 0, "FOFC counted cells with fofc=false"
            else:
                assert nfofc > 0, "FOFC never fired: the test would be vacuous"

        assert counts[0] > 0, (
            "the control run needed no energy floor: the rarefaction is too weak to "
            "test FOFC"
        )
        # NOT zero, unlike the hydro cubed-sphere test: the floor this problem trips is
        # a physical one (the rarefaction really is that cold), not a trial-state
        # artefact, so the first-order fallback can only reduce it.  Measured ~8%.
        assert counts[1] < counts[0], (
            "FOFC did not reduce the energy-floor events: {:} vs {:}".format(
                counts[1], counts[0]
            )
        )

        rel = [(last - first) / first for first, last in mass]
        assert max(abs(r) for r in rel) < 5.0e-2
        assert abs(rel[0] - rel[1]) < 5.0e-3, (
            "fofc changes the mass budget by " + str(rel[0] - rel[1])
        )
    finally:
        cleanup()
        testutils.cleanup()


def test_fofc_cs_vceil_clips_and_flags_cpu():
    """<mhd>/vceil must clip |v| on the real pass and flag the cell for FOFC."""
    try:
        cleanup()
        assert testutils.run(_input, ["mhd/fofc=false", "mhd/vceil=15.0"])
        nfofc, _, nvceil = read_log(_basename + ".log")
        assert nvceil > 0, "vceil never fired"
        assert nfofc == 0
        cleanup()

        assert testutils.run(_input, ["mhd/fofc=true", "mhd/vceil=15.0"])
        nfofc, _, _ = read_log(_basename + ".log")
        assert nfofc > 0, "the velocity ceiling did not flag any cell for FOFC"
    finally:
        cleanup()
        testutils.cleanup()
