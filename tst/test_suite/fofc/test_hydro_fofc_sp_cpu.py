"""
First-order flux correction (FOFC) on the SPHERICAL-POLAR grid.

Runs a strong radial double-rarefaction shock tube in a closed spherical shell
(reflecting radial walls, polar/periodic angles), with and without FOFC.  Checks that
  - the run completes with fofc=true and nothing in it is a NaN,
  - the event log records a non-zero FOFC count (otherwise the test is vacuous),
  - FOFC removes the ENERGY-floor events the same run needs without it.  This is the
    curvilinear trial update doing its job: with the old Cartesian dt/dx divergence the
    trial state on this grid was not the state the real update produces, so the wrong
    cells were flagged.
  - total mass agrees between the two runs, and the closed domain's total energy is
    conserved BETTER with FOFC than without.
"""

# Modules
import os
import glob
import shutil
import pytest  # noqa: F401
import test_suite.testutils as testutils
import athena_read
import numpy as np

input_file = "inputs/fofc_sp.athinput"
_logfile = "fofc_sp.log"
_hstfile = "fofc_sp.hydro.hst"


def cleanup_logs():
    """Remove the event log, history file and binary dumps left by a run."""
    for name in [_logfile, _hstfile]:
        if os.path.exists(name):
            os.remove(name)
    shutil.rmtree("bin", ignore_errors=True)


def read_log(path):
    """Return (max fofc count, summed eos_efloor count) from an event log."""
    nfofc, nefloor = 0, 0
    with open(path, "r") as fp:
        for line in fp:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 8:
                nefloor += int(f[2])
                nfofc = max(nfofc, int(f[7]))
    return nfofc, nefloor


def run_case(fofc):
    """Run the input with fofc on/off; return (nfofc, nefloor, history)."""
    cleanup_logs()
    flag = "true" if fofc else "false"
    assert testutils.run(input_file, ["hydro/fofc=" + flag]), (
        f"spherical-polar FOFC run failed for fofc={flag}"
    )
    nfofc, nefloor = read_log(_logfile)
    hst = athena_read.hst(_hstfile)
    return nfofc, nefloor, hst


def test_fofc_sp_cpu():
    """FOFC must fire on the spherical-polar grid and must remove the floor events."""
    try:
        nfofc_on, nefloor_on, hst_on = run_case(True)
        nfofc_off, nefloor_off, hst_off = run_case(False)

        assert nfofc_on > 0, "FOFC never fired: the test would be vacuous"
        assert nfofc_off == 0, "FOFC counted cells with fofc=false"

        for key, col in hst_on.items():
            assert np.all(np.isfinite(col)), f"non-finite {key} in history (fofc=true)"

        # FOFC exists to keep the conserved-to-primitive inversion off its floors
        assert nefloor_off > 0, (
            "the control run needed no energy floor: the blast is too weak to test FOFC"
        )
        assert nefloor_on == 0, (
            f"FOFC left {nefloor_on} energy-floor events on the spherical-polar grid "
            f"(control run: {nefloor_off})"
        )

        # closed domain: mass must agree between the two runs
        m_on, m_off = hst_on["mass"][-1], hst_off["mass"][-1]
        assert abs(m_on - m_off)/m_off < 1.0e-4, (
            f"final mass differs: fofc={m_on!r} vs no-fofc={m_off!r}"
        )

        # ...and the total energy is conserved exactly unless a floor fires
        e0 = hst_on["tot-E"][0]
        drift_on = abs(hst_on["tot-E"][-1] - e0)
        drift_off = abs(hst_off["tot-E"][-1] - e0)
        assert drift_on <= drift_off, (
            f"FOFC made energy conservation worse: {drift_on!r} vs {drift_off!r}"
        )
    finally:
        cleanup_logs()
        for name in glob.glob("fofc_sp*"):
            if os.path.isfile(name):
                os.remove(name)
        testutils.cleanup()
