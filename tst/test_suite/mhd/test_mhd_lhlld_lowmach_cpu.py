"""
Low-Mach-number regression test for the LHLLD Riemann solver.

Runs the magnetized Gresho vortex (pgen_name = gresho_mhd) at Mach 0.01 with a uniform
x1-field of plasma beta 100.  A uniform field exerts no force, so the vortex is a
stationary solution and every bit of kinetic energy lost over one turnover is either
numerical dissipation or work done against the field.

Checks that
  - both hlld and lhlld run the vortex to one turnover with a finite history,
  - lhlld retains measurably more kinetic energy than hlld over that turnover at 64^2,
    which is the point of the low-dissipation fix of Minoshima & Miyoshi (2021).

Most of the kinetic energy the vortex loses here is wound into the field rather than
dissipated (at Mach 0.01 and beta 100 the Alfven speed is ten times the peak rotation
speed), so the gain from the fix is a factor of order 1.1, not an order of magnitude;
what the thresholds below guard is that the fix is present and pointing the right way.
"""

# Modules
import os
import glob
import pytest
import test_suite.testutils as testutils
import athena_read
import numpy as np

_input = "inputs/gresho_mhd.athinput"
_basename = "gresho"

# Kinetic energy fraction that must survive one turnover at 64^2, per solver.  hlld is
# listed at the level it is known to reach, so a change of the vortex setup that made
# hlld better is caught here too.  Measured: hlld 0.274, lhlld 0.302.
_min_retention = {"hlld": 0.20, "lhlld": 0.22}
# lhlld must beat hlld by at least this factor (measured 1.10).  The margin is thin
# because the chi magnetic floor deliberately gives up the extra low-Mach dissipation
# reduction once the fast speed, not the sound speed, sets the numerical viscosity --
# which it does here, where the Alfven speed is ten times the peak rotation speed.
_min_gain = 1.05


def cleanup():
    """Remove the history files left by a run."""
    for name in glob.glob(_basename + "*.hst"):
        if os.path.exists(name):
            os.remove(name)
    testutils.cleanup()


def ke_retention():
    """Fraction of the initial kinetic energy left at the end of the run."""
    hst = athena_read.hst("{:}.mhd.hst".format(_basename))
    for key, col in hst.items():
        assert np.all(np.isfinite(col)), "non-finite {:} in history".format(key)
    ke = np.asarray(hst["1-KE"]) + np.asarray(hst["2-KE"])
    me = np.asarray(hst["1-ME"]) + np.asarray(hst["2-ME"])
    assert ke[0] > 0.0
    return ke[-1] / ke[0], me[-1] - me[0]


def run_vortex(rsolver):
    """Run one turnover of the vortex with the given solver."""
    cleanup()
    assert testutils.run(
        _input, ["mhd/rsolver=" + rsolver]
    ), "Gresho vortex failed for rsolver={:}".format(rsolver)
    return ke_retention()


def test_lhlld_retains_more_kinetic_energy():
    """The low-Mach fix must keep the vortex alive markedly longer than HLLD."""
    results = {}
    try:
        for rsolver in ("hlld", "lhlld"):
            results[rsolver] = run_vortex(rsolver)
            retention = results[rsolver][0]
            if retention < _min_retention[rsolver]:
                pytest.fail(
                    "{:} kept only {:g} of the initial kinetic energy, "
                    "threshold {:g}".format(
                        rsolver, retention, _min_retention[rsolver]
                    )
                )
        gain = results["lhlld"][0] / results["hlld"][0]
        if gain < _min_gain:
            pytest.fail(
                "lhlld kept {:g} and hlld {:g} of the initial kinetic energy: "
                "a gain of {:g}, threshold {:g}".format(
                    results["lhlld"][0], results["hlld"][0], gain, _min_gain
                )
            )
        # the field is wound up either way; this only asserts the budget is sane
        for rsolver, (_, dme) in results.items():
            assert dme > 0.0, "{:} lost magnetic energy in a wound-up vortex".format(
                rsolver
            )
    finally:
        cleanup()
