"""
Advected-magnetic-vortex regression test for the LHLLD Riemann solver.

Runs the Balsara (2004) MHD vortex in the form used by Leidi et al. (2022, Sect. 5.2) on
a 64^2 periodic box [-5,5]^2, with vortex strength Vtil = 1e-2 (maximum sonic Mach
number 1.55e-2) and beta_K = Btil^2/Vtil^2 = 1, advected once along the grid diagonal.
The vortex is an exact stationary solution and one advective crossing returns it to its
starting point, so the exact solution at t_adv is the initial condition: everything the
run loses is numerical dissipation.

Checks that
  - both hlld and lhlld complete the crossing with finite, bounded L1 errors,
  - lhlld dissipates markedly less of the vortex kinetic energy than hlld.

Measured on this branch at 64^2 (rk2/plm, CFL 0.4):

              L1(v)        Ekin/Ekin0   Emag/Emag0
    hlld      1.137e-03    0.1861       0.7247
    lhlld     1.024e-03    0.3082       0.5549

The clean separation is in the kinetic energy, where lhlld retains 1.66x what hlld does.
The L1 of the velocity separates the two much less (lhlld is only 0.90x hlld): by one
crossing at this resolution both solvers have lost most of the vortex, so L1(v)
saturates near the amplitude of the initial vortex itself and stops discriminating.
The velocity-L1 gate below is therefore set at the measured 0.90 plus margin, NOT at a
larger claimed advantage, and the kinetic-energy gate carries the real signal.

Note also that lhlld ends with LESS magnetic energy than hlld.  That is expected here
rather than a defect: beta_K = 1 puts the Alfven speed at the rotation speed, so the
magnetic floor on chi is active and, as in Leidi et al., the residual O(1/M_Alfven)
pressure diffusion of LHLLD is not switched off.  hlld's larger final magnetic energy
comes from field wound up by a velocity field that its own dissipation has already
smeared out, so it is not a sign of a better-preserved vortex.
"""

# Modules
import os
import pytest
import test_suite.testutils as testutils
import numpy as np

_input = "inputs/balsara_vortex.athinput"
_errfile = "balsara-vortex.dat"

# Columns of <basename>-vortex.dat written by the pgen.
_col_l1v = 12
_col_ekin = 14

# Upper bounds on the L1 error of the velocity, per solver, at the measured value plus
# 15%.  hlld is listed at the level it is known to reach, so a change to the vortex
# setup that made hlld better is caught here too.
_max_l1v = {"hlld": 1.308e-03, "lhlld": 1.177e-03}
# Lower bounds on the surviving fraction of the vortex kinetic energy, measured value
# less 15%.
_min_ekin = {"hlld": 0.158, "lhlld": 0.262}
# lhlld must not be worse than hlld on either metric.  Measured: an L1(v) ratio of 0.900
# (saturated, see the module docstring) and a kinetic-energy gain of 1.656.
_max_l1v_ratio = 0.95
_min_ekin_gain = 1.40


def cleanup():
    """Remove the error file left by a run."""
    if os.path.exists(_errfile):
        os.remove(_errfile)
    testutils.cleanup()


def run_vortex(rsolver):
    """Advect the vortex once with the given solver, and return (L1(v), Ekin/Ekin0)."""
    cleanup()
    assert testutils.run(
        _input, ["mhd/rsolver=" + rsolver]
    ), "Balsara vortex failed for rsolver={:}".format(rsolver)
    data = np.loadtxt(_errfile, ndmin=2)
    assert data.shape[0] == 1, "expected one row in " + _errfile
    row = data[0]
    assert np.all(np.isfinite(row)), "non-finite entry in " + _errfile
    return row[_col_l1v], row[_col_ekin]


def test_lhlld_dissipates_the_vortex_less():
    """The low-Mach fix must carry the advected vortex further than HLLD."""
    results = {}
    try:
        for rsolver in ("hlld", "lhlld"):
            l1v, ekin = run_vortex(rsolver)
            results[rsolver] = (l1v, ekin)
            if l1v > _max_l1v[rsolver]:
                pytest.fail(
                    "{:} velocity L1 is {:g}, threshold {:g}".format(
                        rsolver, l1v, _max_l1v[rsolver]
                    )
                )
            if ekin < _min_ekin[rsolver]:
                pytest.fail(
                    "{:} kept only {:g} of the vortex kinetic energy, "
                    "threshold {:g}".format(rsolver, ekin, _min_ekin[rsolver])
                )
        ratio = results["lhlld"][0] / results["hlld"][0]
        if ratio > _max_l1v_ratio:
            pytest.fail(
                "lhlld velocity L1 is {:g} and hlld {:g}: a ratio of {:g}, "
                "threshold {:g}".format(
                    results["lhlld"][0], results["hlld"][0], ratio, _max_l1v_ratio
                )
            )
        gain = results["lhlld"][1] / results["hlld"][1]
        if gain < _min_ekin_gain:
            pytest.fail(
                "lhlld kept {:g} and hlld {:g} of the vortex kinetic energy: "
                "a gain of {:g}, threshold {:g}".format(
                    results["lhlld"][1], results["hlld"][1], gain, _min_ekin_gain
                )
            )
    finally:
        cleanup()
