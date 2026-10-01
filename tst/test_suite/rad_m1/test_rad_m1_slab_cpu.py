"""
Regression test of the implicit M1 transport: the seeded 2-D He slab (84 x 32, one
MeshBlock, implicit transport at the true c with implicit_vimp and the be levers), 12
cycles with the Eddington closure and 12 with vet_sc (full tensor), compared with stored
numbers from the last history row.

The solver runs at implicit_tol = 1e-12, implicit_lin_tol = 1e-11, so the radiation is
fixed to ~1e-11 (at the loose default 1e-8 / 1e-10 a changed iteration order moves it
more).  RTOL_CONS = 1e-9 for the radiation and the conserved totals, whose spread between
1, 2 and 4 ranks at this tolerance is < 5e-11; RTOL_DYN = 1e-7 for the kinetic energies
and V1max, which the flow amplifies (spread 2e-9 in 12 slab cycles; tests_m1/gates/
README.md, 2).  A deliberate change of the scheme must update REF (the failure message
prints the new numbers).
"""

# Modules
import os
import pytest
import test_suite.rad_m1.m1_common as m1

NLIM = 12
RTOL_CONS = 1.0e-9
RTOL_DYN = 1.0e-7
DYN = ("1-KE", "2-KE", "V1max")
COLS = {"m1slab.hydro.hst": ["time", "mass", "tot-E", "1-KE", "2-KE"],
        "m1slab.user.hst": ["F1top", "F1mid", "F1bot", "V1max", "Etot"]}
# m1-keydefault-0927 (force_reference = wb_arad now defaults to force_reference_work =
# auto -> split; implicit_opac_update and implicit_one_pass are named in the input and
# LEVERS, so only the split acts), gcc 14, CPU serial.  Change from the e89954e2
# numbers: time 9e-6, 1-KE +0.9 % (edd) / -0.8 % (vet), Etot 1.6e-5
# (docs/dev/ke_dt_0926.md, section 8)
# box-wall-truephi (09-28): the wall mass cancellation now takes rho*Phi out with the
# TRUE potential (was Phi_eff under wb_phi_eff); change from the 0927 numbers: time
# 3.7e-9, 1-KE -3.8e-5 (edd) / -6.4e-6 (vet), V1max -8.1e-4 (edd), Etot 1.2e-10.
# rss-default-1001 (58f39f11 + 57fcd2b8): <hydro>/rad_signal_speed is now ON by default
# with implicit M1, so the hydro dt uses the radiation-modified signal speed, which is
# larger than the gas sound speed in the slab: 12 cycles reach t 1.33 (edd) / 1.43
# (vet) instead of 1.94.  57fcd2b8 also fills opac before the first dt on a fresh start
# (before it, step 1 used the gas speed alone): this moved edd again (1.428 -> 1.326);
# vet was the same at both commits.  Every other column follows from the shorter time.
# Regenerated with gcc 14, CPU serial (fixes-1001).
REF = {
    "edd": {"time": 1.3263574735599821e+00, "mass": 1.2738377353803061e+19,
            "tot-E": 8.0476824191865979e+32, "1-KE": 1.2025839135391527e+24,
            "2-KE": 2.7622430561289195e+24, "F1top": 2.4752258138624450e+15,
            "F1mid": 2.4752250100079205e+15, "F1bot": 2.4751827481951940e+15,
            "V1max": 4.4470178967837883e+03, "Etot": 1.0723532628809100e+33},
    "vet": {"time": 1.4280339466865284e+00, "mass": 1.2738377353803041e+19,
            "tot-E": 8.0494317703134749e+32, "1-KE": 8.0414325226132745e+24,
            "2-KE": 5.7081441763546965e+24, "F1top": 2.4750037022798720e+15,
            "F1mid": 2.4750230810029555e+15, "F1bot": 2.4751780139539245e+15,
            "V1max": 1.0426943074024435e+04, "Etot": 1.0738808270535462e+33},
}


def last_row(rundir):
    vals = {}
    for name, cols in COLS.items():
        names, rows = m1.history(rundir, name)
        for c in cols:
            vals[c] = rows[-1][names.index(c)]
    return vals


@pytest.mark.skipif(not m1.HAVE_DATA, reason=m1.NO_DATA)
def test_run():
    """The last history row of each closure must match REF."""
    exe = m1.build(False)
    got = {}
    for key, clo in [("edd", []), ("vet", m1.VET)]:
        rundir = os.path.join(m1.BUILD[False], "slab_" + key)
        m1.run(exe, rundir, ["time/nlim=%d" % NLIM] + m1.LEVERS + m1.TIGHT + clo)
        got[key] = last_row(rundir)
    lines = []
    for k in got:
        vals = ", ".join("\"%s\": %.16e" % (c, v) for c, v in got[k].items())
        lines.append("    \"%s\": {%s}," % (k, vals))
    msg = "new numbers:\n" + "\n".join(lines)
    for key in REF:
        assert REF[key], msg
        for c, v in REF[key].items():
            rtol = RTOL_DYN if c in DYN else RTOL_CONS
            assert abs(got[key][c] - v) <= rtol * abs(v), \
                f"{key} {c}: {got[key][c]:.16e} vs REF {v:.16e}\n{msg}"
