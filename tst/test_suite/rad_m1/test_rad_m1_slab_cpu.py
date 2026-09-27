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
REF = {
    "edd": {"time": 1.9368550497296493e+00, "mass": 1.2738377353803078e+19,
            "tot-E": 8.0476796594703527e+32, "1-KE": 2.0731991679474028e+24,
            "2-KE": 4.0143150826176702e+24, "F1top": 2.4752110390282410e+15,
            "F1mid": 2.4752123125695505e+15, "F1bot": 2.4751767037463765e+15,
            "V1max": 6.0545363318421223e+03, "Etot": 1.0723495548074596e+33},
    "vet": {"time": 1.9366909155541050e+00, "mass": 1.2738377353803092e+19,
            "tot-E": 8.0494410158165634e+32, "1-KE": 1.3699078751047554e+25,
            "2-KE": 4.1156161717895498e+24, "F1top": 2.4751298463142445e+15,
            "F1mid": 2.4751359980914905e+15, "F1bot": 2.4751734924189520e+15,
            "V1max": 1.4451317084849488e+04, "Etot": 1.0739003430020512e+33},
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
