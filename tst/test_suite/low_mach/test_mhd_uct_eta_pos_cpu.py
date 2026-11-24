"""
Regression test of <mhd>/lhlld_emf_eta_pos (LHLLD Step 8, default on) under
<mhd>/emf_scheme = uct_hlld (uctaudit-1103).

The UCT-HLLD corner EMF is built from per-face coefficients (UctFaceStore, MD21), not
from the face EMFs, so a face-flux modifier acts on CT only if it is mirrored in those
coefficients; eta_pos is (ceta/2 added to d^L and d^R).  Two cases:

1. inputs/uct_eta_pos.athinput: Mach 0.1 colliding flows across a B_y jump, 2-D but
   aligned with x1.  There the UCT corner EMF reduces to the x1-face flux, so
   - uct_hlld reproduces gs07 (which reads the face EMFs) to round-off, eta_pos on and off;
   - eta_pos changes the uct_hlld field (it acts);
   - it removes magnetic energy (ME_on <= ME_off at every output: no ME is created by the
     negative numerical resistivity at the converging contact faces), the total energy
     is the same in both runs, and the removed ME reappears as internal energy.
2. inputs/uct_eta_pos_2d.athinput: a 2dx velocity checkerboard (-1)^(i+j) in the
   Gresho ring field (vortex off), periodic: genuinely 2-D corner averaging.  eta_pos
   on vs off under uct_hlld: ME_on < ME_off at every output, the total energy is
   conserved in both runs, and the ME difference is heat.

Measured (uctaudit-1103, CPU): case 1 uct_hlld - gs07 max|dB_y| 6e-16, on - off max|dB_y|
2.1e-3, ME_off - ME_on 7.8e-8 (9e-5 of ME), dE_int +7.6e-8, dKE +1.9e-9, dE_tot 3e-17;
case 2 ME_off - ME_on 8.0e-7, dE_int +8.0e-7, |dE_tot/E| <= 3e-15.
"""

# Modules
import os
import glob
import shutil
import pytest
import numpy as np
import test_suite.testutils as testutils

_in1 = "inputs/uct_eta_pos.athinput"
_in2 = "inputs/uct_eta_pos_2d.athinput"
_cols = ["dens", "velx", "vely", "velz", "eint", "bcc1", "bcc2", "bcc3"]


def cleanup():
    """Remove the files left by a run."""
    for name in glob.glob("uct_eta_pos*.hst"):
        os.remove(name)
    if os.path.isdir("tab"):
        shutil.rmtree("tab")
    testutils.cleanup()


def hst(basename):
    """(E_tot, KE, ME, E_int) time series of a history file."""
    h = np.loadtxt("{:}.mhd.hst".format(basename))
    assert np.all(np.isfinite(h)), "non-finite history"
    e = h[:, 6]
    ke = h[:, 7:10].sum(axis=1)
    me = h[:, 10:13].sum(axis=1)
    return e, ke, me, e - ke - me


def run(inp, basename, scheme, eta_pos, tab=False):
    """One run: its history and (case 1) the final x1 profile."""
    cleanup()
    assert testutils.run(
        inp,
        ["mhd/emf_scheme=" + scheme, "mhd/lhlld_emf_eta_pos=" + eta_pos],
    ), "{:} failed for emf_scheme={:} eta_pos={:}".format(inp, scheme, eta_pos)
    out = {"hst": hst(basename)}
    if tab:
        a = np.loadtxt(sorted(glob.glob("tab/*.tab"))[-1], comments="#")
        out.update({k: a[:, 3 + n] for n, k in enumerate(_cols)})
    return out


def check_budget(on, off, what):
    """eta_pos removes ME at every output and the removed ME is heat."""
    e1, k1, m1, i1 = on["hst"]
    e0, k0, m0, i0 = off["hst"]
    scale = np.max(np.abs(e0))
    assert np.max(np.abs(e1 - e0)) < 1.0e-12 * scale, (
        "{:}: total energy differs between eta_pos on and off".format(what))
    assert np.all(m1 - m0 <= 1.0e-14 * np.max(m0)), (
        "{:}: eta_pos created magnetic energy at some output".format(what))
    dme = m1[-1] - m0[-1]
    assert dme < -1.0e-6 * m0[-1], (
        "{:}: eta_pos removed only {:g} of the magnetic energy".format(what, dme))
    assert (i1[-1] - i0[-1]) > -0.9 * dme, (
        "{:}: the removed magnetic energy {:g} is not heat (dE_int {:g})".format(
            what, -dme, i1[-1] - i0[-1]))


def test_uct_eta_pos_colliding():
    """Aligned colliding flows: UCT-HLLD = face flux with eta_pos on and off."""
    try:
        r = {(s, e): run(_in1, "uct_eta_pos", s, e, tab=True)
             for s in ("uct_hlld", "gs07") for e in ("true", "false")}
        for e in ("true", "false"):
            for k in ("bcc2", "dens"):
                d = np.max(np.abs(r[("uct_hlld", e)][k] - r[("gs07", e)][k]))
                assert d < 1.0e-12, (
                    "eta_pos={:}: uct_hlld and gs07 differ by {:g} in {:}".format(e, d, k))
        d = np.max(np.abs(r[("uct_hlld", "true")]["bcc2"] -
                          r[("uct_hlld", "false")]["bcc2"]))
        assert d > 1.0e-4, "eta_pos changes the uct_hlld B_y by only {:g}".format(d)
        check_budget(r[("uct_hlld", "true")], r[("uct_hlld", "false")], "colliding")
    finally:
        cleanup()


def test_uct_eta_pos_2d():
    """Genuinely 2-D converging checkerboard: eta_pos acts through the UCT corner."""
    try:
        on = run(_in2, "uct_eta_pos_2d", "uct_hlld", "true")
        off = run(_in2, "uct_eta_pos_2d", "uct_hlld", "false")
        for x in (on, off):
            e = x["hst"][0]
            assert abs(e[-1] - e[0]) < 1.0e-12 * abs(e[0]), "total energy not conserved"
        check_budget(on, off, "2-D")
    finally:
        cleanup()
