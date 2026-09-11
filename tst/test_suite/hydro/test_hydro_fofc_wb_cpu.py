"""
First-order flux correction (FOFC) on top of the WELL-BALANCED (WB) schemes.

A 2D polytropic atmosphere in constant gravity (pgen wb_atm) is evolved with both
well-balanced schemes in turn, and every test below is parametrised over the two:

  "static"   <hydro>/wellbalance_static (+ wellbalance_static_reconst).  The pgen hands
             the analytic profile to the hydro module as the stored background
             (u0wb/w0wb/w0facewb), so the unperturbed problem is a machine-exact fixed
             point.  During Hydro::Fluxes the perturbation reconstruction has SUBTRACTED
             that background from w0 (RemoveWbVar), so in that window w0 is NOT a full
             state.  The first-order fallback is a physical Riemann solve and must be
             handed the FULL state; hydro_fofc.cpp adds the background back for the
             flagged cells only.

  "dynamic"  <hydro>/wellbalance_dynamic + wb_x2 + wb_rho (Kappeli & Mishra), with
             wb_option = isentropic, which is exactly this atmosphere's stratification.
             There is no stored background at all: the scheme rebuilds a local
             hydrostatic one per stencil out of the current state and the gravitational
             potential the pgen fills (phicc0, phi0), and the gravitational source term
             is that background's own pressure difference across the cell.  w0 therefore
             stays a full state throughout and FOFC needs no special handling -- which is
             what these cases pin down, so that the static-only fix-up in hydro_fofc.cpp
             cannot silently start applying (or failing to apply) here.

A strong radially diverging velocity kick high up in the atmosphere then drives a wide
region through the density floor, which is what makes FOFC fire.  Checks:

  1. with fofc=true the run completes, FOFC fires, and nothing is a NaN;
  2. with the kick switched off (vamp=0) fofc=true flags nothing and is bit-identical to
     fofc=false, and the balance itself survives -- exactly for the static scheme, to
     round-off for the dynamic one (see test_fofc_wb_balanced_state_untouched);
  3. total mass with fofc=true agrees with the fofc=false control;
  4. THE KEY ONE, see test_fofc_wb_fallback_sees_full_state below.
"""

# Modules
import os
import shutil
import pytest
import test_suite.testutils as testutils
import athena_read
import bin_convert
import numpy as np

input_file = "inputs/fofc_wb.athinput"
_logfile = "fofc_wb.log"
_hstfile = "fofc_wb.hydro.hst"
_tabfile = "tab/fofc_wb.hydro_w.00001.tab"
_binfile = "bin/fofc_wb.hydro_w.00001.bin"

# the two well-balanced schemes every test below is run with
wb_schemes = ["static", "dynamic"]


def wb_flags(wb):
    """Command-line flags selecting one of the two WB schemes, or neither.

    `wb` is "static", "dynamic" or "off".  Note that the dynamic scheme is switched on
    INSTEAD of the static one: they are alternative discretisations of the same balance
    and are not meant to be combined.
    """
    static = (wb == "static")
    dynamic = (wb == "dynamic")
    return [
        "hydro/wellbalance_static=" + ("true" if static else "false"),
        "hydro/wellbalance_static_reconst=" + ("true" if static else "false"),
        "hydro/wellbalance_dynamic=" + ("true" if dynamic else "false"),
        "hydro/wb_x2=" + ("true" if dynamic else "false"),
        "hydro/wb_rho=" + ("true" if dynamic else "false"),
    ]


def cleanup_logs():
    """Remove the event log, history file and dumps left by a run."""
    for name in [_logfile, _hstfile]:
        if os.path.exists(name):
            os.remove(name)
    shutil.rmtree("bin", ignore_errors=True)
    shutil.rmtree("tab", ignore_errors=True)


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
    """Final 2D state from a binary dump, MeshBlocks sorted by logical location so the
    arrays line up between runs with different block-to-rank assignments."""
    fd = bin_convert.read_binary(path)
    ll = fd["mb_logical"]
    order = np.lexsort((ll[:, 0], ll[:, 1], ll[:, 2], ll[:, 3]))
    return {k: np.asarray(v)[order] for k, v in fd["mb_data"].items()}


def run_case(wb, fofc, vamp=None):
    """Run one configuration; return its counters and its final state."""
    cleanup_logs()
    flags = ["hydro/fofc=" + ("true" if fofc else "false")] + wb_flags(wb)
    if vamp is not None:
        flags.append("problem/vamp=" + repr(vamp))
    assert testutils.run(input_file, flags), (
        f"wb atmosphere run failed for wb={wb}, fofc={fofc}, vamp={vamp}"
    )
    # athena_read raises on a NaN (testutils sets check_nan_flag)
    out = {
        "nfofc": max_fofc_count(_logfile),
        "hst": athena_read.hst(_hstfile),
        "tab": athena_read.tab(_tabfile),
        "field": read_field(_binfile),
    }
    cleanup_logs()
    return out


def assert_finite(case, tag):
    """Every history integral, every point of the cut and of the 2D field is finite."""
    for key, col in case["hst"].items():
        assert np.all(np.isfinite(col)), f"non-finite {key} in history of {tag}"
    for key, col in case["tab"].items():
        assert np.all(np.isfinite(col)), f"non-finite {key} in the x1 cut of {tag}"
    for key, arr in case["field"].items():
        assert np.all(np.isfinite(arr)), f"non-finite {key} in the 2D field of {tag}"


def l1_rel(new, ref):
    """L1 norm of new - ref, normalised by the L1 norm of ref."""
    return np.abs(new - ref).sum()/np.abs(ref).sum()


@pytest.mark.parametrize("wb", wb_schemes)
def test_fofc_wb_fires_and_stays_finite(wb):
    """FOFC must fire on the perturbed well-balanced atmosphere, and stay finite."""
    try:
        on = run_case(wb=wb, fofc=True)
        if on["nfofc"] <= 0:
            pytest.fail(f"FOFC never fired with wb={wb}: the test is vacuous")
        assert_finite(on, f"wb={wb}, fofc=on")
    finally:
        cleanup_logs()
        testutils.cleanup()


@pytest.mark.parametrize("wb", wb_schemes)
def test_fofc_wb_balanced_state_untouched(wb):
    """With no perturbation the atmosphere is at rest and must stay there: FOFC must
    flag nothing and change nothing, bit for bit.

    How exactly the two schemes hold the balance differs, and the numbers below are the
    measured ones.  The STATIC scheme is exact by construction: the reconstructed
    deviation is identically zero and the vertical momentum stays at a hard 0.  The
    DYNAMIC scheme rebuilds the background per stencil out of pow()/log() closed forms,
    so it holds the balance only to round-off; measured on this input at t = 0.04
    (11 cycles), the volume-integrated vertical momentum reaches 3.5e-16 (against a
    total mass of 0.982) and the largest velocity anywhere in the domain is 7.2e-16,
    while the mass and the total energy are still bit-constant.  The gates are 1e-12,
    i.e. ~3000x the measured residual and still far below anything physical.
    """
    try:
        off = run_case(wb=wb, fofc=False, vamp=0.0)
        on = run_case(wb=wb, fofc=True, vamp=0.0)
        assert on["nfofc"] == 0, (
            f"FOFC flagged {on['nfofc']} cells in a balanced atmosphere (wb={wb})"
        )
        # the tab cut is written with %24.16e, so this is a bitwise comparison of the
        # final state along the line through the (absent) kick
        for key in ("dens", "velx", "vely", "eint"):
            assert np.array_equal(on["tab"][key], off["tab"][key]), (
                f"fofc=true changed {key} in the balanced atmosphere (wb={wb})"
            )
        # mass and total energy are conserved to the last bit by BOTH schemes here
        for key in ("mass", "tot-E"):
            col = on["hst"][key]
            assert np.array_equal(col, col[0]*np.ones_like(col)), (
                f"{key} drifted in the balanced atmosphere with fofc=true (wb={wb})"
            )
        # the balance itself: no vertical momentum, no motion (see the docstring)
        mom2 = np.abs(on["hst"]["2-mom"]).max()
        assert mom2 < 1.0e-12, (
            f"vertical momentum grew to {mom2:g} in the balanced atmosphere (wb={wb})"
        )
        vmax = max(np.abs(on["field"][key]).max() for key in ("velx", "vely", "velz"))
        assert vmax < 1.0e-12, (
            f"the balanced atmosphere moved at {vmax:g} (wb={wb})"
        )
    finally:
        cleanup_logs()
        testutils.cleanup()


@pytest.mark.parametrize("wb", wb_schemes)
def test_fofc_wb_conserves_mass(wb):
    """FOFC must not change the mass budget beyond the floors' own bookkeeping."""
    try:
        off = run_case(wb=wb, fofc=False)
        on = run_case(wb=wb, fofc=True)
        m_off = off["hst"]["mass"][-1]
        m_on = on["hst"]["mass"][-1]
        # both runs create mass at the density floor; FOFC changes WHICH cells are
        # floored, so the two budgets differ at the level of the floor creation itself
        # (measured 1.0e-4 of the total mass for either scheme), not at the level of a
        # lost flux
        dev = abs(m_on - m_off)/abs(m_off)
        assert dev < 1.0e-3, (
            f"mass differs by {dev:g} between fofc on and off (wb={wb})"
        )
    finally:
        cleanup_logs()
        testutils.cleanup()


@pytest.mark.parametrize("wb", wb_schemes)
def test_fofc_wb_fallback_sees_full_state(wb):
    """THE KEY TEST: the first-order fallback must be handed the FULL state.

    The same physical problem is run twice with fofc=true, once with the WB scheme on
    and once with it off.  The two schemes differ only in how the same fluxes are
    assembled (deviation reconstruction + background flux removal, or the local
    hydrostatic reconstruction of the dynamic scheme, versus the plain scheme), so their
    solutions may differ only at truncation level -- and crucially, by the SAME amount
    whether or not FOFC is on, because FOFC is supposed to be blind to the split.  If
    the fallback used the perturbation state instead, the flagged cells (which here sit
    high in the atmosphere, where rho - rho_background is a small difference of two
    comparable numbers and goes negative) would be solved with a nonsensical density;
    running exactly this test against a deliberately sabotaged hydro_fofc.cpp (wbpert_
    forced to false) gives NaN in the history integrals within 44 cycles and collapses
    the FOFC count from 4474 to 84.  The dynamic case is the guard from the other side:
    there w0 IS a full state, and the same fix-up must NOT be applied to it.

    Tolerances.  Measured on this input, WB-on/WB-off L1 differences with FOFC on and
    off, and their ratio:

                    static                        dynamic
        dens    1.6e-4 / 1.4e-4 = 1.2      1.02e-4 / 0.93e-4 = 1.09
        velx        ratio 1.8              1.15e-3 / 0.84e-3 = 1.37
        vely        ratio 1.4              1.94e-3 / 1.38e-3 = 1.40
        eint        ratio 1.3              2.61e-4 / 2.42e-4 = 1.08

    The gate is ratio < 10, i.e. FOFC is allowed to amplify the legitimate WB-on/WB-off
    difference by an order of magnitude before the test complains, which leaves a factor
    ~5 of headroom while still being far below the order-unity corruption a
    perturbation-state fallback produces.  The absolute gate (L1 < 1e-2, a factor 5 to 60
    above the measured values) catches the case where BOTH runs are corrupted in the same
    way and the ratio stays small.
    """
    try:
        wb_on_fofc = run_case(wb=wb, fofc=True)
        wb_off_fofc = run_case(wb="off", fofc=True)
        wb_on_ctrl = run_case(wb=wb, fofc=False)
        wb_off_ctrl = run_case(wb="off", fofc=False)
        for tag, case in (
            (f"wb={wb}, fofc=on", wb_on_fofc),
            ("wb=off, fofc=on", wb_off_fofc),
        ):
            assert_finite(case, tag)

        # FOFC must flag essentially the same region with and without the WB split
        n_on = wb_on_fofc["nfofc"]
        n_off = wb_off_fofc["nfofc"]
        assert n_on > 0 and n_off > 0, "FOFC did not fire in one of the two runs"
        # measured 4474 (static) and 4436 (dynamic) against 4428 with the scheme off,
        # i.e. 1% and 0.2%; a perturbation-state fallback gave 84
        rel = abs(n_on - n_off)/max(n_on, n_off)
        assert rel < 0.25, (
            f"FOFC flagged {n_on} cells with wb={wb} but {n_off} with it off"
        )

        for key in ("dens", "velx", "vely", "eint"):
            with_fofc = l1_rel(wb_on_fofc["field"][key], wb_off_fofc["field"][key])
            without = l1_rel(wb_on_ctrl["field"][key], wb_off_ctrl["field"][key])
            assert with_fofc < 1.0e-2, (
                f"{key}: WB-on and WB-off disagree by {with_fofc:g} in L1 with FOFC on"
            )
            ratio = with_fofc/max(without, 1.0e-30)
            assert ratio < 10.0, (
                f"{key}: FOFC amplified the WB-on/WB-off difference by {ratio:g} "
                f"({with_fofc:g} with FOFC, {without:g} without)"
            )
    finally:
        cleanup_logs()
        testutils.cleanup()
