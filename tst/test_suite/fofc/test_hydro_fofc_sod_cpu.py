"""
First-order flux correction (FOFC) regression test in 1D.

Runs Sod's shocktube with <hydro>/fofc=true.  Two things are checked:
  - the run completes: before the one_d branch was added to the FOFC stencil widening
    in hydro_fluxes.cpp, a 1D run indexed j/k ghost zones that do not exist;
  - no cell in Sod's tube ever needs a floor, so FOFC flags nothing and the solution
    must be identical, to round-off, to the same run with fofc=false.
"""

# Modules
import pytest
import test_suite.testutils as testutils
import athena_read
import numpy as np

input_file = "inputs/fofc_sod.athinput"
# plm only: FOFC with a PPM/WENO reconstruction demands nghost >= 4, and the point
# here is the 1D stencil widening, not the reconstruction
_recon = ["plm"]
_flux = ["llf", "hlle", "hllc"]

# output 1 is t=0.05 (the shock, contact and rarefaction are all well inside the
# domain but the solution is still smooth on the scale of the floors)
_early = "tab/fofc_sod.hydro_w.00001.tab"
_final = "tab/fofc_sod.hydro_w.00005.tab"


def arguments(rv, fv, fofc):
    """Assemble arguments for the run command"""
    return [
        "hydro/reconstruct=" + rv,
        "hydro/rsolver=" + fv,
        "hydro/fofc=" + fofc,
        "time/integrator=" + ("rk2" if rv == "plm" else "rk3"),
    ]


@pytest.mark.parametrize("rv", _recon)
@pytest.mark.parametrize("fv", _flux)
def test_fofc_sod_matches_no_fofc(rv, fv):
    """FOFC must be a no-op on a problem that never hits a floor."""
    try:
        assert testutils.run(input_file, arguments(rv, fv, "false"))
        off_early = athena_read.tab(_early)
        off_final = athena_read.tab(_final)
        assert testutils.run(input_file, arguments(rv, fv, "true"))
        on_early = athena_read.tab(_early)
        on_final = athena_read.tab(_final)
        for key in ("dens", "velx", "eint"):
            for tag, ref, new in (
                ("t=0.05", off_early[key], on_early[key]),
                ("t=0.25", off_final[key], on_final[key]),
            ):
                dev = np.max(np.abs(new - ref) / (np.abs(ref) + 1.0e-30))
                if dev > 1.0e-12:
                    pytest.fail(
                        f"fofc=true changed {key} at {tag} for {rv}+{fv}: "
                        f"max relative deviation {dev:g}"
                    )
    finally:
        testutils.cleanup()
