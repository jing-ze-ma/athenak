"""
First-order flux correction (FOFC) regression test in 2D.

Runs a strong 2D double-rarefaction blast that opens a near-vacuum, so the trial
conserved state genuinely needs floors and the first-order fallback is exercised.
Checks that
  - the run completes with fofc=true,
  - the event log records a non-zero FOFC count at least once (otherwise the test
    would be silently testing nothing),
  - neither the history integrals nor the final slice contain a NaN.
"""

# Modules
import os
import glob
import shutil
import pytest
import test_suite.testutils as testutils
import athena_read
import numpy as np

input_file = "inputs/fofc_blast.athinput"
_shock_dir = ["1", "2"]  # exercise the x1 and the x2 FOFC sweeps
# "general" is the general-EOS interface in its gamma mode (no table needed); it takes
# the single-state LLF fallback through wder, the derived-variable array
_eos = ["ideal", "general"]
_logfile = "fofc_blast.log"
_hstfile = "fofc_blast.hydro.hst"


def cleanup_logs():
    """Remove the event log, history file and binary dumps left by a run."""
    for name in [_logfile, _hstfile] + glob.glob("fofc_blast*.tab"):
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


@pytest.mark.parametrize("ev", _eos)
@pytest.mark.parametrize("sd", _shock_dir)
def test_fofc_blast_runs_and_fires(sd, ev):
    """FOFC must fire, and the run must stay finite."""
    try:
        cleanup_logs()
        assert testutils.run(
            input_file,
            ["problem/shock_dir=" + sd, "hydro/fofc=true", "hydro/eos=" + ev],
        ), f"FOFC blast run failed for shock_dir={sd}, eos={ev}"

        nfofc = max_fofc_count(_logfile)
        if nfofc <= 0:
            pytest.fail(
                f"FOFC never fired for shock_dir={sd}, eos={ev}: "
                f"the test would be vacuous"
            )

        # athena_read raises on NaN (testutils sets check_nan_flag), and the history
        # file holds volume integrals over every cell in the domain
        hst = athena_read.hst(_hstfile)
        for key, col in hst.items():
            if key == "time":
                continue
            assert np.all(np.isfinite(col)), f"non-finite {key} in history"
        tabs = sorted(glob.glob("tab/fofc_blast.hydro_w.*.tab"))
        assert tabs, "no tab output was written"
        data = athena_read.tab(tabs[-1])
        for key in ("dens", "velx", "eint"):
            assert np.all(np.isfinite(data[key])), f"non-finite {key} in final slice"
    finally:
        cleanup_logs()
        testutils.cleanup()
