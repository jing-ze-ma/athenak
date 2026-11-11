#!/usr/bin/env python3
"""rw_p333-1112: rw_amp1e6 driver (dt_retry off) + tab dumps of m1_vet (D tensor) and m1_vsc
(formal-solution J, F/c) at the hydro cadence, for the closure check of the P 3.33 tau 10 cell."""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "rw_amp1e6"))
import run_radwave as rw  # noqa: E402
_mi = rw.make_input


def make_input(c, bg, w, vec):
    text, period, dt, icfl, nx2, nx3 = _mi(c, bg, w, vec)
    assert "dt_retry = true" in text
    text = text.replace("dt_retry = true", "dt_retry = false")
    sl = "slice_x2 = %.17g\n" % (0.5 / nx2) if nx2 > 1 else ""
    for no, var in ((3, "m1_vet"), (4, "m1_vsc")):
        text += ("\n<output%d>\nfile_type = tab\nvariable = %s\ndata_format = %%.16e\n"
                 "dt = %.17g\n%s" % (no, var, period / rw.NDUMP, sl))
    return text, period, dt, icfl, nx2, nx3


rw.make_input = make_input
if __name__ == "__main__":
    rw.main()
