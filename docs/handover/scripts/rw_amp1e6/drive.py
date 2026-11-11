#!/usr/bin/env python3
"""rw_amp1e6-1111: run_radwave.py with dt_retry forced off (nominal dt per level)."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_radwave as rw  # noqa: E402
_mi = rw.make_input


def make_input(c, bg, w, vec):
    text, period, dt, icfl, nx2, nx3 = _mi(c, bg, w, vec)
    assert "dt_retry = true" in text
    text = text.replace("dt_retry = true", "dt_retry = false")
    return text, period, dt, icfl, nx2, nx3


rw.make_input = make_input
if __name__ == "__main__":
    rw.main()
