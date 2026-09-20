#!/usr/bin/env python3
"""tests_r10: last-8-face F_2s/F_req and emergent-luminosity ratio, mode 3
(rt_col3_skip_sweep=true), from a case's t=0 dump and its 200-cycle-relaxed dump.

Columns of mltfaces*.txt (red_giant, "i r T_f p_f grad grad_ad x F_mlt F_rad F_used
F_req F_raddiff w_blend F_rad_col0 F_conv_res grad_rad D F_2s"):
  i=0 r=1 T=2 Freq=10 Fdiff=11 w=12 F2s=17
L_out/L at a face = F_2s(face) * 4 pi r(face)^2 / L, with L the case's own lstar.
"""
import sys
import numpy as np


def load_case(path, rin, lstar, n=8):
    d = np.loadtxt(path)
    r, freq, f2s, w = d[:, 1], d[:, 10], d[:, 17], d[:, 12]
    meas = f2s/freq
    lout = 4.0*np.pi*r**2*f2s/lstar
    rows = list(zip(r[-n:]/rin, r[-n:], freq[-n:], meas[-n:], w[-n:], lout[-n:]))
    return rows


if __name__ == "__main__":
    # args: tag rin lstar t0file relaxfile
    tag, rin, lstar, f0, f1 = sys.argv[1:6]
    rin, lstar = float(rin), float(lstar)
    print("### %s  (rin=%.3e lstar=%.4e)" % (tag, rin, lstar))
    print("  --- t=0 (mode 0 default sweep, existing dump) ---" if "t0" in f0 else "")
    print("   r/r_in     r[cm]      F_req      F_2s/F_req    w      L_out/L")
    for row in load_case(f0, rin, lstar):
        print("  %7.4f  %.4e  %.4e   %8.4f   %6.3f   %8.4f" % row)
    print("  --- relaxed (mode 3, 200 cycles) ---")
    for row in load_case(f1, rin, lstar):
        print("  %7.4f  %.4e  %.4e   %8.4f   %6.3f   %8.4f" % row)
