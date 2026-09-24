#!/usr/bin/env python3
"""T-S3 Marshak wave: L1(E) of the Cartesian and the thin-shell (r in [R, R+1]) runs vs
the S_N reference of tests_m1/runs_3a/t6_ref_sn.txt (t6_marshak.py --marshak), and the
shell-vs-Cartesian difference.  usage: ts3.py cart.tab sph.tab R"""
import sys

import numpy as np

import s1lib
from s1lib import common

ref = np.loadtxt('/viper/ptmp2/jinma/wt_m1sp/tests_m1/runs_3a/t6_ref_sn.txt')
c = s1lib.load_tab(sys.argv[1])['m1_e']
s = s1lib.load_tab(sys.argv[2])['m1_e']
R = float(sys.argv[3])
n = c.size
xc = (np.arange(n) + 0.5)/n
_, x1v = s1lib.rgrid(R, R + 1.0, n)
xs = x1v - R
lc = common.l1_rel(c, np.interp(xc, ref[:, 0], ref[:, 1]))
ls = common.l1_rel(s, np.interp(xs, ref[:, 0], ref[:, 1]))
print(f"L1(E) vs S_N: Cartesian {lc:.5f}  shell R={R:g} {ls:.5f}  "
      f"shell-vs-Cartesian L1 {common.l1_rel(s, c):.3e}  max rel "
      f"{np.max(np.abs(s - c)/c):.3e}")
