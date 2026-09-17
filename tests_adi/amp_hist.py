"""Amplitude history of the initial harmonic: A(t) = sum V (T-T0) f / sum V f^2.

The volume weight is the exact cubed-sphere cell volume up to a factor common to every
cell of a radial plane (the gnomonic solid angle), which cancels in the ratio only
approximately -- good enough to see a sign change or growth, which is what the stability
arm is looking for.  Prints A(t)/A(0) per dump.
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, "/viper/u2/jinma/ATHENAK/bench/wt_he4_adi/vis/python")
import bin_convert  # noqa: E402

GAMMA = 1.666667
d = sys.argv[1]
files = sorted(glob.glob(os.path.join(d, "bin", "lap.hydro_w.*.bin")))
r0 = bin_convert.read_binary(files[0])
e0 = np.asarray(r0["mb_data"]["eint"])
rho0 = np.asarray(r0["mb_data"]["dens"])
f = e0 * (GAMMA - 1.0) / rho0 - 1.0
den = np.sum(f * f)
out = []
for fn in files:
    r = bin_convert.read_binary(fn)
    e = np.asarray(r["mb_data"]["eint"])
    rho = np.asarray(r["mb_data"]["dens"])
    t = e * (GAMMA - 1.0) / rho - 1.0
    out.append((r["time"], np.sum(t * f) / den, np.abs(t).max() / np.abs(f).max()))
for tm, a, mx in out:
    print("t=%12.6e  A/A0=%14.6e  max|dT|/max|dT0|=%14.6e" % (tm, a, mx))
mono = all(out[i][1] <= out[i - 1][1] + 1e-12 for i in range(1, len(out)))
print("monotone decreasing:", mono, " final A/A0 =", out[-1][1],
      " sign changes:", sum(1 for i in range(1, len(out))
                            if out[i][1] * out[i - 1][1] < 0))
