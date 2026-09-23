#!/usr/bin/env python3
"""40-period radwave stability probe: max and end |Z|/A0 of the box fundamental, the exact
decay, grid noise (rms of every other rho mode at the end)/A0, NON-CONVERGED count, stage
fallbacks.  ABORTED = non-finite state.   python3 stiff_table.py RUNS_DIR"""
import glob
import json
import math
import os
import re
import sys
import numpy as np

print("%-66s %9s %9s %9s %9s %4s %4s" % ("case", "max|Z|/A0", "end|Z|/A0", "exact",
                                          "noise/A0", "NC", "fb"))
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*nper40*"))):
    tag = os.path.basename(f).replace("edd_impl_x1_", "").replace("_enthplm", "")
    tag = tag.replace("_nper40_ny4", "").replace("_vimptrue", "")
    log = open(f + "/log.txt").read()
    nc = log.count("NON-CONVERGED after")
    m = re.search(r"stage fallbacks=(\S+)", log)
    fb = int(float(m.group(1))) if m else -1
    j = f + "/results.json"
    if not os.path.exists(j):
        print("%-66s %s" % (tag, "ABORTED (non-finite state)"))
        continue
    r = json.load(open(j))
    z = np.array(r["absZ"])
    a0 = z[0]
    wr = complex(*r["refs"][r["ref_key"]])
    rho = np.array(r["rho_end"]) - 1.0
    fh = np.fft.rfft(rho) / len(rho)
    noise = math.sqrt(2 * np.sum(abs(fh[2:]) ** 2)) / (2 * a0)
    print("%-66s %9.3g %9.3g %9.3g %9.2g %4d %4d" % (
        tag, z.max() / a0, z[-1] / a0, math.exp(wr.imag * r["t_end"]), noise, nc, fb))
