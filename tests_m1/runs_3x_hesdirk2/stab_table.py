#!/usr/bin/env python3
"""Boundedness of the nper=10 radwave runs: |Z(t)|/A0 of the box fundamental (max over
the 16 dumps and at the end), the exact decay exp(Im w_ref t_end), and the grid noise =
rms of every non-fundamental Fourier mode of rho at the end over A0.
   python3 stab_table.py RUNS_DIR"""
import glob
import json
import math
import os
import sys
import numpy as np

AMP = 1.0e-5
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*nper10*", "results.json"))):
    r = json.load(open(f))
    tag = os.path.basename(os.path.dirname(f))
    z = np.array(r["absZ"])
    a0 = z[0]
    wr = complex(*r["refs"][r["ref_key"]])
    rho = np.array(r["rho_end"]) - 1.0
    fh = np.fft.rfft(rho) / len(rho)
    noise = math.sqrt(2.0 * np.sum(np.abs(fh[2:]) ** 2)) / (2.0 * a0)
    log = open(os.path.join(os.path.dirname(f), "log.txt")).read()
    nc = log.count("NON-CONVERGED after")
    fb = [s for s in log.splitlines() if "stage fallbacks" in s]
    print("%-78s max|Z|/A0 %.3f end %.3e exact %.3e noise %.1e ncyc %s NC %d %s" % (
        tag, z.max() / a0, z[-1] / a0, math.exp(wr.imag * r["t_end"]), noise,
        r.get("ncycle"), nc, fb[0].split(":")[1].strip() if fb else ""))
