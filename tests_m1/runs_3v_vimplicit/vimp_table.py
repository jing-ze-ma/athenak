#!/usr/bin/env python3
"""Tables of the runs_3v_vimplicit radwave gates from results.json (run_radwave.py).
  python3 vimp_table.py <runs dir> <list>...
Per case: fitted omega (meaningful for nper = 1 only), L1(rho), L1(E), the fundamental
amplitude |Z| at the end over its start, the largest non-fundamental Fourier amplitude of
rho at the end (grid noise) in units of the initial amplitude, Picard and inner means."""
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_radwave as rw  # noqa: E402

D = sys.argv[1]
for lst in sys.argv[2:]:
    print("==", os.path.basename(lst))
    for ln in open(lst):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        c = rw.parse_case(ln)
        tag = rw.case_tag(c)
        p = os.path.join(D, tag, "results.json")
        if not os.path.exists(p):
            print("%-72s MISSING" % tag)
            continue
        r = json.load(open(p))
        a = np.abs(np.array(r["absZ"]))
        rho = np.array(r["rho_end"])
        f = np.fft.rfft(rho - rho.mean()) / len(rho) * 2 / rw.AMP
        f[1] = 0.0
        noise = np.abs(f).max()
        pic = inn = float("nan")
        log = open(os.path.join(D, tag, "log.txt")).read()
        m = re.search(r"Picard iterations mean=(\S+)", log)
        if m:
            pic = float(m.group(1))
        m = re.search(r"inner iterations mean=(\S+)", log)
        if m:
            inn = float(m.group(1))
        nc = re.search(r"NON-CONVERGED=(\S+)", log)
        print("%-72s w=(%9.4f,%8.4f) L1rho %.3e L1E %.3e Zend/Z0 %.3e noise %.1e "
              "Picard %.2f inner %.2f NC %s"
              % (tag, r["w_re"], r["w_im"], r["l1_rho"], r["l1_E"], a[-1] / a[0], noise,
                 pic, inn, nc.group(1) if nc else "?"))
