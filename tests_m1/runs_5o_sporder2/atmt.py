#!/usr/bin/env python3
"""Self-convergence of the vet_col atmosphere transient (scripts/atmt.sh): L1_V of
E_n - R E_2n over L1_V(E_n), per arm.  usage: atmt.py [old new]"""
import glob
import os
import sys

import numpy as np

import solib as so

NS = [32, 64, 128, 256]
for a in sys.argv[1:] or ['old', 'new']:  # also o2only, mlin
    E = {}
    for n in NS:
        f = sorted(glob.glob(os.path.join(so.W, 'cpu', 'atmt', '%s%d' % (a, n), 'tab',
                                          '*.m1.*.tab')))[-1]
        E[n] = so.tab(f)['m1_e']
    errs = []
    for c, f in zip(NS[:-1], NS[1:]):
        rfc = so.rfaces(1.0, 5.0, c)
        w = np.diff(rfc**3)
        ef = so.restrict_r(E[f], so.rfaces(1.0, 5.0, f))
        errs.append(np.sum(w*np.abs(E[c] - ef))/np.sum(w*np.abs(E[c])))
    print('%s: L1 self-diff %s  orders %s'
          % (a, ' '.join('%.3e' % e for e in errs),
             ' '.join('%.2f' % p for p in so.order(errs))))
