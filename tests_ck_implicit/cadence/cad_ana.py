"""cad_ana.py <root>: per arm, from <arm>/run.log (rank 0 ck_cadence / ck_implicit lines):
full calls, guard calls, mean fraction of columns refreshed per linearised step, mean
edef (linear deposit minus Q0 dt, over sum |Q0| dt), mean linerr (linear model error at
the refresh), capped cells, passes per call, non-converged calls, mean |ckdesum|."""
import glob
import os
import re
import sys

import numpy as np

ROOT = sys.argv[1]
print('%-6s %6s %6s %6s %9s %9s %9s %9s %6s %5s %4s %9s' % (
    'arm', 'cycles', 'full', 'guard', 'fref', 'edef', '|edef|', 'linerr', 'ncap', 'pass',
    'nc', '|gap|'))
for d in sorted(glob.glob(os.path.join(ROOT, '*', 'run.log'))):
    arm = os.path.basename(os.path.dirname(d))
    txt = open(d).read()
    cl = re.findall(r'### ck_cadence ncycle=\d+ rank=0 mode=(\w) nref=(\d+) ncol=(\d+) '
                    r'fref=(\S+) dvmax=\S+ edef=(\S+) sde=\S+ sq=\S+ ncap=(\d+) '
                    r'linerr=(\S+)', txt)
    ps = [int(p) for p in re.findall(r'passes=(\d+) rank=0', txt)]
    gap = [abs(float(g)) for g in re.findall(r'ckdesum=(\S+)', txt)]
    nc = txt.count('NOT-CONVERGED')
    if cl:
        mode = [c[0] for c in cl]
        lin = [c for c in cl if c[0] != 'F']
        fref = np.mean([float(c[3]) for c in lin]) if lin else 0.0
        ed = [float(c[4]) for c in lin]
        le = [float(c[6]) for c in cl if float(c[6]) >= 0.0]
        ncap = sum(int(c[5]) for c in cl)
        print('%-6s %6d %6d %6d %9.3e %9.2e %9.2e %9.2e %6d %5.2f %4d %9.2e' % (
            arm, len(cl), mode.count('F'), mode.count('G'), fref,
            np.mean(ed) if ed else 0, np.mean(np.abs(ed)) if ed else 0,
            np.mean(le) if le else -1, ncap, np.mean(ps) if ps else 0, nc,
            np.mean(gap) if gap else 0))
    else:
        print('%-6s %6s %6d %6s %9s %9s %9s %9s %6s %5.2f %4d %9.2e' % (
            arm, '-', len(ps), '-', '-', '-', '-', '-', '-', np.mean(ps) if ps else 0, nc,
            np.mean(gap) if gap else 0))
