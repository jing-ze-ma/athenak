"""score_all.py: L1 errors and orders of arm CC1 for every deck (nltek against exact transport)."""
import glob, os, sys
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H + '/ana')
import score
import s2
from nlte_exact import exact
EX = {}
for dk in sorted(os.listdir(H + '/suite')):
    runs = []
    for d in glob.glob('%s/suite/%s/CC1_x*' % (H, dk)):
        log = open(d + '/run.log').read()
        if 'rc=0' not in log:
            continue
        f = float(d.split('_x')[-1])
        if dk.startswith('nltek'):
            m = s2.nlte(d)
            if m['eps'] not in EX:
                EX[m['eps']] = exact(m['eps'], 6000)
            t, J = EX[m['eps']]
            s = (m['tau'] > 1e-3) & (m['tau'] < 1e5)
            e = np.interp(np.log(m['tau'][s]), np.log(t), J)
            v = float(np.abs(m['E'][s] - e).sum()/e.sum())
        else:
            v = score.metric(dk, d)[1]
        runs.append((f, v))
    runs.sort()
    errs = '  '.join('x%g %.3e' % r for r in runs)
    ords = '  '.join('%.2f' % np.log2(runs[n][1]/runs[n+1][1]) for n in range(len(runs) - 1)
                     if runs[n+1][0] == 2*runs[n][0])
    print('%-16s %s | order %s' % (dk, errs, ords))
