"""score7.py ARM...: L1 score per deck for the given arms (score.metric / nltek exact)"""
import glob, os, sys
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H + '/ana')
import score, s2
from nlte_exact import exact
EX = {}
for dk in sorted(os.listdir(H + '/suite')):
    out = []
    for arm in sys.argv[1:]:
        ds = glob.glob('%s/suite/%s/%s_x*' % (H, dk, arm))
        if not ds: continue
        d = ds[0]; log = open(d + '/run.log').read()
        nw = log.count('CLEAN WARNING')
        if 'rc=0' not in log: out.append('%s FAIL' % arm); continue
        if dk.startswith('nltek'):
            m = s2.nlte(d)
            if m['eps'] not in EX: EX[m['eps']] = exact(m['eps'], 6000)
            t, J = EX[m['eps']]
            s = (m['tau'] > 1e-3) & (m['tau'] < 1e5)
            e = np.interp(np.log(m['tau'][s]), np.log(t), J)
            v = float(np.abs(m['E'][s] - e).sum()/e.sum())
        else:
            v = score.metric(dk, d)[1]
        out.append('%s %.4e warn %d' % (arm, v, nw))
    print(dk, ' | '.join(out))
