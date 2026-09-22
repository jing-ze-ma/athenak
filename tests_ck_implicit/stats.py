import sys
import re
import numpy as np
for f in sys.argv[1:]:
    p, r, g, nc = [], [], [], 0
    for ln in open(f + '/run.log'):
        if '### ck_implicit' not in ln:
            continue
        d = dict(re.findall(r'(\w+)=([-\deE.+]+)', ln))
        p.append(int(d['passes']))
        r.append(float(d['res']))
        g.append(abs(float(d['ckdesum'])))
        nc += ('NOT-CONVERGED' in ln)
    t = open(f + '/t.txt').read().strip()
    if not p:
        print('%-14s wall %7s  (semi-implicit, 2 sweeps/step)' % (f, t))
        continue
    print('%-14s wall %7s  calls %3d  passes mean %4.2f max %d  sweeps/step %4.2f  '
          'res max %.2e  gap max %.2e  nonconv %d'
          % (f, t, len(p), np.mean(p), max(p), 2 * np.mean(p), max(r), max(g), nc))
