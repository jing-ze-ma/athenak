import sys
import re
import numpy as np
for f in sys.argv[1:]:
    p, g = [], []
    for ln in open(f + '/run.log'):
        if '### ck_implicit' not in ln:
            continue
        d = dict(re.findall(r'(\w+)=([-\deE.+]+)', ln))
        p.append(int(d['passes']))
        g.append(abs(float(d['ckdesum'])))
    p, g = np.array(p), np.array(g)
    h = len(p) // 2
    print('%-12s calls %3d  passes: all %.2f  last half %.2f  gap last half %.2e'
          % (f, len(p), p.mean(), p[h:].mean(), g[h:].max()))
