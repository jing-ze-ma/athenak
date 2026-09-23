# usage: hstrst.py <base> full rst... : hst rows of the restart run vs the full run, matched by time
import sys, numpy as np
D = sys.argv[1]
for c in sys.argv[2:]:
    for h in ('hydro', 'user'):
        a = np.loadtxt(f'{D}/{c}_full/m1slab.{h}.hst'); b = np.loadtxt(f'{D}/{c}_rst/m1slab.{h}.hst')
        ta = {r[0]: r for r in a}
        m = [(ta[r[0]], r) for r in b if r[0] in ta]
        eq = all((x == y).all() for x, y in m)
        print(c, h, 'rows rst', len(b), 'matched', len(m), 'BITWISE' if eq and len(m) == len(b) else 'DIFF',
              'times', b[0, 0], b[-1, 0])
