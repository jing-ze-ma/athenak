"""cmpdata.py <dirA> <dirB>: compare every output of two runs by DATA (bin: the mesh-block
arrays and time/cycle, ignoring the header, which carries the input file; others: bytes).
Prints '<n files> <n differ> <max rel diff>'."""
import sys
import os
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_thinstab/vis/python')
import bin_convert as bc  # noqa: E402

a, b = sys.argv[1], sys.argv[2]
n = nd = 0
mx = 0.0
for root, _, fs in os.walk(a):
    for f in fs:
        if f == 'log.txt':
            continue
        pa = os.path.join(root, f)
        pb = os.path.join(b, os.path.relpath(pa, a))
        n += 1
        if f.endswith('.bin'):
            x, y = bc.read_binary(pa), bc.read_binary(pb)
            d = 0.0
            for v in x['var_names']:
                for p, q in zip(x['mb_data'][v], y['mb_data'][v]):
                    p, q = np.asarray(p, float), np.asarray(q, float)
                    d = max(d, float(np.max(np.abs(p - q) / (np.abs(p) + 1e-300))))
            if x['time'] != y['time'] or d > 0:
                nd += 1
                mx = max(mx, d)
        else:
            if open(pa, 'rb').read() != open(pb, 'rb').read():
                nd += 1
print('%d files, %d differ, max rel diff %.2e' % (n, nd, mx))
