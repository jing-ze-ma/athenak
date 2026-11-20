# licmp.py RUN_A RUN_B : layout-invariance compare of two runs (hst columns, final bins per output, MeshBlocks matched by logical location)
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import bin_convert as bc
A, B = [os.path.join(x, 'a') for x in sys.argv[1:3]]
def rows(f): return np.array([[float(v) for v in l.split()] for l in open(f) if not l.startswith('#')])
for h in sorted(glob.glob(A + '/*.hst')):
    g = os.path.join(B, os.path.basename(h))
    if not os.path.exists(g): print('hst', os.path.basename(h), 'missing in B'); continue
    a, b = rows(h), rows(g); n = min(len(a), len(b)); a, b = a[:n], b[:n]
    d = a != b; cols = sorted(set(np.where(d)[1] + 1)); r0 = np.where(d.any(1))[0]
    rel = np.abs(a - b) / (np.abs(a) + np.abs(b) + 1e-300)
    print('hst', os.path.basename(h), 'rows', n, 'differing cols (1-based)', cols, 'first diff row', (r0[0] if len(r0) else None), 'max rel %.3e' % rel.max())
last = {}
for f in sorted(glob.glob(A + '/bin/*.bin')): last[f.split('.')[-3]] = f
for k, f in sorted(last.items()):
    g = os.path.join(B, 'bin', os.path.basename(f))
    if not os.path.exists(g): print('bin', os.path.basename(f), 'missing in B'); continue
    da, db = bc.read_binary(f), bc.read_binary(g)
    key = lambda d, m: tuple(np.asarray(d['mb_logical'][m]).tolist()) + tuple(np.round(np.asarray(d['mb_geometry'][m]), 12).tolist())
    ia = sorted(range(da['n_mbs']), key=lambda m: key(da, m)); ib = sorted(range(db['n_mbs']), key=lambda m: key(db, m))
    if [key(da, m) for m in ia] != [key(db, m) for m in ib]: print('bin', os.path.basename(f), 'MeshBlock sets differ'); continue
    for v in da['var_names']:
        x = np.stack([np.asarray(da['mb_data'][v][m]) for m in ia]); y = np.stack([np.asarray(db['mb_data'][v][m]) for m in ib])
        nd = int((x != y).sum()); rel = (np.abs(x - y) / (np.abs(x) + np.abs(y) + 1e-300)).max()
        print('bin', os.path.basename(f), v, 'identical' if nd == 0 else 'DIFF cells %d max rel %.3e' % (nd, rel))
