# usage: python3 bitcmp.py <dirA> <dirB> <athenak>/vis/python : bitwise compare of every bin file (data arrays + time) and the hst
import sys, glob, os
import numpy as np
sys.path.insert(0, sys.argv[3] if len(sys.argv) > 3 else "vis/python")
import bin_convert as bc
a, b = sys.argv[1], sys.argv[2]
fs = sorted(glob.glob(os.path.join(a, 'bin', '*.bin')))
bad = 0
for fa in fs:
    fb = os.path.join(b, 'bin', os.path.basename(fa))
    if not os.path.exists(fb):
        print('MISSING', fb); bad += 1; continue
    da, db = bc.read_binary(fa), bc.read_binary(fb)
    ok = da['time'] == db['time'] and set(da['mb_data']) == set(db['mb_data'])
    nd = 0
    for k in da['mb_data']:
        x, y = np.asarray(da['mb_data'][k]), np.asarray(db['mb_data'].get(k))
        if x.shape != y.shape or not np.array_equal(x, y):
            nd += 1; ok = False
            print('  DIFF', os.path.basename(fa), k, np.nanmax(np.abs(x - y)) if x.shape == y.shape else 'shape')
    print('%-40s cycle %s time %r  %s' % (os.path.basename(fa), da.get('cycle'), da['time'],
                                         'IDENTICAL' if ok else 'DIFFERENT'))
    bad += (not ok)
for ha in glob.glob(os.path.join(a, '*.hst')):
    hb = os.path.join(b, os.path.basename(ha))
    la = [l for l in open(ha) if not l.startswith('#')]
    lb = [l for l in open(hb) if not l.startswith('#')] if os.path.exists(hb) else []
    print(os.path.basename(ha), 'rows', len(la), 'IDENTICAL' if la == lb else 'DIFFERENT')
    bad += (la != lb)
print('BITWISE_PASS' if bad == 0 and fs else 'BITWISE_FAIL (%d, %d files)' % (bad, len(fs)))
