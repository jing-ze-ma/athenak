"""cmpruns.py A B [B2 ...]: hst files (numeric rows, '#' lines ignored) and every bin dump of A vs B (all variables,
exact). B2..: further dirs whose hst rows/dumps fill in (chain links in one dir already; extra dirs optional).
Prints BITWISE / DIFF per file and a final verdict line CMPRUNS BITWISE|DIFF."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "e1caltech-1013"))
import bin_convert as bc
A, B = sys.argv[1], sys.argv[2]
ok = True
for h in sorted(glob.glob(os.path.join(A, '*.hst'))):
    hb = os.path.join(B, os.path.basename(h))
    ra = [l for l in open(h) if not l.startswith('#')]
    rb = [l for l in open(hb) if not l.startswith('#')] if os.path.exists(hb) else []
    da = {l.split()[0]: l for l in ra}; db = {l.split()[0]: l for l in rb}
    com = [k for k in da if k in db]
    nd = sum(1 for k in com if da[k] != db[k])
    print('hst %s: rows A %d B %d common %d differing %d' % (os.path.basename(h), len(ra), len(rb), len(com), nd))
    if nd or len(com) != len(ra): ok = False
fa = sorted(glob.glob(os.path.join(A, 'bin', '*.bin')))
nb = nbd = 0
for f in fa:
    g = os.path.join(B, 'bin', os.path.basename(f))
    if not os.path.exists(g): print('missing in B', os.path.basename(f)); ok = False; continue
    x = bc.read_binary(f); y = bc.read_binary(g)
    same = (x['time'] == y['time']) and all(np.array_equal(np.asarray(x['mb_data'][v]), np.asarray(y['mb_data'][v])) for v in x['var_names'])
    nb += 1
    if not same:
        nbd += 1; ok = False
        print('bin DIFF', os.path.basename(f), x['time'], y['time'])
print('bins compared %d differing %d' % (nb, nbd))
print('CMPRUNS', 'BITWISE' if ok else 'DIFF')
