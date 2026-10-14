#!/usr/bin/env python3
"""chaincmp.py C L1 L2: a continuous run C vs a 2-link chain (L1, then L2 restarted from L1's FINAL rst).
hst: every data row of C (by time) must appear in L1+L2 with identical text, and the chain must have no extra or missing
row times.  bin: every dump of C (by stem and time) must exist in L1 or L2 with identical bytes of every variable."""
import sys, glob, os, re
import numpy as np
sys.path.insert(0, os.environ.get('BINCONV_DIR', os.path.dirname(os.path.abspath(__file__))))
import bin_convert
C, L1, L2 = sys.argv[1:4]; bad = 0
for h in sorted(glob.glob(os.path.join(C, '*.hst'))):
    b = os.path.basename(h)
    rows = lambda f: [l.split() for l in open(f) if l.strip() and not l.startswith('#')] if os.path.exists(f) else []
    rc = {r[0]: r for r in rows(h)}; rch = {}
    for d in (L1, L2):
        for r in rows(os.path.join(d, b)): rch[r[0]] = r
    miss = sorted(set(rc) - set(rch), key=float); extra = sorted(set(rch) - set(rc), key=float)
    diff = [t for t in rc if t in rch and rc[t] != rch[t]]
    print('hst %s: C rows %d, chain rows %d, missing %d %s, extra %d %s, differing %d' % (b, len(rc), len(rch), len(miss), miss[:4], len(extra), extra[:4], len(diff)))
    bad += len(miss) + len(extra) + len(diff)
def dumps(d):
    out = {}
    for f in glob.glob(os.path.join(d, 'bin', '*.bin')):
        stem = re.sub(r'\.\d+\.bin$', '', os.path.basename(f)); x = bin_convert.read_binary(f)
        out[(stem, round(x['time'], 6))] = x
    return out
dc = dumps(C); dch = dumps(L1); dch.update(dumps(L2))
for k in sorted(dc):
    if k not in dch: print('bin %s t=%s MISSING in chain' % k); bad += 1; continue
    for v in dc[k]['var_names']:
        a = np.concatenate([np.ravel(q) for q in dc[k]['mb_data'][v]]); b = np.concatenate([np.ravel(q) for q in dch[k]['mb_data'][v]])
        if a.tobytes() != b.tobytes(): print('bin %s t=%s %s differs' % (k[0], k[1], v)); bad += 1
ex = sorted(set(dch) - set(dc))
print('bin: C dumps %d, chain dumps %d, extra in chain %d %s' % (len(dc), len(dch), len(ex), ex[:4]))
print('CHAINCMP %s' % ('BITWISE' if bad == 0 else '%d problems' % bad))
