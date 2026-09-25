#!/usr/bin/env python3
"""rstcmp.py straight cont: the continuation's rst files (payload after <par_end>) and its
hst rows vs the same files / rows of the straight run (a restart gate; runs_5h_sph2)."""
import glob
import os
import sys


def payload(p):
    d = open(p, 'rb').read()
    k = d.find(b'<par_end>')
    return d[k:] if k >= 0 else d


a, b = sys.argv[1], sys.argv[2]
ok = True
n = 0
for pb in sorted(glob.glob(os.path.join(b, 'rst', '*.rst'))):
    pa = os.path.join(a, 'rst', os.path.basename(pb))
    s = os.path.exists(pa) and payload(pa) == payload(pb)
    ok &= s
    n += 1
    print(('same ' if s else 'DIFF ') + os.path.basename(pb))
for hb in sorted(glob.glob(os.path.join(b, '*.hst'))):
    ha = os.path.join(a, os.path.basename(hb))
    rb = [ln for ln in open(hb) if not ln.startswith('#')]
    ra = set(ln for ln in open(ha) if not ln.startswith('#'))
    common = [ln for ln in rb if ln in ra]
    s = len(common) == len(rb)
    ok &= s
    print(f"{'same ' if s else 'DIFF '}{os.path.basename(hb)} "
          f"rows {len(common)}/{len(rb)}")
print(f"RESTART VERDICT {'BITWISE' if ok and n > 0 else 'DIFFERENT'} rst={n}")
