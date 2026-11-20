#!/usr/bin/env python3
# picktop.py SQLITE C0 NTOP SUBCYC: top NTOP kernel keys (kname.py keys) of the window, one line each:
#   key <TAB> launch-skip <TAB> mangled name of its dominant kernel (ncu exact match with --kernel-name-base mangled;
#   skip = launches of that kernel before the window start, so ncu profiles its first launch inside the window)
import sys, os, collections; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import open_db, key, window
c, S = open_db(sys.argv[1]); c0 = int(sys.argv[2]); NT = int(sys.argv[3]); sub = float(sys.argv[4])
w0, w1, ncy = window(c, S, c0, 99, sub)
kn = c.execute("select start,end,demangledName,mangledName from CUPTI_ACTIVITY_KIND_KERNEL order by start").fetchall()
tot = collections.Counter(); names = collections.defaultdict(collections.Counter)
for st, en, dn, mn in kn:
    if w0 <= st < w1:
        k, _ = key(S.get(dn, '?')); tot[k] += en - st; names[k][mn] += en - st
print(f'# window cycles {ncy:g}; kernel sum {sum(tot.values())/ncy/1e6:.1f} ms/cyc')
for k, v in tot.most_common(NT):
    mn = names[k].most_common(1)[0][0]; skip = sum(1 for st, en, d, m in kn if st < w0 and m == mn)
    print(f'{k} {v/ncy/1e6:.1f}ms\t{skip}\t{S.get(mn, "?")}')
