#!/usr/bin/env python3
"""npass.py LOG...: passes per (cycle, t2s, pos) from 'CLEAN pass' lines; mean/max/sum overall and per t2s; CW count."""
import sys, re, collections
for f in sys.argv[1:]:
    st = collections.OrderedDict(); cw = 0
    for l in open(f, errors='ignore'):
        if 'CLEAN WARNING' in l: cw += 1
        m = re.search(r'CLEAN pass cycle (\d+) t2s (\d+) pos (\d+) pass (\d+)', l)
        if m: st[(m.group(1), m.group(2), m.group(3))] = int(m.group(4)) + 1
    v = list(st.values())
    if not v: print(f, 'no passes'); continue
    by = collections.defaultdict(list)
    for (c, t, p), n in st.items(): by[t].append(n)
    print('%s: stages %d passes mean %.2f max %d sum %d CW %d | ' % (f, len(v), sum(v)/len(v), max(v), sum(v), cw) +
          ' '.join('t2s%s %.2f/%d' % (t, sum(x)/len(x), max(x)) for t, x in sorted(by.items())))
