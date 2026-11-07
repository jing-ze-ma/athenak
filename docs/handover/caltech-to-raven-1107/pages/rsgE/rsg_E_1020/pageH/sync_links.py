#!/usr/bin/env python3
"""pagefast-1102: the link part of pageH/sync.sh in one process (was ~8500 ln/sed forks, ~3.5 min on NFS).
Same link set as sync.sh.pre_fast_1102: every prodH/bin/*.bin, plus prodG/bin/rsg.<t>.NNNNN.bin with NNNNN below the first prodH
index of type <t> (only types present in prodH).  Existing correct links are left alone; links are never removed."""
import os, re, sys
R = sys.argv[1]; P = R + '/pageH'; H = R + '/prodH/bin'; G = R + '/prodG/bin'; B = P + '/bin'
pat = re.compile(r'^rsg\.([a-z0-9_]+)\.([0-9]+)\.bin$')
want = {}
hs = sorted(f for f in os.listdir(H) if f.endswith('.bin'))
for f in hs:
    want[f] = H + '/' + f
first = {}
for f in hs:
    m = pat.match(f)
    if m:
        t, i = m.group(1), m.group(2)
        if t not in first or i < first[t]:
            first[t] = i
for f in sorted(os.listdir(G)):
    m = pat.match(f)
    if m and m.group(1) in first and int(m.group(2), 10) < int(first[m.group(1)], 10):
        want[f] = G + '/' + f
have = {}
with os.scandir(B) as it:
    for e in it:
        if e.is_symlink():
            have[e.name] = os.readlink(e.path)
n = 0
for f, tgt in want.items():
    if have.get(f) != tgt:
        tmp = B + '/.' + f + '.lnk'
        try:
            os.remove(tmp)
        except FileNotFoundError:
            pass
        os.symlink(tgt, tmp); os.replace(tmp, B + '/' + f); n += 1
print('sync_links: %d new/changed of %d' % (n, len(want)))
