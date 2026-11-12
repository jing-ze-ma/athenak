#!/usr/bin/env python3
# mkdeck.py IN OUT block/key=value ...  : insert (or replace) keys in their <block> of an athinput
import sys, re
src, dst, kv = sys.argv[1], sys.argv[2], sys.argv[3:]
L = open(src).read().split('\n')
for x in kv:
    bk, v = x.split('=', 1); b, k = bk.split('/', 1)
    hdr = None
    for i, l in enumerate(L):
        if re.match(r'<\s*%s\s*>' % re.escape(b), l.strip()): hdr = i; break
    if hdr is None: L += ['<%s>' % b, '%s = %s' % (k, v)]; continue
    j = hdr + 1; done = False
    while j < len(L) and not L[j].strip().startswith('<'):
        if re.match(r'\s*%s\s*=' % re.escape(k), L[j]): L[j] = '%s = %s' % (k, v); done = True; break
        j += 1
    if not done: L.insert(hdr + 1, '%s = %s' % (k, v))
open(dst, 'w').write('\n'.join(L))
