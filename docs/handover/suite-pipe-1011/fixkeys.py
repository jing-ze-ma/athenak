#!/usr/bin/env python3
"""fixkeys.py <in> <out> blk/key=val ...: copy <in> to <out>, adding every command-line key missing
from its block (athenak: a command-line key absent from a fresh input is FATAL).  Values of keys already
present are left alone (the command line overrides them).  Prints the added keys."""
import sys
src, dst, kv = sys.argv[1], sys.argv[2], sys.argv[3:]
lines = open(src).read().split('\n')
have, blk = {}, None
for ln in lines:
    s = ln.strip()
    if s.startswith('<') and s.endswith('>'):
        blk = s[1:-1]; have.setdefault(blk, set()); continue
    if blk and '=' in s and not s.startswith('#'):
        have[blk].add(s.split('=')[0].strip())
add = {}
ADDED = {}
for a in kv:
    if '=' not in a or '/' not in a.split('=')[0]:
        continue
    k, v = a.split('=', 1)
    b, key = k.split('/', 1)
    if key not in have.get(b, set()):
        add.setdefault(b, []).append((key, v)); ADDED.setdefault(b, []).append((key, v))
out = []
for ln in lines:
    out.append(ln)
    s = ln.strip()
    if s.startswith('<') and s.endswith('>') and s[1:-1] in add:
        for key, v in add.pop(s[1:-1]):
            out.append('%-24s = %s' % (key, v))
for b, kvs in add.items():   # block absent
    out.append('<%s>' % b)
    out += ['%-24s = %s' % (key, v) for key, v in kvs]
open(dst, 'w').write('\n'.join(out))
print('added:', ' '.join('%s/%s=%s' % (b, k, v) for b, kv2 in ADDED.items() for k, v in kv2))
