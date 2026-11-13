#!/usr/bin/env python3
"""setkeys.py <in> <out> blk/key=val ...: copy, SETTING each key (replace the value if present, else add it to its block)."""
import sys
src, dst, kv = sys.argv[1], sys.argv[2], sys.argv[3:]
add = {}
for a in kv:
    k, v = a.split('=', 1); b, key = k.split('/', 1); add.setdefault(b, {})[key] = v
out, blk = [], None
for ln in open(src).read().split('\n'):
    s = ln.strip()
    if s.startswith('<') and s.endswith('>'):
        blk = s[1:-1]; out.append(ln)
        for key, v in add.get(blk, {}).items():
            out.append('%-24s = %s' % (key, v))
        continue
    if blk in add and '=' in s and not s.startswith('#') and s.split('=')[0].strip() in add[blk]:
        continue
    out.append(ln)
open(dst, 'w').write('\n'.join(out))
