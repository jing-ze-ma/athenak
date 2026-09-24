#!/usr/bin/env python3
# usage: setkeys.py in out block key=val ...  (replace key in <block> if present, else add)
import sys
src, dst, blk = sys.argv[1:4]
kv = dict(x.split('=', 1) for x in sys.argv[4:])
out, cur, done = [], None, set()
lines = open(src).read().splitlines(True)
for ln in lines:
    s = ln.strip()
    if s.startswith('<') and s.endswith('>'):
        if cur == blk:
            out += [f'{k} = {v}\n' for k, v in kv.items() if k not in done]
            done |= set(kv)
        cur = s[1:-1]
        out.append(ln)
        continue
    key = s.split('=')[0].strip() if '=' in s and not s.startswith('#') else None
    if cur == blk and key in kv:
        out.append(f'{key} = {kv[key]}\n')
        done.add(key)
    else:
        out.append(ln)
if cur == blk:
    out += [f'{k} = {v}\n' for k, v in kv.items() if k not in done]
open(dst, 'w').write(''.join(out))
