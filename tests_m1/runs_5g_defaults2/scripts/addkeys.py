#!/usr/bin/env python3
# usage: addkeys.py in out block key=val ...   (insert keys after <block>)
import sys
src, dst, blk = sys.argv[1:4]
kv = sys.argv[4:]
out = []
for ln in open(src):
    out.append(ln)
    if ln.strip() == '<%s>' % blk:
        out += ['%s = %s\n' % tuple(x.split('=', 1)) for x in kv]
open(dst, 'w').write(''.join(out))
