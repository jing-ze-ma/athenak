"""setkeys.py deck block/key=val ... : set (or add) keys in an athinput file in place"""
import sys
fn = sys.argv[1]
kv = [a.split('=', 1) for a in sys.argv[2:]]
lines = open(fn).read().split('\n')
for bk, v in kv:
    b, k = bk.split('/')
    blk, done, ins = '', False, None
    for n, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith('<'):
            if blk == b and ins is None:
                ins = n
            blk = s.strip('<>').strip()
            continue
        if blk == b and '=' in s and s.split('=')[0].strip() == k:
            lines[n] = '%s = %s' % (k, v)
            done = True
    if not done:
        if ins is None:
            ins = len(lines)
        lines.insert(ins, '%s = %s' % (k, v))
open(fn, 'w').write('\n'.join(lines))
