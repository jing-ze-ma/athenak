#!/usr/bin/env python
"""vgdrot-1011 per-run timing: wall time per cycle c = elapsed(c+1) - elapsed(c) from the
'elapsed=... cycle=...' lines.  Cycles K0 .. last-1 (the LAST step, which carries the end-of-run
outputs, is excluded).  steady = c % ROT != 0, rot = c % ROT == 0 (the vet_gd rotation cycles).
usage: ana.py run.log [K0] [ROT]"""
import re, sys
f = sys.argv[1]; k0 = int(sys.argv[2]) if len(sys.argv) > 2 else 2
rot = int(sys.argv[3]) if len(sys.argv) > 3 else 10
t = {}
for ln in open(f, errors='replace'):
    m = re.search(r'elapsed=([0-9.eE+-]+) cycle=(\d+)', ln)
    if m: t[int(m.group(2))] = float(m.group(1))
cs = sorted(t)
if len(cs) < 3: print('ncyc=%d NA' % len(cs)); sys.exit()
last = cs[-1]
# cycle c = the step c -> c+1; the last step (last-1 -> last) excluded
dt = {c: t[c+1] - t[c] for c in cs if c >= k0 and c + 1 < last}
st = [v for c, v in dt.items() if c % rot != 0]
ro = [v for c, v in dt.items() if c % rot == 0]
al = list(dt.values())
mean = lambda a: sum(a)/len(a) if a else float('nan')
print('ncyc=%d steady=%.4f (n=%d) rot=%.4f (n=%d: %s) avg=%.4f (n=%d) init=%.2f' % (
    last, mean(st), len(st), mean(ro), len(ro), ' '.join('%.3f' % v for c, v in sorted(dt.items()) if c % rot == 0),
    mean(al), len(al), t[cs[0]]))
