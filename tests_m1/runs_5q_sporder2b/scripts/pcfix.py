#!/usr/bin/env python3
"""pcfix.py prep <input> <probe input>: a copy with a restart output every cycle.
pcfix.py <probe run dir> <input>: when the probe's parameter dump (rst/bin header)
resolved implicit_precond = mg_gc, name implicit_precond = rbgs_fwd in <input> (the
default before m1-sp-order2b)."""
import glob
import os
import sys

if sys.argv[1] == 'prep':
    s = open(sys.argv[2]).read()
    s += '\n<output99>\nfile_type = rst\ndcycle = 1\n'
    open(sys.argv[3], 'w').write(s)
    sys.exit(0)
d, inp = sys.argv[1], sys.argv[2]
txt = ''
for f in glob.glob(os.path.join(d, '**', '*.rst'), recursive=True) + \
        glob.glob(os.path.join(d, '**', '*.bin'), recursive=True):
    b = open(f, 'rb').read()
    k = b.find(b'<par_end>')
    txt += b[:k].decode('latin-1') if k > 0 else ''
mg = any(ln.split('=')[0].strip() == 'implicit_precond' and 'mg_gc' in ln
         for ln in txt.splitlines())
if not txt:
    print('pcfix: no parameter dump in', d)
if mg:
    out = []
    for ln in open(inp):
        out.append(ln)
        if ln.strip() == '<rad_m1>':
            out.append('implicit_precond = rbgs_fwd\n')
    open(inp, 'w').write(''.join(out))
print('pcfix', os.path.basename(inp), 'mg_gc -> rbgs_fwd' if mg else 'unchanged')
