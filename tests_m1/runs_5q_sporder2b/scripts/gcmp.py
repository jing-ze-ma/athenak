#!/usr/bin/env python3
"""gcmp.py dirA dirB: bitwise compare of every output file (bin/rst from <par_end> on,
the header echoes the input).  One line per pair and a verdict (from runs_5h_sph2)."""
import os
import sys


def payload(p):
    d = open(p, 'rb').read()
    if p.endswith(('.bin', '.rst')):
        k = d.find(b'<par_end>')
        if k >= 0:
            d = d[k:]
    return d


a, b = sys.argv[1], sys.argv[2]
nf = nd = 0
for root, _, files in os.walk(a):
    for f in sorted(files):
        if f == 'log.txt' or f.startswith(('slurm', 'log.')):
            continue
        pa = os.path.join(root, f)
        pb = os.path.join(b, os.path.relpath(pa, a))
        if not os.path.exists(pb):
            print('MISSING', pb)
            nd += 1
            continue
        nf += 1
        same = payload(pa) == payload(pb)
        nd += (not same)
        print(('same ' if same else 'DIFF ') + os.path.relpath(pa, a))
print(f'VERDICT {"BITWISE" if nd == 0 and nf > 0 else "DIFFERENT"} files={nf} diff={nd}')
