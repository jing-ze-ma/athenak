#!/usr/bin/env python3
"""reblock_identity.py (fastfloor-1109) RST DECK OUTDIR NAME : identity remap (same grid) of an AthenaK he_star_m1 + rad_m1 restart into a
problem/he_ic_remap_file (HSREMAP1) fresh start, so the SAME state can start on any MeshBlock layout.
Arrays copied verbatim from the restart (u0 hydro incl. rho*Phi, M1 u0, f0x1), lateral ghosts = periodic images,
x1 ghosts = the restart's.  Writes NAME.remap and NAME.athinput (DECK + start_time + he_ic_remap_file); the
meshblock size is set on the command line of each run.  Uses read_rst/gather of tools/e1remap/remap_e1w.py."""
import hashlib
import os
import re
import struct
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from remap_e1w import read_rst, gather  # noqa: E402

rst, deck, outdir, name = sys.argv[1:5]
R = read_rst(rst)
print('rst t=%.17g cycle %d nmb %d mesh %s mb %s ng %d nh %d counts %s' % (
    R['time'], R['ncycle'], R['nmb'], R['mesh'], R['mb'], R['ng'], R['nh'], R['counts']))
NX1, NX2, NX3 = R['mesh']; ng = R['ng']; par = R['params']
hyd = gather(R, 'hyd', R['nh'])[:5]
m1 = gather(R, 'm1', 4)
f0 = gather(R, 'f0x1', 1, face1=True)[0]
pg = R['blocks'].get(b'PGENST01', b'')


def padlat(x):
    return np.pad(x, ((ng, ng), (ng, ng), (0, 0)), mode='wrap')


os.makedirs(outdir, exist_ok=True)
fn = os.path.join(outdir, name + '.remap')
hdr = b'HSREMAP1' + struct.pack('<8i', 1, ng, NX1, NX2, NX3, 5, 4, 0) + struct.pack(
    '<7d', R['time'], float(par['mesh/x1min']), float(par['mesh/x1max']), float(par['mesh/x2min']),
    float(par['mesh/x2max']), float(par['mesh/x3min']), float(par['mesh/x3max']))
with open(fn, 'wb') as f:
    f.write(hdr); f.write(struct.pack('<Q', len(pg))); f.write(pg)
    for n in range(5):
        f.write(np.ascontiguousarray(padlat(hyd[n]), dtype='<f8').tobytes())
    for n in range(4):
        f.write(np.ascontiguousarray(padlat(m1[n]), dtype='<f8').tobytes())
    f.write(np.ascontiguousarray(padlat(f0), dtype='<f8').tobytes())
print('wrote', fn, os.path.getsize(fn), 'pgen block', len(pg))
d = open(deck).read()


def setkey(txt, blk, key, val):
    lines = txt.split('\n'); bi = None
    for n, ln in enumerate(lines):
        if re.match(r'\s*<\s*%s\s*>' % re.escape(blk), ln):
            bi = n; continue
        if bi is not None and re.match(r'\s*<', ln):
            break
        if bi is not None and re.match(r'\s*%s\s*=' % re.escape(key), ln):
            lines[n] = '%-24s = %s   # reblock_identity.py' % (key, val); return '\n'.join(lines)
    lines.insert(bi + 1, '%-24s = %s   # reblock_identity.py' % (key, val)); return '\n'.join(lines)


d = setkey(d, 'time', 'start_time', repr(R['time']))
d = setkey(d, 'problem', 'he_ic_remap_file', os.path.abspath(fn))
open(os.path.join(outdir, name + '.athinput'), 'w').write(d)
print('wrote deck', os.path.join(outdir, name + '.athinput'))
