#!/usr/bin/env python3
"""CSEDGE (raven 10-19) cs edge-row diagnostic.  usage: csedge_diag.py <athenak vis/python dir> <bin dir> [basename=rsg]
For every <base>.hydro_w.NNNNN.bin: rms of the ALONG-EDGE tangential velocity (contravariant v3 on xi-edges, v2 on eta-edges)
in the edge row (row0), the next row (row1) and a panel-interior line, at radial cells i = 0, 5, 10, 20, 40, plus rms v_r at the
cube-vertex 2x2 cells vs interior at i = 0.  The cs edge instability = row0 >> row1 ~ interior, growing e-fold ~0.5 d (N=128)."""
import sys, glob, os
import numpy as np
sys.path.insert(0, sys.argv[1])
import bin_convert as bc
base = sys.argv[3] if len(sys.argv) > 3 else 'rsg'
I = [0, 5, 10, 20, 40]
def asm(d, var, ni):
    nmb = len(d['mb_logical']); per = nmb // 6
    n3, n2, n1 = d['mb_data'][var][0].shape; nb = int(round(np.sqrt(per))); N = n2 * nb
    P = np.empty((6, N, N, ni))
    for m in range(nmb):
        p = m // per; lx1, lx2, lx3, lev = d['mb_logical'][m]
        P[p, lx3*n3:(lx3+1)*n3, lx2*n2:(lx2+1)*n2, :] = d['mb_data'][var][m][..., :ni]
    return P
rms = lambda a: float(np.sqrt(np.mean(a**2)))
for f in sorted(glob.glob(os.path.join(sys.argv[2], '%s.hydro_w.*.bin' % base))):
    h = bc.read_binary(f); ni = max(I) + 1
    v1 = asm(h, 'velx', ni); v2 = asm(h, 'vely', ni); v3 = asm(h, 'velz', ni)
    N = v2.shape[1]; s = slice(N//8, N - N//8)
    out = []
    for i in I:
        r0 = np.concatenate([v3[:, s, 0, i], v3[:, s, N-1, i], v2[:, 0, s, i], v2[:, N-1, s, i]], axis=1)
        r1 = np.concatenate([v3[:, s, 1, i], v3[:, s, N-2, i], v2[:, 1, s, i], v2[:, N-2, s, i]], axis=1)
        ri = np.concatenate([v3[:, s, N//2, i], v2[:, N//2, s, i]], axis=1)
        out.append('i%d %.3g/%.3g/%.3g' % (i, rms(r0), rms(r1), rms(ri)))
    vt = np.concatenate([v1[:, a, b, 0] for a in (slice(0, 2), slice(N-2, N)) for b in (slice(0, 2), slice(N-2, N))], axis=1)
    print('%s t %.3f d | along-edge rms row0/row1/interior cm/s: %s | vr vtx/int i0 %.3g/%.3g'
          % (os.path.basename(f), float(h['time'])/86400, ' '.join(out), rms(vt), rms(v1[:, s, s, 0])), flush=True)
