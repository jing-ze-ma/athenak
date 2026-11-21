#!/usr/bin/env python3
"""seam.py BINDIR RPH [BINVIS]: periodic-seam artifact metric from AthenaK bins (hydro_w + m1, last output).
For each variable and radial band (r_ph +-5 cells; atmosphere 1.1-1.4 r_ph): ratio = mean |q(seam pair)| jump / median interior |dq| (cell to cell, same direction).
phi seam: k=N3-1 vs k=0; theta seam: j=N2-1 vs j=0. Ratio ~1 = no seam artifact."""
import sys, glob, numpy as np
sys.path.insert(0, sys.argv[3] if len(sys.argv) > 3 else '.')
import bin_convert as bc
D, rph = sys.argv[1], float(sys.argv[2])
def glob_arr(f):
    d = bc.read_binary(f); lg = np.array(d['mb_logical']); n3, n2, n1 = np.shape(d['mb_data'][d['var_names'][0]][0])
    N = [(lg[:, i].max() + 1) * s for i, s in zip(range(3), (n1, n2, n3))]
    x1 = np.zeros(N[0]); out = {}
    for v in d['var_names']:
        A = np.zeros((N[2], N[1], N[0]))
        for b, (l1, l2, l3, _) in enumerate(lg):
            A[l3*n3:(l3+1)*n3, l2*n2:(l2+1)*n2, l1*n1:(l1+1)*n1] = d['mb_data'][v][b]
            x1[l1*n1:(l1+1)*n1] = d['mb_x1v'][b]
        out[v] = A
    return out, x1
for kind in ['hydro_w', 'm1']:
    fs = sorted(glob.glob('%s/*.%s.*.bin' % (D, kind)))
    if not fs: continue
    Q, x1 = glob_arr(fs[-1]); i0 = np.searchsorted(x1, rph)
    bands = {'r_ph': slice(max(i0-5, 0), i0+5), 'atm': slice(np.searchsorted(x1, 1.1*rph), max(np.searchsorted(x1, 1.4*rph), np.searchsorted(x1, 1.1*rph)+1))}
    for v, A in Q.items():
        s = []
        for bn, sl in bands.items():
            B = A[:, :, sl]
            for ax, nm in ((0, 'phi'), (1, 'th')):
                d_in = np.abs(np.diff(B, axis=ax)); med = np.median(d_in) + 1e-300
                seam = np.abs(np.take(B, -1, axis=ax) - np.take(B, 0, axis=ax)).mean()
                s.append('%s/%s %.2f' % (bn, nm, seam / med))
        print('%-8s %-8s %s' % (fs[-1].split('.')[-2], v, '  '.join(s)))
