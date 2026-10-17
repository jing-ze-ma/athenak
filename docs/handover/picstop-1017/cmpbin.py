#!/usr/bin/env python3
"""picstop-1017: solution difference of arms vs a reference at the same time (last bin of each kind).
usage: cmpbin.py REFDIR DIR [DIR ...] [--ib i0:i1]   (needs bin_convert.py next to it or in vis/python)
Per arm and variable (dens, eint (or press), m1_e, m1_f1): max |a-r|/|r| and mean |a-r|/|r| over all cells and over the radial
band i0:i1 (default 190:240, the AG Car tail cells), and the time check (both files at the same t)."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, '/raven/u/jinma/ATHENAK/wt-picstop/vis/python')
import bin_convert as bc
args = [a for a in sys.argv[1:] if not a.startswith('--')]
ib = (190, 240)
for i, a in enumerate(sys.argv):
    if a == '--ib': ib = tuple(int(x) for x in sys.argv[i + 1].split(':'))
args = [a for a in args if ':' not in a]
def load(d, kind):
    fs = sorted(glob.glob(d + '/bin/*.%s.*.bin' % kind))
    if not fs: return None, None
    f = bc.read_binary(fs[-1]); out = {}
    for v in f['var_names']:
        out[v] = np.concatenate([np.asarray(x).ravel() for x in f['mb_data'][v]]) if f['n_mbs'] > 1 else np.asarray(f['mb_data'][v][0])
    return f['time'], out
ref = args[0]
R = {k: load(ref, k) for k in ('hydro_w', 'm1')}
print('ref', ref, 't', R['hydro_w'][0], R['m1'][0], '| band i', ib)
for d in args[1:]:
    line = []
    for k, vs in (('hydro_w', ('dens', 'eint', 'press')), ('m1', ('m1_e', 'm1_f1'))):
        t, A = load(d, k); tr, B = R[k]
        if A is None or B is None: line.append('%s: missing' % k); continue
        if abs(t - tr) > 1e-9 * abs(tr): line.append('%s: t %.6e != ref %.6e' % (k, t, tr)); continue
        for v in vs:
            if v not in A or v not in B: continue
            a, b = A[v], B[v]
            den = np.maximum(np.abs(b), 1e-300); r = np.abs(a - b) / den
            if v == 'm1_f1': r = np.abs(a - b) / np.maximum(np.abs(B['m1_e']) * 3e10 if 'm1_e' in B else den, 1e-300)
            rb = r[..., ib[0]:ib[1]] if r.ndim == 3 else r
            line.append('%s max %.1e mean %.1e | band max %.1e mean %.1e' % (v, r.max(), r.mean(), rb.max(), rb.mean()))
    print('%-12s ' % os.path.basename(d.rstrip('/')) + ' ; '.join(line))
