#!/usr/bin/env python3
"""ord_hebox.py <dir> <job> -- temporal self-convergence and ringing check of the He-box order runs.

Runs <dir>/<scheme><cfl>_<job>/ (scheme b = base, n = implicit_opac_newton; cfl tags 6 3 15 075) all end at
the same tlim and write a double-precision restart there.  For each scheme:
  e(c) = || u_c - u_{c/2} ||_1 / || u_{c/2} ||_1   (interior cells, per variable)
  order = log2(e(2c)/e(c)).
Also || u_b - u_n || at each cfl, relative to e(c) (the time error): Newton should sit far below it.
Ringing: from the per-step hst, the odd-even (period-2) content of the step-to-step increments,
  R = |sum (-1)^n d_n| / sum |d_n|  over the last 3/4 of the run, d_n = x_{n+1} - x_n,
(R -> 1 for a pure 2-step oscillation, ~ 1/sqrt(N) for noise, ~ 0 for a smooth drift),
and the largest |d_{n+1} + d_n| / (|d_{n+1}| + |d_n|) sign-flip fraction.
"""
import os, sys
import numpy as np
sys.path.insert(0, '/u/jma20/ATHENAK/scripts')
import cmp_rst as C

D, J = sys.argv[1], sys.argv[2]
CFL = ['6', '3', '15', '075']
VARS = None


def load(tag):
    d = f'{D}/{tag}_{J}/rst'
    f = sorted(x for x in os.listdir(d) if x.endswith('.rst'))[-1]
    R = C.Rst(os.path.join(d, f), True, None)
    segs = C.build_layout(R.info, R.D, None)
    fa = np.frombuffer(R.payload, dtype=R.fdt).reshape(R.nmb, -1)
    out, off = {}, 0
    ng = int(R.info['mb_indcs'][0])
    for name, varnames, shape in segs:
        n = int(np.prod(shape))
        if name in ('hydro_u0', 'rad_m1_u0'):
            a = fa[:, off:off + n].reshape((R.nmb,) + shape).astype(np.float64)
            a = a[:, :, ng:-ng, ng:-ng, ng:-ng]
            for v, vn in enumerate(varnames):
                out[vn.replace('hydro_u0.', '').replace('rad_m1_u0', 'm1')] = a[:, v]
        off += n
    return out, R.info['time']


def rel(a, b):
    return np.abs(a - b).sum() / max(np.abs(b).sum(), 1e-300)


SCH = [s for s in 'bn' if os.path.isdir(f'{D}/{s}{CFL[0]}_{J}')]
S = {}
for s in SCH:
    for c in CFL:
        S[s + c] = load(s + c)
t = {k: v[1] for k, v in S.items()}
print('final times:', ', '.join(f'{k}={v:.12g}' for k, v in t.items()))
names = list(S['b6'][0].keys())
for s, lab in [x for x in (('b', 'base'), ('n', 'newton')) if x[0] in SCH]:
    print(f'\n== {lab}: e(cfl) = ||u_c - u_c/2||_1/||u_c/2||_1   [order]')
    print(f"{'var':10s}" + ''.join(f'{"e(" + a + "/" + b + ")":>16s}' for a, b in zip(CFL, CFL[1:])) +
          f"{'p1':>7s}{'p2':>7s}")
    for vn in names:
        e = [rel(S[s + a][0][vn], S[s + b][0][vn]) for a, b in zip(CFL, CFL[1:])]
        p = [np.log2(e[i] / e[i + 1]) if e[i + 1] > 0 else float('nan') for i in range(len(e) - 1)]
        print(f'{vn:10s}' + ''.join(f'{x:16.3e}' for x in e) + ''.join(f'{x:7.2f}' for x in p))
if len(SCH) == 2:
  print('\n== ||u_base - u_newton||_1 / ||u||_1 at each cfl, and its ratio to the time error e(cfl)')
for vn in (names if len(SCH) == 2 else []):
    row = []
    for i, c in enumerate(CFL):
        d = rel(S['n' + c][0][vn], S['b' + c][0][vn])
        et = rel(S['b' + CFL[i]][0][vn], S['b' + CFL[i + 1]][0][vn]) if i + 1 < len(CFL) else float('nan')
        row.append(f'{d:10.2e} ({d / et if et == et and et > 0 else float("nan"):7.1e})')
    print(f'{vn:10s} ' + '  '.join(row))

print('\n== ringing (per-step hst): R = |sum (-1)^n d_n| / sum |d_n|, flip = fraction of sign flips of d_n')
for f, cols in (('m1slab.hydro.hst', {7: 'tot-E', 8: '1-KE', 9: '2-KE'}),
                ('m1slab.user.hst', {3: 'F1top', 5: 'F1bot', 6: 'V1max', 7: 'Etot', 8: 'Fres'})):
    for s in SCH:
        for c in CFL:
            h = np.loadtxt(f'{D}/{s}{c}_{J}/{f}')
            n0 = len(h) // 4
            parts = []
            for ci, nm in cols.items():
                x = h[n0:, ci - 1]
                d = np.diff(x)
                sgn = (-1.0) ** np.arange(len(d))
                R = abs((sgn * d).sum()) / max(np.abs(d).sum(), 1e-300)
                flip = np.mean(np.sign(d[1:]) != np.sign(d[:-1])) if len(d) > 1 else 0
                parts.append(f'{nm}:R={R:.3f},flip={flip:.2f}')
            print(f'{s}{c:4s} n={len(h):5d} ' + '  '.join(parts))
