#!/usr/bin/env python3
"""idort-1019: cs + sp diffusion-mode order tests (csport W1 harness): space (N 8..64, dt ~ dx) and time (N16, dt/q) self-convergence,
error classes all / interior / seam / vertex (cs).  arms r19 (harness keys) and p19 (PROD).  -> csconv.json, fig/conv_cs.png"""
import os, sys, json, io, contextlib
sys.path.insert(0, '/raven/ptmp/jinma/idort1020/cs/py'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, csconv as CC
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; import style as S
def ev(mesh, case, tag, kind):
    os.environ['CS_MESH'] = mesh; os.environ['CS_TAG'] = tag
    runs = [(s, None) for s in (1, 2, 4, 8) if os.path.exists(CC.deck(case, s, None)[0] + '/rc.txt') and open(CC.deck(case, s, None)[0] + '/rc.txt').read().strip() == '0'] \
        if kind == 'space' else [(2, q) for q in (1, 2, 4, 8) if os.path.exists(CC.deck(case, 2, q)[0] + '/rc.txt') and open(CC.deck(case, 2, q)[0] + '/rc.txt').read().strip() == '0']
    # keep the contiguous chain
    if len(runs) < 3: return None
    E, V, C = {}, {}, {}
    for s, q in runs:
        d = CC.deck(case, s, q)[0]
        E[(s, q)], _ = (CC.load_sp if mesh == 'sp' else CC.load)(d)
        V[(s, q)], C[(s, q)] = (CC.geom_sp if mesh == 'sp' else CC.geom)(8 * s, 16 * s)
    names = ['all', 'interior', 'seam', 'vertex']; res = {n: [] for n in names}
    for a, b in zip(runs[:-1], runs[1:]):
        Ea = E[a]; Eb = CC.restrict(E[b], V[b], True) if kind == 'space' else E[b]; w = V[a]
        mean = np.sum(w * Ea) / np.sum(w); den = np.sum(w * np.abs(Ea - mean)) / np.sum(w); dif = np.abs(Ea - Eb)
        for n, c in zip(names, [None, 0, 1, 2]):
            msk = np.ones_like(dif, bool) if c is None else (C[a] == c)
            if msk.sum() == 0: continue
            res[n].append(float(np.sum((w * dif)[msk]) / np.sum(w[msk]) / den))
    res = {k: v for k, v in res.items() if v}
    x = [8 * s for s, q in runs[:-1]] if kind == 'space' else [q for s, q in runs[:-1]]
    return dict(x=x, err=res, order={k: [float(np.log2(a / b)) for a, b in zip(v[:-1], v[1:])] for k, v in res.items()})
out = {}
for mesh in ('cs', 'sp'):
    for case in ('gd', 'gda', 'edd'):
        for tag in ('p19', 'r19', 'pnr19'):
            for kind in ('space', 'time'):
                try:
                    with contextlib.redirect_stdout(io.StringIO()):
                        r = ev(mesh, case, tag, kind)
                except Exception as ex:
                    print('ERR', mesh, case, tag, kind, ex); r = None
                if r:
                    out['%s|%s|%s|%s' % (mesh, case, tag, kind)] = r
                    print('%s|%s|%s|%-5s x=%s ' % (mesh, case, tag, kind, r['x']) + ' ; '.join('%s %s (p %s)' % (k, ' '.join('%.2e' % e for e in v), ' '.join('%.2f' % o for o in r['order'][k])) for k, v in r['err'].items()))
json.dump(out, open('/raven/ptmp/jinma/idort1020/page/csconv.json', 'w'), indent=1)
fig, axs = plt.subplots(2, 2, figsize=(15, 11))
LAB = {'gd': 'thick diffusion mode (rho kappa_s 1000)', 'gda': 'diffusion mode in an atmosphere (tau 20, Marshak top)'}
for j, case in enumerate(('gd', 'gda')):
    for i, kind in enumerate(('space', 'time')):
        ax = axs[i, j]
        for mesh, cls, col, m in (('cs', 'all', S.C1, 'o'), ('cs', 'seam', S.C2, 'o'), ('cs', 'vertex', S.C4, 'o'), ('sp', 'all', S.C3, 's')):
            r = out.get('%s|%s|p19|%s' % (mesh, case, kind))
            if not r or cls not in r['err']: continue
            x = r['x'] if kind == 'time' else r['x']
            ax.loglog(x, r['err'][cls], marker=m, ms=8, color=col, label='%s %s cells  p = %.2f' % ('cubed sphere' if mesh == 'cs' else 'sp wedge', cls, r['order'][cls][-1]))
        if ax.get_lines(): S.slope2(ax)
        ax.set_title('%s: %s' % (LAB[case], 'space (dt ~ dx)' if kind == 'space' else 'time (N 16, dt/q)'), fontsize=13)
        ax.set_xlabel('N panel cells (coarse level of the pair)' if kind == 'space' else 'q (coarse level of the pair)'); ax.set_ylabel('L1 self-convergence difference'); ax.legend(fontsize=11)
fig.suptitle('Production closure (vet_col + vet_gd) on the cubed sphere and the sp wedge, PROD keys: l = 2 radiation diffusion mode')
fig.tight_layout(); fig.savefig('/raven/ptmp/jinma/idort1020/page/fig/conv_cs.png', dpi=150)
