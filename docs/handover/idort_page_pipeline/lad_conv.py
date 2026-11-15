#!/usr/bin/env python3
"""idort-1019: Cartesian self-convergence ladders (suite lad/: rw_*, pulse_k*, atm) from od.py output eval/out_<arm>/lad.txt
-> ladconv.json + fig/conv_lad_cart.png (X = space at fixed dt, C = dt ~ dx, T = dt refinement at fixed N)"""
import re, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; import style as S
E = '/raven/ptmp/jinma/idort1020/eval'
def parse(fn):
    out = {}
    for ln in open(fn):
        m = re.match(r'^(\S+)\s+(\S+)\s+(\S+)\s+([XCT])\s+(\S+)\s+lev (\S+)\s+L1 ([^|]+)\|\s*ord(.*)$', ln)
        if not m: continue
        case, arm, ts, kind, var, lev, l1, od = m.groups()
        out['%s|%s|%s|%s' % (case, arm, kind, var)] = dict(lev=[int(x) for x in lev.split('-')], l1=[float(x) for x in l1.split()], ord=[float(x) for x in od.split()])
    return out
R = {a: parse('%s/out_%s/lad.txt' % (E, a)) for a in ('p', 'ref') if os.path.exists('%s/out_%s/lad.txt' % (E, a))}
V10 = parse('/raven/ptmp/jinma/idort_page/cache/fin3ana/lad_FIN.txt') if os.path.exists('/raven/ptmp/jinma/idort_page/cache/fin3ana/lad_FIN.txt') else {}
json.dump(dict(p=R.get('p', {}), ref=R.get('ref', {}), v10=V10), open('/raven/ptmp/jinma/idort1020/page/ladconv.json', 'w'), indent=0)
cases = sorted({k.split('|')[0] for k in R.get('p', {})})
fig, axs = plt.subplots(len(cases), 3, figsize=(14, 3.1 * len(cases)), squeeze=False)
for i, c in enumerate(cases):
    for j, kind in enumerate('XCT'):
        ax = axs[i, j]
        for arm, mf in (('p', True), ('ref', False)):
            for var, col in (('e', S.C1), ('dens', S.C2), ('f1', S.C3)):
                ks = [k for k in R.get(arm, {}) if k.split('|')[0] == c and k.split('|')[1] in ('FIN', 'FINs') and k.split('|')[2] == kind and k.split('|')[3] == var]
                if not ks: continue
                r = R[arm][ks[0]]; x = r['lev'][:-1]
                ax.loglog(x, r['l1'], marker='o', ls='-' if mf else ':', color=col, mfc=col if mf else 'white',
                          label='%s %s p=%.2f' % ('PROD' if mf else 'ref', var, r['ord'][-1]))
        if ax.get_lines(): S.slope2(ax); ax.legend(fontsize=7)
        ax.set_title('%s %s' % (c, {'X': 'X space (fixed dt)', 'C': 'C dt ~ dx', 'T': 'T dt / l (fixed N)'}[kind]), fontsize=9.5)
fig.suptitle('Cartesian vet_sc self-convergence ladders (hesdirk2): L1 of successive-level differences; p = finest-pair order')
fig.tight_layout(); fig.savefig('/raven/ptmp/jinma/idort1020/page/fig/conv_lad_cart.png', dpi=150)
for a in R:
    for k, r in sorted(R[a].items()):
        print(a, k, ' '.join('%.2e' % x for x in r['l1']), '| ord', ' '.join('%.2f' % x for x in r['ord']))
