#!/usr/bin/env python3
"""idort-1019: sp (vet_col + vet_gd wedge) >= 4-level self-convergence of the suite's linear-mode ladders, arms p (PROD) / ref.
Error e_N = L1(u_N - R u_2N) / mean|u_finest| (R = volume restriction; radial for rsw*/marsh/patmT, phi for rlat*), as cmpmk.xord.
Writes spxconv.json and fig/conv_sp_*.png"""
import os, sys, json, glob
os.environ['SPT_RUNS'] = '/raven/ptmp/jinma/idort1020/page/spxtree'
sys.path.insert(0, '/raven/ptmp/jinma/suiteF_1012/sp/ana'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cmpmk, spl, spt
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; import style as S
R = '/raven/ptmp/jinma/idort1020/runs/spx'; T = os.environ['SPT_RUNS']
for d in glob.glob(R + '/*/L*/*'):
    if not os.path.exists(d + '/rc.txt') or open(d + '/rc.txt').read().strip() != '0': continue
    case, lev, arm = d.split('/')[-3:]
    dst = '%s/%s/%s_%s' % (T, case, {'p': 'P', 'ref': 'R'}[arm], lev[1:])
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if not os.path.lexists(dst): os.symlink(d, dst)
    if not os.path.lexists(dst + '/log') and os.path.exists(d + '/run.log'): pass
def xerr(case, arm, levs):
    rows = [cmpmk.lastdata(case, arm, n) for n in levs]
    keys = [k for k in rows[0] if k != 'time' and k in rows[-1]]
    c = spt.case_cfg(case); out = {}
    for k in keys:
        errs = []
        for q in range(len(levs) - 1):
            a, b = rows[q][k], rows[q + 1][k]
            if c['test'] in ('drot', 'rlat', 'rlatt'):
                bb = 0.5 * (b[0::2] + b[1::2]); errs.append(float(np.mean(np.abs(a - bb))))
            else:
                rf = spl.rfaces(c['r'][0], c['r'][1], levs[q + 1], c['poly'])
                bb = spl.restrict_r(b, rf, c['geo'] != 'cart'); errs.append(float(np.mean(np.abs(a - bb))))
        sc = float(np.mean(np.abs(rows[-1][k]))) or 1.0
        out[k] = [e / sc for e in errs]
    return out
CASES = [('rsw', 'X'), ('rswt', 'X'), ('rswv', 'X'), ('rswtv', 'X'), ('rsw_C', 'C'), ('rswt_C', 'C'), ('rswv_C', 'C'), ('rswtv_C', 'C'),
         ('marsh_r1_u', 'X'), ('patmT_L16_s', 'C'), ('rlat', 'L'), ('rlatt', 'L'), ('rlat_om5', 'L'), ('rlatt_om5', 'L'),
         ('rswv_F', 'X'), ('rswtv_F', 'X'), ('rswv_F_C', 'C'), ('rswtv_F_C', 'C')]
res = {}
for case, kind in CASES:
    for arm in ('P', 'R'):
        levs = sorted(int(x.split('_')[-1]) for x in os.listdir(T + '/' + case) if x.startswith(arm + '_')) if os.path.isdir(T + '/' + case) else []
        # contiguous factor-2 chain ending at the finest level
        if not levs: continue
        ch = [levs[-1]]
        for x in reversed(levs[:-1]):
            if x * 2 == ch[0]: ch.insert(0, x)
            else: break
        levs = ch
        if len(levs) < 3: continue
        try:
            e = xerr(case, arm, levs)
        except Exception as ex:
            print('ERR', case, arm, ex); continue
        Ns = levs[:-1]
        res['%s|%s' % (case, arm)] = dict(kind=kind, N=Ns, levs=levs, err=e,
            order={k: [float(x) for x in np.log2(np.array(v[:-1]) / np.array(v[1:]))] for k, v in e.items()},
            fit={k: float(-np.polyfit(np.log(Ns), np.log(v), 1)[0]) for k, v in e.items()})
# dt-only sweep (rswt at N 256, dt / q): same-grid differences
try:
    qs = [1, 2, 4, 8, 16]; rows = [cmpmk.lastdata('rswt_dt', 'P', q) for q in qs]
    e = {k: [float(np.mean(np.abs(rows[i][k] - rows[i + 1][k]))) / (float(np.mean(np.abs(rows[-1][k]))) or 1.0) for i in range(len(qs) - 1)] for k in rows[0] if k != 'time'}
    res['rswt_dt|P'] = dict(kind='T', N=qs[:-1], levs=qs, err=e, order={k: [float(x) for x in np.log2(np.array(v[:-1]) / np.array(v[1:]))] for k, v in e.items()},
                            fit={k: float(-np.polyfit(np.log(qs[:-1]), np.log(v), 1)[0]) for k, v in e.items()})
except Exception as ex:
    print('ERR rswt_dt', ex)
json.dump(res, open('/raven/ptmp/jinma/idort1020/page/spxconv.json', 'w'), indent=1)
for k, r in res.items():
    for v in ('m1_e', 'dens'):
        if v in r['err']: print('%-18s %-5s %s N=%s e=%s pairs=%s fit=%.2f' % (k, v, r['kind'], r['N'], ' '.join('%.2e' % x for x in r['err'][v]), ' '.join('%.2f' % x for x in r['order'][v]), r['fit'][v]))
# ---- figure (page): 4 sp wave ladders, PROD
cases = ['rsw_C', 'rswt_C', 'rswv_F_C', 'rlat']
TT = {'rsw_C': 'radiation-acoustic pulse, rho kappa_s 1000', 'rswt_C': 'same, rho kappa_s 10 (semi-thin)', 'rswv_F_C': 'moving shell v_r 0.01 (approved deck)', 'rlat': 'lateral (phi) wave, rho kappa_s 1000'}
fig, axs = plt.subplots(2, 2, figsize=(15, 11))
for ax, c in zip(axs.flat, cases):
    r = res.get(c + '|P')
    if not r: continue
    for v, col, lab in (('m1_e', S.C1, 'E'), ('dens', S.C2, 'rho'), ('m1_f1', S.C3, 'F_r')):
        if v in r['err']: ax.loglog(r['N'], r['err'][v], marker='o', ms=8, color=col, label='%s  p = %.2f' % (lab, r['order'][v][-1]))
    S.slope2(ax); ax.set_title(TT[c], fontsize=13); ax.legend(fontsize=11); ax.set_xlabel('N (coarse level of the pair), dt ~ dx'); ax.set_ylabel('L1 self-convergence difference / mean')
fig.suptitle('Spherical wedge (vet_col + vet_gd), PROD keys: linear radiation-acoustic modes')
fig.tight_layout(); fig.savefig('/raven/ptmp/jinma/idort1020/page/fig/conv_sp.png', dpi=150)
