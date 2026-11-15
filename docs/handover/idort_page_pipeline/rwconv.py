#!/usr/bin/env python3
"""idort-1019: Cartesian linear radiation-modified acoustic wave convergence (exact damped eigenmode, one period, dt ~ dx).
Reads results.json of idort1019/rw/runs (PROD19, TIGHT, xy) and v10 (idort_page/reruns/fin3/rw, arm FIN). Writes
fig/conv_rw_cart.png and rwconv.json {group: {N: [...], l1_rho: [...], l1_E: [...], order_fit, order_fine}}."""
import json, glob, os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import style as S
P = '/raven/ptmp/jinma/idort1020'
src = [(P + '/rw/runs', None)]
G = {}
for root, lab in src:
    for f in glob.glob(root + '/*/results.json'):
        r = json.load(open(f)); c = r['case']; arm = lab or c['extra'].get('arm')
        if lab == 'v10' and c['extra'].get('arm') != 'FIN': continue
        if c['nt'] != 2 * c['N'] or any(x in c['extra'] for x in ('nmu', 'cfl')) or (lab is None and arm == 'FIN'): continue
        k = '%s|%s|P%g|t%g' % (arm, c['dir'], c['prat'], c['tau'])
        G.setdefault(k, {})[c['N']] = (r['l1_rho'], r['l1_E'], r.get('ncycle'))
def order(N, e):
    N, e = np.array(N, float), np.array(e, float)
    fit = -np.polyfit(np.log(N), np.log(e), 1)[0] if len(N) >= 2 else np.nan
    fine = np.log(e[-2] / e[-1]) / np.log(N[-1] / N[-2]) if len(N) >= 2 else np.nan
    pairs = [np.log(e[q] / e[q + 1]) / np.log(N[q + 1] / N[q]) for q in range(len(N) - 1)]
    return float(fit), float(fine), [float(x) for x in pairs]
out = {}
for k, d in sorted(G.items()):
    N = sorted(d); er = [d[n][0] for n in N]; ee = [d[n][1] for n in N]
    out[k] = dict(N=N, l1_rho=er, l1_E=ee, rho=order(N, er), E=order(N, ee))
json.dump(out, open(P + '/page/rwconv.json', 'w'), indent=1)
Ps, Ts = (0.00333333, 3.33333), (0.01, 1.0, 100.0)
sty = {'PROD19|x1': (S.C1, 'o', '-', 'PROD 1-D'), 'PROD19|xy': (S.C4, '^', '-', 'PROD 2-D (xy)')}
for p, tag in zip(Ps, ('P001', 'P10')):
    fig, axs = plt.subplots(2, 3, figsize=(15, 10), sharex=True)
    for i, VAR in enumerate(('rho', 'E')):
        for j, t in enumerate(Ts):
            ax = axs[i, j]
            for sk, (col, m, ls, lab) in sty.items():
                arm, dr = sk.split('|')
                kk = [k for k in out if k.startswith('%s|%s|' % (arm, dr)) and abs(float(k.split('|')[2][1:]) / p - 1) < 1e-3 and abs(float(k.split('|')[3][1:]) / t - 1) < 1e-6]
                if not kk: continue
                o = out[kk[0]]
                ax.loglog(o['N'], o['l1_' + VAR], marker=m, ms=8, ls=ls, color=col, label='%s  p = %.2f' % (lab, o[VAR][1]))
            S.slope2(ax)
            ax.set_xticks([16, 32, 64, 128, 256]); ax.set_xticklabels(['16', '32', '64', '128', '256']); ax.minorticks_off()
            ax.set_title('tau_lambda = %g' % t)
            ax.legend(fontsize=11, loc='lower left')
            if i == 1: ax.set_xlabel('N cells per wavelength (dt ~ dx)')
            if j == 0: ax.set_ylabel('L1(%s - exact) / amplitude' % VAR)
    fig.suptitle('Radiation-modified acoustic waves, P_rad/P_gas = %s (cf. Ma+ Sect. 4.3.2)' % ('3.3e-3' if p < 1 else '3.33'), fontsize=17)
    fig.tight_layout(); fig.savefig(P + '/page/fig/conv_rw_%s.png' % tag, dpi=150)
print('\n'.join('%-28s N=%s rho fit %.2f fine %.2f | E fit %.2f fine %.2f' % (k, o['N'], o['rho'][0], o['rho'][1], o['E'][0], o['E'][1]) for k, o in out.items()))
