#!/usr/bin/env python3
"""10-21 USER: IDORT linear waves, (P_rad/P_gas, tau_lambda) grid with the production method (PROD19, x1, dt ~ dx).
Reads rwconv.json (run rwconv.py first). Writes fig/conv_rw_grid.png (error vs N per (P, tau)) and
fig/conv_rw_ordermap.png (finest-pair order map for rho and E)."""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import style as S
P = '/raven/ptmp/jinma/idort1020'
out = json.load(open(P + '/page/rwconv.json'))
Ps = [3.33333e-3, 3.33333e-2, 0.333333, 3.33333]; Ts = [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0]
def get(p, t):
    for k, o in out.items():
        a, d, pp, tt = k.split('|')
        if a == 'PROD19' and d == 'x1' and abs(float(pp[1:]) / p - 1) < 1e-3 and abs(float(tt[1:]) / t - 1) < 1e-6: return o
fig, axs = plt.subplots(len(Ps), len(Ts), figsize=(4.0 * len(Ts), 3.4 * len(Ps)), sharex=True)
M = {v: np.full((len(Ps), len(Ts)), np.nan) for v in ('rho', 'E')}
for i, p in enumerate(Ps):
    for j, t in enumerate(Ts):
        ax = axs[i, j]; o = get(p, t)
        if o:
            for v, col, m in (('rho', S.C1, 'o'), ('E', S.C4, '^')):
                ax.loglog(o['N'], o['l1_' + v], marker=m, ms=6, color=col, label='%s p = %.2f' % (v, o[v][1])); M[v][i, j] = o[v][1]
            S.slope2(ax); ax.legend(fontsize=9, loc='lower left')
        else: ax.text(0.5, 0.5, 'pending', ha='center', transform=ax.transAxes)
        ax.set_xticks([16, 32, 64, 128, 256]); ax.set_xticklabels(['16', '32', '64', '128', '256']); ax.minorticks_off()
        if i == 0: ax.set_title('tau_lambda = %g' % t)
        if j == 0: ax.set_ylabel('P_rad/P_gas = %.3g\nL1 / amplitude' % p)
        if i == len(Ps) - 1: ax.set_xlabel('N per wavelength (dt ~ dx)')
fig.suptitle('Radiation-modified acoustic waves, production method (Cartesian vet_sc, hesdirk2), P_rad/P_gas 3.3e-3 - 3.33 x tau 0.01 - 1000 (33.3 omitted: E/(rho c^2) = 1 with c_light = 10, outside the Newtonian regime); dashed slope 2', fontsize=16)
fig.tight_layout(); fig.savefig(P + '/page/fig/conv_rw_grid.png', dpi=110)
fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
for ax, v in zip(axs, ('rho', 'E')):
    im = ax.imshow(M[v], cmap='RdYlGn', vmin=0.5, vmax=2.5, origin='lower', aspect='auto')
    for i in range(len(Ps)):
        for j in range(len(Ts)):
            if np.isfinite(M[v][i, j]): ax.text(j, i, '%.2f' % M[v][i, j], ha='center', va='center', fontsize=11)
    ax.set_xticks(range(len(Ts))); ax.set_xticklabels(['%g' % t for t in Ts]); ax.set_yticks(range(len(Ps))); ax.set_yticklabels(['%.3g' % p for p in Ps])
    ax.set_xlabel('tau_lambda'); ax.set_ylabel('P_rad/P_gas'); ax.set_title('finest-pair order (N 128 -> 256), %s' % v)
fig.colorbar(im, ax=axs, shrink=0.8, label='order')
fig.savefig(P + '/page/fig/conv_rw_ordermap.png', dpi=130, bbox_inches='tight')
print({v: np.round(M[v], 2).tolist() for v in M})
