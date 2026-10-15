#!/usr/bin/env python3
"""bsg_fr2d_1015 compare.py: cross-run summary figure + table (reads out/<run>.json and out/<run>_prof.npz)."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.argv = ['x']
exec(open('/raven/ptmp/jinma/bsg_fr2d_1015/fr2d.py').read().split('# ---------------------------------------------------------------------------------------------- opacity')[0])
RUNL = ['P', 'repro', 'sk67', 'agcE1', 'agcL64', 'agcS']
COL = dict(P='k', repro='tab:blue', sk67='tab:cyan', agcE1='tab:red', agcL64='tab:orange', agcS='tab:purple')
fig, axs = plt.subplots(1, 2, figsize=(15, 5.8), constrained_layout=True)
rows = []
for r in RUNL:
    d = json.load(open(os.path.join(OUT, r + '.json')))
    p = np.load(os.path.join(OUT, r + '_prof.npz'))
    q = p['qph'][np.isfinite(p['qph'])].ravel()
    x2f = p['x2f']
    w = (np.ones(p['qph'].shape) * (np.cos(x2f[:-1]) - np.cos(x2f[1:]))[None, :])[np.isfinite(p['qph'])].ravel()
    axs[0].hist(q, bins=np.linspace(-0.5, 3.0, 141), weights=w / w.sum(), histtype='step', lw=1.8, color=COL[r],
                label='%s (t %.0f d): p2/50/98 %.2f/%.2f/%.2f' % (r, d['t_d'], *d['q_ph_p2_50_98']))
    nt = d['neg_frac_by_local_tau']
    tc = [np.sqrt(a * b) for a, b, n, f in nt]
    fr = [max(f, 1e-5) if f is not None else np.nan for a, b, n, f in nt]
    axs[1].loglog(tc, fr, 'o-', color=COL[r], label=r)
    R = RUNS[r]
    x1f = faces_from_file(R['gridfile']) if 'gridfile' in R else stretch_faces(
        R['grid']['nx1'], R['grid']['x1min'], R['grid']['x1max'], R['grid'].get('c', ()), R['grid']['p'])
    rph = d['rph_Rsun_mean'] * RSUN
    i = np.searchsorted(x1f, rph) - 1
    drph = x1f[i + 1] - x1f[i]
    s = d['surf']
    ncell = (s['0.1']['r_Rsun'] - s['10']['r_Rsun']) * RSUN / drph
    rows.append((r, d['t_d'], d['nx'], d['rph_Rsun_mean'], drph / RSUN, ncell, d['frac_Fph_le0_cols'], d['q_ph_p2_50_98'],
                 d['q_box3_p2_50_98'], d['L_top_over_L'], s['1']['Fconv_over_Frad'], s['1']['Fradadv'] / s['1']['Frad'],
                 d['Teff_local_p2_50_98'], d['Tph_K_p2_50_98']))
axs[0].axvline(0, color='0.6', lw=0.8)
axs[0].set_xlabel(r'$F_r r^2/(L_*/4\pi)$ at $\tau_R = 2/3$ (lab frame, per column, area-weighted)')
axs[0].set_ylabel('fraction of solid angle per bin')
axs[0].legend(fontsize=8)
axs[0].set_title('distribution of the photospheric flux (no smoothing)')
axs[1].axvspan(0.1, 10, color='0.9')
axs[1].set_xlabel(r'local $\tau_R$ (cell centre, integrated inward from the outer boundary)')
axs[1].set_ylabel(r'fraction of cells with $F_r \leq 0$ (floor $10^{-5}$ = none)')
axs[1].set_ylim(5e-6, 1)
axs[1].legend(fontsize=9)
axs[1].set_title(r'where does $F_r$ change sign? (grey: $0.1 < \tau_R < 10$)')
fig.savefig(os.path.join(OUT, 'compare_runs.png'), dpi=110)
hdr = '| run | t [d] | nx1 x nx2 x nx3 | <r_ph> [Rsun] | dr(r_ph) [Rsun] | cells tau 0.1->10 | F_r(r_ph)<=0 | q p2/50/98 | q 3x3 p2/50/98 | L_top/L | F_conv/F_rad (tau 1) | v(E+P_rr)/F_r (tau 1) | T_eff local p2/50/98 [K] | T_gas(tau 2/3) p2/50/98 [K] |'
print(hdr)
print('|' + '---|' * 14)
for r in rows:
    print('| %s | %.1f | %dx%dx%d | %.2f | %.3f | %.1f | %.2f %% | %.2f / %.2f / %.2f | %.2f / %.2f / %.2f | %.3f | %.1e | %.1e | %.0f / %.0f / %.0f | %.0f / %.0f / %.0f |' % (
        r[0], r[1], *r[2], r[3], r[4], r[5], 100 * r[6], *r[7], *r[8], r[9], r[10], r[11], *r[12], *r[13]))
