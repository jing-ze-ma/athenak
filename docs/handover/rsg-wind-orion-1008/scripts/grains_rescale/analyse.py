"""Wind Mdot / v_inf integrated over the FT launch distribution. Each v0 run represents the
bin [v0 - 1.25, v0 + 1.25) km/s (last bin to infinity); bin mass flux = mdot0 [exp(-v_lo/vc)
- exp(-v_hi/vc)]. usage: ../venv/bin/python analyse.py > out/tables.md"""
import sys
import pickle
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa
sys.path.insert(0, '/orion/ptmp/jinma/rsg_wind_1008/grains')
import grainlib as gl  # noqa
from rescale import SCAL, OPT, ALPHA, XSEED, GGAS, TGAS, O  # noqa

PL = O.replace('out/', 'plots/')
res = pickle.load(open(O + 'grid.pkl', 'rb'))
YR = gl.YR / gl.Msun
grp = {}
for r in res:
    grp.setdefault((r['star'], r['sc'], r['var'], r['alpha'], r['xseed'], r['Ggas'], r['Tg'],
                    r['dens']), []).append(r)
S = {}
for k, rs in grp.items():
    rs = sorted(rs, key=lambda x: x['v0'])
    vc, m0, ve = rs[0]['v_con'] / 1e5, rs[0]['mdot0'], rs[0]['vesc'] / 1e5
    md = {'dust': 0., 'ball': 0.}
    vw = {'dust': 0., 'ball': 0.}
    vlist = []
    for i, r in enumerate(rs):
        lo = r['v0'] - 1.25
        hi = r['v0'] + 1.25 if i < len(rs) - 1 else np.inf
        dm = m0 * (np.exp(-lo / vc) - np.exp(-hi / vc))
        r['dm'] = dm
        if r['outcome'] == 'wind' and r['vinf'] > 0:
            c = 'dust' if r['rmax0'] < 20 else 'ball'
            md[c] += dm
            vw[c] += dm * r['vinf'] / 1e5
            if c == 'dust':
                vlist.append(r['vinf'] / 1e5)
    tot = md['dust'] + md['ball']
    S[k] = dict(md_dust=md['dust'] * YR, md_ball=md['ball'] * YR, md=tot * YR,
                v_dust=vw['dust'] / md['dust'] if md['dust'] else np.nan,
                v_tot=(vw['dust'] + vw['ball']) / tot if tot else np.nan,
                v_dust_med=np.median(vlist) if vlist else np.nan,
                ball_an=m0 * np.exp(-ve * np.sqrt(1 - k[5]) / vc) * YR, vc=vc, ve=ve)
    S[k]['obs'] = (1e-7 <= S[k]['md'] <= 1e-5) and (10 <= S[k]['v_tot'] <= 40)


def scl(sc):
    return 'f_rho=%g' % sc[1] if sc[0] == 'f' else 'vesc/vcon=%g' % sc[1]


print('# rescale scan (dens=ft unless noted)\n')
print('## Per (star, scaling): ballistic Mdot (analytic, Ggas 0 / 0.1), dust-driven Mdot and '
      'mass-weighted v_inf: median [min-max] over optics x alpha x seed x Ggas x Tg (64)\n')
print('| star | scaling | v_esc/v_con | ballistic Mdot G0 / G0.1 | dust Mdot med [range] | '
      'v_inf dust (mass-wt) med [range] | total Mdot med | v_inf total med | grains form r/R '
      '(pure / Fe1e-3) | a=0.1um at r/R (pure / Fe1e-3) |')
print('|---' * 10 + '|')
for st in gl.STARS:
    for sc in SCAL:
        ks = [k for k in S if k[0] == st and k[1] == sc and k[7] == 'ft']
        b = [S[k]['ball_an'] for k in ks if k[5] == 0][0], [S[k]['ball_an'] for k in ks
                                                             if k[5] == 0.1][0]
        md = np.array([S[k]['md_dust'] for k in ks])
        vd = np.array([S[k]['v_dust'] for k in ks])
        mt = np.array([S[k]['md'] for k in ks])
        vt = np.array([S[k]['v_tot'] for k in ks])
        pos = md[md > 0]
        rr = {}
        for var in ['pure', 'Fe1e-3']:
            q = [r for r in res if r['star'] == st and r['sc'] == sc and r['var'] == var
                 and r['dens'] == 'ft']
            rr[var] = (np.nanmedian([r.get('r_seed', np.nan) for r in q]) if any(np.isfinite(r.get('r_seed', np.nan))
                       for r in q) else np.nan,
                       np.nanmedian([r.get('r_01', np.nan) for r in q]) if any(np.isfinite(r.get('r_01', np.nan))
                       for r in q) else np.nan)
        print('| %s | %s | %.1f | %.1e / %.1e | %s | %s | %.1e | %.0f | %.2f / %.2f | '
              '%.2f / %.2f |' % (
                  st, scl(sc), S[ks[0]]['ve'] / S[ks[0]]['vc'], b[0], b[1],
                  '%.1e [%.1e-%.1e] (%d/64 >0)' % (np.median(pos), pos.min(), pos.max(),
                                                   len(pos)) if len(pos) else '0',
                  '%.0f [%.0f-%.0f]' % (np.nanmedian(vd), np.nanmin(vd), np.nanmax(vd))
                  if np.isfinite(vd).any() else '-', np.median(mt), np.nanmedian(vt),
                  rr['pure'][0], rr['Fe1e-3'][0], rr['pure'][1], rr['Fe1e-3'][1]))

print('\n## Observed-box map: cases (of 8 = seed x Ggas x Tg) with TOTAL wind Mdot in '
      '1e-7..1e-5 AND mass-weighted v_inf 10-40 km/s; [dust-only box count in brackets]\n')
hdr = [(o, a) for o in OPT for a in ALPHA]
print('| star | scaling | ' + ' | '.join('%s a=%g' % h for h in hdr) + ' |')
print('|---|---' + '|---' * len(hdr) + '|')
for st in gl.STARS:
    for sc in SCAL:
        cells = []
        for o, a in hdr:
            ks = [k for k in S if k[:4] == (st, sc, o, a) and k[7] == 'ft']
            n = sum(S[k]['obs'] for k in ks)
            nd = sum((1e-7 <= S[k]['md_dust'] <= 1e-5) and (10 <= S[k]['v_dust'] <= 40)
                     for k in ks)
            cells.append('%d [%d]' % (n, nd))
        print('| %s | %s | ' % (st, scl(sc)) + ' | '.join(cells) + ' |')

print('\n## mc density check (pure): total Mdot median / v_inf median, ft vs mc\n')
for st in gl.STARS:
    line = []
    for sc in SCAL:
        v = []
        for d in ['ft', 'mc']:
            ks = [k for k in S if k[0] == st and k[1] == sc and k[2] == 'pure' and k[7] == d]
            v.append('%.1e/%.0f' % (np.median([S[k]['md'] for k in ks]),
                                    np.nanmedian([S[k]['v_tot'] for k in ks])))
        line.append('%s: %s vs %s' % (scl(sc), v[0], v[1]))
    print('- %s: %s' % (st, '; '.join(line)))

# v_inf vs v0 and Gamma_d, dust-driven runs
dw = [r for r in res if r['outcome'] == 'wind' and r['vinf'] > 0 and r['rmax0'] < 20]
print('\n## Dust-driven runs: v_inf distribution; runs with v_inf in 10-40 km/s\n')
for st in gl.STARS:
    d = [r for r in dw if r['star'] == st]
    v = np.array([r['vinf'] for r in d]) / 1e5
    lo = [r for r in d if 10 <= r['vinf'] / 1e5 <= 40]
    print('- %s: %d dust-driven runs, v_inf min %.1f, 10/50/90 %% = %.0f/%.0f/%.0f km/s; %d '
          'with 10-40 km/s' % (st, len(d), v.min(), *np.percentile(v, [10, 50, 90]), len(lo)))
    if lo:
        g = np.array([r['gdmax'] for r in lo])
        print('  those: v0 %s km/s; Gd_max med %.2f [%.2f-%.2f]; r(Gamma>1) med %.2f; v(at '
              'Gamma>1) med %.1f km/s (neg = falling); optics %s; scalings %s' % (
                  sorted(set(r['v0'] for r in lo)), np.median(g), g.min(), g.max(),
                  np.nanmedian([r['r_G1'] for r in lo]),
                  np.nanmedian([r['v_G1'] for r in lo]) / 1e5,
                  sorted(set(r['var'] for r in lo)), sorted(set(scl(r['sc']) for r in lo))))
        fb = [r for r in lo if np.isfinite(r['v_G1']) and r['v_G1'] < 0]
        print('  of which Gamma>1 reached while falling back: %d; their mass share of the '
              'wind bins: see map' % len(fb))
    # binned v_inf vs Gd_max
    g = np.array([r['gdmax'] for r in d])
    for a_, b_ in [(0, 1.5), (1.5, 3), (3, 6), (6, 100)]:
        m = (g >= a_) & (g < b_)
        if m.any():
            print('  Gd_max %g-%g: n=%d v_inf median %.0f [%.0f-%.0f]' %
                  (a_, b_, m.sum(), np.median(v[m]), v[m].min(), v[m].max()))
    vv0 = np.array([r['v0'] for r in d])
    print('  v_inf median by v0: ' + ', '.join('%g:%.0f' % (x, np.median(v[vv0 == x]))
                                               for x in sorted(set(vv0))))

# ---------- plot ----------
cols = dict(zip(OPT, ['C0', 'C2', 'C1', 'C3']))
fig, ax = plt.subplots(2, 2, figsize=(14, 8.5))
xs = np.arange(len(SCAL))
for i, st in enumerate(gl.STARS):
    for o, a in hdr:
        mm, mlo, mhi, vv, vlo, vhi = [], [], [], [], [], []
        for sc in SCAL:
            ks = [k for k in S if k[:4] == (st, sc, o, a) and k[7] == 'ft']
            m = np.array([S[k]['md'] for k in ks])
            v = np.array([S[k]['v_tot'] for k in ks])
            mm.append(np.median(m)); mlo.append(m.min()); mhi.append(m.max())
            vv.append(np.nanmedian(v) if np.isfinite(v).any() else np.nan)
            vlo.append(np.nanmin(v) if np.isfinite(v).any() else np.nan)
            vhi.append(np.nanmax(v) if np.isfinite(v).any() else np.nan)
        ls = '-' if a == 1 else '--'
        off = 0.04 * (OPT.index(o) - 1.5) + (0.02 if a == 1 else -0.02)
        ax[0, i].errorbar(xs + off, mm, yerr=[np.array(mm) - mlo, np.array(mhi) - mm],
                          c=cols[o], ls=ls, marker='o', ms=3, capsize=2, label='%s a=%g' % (o, a))
        ax[1, i].errorbar(xs + off, vv, yerr=[np.array(vv) - vlo, np.array(vhi) - vv],
                          c=cols[o], ls=ls, marker='o', ms=3, capsize=2)
    ax[0, i].set_yscale('log')
    ax[0, i].axhspan(1e-7, 1e-5, color='0.85', zorder=0)
    ax[0, i].axhspan(2e-6, 4e-6, color='0.6', zorder=0)
    ax[1, i].axhspan(10, 40, color='0.85', zorder=0)
    ax[1, i].axhspan(10, 15, color='0.6', zorder=0)
    ax[0, i].set_title('%s: total wind Mdot (median, range over seed x Ggas x Tg)' % st,
                       fontsize=9)
    ax[1, i].set_title('%s: mass-weighted v_inf; grey = observed RSG, dark = Betelgeuse' % st,
                       fontsize=9)
    ax[0, i].set_ylabel('Mdot [Msun/yr]')
    ax[1, i].set_ylabel('v_inf [km/s]')
    for a_ in ax[:, i]:
        a_.set_xticks(xs)
        a_.set_xticklabels([scl(s) for s in SCAL], rotation=30, fontsize=7)
ax[0, 0].legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(PL + 'mdot_vinf_vs_scaling.png', dpi=110)

# v_inf vs Gd_max scatter
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
for i, st in enumerate(gl.STARS):
    d = [r for r in dw if r['star'] == st]
    ax[i].scatter([r['gdmax'] for r in d], [r['vinf'] / 1e5 for r in d],
                  c=[r['v0'] for r in d], s=6, cmap='viridis')
    ax[i].axhspan(10, 40, color='0.85', zorder=0)
    ax[i].set_xscale('log')
    ax[i].set_xlabel('max Gamma_dust')
    ax[i].set_ylabel('v_inf [km/s] (colour: v0)')
    ax[i].set_title('%s dust-driven runs, all scalings' % st, fontsize=9)
fig.tight_layout()
fig.savefig(PL + 'vinf_vs_gammad.png', dpi=110)
