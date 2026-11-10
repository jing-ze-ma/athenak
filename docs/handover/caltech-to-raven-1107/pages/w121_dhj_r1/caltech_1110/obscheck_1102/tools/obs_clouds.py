"""obs_clouds.py: post-processed cloud grid (em3d clouds=, emcl_*.npz) -> band metrics per case,
page/obs_clouds_w121.png, obs_clouds_w76.png, obs_clouds_summary.png, page/obs_clouds.json."""
import glob, json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/obscheck_1102/tools')
import bands  # noqa: E402

O = '/resnick/groups/carnegie_poc/jingze/obscheck_1102/'
PG = O + 'page/'
MOD, D1 = '#2a78d6', '#1b1b1b'
EXTC = {'0.43': '#1baf7a', '1.3': '#eb6834', 'inf': '#6250d6'}
SETM = {'A': '^', 'C': 'o', 'B': 's'}
NOTE = 'post-processing only: clouds added to the frozen hydro T(p); no radiative feedback'
# observed (day, night, offset east) per band: W121 Mikal-Evans+23 NIRSpec, Morello+23 Spitzer; W76 May & Stevenson 21
OBSV = {'w121': {'G395H_NRS2': ((4924, 9), (630, 10), (2.66, 0.12), 'NIRSpec NRS2 (Mikal-Evans+23)'),
                 'G395H_NRS1': ((3924, 7), (136, 8), (3.36, 0.11), 'NIRSpec NRS1 (Mikal-Evans+23)'),
                 'IRAC2': ((5080, 90), (720, 270), (-5.0, 3.2), 'Spitzer 4.5 um (Morello+23)')},
        'w76': {'IRAC2': ((3729, 52), (801, 89), (0.67, 0.2), 'Spitzer 4.5 um (May & Stevenson 21)'),
                'IRAC1': ((2539, 30), (931, 96), (0.68, 0.48), 'Spitzer 3.6 um (May & Stevenson 21)')}}
plt.rcParams.update({'font.size': 9.5, 'axes.titlesize': 10, 'axes.grid': True, 'grid.color': '#e3e3e3',
                     'grid.linewidth': 0.6, 'axes.edgecolor': '#888', 'legend.fontsize': 8, 'legend.frameon': False})


def fourier(ph, c):
    A = lambda p: np.stack([np.ones_like(p), np.cos(2*np.pi*p), np.sin(2*np.pi*p), np.cos(4*np.pi*p),   # noqa: E731
                            np.sin(4*np.pi*p)], 1)
    co = np.linalg.lstsq(A(ph), c, rcond=None)[0]
    pf = np.linspace(0, 1, 2001)
    cf = A(pf)@co
    return dict(day=np.interp(0.5, pf, cf)*1e6, night=cf[0]*1e6, amp=(cf.max() - cf.min())*1e6,
                off=(0.5 - pf[np.argmax(cf)])*360, curve=(pf, cf*1e6))


def metrics(key):
    S = bands.SYS[key]
    F = sorted(glob.glob(O + 'out/emcl_%s_0*.npz' % key))
    if not F:
        return None
    E = [np.load(f) for f in F]
    cases = json.loads(str(E[0]['cases']))
    lam, ph = E[0]['lam_um'], E[0]['phases']
    tab = np.load(str(E[0]['table']))
    edges = tab['lam_edges_um']
    Fs = bands.star_flux(S['Teff'], lam, edges)
    R = np.array([e['Fp'] for e in E])/(Fs*(S['Rs']*bands.RSUN)**2)        # (dump, case, phase, w)
    dl = edges[:, 1] - edges[:, 0]
    out = {'cases': cases, 'files': F, 'frac': json.loads(str(E[0]['cell_frac'])), 'band': {}}
    clear = np.load(O + 'out/bands_%s.npz' % key, allow_pickle=True)
    for bn in OBSV[key]:
        a, b = bands.BANDS[bn]
        w = ((lam >= a) & (lam < b))*lam*dl
        c = (R*Fs*w).sum(-1)/(Fs*w).sum()                                   # (dump, case, phase)
        res = []
        for k in range(len(cases)):
            per = [fourier(ph, c[d, k]) for d in range(c.shape[0])]
            mm = {q: float(np.mean([p[q] for p in per])) for q in ('day', 'night', 'amp', 'off')}
            mm.update({q + '_sd': float(np.std([p[q] for p in per])) for q in ('day', 'night', 'amp', 'off')})
            mm['curve'] = fourier(ph, c[:, k].mean(0))['curve']
            res.append(mm)
        st = clear['st_' + bn].item()
        cl = dict(day=st['day']*1e6, night=st['night']*1e6, amp=st['amp']*1e6, off=st['off_max_deg'],
                  curve=(np.r_[clear['phases'], 1.0], np.r_[clear['pc_' + bn], clear['pc_' + bn][0]]*1e6))
        (od, ode), (on, one), (oo, ooe), lab = OBSV[key][bn]
        oa = od - on
        for r in res:
            r['score'] = abs(r['night'] - on)/max(on, 100) + abs(r['amp'] - oa)/oa + abs(r['off'] - oo)/10.0
        out['band'][bn] = dict(res=res, clear=cl, obs=dict(day=(od, ode), night=(on, one), amp=(oa, np.hypot(ode, one)),
                                                            off=(oo, ooe)), label=lab)
    return out


def fig_planet(key, M, name):
    bn = list(OBSV[key])[0]
    B = M['band'][bn]
    cases = M['cases']
    fig, axs = plt.subplots(1, 4, figsize=(14, 4.2))
    for ax, q, t in zip(axs, ('night', 'amp', 'off', 'day'), ('nightside Fp/Fs [ppm]', 'peak-to-trough amplitude [ppm]',
                                                            'hotspot offset [deg east]', 'dayside Fp/Fs [ppm]')):
        v, e = B['obs'][q]
        ax.axhspan(v - e, v + e, color='#999', alpha=0.35, lw=0)
        ax.axhline(v, color=D1, lw=1.5, label='observed')
        ax.axhline(B['clear'][q], color=MOD, lw=1.5, ls='--', label='model, no clouds')
        for k, c in enumerate(cases):
            ex = str(c['ext'])
            x = np.log10(c['R_um']) + {'A': -0.08, 'C': 0.0, 'B': 0.08}[c['set']] + (0.15 if c['mode'] == 'abs' else 0) \
                + (0.25 if c.get('f', 1.0) != 1.0 else 0)
            r = B['res'][k]
            ax.errorbar([x], [r[q]], yerr=[r[q + '_sd']], fmt=SETM[c['set']], color=EXTC[ex], ms=6,
                        mfc=EXTC[ex] if c['mode'] == 'ext' else 'white', mec=EXTC[ex], mew=1.4, elinewidth=1)
        ax.set_xticks([-1, 0, np.log10(5)])
        ax.set_xticklabels(['0.1', '1', '5'])
        ax.set_xlabel('particle radius [um]')
        ax.set_title(t)
    h = [plt.Line2D([], [], color=EXTC[e], marker='o', ls='', label='extent %s' % ({'0.43': '1 H (0.43 dex)', '1.3': '3 H (1.3 dex)', 'inf': 'to model top'}[e]))
         for e in EXTC]
    h += [plt.Line2D([], [], color='#555', marker=SETM[s], ls='', label={'A': 'Al2O3+CaTiO3', 'C': '+ silicates', 'B': '+ silicates + Fe'}[s]) for s in SETM]
    h += [plt.Line2D([], [], color='#555', marker='o', mfc='white', ls='', label='absorption only (open, +0.15); f=0.1 mass at +0.25'),
          plt.Line2D([], [], color=D1, label='observed +-1 sigma'), plt.Line2D([], [], color=MOD, ls='--', label='model, no clouds')]
    fig.legend(handles=h, loc='lower center', ncol=5, fontsize=8, frameon=False)
    fig.suptitle('%s cloud-deck test, %s; filled = extinction (abs + (1-g) scat) as absorption. %s' % (name, B['label'], NOTE), y=1.0, fontsize=10)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(PG + 'obs_clouds_%s.png' % key, dpi=90)
    plt.close(fig)


def fig_summary(MS):
    fig, axs = plt.subplots(1, len(MS), figsize=(6.2*len(MS), 4.2))
    best = {}
    for ax, (key, name, M) in zip(np.atleast_1d(axs), MS):
        bn = list(OBSV[key])[0]
        B = M['band'][bn]
        k = int(np.argmin([r['score'] for r in B['res']]))
        best[key] = (k, M['cases'][k], B['res'][k], B['clear'], B['obs'])
        pf, cf = B['clear']['curve']
        ax.plot(pf, cf, color=MOD, ls='--', lw=1.5, label='model, no clouds')
        pf, cf = B['res'][k]['curve']
        ax.plot(pf, cf, color='#eb6834', lw=2, label='best cloud case: %s' % M['cases'][k]['name'])
        o = B['obs']
        ax.errorbar([0.5 - o['off'][0]/360], [o['day'][0]], yerr=[o['day'][1]], fmt='o', color=D1, label='observed day (at peak phase)')
        ax.errorbar([0.0, 1.0], [o['night'][0]]*2, yerr=[o['night'][1]]*2, fmt='s', color=D1, label='observed night')
        if key == 'w121':
            es = [e for e in json.load(open(O + 'obsdata.json'))['WASP-121b']['emission'] if e['id'].startswith('W121_NIRSpec_G395H_EvansSoma')][0]
            p = es['band_avg_phase_curve']['points']
            ax.errorbar([q['phase_center'] for q in p], [q['NRS2_ppm']['v'] for q in p], yerr=[q['NRS2_ppm']['e'] for q in p],
                        fmt='.', color='#777', ms=4, label='Evans-Soma+25 NRS2')
        ax.set_xlim(0, 1)
        ax.set_ylim(bottom=0)
        ax.set_xlabel('orbital phase (0.5 = eclipse)')
        ax.set_ylabel('Fp/Fs [ppm]')
        ax.set_title('%s, %s' % (name, B['label']))
        ax.legend(loc='upper left', fontsize=7.5)
    fig.suptitle('Best-matching cloud deck vs clear model vs data (%s)' % NOTE, y=1.0, fontsize=10)
    fig.tight_layout()
    fig.savefig(PG + 'obs_clouds_summary.png', dpi=90)
    plt.close(fig)
    return best


def main():
    MS = []
    J = {}
    for key, name in (('w121', 'WASP-121b'), ('w76', 'WASP-76b')):
        M = metrics(key)
        if M is None:
            continue
        fig_planet(key, M, name)
        MS.append((key, name, M))
        J[key] = {bn: dict(clear={q: B['clear'][q] for q in ('day', 'night', 'amp', 'off')}, obs=B['obs'],
                           cases=[dict(c, **{q: r[q] for q in ('day', 'night', 'amp', 'off', 'score')}) for c, r in zip(M['cases'], B['res'])])
                  for bn, B in M['band'].items()}
        J[key]['frac'] = M['frac']
    if MS:
        best = fig_summary(MS)
        J['best'] = {k: dict(case=v[1], night=v[2]['night'], amp=v[2]['amp'], off=v[2]['off'], day=v[2]['day']) for k, v in best.items()}
    json.dump(J, open(PG + 'obs_clouds.json', 'w'), indent=1, default=float)
    for k, v in J.get('best', {}).items():
        print(k, 'best', v)


if __name__ == '__main__':
    main()
