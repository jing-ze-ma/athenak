"""obs_plots.py: model (em3d + bands.py, fe_hires + tr_fe.py) vs observations (obsdata.json) ->
page/obs_*.png + page/obs_takeaways.json (one plain line per plot) + page/obs_numbers.txt."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

O = '/resnick/groups/carnegie_poc/jingze/obscheck_1102/'
PG = O + 'page/'
MOD, D1, D2, D3, D4 = '#2a78d6', '#1b1b1b', '#eb6834', '#1baf7a', '#eda100'
plt.rcParams.update({'font.size': 9.5, 'axes.titlesize': 10, 'figure.facecolor': 'white', 'axes.facecolor': 'white',
                     'savefig.facecolor': 'white', 'axes.grid': True, 'grid.color': '#e3e3e3', 'grid.linewidth': 0.6,
                     'axes.edgecolor': '#888', 'legend.fontsize': 8, 'legend.frameon': False})
OBS = json.load(open(O + 'obsdata.json'))
TK, NUM = {}, []


def ent(pl, sec, idp):
    return [e for e in OBS[pl][sec] if e['id'].startswith(idp)][0]


def band(e, name):
    return [b for b in e['bands'] if b['name'].startswith(name)][0]


def ve(x):
    if x is None:
        return np.nan, np.nan
    if 'e' in x:
        return x['v'], x['e']
    return x['v'], 0.5*(x.get('ep', 0) + x.get('em', 0))


def ebar(ax, x, y, e, c=D1, m='o', lab=None, ms=5, **k):
    ax.errorbar(x, y, yerr=e, fmt=m, color=c, ms=ms, mfc=c, mec='white', mew=0.6, elinewidth=1, capsize=0, label=lab,
                zorder=5, **k)


def model_pc(ax, B, bn, lab='model (3-D RT)'):
    ph = np.r_[B['phases'], 1.0]
    y = np.r_[B['pc_' + bn], B['pc_' + bn][0]]*1e6
    s = np.r_[B['pcsd_' + bn], B['pcsd_' + bn][0]]*1e6
    ax.fill_between(ph, y - s, y + s, color=MOD, alpha=0.2, lw=0)
    ax.plot(ph, y, color=MOD, lw=2, label=lab)
    st = B['st_' + bn].item()
    return st


def binspec(wl, f, e, nb):
    wl, f, e = map(np.asarray, (wl, f, e))
    edges = np.linspace(wl.min(), wl.max() + 1e-9, nb + 1)
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (wl >= a) & (wl < b) & np.isfinite(f)
        if m.sum():
            w = 1/e[m]**2
            out.append(((wl[m]*w).sum()/w.sum(), (f[m]*w).sum()/w.sum(), 1/np.sqrt(w.sum())))
    return np.array(out).T


def mspec(B, ph, R=None):
    i = np.argmin(abs(B['phases'] - ph))
    return B['lam'], B['ratio'][i]*1e6


# ------------------------------------------------------------------ WASP-121b phase curves
def w121_pc(B):
    pl = 'WASP-121b'
    fig, axs = plt.subplots(2, 3, figsize=(12, 6.6))
    es = ent(pl, 'emission', 'W121_NIRSpec_G395H_EvansSoma')
    me = ent(pl, 'emission', 'W121_NIRSpec_G395H_MikalEvans')
    pts = es['band_avg_phase_curve']['points']
    rows = []
    for k, (ax, det, bn) in enumerate(((axs[0, 0], 'NRS1', 'G395H_NRS1'), (axs[0, 1], 'NRS2', 'G395H_NRS2'))):
        st = model_pc(ax, B, bn)
        x = [p['phase_center'] for p in pts]
        y = [p[det + '_ppm']['v'] for p in pts]
        e = [p[det + '_ppm']['e'] for p in pts]
        ebar(ax, x, y, e, lab='Evans-Soma+25 (band mean)', ms=4)
        b = band(me, det)
        ax.set_title('JWST NIRSpec G395H %s (%.2f-%.2f um)' % (det, *b['band_um']))
        rows.append((det, st, ve(b['dayside_ppm']), ve(b['nightside_ppm']), ve(b['hotspot_offset_deg_east']),
                     ve(b['Tday_K']), ve(b['Tnight_K']), 'Mikal-Evans+23'))
    mo = ent(pl, 'emission', 'W121_Spitzer_Morello')
    dv = ent(pl, 'emission', 'W121_Spitzer_Davenport')
    for ax, bn in ((axs[0, 2], 'IRAC1'), (axs[1, 0], 'IRAC2')):
        st = model_pc(ax, B, bn)
        b = band(mo, bn)
        ebar(ax, [0.5, 0.0, 1.0], [b['F_phase0.5_ppm']['v'], b['F_phase0_ppm']['v'], b['F_phase0_ppm']['v']],
             [b['F_phase0.5_ppm']['e'], b['F_phase0_ppm']['e'], b['F_phase0_ppm']['e']], lab='Morello+23 (phase 0.5 / 0)')
        b2 = band(dv, bn)
        ebar(ax, [0.52, 0.02, 1.02], [b2['dayside_ppm']['v'], b2['nightside_ppm']['v'], b2['nightside_ppm']['v']],
             [b2['dayside_ppm']['e'], b2['nightside_ppm']['e'], b2['nightside_ppm']['e']], c=D2, m='s',
             lab='Davenport+25 (day / night)')
        ax.set_title('Spitzer IRAC %s um' % ('3.6' if bn == 'IRAC1' else '4.5'))
        rows.append((bn, st, ve(b['F_phase0.5_ppm']), ve(b['F_phase0_ppm']), ve(b['hotspot_offset_deg_east']),
                     ve(b['Tday_K']), ve(b['Tnight_K']), 'Morello+23'))
    ax = axs[1, 1]
    st = model_pc(ax, B, 'WFC3')
    h = ent(pl, 'emission', 'W121_HST_WFC3_G141_MikalEvans2022')
    p = [q for q in h['broadband_phase_curve_binned']['points'] if not q['in_or_near_transit_eclipse']]
    ebar(ax, [q['phase'] for q in p], [q['fp_ppm'] for q in p], [q['err_ppm'] for q in p], ms=3.5,
         lab='Mikal-Evans+22 (binned)')
    ax.set_title('HST WFC3 G141 (1.12-1.64 um)')
    hb = band(h, 'G141')
    e20 = band(ent(pl, 'emission', 'W121_HST_WFC3_G141_eclipse_MikalEvans2020'), 'G141')
    rows.append(('WFC3', st, ve(e20['dayside_ppm']), (np.nan, np.nan), ve(hb['hotspot_offset_deg_east']),
                 (np.nan, np.nan), (np.nan, np.nan), 'Mikal-Evans+20/22'))
    ax = axs[1, 2]
    st = model_pc(ax, B, 'SOSS')
    sp = band(ent(pl, 'emission', 'W121_NIRISS_SOSS'), 'Order 1')
    ebar(ax, [0.5], [sp['dayside_ppm']['v']], [sp['dayside_ppm']['ep']], lab='Splinter+25 (dayside)')
    ax.set_title('JWST NIRISS SOSS order 1 (0.85-2.85 um)')
    rows.append(('SOSS', st, ve(sp['dayside_ppm']), (np.nan, np.nan), ve(sp['hotspot_offset_deg_east']),
                 (2717, 17), (1562, 18), 'Splinter+25 (T = bolometric)'))
    for ax in axs.flat:
        ax.set_xlim(0, 1)
        ax.set_ylim(bottom=min(0, ax.get_ylim()[0]))
        ax.set_xlabel('orbital phase (0.5 = eclipse)')
        ax.set_ylabel('Fp/Fs [ppm]')
        ax.legend(loc='upper left')
    fig.suptitle('WASP-121b phase curves: hydro p1 (rot 298-300) through 3-D emission RT vs observations', y=0.995)
    fig.tight_layout()
    fig.savefig(PG + 'obs_w121_phasecurves.png', dpi=90)
    plt.close(fig)
    return rows


def w121_spec(B):
    pl = 'WASP-121b'
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
    k26 = ent(pl, 'emission', 'W121_JWST_panchromatic_Kahle')['dayside_spectrum']
    es = ent(pl, 'emission', 'W121_NIRSpec_G395H_EvansSoma')
    ax = axs[0]
    for ph, ls in ((0.5, '-'), (0.425, '--')):
        l, f = mspec(B, ph)
        ax.plot(l, f, color=MOD, lw=1.8 if ls == '-' else 1.2, ls=ls, label='model, phase %.3f' % ph)
    for key, nb, c, lab in (('NIRISS_order2', 6, D3, 'NIRISS SOSS o2 (Kahle+26)'), ('NIRISS_order1', 20, D3, 'NIRISS SOSS o1 (Kahle+26)'),
                            ('NIRSpec_G395H', 24, D1, 'NIRSpec G395H at eclipse (Kahle+26 from Evans-Soma+25)')):
        s = k26[key]
        x, y, e = binspec(s['wl_um'], s['fpfs_ppm'], s['err_ppm'], nb)
        ebar(ax, x, y, e, c=c, lab=lab, ms=4)
    w20 = ent(pl, 'emission', 'W121_HST_WFC3_G141_eclipse_MikalEvans2020')['dayside_spectrum']
    ebar(ax, w20['wl_um'], w20['fpfs_ppm'], 0.5*(np.array(w20['err_ppm_plus']) + np.array(w20['err_ppm_minus'])), c=D2,
         m='s', ms=3.5, lab='HST WFC3 (Mikal-Evans+20)')
    mo = ent(pl, 'emission', 'W121_Spitzer_Morello')
    for bn, x in (('IRAC1', 3.55), ('IRAC2', 4.5)):
        b = band(mo, bn)
        ebar(ax, [x], [b['F_phase0.5_ppm']['v']], [b['F_phase0.5_ppm']['e']], c=D4, m='D', ms=6,
             lab='Spitzer (Morello+23)' if bn == 'IRAC1' else None)
    ax.set_title('Dayside (eclipse) emission spectrum')
    ax = axs[1]
    l, f = mspec(B, 0.9375)
    ax.plot(l, f, color=MOD, lw=1.8, label='model, phase 0.9375')
    s = es['nightside_spectrum']
    x, y, e = binspec(s['wl_um'], s['fpfs_ppm'], np.maximum(np.array(s['err_ppm']), 1), 16)
    ebar(ax, x, y, e, lab='NIRSpec G395H phase 0.925-0.95 (Evans-Soma+25)', ms=4)
    hn = ent(pl, 'emission', 'W121_HST_WFC3_G141_MikalEvans2022')['nightside_spectrum']
    ebar(ax, hn['wl_um'], hn['fpfs_ppm'], hn['err_ppm'], c=D2, m='s', ms=3.5, lab='HST WFC3 phase 0.9-1.0 (Mikal-Evans+22)')
    for bn, x in (('IRAC1', 3.55), ('IRAC2', 4.5)):
        b = band(mo, bn)
        ebar(ax, [x], [b['F_phase0_ppm']['v']], [b['F_phase0_ppm']['e']], c=D4, m='D', ms=6,
             lab='Spitzer phase 0 (Morello+23)' if bn == 'IRAC1' else None)
    ax.set_title('Nightside emission spectrum')
    for ax in axs:
        ax.set_xlabel('wavelength [um]')
        ax.set_ylabel('Fp/Fs [ppm]')
        ax.set_xlim(0.55, 5.3)
        ax.legend(loc='upper left')
    fig.suptitle('WASP-121b emission spectra: model vs JWST / HST / Spitzer', y=1.0)
    fig.tight_layout()
    fig.savefig(PG + 'obs_w121_spectra.png', dpi=90)
    plt.close(fig)
    # spectral shape diagnostics: dayside 1.4 um water band vs continuum; NIRSpec 4.3-4.6 vs 3.5 (CO emission)
    l, f = mspec(B, 0.5)

    def mean(a, b):
        m = (l > a) & (l < b)
        return f[m].mean()
    return dict(model_h2o_14_over_12=mean(1.35, 1.45)/mean(1.2, 1.3), model_co_45_over_37=mean(4.3, 4.7)/mean(3.5, 3.8))


def w76_pc(B):
    pl = 'WASP-76b'
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.0))
    ms = ent(pl, 'emission', 'W76_Spitzer_MayStevenson')
    dg = ent(pl, 'emission', 'W76_Spitzer45_Dang')
    rows = []
    for ax, bn in ((axs[0], 'IRAC1'), (axs[1], 'IRAC2')):
        st = model_pc(ax, B, bn)
        b = band(ms, bn)
        d, de = ve(b['eclipse_depth_ppm'])
        n, ne = ve(b['nightside_ppm_derived'])
        ebar(ax, [0.5, 0.0, 1.0], [d, n, n], [de, ne, ne], lab='May & Stevenson 21 (eclipse / derived night)')
        if bn == 'IRAC2':
            b2 = band(dg, bn)
            ebar(ax, [0.52, 0.02, 1.02], [b2['eclipse_depth_ppm']['v'], b2['min_ppm']['v'], b2['min_ppm']['v']],
                 [b2['eclipse_depth_ppm']['e'], b2['min_ppm']['e'], b2['min_ppm']['e']], c=D2, m='s', lab='Dang+25 (eclipse / min)')
        ax.set_title('Spitzer IRAC %s um' % ('3.6' if bn == 'IRAC1' else '4.5'))
        ax.set_xlabel('orbital phase (0.5 = eclipse)')
        ax.set_ylabel('Fp/Fs [ppm]')
        ax.set_xlim(0, 1)
        ax.set_ylim(bottom=0)
        ax.legend(loc='upper left')
        rows.append((bn, st, (d, de), (n, ne), ve(b['hotspot_offset_deg_east']), ve(b['Tday_K']), ve(b['Tnight_K']),
                     'May & Stevenson 21'))
    ax = axs[2]
    l, f = mspec(B, 0.5)
    ax.plot(l, f, color=MOD, lw=1.8, label='model, phase 0.5')
    for idp, c, m, lab in (('W76_HST_WFC3_Fu', D1, 'o', 'HST WFC3 (Fu+21)'), ('W76_HST_WFC3_Edwards', D2, 's', 'HST WFC3 (Edwards+20)')):
        s = ent(pl, 'emission', idp)['dayside_spectrum']
        ebar(ax, s['wl_um'], s['fpfs_ppm'], s['err_ppm'], c=c, m=m, ms=3.5, lab=lab)
    for bn, x in (('IRAC1', 3.55), ('IRAC2', 4.5)):
        b = band(ms, bn)
        ebar(ax, [x], [b['eclipse_depth_ppm']['v']], [b['eclipse_depth_ppm']['e']], c=D4, m='D', ms=6,
             lab='Spitzer (May & Stevenson 21)' if bn == 'IRAC1' else None)
    ax.set_xlim(0.55, 5.3)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('wavelength [um]')
    ax.set_ylabel('Fp/Fs [ppm]')
    ax.set_title('Dayside emission spectrum')
    ax.legend(loc='upper left')
    fig.suptitle('WASP-76b: hydro p1 (rot 298-300) through 3-D emission RT vs observations (no JWST phase curve published yet)', y=1.0)
    fig.tight_layout()
    fig.savefig(PG + 'obs_w76_phasecurves.png', dpi=90)
    plt.close(fig)
    l, f = mspec(B, 0.5)
    m1, m2 = (l > 1.35) & (l < 1.45), (l > 1.2) & (l < 1.3)
    return rows, f[m1].mean()/f[m2].mean()


def summary(rows121, rows76):
    """dot plot: per band, model minus observed for day flux, night flux, peak offset, T_day, T_night."""
    fig, axs = plt.subplots(1, 4, figsize=(13, 4.4), sharey=True)
    labels, k = [], 0
    for pl, rows in (('W121', rows121), ('W76', rows76)):
        for (bn, st, day, night, off, Td, Tn, ref) in rows:
            y = -k
            labels.append('%s %s' % (pl, bn))
            for ax, mv, ov in ((axs[0], st['day']*1e6, day), (axs[1], st['night']*1e6, night),
                               (axs[2], st['off_max_deg'], off)):
                if np.isfinite(ov[0]):
                    ax.errorbar([ov[0]], [y], xerr=[ov[1]], fmt='o', color=D1, ms=5, mec='white', zorder=5)
                ax.plot([mv], [y], 'D', color=MOD, ms=6, zorder=6)
                if np.isfinite(ov[0]):
                    ax.plot([ov[0], mv], [y, y], color='#bbb', lw=1, zorder=1)
            ax = axs[3]
            for mv, ov, c in ((st['Tday'], Td, D2), (st['Tnight'], Tn, D3)):
                if np.isfinite(ov[0]):
                    ax.errorbar([ov[0]], [y], xerr=[ov[1]], fmt='o', color=c, ms=5, mec='white', zorder=5)
                    ax.plot([ov[0], mv], [y, y], color='#bbb', lw=1, zorder=1)
                ax.plot([mv], [y], 'D', color=c, ms=6, mfc='white', mew=1.5, zorder=6)
            NUM.append('%-5s %-11s model day %5.0f night %5.0f ppm off %+5.1f deg Tday %4.0f Tnight %4.0f | obs (%s) day %s night %s off %s Tday %s Tnight %s' % (
                pl, bn, st['day']*1e6, st['night']*1e6, st['off_max_deg'], st['Tday'], st['Tnight'], ref,
                '%.0f+-%.0f' % day, '%.0f+-%.0f' % night, '%.1f+-%.1f' % off, '%.0f+-%.0f' % Td, '%.0f+-%.0f' % Tn))
            k += 1
        k += 0.5
    ys = [-i for i in range(len(rows121))] + [-(len(rows121) + 0.5 + i) for i in range(len(rows76))]
    axs[0].set_yticks(ys)
    axs[0].set_yticklabels(labels)
    for ax, t in zip(axs, ('dayside Fp/Fs [ppm]', 'nightside Fp/Fs [ppm]', 'hotspot offset [deg east]',
                          'brightness T [K] (orange day, green night)')):
        ax.set_xlabel(t)
    axs[0].plot([], [], 'D', color=MOD, label='model')
    axs[0].plot([], [], 'o', color=D1, label='observed (+-1 sigma)')
    axs[0].legend(loc='lower right')
    axs[1].set_xscale('symlog', linthresh=300)
    fig.suptitle('Summary: model (diamonds) vs observed (dots), per band', y=1.0)
    fig.tight_layout()
    fig.savefig(PG + 'obs_summary.png', dpi=90)
    plt.close(fig)


def fe_tr():
    fig, axs = plt.subplots(2, 2, figsize=(12, 7.4))
    res = {}
    for r, (key, pl) in enumerate((('w76', 'WASP-76b'), ('w121', 'WASP-121b'))):
        z = np.load(O + 'out/trfe_%s.npz' % key, allow_pickle=True)
        ph, m, sd = z['phases'], z['mean'], z['sd']
        inn = abs(ph) < z['PH14'] - 1e-6
        ax = axs[r, 0]
        for k, c, ls, lab in ((0, MOD, '-', 'model, both limbs (as observed)'), (1, D3, '--', 'model, morning limb only'),
                              (2, D2, '--', 'model, evening limb only')):
            ax.plot(ph[inn], m[inn, k, 1], color=c, ls=ls, lw=2 if k == 0 else 1.3, label=lab)
        ax2 = axs[r, 1]
        ax2.plot(ph[inn], m[inn, 0, 0], color=MOD, lw=2, label='model CCF amplitude (both limbs)')
        F = {e['id']: e for e in OBS[pl]['Fe_transmission_hires']}
        if key == 'w76':
            e = F['W76_FeCCF_Ehrenreich2020']['time_series_digitised']
            ebar(ax, e['phase'], e['v_kms'], e['v_err_kms'], lab='Ehrenreich+20 (ESPRESSO, digitised)')
            ebar(ax2, e['phase'], e['depth_ppm'], e['depth_err_ppm'], lab='Ehrenreich+20 depth (digitised)')
            segs = [('W76_FeCCF_Kesseli2022', 'v_kms', D4, 'Kesseli+22'), ('W76_FeCCF_KesseliSnellen2021', 'v_kms', '#6250d6', 'Kesseli & Snellen 21'),
                    ('W76_FeRetrieval_Gandhi2022', 'dVsys_kms', '#e87ba4', 'Gandhi+22')]
            obs_v = np.array(e['v_kms']); obs_p = np.array(e['phase'])
        else:
            e = F['W121_FeCCF_Seidel2025']['time_series_digitised']
            ebar(ax, e['phase'], e['v_kms'], e['v_err_kms'], lab='Seidel+25 (ESPRESSO 4-UT, digitised)')
            segs = [('W121_FeCCF_Seidel2025', 'v_kms', D4, 'Seidel+25 segments'), ('W121_FeRetrieval_Gandhi2023', 'dVsys_kms', '#e87ba4', 'Gandhi+23 (2D)')]
            obs_v = np.array(e['v_kms']); obs_p = np.array(e['phase'])
            for s in F['W121_CCF_Borsa2021']['segments']:
                pass
        for idp, vk, c, lab in segs:
            first = True
            for s in F[idp]['segments']:
                pr = s.get('phase_range')
                if pr is None or isinstance(pr[0], list):
                    continue
                v, ev = ve(s[vk])
                ax.fill_between(pr, v - ev, v + ev, color=c, alpha=0.35, lw=0)
                ax.plot(pr, [v, v], color=c, lw=2.5, label=lab if first else None)
                first = False
            if key == 'w76' and idp == 'W76_FeCCF_Kesseli2022':
                for s in F[idp]['segments']:
                    a, ea = ve(s['amplitude_ppm'])
                    ax2.fill_between(s['phase_range'], a - ea, a + ea, color=c, alpha=0.35, lw=0)
                    ax2.plot(s['phase_range'], [a, a], color=c, lw=2.5, label='Kesseli+22 amplitude' if s is F[idp]['segments'][0] else None)
        if key == 'w121':
            for s in F['W121_CCF_Borsa2021']['segments']:
                a, ea = ve(s['contrast_ppm'])
                pr = [-z['PH14'], 0] if 'first' in s['label'] else [0, z['PH14']]
                ax2.fill_between(pr, a - ea, a + ea, color=D4, alpha=0.35, lw=0)
                ax2.plot(pr, [a, a], color=D4, lw=2.5, label='Borsa+21 contrast (stellar mask CCF)' if 'first' in s['label'] else None)
            ax2.text(0.02, 0.04, 'Seidel+25 give no Fe I depth', transform=ax2.transAxes, fontsize=8, color='#555')
        # model vs digitised series: rms difference and the egress-side mean
        mv = np.interp(obs_p, ph[inn], m[inn, 0, 1])
        res[key] = dict(rms=float(np.sqrt(np.nanmean((mv - obs_v)**2))), bias=float(np.nanmean(mv - obs_v)),
                        mod_first=float(np.nanmean(m[inn & (ph < -0.017), 0, 1])), mod_last=float(np.nanmean(m[inn & (ph > 0.017), 0, 1])),
                        obs_first=float(np.nanmean(obs_v[obs_p < -0.017])), obs_last=float(np.nanmean(obs_v[obs_p > 0.017])),
                        mod_A=float(np.nanmean(m[inn, 0, 0])), mod_A_first=float(np.nanmean(m[inn & (ph < 0), 0, 0])),
                        mod_A_last=float(np.nanmean(m[inn & (ph > 0), 0, 0])))
        for a in (ax, ax2):
            a.axvline(0, color='#999', lw=0.8)
            for x in (-z['PH14'], z['PH14'], -z['PH23'], z['PH23']):
                a.axvline(x, color='#ccc', lw=0.8, ls=':')
            a.set_xlabel('orbital phase (0 = mid-transit; dotted = contacts)')
            a.legend(loc='lower left' if a is ax else 'upper left', fontsize=7.5)
        ax.set_ylabel('Fe I line centre [km/s] (- = blueshift)')
        ax2.set_ylabel('Fe I absorption [ppm]')
        ax.set_title('%s: Fe I velocity vs phase (LSF %.2f km/s)' % (pl, float(z['LSF'])))
        ax2.set_title('%s: Fe I depth vs phase' % pl)
        ax2.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(PG + 'obs_fe_timeresolved.png', dpi=90)
    plt.close(fig)
    return res


if __name__ == '__main__':
    import os
    os.makedirs(PG, exist_ok=True)
    B1 = np.load(O + 'out/bands_w121.npz', allow_pickle=True)
    B2 = np.load(O + 'out/bands_w76.npz', allow_pickle=True)
    r121 = w121_pc(B1)
    sp = w121_spec(B1)
    r76, h2o76 = w76_pc(B2)
    summary(r121, r76)
    fe = fe_tr()
    json.dump(dict(spec121=sp, h2o76=h2o76, fe=fe), open(PG + 'obs_metrics.json', 'w'), indent=1)
    open(PG + 'obs_numbers.txt', 'w').write('\n'.join(NUM) + '\n' + json.dumps(dict(spec121=sp, h2o76=h2o76, fe=fe), indent=1) + '\n')
    import obs_clouds
    obs_clouds.main()
    print('\n'.join(NUM))
    print(json.dumps(dict(spec121=sp, h2o76=h2o76, fe=fe), indent=1))
