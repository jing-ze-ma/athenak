"""compare57.py PAGEDIR: the 57.4 d comparison of the BSG reproduction run with Ma+2026
(10-08).  Reads (read only) the merged run view's user.hst, the bsg_figs dump reductions
cache/repro_sf_4n/red_NNNNN.npz, cache/repro_sf_4n/numbers.json (Fig 6 statistics) and the
per-dump comparison reductions cache/cmp57/c_NNNNN.npz (cmp57_dumps.py).  Writes
cache/cmp57.json, out/ours_cmp57.png (paper Figs 5 and 7 analogues + time-averaged L(r))
and adds the keys cmp-body / w-cmp / n-cmp to PAGEDIR/status.json (after gen_status.py).
Averaging window: the paper's steady window, 30 d to the end (57.4 d), all dumps in it
(every ~0.5 d) and all hst rows (every 1000 s)."""
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

S = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, S)
import bsg_figs as B  # noqa: E402

RUN = os.environ.get('BSG_REPRO_RUN', '/viper/ptmp2/jinma/bsg_viper_1006/merged_repro')
RSUN, DAY, YR, MSUN = B.RSUN, B.DAY, 3.15576e7, 1.989e33
C_LIGHT, SIGMA = 2.99792458e10, 5.670374419e-5
L_IC = B.L_STAR
T0, T1 = 30.0, 57.41
RS = '&#9737;'


def rph(r, tm):
    return B.r_photo(r, tm)


def main():
    page = sys.argv[1]
    hst = B.read_hst(os.path.join(RUN, 'bsg3d.user.hst'))
    th = hst['time'] / DAY
    om = np.pi / 3
    ws = 4 * np.pi / om
    hw = (th >= T0) & (th <= T1)
    # --- bsg_figs reductions of every dump (time series)
    reds = [dict(np.load(f)) for f in sorted(glob.glob(os.path.join(
        B.cache_dir(RUN), 'red_*.npz')))]
    reds = [d for d in reds if 'tau_mean' in d]
    td = np.array([float(d['time']) / DAY for d in reds])
    r = reds[0]['r']
    rr = r / RSUN
    rph3 = np.array([rph(r, d['tau_mean']) / RSUN for d in reds])
    dw = (td >= T0 - 1e-3) & (td <= T1)
    W = [d for d, k in zip(reds, dw) if k]
    i40 = int(np.argmin(np.abs(td - 40.0)))
    # --- comparison reductions
    cmp = {int(os.path.basename(f)[2:7]): dict(np.load(f)) for f in
           glob.glob(os.path.join(S, 'cache', 'cmp57', 'c_*.npz'))}
    CW = [cmp[k] for k in sorted(cmp) if T0 - 1e-3 <= float(cmp[k]['time']) / DAY <= T1]
    tcw = np.array([float(c['time']) / DAY for c in CW])
    c40 = CW[int(np.argmin(np.abs(tcw - 40.0)))]
    cend = CW[-1]
    res = dict(window_d=[float(tcw.min()), float(tcw.max())], n_dumps=len(CW),
               n_hst=int(hw.sum()), t40_d=float(c40['time']) / DAY,
               t_end_d=float(cend['time']) / DAY)
    # 1. luminosity
    Ltop = hst['L_top'][hw] * ws
    res['L_top_mean'] = float(Ltop.mean())
    res['L_top_over_LIC'] = float(Ltop.mean() / L_IC)
    Lrad = np.mean([c['L_rad'] for c in CW], 0)
    Leg = np.mean([c['L_enth_gas'] for c in CW], 0)
    Ler = np.mean([c['L_enth_rad'] for c in CW], 0)
    Lk = np.mean([c['L_kin'] for c in CW], 0)
    Ltot = Lrad + Leg + Lk                  # lab radiative flux already holds 4/3 E v
    Lconv = Leg + Ler + Lk
    rph_w = float(np.nanmean(rph3[dw]))
    ib = (rr >= 30) & (rr <= rph_w - 1.0)
    for k in ('L_in', 'L_bot', 'L_mid', 'L_top'):
        res['hst_' + k] = float(hst[k][hw].mean() * ws / L_IC)
    res['hst_L_top_rms'] = float(hst['L_top'][hw].std() * ws / L_IC)
    res['r_mid'] = float(rr[len(rr) // 2])
    res['dEtot_dt_over_L'] = float(np.polyfit(hst['time'][hw], hst['Etot'][hw], 1)[0] * ws
                                   / L_IC)
    res['Ltot_over_LIC_25_rph'] = [float((Ltot[ib] / L_IC).min()),
                                   float((Ltot[ib] / L_IC).max())]
    fc = Lconv / Ltot
    j = np.where(ib)[0][np.argmax(fc[ib])]
    res['Fconv_frac_max'] = float(fc[j])
    res['Fconv_frac_max_r'] = float(rr[j])
    # 2. photosphere, Teff
    res['rph3D_mean'] = rph_w
    res['rph3D_range'] = [float(np.nanmin(rph3[dw])), float(np.nanmax(rph3[dw]))]
    res['rph3D_40'] = float(rph3[i40])
    res['rph3D_end'] = float(rph3[-1])
    res['rph3D_30'] = float(rph3[np.argmin(np.abs(td - 30.0))])
    ro40 = c40['r_own'].ravel() / RSUN
    roW = np.concatenate([c['r_own'].ravel() for c in CW]) / RSUN
    res['r_own_95_40'] = [float(x) for x in np.percentile(ro40, [2.5, 97.5])]
    res['r_own_95_W'] = [float(x) for x in np.percentile(roW, [2.5, 97.5])]
    res['Teff'] = float((Ltop.mean() / (4 * np.pi * SIGMA * (rph_w * RSUN) ** 2)) ** 0.25)
    res['Teff_paper_implied'] = float((L_IC / (4 * np.pi * SIGMA * (43 * RSUN) ** 2)) ** 0.25)
    # 3. opacity peaks, ram pressure, tau = c / v
    kap = np.mean([d['kap_sh'] for d in W], 0)
    m = (rr > 30) & (rr < 46)
    jfe = np.where(m)[0][np.argmax(kap[m])]
    res['fe_peak_r'] = float(rr[jfe])
    res['fe_peak_kappa'] = float(kap[jfe])
    res['kappa_edd'] = float(B.KAPPA_EDD)
    m2 = (rr > rr[jfe] + 3) & (rr < rph_w)
    if m2.any():
        jhe = np.where(m2)[0][np.argmax(kap[m2])]
        res['he_bump_r'] = float(rr[jhe])
        res['he_bump_kappa'] = float(kap[jhe])
    ram = []
    for d in W:
        pr = d['p0_d11'] * d['p0_erad'] if 'p0_d11' in d else d['p0_erad'] / 3
        ram.append(np.mean(0.5 * d['p0_rho'] * d['p0_v2'] / (d['p0_pg'] + pr), 0))
    ram = np.mean(ram, 0)
    res['ram_max_below_rph'] = float(ram[rr < rph_w].max())
    for thr in (0.05, 0.1):
        k = np.where((ram >= thr) & (rr < rph_w))[0]
        res['ram_r_%g' % thr] = float(rr[k[0]]) if len(k) else None
    rcv = []
    for c in CW:
        x = c['tau_mean'] * c['vrms_sh'] / C_LIGHT     # = 1 where tau = c / v
        k = np.where((x[:-1] >= 1) & (x[1:] < 1))[0]
        rcv.append(rr[k[-1]] if len(k) else np.nan)
    res['r_tau_cv_mean'] = float(np.nanmean(rcv))
    res['r_tau_cv_range'] = [float(np.nanmin(rcv)), float(np.nanmax(rcv))]
    # 4. velocities
    vW = np.concatenate([c['v_own'].ravel() for c in CW])
    v40 = c40['v_own'].ravel()
    res['v_own_W'] = [float(x) for x in np.percentile(vW, [2.5, 10, 50, 90, 97.5])]
    res['v_own_40'] = [float(x) for x in np.percentile(v40, [2.5, 10, 50, 90, 97.5])]
    j36 = int(np.argmin(np.abs(rr - 36.0)))
    res['vr_rms_36'] = float(np.mean([c['vr_rms_sh'][j36] for c in CW]) / 1e5)
    res['vr_rms_fe'] = float(np.mean([c['vr_rms_sh'][jfe] for c in CW]) / 1e5)
    s36 = [d for d in W]
    res['vr_shell406_40'] = [float(reds[i40]['s36_vr'].min() / 1e5),
                             float(reds[i40]['s36_vr'].max() / 1e5)]
    res['vr_shell406_W'] = [float(min(d['s36_vr'].min() for d in s36) / 1e5),
                            float(max(d['s36_vr'].max() for d in s36) / 1e5)]
    # 5. energies, mass flow
    ke = hst['KE_int'][hw]
    res['KE_mean'] = float(ke.mean())
    res['KE_range'] = [float(ke.min()), float(ke.max())]
    Eg0, U0 = float(reds[0]['E_grav']), float(reds[0]['U_gas'] + reds[0]['U_rad'])
    for tt in (25.0, T1):
        k = int(np.argmin(np.abs(td - tt)))
        res['dEgrav_%g' % tt] = float(reds[k]['E_grav']) - Eg0
        res['dU_%g' % tt] = float(reds[k]['U_gas'] + reds[k]['U_rad']) - U0
    menv = float(reds[0]['M_wedge']) * ws
    res['M_env_Msun'] = menv / MSUN
    inner = rr < 40
    f10 = [np.abs(d['mdot_r'][inner]).max() / menv * YR for d, t in zip(reds, td) if t <= 10]
    fW = [np.abs(d['mdot_r'][inner]).max() / menv * YR for d in W]
    res['mdot_env_max_0_10'] = float(max(f10))
    res['mdot_env_max_W'] = float(max(fW))
    res['mdot_env_median_W'] = float(np.median(fW))
    res['Mdot_top_Msun_yr'] = float(hst['Mdot_top'][hw].mean() * ws / MSUN * YR)
    res['M_tot_change_frac_W'] = float((hst['M_tot'][hw][-1] - hst['M_tot'][hw][0]) /
                                       hst['M_tot'][hw][0])
    # 6. density structure
    rho = np.mean([c['rho_sh'] for c in CW], 0)
    clump = np.mean([c['rho2_sh'] / c['rho_sh'] ** 2 for c in CW], 0)
    mi = (rr > 30) & (rr < rph_w)
    dl = np.diff(np.log(rho))
    inv = np.where(mi[1:] & (dl > 0))[0]
    res['rho_inversion_r'] = ([float(rr[inv[0] + 1]), float(rr[inv[-1] + 1])]
                              if len(inv) else None)
    res['rho_inversion_frac_dumps'] = float(np.mean(
        [np.any(np.diff(np.log(c['rho_sh']))[mi[1:]] > 0) for c in CW]))
    res['clump_max_below_rph'] = float(clump[mi].max())
    res['clump_max_r'] = float(rr[mi][np.argmax(clump[mi])])
    jp = int(np.argmin(np.abs(rr - rph_w)))
    res['clump_at_rph'] = float(clump[jp])
    # 7. light curve
    num = json.load(open(os.path.join(B.cache_dir(RUN), 'numbers.json')))
    f6 = num['fig6']
    res['lc_window'] = f6['t_range_d']
    res['lc_flux_range'] = [f6['flux_min'], f6['flux_max']]
    res['lc_rms'] = f6['rms_eq8']
    ps = f6.get('psd_eq8', {})
    res['psd_f50'] = ps.get('f50_muHz')
    res['psd_f80'] = ps.get('f80_muHz')
    res['psd_peak_binned'] = ps.get('peak_binned_noMode_muHz')
    res['psd_peak_fP'] = ps.get('peak_fP_noMode_muHz')
    res['psd_nu_char'] = ps.get('harvey', {}).get('nu_char_muHz')
    res['psd_slope'] = ps.get('slope_10_100')
    res['mode'] = f6.get('mode')
    json.dump(res, open(os.path.join(S, 'cache', 'cmp57.json'), 'w'), indent=1)
    print(json.dumps(res, indent=1))
    # --- figure: (a) tau(r) as paper Fig 5, (b) |v| at own photosphere as Fig 7,
    # (c) time-averaged L(r) components
    fig, ax = plt.subplots(1, 3, figsize=(14.4, 4.6))
    a = ax[0]
    tp = np.mean([c['tau_p'] for c in CW], 0)
    tmw = np.mean([c['tau_mean'] for c in CW], 0)
    cv = C_LIGHT / np.mean([c['vrms_sh'] for c in CW], 0)
    a.fill_between(rr, tp[0], tp[2], color='C0', alpha=0.3, label='95 % of columns')
    a.semilogy(rr, tmw, 'C0-', lw=2, label=r'$\bar\tau(r)$, angle mean')
    a.semilogy(rr, cv, 'k:', lw=1.5, label=r'$c/\bar v$')
    a.axvline(res['fe_peak_r'], color='0.4', ls='--', lw=1, label='Fe opacity peak')
    a.plot([rph_w], [1.0], 'bo', label=r'photosphere, 3D ($\bar\tau$ = 1)')
    a.set_xlim(32, 56)
    a.set_ylim(1e-1, 1e7)
    a.set_xlabel(r'r [$R_\odot$]')
    a.set_ylabel(r'optical depth $\tau$ from r')
    a.set_title('(a) as paper Fig. 5, mean over %.0f-%.1f d' % (tcw.min(), tcw.max()),
                fontsize=10)
    a.legend(fontsize=8, loc='upper right')
    a = ax[1]
    vb = np.arange(0, 151, 1.0)
    hW = np.mean([c['v_own_hist'] for c in CW], 0)
    a.step(vb[:-1], hW, where='post', color='C0', lw=1.5,
           label='ours, mean over %.0f-%.1f d' % (tcw.min(), tcw.max()))
    a.step(vb[:-1], c40['v_own_hist'], where='post', color='C1', lw=1,
           label='ours, t = %.1f d' % (float(c40['time']) / DAY))
    a.axvline(40, color='k', ls='--', lw=1, label=r'SK -67 133 $\Theta$ = 40 km/s')
    a.axvspan(30, 40, color='0.85', alpha=0.6, lw=0, label='paper median 30-40 km/s')
    a.set_xlim(0, 100)
    a.set_xlabel('|v| at each column\'s own photosphere [km/s]')
    a.set_ylabel('relative fraction (per 1 km/s, solid-angle weighted)')
    a.set_title('(b) as paper Fig. 7', fontsize=10)
    a.legend(fontsize=8)
    a = ax[2]
    for y, lab, col in ((Ltot, 'total', 'k'), (Lrad - Ler, 'radiative diffusion', 'C0'),
                        (Lconv, 'convective (gas + rad. enthalpy + kinetic)', 'C3'),
                        (Lk, 'kinetic', 'C2')):
        a.plot(rr, y / L_IC, color=col, lw=1.5, label=lab)
    a.axhline(1, color='0.5', lw=0.8)
    a.axvline(rph_w, color='b', ls='--', lw=1)
    a.set_xlim(20, 52)
    a.set_ylim(-0.6, 1.6)
    a.set_xlabel(r'r [$R_\odot$]')
    a.set_ylabel(r'$L(r)/L_{\rm IC}$ (whole-star equivalent)')
    a.set_title('(c) time-averaged L(r), %.0f-%.1f d' % (tcw.min(), tcw.max()), fontsize=10)
    a.legend(fontsize=8, loc='lower left')
    fig.tight_layout()
    fig.savefig(os.path.join(S, 'out', 'ours_cmp57.png'), dpi=100)
    plt.close(fig)
    # --- status keys
    st_path = os.path.join(page, 'status.json')
    st = json.load(open(st_path))
    st.update(status_keys(res))
    json.dump(st, open(st_path, 'w'), indent=1)


def e(x, nd=1):
    ex = int(np.floor(np.log10(abs(x))))
    return f'{x / 10 ** ex:.{nd}f}&times;10<sup>{ex}</sup>'


def status_keys(r):
    R = f' R<sub>{RS}</sub>'
    win = f'{r["window_d"][0]:.1f}&ndash;{r["window_d"][1]:.1f} d'
    vw, v4 = r['v_own_W'], r['v_own_40']
    rows = [
        ('Luminosity', 'input log L/L<sub>&#9737;</sub> = 5.161 (L = 5.58&times;10<sup>38</sup> '
         'erg/s); 3-D L<sub>bol</sub> not quoted',
         f'mean over {win}, units of the injected L<sub>in</sub> = L<sub>IC</sub>: '
         f'L<sub>top</sub> = {r["L_top_over_LIC"]:.3f} (rms {r["hst_L_top_rms"]:.3f}; '
         f'{e(r["L_top_mean"], 2)} erg/s), L at {r["r_mid"]:.1f}{R} = {r["hst_L_mid"]:.3f}; '
         f'L(r) of the dumps {r["Ltot_over_LIC_25_rph"][0]:.2f}&ndash;'
         f'{r["Ltot_over_LIC_25_rph"][1]:.2f} from 30{R} to 1{R} below the photosphere '
         f'(panel c). The domain still gains energy at {r["dEtot_dt_over_L"]:.3f} L',
         'close (4.5 % below the input L: the envelope is not yet in thermal balance)'),
        ('3-D photosphere', '43' + R + ' at 40 d (&tau;&#772; = 1), down from the 1-D MESA ~51'
         + R + '; "slight shrinking" to 57.4 d',
         f'{r["rph3D_mean"]:.1f}{R} (mean, {win}; range {r["rph3D_range"][0]:.1f}&ndash;'
         f'{r["rph3D_range"][1]:.1f}); {r["rph3D_40"]:.1f} at 40 d, {r["rph3D_end"]:.1f} at '
         f'{r["t_end_d"]:.1f} d (IC, our 1-D integration: 53.0)', '<b>differs</b>: ours stays near the 1-D radius, '
         f'{r["rph3D_mean"] - 43:.1f}{R} above the paper\'s'),
        ('Local photosphere spread', '42&ndash;44' + R + ' (95 % of angles, 40 d), '
         '&asymp;2' + R + ' or 5 %',
         f'{r["r_own_95_40"][0]:.1f}&ndash;{r["r_own_95_40"][1]:.1f}{R} at 40 d; '
         f'{r["r_own_95_W"][0]:.1f}&ndash;{r["r_own_95_W"][1]:.1f}{R} pooled over {win} '
         f'(each column\'s own &tau;<sub>R</sub> = 1)', 'similar width (~2' + R + '), '
         'shifted outward with the photosphere'),
        ('T<sub>eff</sub>', '15835 K quoted (MESA, 1-D); A12 with 43' + R + ' and the input L '
         f'gives {r["Teff_paper_implied"]:.0f} K (not quoted); SK&nbsp;&minus;67&nbsp;133 '
         '16630 &plusmn; 850 K',
         f'{r["Teff"]:.0f} K (A12 with mean L<sub>top</sub> and mean r<sub>ph,3D</sub>)',
         '<b>differs</b> from the paper\'s 3-D value (radius); close to the 1-D model, '
         '1.4&sigma; below SK&nbsp;&minus;67&nbsp;133'),
        ('Fe opacity peak', '&asymp;35' + R + ' (3-D, 40 d); MESA ~40' + R,
         f'{r["fe_peak_r"]:.1f}{R}, shell-mean &kappa;<sub>R</sub> = '
         f'{r["fe_peak_kappa"]:.2f} cm<sup>2</sup>/g (&kappa;<sub>Edd</sub> = '
         f'{r["kappa_edd"]:.2f}), mean over {win}',
         f'<b>differs</b> by {r["fe_peak_r"] - 35:.0f}{R} (41.4{R} at t = 0 with the '
         'paper\'s own IC and tables)'),
        ('He opacity bump', '&asymp;42&ndash;43' + R + ', merged with the Fe zone',
         (f'{r["he_bump_r"]:.1f}{R} (&kappa;<sub>R</sub> = {r["he_bump_kappa"]:.2f})'
          if 'he_bump_r' in r else 'no separate maximum'), 'shifted outward with the envelope'),
        ('Ram / total pressure', 'significant above &asymp;33' + R + ' (colour scale to 0.5)',
         f'&phi; = 0 plane, &theta;-mean over {win}: max {r["ram_max_below_rph"]:.2f} below '
         f'the photosphere; &ge; 0.05 from {r["ram_r_0.05"]:.1f}{R}, &ge; 0.1 from '
         f'{r["ram_r_0.1"]:.1f}{R}',
         'convective zone sits further out, around our Fe peak'),
        ('&tau;&#772; = c/v&#772; depth', 'just above the Fe peak (&asymp;35&ndash;36' + R + ')',
         f'{r["r_tau_cv_mean"]:.1f}{R} (mean; {r["r_tau_cv_range"][0]:.1f}&ndash;'
         f'{r["r_tau_cv_range"][1]:.1f})', 'same position relative to the Fe peak'),
        ('Convective flux', 'not quantified (only the &tau;&#772; = c/v&#772; criterion)',
         f'L<sub>conv</sub>/L at most {r["Fconv_frac_max"]:.2f} (at '
         f'{r["Fconv_frac_max_r"]:.1f}{R}), gas + radiation enthalpy + kinetic flux, '
         f'mean over {win} (panel c)', 'ours only'),
        ('|v| at the local photosphere', 'broad 10&ndash;70 km/s, median 30&ndash;40 km/s '
         '(40 d); SK&nbsp;&minus;67&nbsp;133 &Theta; &asymp; 40 km/s',
         f'median {vw[2]:.0f} km/s, 10&ndash;90 % {vw[1]:.0f}&ndash;{vw[3]:.0f}, 95 % '
         f'{vw[0]:.0f}&ndash;{vw[4]:.0f} km/s pooled over {win}; at 40 d median {v4[2]:.0f}, '
         f'10&ndash;90 % {v4[1]:.0f}&ndash;{v4[3]:.0f} km/s',
         ('consistent' if 30 <= vw[2] <= 40 else
          ('slightly lower median, same broad range' if 20 <= vw[2] < 30 else
           ('<b>lower</b>' if vw[2] < 30 else '<b>higher</b>')))
         + ' (panel b)'),
        ('v<sub>r</sub> on the Fig. 4 shell', '|v<sub>r</sub>| up to ~60 km/s at 36' + R,
         f'at 40.6{R} (just below our Fe peak): {r["vr_shell406_40"][0]:+.0f} to '
         f'{r["vr_shell406_40"][1]:+.0f} km/s at 40 d; rms v<sub>r</sub> at the Fe peak '
         f'{r["vr_rms_fe"]:.0f} km/s', 'comparable; our narrow downflows reach '
         f'{-r["vr_shell406_W"][0]:.0f} km/s'),
        ('Kinetic energy', '&asymp;10<sup>43</sup> erg at 20 d, steady after ~25 d',
         f'KE<sub>int</sub> {e(r["KE_mean"])} erg (mean, {win}; range '
         f'{e(r["KE_range"][0])}&ndash;{e(r["KE_range"][1])})',
         ('consistent' if r['KE_mean'] > 3e42 else '<b>lower</b> by '
          f'{1e43 / r["KE_mean"]:.0f}&times;')),
        ('Energy release', '&Delta;E<sub>grav</sub>, &Delta;U &asymp; &minus;5&times;10<sup>44</sup> '
         'erg by 20&ndash;25 d',
         f'&Delta;E<sub>grav</sub> = {e(r["dEgrav_25"])}, &Delta;U = {e(r["dU_25"])} erg at 25 d; '
         f'{e(r["dEgrav_57.41"])} and {e(r["dU_57.41"])} erg at the end',
         '<b>differs</b>: ~10&times; smaller &Delta;E<sub>grav</sub> and &Delta;U of the '
         'opposite sign; our envelope does not release the paper\'s 5&times;10<sup>44</sup> '
         'erg'),
        ('Radial mass flow', '&#7744;<sub>r</sub>/M<sub>env</sub> &gt; 4 yr<sup>&minus;1</sup> '
         'below 40' + R + ' in the first 10 d; damped after 30 d',
         f'max |&#7744;<sub>r</sub>|/M<sub>env</sub> below 40{R}: '
         f'{r["mdot_env_max_0_10"]:.1f} yr<sup>&minus;1</sup> in 0&ndash;10 d, median '
         f'{r["mdot_env_median_W"]:.2f} (max {r["mdot_env_max_W"]:.2f}) in {win}',
         'similar start-up transient (slightly weaker); quiet in the window, as the paper'),
        ('Mass loss', 'not quoted (outflow outer boundary)',
         f'{r["Mdot_top_Msun_yr"]:.1e} M<sub>{RS}</sub>/yr through the top (whole-star '
         f'equivalent, mean over {win}: no outflow); domain mass change '
         f'{100 * r["M_tot_change_frac_W"]:+.1g} %', 'ours only (no wind, as expected for '
         'an outflow boundary far above the photosphere)'),
        ('Density inversion, clumping', 'not quantified ("large density variations" '
         'near the surface)',
         (f'shell-mean &rho; inversion at {r["rho_inversion_r"][0]:.1f}&ndash;'
          f'{r["rho_inversion_r"][1]:.1f}{R} in the time mean'
          if r['rho_inversion_r'] else 'no shell-mean &rho; inversion in the time mean') +
         f' (in {100 * r["rho_inversion_frac_dumps"]:.0f} % of the dumps); '
         f'&lang;&rho;<sup>2</sup>&rang;/&lang;&rho;&rang;<sup>2</sup> up to '
         f'{r["clump_max_below_rph"]:.2f} (at {r["clump_max_r"]:.1f}{R}), '
         f'{r["clump_at_rph"]:.2f} at the photosphere', 'ours only'),
        ('Light curve (eq. 8)', 'up to ~10 % (bolometric), 30&ndash;57.4 d',
         f'{r["lc_flux_range"][0]:.3f}&ndash;{r["lc_flux_range"][1]:.3f}, rms '
         f'{r["lc_rms"]:.3f} over {r["lc_window"][0]:.1f}&ndash;{r["lc_window"][1]:.1f} d',
         ('consistent' if max(1 - r['lc_flux_range'][0], r['lc_flux_range'][1] - 1) > 0.05
          else f'<b>weaker</b> by ~{0.1 / max(1 - r["lc_flux_range"][0], r["lc_flux_range"][1] - 1):.0f}&times;')),
        ('PSD', 'red noise, peak &asymp;10&ndash;20 &micro;Hz; 80 % power below &asymp;20 '
         '&micro;Hz (TESS SK&nbsp;&minus;67&nbsp;133: 6&ndash;7)',
         f'80 % below {r["psd_f80"]:.1f} &micro;Hz, 50 % below {r["psd_f50"]:.1f}; peak '
         f'(log-binned, f&middot;P) {r["psd_peak_fP"]:.1f} &micro;Hz; Harvey '
         f'&nu;<sub>char</sub> {r["psd_nu_char"]:.1f} &micro;Hz' if r['psd_f80'] else 'n/a',
         '<b>differs</b>: our power sits at 3&ndash;4 &micro;Hz and is dominated by one '
         f'coherent radial mode (P = {r["mode"]["P_d"]:.2f} d, '
         f'{100 * r["mode"]["var_frac"]:.0f} % of the variance), not convective red noise; '
         'close to TESS\'s 2&ndash;4 &micro;Hz for the wrong reason'),
    ]
    body = ''.join(f'<tr><th scope="row">{a}</th><td>{b}</td><td>{c}</td><td>{d}</td></tr>'
                   for a, b, c, d in rows)
    return {'cmp-body': body, 'w-cmp': win,
            'w-cmp57': f'mean over {win} ({r["n_dumps"]} dumps)'}


if __name__ == '__main__':
    main()
