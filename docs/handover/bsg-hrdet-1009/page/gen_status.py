"""gen_status.py (repro copy): write status.json for the BSG REPRODUCTION page from
scripts_repro/cache/<run basename>/numbers.json (bsg_figs.py), scripts_repro/cache/cost.json
(cost_repro.py) and scripts_repro/cache/ke.json (ke_compare.py).  Run dir: $BSG_REPRO_RUN
(default repro_sf_4n since 10-06).  The page fills its cost banner, run status, figure time labels and the
"Definitions and differences" notes from this file.

usage: python3 gen_status.py PAGEDIR"""
import json
import os
import sys
import time

S = os.environ.get('BSG_WORK', os.path.dirname(os.path.abspath(__file__))).rstrip('/') + '/'
RS = '&#9737;'
RUNDIR = os.environ.get('BSG_REPRO_RUN', '/raven/ptmp/jinma/bsg_raven_1002/repro_sf_4n')
RNAME = os.path.basename(os.path.normpath(RUNDIR))
# 10-07: the run continued on viper (bsg_viper_1006/repro) from the Raven rst; RUNDIR is then
# the merged view (Raven t <= t_rst + viper t > t_rst, build_merged.py), cache/merged_repro is
# a symlink to cache/repro_sf_4n (same dumps, same reductions).
MERGED = os.path.exists(os.path.join(RUNDIR, 'MERGED.txt'))
if MERGED:
    T_RST_D = float(open(os.path.join(RUNDIR, 'MERGED.txt')).readlines()[1].split()[5])
    _ml = open(os.path.join(RUNDIR, 'MERGED.txt')).readlines()
    _t2 = [float(x.split()[5]) for x in _ml if x.startswith('t_rst2')]
    CONT = (f'continued on viper from {T_RST_D:.1f} d (Raven 16 A100 &rarr; viper 16 '
            f'MI300A, same settings and binary commit d6164a6a)')
    if _t2:   # 10-08: third part, the Raven twin repro_sf_cont from the viper rst
        CONT += (f' and finished on Raven from {_t2[0]:.1f} d (16 A100, restarted from the '
                 f'viper restart file)')
    RLABEL = 'repro_sf_4n'
    HW = f'16 blocks, 1 per GPU; {CONT}'
else:
    CONT, RLABEL = '', os.environ.get('BSG_RUN_LABEL', RNAME)
    HW = os.environ.get('BSG_HW_LABEL', '4 Raven nodes, 16 A100, 16 blocks')
d = json.load(open(S + 'cache/' + RNAME + '/numbers.json'))
c = json.load(open(S + 'cache/cost.json'))
ke = json.load(open(S + 'cache/ke.json'))
# RUN before 10-07: ... vet_gd closure; 4 Raven nodes, 16 A100, 16 blocks')
RUN = (f'seeded reproduction run ({RLABEL}): the paper\'s 1-D IC as is, 20&ndash;80 '
       f'R<sub>{RS}</sub>, 1 % temperature seed, force reference work = split, no sponge, '
       f'vet_gd closure; {HW}')
r_bot, r_top = d.get('r_bot_Rsun', 20.0), d.get('r_top_Rsun', 80.0)
grid = d.get('grid', [267, 256, 256])
f2, f3, f4, f6 = d['fig2'], d['fig3'], d['fig4'], d['fig6']
fi = d.get('fig_int')
t_end = f2['t_hst_max_d']
tdump = f2['t_dumps_d']
vr = f4['maps']['Radial Velocity']
rph = f2['r_photo_3D_Rsun']
rhom = f4['maps']['Gas Density']
Tm = f4['maps']['Gas Temperature']


def sci(x):
    e = int(f'{abs(x):e}'.split('e')[1])
    m = x / 10 ** e
    return f'{m:.1f}&times;10<sup>{e}</sup>'


vspread = vr['max'] - vr['min']
rel_rho = (rhom['max'] - rhom['min']) / abs(rhom['mean_physical'])
steady = t_end >= 30.0
if vspread < 0.05 and rel_rho < 1e-4:
    stage = (f'the shell is still laterally uniform (v<sub>r</sub> = {vr["mean_physical"]:+.2f} '
             f'km/s everywhere, &rho; uniform to {rel_rho:.0e}): lateral structure has not '
             f'grown yet')
elif vspread < 5.0:
    stage = (f'lateral structure is growing (v<sub>r</sub> spans {vr["min"]:+.2f} to '
             f'{vr["max"]:+.2f} km/s), but convection has not developed yet')
else:
    stage = (f'convective flows have developed (v<sub>r</sub> {vr["min"]:+.1f} to '
             f'{vr["max"]:+.1f} km/s)')
ke_now, ke_lo, ke_hi = ke['ke_now'], ke['ke_seed_min_after5'], ke['ke_seed_max_after5']
if ke_lo is not None and ke_lo > 3e41:
    conv = (f'Convection has developed: the kinetic energy KE<sub>int</sub> rose from the seed '
            f'and has stayed at {sci(ke_lo)}&ndash;{sci(ke_hi)} erg since about day 5 (now '
            f'{sci(ke_now)} erg), about {ke_now / ke["paper_ke"]:.2f} of the paper\'s '
            f'~10<sup>43</sup> erg.')
else:
    conv = (f'Convection has not settled yet: KE<sub>int</sub> = {sci(ke_now)} erg '
            f'(paper ~10<sup>43</sup> erg).')
FINISHED = t_end >= 57.3
banner = ((f'The run finished at day {t_end:.1f}, the paper\'s end time. ' if FINISHED else
           f'The run is at day {t_end:.1f} of 57.4. ') +
          ('It has passed the 30 d after which the paper calls its envelope steady. '
           if steady else
           'It is in its start-up transient: it starts, like the paper, from the '
           'radiative-only 1-D profile, which is out of thermal and hydrostatic balance, and '
           'the paper calls its envelope steady only after 30 d, so every panel on our side '
           'is still transient. ') + conv +
          f' At {f4["r_cell_Rsun"]:.1f} R<sub>{RS}</sub> {stage}. The 3D photosphere has moved '
          f'from {rph[0]:.1f} to {rph[-1]:.1f} R<sub>{RS}</sub>.')
parts = c.get('parts', [])
part_txt = '; '.join(
    f'{p["name"]} ({p["layout"]}): {p["wall_h"]:.2f} h wall ({p["wall_src"]}, job(s) '
    f'{", ".join(x["job"] for x in p["sacct"]) or ", ".join(p["jobs"])}), '
    f'{p["gpu_h"]:.0f} {p["gpu"]} GPU-hours, t = {p["t0_d"]:.1f}&ndash;{p["t1_d"]:.1f} d'
    for p in parts)
_pace = '; '.join(f'{p["name"]} {p["s_per_cycle_mean"]:.2f} s per step' for p in parts)
cost_note = (f'Measured from the run logs and sacct of all {len(parts)} parts ({part_txt}): '
             f'{c["wall_h"]:.2f} h wall-clock for {c["t_now_d"]:.2f} model days (cycle '
             f'{c["cycle"]}). The run {CONT}. Mean pace per part: {_pace}; '
             f'dt = {c["dt_recent_s"]:.1f} s at the end (it started at 87 s and did not '
             f'shrink with the convection), {grid[0] * grid[1] * grid[2]:,} cells. ' +
             ('' if FINISHED else
              f'The projection assumes the current pace holds to {c["tlim_d"]:.1f} d. ') +
             f'GPU-hours add A100 and MI300A hours (not '
             f'the same unit of work). Ma+2026 do not state their cost.')
st = {
    'c-wall': f'{c["wall_h"]:.1f} h',
    'c-wallcap': (f'wall-clock, all {len(parts)} parts ({c["t_now_d"]:.1f} d of model time)'
                  if FINISHED else
                  f'wall-clock so far ({c["t_now_d"]:.1f} d of model time)'),
    'c-proj': (f'{c["wall_h"]:.1f} h' if FINISHED else f'~{c["proj_total_h"]:.0f} h'),
    'c-projcap': ('total wall-clock for the full 57.4 d (finished)' if FINISHED else
                  'projected wall-clock for 57.4 d at the current dt'),
    'c-gpuh': f'{c["gpu_h"]:.0f}',
    'c-gpuhcap': (('GPU-hours in total, ' if FINISHED else 'GPU-hours so far, ') + ' + '.join(f'{p["gpu_h"]:.0f} {p["gpu"]}'
                                                    for p in parts) +
                  ('' if FINISHED else f' (~{c["proj_gpu_h"]:.0f} projected)')),
    'c-note': cost_note,
    's-time': f'{t_end:.2f} d',
    's-run': RUN,
    's-grid': f'{grid[0]} &times; {grid[1]} &times; {grid[2]}',
    's-snaps': f'0&ndash;{tdump[-1]:.2f} d, every 0.5 d ({len(tdump)})',
    's-updated': time.strftime('%Y-%m-%d %H:%M CEST'),
    's-banner': banner,
    'w-fig2': f't = 0&ndash;{t_end:.2f} d',
    'w-fig3': f't = {f3["t_d"]:.2f} d',
    'w-fig4': f't = {f4["t_d"]:.2f} d, r = {f4["r_cell_Rsun"]:.2f} R<sub>{RS}</sub>',
    'w-fig6': f't = {f6["t_range_d"][0]:.2f}&ndash;{f6["t_range_d"][1]:.2f} d',
    'n-fig2': (f'Energies are sums over our wedge (&pi;/3 sr), as in the paper\'s appendix. '
               f'M<sub>env</sub> is the whole-star equivalent of the wedge mass, '
               f'{f2["M_envelope_whole_star_Msun"]:.3f} M<sub>{RS}</sub> for '
               f'{r_bot:.0f}&ndash;{r_top:.0f} R<sub>{RS}</sub> (the paper\'s domain). Colour maps '
               f'have one column per 0.5-d snapshot and are drawn to 65 R<sub>{RS}</sub> like '
               f'the paper\'s. So far &Delta;E<sub>grav</sub> = '
               f'{sci(f2["dE_grav_1e44"] * 1e44)} erg, &Delta;U = '
               f'{sci(f2["dU_int_1e44"] * 1e44)} erg and E<sub>k</sub> = '
               f'{sci(f2["E_k_1e44"] * 1e44)} erg (paper: about 5&times;10<sup>44</sup> and '
               f'10<sup>43</sup> erg in its transient). The 3D photosphere is at '
               f'{min(rph):.1f}&ndash;{max(rph):.1f} R<sub>{RS}</sub>; the IC\'s is at '
               f'{f3["r_photo_IC_Rsun"]:.1f} R<sub>{RS}</sub>.'),
    'n-fig3': (f'P<sub>rad</sub> = D<sub>rr</sub>E, the Eddington tensor the vet_gd closure '
               f'used, from the m1_vet dump of the same time. The right panel is the paper\'s '
               f'IC as run: radiative only (F<sub>MLT</sub> = 0), so it has no convective '
               f'shell to shade, where the paper\'s panel marks the two MESA convection zones. '
               f'Its Fe opacity peak is at {f3["fe_peak_ic_Rsun"]:.1f} R<sub>{RS}</sub> '
               f'(&kappa;<sub>R</sub> = {f3["fe_peak_ic_kappa"]:.2f} cm<sup>2</sup> g<sup>&minus;1</sup>, '
               f'&kappa;/&kappa;<sub>Edd</sub> = {f3["fe_peak_ic_kappa"] / f3["kappa_edd"]:.2f} '
               f'with &kappa;<sub>Edd</sub> = 4&pi;cGM/L = {f3["kappa_edd"]:.2f}), not at the '
               f'35 R<sub>{RS}</sub> the paper quotes. At {f3["t_d"]:.2f} d the shell-mean Fe '
               f'peak is at {f3["fe_peak_3d_Rsun"]:.1f} R<sub>{RS}</sub> '
               f'(&kappa;<sub>R</sub> = {f3["fe_peak_3d_kappa"]:.2f}); the He&nbsp;II bump near '
               f'{f3["heii_bump_3d_Rsun"]:.1f} R<sub>{RS}</sub> (T &asymp; 4&times;10<sup>4</sup> K) '
               f'reaches &kappa;<sub>R</sub> = {f3["heii_bump_3d_kappa"]:.2f}. Largest '
               f'&kappa;/&kappa;<sub>Edd</sub> in the plane: {f3["kap_over_edd_max_3d"]:.2f}. '
               f'Largest ram/total pressure below the photosphere: '
               f'{f3["ram_ratio_max_below_rphoto"]:.3g} (paper: about 0.5 above 33 '
               f'R<sub>{RS}</sub> in steady state). Opacities: Linhao Ma\'s TOPS tables '
               f'(padded with the edge values), the tables the run uses.'),
    'n-fig4': (f'The shell stays at 40.6 R<sub>{RS}</sub> (nearest cell {f4["r_cell_Rsun"]:.2f}), '
               f'as on the earlier page. The paper puts it 1 R<sub>{RS}</sub> above its Fe peak '
               f'(35 &rarr; 36 R<sub>{RS}</sub>); in this run the Fe peak of the same IC is at '
               f'{f3["fe_peak_3d_Rsun"]:.1f} R<sub>{RS}</sub>, so 40.6 is just inside it. At '
               f'{f4["t_d"]:.2f} d: v<sub>r</sub> from {vr["min"]:+.3f} to {vr["max"]:+.3f} km/s; '
               f'shell means &rho; = {rhom["mean_physical"]:.2f}&times;10<sup>&minus;8</sup> '
               f'g cm<sup>&minus;3</sup> and T = {Tm["mean_physical"]:.3f}&times;10<sup>5</sup> K '
               f'(paper, steady: 5&ndash;8&times;10<sup>&minus;8</sup>, '
               f'1.5&ndash;1.7&times;10<sup>5</sup> K). v<sub>r</sub> uses the paper\'s range '
               f'(&plusmn;60 km/s) once our flows reach 5 % of it, else our own; &rho; and T '
               f'are absolute values with a fixed range = shell mean of the newest dump '
               f'&plusmn; a relative amplitude (1&ndash;99th percentile, at least 0.1 %, so '
               f'round-off is not painted as structure).'),
    'n-fig6': (f'Normalised as the paper\'s eq. 8: the wedge luminosity (L<sub>top</sub>, the '
               f'Marshak top face) is linearly detrended and its fluctuation scaled by '
               f'&radic;(&Omega;/4&pi;) = 0.289 to a whole-star equivalent. Both panels use '
               f'(&Delta;F/F)<sup>2</sup>/&micro;Hz and the paper\'s axis range. Our flux spans '
               f'{f6["flux_min"]:.3f}&ndash;{f6["flux_max"]:.3f} over '
               f'{f6["t_range_d"][0]:.2f}&ndash;{f6["t_range_d"][1]:.2f} d; mean whole-star '
               f'L<sub>top</sub> = {sci(f6["L_top_mean_whole_star_erg_s"])} erg/s (IC: '
               f'5.58&times;10<sup>38</sup>). The dotted line marks our Nyquist frequency '
               f'(hst every {f6.get("cadence_s", 1000):.0f} s). ' +
               ('Only the window after 30 d compares with the paper.' if steady else
                'Before 30 d the whole curve is the start-up transient (the window falls back '
                'to the full record), so the orange "radial mode" fit is not meaningful yet.')),
    **({'w-int': (f't = {fi["t_d"]:.2f} d, own &tau;<sub>R</sub> = 1, r = '
                  f'{fi["r_own_p1_p99_Rsun"][0]:.1f}&ndash;{fi["r_own_p1_p99_Rsun"][1]:.1f} '
                  f'R<sub>{RS}</sub>'),
        'n-int': (f'At {fi["t_d"]:.2f} d, F<sub>r</sub> is taken at each column\'s own '
                  f'&tau;<sub>R</sub> = 1 cell, at {fi["r_own_p1_p99_Rsun"][0]:.2f}&ndash;'
                  f'{fi["r_own_p1_p99_Rsun"][1]:.2f} R<sub>{RS}</sub> (1st&ndash;99th '
                  f'percentile; median {fi["r_ph_Rsun"]:.2f}). There '
                  f'F<sub>r</sub>/&lang;F<sub>r</sub>&rang; spans {fi["q_min"]:.3f}&ndash;'
                  f'{fi["q_max"]:.3f} (rms {fi["q_sd"]:.3f}; F<sub>r</sub> &lt; 0 in '
                  f'{100 * fi["fneg"]:.2f} % of the columns), &lang;T<sub>eff</sub>&rang; = '
                  f'{fi["Teff_mean"]:.0f} K and 4&pi;&lang;r<sup>2</sup>F<sub>r</sub>&rang; = '
                  f'{sci(fi["L_4pir2F"])} erg/s = {fi["L_over_Lstar"]:.3f} L.')}
       if fi else {}),
    'w-ke': f't = 0&ndash;{ke["t_end_d"]:.2f} d (unseeded full run: 0&ndash;{ke["t_end_full_d"]:.2f} d)',
    'n-ke': (f'KE<sub>int</sub> is column 11 of user.hst (gas kinetic energy inside r<sub>int</sub>, '
             f'summed over the wedge). Unseeded full run: peak {sci(ke["ke_full_max"])} erg at '
             f'{ke["t_full_max_d"]:.1f} d (the start-up adjustment), then a decline to '
             f'{sci(ke["ke_full_7d"])} erg at 7 d and {sci(ke["ke_full_end"])} erg when it '
             f'stopped at {ke["t_end_full_d"]:.1f} d. Seeded split run: peak '
             f'{sci(ke["ke_seed_max"])} erg at {ke["t_seed_max_d"]:.1f} d, '
             f'{sci(ke["ke_seed_7d"])} erg at 7 d, {sci(ke["ke_now"])} erg now. The dashed '
             f'line is the paper\'s kinetic-energy level read from its Fig. 2 '
             f'(~10<sup>43</sup> erg), a reference only.'),
    '__v': str(int(time.time())),
}
out = sys.argv[1]
json.dump(st, open(out.rstrip('/') + '/status.json', 'w'), indent=1)
print('wrote', out.rstrip('/') + '/status.json')
