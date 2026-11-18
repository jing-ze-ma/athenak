#!/usr/bin/env python3
"""accsp_prod_page (11-02): publish/index.html of "RY Per Accretor 3-D" for the sp full-sphere PRODUCTION (accsp_1102/prod/fresh).
New figures/movies have an sp_ prefix (publish/img/sp_*, publish/movie/sp_*). The earlier runs stay on the artifact and are only
REFERENCED here, never re-sent: cs N=128 test img/fig_*.{png,jpg} + movie/cs128_*, RUN12 img/run12_* + movie/acc3d_*.
CSS = acc3d_run12_page/bundle/page.css (the page's existing theme tokens). Text from img/status.json; takeaways computed."""
import os
import json
import html
import shutil
import time

P = '/resnick/groups/carnegie_poc/jingze/accsp_prod_page'
R12 = '/resnick/groups/carnegie_poc/jingze/acc3d_run12_page'
IMG, PUB = P + '/img', P + '/publish'
CSS = open(R12 + '/bundle/page.css').read() + '\n.take { font-weight: 600; color: var(--fg); max-width: 86ch; margin: 0; }\n'
S = json.load(open(IMG + '/status.json'))
PORB, DAYS, TS = 0.473169, 8.0521, 0.01
esc = html.escape
os.makedirs(PUB + '/img', exist_ok=True)
NEW = ('sp_hist.png', 'sp_keband.png', 'sp_slices_latest.jpg', 'sp_phot_latest.jpg', 'sp_tip.png', 'sp_trc_budget.png', 'sp_trc_maps.jpg', 'sp_trc_prof.png')
TRC = json.load(open(IMG + '/trc.json')) if os.path.exists(IMG + '/trc.json') else {}   # trc-1102
for f in NEW:
    if os.path.exists(f'{IMG}/{f}'):
        shutil.copy2(f'{IMG}/{f}', f'{PUB}/img/{f}')


def fig(fn, alt, take, cap, kept=False):
    if not kept and not os.path.exists(f'{PUB}/img/{fn}'):
        return f'<p class="status">[figure {esc(fn)} not built yet]</p>'
    tk = f'<p class="take">{take}</p>' if take else ''
    return f'{tk}<figure><img alt="{esc(alt)}" loading="lazy" src="img/{fn}"><figcaption>{cap}</figcaption></figure>'


def movie(name, cap, kept=False):
    if not kept and not os.path.exists(f'{PUB}/movie/{name}.mp4'):
        return '<p class="status">movie not built yet</p>'
    n = ''
    try:
        j = json.load(open(f'{PUB}/movie/{name}_frames.json')); n = f" {j['n']} frames at {j['fps']} fps."
    except Exception:
        pass
    return (f'<figure><video controls muted playsinline preload="metadata" poster="movie/{name}_poster.png" src="movie/{name}.mp4"></video>'
            f'<figcaption>{cap}{n}</figcaption></figure>')


def g(x, f='%.3g'):
    return '–' if x is None else (f % x)


t = S.get('t', 0.0); sr, ss = S['stream'], S['static']
ms_last, ms_all = sr['ms'] if sr['ms'] else (None, None)
on = t - TS
reds = S.get('reds', [])
lastred = reds[-1] if reds else None
health_bad = bool(S.get('hst_nan') or sr['nan'] or sr['fatal'] or (lastred and lastred['nonfinite']))
health = ('<b>flag</b>: ' if health_bad else 'clean: ') + (f"hst NaN {S.get('hst_nan')}, run.log NaN/FATAL lines {sr['nan']}/{sr['fatal']}, "
                                                            f"non-finite cells in the latest dump {lastred['nonfinite'] if lastred else '–'}")
rate = None
if ms_all and S.get('dt'):
    rate = S['dt'] / (ms_all / 1e3) * 3600          # code time per wall hour at the current dt
kv = [
    ('status', f"t = {t:.5f} code = {t * DAYS * 24:.2f} h = {t / PORB:.4f} orbits; stream on since t 0.0100 ({on * DAYS * 24:.2f} h, "
               f"{on / PORB:.4f} orbits); cycle {sr['cycle']}, dt {S.get('dt', 0):.3g}; run {esc(S['state'])}"),
    ('speed', f"{g(ms_last, '%.0f')} ms/cycle over the last 50 cycles, {g(ms_all, '%.0f')} ms/cycle since the restart"
              + (f"; ≈ {rate * 24 / PORB:.2f} orbits per wall day at this dt" if rate else '')
              + '; 1 node × 4 H200 (hpc-sm-01-05, alloc 5719508), 8 MeshBlocks 500 × 116 × 64'),
    ('health', health + f"; stream: {sr['posfix']} time2_posfix re-solves, {sr['BE']} backward-Euler step redos"),
    ('binary', 'athena_gpu_ry_per_accretor_r2 (branch ringsmooth-1102 e67c9e967 = accsp-1102 167c05ce4 + the polar ring-averaging fix; md5 c9163a55)'),
    ('keys', '<code>rad_m1/polar_ring_avg = true</code>, <code>polar_ring_avg_ke = rad</code>, <b style="display:inline;color:var(--accent)">'
             '<code>polar_ring_avg_gz = true</code></b> (new: the coarse-group tangential radiation force deposited smoothly across the phi faces), '
             '<code>vet_gd = true</code> (twin off), <code>problem/roche_wb_tan = true</code> (WB tide), <code>floor_m1_vcap = true</code>, '
             '<code>time2_floor_eint = 1e-7</code>'),
    ('decks', 'static <code>ringsmooth_1102/decks/sp_static_gz.athinput</code> t 0 – 0.01 (run ringsmooth_1102/long/L_gz, same binary); stream '
              '<code>sp_stream_gz.athinput</code> from its rst 00002, tlim 1e9; outputs hydro_w every 0.005, m1 every 0.01, rst every 0.05'),
    ('grid', 'spherical polar, full sphere: nx 500 (r 6.3 – 24.34 R<sub>⊙</sub> stretched) × 116 θ × 512 φ; θ plateau |lat| ≤ 20.39° at 0.703°, '
             'stretched caps (q 1.073) to 5.45° polar cells; φ 0.703°; the 7 polar rows ring-averaged in φ groups 16, 4, 4, 2, 2, 2, 2'),
    ('now', f"L<sub>meas</sub> {g(S.get('Lmeas'), '%.4g')} L<sub>⊙</sub> · L<sub>top</sub> {g(S.get('Ltop'), '%.4g')} L<sub>⊙</sub> · "
            f"M<sub>in</sub> {g(S.get('Min'), '%.4g')} · KE {g(S.get('KE'), '%.4g')} (lateral {g(S.get('KElat'), '%.4g')}) code · "
            f"M<sub>fl</sub> {g(S.get('Mfl'), '%.3g')} · Picard {g(S.get('Pic'), '%.2f')} (max {g(S.get('Picmax'), '%.0f')})"),
    ('system', 'M<sub>acc</sub> 16, M<sub>don</sub> 18 M<sub>⊙</sub>, a 33.26 R<sub>⊙</sub>, P 3.81 d (1 code time = 8.05 d); R<sub>acc</sub> ≈ 9.0 R<sub>⊙</sub>, '
               'T<sub>eff</sub> 28.3 kK; d<sub>L1</sub> 16.23; stream 1e-4 M<sub>⊙</sub>/yr from a nozzle at r 14.5'),
    ('physics', '3-D Roche gravity + WB tangential tide, Coriolis + centrifugal; implicit grey M1 (hesdirk2) with the global-direction VET (vet_gd); '
                'MESA-kap Rosseland + Planck tables; general EOS (Saha H/He)'),
]
kvh = '\n'.join(f'<div><b>{k}</b>{v}</div>' for k, v in kv)

# ---- takeaways (computed)
tk_hist = (f"dt {S.get('dt', 0):.3g} now, smallest in the stream phase {g(S.get('dt_min_stream'))} (cs N=128 test: {g(S.get('cs_dt_min_stream'))}); "
           f"L<sub>meas</sub> {g(S.get('Lmeas'), '%.4g')} L<sub>⊙</sub>; {'no NaN/FATAL so far' if not health_bad else 'see the health flag'}.")
kb = S.get('keband')
if kb:
    rat = kb['sp_tot'] / kb['cs_tot'] if kb['cs_tot'] else float('nan')
    ratp = kb['sp_pole'] / kb['cs_pole'] if kb['cs_pole'] else float('nan')
    tk_ke = (f"At t {kb['t_sp']:.4f} the sp lateral KE is {rat:.2f}× the cs N=128 test (|lat| &gt; 60°: {ratp:.2f}×)"
             + (' – the polar bands match the cs test.' if ratp < 1.5 else ' – the polar bands still carry more lateral KE than on the cubed sphere.')
             + (f" Latest (t {kb['t_latest']:.4f}): {kb['latest_tot']:.3g} code; the cs128 test stopped at t {kb['t_cs']:.3f}." if kb['t_latest'] > kb['t_sp'] + 1e-6 else ''))
else:
    tk_ke = 'Not built yet (needs a reduced dump).'
tip = S.get('tip')
tk_tip = (f"Median R<sub>phot</sub> {lastred['Rph_med']:.3f} R<sub>⊙</sub> (2nd percentile {lastred['Rph_p2']:.3f}); at the L1 tip v<sub>r</sub> is "
          f"{tip['vr_min']:.1f} to {tip['vr_max']:.1f} km/s and ρ up to {tip['rho_ratio_max']:.3g}× the t = 0 value." if (lastred and tip) else '')
sd, pd = S.get('slice_dump', {}), S.get('phot_dump', {})
built = time.strftime('%Y-%m-%d %H:%M %Z')
warn = ', '.join(esc(w) for w in sr.get('warn', []))

page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>RY Per Accretor 3-D</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Source+Sans+3:wght@400;600&family=Source+Code+Pro:wght@400;500&display=swap">
<style>{CSS}</style></head><body>
<main>
<header>
<p class="note"><b>Since 11-02 evening this page shows the spherical-polar full-sphere PRODUCTION run</b> (launched 11-02 18:44, stream phase from t 0.01).
The cubed-sphere N = 128 test (shown here earlier on 11-02) and RUN12 (10-24 – 10-30) are kept at the bottom.</p>
<h1>RY Per accretor in 3-D: full-sphere production</h1>
<p class="lede">AthenaK radiation hydrodynamics of the 16 M<sub>⊙</sub> mass gainer of RY Per on the full sphere (spherical polar, 500 × 116 × 512, stretched θ with
a fine ±20° core), with the 3-D Roche potential, the WB tangential tide, Coriolis and centrifugal forces, in the frame co-rotating with the binary
(donor towards +x, longitude 0). Static phase to t 0.01, then the 1e-4 M<sub>⊙</sub>/yr L1 stream; the polar ring averaging uses the new
<code>polar_ring_avg_gz</code> fix that removed the grid-scale v<sub>φ</sub> noise in the polar rows.</p>
<div class="kv">
{kvh}
</div>
<p class="status">page built {built} on Caltech by accsp_prod_page/refresh.sh from the hst, logs and dumps present then (read only){'; unread input keys: ' + warn if warn else ''}</p>
</header>
<section>
<h2>History</h2>
{fig('sp_hist.png', 'history', tk_hist, 'Static phase (t &lt; 0.01, the ringsmooth_1102/long/L_gz run that the production restarts from) and the stream phase; dotted line = stream on; dashed = the cs N = 128 test (test128wb) over its t 0 – 0.03. Row 1: dt, wall ms per cycle (run.log every 50 cycles), L<sub>meas</sub> and L<sub>top</sub>. Row 2: nozzle mass M<sub>in</sub>/M<sub>nzf</sub>, total KE, lateral and radial KE from hydro.hst (cs: panel components, no metric cross term). Row 3: booked floor mass M<sub>fl</sub> with the cumulative M1 fallbacks (time2_posfix re-solves, backward-Euler redos), Picard iterations, WB torque and the stream tracer mass.')}
</section>
<section>
<h2>Lateral KE by latitude: sp vs the cs N = 128 test</h2>
{fig('sp_keband.png', 'KE by latitude band', tk_ke, 'Lateral KE ½ρ(v<sub>θ</sub>² + v<sub>φ</sub>²) dV by |latitude| band (accsp_1102/tools/keband.py: sp orthonormal components, cs exact gnomonic metric) at the latest sp dump that the cs N = 128 test also has (it ran to t 0.03), the deep (r &lt; 8) and envelope (r 8 – 9.3) parts, both totals and the |lat| &gt; 60° part vs time, and the φ-rms of v<sub>φ</sub> in the 10 rows next to each pole at r 6.8 (the 11-02 source of the excess lateral KE before the gz fix).')}
</section>
<section>
<h2>Slices</h2>
{fig('sp_slices_latest.jpg', 'slices', '', f"Orbital plane (top, the two θ rows next to the equator averaged) and the x–z plane through the donor axis (bottom, φ = 0 and φ = π half planes, pole to pole), t = {sd.get('t', 0):.4f} (dump {sd.get('idx', '–')}, latest with an m1 dump). Columns: log ρ, log T<sub>gas</sub> (run EOS), log P<sub>rad</sub>/P<sub>gas</sub> (E/3), v<sub>r</sub>, v<sub>φ</sub> (rotating frame). T<sub>gas</sub> and P<sub>rad</sub>/P<sub>gas</sub> are shown everywhere, including the floor-class gas (ρ &lt; 1e-7 code) around the star, whose T is EOS noise on the H-ionisation step and which holds no mass. Dashed circle = d<sub>L1</sub>; dark disc = the inner boundary r 6.3.")}
</section>
<section>
<h2>Photosphere and the L1 tip</h2>
{fig('sp_phot_latest.jpg', 'photosphere', tk_tip, f"t = {pd.get('t', 0):.4f} (dump {pd.get('idx', '–')}). Top: Rosseland photosphere R(τ<sub>R</sub> = 2/3) per (φ, θ) column (τ inward from x1max with the deck's MESA-kap table, T<sub>gas</sub> from the run EOS; longitude 0 = towards the donor) and T<sub>gas</sub> there; native cells, so the polar caps show the 5.45° cells. Bottom: the L1-facing tip, log ρ and v<sub>r</sub> in the orbital and x–z planes (dashed = median R<sub>phot</sub>).")}
{fig('sp_tip.png', 'tip', '', 'L1-tip column (φ = 0, the 2 × 2 cells at the equator): ρ/ρ<sub>t=0</sub> and v<sub>r</sub> (up to 12 dumps; dotted = median R<sub>phot</sub>), and the solid-angle-weighted median and 2nd percentile of R<sub>phot</sub> vs time.')}
</section>
<section>
<h2>Stream tracer</h2>
{(fig('sp_trc_budget.png', 'stream tracer budget', esc(TRC.get('take_budget', '')), 'Passive scalars s0 (gas written by the nozzle, the stream) and s1 (gas that entered through the other user boundaries), from ryper.trc.hst since the stream was switched on: injected Nz0, in the domain Ms0, removed by the sinks Sk0, out through r<sub>in</sub> Fi0 and through r<sub>out</sub> + the θ walls Fo0 (cumulative registers of the pgen, acctrc-1019). Middle: budget closure ΔMs − (sources − sinks − outflows), relative to the injected s0 / the s1 inflow. Right: the s1 budget.') + fig('sp_trc_maps.jpg', 'stream tracer maps', esc(TRC.get('take_maps', '')), 'Stream-gas fraction s0 (log scale, dark = below 1e-4) in the orbital plane (the two θ rows next to the equator), the x–z plane through the donor axis, and in the cell of each (φ, θ) column where the Rosseland optical depth reaches 2/3 (same photosphere as above). Dashed circle = d<sub>L1</sub>, dark disc = r<sub>in</sub>.') + fig('sp_trc_prof.png', 'stream tracer radial profile', esc(TRC.get('take_prof', '')), 'Left: stream-gas mass per unit radius, ρ s0 summed over each spherical shell, for up to 10 stream dumps (colour = time). Right: the stream-gas fraction of each shell mass: how deep the accreted gas has mixed. Dotted: median R<sub>phot</sub> of the latest dump and r<sub>in</sub>.')) if TRC.get('present') else '<p class="status">tracer section not built yet (' + esc(TRC.get('why', 'no trc.json')) + ')</p>'}
</section>
<section>
<h2>Movies</h2>
{movie('sp_eq', 'Orbital plane and x–z plane, same panels as the slice figure, one frame per hydro dump (static every 0.0025, stream every 0.005 code); dumps without an m1 dump at the same time leave the P<sub>rad</sub>/P<sub>gas</sub> panels empty.')}
{movie('sp_phot', 'Photosphere and L1 tip, same panels as the photosphere figure, one frame per hydro dump.')}
</section>
<section>
<h2>Earlier runs</h2>
<h3>cs N = 128 test (test128wb, 11-01 – 11-02)</h3>
<p>Cubed sphere, 6 × 128 × 128 columns × 500 radial cells, same physics and keys except the sp pole handling; static to t 0.01, stream to
t 0.03 (5.8 h, 0.063 orbits), 2.4 s per cycle on 4 H200. It is the reference for the lateral-KE check above. Figures and movies are its final ones.</p>
{fig('fig_hist.png', 'cs128 history', '', 'cs N = 128: hst history (static + stream; dotted = stream on at t 0.01).', kept=True)}
{fig('fig_ke_static.png', 'cs128 static KE', '', 'cs N = 128 static phase: KE with and without the WB tide, and d ln KE/dt.', kept=True)}
{fig('fig_cmp96.png', 'cs128 vs cs96', '', 'cs N = 128 stream phase against the cs96 long tests long_base and long_F.', kept=True)}
{fig('fig_slices_latest.jpg', 'cs128 slices', '', 'cs N = 128 slices at t 0.03 (orbital and x–z planes).', kept=True)}
{fig('fig_phot_latest.jpg', 'cs128 photosphere', '', 'cs N = 128 photosphere and L1 tip at t 0.03.', kept=True)}
{fig('fig_tip.png', 'cs128 tip', '', 'cs N = 128 L1-tip column and R<sub>phot</sub> vs time.', kept=True)}
{movie('cs128_eq', 'cs N = 128: orbital and x–z planes, one frame per hydro dump.', kept=True)}
{movie('cs128_phot', 'cs N = 128: photosphere and L1 tip, one frame per hydro dump.', kept=True)}
<h3>RUN12 (spherical-polar ±20° wedge, 10-24 – 10-30)</h3>
<p>±20.4° wedge about the orbital plane (420 × 58 × 512, closed θ walls), static to t 0.05 then the stream; last shown at t = 2.0244 code
= 16.3 d = 4.28 orbits (405 hydro dumps).</p>
{fig('run12_hist.png', 'run12 history', '', 'RUN12 hst history to t 2.02 (10-30).', kept=True)}
{fig('run12_map.jpg', 'run12 map', '', 'RUN12 maps at t = 2.020 (hydro dump 404); grey = cells the code resets every stage.', kept=True)}
{movie('acc3d_3d', 'RUN12 3-D view (r ≤ d<sub>L1</sub>): nested density isosurfaces coloured by T<sub>gas</sub>, to t 2.02.', kept=True)}
{movie('acc3d_eq', 'RUN12 orbital plane and edge-on cut to t 2.02.', kept=True)}
</section>
<section>
<h2>Caveats</h2>
<ul>
<li>P<sub>rad</sub> = E/3 is the isotropic value; in the optically thin atmosphere the VET P<sub>rr</sub> differs.</li>
<li>The photosphere is a radial ray trace per column (no oblique rays); the stream columns near longitude 0 reach τ = 2/3 at the nozzle.</li>
<li>The floor-class cells (ρ &lt; 1e-7 code) are mostly the ambient around the star, not failures; M<sub>fl</sub> is the floor mass the code books.</li>
<li>The cs N = 128 comparison ends at t 0.03; later sp times are compared with that end state only.</li>
</ul>
</section>
<footer>Page built by /resnick/groups/carnegie_poc/jingze/accsp_prod_page/refresh.sh (reduce_sp.py → mk_figs.py → ffmpeg → build_page.py, adapted
from acc3d_cs128_page) from accsp_1102/prod/fresh/stream and ringsmooth_1102/long/L_gz (static phase); comparison csacc_1031/test128wb.</footer>
</main>
</body></html>
'''
open(PUB + '/index.html', 'w').write(page)
print('index.html', os.path.getsize(PUB + '/index.html'), 'B')
