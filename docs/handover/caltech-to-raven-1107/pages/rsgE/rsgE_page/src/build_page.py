#!/usr/bin/env python3
"""RSGPAGE 10-19: assemble publish/index.html for the RSG g25 cs production (figures inlined base64, one file).
usage: build_page.py <imgdir> <out.html> <status>"""
import sys
import os
import json
import base64
import html

IMG, OUTF, STATUS = sys.argv[1], sys.argv[2], sys.argv[3]
T = json.load(open(IMG + '/table.json'))
HERE = os.path.dirname(os.path.abspath(__file__))
CSS = open(os.path.join(HERE, 'page.css')).read()
SITE = os.environ.get('SITE', 'Raven')
LAYOUT = os.environ.get('LAYOUT', '6 nodes × 4 A100, 24 MeshBlocks 591×64×64 (1 per GPU); chain of 24 h links')
PUB = os.path.dirname(os.path.abspath(OUTF))


def img(fn, alt):
    p = os.path.join(IMG, fn)
    if not os.path.exists(p):
        return f'<p class="status">[missing figure {fn}]</p>'
    mt = 'image/jpeg' if fn.endswith('.jpg') else 'image/png'
    b = base64.b64encode(open(p, 'rb').read()).decode()
    return f'<img alt="{html.escape(alt)}" loading="lazy" src="data:{mt};base64,{b}">'


def fig(fn, alt, cap):
    return f'<figure>{img(fn, alt)}<figcaption>{cap}</figcaption></figure>'


jobs = T['jobs']
running = [j for j in jobs if len(j) > 1 and j[1] == 'RUNNING']
queued = [j for j in jobs if len(j) > 1 and j[1] == 'PENDING']
jobtxt = '; '.join(f'{j[0]} {j[1].lower()} {j[2]}' for j in jobs) or f'squeue: no {T.get("jobname", "")} jobs'
health_bad = T['fatal'] + T['div'] + T['nc'] + T['nonfinite'] + T['nan']
health = (f'<span class="{"bad" if health_bad else "ok"}">FATAL {T["fatal"]}, DIVERGED {T["div"]}, summary NC {T["nc"]}, '
          f'NaN in logs {T["nan"]}, non-finite cells in dumps {T["nonfinite"]}</span>')

slices = ''.join(fig(fn, 'slices t %.1f d' % t, f'Slices at t = {t:.1f} d (latest reduced dump).') for fn, t in T['slices'])
phot = ''.join(fig(fn, 'photosphere band t %.1f d' % t, f'Photosphere band 0.9–1.2 R<sub>ph</sub>, equatorial ring unwrapped, t = {t:.1f} d.') for fn, t in T.get('phot', []))
maps = ''.join(fig(fn, 'maps t %.1f d' % t, f'Mollweide maps at t = {t:.1f} d.') for fn, t in T['maps'])
MJ = json.load(open(IMG + '/movies.json')) if os.path.exists(IMG + '/movies.json') else None
EV = T.get('ev', {})
evtxt = ('; '.join(f'{k} {EV[k]}' for k in ('eos_dfloor', 'eos_tfloor', 'eos_vceil', 'eos_fail', 'fofc', 'm1_pos_flr') if k in EV)
         + f' (cell-stages in the 100 cycles before cycle {EV.get("cycle", "?")}, <code>rsg.log</code>)') if EV else 'rsg.log not read'

kv = [
    ('status', f'{html.escape(SITE)}: {len(running)} link running, {len(queued)} queued — {html.escape(jobtxt)}'),
    ('now', f't = {T["t_s"]:.4g} s = <b>{T["t_d"]:.2f} d</b>, cycle {T["cyc"]}, dt {T["dt"]:.3g} s, {T["spc"]:.2f} s/cycle '
            f'(≈ {86400 / max(T["dt"], 1e-30) * T["spc"] / 3600:.2f} h wall per sim day); Picard mean {T["pic"]:.2f}, max {T["picmax"]:.0f}'),
    ('health', health),
    ('floors / limiters', evtxt + (f'; T<sub>gas</sub> at the 2000 K floor in the last dump: {T["nfloorT"]:,} cells '
                                   f'({100 * T["nfloorT"] / max(T["ncell"], 1):.1f} % of the mesh)' if T.get('nfloorT', -1) >= 0 else '')),
    ('last dump', f'{T["t_dump"]:.1f} d ({T["ndump"]} dumps reduced): median R<sub>ph</sub> {T["rph"]:.4f} R<sub>ph,IC</sub>, '
                  f'T<sub>eff</sub> {T["teff"]:.0f} K, max |v<sub>r</sub>| {T["vmax"]:.1f} km/s'),
    ('budget', f'L<sub>top</sub>/L<sub>in</sub> {T["LL"]:.4f}, L<sub>bot</sub>/L<sub>in</sub> {T["Lbot"]:.4f}; '
               f'ΔE<sub>tot</sub>/|E<sub>0</sub>| {T["dE"]:.2e} vs ∫(L<sub>in</sub>−L<sub>top</sub>)dt/|E<sub>0</sub>| {T["dEb"]:.2e}; '
               f'ΔM/M {T["dM"]:.2e}; KE<sub>int</sub> {T["ke"]:.3e} erg'),
    ('star', 'RSG “golden25” MESA model (profile438): R<sub>ph</sub> 5.800×10<sup>13</sup> cm = 834 R<sub>⊙</sub> (IC column), '
             'L 6.85×10<sup>38</sup> erg/s = 1.79×10<sup>5</sup> L<sub>⊙</sub>, T<sub>eff</sub> 4.11 kK; X 0.654, Z 0.0198; '
             'deep IC column <code>ic_rsg_g25_deep.txt</code> (bottom 11 R<sub>⊙</sub>)'),
    ('grid', 'Global cubed sphere (gnomonic equiangular), N = 128 per face (0.70°), 6 panels; 548 stretched radial cells '
             'r = 1.761×10<sup>12</sup> – 1.761×10<sup>14</sup> cm (0.03 – 3 R<sub>ph</sub>; geometric head below the old 591-cell face 110); '
             'lateral cell merging ×8 / ×4 / ×2 below 2.06 / 3.43 / 5.72×10<sup>12</sup> cm (radiation unmerged); 54 M cells'),
    ('layout', html.escape(LAYOUT)),
    ('physics', 'Implicit grey M1 radiation with the VET closure (vet_gd), general EOS (H/He ionisation), MESA kap Rosseland table, '
                'adaptive MLT closure with a cosine ramp: start 43.0 d, length 51.4 d = 7/σ<sub>max</sub> from the IC N<sup>2</sup>, MLT off at 94.4 d; α<sub>MLT</sub> 3; angular-momentum transport on; shared seam flux; '
                'floors ρ 10<sup>−20</sup>, T 2000 K, sponge on ρ &lt; 2×10<sup>−17</sup> at r &gt; 1.58×10<sup>14</sup> cm; '
                'dt_retry (cut 0.9, regrow 1.1 after 20 good steps)'),
    ('seed', 'N<sup>2</sup>-placed surface seed: real spherical harmonics Y<sub>lm</sub>, l = 16–32 (833 terms, all m), amplitude 10<sup>−3</sup>, rng 1234, '
             'radial band 0.97–0.999 R<sub>ph</sub> (5.626–5.794×10<sup>13</sup> cm), where the IC N<sup>2</sup> growth is fastest (0.993–0.996 R<sub>ph</sub>, 1/σ = 7.4 d, H recombination)'),
    ('binary / deck', 'fin51-r-fix ff1f52bc <code>athena_gpu_he_star_m1_ff1f52bc_nofma</code> md5 cd998ac0; deck <code>rsg_E_1022_n2seed.athinput</code> = the production deck + '
                      '<code>implicit_face_opac = harmonic</code>, <code>time2_vet_col = extrap2</code>, <code>time2_accept_fallback = true</code>, <code>time2_posfix = 4</code> '
                      '(the dt-collapse fix: it passed the 107 d and 129.8 d collapse points of the previous run in restart tests) + the N<sup>2</sup> seed and MLT ramp; '
                      'started FRESH 2026-10-22 09:29 PDT. The previous run (deep l 1–8 seed, harmonic only) collapsed at 107–130 d and is archived.'),
]
kvh = ''.join(f'<div><b>{k}</b>{v}</div>' for k, v in kv)

TG_NOTE = ('T<sub>gas</sub> is inverted per cell from (ρ, e<sub>int</sub>) with the run’s own general EOS (H/He Saha, '
           'X 0.6538, Y 0.3264, no H<sub>2</sub>, no radiation term; numpy port of eos_composition.hpp on the code’s 1001 × 876 '
           'table grid); grey = cells at the 2000 K temperature floor.')


def _video1(fn, cap, poster=''):
    pa = f' poster="movies/{poster}"' if poster else ''
    return (f'<figure><video controls playsinline muted preload="metadata"{pa} src="movies/{fn}"></video>'
            f'<figcaption>{cap} (<a href="movies/{fn}">mp4</a>)</figcaption></figure>')


def video(fn, cap, poster=''):
    """movfix-1102 (USER): a 2-D movie is split into parts of <= ~14 MB (movies.json 'parts', oldest first): the newest part
    plays here, the earlier parts sit below it, newest first, in a collapsed list."""
    P = (MJ or {}).get('parts', {}).get(fn) if 'MJ' in globals() else None
    if not P:
        return _video1(fn, cap, poster)
    rng = lambda q: f"dumps {q['i0']}–{q['i1']}, t = {q['t0']:.1f}–{q['t1']:.1f} d, {q['n']} frames"  # noqa: E731
    k = len(P)
    h = _video1(P[-1]['file'], f"{cap} <br>Part {k} of {k} (newest; {rng(P[-1])}).")
    if k > 1:
        h += (f'<details><summary>Earlier parts of this movie ({k - 1}, newest first; every dump is one frame)</summary>'
              + ''.join(f'<figure><video controls playsinline muted preload="none" src="movies/{q["file"]}"></video>'
                        f'<figcaption>Part {j + 1} of {k}: {rng(q)} (<a href="movies/{q["file"]}">mp4</a>)</figcaption></figure>'
                        for j, q in reversed(list(enumerate(P[:-1]))))
              + '</details>')
    return h


def section3d(name):
    """optional 3-D view (10-21): movies/<name>_3d.mp4 (+ poster movies/<name>_3d_poster.png) when present in publish/movies"""
    mv = os.path.join(PUB, 'movies')
    if not os.path.exists(os.path.join(mv, name + '_3d.mp4')):
        return ''
    def one(f, cap):
        if not os.path.exists(os.path.join(mv, f + '.mp4')):
            return ''
        po = f + '_poster.png' if os.path.exists(os.path.join(mv, f + '_poster.png')) else ''
        return video(f + '.mp4', cap, po)
    return ('<section>\n<h2>3-D view</h2>\n<p>The star as a full sphere: the corrugated photosphere (the true &tau;<sub>R</sub> = 2/3 surface of every column), '
            'coloured by T<sub>eff</sub> and drawn at its true radius everywhere (no exaggeration; lightly smoothed over neighbouring columns); wedges are cut out and every cut face shows the radial velocity v<sub>r</sub>. Colour ranges are 1–99 % of the newest dump. '
            'Made by the 3-D movie pipeline (<code>/resnick/groups/carnegie_poc/jingze/vis3d</code>).</p>\n'
            + one(name + '_3d', 'Single cut: one 90° × 90° wedge (longitude 0–90°, latitude 0–90°) removed; the three cut faces show v<sub>r</sub> out to the local photosphere.')
            + '\n' + one(name + '_3d_layers', 'Layered: the photosphere cut by 90° × 180°, a deeper shell cut by 90° × 135° and a deeper one by 90° × 90°, so each layer shows through the one above it. The deeper shells show the relative temperature fluctuation δT/⟨T⟩ on that shell (T<sub>gas</sub> from the run’s EOS).')
            + '\n</section>\n')


if MJ:
    L = MJ['lim']
    movies = f"""
<section>
<h2>Slices and movies</h2>
<p>One frame per 1 d dump ({MJ['n']} frames, t = {MJ['t_first']:.0f}–{MJ['t_last']:.0f} d), every cubed-sphere cell drawn with its true edges (no resampling). The planar slices are the ring of panels 4, 0, 1, 2 at η = 0 (mean of the two cell rows next to x<sub>3</sub> = 0): 512 angles × the native radial cells. Colour scales are fixed per movie from percentiles pooled over every dump. {TG_NOTE}</p>
{video('rsgE_eq.mp4', f"Whole star, equatorial plane to ±{MJ['rmax']} R<sub>ph,IC</sub>: log ρ, v<sub>r</sub> (±{L['eq']['vr']:.3g} km/s = 99.5 % of |v<sub>r</sub>| inside R<sub>ph</sub>; faster gas outside saturates) and T<sub>gas</sub>.")}
{fig('fig_eq_latest.jpg', 'equatorial slice, latest dump', 'Still: the latest whole-star frame (t = %.1f d).' % MJ['t_last'])}
{video('rsgE_inner.mp4', f"Inner zoom r &lt; {MJ['zr']} R<sub>ph,IC</sub> (dashed: MESA CZ base 0.049 R<sub>ph</sub>), where the fastest flows are: log ρ − ⟨log ρ⟩<sub>shell</sub> (±{L['inner']['dlr']:.2g}), v<sub>r</sub> (±{L['inner']['vr']:.3g} km/s), δT<sub>gas</sub>/⟨T<sub>gas</sub>⟩ (±{L['inner']['dT']:.2g} %). The lateral cell merging (×8/×4/×2 inside 0.036/0.059/0.099 R<sub>ph</sub>) shows as rings.")}
{video('rsgE_phot.mp4', f"Photosphere band {MJ['pb'][0]}–{MJ['pb'][1]} R<sub>ph,IC</sub>, equatorial ring unwrapped (φ across, r up; every radial cell resolved): log ρ, v<sub>r</sub> (±{L['phot']['vr']:.3g} km/s) and δT<sub>gas</sub>/⟨T<sub>gas</sub>⟩ (±{L['phot']['dT']:.2g} %). Shows the infalling shell near 1.08 R<sub>ph</sub> and the outflow band just above R<sub>ph</sub>.")}
{video('rsgE_map.mp4', f"Mollweide maps, every 0.70° cell drawn, shells 0.90 and 0.96 R<sub>ph,IC</sub>: v<sub>r</sub> − ⟨v<sub>r</sub>⟩<sub>shell</sub> (±{L['map']['mvr_090']:.2g} / ±{L['map']['mvr_096']:.2g} km/s) and δT<sub>gas</sub>/⟨T<sub>gas</sub>⟩ (±{L['map']['mdT_090']:.2g} / ±{L['map']['mdT_096']:.2g} %). Longitude 0 is +x (centre of panel 0); graticule every 45°, so the vertical cube edges are the meridians at ±45°, ±135°.")}
</section>
"""
else:
    movies = ''
SEC3D = section3d('rsgE')


def section_ls():
    """11-01: Lomb-Scargle spectrum of the light curve, relaxed phase (mk_figs.py + lomb.py, every refresh)"""
    S = T.get('ls') or {}
    if not S or not S.get('peaks'):
        return ''
    D = S['det']
    fap = lambda p: '&lt; 10<sup>−10</sup>' if p < 1e-10 else f'{p:.1e}'
    pks = '; '.join(f'<b>{p["P"]:.0f} ± {p["dP"] / 2:.0f} d</b> (power {p["pow"]:.3f}, FAP {fap(p["fap"])})' for p in S['peaks'])
    p1 = S['peaks'][0]['P']
    relax = (f'<b>Relaxed phase: t ≥ {S["ts"]:.0f} d</b> — the first time (on a 50 d grid) after which the {S["wrun"]}-d running means of '
             f'L<sub>top</sub>, the median R<sub>ph</sub> and KE<sub>int</sub> all stay within {S["ksig"]:.0f}σ of their level in the second half of the run '
             f'(last excursions end at {D["L"]["last"]:.0f}, {D["R"]["last"]:.0f} and {D["KE"]["last"]:.0f} d: the MLT ramp, the R<sub>ph</sub> expansion and the 875–895 d ejection with its luminosity overshoot are excluded); '
             f'a linear fit to the 100-d block means over the window gives trend t-values {S["trend"][0]:.1f} / {S["trend"][1]:.1f} / {S["trend"][2]:.1f} (no significant secular trend). '
             f'Window {S["ts"]:.0f}–{S["tend"]:.0f} d = {S["T"]:.0f} d ≈ {S["T"] / p1:.1f} cycles of the dominant period.')
    cap = (f'Series: L<sub>top</sub>/L<sub>in</sub>, the radiative luminosity through the outer boundary at r = 1.761×10<sup>14</sup> cm (3 R<sub>ph,IC</sub>; the same series as the light curve in the History panel; '
           f'the light-travel time from R<sub>ph</sub> out to there is {S["tlt"]:.2f} d and the gas above R<sub>ph</sub> is optically thin, so delays do not matter for periods of tens–hundreds of days), '
           f'binned to daily means of the hst rows (prodG before 1675 d, prodH from 1675 d on: in the 1675–1895 d overlap only prodH is used), mean and a linear trend '
           f'({S["slope"]:+.4f} L<sub>in</sub> per 1000 d) removed; generalised Lomb–Scargle (floating mean, standard normalisation, numpy port of astropy’s), periods 20 d – T/2. '
           f'Top peaks: {pks}; ± = half the resolution P²/T. The 1 % false-alarm level ({S["fal"]:.3f}, Baluev 2008) assumes white noise; '
           f'the convective light curve is red, so it overstates significance. Blue: prodH alone (T = {S["Th"]:.0f} d), top peak {S["Ph1"]:.0f} d. '
           f'Expected radial fundamental for M = {S["M"]:.1f} M<sub>⊙</sub> (core point mass + mesh mass) and the relaxed mean R<sub>ph</sub> = {S["R"]:.0f} R<sub>⊙</sub>: '
           f'P = Q (M/M<sub>⊙</sub>)<sup>−1/2</sup>(R/R<sub>⊙</sub>)<sup>3/2</sup> = {S["Pq"][0]:.0f}–{S["Pq"][1]:.0f} d for Q = 0.04–0.08 d (green band); '
           f'Betelgeuse’s 416 d fundamental (Joyce et al. 2020) scaled to this M, R: {S["Pbg"]:.0f} d. The dominant {p1:.0f} d peak corresponds to Q = {p1 / S["M"] ** -0.5 / S["R"] ** 1.5:.3f} d.')
    return ('<section>\n<h2>Light-curve power spectrum</h2>\n'
            f'<p>{relax}</p>\n'
            + fig('fig_lombscargle.png', 'Lomb-Scargle power spectrum of the light curve', cap) + '\n</section>\n')


def section_sg():
    """11-01: supergranulation analysis (static, rsgE_supergran_1101; relaxed phase 1150-3061 d)"""
    D = '/resnick/groups/carnegie_poc/jingze/rsgE_supergran_1101/'
    if not os.path.exists(D + 'fig1_spectra.png'):
        return ''
    intro = ('<p><b>Is there supergranulation?</b> No distinct supergranular scale. The surface flow has deep, long-lived giant cells at '
             'l = 1–2 and, on top of them, one continuous spectrum of cells from l ≈ 6 down to the grid cutoff at l ≈ 45–60: no separate '
             'intermediate-scale hump in power, no break in lifetime and no separate divergence network, at any depth.</p>\n'
             '<div class="tbl"><table><tr><th>scale</th><th>size</th><th>lifetime</th></tr>'
             '<tr><td>giant cells, l ≈ 2</td><td>whole star (dipole/quadrupole), deep to surface</td><td>≈ 450–600 d (deep shells)</td></tr>'
             '<tr><td>surface cells, l ≈ 10</td><td>0.5–0.6 R ≈ 500–600 R<sub>⊙</sub></td><td>≈ 100 d</td></tr>'
             '<tr><td>smallest resolved, l ≈ 30–35</td><td>≈ 0.2 R ≈ 200 R<sub>⊙</sub></td><td>≈ 35 d</td></tr></table></div>\n'
             '<p>Method: real spherical-harmonic projection of the native cubed-sphere shell maps (v<sub>r</sub>, T, T<sub>eff</sub>, horizontal divergence) '
             'from 1912 daily reductions plus 40 raw dumps for 0.95–0.99 R<sub>ph</sub>; lifetimes = 1/e time of the per-l temporal autocorrelation. '
             'Compared with CO5BOLD RSG models (Freytag et al. 2002; Chiavassa et al. 2009, 2011): their large cells match our l ≈ 10 cells in size '
             '(ours live ≈ 100 d, not years; what lasts years here is the deep l = 1–2 flow). <b>Resolution caveat:</b> 128 cells per panel '
             '(≈ 12 R<sub>⊙</sub> at the surface); power rolls off from l ≈ 45–50, so a separate small-granule population (50–100 R<sub>⊙</sub>, l ≈ 60–125) '
             'cannot be tested; that needs 2–4× the angular resolution.</p>\n'
             + '<p><b>Blob-suppressed check.</b> Small, short-lived blobs could mask a larger pattern, as granulation masks solar supergranulation, so the search was repeated with 30–200 d time averages (pulsation removed), the divergence and vorticity of the mean flow, 600-d Lagrangian corks at 0.97 R<sub>ph</sub> (2000–2599 d), sign(v<sub>r</sub>) spectra and watershed cell segmentation. All recover the same single population: cells ≈ 0.4 R<sub>ph</sub> ≈ 400 R<sub>⊙</sub> across (l ≈ 10–15, ≈ 15–20 H<sub>p</sub> at the surface) living ≈ 100–200 d, on the giant l ≤ 4 flow (&gt; 1 yr). The corks never stall at an intermediate scale (network 0.23 → 0.36 → 0.75 R<sub>ph</sub> after 15 / 100 / 200 d) and end in 4–8 sinks ≈ 0.6–0.8 R<sub>ph</sub> apart; the cell-size distribution is single-peaked and barely changes with averaging (0.38 → 0.42 R<sub>ph</sub>). The blobs are substructure inside these cells, not a mask over a larger pattern.</p>\n'
             + '<p><b>Granule-like (localized) decomposition.</b> Global Y<sub>lm</sub> modes do not look like cells, so the maps were also split with dyadic spherical needlets (localized wavelets, l ≈ 2<sup>j</sup>) and a watershed catalogue of individual cells (v<sub>r</sub> at 0.97 R<sub>ph</sub>, 600 daily maps). Every band looks like granules, but it is one self-similar, continuously fragmenting hierarchy: blob size just follows the band, the energy per scale is a broad plateau with no separate peak, lifetime falls smoothly with scale, and the cell catalogue is single-peaked with size, contrast and lifetime changing continuously (3 d below 0.15 R<sub>ph</sub>, 69 d at 0.3–0.6 R<sub>ph</sub>, ≈ 200 d above 0.6 R<sub>ph</sub>). The dominant structures are hot 0.3–0.5 R<sub>ph</sub> upflow cells (l ≈ 8–16, 60–110 d); no second (supergranular or granular) population appears at l ≈ 4–64.</p>\n')
    figs = [('fig1_spectra.png', 'Angular power spectra', 'Compensated angular power spectra of v<sub>r</sub> at 9 depths, and of T<sub>eff</sub>, the R<sub>ph</sub> corrugation and the horizontal divergence: flat from l ≈ 6 to ≈ 45, no separate hump.'),
            ('fig2_lifetimes.png', 'Lifetimes vs degree', 'Lifetime vs degree l and band autocorrelations: a smooth power law τ ∝ l<sup>−0.8</sup>, no break. Near the surface the low-l autocorrelation oscillates with the ≈ 500 d pulsation.'),
            ('fig3_maps_div.jpg', 'Maps at 2400 d', 'Mollweide maps at t = 2400 d: T<sub>eff</sub>, v<sub>r</sub> and the native horizontal divergence at 0.97 R<sub>ph</sub>, and the divergence filtered to l = 1–4, 5–16 and 17–64. The network outlines only the smallest resolved cells.'),
            ('fig4_cellsize_depth.png', 'Cell size vs depth', 'Angular autocorrelation and cell size vs depth, against the pressure scale height and the resolution limit.'),
            ('fig5_timeavg_spectra.png', 'Time-averaged spectra', '<b>Blob-suppressed search.</b> Spectra of 30–200 d time-averaged v<sub>r</sub>, T, horizontal divergence, vorticity and sign(v<sub>r</sub>) at 0.97 R<sub>ph</sub> (pulsation l ≤ 2 removed), plus 1912-d averages of v<sub>r</sub> at 0.90 R<sub>ph</sub> and T<sub>eff</sub>: averaging removes the l ≳ 30 blobs and leaves the same plateau at l ≈ 6–30, no new hump.'),
            ('fig6_timeavg_corks_maps.jpg', 'Mean-flow maps and corks', 'Daily v<sub>r</sub> vs the 100-d mean with watershed cells, the 100-d mean divergence and vorticity, and 300 000 corks advected in the horizontal flow after 60 d (network cells) and 300 d (sinks at the giant-cell downflow junctions).'),
            ('fig7_cellsizes_persistence.png', 'Cell sizes and persistence', 'Watershed and cork-network cell-size distributions, and persistence of the window-mean pattern vs l: one single-peaked population of ≈ 0.4 R<sub>ph</sub> cells living ≈ 100–200 d; the cork network coarsens continuously without stalling at an intermediate scale.'),
            ('fig8_needlet_snapshot.jpg', 'Needlet decomposition at 2400 d', '<b>Granule-like (localized) decomposition.</b> T<sub>eff</sub> and v<sub>r</sub> at 0.97 R<sub>ph</sub> at t = 2400 d, full field and split into dyadic spherical needlet scales j = 1–6 (l ≈ 2<sup>j</sup>): every band looks like granules, with blob size following the band.'),
            ('fig9_needlet_spectrum_blobs.png', 'Needlet energy and blobs', 'Needlet energy per scale (a broad plateau at l ≈ 8–32, no separate bump), blob sizes per scale (0.75 → 0.07 R<sub>ph</sub> for l ≈ 4–64) and band lifetimes (≈ 110 → 21 d), falling smoothly with scale.'),
            ('fig10_granule_catalogue.png', 'Cell catalogue', 'Watershed catalogue of individual upflow cells (413 per day, 248k in total, tracked by overlap): one skewed size distribution (median 0.15 R<sub>ph</sub>, area-weighted peak 0.27 R<sub>ph</sub>); lifetime and v<sub>r</sub> / T contrast rise continuously with size; ≈ 70 merges and 70 splits per day.')]
    return ('<section>\n<h2>Supergranulation</h2>\n' + intro
            + '\n'.join(fig(D + f, a, c) for f, a, c in figs) + '\n</section>\n')


def section_modes():
    # USER 11-06: mode analysis (rsgmode_1106, static snapshot to 6269 d; files in the page's modes/, not inlined)
    figs = [('vr_spacetime.png', 'Shell-averaged radial velocity vs radius and time', 'Shell-averaged v<sub>r</sub>(r, t). Top: mass-weighted &lt;&rho; v<sub>r</sub>&gt;/&lt;&rho;&gt;. Middle: divided by the linear F-mode shape. Bottom: area-weighted, biased outward because upflows cover most of the area. Black line: median R<sub>ph</sub>. Dotted lines: code or key changes.'),
            ('modal_energy.png', 'Modal kinetic energy vs time', 'Peak KE in each linear mode band, F share of F + 1O, F phase against a 456 d clock, and mass-mean v<sub>r</sub> inside 0.9 R<sub>0</sub>. The drop in the last ~300 d is a filter edge effect, not real.'),
            ('radial_structure.png', 'Radial structure of the modes', 'Amplitude and phase against radius. F has flat phase with no node; the 211-228 d signal flips by about 180&deg; near 0.9-1.0 R<sub>0</sub>. Ignore the "348 d" 1O-band peaks, which are F leakage.'),
            ('filtered.png', 'Band-filtered time series', 'Band-filtered velocities.'),
            ('scalogram.png', 'Scalogram', 'Time-period power. Below about 1100 d the star&apos;s expansion contaminates the F band.'),
            ('pod.png', 'POD modes', 'Mass-weighted POD for t &gt; 1100 d: POD0 83 % (0 nodes, 460 d), POD1 7 % (1 node, 227-251 d). POD2 (6 %, 172 d) has its power mostly near and above R<sub>ph</sub>.')]
    tbl = ('<div class="tbl"><table><tr><th>t [d]</th><th>F period [d]</th><th>v<sub>F</sub> at R<sub>ph</sub> [km/s]</th><th>KE F / 1O [10<sup>43</sup> erg]</th><th>Dominant</th></tr>'
           '<tr><td>1100-2056</td><td>400-530, phase poorly defined</td><td>2.9 to 1.0</td><td>3 to 0.45 / 0.1-0.36</td><td>F at the surface; ~215 d in deeper mass-mean v</td></tr>'
           '<tr><td>2056-2396</td><td>nearly gone</td><td>~0.6</td><td>0.16 / 0.40</td><td>1O</td></tr>'
           '<tr><td>2400-2800</td><td>growing</td><td>0.6 to 3.7</td><td>0.16 to 5.3 / 0.5 to 0.17</td><td>transition</td></tr>'
           '<tr><td>2800-4500</td><td>455.5, coherent</td><td>3.3-4.4</td><td>4-7 / 0.07-0.48</td><td>F</td></tr>'
           '<tr><td>4500-5100</td><td>weakened</td><td>1.7-2.2</td><td>1.2-1.9 / ~0.5</td><td>F, narrowly</td></tr>'
           '<tr><td>5100-6269</td><td>~470</td><td>3.2 to 6.9</td><td>4-18 / 0.04-0.33</td><td>F</td></tr></table></div>')
    intro = ('<p>Analysis of the dumps to 6269 d (snapshot from 2026-11-06, not updated with each refresh). '
             'A one-node ~215-245 d first overtone (1O) led around 2056-2396 d. A nodeless ~455 d fundamental (F) took over by about 2570 d and has dominated since, with one dip near 4850 d. '
             'It looks more like F fading and regrowing over an always-present 1O than a sharp switch: F&apos;s phase slips ~0.4 cycles near 2240 d, when its amplitude is near zero. '
             'Linear adiabatic periods from the mean profile: F 477 d (431-482 d with the outer boundary), 1O 261 d (216-273 d). '
             'No code change lines up with the 2570 d turn, which suggests a physical origin, but no control run exists. F is still growing, so saturation is not established. '
             'Total KE<sub>r</sub> is about 10<sup>46</sup> erg, so most radial motion is not modal.</p>')
    # USER 11-06: nu-r and l-nu diagrams (rsgmode_1106/lnu, scripts in hgmode_1106/scripts; static snapshot, files in modes/)
    figs2 = [('nur_rsg.png', 'Power against frequency and radius', 'Mass-weighted shell-mean v<sub>r</sub> to 6269 d. Top: each radius normalized by its own maximum. Bottom: weighted by &rho;r<sup>2</sup>. Only 2800-4500 d resolves F from 1O (1/T = 0.00059 /d). In 2056-2396 d (T = 340 d, 1/T = 0.00294 /d) F and 1O merge into one blob; the dark line near 0.92 R<sub>0</sub> is a phase flip, consistent with the one-node 1O.'),
             ('lnu_rsg_B.png', 'Power against degree and frequency, 2800-4500 d', 'Shell maps of v<sub>r</sub>, 2800-4500 d, about 3.7 F cycles. Convection dominates at all l: a smooth dome, no ridge.'),
             ('lnu_rsg_A.png', 'Power against degree and frequency, 2056-2396 d', 'Same for 2056-2396 d, 0.7 F cycles: the F frequency lies below 1/T, so nothing is resolved near F or 1O here. Shown for completeness.'),
             ('lnu_rsg_check.png', 'Split-window check', 'Spectra at l = 0, 1, 2, 3, 4, 6 in 2800-4500 d, its two halves and the later 5100-6292 d. The l = 0 F peak and the l = 1 peak near it repeat in every window at 0.6, 0.9 and 1.085 R<sub>0</sub>; l &ge; 2 shows no repeating peak of comparable contrast.'),
             ('lnu_rsg_l1.png', 'The l = 1 signal against l = 0', 'l = 0, 1, 2 over 2800-6292 d (1/T = 0.00029 /d), then F-band envelopes (330-700 d filter), dipole axis and dipole phase against a<sub>00</sub>. The filter smears about a period and the ends are unreliable.')]
    sub2 = ('<h3>Frequency diagrams: radius and spherical harmonic degree</h3>\n'
            '<p>Static snapshot from 2026-11-06: nu-r to 6269 d, l-nu to 6292 d; not updated with each refresh. Native shell maps at 0.30, 0.60, 0.90, 0.98 and 1.085 R<sub>0</sub> (R<sub>ph</sub> &asymp; 1.20 R<sub>0</sub>, so the top shell is 0.90 R<sub>ph</sub>), real Y<sub>lm</sub> to l = 100, daily dumps resampled to 1 d, line removed, Hann window. Nyquist 0.5 /d; power there is about 10<sup>-8</sup> of the peak, so aliasing is negligible. Three stale reductions (dumps 3721-3723) are dropped and bridged. The l = 0 row is area-weighted, unlike the mass-weighted nu-r plot. The roll-off marked on the plots (l E<sub>l</sub> below half its l = 6-30 plateau) is l &asymp; 52, 61, 76, 89 at 0.30, 0.60, 0.90, 1.085 R<sub>0</sub>.</p>\n'
            '<p><b>Measured.</b> l = 0: F at 454-469 d at all shells, in both halves of 2800-4500 d and in 5100-6292 d. l = 1: a peak at 530-578 d in 2800-4500, 5100-6292 and 2800-6292 d (545-648 d in the coarser halves), separate from F over 2800-6292 d (458 d against 553 d, difference 1.3/T, just under the Hann width 1.44/T), with an amplitude comparable to or above l = 0 at 0.6-1.085 R<sub>0</sub>. It keeps going through 2000-2300 d when F nearly vanishes, its phase drifts steadily against a<sub>00</sub>, and it is close to linearly polarized (linear fraction about 0.9) with an axis that turns slowly, about 40&deg; over 4000 d. It rises outward with no phase flip from 0.6 to 1.085 R<sub>0</sub>. Weaker l = 1 peaks repeat near 250 d and 165 d. No ridge at l &ge; 2 at any radius: a persistence test on l &le; 12 finds only marginal, scattered peaks.</p>\n'
            '<p><b>Inferred, not established.</b> The l = 1 signal behaves like its own standing oscillation, not F modulated by asphericity. Its identity is open: an l = 1 acoustic mode, or the envelope sloshing about the fixed central gravity (the code&apos;s gravity is a frozen monopole, so that motion need not exist in a real star). A third option is not ruled out: the deep l = 1-2 giant cells, whose lifetimes of about 450-600 d (supergranulation section) overlap this period. No non-radial eigenmodes were computed and the momentum budget was not checked. l = 0 at 136 d repeats across shells and windows and is near the linear 3O (134 d), but the identification rests on that match alone. LAWE values are linear adiabatic estimates from a shell mean. The maps carry no rotation removal.</p>\n')
    return ('<section>\n<h2>Radial pulsation modes</h2>\n' + intro + tbl + '\n'
            + '\n'.join(f'<figure><img alt="{html.escape(a)}" loading="lazy" src="modes/{f}"><figcaption>{c}</figcaption></figure>' for f, a, c in figs) + '\n' + sub2
            + '\n'.join(f'<figure><img alt="{html.escape(a)}" loading="lazy" src="modes/{f}"><figcaption>{c}</figcaption></figure>' for f, a, c in figs2) + '\n</section>\n')


SECSG = section_sg()
SECMODES = section_modes()
SECLS = section_ls()
NOTES = open(os.path.join(HERE, 'notes.html')).read() if os.path.exists(os.path.join(HERE, 'notes.html')) else ''

page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>RSG E Caltech</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Source+Sans+3:wght@400;600&family=Source+Code+Pro:wght@400;500&display=swap">
<style>{CSS}</style></head><body>
<main>
<header>
<h1>Red supergiant g25, run E on Caltech</h1>
<p class="lede">AthenaK radiation hydrodynamics of the full convective envelope of a red supergiant on a global cubed-sphere mesh,
from 0.03 R<sub>ph</sub> (inside the radiative layer below the convection zone) to 3 R<sub>ph</sub>. This is the Caltech copy of RSG E (now the only one running),
restarted from t = 0 on 10-22 with the dt-collapse fix, an N<sup>2</sup>-placed surface seed and a growth-rate MLT ramp; it runs in short debug-QOS links on 1 node × 4 H200 until a 7-d allocation starts. This page reads the run’s history, logs and 1 d dumps.</p>
<div class="kv">{kvh}</div>
<p class="status">{html.escape(STATUS)}</p>
</header>

<section>
<h2>History</h2>
{fig('fig_hist.png', 'history panels', 'From <code>rsg.user.hst</code> (every ~cycle; rows rewritten after a restart are dropped). L<sub>in</sub> is the luminosity imposed at the inner boundary; L<sub>bot</sub>, L<sub>mid</sub> and L<sub>top</sub> are the code’s radiative luminosities at the bottom, middle and top of the star region; L<sub>MLT</sub> is the sub-grid convective luminosity of the adaptive closure. The energy panel compares the change of the total energy E<sub>tot</sub> with the time integral of L<sub>in</sub> − L<sub>top</sub>; this ignores the flux between the top of the star region and the outer boundary and the sponge / source terms, so it is a consistency check, not a strict conservation test. R<sub>ph</sub> and T<sub>eff</sub> come from the dumps: per column, R<sub>ph</sub> is where σT<sub>rad</sub><sup>4</sup>r<sup>2</sup> equals the column’s outgoing r<sup>2</sup>F<sub>r</sub> at 1.2–1.4 R<sub>ph</sub> (grey Eddington).')}
</section>

{SECLS}
{SECMODES}
{SECSG}

<section>
<h2>Radial profiles</h2>
{fig('fig_profiles.png', 'radial profiles', 'Solid-angle-weighted shell means over the full sphere at a few dump times (t = 0 is the first dump, after the in-code IC relax). Temperature: T<sub>gas</sub> from the run’s EOS (solid) and T<sub>rad</sub> = (E/a)<sup>1/4</sup> (dashed). Dotted white: the IC column (ρ, T<sub>gas</sub>, scaffold F<sub>MLT</sub>/F). Velocities are ρ-weighted rms; the tangential speed uses the cubed-sphere’s non-orthogonal unit basis (e<sub>1</sub>·e<sub>2</sub> = −sin ξ sin η). The resolved convective luminosity is 4πr²(⟨ρv<sub>r</sub>h⟩ − ⟨ρv<sub>r</sub>⟩⟨h⟩) with the total enthalpy h = (e<sub>int</sub> + P<sub>gas</sub> + E + P<sub>rr</sub>)/ρ (gas terms from the run’s EOS table incl. ionisation; P<sub>rr</sub> from the variable Eddington tensor on the dumps that carry it, E/3 otherwise); dashed: the kinetic-energy flux 4πr²(⟨ρv<sub>r</sub>v²/2⟩ − ⟨ρv<sub>r</sub>⟩⟨v²/2⟩).')}
{fig('fig_fluxbudget.png', 'energy-flux budget', 'Shell-integrated luminosities in units of L<sub>in</sub>, averaged over up to 10 of the dumps that carry the code’s Eddington tensor (m1_vet, every 5th dump) at mid-run (left) and at the latest times (right): comoving radiative L<sub>rad</sub> = 4πr²⟨F<sub>r</sub> − (E + P<sub>rr</sub>)v<sub>r</sub>⟩ (the lab-frame M1 flux minus the advected radiation enthalpy, as the code’s own split A = vE + v·P, which is counted in the enthalpy flux instead; P<sub>rr</sub> = D<sub>rr</sub>E from the run’s variable Eddington tensor; the off-diagonal P<sub>rθ</sub>v<sub>θ</sub> + P<sub>rφ</sub>v<sub>φ</sub> is not dumped), resolved enthalpy flux with the total h = (e<sub>int</sub> + P<sub>gas</sub> + E + P<sub>rr</sub>)/ρ (dashed orange: gas part only), kinetic-energy flux, and their sum (white; 1 in a thermally relaxed state). The enthalpy and kinetic fluxes are the parts carried by the fluctuations about the shell-mean mass flux; what the mean radial flow (breathing / pulsation) carries, 4πr²⟨ρv<sub>r</sub>⟩⟨h + v²/2⟩ including the radiation enthalpy, is shown separately (grey dashed), and the dotted white curve adds it back (the lab-frame total). Dotted: the old ideal-gas enthalpy flux h = (5/3)e<sub>int</sub>/ρ. Cell-centred F<sub>r</sub> is noisy (odd–even) in the very optically thick deep shells; the history’s face-based L<sub>bot</sub> is the reliable number there.')}
</section>

{movies}

{SEC3D}

<section>
<h2>Slices</h2>
<p>Static slices of the latest dump on native cells (PNG, 120 dpi). Rows 1–2: equatorial (z = 0, panels 4,0,1,2 at η = 0) and meridional (y = 0, panels 5,0,3,2 at ξ = 0) rings through the whole star and the CSM out to 2.5 R<sub>ph,IC</sub> (sponge beyond 2.72), colour ranges from 0.25–1 R<sub>ph</sub>. Row 3: zoom r &lt; 0.25 R<sub>ph</sub> with its own ranges (dashed: MESA CZ base). Columns: log ρ; log ρ − ⟨log ρ⟩<sub>shell</sub>; T<sub>gas</sub> (grey = 2000 K floor); v<sub>r</sub>. Below: the photosphere band unwrapped (log ρ, v<sub>r</sub>, T<sub>gas</sub>, δT<sub>gas</sub>/⟨T<sub>gas</sub>⟩).</p>
{slices}
{phot}
</section>

<section>
<h2>Maps</h2>
<p>Whole-sphere Mollweide maps of the latest dump, each 0.70° cubed-sphere cell drawn (every screen pixel takes its own cell): v<sub>r</sub> and δT<sub>gas</sub>/⟨T<sub>gas</sub>⟩ on the shells 0.11 (fastest deep flows), 0.3, 0.9, 0.96 (below the sharp photospheric T front at ~0.975) and 1.085 R<sub>ph,IC</sub> (the infalling shell); a shell that moves as a whole is shown as v<sub>r</sub> − ⟨v<sub>r</sub>⟩, with the mean in the label. Also R<sub>ph</sub>/R<sub>ph,IC</sub> and T<sub>eff</sub> per column (grey Eddington, from T<sub>rad</sub> and the flux).</p>
{maps}
</section>

{NOTES}


<footer>Run: Caltech {html.escape(os.environ.get('RUN', ''))} (READ ONLY). Page rebuilt by /resnick/groups/carnegie_poc/jingze/rsgE_page/update.sh (reduce.py → mk_figs.py → mk_movies.py → build_page.py; base: Raven RSGPAGE 10-19, athenak c542d6b0 docs/handover/rsg-page). Movies are separate files under movies/.</footer>
</main>
</body></html>
'''
open(OUTF, 'w').write(page)
print(OUTF, os.path.getsize(OUTF))
