#!/usr/bin/env python3
"""hgglob page 10-21 (from RSGPAGE 10-19): assemble publish/index.html for the He giant GLOBAL cs run (figures inlined base64, one file;
movies are separate files under publish/movies, incl. the optional 3-D movie hggl_3d.mp4 + hggl_3d_poster.png).
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
SITE = os.environ.get('SITE', 'Caltech')
LAYOUT = os.environ.get('LAYOUT', '2 nodes × 4 H200 = 8 ranks, 24 MeshBlocks 704×64×64 (3 per GPU); 7 d links')
BASE = os.environ.get('BASENAME', 'hegiant')
TFL = float(os.environ.get('TFLOOR', '2400'))
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
         + f' (cell-stages in the 80 cycles before cycle {EV.get("cycle", "?")}, <code>{BASE}.log</code>)') if EV else f'{BASE}.log not read'

kv = [
    ('status', f'{html.escape(SITE)}: {len(running)} link running, {len(queued)} queued — {html.escape(jobtxt)}'),
    ('now', f't = {T["t_s"]:.4g} s = <b>{T["t_d"]:.2f} d</b>, cycle {T["cyc"]}, dt {T["dt"]:.3g} s, {T["spc"]:.2f} s/cycle '
            f'(≈ {86400 / max(T["dt"], 1e-30) * T["spc"] / 3600:.2f} h wall per sim day); Picard mean {T["pic"]:.2f}, max {T["picmax"]:.0f}'),
    ('health', health),
    ('floors / limiters', evtxt + (f'; T<sub>gas</sub> at the {TFL:.0f} K floor in the last dump: {T["nfloorT"]:,} cells '
                                   f'({100 * T["nfloorT"] / max(T["ncell"], 1):.1f} % of the mesh)' if T.get('nfloorT', -1) >= 0 else '')),
    ('last dump', f'{T["t_dump"]:.1f} d ({T["ndump"]} dumps reduced): median R<sub>ph</sub> {T["rph"]:.4f} R<sub>ph,IC</sub>, '
                  f'T<sub>eff</sub> {T["teff"]:.0f} K, max |v<sub>r</sub>| {T["vmax"]:.1f} km/s'),
    ('budget', f'L<sub>top</sub>/L<sub>in</sub> {T["LL"]:.4f}, L<sub>bot</sub>/L<sub>in</sub> {T["Lbot"]:.4f}; '
               f'ΔE<sub>tot</sub>/|E<sub>0</sub>| {T["dE"]:.2e} vs ∫(L<sub>in</sub>−L<sub>top</sub>)dt/|E<sub>0</sub>| {T["dEb"]:.2e}; '
               f'ΔM/M {T["dM"]:.2e}; KE<sub>int</sub> {T["ke"]:.3e} erg'),
    ('star', 'He giant from the MESA model 14p5M_widebinary_profile14 (S. Wu; He star 2.754 M<sub>⊙</sub>): IC column '
             '<code>ic_hegiant_14p5M_mk.txt</code> (MESA T(P) re-integrated in HSE with the code’s EOS and gravity) with '
             'R<sub>ph</sub> 7.947×10<sup>12</sup> cm = 114.2 R<sub>⊙</sub>, L 1.421×10<sup>38</sup> erg/s = 3.71×10<sup>4</sup> L<sub>⊙</sub> (imposed at the inner wall: L<sub>in</sub> 1.404×10<sup>38</sup> erg/s = 3.67×10<sup>4</sup> L<sub>⊙</sub>), '
             'T<sub>eff</sub> 7.50 kK; X 0.00095, Y 0.9775, Z 0.0215; convective (scaffold F<sub>MLT</sub>/F &gt; 1 %) from the inner wall to ≈ R<sub>ph</sub>'),
    ('grid', 'Global cubed sphere (gnomonic equiangular), N = 128 per face (0.70°), 6 panels; 704 stretched radial cells '
             'r = 5.566×10<sup>11</sup> – 2.474×10<sup>13</sup> cm (0.070 – 3.11 R<sub>ph</sub>; radial grid fit145b_704, fine band at the photosphere); '
             'no lateral cell merging; 69 M cells'),
    ('layout', html.escape(LAYOUT)),
    ('physics', 'Implicit grey M1 radiation with the VET closure (vet_col + vet_gd, half-range flux blend), hesdirk2; '
                'general EOS (H/He ionisation, metals inert with a<sub>metal</sub> 1.80, no H<sub>2</sub>, no radiation term), '
                'MESA kap Rosseland + ext2 Planck tables; adaptive MLT closure (α<sub>MLT</sub> 2) with the automatic growth-rate ramp (start 0, length t<sub>sat</sub> = 10.5/(0.61 σ<sub>max</sub>) = 13.3 d from the IC N<sup>2</sup> on the code cells; MLT off at 13.3 d); '
                'angular-momentum transport on; shared seam flux; inner wall closed with the luminosity imposed as a flux, outer outflow; '
                'floors ρ 10<sup>−18</sup>, T 2400 K, sponge (rate 10<sup>−3</sup>) on ρ &lt; 2×10<sup>−14</sup> at r &gt; 2.23×10<sup>13</sup> cm (2.80 R<sub>ph</sub>); '
                'dt_retry (cut 0.9, regrow 1.1 after 20 good steps; non-admissible hesdirk2 stage → retry)'),
    ('seed', 'Automatic N<sup>2</sup> seed (he_seed_band = n2): real spherical harmonics Y<sub>lm</sub>, l = 16–64 (3969 terms, all m; window l<sub>hi</sub> = min(π r/H<sub>p</sub>, 8 cells per wavelength), l<sub>lo</sub> = l<sub>hi</sub>/4), '
             'amplitude 10<sup>−3</sup>, rng 22102, radial band 7.753–7.929×10<sup>12</sup> cm = 0.976–0.998 R<sub>ph</sub> where the IC growth rate is ≥ 0.3 σ<sub>max</sub> (1/σ<sub>max</sub> = 0.78 d)'),
    ('binary / deck', 'fin51-r-fix d01ce628c <code>athena_gpu_he_star_m1_d01ce628_nofma</code> md5 83b6545e (incremental build); deck <code>hgglob_1022_n2seed.athinput</code> = the production deck + '
                      '<code>implicit_face_opac = harmonic</code>, <code>time2_vet_col = extrap2</code>, <code>time2_accept_fallback = true</code>, <code>time2_posfix = 4</code> (the dt-collapse fix) + '
                      'the N<sup>2</sup> seed, the automatic MLT ramp and the GPU-memory keys (dt_retry_host, time2_posfix_host, restart_pack_per_mb); started FRESH 2026-10-22 15:32 PDT. '
                      'The previous run (deep l 1–8 seed, harmonic only) collapsed at 18.8–20 d and is archived.'),
]
kvh = ''.join(f'<div><b>{k}</b>{v}</div>' for k, v in kv)

TG_NOTE = ('T<sub>gas</sub> is inverted per cell from (ρ, e<sub>int</sub>) with the run’s own general EOS (H/He Saha, '
           'X 0.00095, Y 0.9775, Z 0.0215 inert with a<sub>metal</sub> 1.80, no H<sub>2</sub>, no radiation term; numpy port of eos_composition.hpp on the code’s 1001 × 876 '
           'table grid; it reproduces the IC column’s T to 1.5×10<sup>−4</sup>); grey = cells at the %.0f K temperature floor.' % TFL)


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
    """optional 3-D view: movies/<name>_3d.mp4 (+ poster movies/<name>_3d_poster.png) when present in publish/movies"""
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
<p>One frame per 1 d dump ({MJ['n']} frame{'s' if MJ['n'] != 1 else ''}, t = {MJ['t_first']:.0f}–{MJ['t_last']:.0f} d), every cubed-sphere cell drawn with its true edges (no resampling). The planar slices are the ring of panels 4, 0, 1, 2 at η = 0 (mean of the two cell rows next to x<sub>3</sub> = 0): 512 angles × the native radial cells. Colour scales are fixed per movie from percentiles pooled over every dump. {TG_NOTE}</p>
{video('hggl_eq.mp4', f"Whole star, equatorial plane to ±{MJ['rmax']} R<sub>ph,IC</sub>: log ρ, v<sub>r</sub> (±{L['eq']['vr']:.3g} km/s = 99.5 % of |v<sub>r</sub>| inside R<sub>ph</sub>; faster gas outside saturates) and T<sub>gas</sub>.")}
{fig('fig_eq_latest.jpg', 'equatorial slice, latest dump', 'Still: the latest whole-star frame (t = %.1f d).' % MJ['t_last'])}
{video('hggl_inner.mp4', f"Inner zoom r &lt; {MJ['zr']} R<sub>ph,IC</sub> (dashed: the inner wall r<sub>in</sub> = 0.070 R<sub>ph</sub>): log ρ − ⟨log ρ⟩<sub>shell</sub> (±{L['inner']['dlr']:.2g}), v<sub>r</sub> (±{L['inner']['vr']:.3g} km/s), δT<sub>gas</sub>/⟨T<sub>gas</sub>⟩ (±{L['inner']['dT']:.2g} %).")}
{video('hggl_phot.mp4', f"Photosphere band {MJ['pb'][0]}–{MJ['pb'][1]} R<sub>ph,IC</sub>, equatorial ring unwrapped (φ across, r up; every radial cell resolved): log ρ, v<sub>r</sub> (±{L['phot']['vr']:.3g} km/s) and δT<sub>gas</sub>/⟨T<sub>gas</sub>⟩ (±{L['phot']['dT']:.2g} %).")}
{video('hggl_map.mp4', f"Mollweide maps, every 0.70° cell drawn, shells 0.90 and 0.98 R<sub>ph,IC</sub>: v<sub>r</sub> − ⟨v<sub>r</sub>⟩<sub>shell</sub> (±{L['map']['mvr_090']:.2g} / ±{L['map']['mvr_098']:.2g} km/s) and δT<sub>gas</sub>/⟨T<sub>gas</sub>⟩ (±{L['map']['mdT_090']:.2g} / ±{L['map']['mdT_098']:.2g} %). Longitude 0 is +x (centre of panel 0); graticule every 45°, so the vertical cube edges are the meridians at ±45°, ±135°.")}
</section>
"""
else:
    movies = ''
SEC3D = section3d('hggl')
NOTES = open(os.path.join(HERE, 'notes.html')).read() if os.path.exists(os.path.join(HERE, 'notes.html')) else ''

page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>He Giant Global Caltech</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Source+Sans+3:wght@400;600&family=Source+Code+Pro:wght@400;500&display=swap">
<style>{CSS}</style></head><body>
<main>
<header>
<h1>He giant, global run on Caltech</h1>
<p class="lede">AthenaK radiation hydrodynamics of the convective envelope of a 14.5 M<sub>⊙</sub>-progenitor He giant (2.75 M<sub>⊙</sub> He star, R<sub>ph</sub> 114 R<sub>⊙</sub>)
on a global cubed-sphere mesh (full sphere), from a closed inner wall at 0.07 R<sub>ph</sub> to 3.1 R<sub>ph</sub>. Restarted from t = 0 on 10-22 with the dt-collapse fix, an automatic N<sup>2</sup>-placed surface seed and a growth-rate MLT ramp; it runs in debug-QOS links on 1 node until a larger allocation starts.
This page reads the run’s history, logs and 1 d dumps.</p>
<div class="kv">{kvh}</div>
<p class="status">{html.escape(STATUS)}</p>
</header>

<section>
<h2>History</h2>
{fig('fig_hist.png', 'history panels', f'From <code>{BASE}.user.hst</code> (every ~cycle; rows rewritten after a restart are dropped). L<sub>in</sub> is the luminosity imposed at the inner boundary; L<sub>bot</sub>, L<sub>mid</sub> and L<sub>top</sub> are the code’s radiative luminosities at the bottom, middle and top of the star region; L<sub>MLT</sub> is the sub-grid convective luminosity of the adaptive closure. The energy panel compares the change of the total energy E<sub>tot</sub> with the time integral of L<sub>in</sub> − L<sub>top</sub>; this ignores the flux between the top of the star region and the outer boundary and the sponge / source terms, so it is a consistency check, not a strict conservation test. R<sub>ph</sub> and T<sub>eff</sub> come from the dumps: per column, R<sub>ph</sub> is where σT<sub>rad</sub><sup>4</sup>r<sup>2</sup> equals the column’s outgoing r<sup>2</sup>F<sub>r</sub> at 1.2–1.4 R<sub>ph</sub> (grey Eddington).')}
</section>

<section>
<h2>Radial profiles</h2>
{fig('fig_profiles.png', 'radial profiles', 'Solid-angle-weighted shell means over the full sphere at a few dump times (t = 0 is the first dump, after the in-code IC relax). Temperature: T<sub>gas</sub> from the run’s EOS (solid) and T<sub>rad</sub> = (E/a)<sup>1/4</sup> (dashed). Dotted white: the IC column (ρ, T<sub>gas</sub>, scaffold F<sub>MLT</sub>/F). Velocities are ρ-weighted rms; the tangential speed uses the cubed-sphere’s non-orthogonal unit basis (e<sub>1</sub>·e<sub>2</sub> = −sin ξ sin η). The resolved convective luminosity is 4πr²(⟨ρv<sub>r</sub>h⟩ − ⟨ρv<sub>r</sub>⟩⟨h⟩) with the total enthalpy h = (e<sub>int</sub> + P<sub>gas</sub> + E + P<sub>rr</sub>)/ρ (gas terms from the run’s EOS table incl. ionisation; P<sub>rr</sub> from the variable Eddington tensor on the dumps that carry it, E/3 otherwise); dashed: the kinetic-energy flux 4πr²(⟨ρv<sub>r</sub>v²/2⟩ − ⟨ρv<sub>r</sub>⟩⟨v²/2⟩).')}
{fig('fig_fluxbudget.png', 'energy-flux budget', 'Shell-integrated luminosities in units of L<sub>in</sub>, averaged over up to 10 of the dumps that carry the code’s Eddington tensor (m1_vet, every 5th dump) at mid-run (left) and at the latest times (right): comoving radiative L<sub>rad</sub> = 4πr²⟨F<sub>r</sub> − (E + P<sub>rr</sub>)v<sub>r</sub>⟩ (the lab-frame M1 flux minus the advected radiation enthalpy, as the code’s own split A = vE + v·P, which is counted in the enthalpy flux instead; P<sub>rr</sub> = D<sub>rr</sub>E from the run’s variable Eddington tensor; the off-diagonal P<sub>rθ</sub>v<sub>θ</sub> + P<sub>rφ</sub>v<sub>φ</sub> is not dumped), resolved enthalpy flux with the total h = (e<sub>int</sub> + P<sub>gas</sub> + E + P<sub>rr</sub>)/ρ (dashed orange: gas part only), kinetic-energy flux, and their sum (white; 1 in a thermally relaxed state). The enthalpy and kinetic fluxes are the parts carried by the fluctuations about the shell-mean mass flux; what the mean radial flow (breathing / pulsation) carries, 4πr²⟨ρv<sub>r</sub>⟩⟨h + v²/2⟩ including the radiation enthalpy, is shown separately (grey dashed), and the dotted white curve adds it back (the lab-frame total). Dotted: the old ideal-gas enthalpy flux h = (5/3)e<sub>int</sub>/ρ. Cell-centred F<sub>r</sub> is noisy (odd–even) in the very optically thick deep shells; the history’s face-based L<sub>bot</sub> is the reliable number there.')}
</section>

{movies}

{SEC3D}

<section>
<h2>Slices</h2>
<p>Static slices of the latest dump on native cells (PNG, 120 dpi). Rows 1–2: equatorial (z = 0, panels 4,0,1,2 at η = 0) and meridional (y = 0, panels 5,0,3,2 at ξ = 0) rings through the whole star and the ejecta out to 3.0 R<sub>ph,IC</sub> (sponge beyond 2.77, domain edge 3.07), colour ranges from 0.25–1 R<sub>ph</sub>. Row 3: zoom r &lt; 0.25 R<sub>ph</sub> with its own ranges (dashed: the inner wall). Columns: log ρ; log ρ − ⟨log ρ⟩<sub>shell</sub>; T<sub>gas</sub> (grey = {TFL:.0f} K floor); v<sub>r</sub>. Below: the photosphere band unwrapped (log ρ, v<sub>r</sub>, T<sub>gas</sub>, δT<sub>gas</sub>/⟨T<sub>gas</sub>⟩).</p>
{slices}
{phot}
</section>

<section>
<h2>Maps</h2>
<p>Whole-sphere Mollweide maps of the latest dump, each 0.70° cubed-sphere cell drawn (every screen pixel takes its own cell): v<sub>r</sub> and δT<sub>gas</sub>/⟨T<sub>gas</sub>⟩ on the shells 0.1, 0.3, 0.6, 0.9 and 0.98 R<sub>ph,IC</sub>; a shell that moves as a whole is shown as v<sub>r</sub> − ⟨v<sub>r</sub>⟩, with the mean in the label. Also R<sub>ph</sub>/R<sub>ph,IC</sub> and T<sub>eff</sub> per column (grey Eddington, from T<sub>rad</sub> and the flux).</p>
{maps}
</section>

{NOTES}


<footer>Run: Caltech {html.escape(os.environ.get('RUN', ''))} (READ ONLY). Page rebuilt by /resnick/groups/carnegie_poc/jingze/hgglob_page/update.sh (reduce.py → mk_figs.py → mk_movies.py → build_page.py; copy of the RSG E page pipeline rsgE_page, base: Raven RSGPAGE 10-19). Movies are separate files under movies/.</footer>
</main>
</body></html>
'''
open(OUTF, 'w').write(page)
print(OUTF, os.path.getsize(OUTF))
