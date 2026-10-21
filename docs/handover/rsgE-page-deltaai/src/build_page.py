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
    ('now', f't = {T["t_s"]:.4g} s = <b>{T["t_d"]:.2f} d</b>, cycle {T["cyc"]}, dt {T["dt"]:.1f} s, {T["spc"]:.2f} s/cycle '
            f'(≈ {86400 / max(T["dt"], 1) * T["spc"] / 3600:.2f} h wall per sim day); Picard mean {T["pic"]:.2f}, max {T["picmax"]:.0f}'),
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
                'adaptive MLT closure (no ramp), α<sub>MLT</sub> 3; angular-momentum transport on; shared seam flux; '
                'floors ρ 10<sup>−20</sup>, T 2000 K, sponge on ρ &lt; 2×10<sup>−17</sup> at r &gt; 1.58×10<sup>14</sup> cm; '
                'dt_retry (cut 0.9, regrow 1.1 after 20 good steps)'),
    ('seed', 'Real spherical harmonics Y<sub>lm</sub>, l = 1–8 (all m), amplitude 10<sup>−3</sup>, random signs, '
             'radial band 2.95×10<sup>12</sup> – 5.74×10<sup>13</sup> cm'),
    ('binary / deck', 'fin51-r-fix @ 4de8168d: <code>athena_hesnf_gpu_4de8168d2b23</code> md5 c56693cb; deck '
                      '<code>rsg_E_clip.athinput</code> md5 9bd90b98 (Raven RSG E template 4c442452 + DeltaAI paths); '
                      'started fresh 2026-10-20 18:56 CDT, dt 10<sup>3</sup> s; outputs hydro_w + m1 every 1 d, m1_vet 5 d, rst 0.5 d in link 1 then 5 d; dumps move from /work/hdd to /projects after link 1. '
                      'Head-start copy of Raven’s RSG E production'),
]
kvh = ''.join(f'<div><b>{k}</b>{v}</div>' for k, v in kv)

TG_NOTE = ('T<sub>gas</sub> is inverted per cell from (ρ, e<sub>int</sub>) with the run’s own general EOS (H/He Saha, '
           'X 0.6538, Y 0.3264, no H<sub>2</sub>, no radiation term; numpy port of eos_composition.hpp on the code’s 1001 × 876 '
           'table grid); grey = cells at the 2000 K temperature floor.')


def video(fn, cap):
    return (f'<figure><video controls playsinline muted loop preload="metadata" src="movies/{fn}"></video>'
            f'<figcaption>{cap} (<a href="movies/{fn}">mp4</a>)</figcaption></figure>')


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
{video('rsgE_map.mp4', f"Mollweide maps, every 0.70° cell drawn, shells 0.90 and 0.98 R<sub>ph,IC</sub>: v<sub>r</sub> − ⟨v<sub>r</sub>⟩<sub>shell</sub> (±{L['map']['mvr_090']:.2g} / ±{L['map']['mvr_098']:.2g} km/s) and δT<sub>gas</sub>/⟨T<sub>gas</sub>⟩ (±{L['map']['mdT_090']:.2g} / ±{L['map']['mdT_098']:.2g} %). Longitude 0 is +x (centre of panel 0); graticule every 45°, so the vertical cube edges are the meridians at ±45°, ±135°.")}
</section>
"""
else:
    movies = ''
NOTES = open(os.path.join(HERE, 'notes.html')).read() if os.path.exists(os.path.join(HERE, 'notes.html')) else ''

page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>RSG E DeltaAI</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Source+Sans+3:wght@400;600&family=Source+Code+Pro:wght@400;500&display=swap">
<style>{CSS}</style></head><body>
<main>
<header>
<h1>Red supergiant g25, run E on DeltaAI</h1>
<p class="lede">AthenaK radiation hydrodynamics of the full convective envelope of a red supergiant on a global cubed-sphere mesh,
from 0.03 R<sub>ph</sub> (inside the radiative layer below the convection zone) to 3 R<sub>ph</sub>. This is the DeltaAI copy of RSG E,
started from t = 0 on 10-20 and run in 2 h interactive plus 48 h links. This page reads the run’s history, logs and 1 d dumps.</p>
<div class="kv">{kvh}</div>
<p class="status">{html.escape(STATUS)}</p>
</header>

<section>
<h2>History</h2>
{fig('fig_hist.png', 'history panels', 'From <code>rsg.user.hst</code> (every ~cycle; rows rewritten after a restart are dropped). L<sub>in</sub> is the luminosity imposed at the inner boundary; L<sub>bot</sub>, L<sub>mid</sub> and L<sub>top</sub> are the code’s radiative luminosities at the bottom, middle and top of the star region; L<sub>MLT</sub> is the sub-grid convective luminosity of the adaptive closure. The energy panel compares the change of the total energy E<sub>tot</sub> with the time integral of L<sub>in</sub> − L<sub>top</sub>; this ignores the flux between the top of the star region and the outer boundary and the sponge / source terms, so it is a consistency check, not a strict conservation test. R<sub>ph</sub> and T<sub>eff</sub> come from the dumps: per column, R<sub>ph</sub> is where σT<sub>rad</sub><sup>4</sup>r<sup>2</sup> equals the column’s outgoing r<sup>2</sup>F<sub>r</sub> at 1.2–1.4 R<sub>ph</sub> (grey Eddington).')}
</section>

<section>
<h2>Radial profiles</h2>
{fig('fig_profiles.png', 'radial profiles', 'Solid-angle-weighted shell means over the full sphere at a few dump times (t = 0 is the first dump, after the in-code IC relax). Temperature: T<sub>gas</sub> from the run’s EOS (solid) and T<sub>rad</sub> = (E/a)<sup>1/4</sup> (dashed). Dotted white: the IC column (ρ, T<sub>gas</sub>, scaffold F<sub>MLT</sub>/F). Velocities are ρ-weighted rms; the tangential speed uses the cubed-sphere’s non-orthogonal unit basis (e<sub>1</sub>·e<sub>2</sub> = −sin ξ sin η). The resolved convective luminosity is 4πr²(⟨ρv<sub>r</sub>h⟩ − ⟨ρv<sub>r</sub>⟩⟨h⟩) with an ideal-gas h = (5/3)e<sub>int</sub>/ρ (no ionisation terms), so it underestimates the true enthalpy flux in the H/He ionisation zones.')}
</section>

{movies}

<section>
<h2>Slices</h2>
<p>Static slices of the latest dump on native cells (PNG, 120 dpi). Rows 1–2: equatorial (z = 0, panels 4,0,1,2 at η = 0) and meridional (y = 0, panels 5,0,3,2 at ξ = 0) rings through the whole star out to 1.5 R<sub>ph,IC</sub>, colour ranges from 0.25–1 R<sub>ph</sub>. Row 3: zoom r &lt; 0.25 R<sub>ph</sub> with its own ranges (dashed: MESA CZ base). Columns: log ρ; log ρ − ⟨log ρ⟩<sub>shell</sub>; T<sub>gas</sub> (grey = 2000 K floor); v<sub>r</sub>. Below: the photosphere band unwrapped (log ρ, v<sub>r</sub>, T<sub>gas</sub>, δT<sub>gas</sub>/⟨T<sub>gas</sub>⟩).</p>
{slices}
{phot}
</section>

<section>
<h2>Maps</h2>
<p>Whole-sphere Mollweide maps of the latest dump, each 0.70° cubed-sphere cell drawn (every screen pixel takes its own cell): v<sub>r</sub> and δT<sub>gas</sub>/⟨T<sub>gas</sub>⟩ on the shells 0.11 (fastest deep flows), 0.3, 0.9, 0.98 and 1.085 R<sub>ph,IC</sub> (the infalling shell); a shell that moves as a whole is shown as v<sub>r</sub> − ⟨v<sub>r</sub>⟩, with the mean in the label. Also R<sub>ph</sub>/R<sub>ph,IC</sub> and T<sub>eff</sub> per column (grey Eddington, from T<sub>rad</sub> and the flux).</p>
{maps}
</section>

{NOTES}


<footer>Run: DeltaAI {html.escape(os.environ.get('RUN', ''))} (READ ONLY). Page rebuilt by /work/nvme/bivj/jma20/rsgE_page/update.sh (reduce.py → mk_figs.py → mk_movies.py → build_page.py; base: Raven RSGPAGE 10-19, athenak c542d6b0 docs/handover/rsg-page). Movies are separate files under movies/.</footer>
</main>
</body></html>
'''
open(OUTF, 'w').write(page)
print(OUTF, os.path.getsize(OUTF))
