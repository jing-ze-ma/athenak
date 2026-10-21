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
jobtxt = '; '.join(f'{j[0]} {j[1].lower()} {j[2]}' for j in jobs) or 'squeue: no rsgP25 jobs'
health_bad = T['fatal'] + T['div'] + T['nc'] + T['nonfinite'] + T['nan']
health = (f'<span class="{"bad" if health_bad else "ok"}">FATAL {T["fatal"]}, DIVERGED {T["div"]}, summary NC {T["nc"]}, '
          f'NaN in logs {T["nan"]}, non-finite cells in dumps {T["nonfinite"]}</span>')

slices = ''.join(fig(fn, 'slices t %.1f d' % t, f'Slices at t = {t:.1f} d.') for fn, t in T['slices'])
maps = ''.join(fig(fn, 'maps t %.1f d' % t, f'Mollweide maps at t = {t:.1f} d.') for fn, t in T['maps'])

kv = [
    ('status', f'{html.escape(SITE)}: {len(running)} link running, {len(queued)} queued — {html.escape(jobtxt)}'),
    ('now', f't = {T["t_s"]:.4g} s = <b>{T["t_d"]:.2f} d</b>, cycle {T["cyc"]}, dt {T["dt"]:.1f} s, {T["spc"]:.2f} s/cycle '
            f'(≈ {86400 / max(T["dt"], 1) * T["spc"] / 3600:.2f} h wall per sim day); Picard mean {T["pic"]:.2f}, max {T["picmax"]:.0f}'),
    ('health', health),
    ('last dump', f'{T["t_dump"]:.1f} d ({T["ndump"]} dumps reduced): median R<sub>ph</sub> {T["rph"]:.4f} R<sub>ph,IC</sub>, '
                  f'T<sub>eff</sub> {T["teff"]:.0f} K, max |v<sub>r</sub>| {T["vmax"]:.1f} km/s'),
    ('budget', f'L<sub>top</sub>/L<sub>in</sub> {T["LL"]:.4f}, L<sub>bot</sub>/L<sub>in</sub> {T["Lbot"]:.4f}; '
               f'ΔE<sub>tot</sub>/|E<sub>0</sub>| {T["dE"]:.2e} vs ∫(L<sub>in</sub>−L<sub>top</sub>)dt/|E<sub>0</sub>| {T["dEb"]:.2e}; '
               f'ΔM/M {T["dM"]:.2e}; KE<sub>int</sub> {T["ke"]:.3e} erg'),
    ('star', 'RSG “golden25” MESA model (profile438): R<sub>ph</sub> 5.805×10<sup>13</sup> cm = 834 R<sub>⊙</sub> (IC), '
             'L 6.85×10<sup>38</sup> erg/s = 1.79×10<sup>5</sup> L<sub>⊙</sub>, T<sub>eff</sub> 4.11 kK; X 0.654, Z 0.0198'),
    ('grid', 'Global cubed sphere (gnomonic equiangular), N = 128 per face (0.70°), 6 panels; 591 stretched radial cells '
             'r = 1.761×10<sup>12</sup> – 1.761×10<sup>14</sup> cm (0.03 – 3 R<sub>ph</sub>; dr 3×10<sup>10</sup> cm at the photosphere); 58 M cells'),
    ('layout', html.escape(LAYOUT)),
    ('physics', 'Implicit grey M1 radiation with the VET closure (vet_gd), general EOS, MESA kap Rosseland table '
                '(gs98 + lowT_fa05), adaptive MLT closure (no ramp), α<sub>MLT</sub> 3; floors re-derived for the RSG '
                '(ρ 10<sup>−20</sup>, T 2000 K, sponge on ρ &lt; 2×10<sup>−17</sup>)'),
    ('seed', 'Real spherical harmonics Y<sub>lm</sub>, l = 1–8 (all m, 80 terms) in the global angles, amplitude 10<sup>−3</sup>, '
             'random signs, radial band 2.94×10<sup>12</sup> – 5.74×10<sup>13</sup> cm (he_seed_ylm)'),
    ('binary / deck', 'ylmseed-1019 @ a1c4c1e3 (prodmerge-1019 74811de1 + Y<sub>lm</sub> seed): '
                      '<code>athena_he_a100_a1c4c1e3</code> md5 56387536; deck <code>rsg_g25_prod_ylm.athinput</code> md5 b4795511; '
                      'IC <code>ic/pass3b/ic_rsg_g25.txt</code> md5 87918742; outputs hydro_w + m1 every 1 d, m1_vet 5 d, rst 0.5 d'),
]
kvh = ''.join(f'<div><b>{k}</b>{v}</div>' for k, v in kv)

page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>RSG g25 Production</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Source+Sans+3:wght@400;600&family=Source+Code+Pro:wght@400;500&display=swap">
<style>{CSS}</style></head><body>
<main>
<header>
<h1>Red supergiant g25: global 3-D production run</h1>
<p class="lede">AthenaK radiation hydrodynamics of the full convective envelope of a red supergiant on a global cubed-sphere mesh,
from 0.03 R<sub>ph</sub> (inside the radiative layer below the convection zone) to 3 R<sub>ph</sub>. The run started from t = 0
with the spherical-harmonic seed on 10-19 and goes on Raven in 24 h links. This page reads the run’s history, logs and 1 d dumps.</p>
<div class="kv">{kvh}</div>
<p class="status">{html.escape(STATUS)}</p>
</header>

<section>
<h2>History</h2>
{fig('fig_hist.png', 'history panels', 'From <code>rsg.user.hst</code> (every ~cycle; rows rewritten after a restart are dropped). L<sub>in</sub> is the luminosity imposed at the inner boundary; L<sub>bot</sub>, L<sub>mid</sub> and L<sub>top</sub> are the code’s radiative luminosities at the bottom, middle and top of the star region; L<sub>MLT</sub> is the sub-grid convective luminosity of the adaptive closure. The energy panel compares the change of the total energy E<sub>tot</sub> with the time integral of L<sub>in</sub> − L<sub>top</sub>; this ignores the flux between the top of the star region and the outer boundary and the sponge / source terms, so it is a consistency check, not a strict conservation test. The dashed grey KE curve is the archived first ~3 h link with the old (4-fold symmetric) seed. R<sub>ph</sub> and T<sub>eff</sub> come from the dumps: per column, R<sub>ph</sub> is where σT<sub>rad</sub><sup>4</sup>r<sup>2</sup> equals the column’s outgoing r<sup>2</sup>F<sub>r</sub> at 1.2–1.4 R<sub>ph</sub> (grey Eddington).')}
</section>

<section>
<h2>Radial profiles</h2>
{fig('fig_profiles.png', 'radial profiles', 'Solid-angle-weighted shell means over the full sphere at a few dump times (t = 0 is the first dump, after the in-code IC relax). Dotted white: the IC column (ρ, T<sub>gas</sub>, scaffold F<sub>MLT</sub>/F). Velocities are ρ-weighted rms; the tangential speed uses the cubed-sphere’s non-orthogonal unit basis (e<sub>1</sub>·e<sub>2</sub> = −sin ξ sin η). The resolved convective luminosity is 4πr²(⟨ρv<sub>r</sub>h⟩ − ⟨ρv<sub>r</sub>⟩⟨h⟩) with an ideal-gas h = (5/3)e<sub>int</sub>/ρ (no ionisation terms), so it underestimates the true enthalpy flux in the H/He ionisation zones.')}
</section>

<section>
<h2>Slices</h2>
<p>Rows 1–2: equatorial (z = 0) and meridional (y = 0) planes through the full sphere, out to 1.5 R<sub>ph,IC</sub>, with colour ranges set over 0.2–1 R<sub>ph</sub> (the deep region saturates). Row 3: zoom on r &lt; 0.2 R<sub>ph</sub> of the equatorial plane with its own ranges; dashed circle = MESA CZ base. All sampled with the nearest cubed-sphere cell (720 angles). Columns: log ρ; log ρ relative to its shell mean (the convective pattern); T<sub>rad</sub> = (E/a)<sup>1/4</sup>; v<sub>r</sub> (symmetric range from the 99.5 % percentile of |v<sub>r</sub>|).</p>
{slices}
</section>

<section>
<h2>Maps</h2>
<p>Whole-sphere Mollweide maps: v<sub>r</sub> on four shells (0.07, 0.3, 0.9, 0.98 R<sub>ph,IC</sub>; at 0.07 R<sub>ph</sub> look for the cube edges), the photosphere radius relative to the IC, and T<sub>eff</sub> per column.</p>
{maps}
</section>

<section>
<h2>What the run shows so far (read at t ≈ 31 d, 10-19)</h2>
<p class="note">This box is written by hand at the page build of 10-19 and is not regenerated; the figures and the header numbers are live.</p>
<ul>
<li><b>Healthy:</b> no FATAL / DIVERGED / non-converged summaries, no NaN, no non-finite cells in any dump. Energy bookkeeping closes: ΔE<sub>tot</sub>/|E<sub>0</sub>| = −3.67×10<sup>−7</sup> against ∫(L<sub>in</sub> − L<sub>top</sub>)dt/|E<sub>0</sub>| = −3.69×10<sup>−7</sup>; mass drift 3×10<sup>−12</sup>.</li>
<li><b>Surface steady:</b> L<sub>top</sub>/L<sub>in</sub> = 1.001; median R<sub>ph</sub> within 0.05 % of the IC, T<sub>eff</sub> 4110 K. No envelope collapse or runaway outer flux; the sponge power (≤ 10<sup>31</sup> erg/s, with spikes that recur more often after ~20 d) is 10<sup>−8</sup> of L.</li>
<li class="bad-li"><b>ANOMALY: flows along the cubed-sphere panel edges near the inner boundary.</b> KE<sub>int</sub> grew from 10<sup>41</sup> to 1.3×10<sup>45</sup> erg (growth slowing after ~15 d), almost all of it at r ≲ 0.12 R<sub>ph</sub>: ρ-weighted rms v<sub>r</sub> 6–10 km/s at 0.04–0.08 R<sub>ph</sub>, single cells up to ~65 km/s, and ~22 km/s rms tangential flow right at the inner boundary (0.03 R<sub>ph</sub>). In both slice planes the largest |v<sub>r</sub>| at 0.04–0.12 R<sub>ph</sub> sits at angles 45°, 135°, 225°, 315° (±3°) at 10, 20 and 31 d, i.e. exactly where the plane crosses cube edges (zoom row of the slices); the Mollweide v<sub>r</sub> map at 0.07 R<sub>ph</sub> shows ±40 km/s streaks tracing all 12 cube edges (the four vertical edges at longitude ±45°, ±135° and the eight polar-cap edges at latitude ±35°), with almost nothing in the panel interiors. The layer below the MESA CZ base (0.049 R<sub>ph</sub>) was meant to be a quiet radiative buffer, so this looks like a numerical cs edge/vertex instability seeded at the small inner radii, not convection. The time step fell from 230 s to ~148 s at 10–12 d as these flows grew.</li>
<li><b>Deep region out of thermal balance:</b> the lab-frame radiative luminosity 4πr²⟨F<sub>r</sub>⟩ reaches 2–3.6 L<sub>in</sub> at 0.04–0.15 R<sub>ph</sub>, plus a resolved enthalpy flux up to 0.2 L<sub>in</sub>. This is probably radiation enthalpy carried by the convective flows (the lab flux includes 4/3 E v), and the excess is drawn from the deep layer itself: T<sub>rad</sub> at 0.04–0.07 R<sub>ph</sub> has dropped by ~0.16 % and ρ at 0.07 R<sub>ph</sub> by 0.9 %, while 0.1–0.2 R<sub>ph</sub> warms by 0.07 %. L<sub>bot</sub>/L<sub>in</sub> drifted from 1.00 to 0.95–0.97 after 15 d.</li>
<li><b>Bulk envelope seed barely growing:</b> at 0.3, 0.6, 0.9 and 0.98 R<sub>ph</sub> the v<sub>r</sub> maps still show the smooth l ≤ 8 seed pattern at 0.3, 0.1, 0.04 and 0.06 km/s, against MLT speeds of ~5–10 km/s. The adaptive closure still carries ~100 % of the flux sub-grid (L<sub>MLT</sub>/L<sub>in</sub> = 1.00). With a turnover time of H<sub>p</sub>/v<sub>MLT</sub> ~ 100 d in the mid-envelope this may simply be early, but it should be watched over the next links.</li>
<li>The archived old-seed link (dashed grey KE) followed the same KE curve to ~17 d, so the deep growth does not come from the seed (consistent with a mesh-edge origin).</li>
</ul>
</section>

<section>
<h2>Caveats</h2>
<ul>
<li>T<sub>rad</sub> is the radiation temperature. It equals T<sub>gas</sub> in the optically thick envelope but not in the thin atmosphere; T<sub>gas</sub> needs the general-EOS table.</li>
<li>R<sub>ph</sub> and T<sub>eff</sub> use the grey Eddington relation per column, not a τ integral.</li>
<li>Slices and maps are nearest-cell samples of the cubed-sphere grid, no interpolation.</li>
<li>The energy-budget line is approximate (see the History caption).</li>
</ul>
</section>
<footer>Run: Raven /raven/ptmp/jinma/rsg_g25_1019/prod/R. Notes: rsg_g25_1019/NOTE.md. Page rebuilt by rsg_g25_1019/page/upd_page.sh (reduce.py → mk_figs.py → build_page.py).</footer>
</main>
</body></html>
'''
open(OUTF, 'w').write(page)
print(OUTF, os.path.getsize(OUTF))
