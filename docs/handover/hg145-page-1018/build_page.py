#!/usr/bin/env python3
"""HGPAGE 10-18: assemble publish/index.html (images inlined base64; the movie stays a separate file movie3d/hg145.mp4).
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
MOV = os.path.join(os.path.dirname(OUTF), 'movie3d')
nfr = 0
MC = json.load(open(os.path.join(MOV, 'hg145.json')))
try:
    MC.update(json.load(open(os.path.join(MOV, 'hg145_scale.json'))))   # the resolved fixed scale of the frames
except Exception:
    pass
SITE = os.environ.get('SITE', 'Raven')
try:
    nfr = json.load(open(os.path.join(MOV, 'hg145_frames.json')))['n']
except Exception:
    pass


def img(fn, alt):
    p = os.path.join(IMG, fn)
    if not os.path.exists(p):
        return f'<p class="status">[missing figure {fn}]</p>'
    mt = 'image/jpeg' if fn.endswith('.jpg') else 'image/png'
    b = base64.b64encode(open(p, 'rb').read()).decode()
    return f'<img alt="{html.escape(alt)}" loading="lazy" src="data:{mt};base64,{b}">'


def fig(fn, alt, cap):
    return f'<figure>{img(fn, alt)}<figcaption>{cap}</figcaption></figure>'


def n(x, f='%.3f'):
    try:
        return f % x
    except Exception:
        return str(x)


jobs = '; '.join(f'{j[0]} {j[1].lower()} {j[2]}' for j in T['jobs'] if len(j) >= 3) or 'sacct unavailable'
LIVE = (f'{SITE} (this page, live)', os.environ.get('LAYOUT', '1 node × 4 GPUs, 4 MeshBlocks 704×64×64'),
        f'{T["t_d"]:.2f} d', f'{T["cyc"]}', f'{T["pic"]:.2f} / {T["picmax"]:.0f}', f'{T["spc"]:.2f}',
        f'{T["LL"]:.3f}', f'FATAL {T["fatal"]}, DIVERGED {T["div"]}, summary NC {T["nc"]}; jobs: {jobs}')
STATIC = {
    'Raven': ('Raven copy', '1 node × 4 A100, chain of 3 × 24 h links', '9.1 d (bundle time)', '–', '3.5 at 9 d', '0.73', '0.998',
              'running (job 31263625, link 1 of 3); separate Raven page build'),
    'DeltaAI': ('DeltaAI copy', '1 node × 4 GH200 (ghx4), deck d70305a7 (@FILES@ only)', '28.0 d', '24150',
                '3 → 5 by ~20 d', '0.49 → 0.68', '0.95 at ~27 d',
                'link 1 (3403162) ended at 28.0 d on an rst write (nvme quota full; I/O, not physics); restarted as 3405272 from 27.5 d'),
}
rows = [LIVE] + [v for k, v in STATIC.items() if k != SITE] + [
    ('Caltech copy', '4 H200 (TASK drafted)', '–', '–', '–', '–', '–', 'TASK drafted, HELD for the USER go (not pushed)')]
rowh = ''.join('<tr>' + ''.join(f'<td{" class=n" if 1 < i < 7 else ""}>{html.escape(c) if i != 0 else "<b>" + html.escape(c) + "</b>"}</td>'
                                for i, c in enumerate(r)) + '</tr>' for r in rows)

slices = ''.join(fig(fn, 'slices t %.2f d' % t, f'Slices at t = {t:.2f} d.') for fn, t in T['slices'])

if os.path.exists(os.path.join(MOV, 'hg145.mp4')):
    video = ('<figure><video controls playsinline preload="metadata" poster="movie3d/hg145_poster.png">'
             '<source src="movie3d/hg145.mp4" type="video/mp4">Your browser does not play MP4 video.</video>'
             f'<figcaption>3-D view of the wedge, one frame per 0.25 d dump ({nfr} frames so far, 22 frames/s). '
             'The cap is the photosphere R<sub>ph</sub>(θ, φ), coloured by the local T<sub>eff</sub>; the four side faces of the '
             'wedge (φ = 0 and 60°, θ = 60° and 120°) are cut from the inner boundary (8 R<sub>⊙</sub>) up to the local photosphere '
             f'and coloured by the radial velocity (±{MC["vr"]:.3g} km/s; T<sub>eff</sub> {MC["teff"][0]:.4g}–{MC["teff"][1]:.4g} kK; both scales are fixed for the whole movie, from the 99.5 % |v<sub>r</sub>| and the 0.5–99.5 % T<sub>eff</sub> percentiles over all dumps so far). The title gives t and the median R<sub>ph</sub>. '
             'The movie is rebuilt as the run goes; only new dumps are rendered.</figcaption></figure>')
else:
    video = '<p class="status">[movie not built yet]</p>'

page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>He Giant 2.75M Production</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Source+Sans+3:wght@400;600&family=Source+Code+Pro:wght@400;500&display=swap">
<style>{CSS}</style></head><body>
<main>
<header>
<h1>He giant (2.75 M<sub>⊙</sub> He star from a 14.5 M<sub>⊙</sub> progenitor): 3-D production run</h1>
<p class="lede">AthenaK radiation hydrodynamics of the extended envelope of a stripped He star (MESA r15140 <code>14p5M_widebinary_profile14</code>, S. Wu) on a π/3 × π/3 spherical wedge. The envelope runs from 8 R<sub>⊙</sub> to well beyond the photosphere at ~114 R<sub>⊙</sub>. The setup is the screened one: MESA-consistent opacity, a MESA-anchored starting column and the adaptive MLT closure. Copies of the same run go on Raven and DeltaAI in walltime-limited links.</p>
<div class="kv">
<div><b>star</b>He star 2.754 M<sub>⊙</sub>, L 3.71×10<sup>4</sup> L<sub>⊙</sub>, R<sub>ph</sub> 114.4 R<sub>⊙</sub> (IC), T<sub>eff</sub> 7.5 kK; envelope X 9.5×10<sup>−4</sup>, Y 0.978, Z 0.0215 (uniform 1–110 R<sub>⊙</sub>)</div>
<div><b>initial condition</b>tools/ic_star spec <code>hegiant_14p5M_mk</code> (ic-mesascaf-1017): MESA T(P) kept below τ<sub>R</sub> 3, radiative equilibrium re-solved above it (Hopf), plus the convective scaffold F<sub>MLT</sub>/F = 1 − ∇<sub>MESA</sub>/∇<sub>rad</sub>. In-code he_ic_equil relax at startup.</div>
<div><b>grid</b>704 × 128 × 128 (r, θ, φ); θ = π/3 – 2π/3, φ = 0 – π/3 (solid angle π/3 sr = 1/12 of the sphere); stretched r = 8 – 355.6 R<sub>⊙</sub> (fit145b_704); 4 MeshBlocks 704×64×64</div>
<div><b>physics</b>Implicit grey M1 radiation with the VET closure (vet_gd twin), general EOS, Rosseland table <code>opac_mk_afe2dd8da2bd</code>: MESA r15140 default kap below log T 4.6–4.8, spliced into the TOPS/ext2 recipe. <code>problem/mlt_closure = adaptive</code>: a live shell-mean closure for the sub-grid convective flux.</div>
<div><b>boundaries</b>Bottom: imposed radiative flux (L<sub>in</sub>), closed hydro. Top: Marshak radiation plus a radial sponge above 320 R<sub>⊙</sub> (rate 10<sup>−3</sup>, ρ<sub>max</sub> 2×10<sup>−14</sup>). Floors: ρ 10<sup>−18</sup>, p 10<sup>−11</sup>, T 2400 K (baseline, kept after the fl15 test).</div>
<div><b>binary / commit</b>prodmerge-1018-hg @ 4c494982 (fin51-r-fix bb360d8f + mltp-1017 + ic-mesascaf-1017). Raven A100 <code>athena_he_a100_4c494982</code> md5 a6fa347d. MI300A strict-md5 gate vs bb360d8f: bsg/e1/hg bitwise.</div>
<div><b>deck</b><code>hg145_prod.athinput</code> md5 331f76a4 (Raven). IC file md5 f9678c2a; output cadence: hst every cycle; hydro_w, m1 and m1_vet every 0.25 d; rst every 0.5 d. tlim 10<sup>9</sup> (walltime-limited links).</div>
<div><b>now</b>t = {T["t_d"]:.2f} d ({T["cyc"]} cycles, dt {T["dt"]:.0f} s); last dump {T["t_dump"]:.2f} d: median R<sub>ph</sub> {T["rph"]:.2f} R<sub>⊙</sub>, T<sub>eff</sub> {T["teff"]:.0f} K; L<sub>top</sub>/L<sub>in</sub> {T["LL"]:.3f}</div>
</div>
<p class="status">{html.escape(STATUS)}</p>
</header>

<section>
<h2>3-D movie</h2>
{video}
</section>

<section>
<h2>Runs and job status</h2>
<p>All copies use the same commit and deck and start from t = 0 with the same IC. The live row is read from this copy's run directory at build time; the other rows are copied from handover NOTEs and are not live.</p>
<div class="tbl"><table><thead><tr><th>copy</th><th>layout</th><th>t reached</th><th>cycle</th><th>Picard mean / max</th><th>s / cycle</th><th>L<sub>top</sub>/L<sub>in</sub> (last)</th><th>status</th></tr></thead><tbody>
{rowh}
</tbody></table></div>
<p class="note">The DeltaAI copy is ahead of Raven. Its L<sub>top</sub>/L<sub>in</sub> slid from 1.00 to 0.95 by ~27 d while Picard rose from 3 to 5 per solve. Raven will reach that time in its first link, so the history panel below will show whether the drift reproduces.</p>
</section>

<section>
<h2>History</h2>
{fig('fig_hist.png', 'history panels', 'From <code>hegiant.user.hst</code> (every cycle; rows rewritten after a restart are dropped). L<sub>top</sub> is the luminosity leaving at the top of the star region and L<sub>in</sub> the imposed bottom luminosity, both for the wedge (1/12 of the star). L<sub>MLT</sub> is the sub-grid convective luminosity of the adaptive closure. R<sub>ph</sub> and T<sub>eff</sub> come from the dumps. Per column, R<sub>ph</sub> is where σT<sub>rad</sub><sup>4</sup> r<sup>2</sup> equals the column’s outgoing flux × r<sup>2</sup>, measured at 1.2–1.4 R<sub>ph</sub> (the grey Eddington T(τ=2/3) = T<sub>eff</sub>), and T<sub>eff</sub> follows from that flux at R<sub>ph</sub>.')}
</section>

<section>
<h2>Radial profiles</h2>
{fig('fig_profiles.png', 'radial profiles', 'Shell means over the wedge (solid-angle weighted) at a few dump times; t = 0 is the first dump, written after the in-code IC relax. The dashed white lines are the MESA-anchored IC column (ρ, T<sub>gas</sub> and the scaffold F<sub>MLT</sub>/F). Top row: full star. Bottom row: zoom on 0.85–1.12 R<sub>ph</sub>, with ρ and T<sub>rad</sub> shown relative to t = 0. The convective flux fraction is the <em>resolved</em> enthalpy flux, ⟨ρ v<sub>r</sub> h⟩ − ⟨ρ v<sub>r</sub>⟩⟨h⟩ with h = (5/3) e<sub>int</sub>/ρ (ideal gas, no ionisation terms), over L<sub>in</sub>/(Ω r<sup>2</sup>). The sub-grid flux of the MLT closure is not in the dumps; its total is L<sub>MLT</sub> in the history.')}
</section>

<section>
<h2>Photosphere</h2>
{fig('fig_rph.png', 'photosphere maps', 'R<sub>ph</sub>(θ, φ) relative to its median, and T<sub>eff</sub>, at about half the current run time (top) and at the last dump (bottom), on a latitude–longitude grid covering the wedge. Definitions as in the history caption.')}
</section>

<section>
<h2>Slices</h2>
<p>Cartesian slices with lengths in R<sub>⊙</sub>. Rows 1–2: equatorial plane (mean of the two θ rows next to π/2), x = r cos φ, y = r sin φ, with a 36 R<sub>⊙</sub> zoom box on the photosphere at φ = 30°. Rows 3–4: meridional plane at φ = 30°, x = r sin θ, z = r cos θ; the straight edges are the θ = 60° and 120° walls. Slices run out to 1.6 R<sub>ph</sub>. Columns: log ρ, log ρ/ρ(t = 0), radiation temperature T<sub>rad</sub> = (E/a)<sup>1/4</sup> (log scale), and v<sub>r</sub> (symmetric range from the 99.5 % percentile of |v<sub>r</sub>|).</p>
{slices}
</section>

<section>
<h2>Pre-production checks (10-18)</h2>
<ul>
<li><b>Gate (MI300A, Viper):</b> strict-md5 gate of 4c494982 vs bb360d8f: bsg, e1 and hg are bitwise (job 12264965).</li>
<li><b>Viper background check</b> (production deck, 3 links on apudev, 1 node × 2 MI300A, 1.27 d): clean. The interior is quiet; the 0.9–0.99 R<sub>ph</sub> band has −0.37 km/s; in the 1.03–1.10 R<sub>ph</sub> layer above the photosphere ρ rose 32 %, and the rise is slowing.</li>
<li><b>DeltaAI smoke, Test A</b> (4 × GH200, production deck, 300 cycles): clean. Restart vs continuous is bitwise in every hst row and the final bins; 0.52 s/cycle; 22 GB per GPU; Picard 3.06; L<sub>top</sub>/L<sub>in</sub> 0.847 → 1.004 by 0.17 d.</li>
<li><b>DeltaAI fl15 screen (B) vs baseline twin (B0)</b>, 1 GH200, 3 d: both clean (NC 0, DIVERGED 0). KE grows in both runs, in the same bands, and B0 saturates first. Decision: keep the baseline floor and sponge (the production deck).</li>
</ul>
<div class="tbl"><table><thead><tr><th>t [d]</th><th>B (fl15) KE<sub>int</sub> [erg]</th><th>B0 (baseline) KE<sub>int</sub> [erg]</th></tr></thead><tbody>
<tr><td class=n>0.25</td><td class=n>2.80e37</td><td class=n>5.03e37</td></tr>
<tr><td class=n>1.00</td><td class=n>5.89e37</td><td class=n>1.04e38</td></tr>
<tr><td class=n>2.00</td><td class=n>8.83e37</td><td class=n>3.33e38</td></tr>
<tr><td class=n>2.50</td><td class=n>2.22e38</td><td class=n>7.93e38</td></tr>
<tr><td class=n>3.00</td><td class=n>6.25e38</td><td class=n>8.42e38 (saturated)</td></tr>
</tbody></table></div>
<div class="tbl"><table><thead><tr><th>band r/R<sub>ph</sub></th><th>B: KE 2 → 3 d (share of ΔKE), rms |v| [km/s]</th><th>B0: KE 2 → 3 d (share), rms |v| [km/s]</th></tr></thead><tbody>
<tr><td class=n>0.00–0.90</td><td>~1e33, 0 %</td><td>~1e33, 0 %</td></tr>
<tr><td class=n>0.90–0.99</td><td>6.3e37 → 3.9e38 (60 %), 0.19 → 0.48</td><td>1.3e38 → 4.8e38 (64 %), 0.27 → 0.53</td></tr>
<tr><td class=n>0.99–1.01</td><td>4.5e37 → 2.6e38 (39 %), 0.29 → 0.70</td><td>2.4e38 → 4.3e38 (35 %), 0.66 → 0.90</td></tr>
<tr><td class=n>1.01–1.03</td><td>5.0e35 → 4.1e36 (1 %), 0.39 → 0.96</td><td>3.2e35 → 6.0e36 (1 %), 0.32 → 1.13</td></tr>
<tr><td class=n>1.03–1.10</td><td>5.7e33 flat, 1.36 → 1.70; max |Δρ/ρ| 1.58</td><td>1.2e32 → 2.9e32, 0.18 → 0.27; max |Δρ/ρ| 0.29</td></tr>
</tbody></table></div>
<p>So the KE growth comes from the near-surface convective layer at 0.90–1.01 R<sub>ph</sub> spinning up, not from the floor. The fl15 floor only delays its onset by ~0.5 d. It also makes the thin 1.03–1.10 R<sub>ph</sub> layer worse.</p>
</section>

<section>
<h2>Caveats</h2>
<ul>
<li>T<sub>rad</sub> is the radiation temperature from E = aT<sup>4</sup>. It equals T<sub>gas</sub> in the coupled envelope but not in the thin atmosphere. The gas temperature needs the general-EOS table.</li>
<li>R<sub>ph</sub> and T<sub>eff</sub> use the grey Eddington relation per column, not a τ integral (the dumps carry no opacity).</li>
<li>The convective flux fraction is only the resolved part, with an ideal-gas enthalpy. The adaptive closure carries the rest as sub-grid flux.</li>
<li>The DeltaAI and Caltech rows are copied from NOTEs and are not live.</li>
</ul>
</section>
<footer>Run: Raven /raven/ptmp/jinma/hg145_prod_1018/prod/R. Notes: hg145_prod_1018/NOTE.md, hegiant145_1015/NOTE.md sect. 9–11, and docs/handover/NOTE-2026-10-18-deltaai-hg145-*.md (bsg-files-1009). Page rebuilt by hg145_prod_1018/page/upd_page.sh (reduce.py → render3d.py → mk_figs.py → build_page.py).</footer>
</main>
</body></html>
'''
open(OUTF, 'w').write(page)
print(OUTF, os.path.getsize(OUTF))
