#!/usr/bin/env python3
"""ACC3D page 10-20: assemble publish/index.html ("RY Per Accretor 3-D"). Figures are inlined (base64); the movies stay
separate files next to it: movie/acc3d_eq.mp4, movie/acc3d_3d.mp4 (+ _poster.png). usage: build_page.py <imgdir> <out.html> <status>
Text facts come from img/prod.json, img/trc.json, img/omega.json and the env (SITE, LAYOUT, BIN_LABEL, DECK_LABEL, RUN_LABEL)."""
import os
import sys
import json
import base64
import html

IMG, OUTF, STATUS = sys.argv[1], sys.argv[2], sys.argv[3]
HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, 'ref_img')
PUB = os.path.dirname(os.path.abspath(OUTF))
E = os.environ.get
SITE, LAYOUT = E('SITE', 'Raven'), E('LAYOUT', '2 nodes x 4 A100, 8 MeshBlocks 420x58x64')
BIN_LABEL = E('BIN_LABEL', 'prodmerge-1019b @ 13206541 (A100 nofma md5 29caa1ab)')
DECK_LABEL = E('DECK_LABEL', 'acc3d_prod_v11amwtr')
RUN_LABEL = E('RUN_LABEL', 'fresh tracer run (run12)')
CHAIN_LABEL = E('CHAIN_LABEL', 'acc3d_chain.sh links, fresh from t = 0')
PORB, DAYS, TS, MDOT = 0.473169, 8.0521, float(E('TS', '0.05')), float(E('STREAM_MDOT', '254.843'))
P = json.load(open(IMG + '/prod.json'))
TR = json.load(open(IMG + '/trc.json')) if os.path.exists(IMG + '/trc.json') else dict(present=False, why='mk_tracer.py not run')
esc = html.escape


def img(fn, alt, d=IMG):
    p = os.path.join(d, fn)
    if not os.path.exists(p):
        return f'<p class="status">[missing figure {esc(fn)}]</p>'
    mt = 'image/jpeg' if fn.endswith('.jpg') else 'image/png'
    return f'<img alt="{esc(alt)}" loading="lazy" src="data:{mt};base64,{base64.b64encode(open(p, "rb").read()).decode()}">'


def fig(fn, alt, cap, d=IMG):
    if not os.path.exists(os.path.join(d, fn)):
        return ''
    return f'<figure>{img(fn, alt, d)}<figcaption>{cap}</figcaption></figure>'


def movies():
    out = ''
    for mp4, cap in [('acc3d_3d', '3-D view (r &le; d<sub>L1</sub>): nested density isosurfaces coloured by the gas temperature '
                      '(EOS-table inversion of &rho;, e<sub>int</sub>); one frame per hydro dump (&Delta;t 0.005 code).'),
                     ('acc3d_eq', 'Orbital plane (top) and edge-on cut through the donor axis (bottom, &phi; 0 on +x, 180&deg; on &minus;x), '
                      'inside the L1 radius (dashed). Columns: log &rho;, log T<sub>gas</sub> (EOS table), log P<sub>rad</sub>/P<sub>gas</sub> '
                      '(P<sub>rad</sub> = E/3 from the m1 dump, P<sub>gas</sub> from the EOS table at T<sub>gas</sub>), v<sub>&phi;</sub> '
                      '(rotating frame). One frame per hydro dump that has an m1 dump at the same time (every 0.01 code).')]:
        if os.path.exists(f'{PUB}/movie/{mp4}.mp4'):
            n = ''
            try:
                n = ' %d frames.' % json.load(open(f'{PUB}/movie/{mp4}_frames.json'))['n']
            except Exception:
                pass
            out += (f'<figure><video controls loop muted playsinline preload="metadata" poster="movie/{mp4}_poster.png" '
                    f'src="movie/{mp4}.mp4"></video><figcaption>{cap}{n}</figcaption></figure>\n')
    return out or '<p class="status">no movie built yet</p>'


def anomalies():
    a = []
    a.append(f"hst NaN entries: {P['hst_nan']}; run.log NaN / FATAL|DIVERGED lines: {P['nan_lines']} / {P['fatal_lines']}; "
             f"largest NON-CONVERGED in a solve summary: {P['nc_max']}" + (' <b>(flag)</b>' if P['nc_max'] or P['fatal_lines'] or P['hst_nan'] else ' (clean)'))
    rt = ', '.join(f'{x:.3f}' for x in P['retry_t'][:12]) + (' …' if len(P['retry_t']) > 12 else '')
    a.append(f"dt retries (rollbacks, eint &le; 0 at the stream onset / contact): {P['retries']}" + (f", at t ≈ {rt}" if P['retries'] else ''))
    a.append(f"dt now {P['dt']:.2e}; smallest dt in the last 200 hst rows {P['dt_min_late']:.2e}"
             + (' <b>(below 1e-6)</b>' if P['dt_min_late'] < 1e-6 else ' (the impact region limits dt to ~3e-6 after contact, as in run11)'))
    j = P['jumps']
    a.append(f"largest single-row relative jumps: M<sub>dom</sub> {j['Mdom'][0]:.1e} (t {j['Mdom'][1]:.4f}), J<sub>dom</sub> {j['Jdom'][0]:.1e} "
             f"(t {j['Jdom'][1]:.4f}), E<sub>tot</sub> {j['Etot'][0]:.1e} (t {j['Etot'][1]:.4f}), L<sub>meas</sub> {j['Lmeas'][0]:.2f} (t {j['Lmeas'][1]:.4f})")
    if 'wall' in P:
        w = P['wall']
        a.append(f"θ-wall screen (t {w['t']:.3f}): mass at r 9.5–16.2 in the 3 rows next to each wall = {100 * w['frac']:.3f} % "
                 f"(uniform {100 * w['uniform']:.0f} %); " + ("<b>ABOVE 1 %</b>" if w['frac'] > 0.01 else "below the 1 % flag"))
    a.append(f"stream rst files kept: {P['nrst']}; hydro_w dumps {P['ndumps']}, m1 dumps {P['nm1']}"
             + ('; <b>chain.txt says STOPPED</b>' if P.get('stopped') else ''))
    return ''.join(f'<li>{x}</li>' for x in a)


def budgets():
    b = P['budget']
    rel = lambda x, y: abs(x) / max(abs(y), 1e-300)   # noqa: E731
    li = [f"<b>Mass</b>: ΔM<sub>dom</sub> {b['dM']:.4g} = M<sub>in</sub> {b['Min']:.4g} − M<sub>out</sub> ({b['Mout']:.3g}) − M<sub>msk</sub> "
          f"({b['Mmsk']:.2g}) + M<sub>fl</sub> ({b['Mfl']:.2g}) + residual {b['rM']:.2e} (max |r| {b['rM_max']:.1e}, "
          f"{rel(b['rM_max'], P['Mdom0']):.0e} of M<sub>dom</sub>). M<sub>out</sub> &lt; 0 is net inflow."]
    li.append(f"<b>Energy</b>: ΔE<sub>tot</sub> {b['dE']:.4g}, E<sub>bnd</sub> {b['Ebnd']:.4g}, E<sub>rel</sub> {b['Erel']:.3g}, "
              f"E<sub>msk</sub> {b['Emsk']:.3g}; unbooked r<sub>E</sub> = ΔE<sub>tot</sub> − E<sub>bnd</sub> − E<sub>rel</sub> + E<sub>msk</sub> = "
              f"{b['rE']:+.3e} = {100 * rel(b['rE'], b['Ebnd']):.2g} % of the boundary throughput ({rel(b['rE'], b['Etot0']):.1e} of |E<sub>tot</sub>|)"
              + (f"; stream phase only: {b['rE_stream']:+.3e} = {100 * rel(b['rE_stream'], b['Ebnd_stream']):.2g} % of its E<sub>bnd</sub>" if 'rE_stream' in b else '')
              + ". Reference: the accwall-1019 gate (0.0055 code from the run11 t 0.42 state) closed it to 1.0e-3 of |E<sub>bnd</sub>| (23 % before the booking fix); short early intervals (the static start, the nozzle filling) show larger fractions.")
    li.append(f"<b>Angular momentum</b> (inertial L<sub>z</sub> about the accretor axis, <code>am_transport = true</code>: exact advective L<sub>z</sub>): "
              f"ΔJ<sub>dom</sub> {b['dJ']:.4g}, nozzle J<sub>in</sub> {b['Jin']:.4g}, J<sub>out</sub> {b['Jout']:.3g}; "
              f"r<sub>J</sub> = ΔJ<sub>dom</sub> − (J<sub>in</sub> − J<sub>out</sub>) = {b['rJ']:+.4g}"
              + (f" ({b['rJ_rate']:.4g} per code time in the stream phase vs J<sub>in</sub> {b['Jin_rate']:.4g})" if 'rJ_rate' in b else '')
              + ". r<sub>J</sub> is the net external torque (donor tide and the non-inertial frame of the orbiting accretor), which the hst does not book "
              "(the accwall gate had r<sub>J</sub> ≈ 3.9× J<sub>in</sub> over 0.0055 code time).")
    return ''.join(f'<li>{x}</li>' for x in li)


def tracer_section():
    if not TR.get('present'):
        return f"""<section><h2>Tracers</h2><p class="status">Tracer section skipped: {esc(TR.get('why', ''))}.</p></section>"""
    T = TR
    rows = ''
    for d in T.get('dumps', []):
        a = d['acc']
        if not a['M'] > 0:
            rows += (f"<tr><td class='n'>{d['t']:.3f}</td><td class='n'>{(d['t'] - TS) * DAYS:.2f}</td><td class='n'>{a['Mtot']:.4g}</td>"
                     f"<td colspan='8' class='status'>no s0 gas at r &lt; 10.3 yet (stream not on the star)</td><td class='n'>{d['phi_impact']:.0f}</td></tr>")
            continue
        rows += (f"<tr><td class='n'>{d['t']:.3f}</td><td class='n'>{(d['t'] - TS) * DAYS:.2f}</td><td class='n'>{a['Mtot']:.4g}</td>"
                 f"<td class='n'>{a['M']:.4g}</td><td class='n'>{a['M_rmeas']:.4g}</td><td class='n'>{a['s']:.2f}</td>"
                 f"<td class='n'>{a['pct'][0]:.2f}–{a['pct'][2]:.2f}</td><td class='n'>{a['srad']:.2f}</td>"
                 f"<td class='n'>{a['rpct'][0]:.2f}–{a['rpct'][2]:.2f}</td>"
                 f"<td class='n'>{('%.3g' % a['tau']) if a.get('tau') else '–'}</td><td class='n'>{d['stream']['s']:.2f}</td>"
                 f"<td class='n'>{d['phi_impact']:.0f}</td></tr>")
    rt = (f"dNz0/dt over the last 0.05 = {T['Nz0_rate']:.2f} (nozzle nominal {MDOT:.2f}, ratio {T['Nz0_rate'] / MDOT:.4f}); "
          f"dMs0/dt = {T['Ms0_rate']:.2f}") if 'Nz0_rate' in T else 'too few rows for a rate'
    nb = ''
    if 'r_nb' in T:
        nb = (f"<li><b>Buoyancy:</b> the accreted gas's s<sub>gas</sub> ({T['dumps'][-1]['acc']['s']:.2f}) equals the pre-stream star's equatorial "
              f"shell entropy at r ≈ {T['r_nb']:.3f} R<sub>⊙</sub>" + (f" (τ ≈ {T['tau_nb']:.3g})" if 'tau_nb' in T else '') + '.</li>')
    sref = ''.join(f" τ {k[9:]}: r {v[0]:.3f}, s<sub>gas</sub> {v[1]:.2f};" for k, v in T.items() if k.startswith('star_s_tau'))
    why = f"<p class='status'>{esc(T['why_dumps'])}</p>" if T.get('why_dumps') else ''
    return f"""<section>
<h2>Tracers: where the stream gas goes</h2>
<p>Two passive scalars (<code>&lt;hydro&gt;/nscalars = 2</code>): <b>s0</b> = stream gas (the nozzle reset writes s0 = 1), <b>s1</b> = gas that
entered through the other user boundaries. <code>ryper.trc.hst</code> books the s0/s1 mass in the domain and the cumulative nozzle (Nz),
sink (Sk), inner-boundary (Fi), outer + θ-wall (Fo) and boundary-inflow (Bin) terms; the registers continue across restarts (1019c).</p>
<ul>
<li><b>Injection:</b> {rt}. Since t<sub>on</sub> = {T['t_on']:.4f}: Nz0 = {T['Nz0']:.4g} (the {T['fill']:.3g} that filled the nozzle at switch-on + {MDOT:.2f} × Δt = {T['nominal']:.4g}),
Ms0 = {T['Ms0']:.4g} in the domain, Fi0 = {T['Fi0']:.3g} through r<sub>in</sub>, Fo0 = {T['Fo0']:.3g} through r<sub>out</sub> + walls, Sk0 = {T['Sk0']:.3g} in the sinks.</li>
<li><b>Budgets:</b> max |ΔMs0 − (Nz0 − Sk0 − Fi0 − Fo0)| = {T['bud0']:.1e}, s1 {T['bud1']:.1e} (round-off of the c2p clip of negative ρs);
s0 range [{T['s0min']:.2g}, {T['s0max']:.4g}], max(s0 + s1) {T['smax']:.4g}.</li>
<li><b>Boundary-inflow tracer:</b> Ms1 = {T['Ms1']:.3g}, Bin = {T['Bin']:.3g} (closed θ walls: only r<sub>out</sub> inflow).</li>
{nb}
<li><b>Original star</b> (end of the static phase, equatorial shell means):{sref}</li>
</ul>
{why}
{fig('fig_trc_hst.png', 'tracer history', 'Tracer history. Left: s0 mass in the domain, the cumulative nozzle term and the nominal nozzle fill + %.2f × (t − t<sub>on</sub>), and the s0 leaving through r<sub>in</sub>, r<sub>out</sub>/walls and the sinks; dotted: s1. Middle: the rates. Right: the s0 and s1 budget residuals.' % MDOT)}
<div class="tbl"><table><thead><tr><th>t</th><th>stream [d]</th><th>s0 in domain</th><th>s0 at r &lt; 10.3</th><th>s0 at r &lt; r<sub>meas</sub></th>
<th>s<sub>gas</sub> (s0-wtd, r &lt; 10.3)</th><th>p10–p90</th><th>s<sub>rad</sub></th><th>r p10–p90</th><th>τ (geo-mean)</th><th>s<sub>gas</sub> stream r &gt; 10.3</th><th>impact φ [°]</th></tr></thead>
<tbody>{rows}</tbody></table></div>
<p class="status">s in k<sub>B</sub>/baryon from the run's own EOS table (thermodynamic integration, matches Saha H/He to ≤ 0.03); T<sub>gas</sub> from the
table inversion; s<sub>rad</sub> = 4aT<sup>3</sup>/3ρ; τ = radial Rosseland depth from r<sub>out</sub> (MESA-kap table of the deck). True stretched faces
(format-1.2 mb_x1f/mb_x2f) for every volume. "On the star" = r &lt; 10.3 R<sub>⊙</sub> (as in the run11 analysis).</p>
{fig('fig_trc_maps.png', 's0 maps', 'Stream-gas fraction s0 in the r–θ plane (R<sub>cyl</sub>, z) at the analysed times, r ≤ d<sub>L1</sub>: left, φ-averaged and mass-weighted (Σ ρ s0 dV / Σ ρ dV over φ); right, the meridional cut through the impact azimuth (largest s0 mass at r 9.3–10.5). Log scale 1e-4 to 1.')}
{fig('fig_trc_prof.png', 's0 radial profiles', 'Mass-weighted s0 profiles: s0 mass per unit radius, the s0 mass fraction of each spherical shell, and the cumulative fraction of the domain s0 mass inside r.')}
{fig('fig_trc_entropy.png', 's0 entropy', 'Left: s0-weighted gas entropy of the accreted gas (r &lt; 10.3) against the original surface layers (t = 0, r 8.8–9.5). Middle: the accreted gas (s0-weighted mean s<sub>gas</sub> at its median radius, bars p10–p90) on the pre-stream star’s equatorial entropy profile. Right: the radial optical depth where the accreted s0 gas sits.')}
{fig('fig_omega_eq.png', 'Omega(r)', 'Equatorial Ω(r) in the inertial frame, in units of Ω<sub>orb</sub> (P 3.81 d): mass-weighted over |lat| &lt; 5° and all φ, Ω = Σ dm R v<sub>φ</sub> / Σ dm R² + Ω<sub>orb</sub>, at the end of the static phase and at up to 7 stream-phase times, out to d<sub>L1</sub>; red dashed: the s0-weighted Ω of the stream gas at the latest time. Dashed: corotation; dotted: Keplerian (16 M<sub>⊙</sub>); dash-dot: the accretor surface; grey dotted: r<sub>meas</sub>. True stretched faces (mk_omega.py).')}
</section>"""


def previous():
    return f"""<section>
<h2>Previous run: Raven run11 (stopped at t 0.73)</h2>
<p>Raven's first 3-D production run (deck v11 → v11amw, 2 × 4 A100, no tracers) ran 0 → 0.731 (5.5 d of stream) and was stopped on 10-19 for this fresh run:
its hse θ walls leaked ~1.8× the stream into the domain (fixed by <code>theta_wall_closed</code>), and its energy budget missed half the injected light
(booking fixed). The figures below are run11's own (fixed images, not recomputed).</p>
<ul>
<li>Accreted-gas entropy (run11, top-ΔM estimate, no tracers): settled accreted gas s<sub>gas</sub> 38.5 (p10–90 38.0–39.1); stream at the nozzle 39.5, cooling to
38.6–38.8 by r 11–12.5, reheated to 40.4 before impact; a hot impact layer (T ~1e5 K) s<sub>gas</sub> 41.0–41.4; original star 37.6 at τ ~6000, 38.8 at τ 100,
40.4 at the photosphere. By gas entropy the accreted gas is neutrally buoyant at τ ≈ 110–160 (0.1–0.2 R<sub>⊙</sub> below the photosphere). The tracer section above replaces this estimate.</li>
<li>Ω(r): a spun-up shear layer that deepened slowly to r ≈ 8.46 at +5.6 d; outside r<sub>meas</sub> the infalling stream (Ω ∝ r<sup>−2</sup>, j ≈ 3150).</li>
</ul>
{fig('run11_prod_hist.png', 'run11 history', 'run11 hst history (Raven, to t 0.73; per-link counters chained).', REF)}
{fig('run11_omega_eq.png', 'run11 Omega', 'run11 equatorial Ω(r) (true faces).', REF)}
{fig('run11_ent_profiles.png', 'run11 entropy profiles', 'run11 entropy vs radius in the impact band (φ 0.15–0.45, |lat| &lt; 1.4°) at t = 0 and t 0.74.', REF)}
{fig('run11_ent_hist.png', 'run11 entropy histograms', 'run11 mass distribution in s<sub>gas</sub> (surface layers, top-ΔM accreted estimate, impact band, stream).', REF)}
</section>"""


maps = ''.join(fig(m['file'], f"t {m['t']:.3f}", f"t = {m['t']:.3f} ({m['t'] / PORB:.3f} P<sub>orb</sub>; hydro dump {m['nh']}, m1 dump {m['nm']})."
                   + (' Latest paired dump.' if m is P['maps'][-1] else '')) for m in P['maps'][::-1])
CSS = open(os.path.join(HERE, 'page.css')).read()
q = '<br>'.join(esc(x) for x in P['queue']) or 'not listed'
chain = '<br>'.join(esc(x) for x in P['chain'][-2:]) or 'no chain.txt'
phase = 'static phase (stream off)' if P['t'] < TS else f"stream on for {(P['t'] - TS) * DAYS:.2f} d ({(P['t'] - TS) / PORB:.3f} P<sub>orb</sub>)"
page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>RY Per Accretor 3-D</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Source+Sans+3:wght@400;600&family=Source+Code+Pro:wght@400;500&display=swap">
<style>{CSS}</style></head><body>
<main>
<header>
<h1>RY Per accretor in 3-D: the fresh tracer run</h1>
<p class="lede">AthenaK radiation hydrodynamics of the 16 M<sub>⊙</sub> mass gainer of RY Per on a 3-D spherical-polar wedge (±20.4° about the orbital
plane, full 2π in φ) with the full 3-D Roche potential, Coriolis and centrifugal forces, in the frame co-rotating with the binary (donor towards +x).
This page follows the {esc(RUN_LABEL)} on {esc(SITE)}: started fresh from t = 0 with closed θ walls and two passive tracers, a static phase to t {TS}, then the
full 1e-4 M<sub>⊙</sub>/yr L1 stream.</p>
<div class="kv">
<div><b>status</b>t = {P['t']:.4f} code = {P['t'] * DAYS:.2f} d = {P['t'] / PORB:.3f} P<sub>orb</sub>; {phase}; cycle {P['cycle']}, dt {P['dt']:.2e}</div>
<div><b>site / layout</b>{esc(SITE)}: {esc(LAYOUT)}</div>
<div><b>binary</b>{esc(BIN_LABEL)}</div>
<div><b>deck</b><code>{esc(DECK_LABEL)}</code> (v11am + <code>theta_wall_closed = true</code> + <code>nscalars = 2</code>; <code>am_transport = true</code>)</div>
<div><b>chain</b>{esc(CHAIN_LABEL)}; {P['segments']} run segments, {P['wall_h']:.1f} h wall so far<br><span class="status">{chain}</span></div>
<div><b>queue</b><span class="status">{q}</span></div>
<div><b>now</b>L<sub>meas</sub> {P['Lmeas']:.3g} L<sub>⊙</sub> · L<sub>top</sub> {P['Ltop']:.3g} L<sub>⊙</sub> · Ṁ<sub>in</sub> {P['Mdot_in']:.1f} · Ṁ<sub>R</sub> {P['Mdot_R']:.1f} code · Picard mean {P['Pic']:.2f} (max {P['Picmax']:.0f})</div>
<div><b>grid</b>420 × 58 × 512 (r, θ, φ), dθ = dφ = 2π/512, θ = π/2 ± 20.39°; stretched r 6.3 – 24.34 R<sub>⊙</sub> (1.5 d<sub>L1</sub>); 8 MeshBlocks 420×58×64</div>
<div><b>system</b>M<sub>acc</sub> 16, M<sub>don</sub> 18 M<sub>⊙</sub>, a 33.26 R<sub>⊙</sub>, P 3.81 d; R<sub>acc</sub> ≈ 9.0 R<sub>⊙</sub>, T<sub>eff</sub> 28.3 kK; d<sub>L1</sub> 16.23</div>
<div><b>stream</b>stream_mdot {MDOT} code (1e-4 M<sub>⊙</sub>/yr), 2-cell nozzle at r 14.5, σ<sub>z</sub> 1.0 R<sub>⊙</sub>; on at t {TS}</div>
<div><b>physics</b>3-D Roche gravity, Coriolis + centrifugal; implicit grey M1 + vet_col; MESA-kap Rosseland + Planck tables; general EOS (Saha H/He); von Zeipel IC; held light neck, donor_irr 0.35</div>
<div><b>θ walls</b><code>theta_bc = hse</code> + <code>theta_wall_closed = true</code> (exactly closed: no mass, energy or x1/x3-momentum flux); radiation <code>implicit_bc_x2 = flux</code> with the von Zeipel factor</div>
</div>
<p class="status">{esc(STATUS)}</p>
</header>
<section>
<h2>Run checks</h2>
<ul>{anomalies()}</ul>
<h3>Conservation budgets</h3>
<ul>{budgets()}</ul>
{fig('fig_budget.png', 'budgets', 'Budget residuals from ryper.user.hst vs time (code units): mass, energy (with the boundary throughput E<sub>bnd</sub>), and L<sub>z</sub> (with the nozzle input and ΔJ<sub>dom</sub>). Dotted: stream switch-on.')}
{fig('fig_prod_hist.png', 'history', 'hst history (static + stream). Row 1: dt (red = dt retries), Picard mean and max, L<sub>meas</sub> and L<sub>top</sub>. Row 2: domain mass change with the cumulative nozzle inflow, boundary outflow (negative = inflow) and floor mass; mass rates (5-row mean; dashed = the nominal stream); relative change of J<sub>z</sub>. Row 3: total energy (user and hydro hst), kinetic energy, radiation energy.')}
</section>
<section>
<h2>Movies</h2>
{movies()}
</section>
{tracer_section()}
<section>
<h2>Maps</h2>
<p>Cartesian, R<sub>⊙</sub>, donor along +x. Rows: equator (mean of the two θ rows next to it, whole domain), impact zoom, meridional φ 0/180, meridional through the
impact azimuth. Columns: log ρ; log T<sub>gas</sub> from the run's EOS table (ρ, e<sub>int</sub>); log P<sub>rad</sub>/P<sub>gas</sub> with P<sub>rad</sub> = E/3 from the
m1 dump at the same time (the dump carries E and F, not the VET tensor) and P<sub>gas</sub> from the EOS table at T<sub>gas</sub>; v<sub>φ</sub> in the rotating frame.</p>
{maps or '<p class="status">no hydro dump with a matching m1 dump yet</p>'}
</section>
{previous()}
<section>
<h2>Caveats</h2>
<ul>
<li>Half the radial resolution of the 2-D bigdom grid and a ±20° wedge: the closed θ walls stand in for the gas above and below.</li>
<li>P<sub>rad</sub> = E/3 is the isotropic (Eddington) value; in the optically thin atmosphere the VET P<sub>rr</sub> differs.</li>
<li>"On the star" (r &lt; 10.3) includes the hot impact layer; the s0 maps and profiles show where inside it the gas sits.</li>
<li>The L<sub>z</sub> residual is the external torque, not booked in the hst.</li>
</ul>
</section>
<footer>Design note: docs/dev/accretor_rhd_design.md §16. Page built by docs/handover/acc3d-page-1020/upd_page.sh (bsg-files-1009) from the run's hst, logs and dumps.
Run dirs: {esc(P['dirs'])} under the run parent on {esc(SITE)}.</footer>
</main>
</body></html>
'''
open(OUTF, 'w').write(page)
print(OUTF, os.path.getsize(OUTF))
