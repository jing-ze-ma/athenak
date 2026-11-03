#!/usr/bin/env python3
"""lowbeta-tests-1102: write pages/lowbeta_cart/index.html from results.json."""
import json
B = '/resnick/groups/carnegie_poc/jingze/lowbeta_1102'
OUT = '/resnick/home/jingze/ATHENAK/pages/lowbeta_cart/index.html'
J = json.load(open(f'{B}/results.json'))
S = J['summary']; R = J['runs']
css = open(f'{B}/page_style.css').read()
NAME = {'hlle': 'HLLE', 'hlld': 'HLLD', 'lhlld': 'LHLLD', 'gs': 'GS07', 'sm': 'emf_smooth',
        'uct': 'UCT-HLLD'}
def lab(s):
    a = s.split('_'); t = NAME[a[0]] + ' + ' + NAME[a[1]]
    return t + (' + FOFC' if len(a) > 2 else '')
def cls(o):
    return 'good' if o >= 1.8 else ('mid' if o >= 1.2 else 'bad')
def ordtab(d, betas, title_n=True):
    """d[scheme][beta] = {N, E, order}"""
    h = '<div class="tw"><table><thead><tr><th>scheme</th>'
    for b in betas: h += f'<th>beta {b}: error at finest N</th><th>orders</th>'
    h += '</tr></thead><tbody>'
    for sc, v in d.items():
        h += f'<tr><td>{lab(sc)}</td>'
        for b in betas:
            x = v.get(b)
            if not x: h += '<td class="ro">-</td><td class="ro">-</td>'; continue
            os_ = ' '.join(f'<span class="{cls(o)}">{o:.2f}</span>' for o in x['order'])
            h += f'<td class="num">N={x["N"][-1]}: {x["E"][-1]:.2e}</td><td class="num">{os_}</td>'
        h += '</tr>'
    return h + '</tbody></table></div>'
def fig(name, alt, take, cap=''):
    return (f'<div class="card"><p class="take">{take}</p><figure><img src="img/{name}.png" alt="{alt}" '
            f'loading="lazy"></figure>' + (f'<figcaption>{cap}</figcaption>' if cap else '') + '</div>')

CB = ['0.2', '0.01', '0.001', '0.0001', '1e-05']
p = []
p.append('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" '
         'content="width=device-width,initial-scale=1"><title>Low-Beta MHD Tests</title>' + css +
         '</head><body><main>')
p.append('<h1>Low-Beta MHD Tests</h1><p class="sub">Cartesian AthenaK (prod-1102 + ucthlld-1102 401c1dbee, CPU serial, '
         'RK2, PLM unless noted, cfl 0.4). Test setups only, no method changed. 921 runs, 2026-11-02. '
         'Branch <code>lowbeta-tests-1102</code> (c1c2a44ad+), data <code>/resnick/groups/carnegie_poc/jingze/lowbeta_1102/results.json</code>.</p>')
p.append('<nav><a href="#verdict">Verdict</a><a href="#cpaw">Alfven wave</a><a href="#lw">Linear waves</a>'
         '<a href="#vortex">Vortex</a><a href="#rd">Rotational disc.</a><a href="#rotor">Rotor</a>'
         '<a href="#blast">Blasts</a><a href="#setup">Setup</a></nav>')
p.append('''<div class="verdict" id="verdict"><h3>Verdict</h3><ul>
<li><b>Smooth problems with PLM stay second order down to beta = 1e-5 for every Riemann solver and corner EMF</b>
(CPAW 1-D/2-D, fast and Alfven linear waves, Christlieb vortex). The pressure error does not get worse in absolute terms as beta falls,
so relative to p it grows as 1/beta: at beta = 1e-5 the 2-D Alfven wave drives p to the floor at N &le; 128 and still dips to 6-14 % of p0 at N = 256.</li>
<li><b>The one real failure in smooth flow: HLLD and LHLLD with the default GS07 corner EMF on a slow wave</b>
(amplitude 1e-6). Order drops to 0.2-0.9 at beta = 1e-2, and at beta &le; 1e-3 (N &ge; 64) or 1e-4 (all N) the wave blows up:
p sits at the floor every step and dt falls by 1e4. HLLE + GS07 is fine; <b>emf_smooth and UCT-HLLD both cure it</b> and give second order.</li>
<li><b>PPM4 does not converge on the 2-D Alfven wave at any beta</b> (B_perp order about 0 between N = 128 and 256, even at beta = 0.2).
With HLLD/LHLLD + GS07 or emf_smooth, its pressure error grows with N below beta = 1e-2 (negative order, up to 100 p0); UCT-HLLD keeps it 30-40x smaller.</li>
<li><b>Blasts: everything finishes at beta 0.016 and 2.5e-4 (2-D and 3-D), but 0.3-1 % of cell updates hit the pressure floor at 2.5e-4.</b>
At beta 3e-6, HLLD/LHLLD with GS07 or emf_smooth (with or without FOFC) stall: dt collapses 1000x by 0.38 t_end with about 10 % of updates at the floor.
UCT-HLLD and HLLE finish; HLLE has 20x fewer floor hits than UCT-HLLD. The low-beta rotor (beta 0.1) is clean for all.</li>
<li><b>UCT-HLLD vs emf_smooth:</b> identical on CPAW-PLM, vortex, fast/Alfven waves; both cure the HLLD slow-wave blow-up; only UCT-HLLD
survives the beta 3e-6 blast and tames the PPM pressure error. Side finding: the built-in <code>cpaw</code> pgen is 3-D only (a 2-D or 1-D run starts from a garbage field).</li>
</ul></div>''')

# ---------------- CPAW
p.append('<h2 id="cpaw">1. Circularly polarised Alfven wave</h2><p>B_par = 1, B_perp = v_perp = 0.1, rho = 1, p lowered to beta = 0.2 ... 1e-5; '
         'oblique 2-D (2N x N box, wave along the diagonal) and 1-D, one period. The exact solution keeps p and |B| uniform, '
         'so any p deviation is error. Test pgen <code>lowbeta_tests iprob=3</code>, because the built-in <code>cpaw</code> builds '
         'its vector potential with a single x3 cell in 2-D but reads index k+1 (2-D runs started with 40x the magnetic energy).</p>')
p.append('<div class="stack">')
p.append(fig('cpaw2d_L1_Bperp_plm', 'CPAW 2-D PLM B_perp error vs N',
             'PLM: B_perp converges at second order at every beta, and all eight scheme pairs lie on top of each other.'))
p.append(fig('cpaw2d_L1_p_plm', 'CPAW 2-D PLM pressure error vs N',
             'PLM: the pressure error (|p - p0|) is about 4e-7 at N = 256 whatever beta is, and falls at order 2-3; low beta does not break convergence.',
             'Because the absolute error does not shrink with p0, the relative error is about 1/beta: 8 % of p0 at beta = 1e-5, N = 256.'))
p.append(fig('cpaw2d_pmin_t', 'CPAW min p vs time',
             'But at beta = 1e-5 the pressure hits the floor almost at once on N = 128 (all schemes); at beta = 1e-4 the dip is about 30 %.',
             'Min over the mesh of p, divided by the exact p, from 40 dumps per period. Floor-hit counts are in the table below.'))
p.append(fig('cpaw2d_L1_Bperp_ppm4', 'CPAW 2-D PPM4 B_perp error vs N',
             'PPM4: B_perp stops converging between N = 128 and 256 at every beta, including the standard beta = 0.2, so this is not a low-beta effect.'))
p.append(fig('cpaw2d_L1_p_ppm4', 'CPAW 2-D PPM4 pressure error vs N',
             'PPM4: with HLLD/LHLLD and GS07 or emf_smooth the pressure error grows with N below beta = 1e-2 (negative order); UCT-HLLD and HLLE stay about 30x lower.'))
p.append(fig('cpaw2d_order_L1_p', 'order vs beta', 'Measured pressure order (N = 128 to 256) against beta: PLM stays at or above 2; PPM4 turns negative.'))
p.append(fig('cpaw1d_L1_p_plm', 'CPAW 1-D', '1-D: the same picture, with far fewer floor hits (only N &le; 64 at beta = 1e-5, and only inside RK stages).'))
p.append('</div>')
# floor table
F = S['cpaw2d']['floors']
h = '<details><summary>numbers: 2-D floor hits / min p over p0 (PLM and PPM4)</summary><div class="tw"><table><thead><tr><th>recon</th><th>scheme</th>'
for b in CB: h += f'<th>beta {b} (N=32 64 128 256)</th>'
h += '</tr></thead><tbody>'
for rec in ('plm', 'ppm4'):
    for sc in ['hlle_gs', 'hlld_gs', 'hlld_sm', 'hlld_uct', 'lhlld_gs']:
        h += f'<tr><td>{rec}</td><td>{lab(sc)}</td>'
        for b in CB:
            cells = []
            for N in (32, 64, 128, 256):
                x = F.get(f'{rec}|{sc}|{b}|{N}')
                if x:
                    hits = x['efloor'] + x['tfloor']
                    cells.append(f'{hits}/{x["pmin_over_p0"]:.2g}' if hits else f'0/{x["pmin_over_p0"]:.2g}')
            h += f'<td class="num">{" ".join(cells)}</td>'
        h += '</tr>'
p.append(h + '</tbody></table></div><p class="note">hits = efloor + tfloor counts summed over all RK stages; min p/p0 from dumps (1e-33 = at the floor). emf_smooth gives the same errors as GS07 to 3 digits for the CPAW.</p></details>')
p.append('<details><summary>numbers: 2-D PLM orders, B_perp and p</summary>' + ordtab(S['cpaw2d']['plm']['L1_Bperp'], CB) +
         ordtab(S['cpaw2d']['plm']['L1_p'], CB) + '</details>')
p.append('<details><summary>numbers: 2-D PPM4 orders, B_perp and p</summary>' + ordtab(S['cpaw2d']['ppm4']['L1_Bperp'], CB) +
         ordtab(S['cpaw2d']['ppm4']['L1_p'], CB) + '</details>')

# ---------------- LW
LB = ['0.369', '0.01', '0.001', '0.0001']
p.append('<h2 id="lw">2. Oblique linear waves (amplitude 1e-6)</h2><p>AthenaK <code>linear_wave</code> (2-D, 2N x N), '
         'background B = (1, sqrt 2, 0.5), p lowered to beta 0.37 / 1e-2 / 1e-3 / 1e-4. The pgen recomputes the eigenvectors '
         'from the background state itself (MHDEigensystemPrim on dens, pgas, b0), so no deck change is needed. '
         'Slow waves at low beta are long (period up to 150), so the slow ladder stops at N = 128 (beta 1e-2, and 3 schemes at 1e-3) or 64.</p>')
p.append('<div class="stack">')
p.append(fig('lw_order', 'linear wave order vs beta',
             'Fast and Alfven waves keep second order down to beta = 1e-4 with every scheme (1.7-1.9 at the finest pair at 1e-4).'))
p.append(fig('lw_slow', 'slow wave error vs N',
             'Slow wave: HLLD/LHLLD + GS07 loses order at beta = 1e-2 and blows up from round-off at beta &le; 1e-3; emf_smooth and UCT-HLLD restore clean second order.',
             'The blown-up runs hit the pressure floor first at cycle 6400-8100 and then every step; dt fell from 1e-2 to 1e-6 and they were stopped.'))
p.append(fig('lw_fast', 'fast wave', 'Fast wave: all schemes identical and second order.'))
p.append(fig('lw_alf', 'Alfven wave', 'Alfven wave: all schemes identical and second order.'))
p.append('</div>')
p.append('<details><summary>numbers: slow-wave orders</summary>' + ordtab(S['lw']['slow'], LB) + '</details>')

# ---------------- vortex
V = S['vortex']
p.append('<h2 id="vortex">3. Christlieb low-pressure magnetised vortex</h2><p>Christlieb et al. 2015 sect. 5.2.1: mu = 5.389489439, kappa = sqrt 2 mu, '
         'centre p = 5.3e-12, domain [-10, 10]^2, t = 0.05, built-in <code>balsara_vortex</code> (same field and pressure formula; '
         'its advection velocity is (0.86, 0.86) instead of (1, 1)). Errors against the shifted exact solution, cell averaged.</p>')
p.append('<div class="stack">' + fig('vortex', 'vortex errors and min p',
         'Second order in v and p from N = 32 to 1024 for all schemes; min p stays positive in every output and follows the exact cell-averaged minimum.',
         'The floor still fires inside RK stages at every N &ge; 64 (3 hits at 64, 463 at 1024: more cells sit near the centre); the schemes are indistinguishable.') + '</div>')
h = '<details><summary>numbers: floor hits and min p vs N (HLLD + GS07; all schemes equal)</summary><div class="tw"><table><thead><tr><th>N</th><th>L1 p</th><th>order</th><th>min p</th><th>floor hits</th></tr></thead><tbody>'
d = V['hlld_gs']
for i, N in enumerate(d['N']):
    o = d['order_p'][i-1] if i else None
    h += f'<tr><td>{N}</td><td class="num">{d["L1_p"][i]:.2e}</td><td class="num">{"" if o is None else f"{o:.2f}"}</td><td class="num">{d["pmin"][i]:.1e}</td><td class="num">{d["floor_hits"][i]}</td></tr>'
p.append(h + '</tbody></table></div></details>')

# ---------------- RD
p.append('<h2 id="rd">4. Low-beta rotational discontinuity (1-D)</h2><p>Zhang, Lani &amp; Poedts 2026 sect. 6.3: rho = p = 0.1, Bx = 5, By = &plusmn;0.1, gamma = 2, t = 0.01, '
         'N = 200 / 800 / 3200 plus 25600 as reference. No floor hits anywhere.</p>')
rd = S['rd']
p.append('<div class="stack">' + fig('rd_profiles', 'RD density and internal energy profiles',
         'HLLD and LHLLD leave a density dip at x = 0.5 of 1.5 % that does not shrink with N (0.09927, 0.09852, 0.09847, 0.09844); HLLE leaves 0.1-0.3 %.',
         'The internal-energy rise behind the fast waves (about 0.5 %) is the same for all three and narrows with N. LHLLD = HLLD to 1e-8 here.') + '</div>')

# ---------------- rotor
p.append('<h2 id="rotor">5. Low-beta rotor</h2><p>Zhang 2026 sect. 6.5: ambient rho = 0.1, rotor rho = 1, p = 0.05, B = 2.5/sqrt(4 pi) (1, 1) '
         '(ambient beta 0.1), gamma = 5/3, t = 0.12, N = 200 / 400 / 800, outflow. Test pgen <code>lowbeta_tests iprob=1</code>.</p>')
p.append('<div class="stack">' + fig('rotor_pmin_t', 'rotor min p vs time',
         'No negative internal energy, no floor hit and no failure for any scheme at any N; min p stays at 5e-3 or more (10 % of the initial p).',
         'Zhang et al. see negative e with HLLD in a first-order implicit code without divergence cleaning; with CT it does not occur here.') +
         fig('rotor_maps', 'rotor pressure maps', 'At N = 400 the pressure maps of GS07, emf_smooth and UCT-HLLD look the same.') + '</div>')

# ---------------- blasts
p.append('<h2 id="blast">6. Blast ladder (robustness only)</h2><p>beta 0.016: Tremblin 2024 sect. 8.2.3 (= GS05 blast; [0,1]^2 periodic, p 10/0.1, gamma 5/3, t = 0.2). '
         'beta 2.5e-4: Balsara &amp; Spicer / MD21 sect. 6.4 (B0 = 100/sqrt(4 pi) at 45 deg, p 1000/0.1, gamma 1.4, t = 0.01, outflow), 2-D N = 256 and 3-D N = 64, 128. '
         'beta 3e-6: Tremblin sect. 8.2.4 (B = 250/sqrt2 (1, 1), p 1000/0.1, gamma 1.4, t = 0.02). Test pgen <code>lowbeta_tests iprob=2</code>. '
         'FOFC is off by default for these pgens; two FOFC arms were added. pfloor = tfloor = FLT_MIN (the defaults), so a floor hit means p would have gone to zero or below.</p>')
p.append('<div class="stack">' + fig('blast_floors', 'blast floor fractions',
         'All arms finish at beta 0.016 and 2.5e-4. At 2.5e-4 about 0.3 % (2-D) to 1 % (3-D, N = 64) of cell updates need the pressure floor, whatever the scheme; FOFC does not change that.',
         'Fraction = efloor count / (2 stages x cycles x cells).') +
         fig('blast_dt', 'beta 3e-6 blast dt', 'beta 3e-6: HLLD/LHLLD with GS07 or emf_smooth (and HLLD + FOFC) lose dt 1000x by 0.38 t_end; UCT-HLLD and HLLE run through.') +
         fig('blast_emin_t', 'blast min p vs time', 'At beta 2.5e-4 and 3e-6 some cell sits at the pressure floor in almost every output from the first step on, for every scheme (HLLE recovers now and then late at 3e-6); at beta 0.016 the floor is touched only briefly, early and near 0.8 t_end.') + '</div>')
h = '<details><summary>numbers: blasts</summary><div class="tw"><table><thead><tr><th>run</th><th>finished</th><th>cycles</th><th>floor fraction</th><th>efloor+tfloor</th><th>FOFC cells</th><th>wall s</th></tr></thead><tbody>'
for k, v in S['blast'].items():
    h += (f'<tr><td>{k}</td><td class="{"good" if v["completed"] else "bad"}">{"yes" if v["completed"] else "stalled"}</td>'
          f'<td class="num">{v["cycles"]}</td><td class="num">{v["floor_fraction"]:.1e}</td><td class="num">{v["floor_hits"]}</td>'
          f'<td class="num">{v["fofc"]}</td><td class="num">{v["wall"]}</td></tr>')
p.append(h + '</tbody></table></div></details>')

# ---------------- setup
p.append('''<h2 id="setup">Setup and what fired</h2><ul>
<li>Schemes: rsolver hlle / hlld / lhlld (with <code>lhlld_cb_cont = true</code>); corner EMF GS07 (default), <code>emf_smooth = 1e-4</code>, <code>emf_scheme = uct_hlld</code> (UCT-HLLD needs hlld/lhlld, 2-D/3-D, no FOFC). PLM everywhere, PPM4 also for the CPAW. RK2, cfl 0.4, nghost 4.</li>
<li>Defaults: FOFC off (no fofc_pgen_default in these pgens), pfloor = dfloor = tfloor = FLT_MIN, no velocity ceiling. Only the energy/temperature floor ever fired (dfloor, vceil and c2p failures were 0 in every run); counts come from the event log (all RK stages).</li>
<li>Binaries: CPU serial builds of 401c1dbee (built-in pgens; and <code>-D PROBLEM=lowbeta_tests</code>). Runs on hpc-sm-01-11 (alloc 5274196), 10-12 cores.</li>
<li>Scripts: <code>inputs/tests/lowbeta/</code> on <code>lowbeta-tests-1102</code> (gen.py decks, analyze.py metrics, plot.py figures).</li></ul>''')
p.append('</main></body></html>')
open(OUT, 'w').write('\n'.join(p))
print('wrote', OUT)
