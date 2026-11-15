#!/usr/bin/env python3
"""idort-1019: content.json for build19.py"""
import json, os, sys
D = '/raven/ptmp/jinma/idort1020/page'
conv = json.load(open(D + '/conv_rows.json'))
F = json.load(open(D + '/fig/figures.json'))
cs = json.load(open(D + '/csconv.json')) if os.path.exists(D + '/csconv.json') else {}
EXTRA = json.load(open(D + '/extra_notes.json')) if os.path.exists(D + '/extra_notes.json') else {}
def fo(k, cls='all', i=-1):
    v = cs.get(k); return ('%.2f' % v['order'][cls][i]) if v and cls in v['order'] else 'pending'
S = []
def row(test, geom, metric, value, a=None, b=None, ref=None):
    ref, verdict = (ref, a) if b is None else (a, b)   # row(..., value, ref, verdict) or row(..., *pv(k), ref=...)
    S.append(dict(test=test, geom=geom, metric=metric, value=value, ref=ref, verdict=verdict, v10=''))
# idort-1020: single PROD arm on 563ab474 (+ deck rsg_E_start1020r keys). Rows that FATAL'd on 563ab474 with
# time2_reject = retry (no BE fallback) are rerun on the train bin with the fallback (T2BIN); their bin is named in the row.
T2BIN = os.environ.get('T2BIN', 'athena_cpu_82275582 (md5 c55079c1, train 82275582 = + t2befb-1020 BE fallback in time2_reject = retry)')
PEND = json.load(open(D + '/t2rows.json')) if os.path.exists(D + '/t2rows.json') else {}
def pv(k, d='rerun pending (time2_reject = retry FATAL on 563ab474)'): return PEND.get(k, (d, 'PENDING'))
row('Crossing beams r4/r12/c0.5/c32 (IDORT 4.1.1)', 'Cartesian', 'median row energy / peak at y 0.65', *pv('cb'), ref='row >= 0.9, peak >= 0.3 E_b')
row('Crossing beams 4x resolution (r4n4x, 256x1024)', 'Cartesian', 'median row energy / peak', *pv('cbn4x'), ref='row >= 0.9, peak >= 0.3 E_b')
row('Crossing beams, Gauss 8x32 (cbg)', 'Cartesian', 'row energy', *pv('cbg'), ref='completes; judge vs J_sc')
row('Crossing beams cbh variant', 'Cartesian', 'row / peak, cost', *pv('cbh'), ref='completes; row >= 0.9')
row('Homogeneous sphere/cylinder', 'Cartesian', 'E, F outside', 'E L1 1.4 %, F 1.1 % (256^2); f_xx(0.5) 0.3331', 'within 10 %', 'PASS')
row('Homogeneous sphere T6b / thin sphere T6a', 'sp wedge', 'E/J-1', 'T6b -0.06 % inside, -1.2 % at r > 3 (N_r 256); T6a 4e-6', 'converging', 'PASS')
row('Non-LTE scattering atmosphere', 'Cartesian', 'L1 vs Eddington analytic', '1.11 % (eps 1e-1) ... 0.01 % (1e-8)', 'L1 <= 2 %', 'PASS')
row('Non-LTE scattering atmosphere, eps 1e-2 (nlte2 r0/w 1000)', 'sp wedge', 'E vs exact per tau band', '-0.28 % / -0.03 % in tau < 0.1 / 0.1-1, E_top sp/cart 0.994 (vet_source = relaxed_esc)', '<= 1 %', 'PASS')
row('Static shadow ring 12x2 / Gauss (shd, shdg)', 'Cartesian', 'L1 of x = 0.4 cut', '0.046 / 0.137 (plotted cut)', 'L1 <= 2 % (suite metric)', 'OPEN')
row('Spherical shadow T7', 'sp wedge', 'lit L1 r 4.05/4.99/5.77', *pv('t7'), ref='lit L1 <= 0.02')
row('Static diffusion', 'Cartesian', 'L1 vs analytic', *pv('sdif'), ref='L1 <= 1 %')
row('Spherical diffusion T4', 'sp wedge', 'L1 vs exact', '7.2e-3 / 1.8e-3 / 4.6e-4 (order 1.99)', '2nd order', 'PASS')
row('Dynamic diffusion ddf / ddn', 'Cartesian', 'peak, centroid, D_eff', *pv('ddf', 'pending (figure extract waits for sdifk)'), ref='peak within 5 %')
row('Thermal equilibrium 1 / 2', 'Cartesian', 'final E, T', 'E 96.79505 T 3.13663 / E 145.78780 T 3.47480', 'targets', 'PASS')
row('Thermal equilibrium (relax1/2)', 'sp wedge', 'error vs dt', '2.7e-9 at dt/1, floor by dt/8', 'floor', 'PASS')
row('Moving thermal equilibrium', 'Cartesian', 'E, v, P_tot', 'E 1.13260, v 2.95538, P_tot 3.000000', 'E 1.1326, v 2.9554', 'PASS')
row('Radiation shocks M 1.2 / 2 / 3', 'Cartesian', 'L1 rho/Tg/Tr vs IDORT transport', 'M1.2 0.0003; M2 (512) 0.0003-0.0004, (1024) 0.0006-0.0009; M3 0.0010-0.0018', 'L1 <= 2 %', 'PASS')
row('Radiative shocks rsx (LE08) subset', 'Cartesian + sp', 'completion', *pv('rsx', 'pending (n1024 / n2048 running)'), ref='completion')
row('Irradiated cloud (IDORT beam)', 'Cartesian', 'halo edge, dense cloud', 'halo -0.208/-0.280/-0.392, cloud -0.067..0.104 / ...', 'qualitative vs IDORT', 'PASS')
row('Marshak wave', 'Cartesian', 'L1(E) vs S_N', *pv('marsh'), ref='L1 <= 0.02')
row('Marshak wave T5', 'sp wedge', 'self-convergence', 'X rho 1.09, E 0.56; T 0.2-0.4 (front)', '-', 'OPEN')
row('Crooked pipe', 'Cartesian', 'probes t=10', '0.469/0.428/0.386/0.343/0.275 keV (c dt/dx 220), 2200 within 0.3 %', 'dt self-consistency', 'PASS')
row('G5 2-D pulse', 'Cartesian', 'd sigma^2/dt / 2D', '0.99957 / 0.99987 / 1.00000 / 1.00000', '1 +- 0.02', 'PASS')
row('P1 / P1T grey atmospheres', 'sp wedge', 'E/J-1 at R_ph', '-0.02 / -0.06 / -0.03 % (N 64/128/256); patmT X E 1.86, rho 1.98', 'converging', 'PASS')
row('P2 thermal channel beam (pinc)', 'sp wedge', 'L/L_MC', *pv('pinc'), ref='within 5 %')
row('P2 n40 decks (HFfZR0G2)', 'sp wedge', 'run', 'FATAL "vet_gd: too many directions" (the deck asks for 20 x 24 directions, more than the build allows); same on 10-19', '-', 'OPEN')
row('Large-radius controls (pulse, rsw/rswt r1000)', 'sp vs Cartesian', 'L1 sp-cart', 'pulse 8.8e-2 (r0/w 1) / 3.7e-4 (1000); rsw 1.4e-4, rswt 2.2e-4', 'r0/w 1000 -> 0', 'PASS')
row('Lateral waves rlat/rlatt (+om5)', 'sp wedge', 'X / T orders, L_z', 'X 1.9-2.4 (rlatt F1 0.77, v 0.90); T 1.66-1.75; rlatt_om5_T E 1.27', '2nd order', 'OPEN')
row('Ladders rw_t10a / rw_p10t1 (time refinement)', 'Cartesian', 'T ladder finest pair', 'error JUMPS at the smallest dt (rw_t10a T E pairs 1.23 3.83 2.77 -8.59, F1 -16.2; rw_p10t1 F1 -6.46, rho -3.71); rw_t10a C ladder E 0.61 / 0.60. Cause on 10-19: the production tolerance floor (implicit_tol 1e-8, lin_tol 1e-10, lin_cnorm 1e-10) over many small steps', '2nd order', 'FAIL')
row('Ladders rw_t1000a, rw_p10t100, pulse k*', 'Cartesian', 'X / C / T orders', 'X, C, T 1.8-2.4 every variable', '2nd order', 'PASS')
summary_note = ('PROD = the suite decks with every <rad_m1> method key and the <time> dt_retry keys of the production RSG cs deck '
                '(rsg_g25_1019/prod/R2/rsg_E_start1020r.athinput) set (keys table below), plus <hydro> am_transport on cs / multi-D sp decks and '
                '<mesh> cs_seam_flux = shared on cs decks. Single arm, one binary (rows rerun on the BE-fallback bin say so). Verdicts use the v10 criteria.')
regress = []
flags = [
 dict(title='cs / sp diffusion mode in an atmosphere (gda): order -1.6 then 1.0 with PROD', text='Space pairs on cs -1.60 / 1.03 (all, seam and vertex cells alike, cs = sp), time-only at N16 1.74 / 2.83. Earlier bisect (10-19, sp): only the PROD half-range / blend / thin face-flux keys (implicit_flux*, implicit_blend*, implicit_thin_*) reproduce it; not a seam / cs issue. The thick diffusion mode (gd) is 2.0 on cs and sp.'),
 dict(title='sp moving shell rswv / rswtv: 0.4-0.8', text='The suite decks are the OLD rswv decks. With the user-approved rswv deck change (implicit_f0_init = zero, ox1_bc = outflow) rswv_F is 1.86-1.91 (E) / 1.54-1.87 (rho) at N 32-256 (dt ~ dx). rswtv (semi-thin) stays 1.05-1.33 with the same fix.'),
 dict(title='sp rswt (semi-thin wave): E order 1.97 -> 1.77 -> 1.37 at N 128-512', text='Spatial, not temporal: the dt-only sweep at N 256 gives 1.83 / 1.82 down to the solver floor (1e-8).'),
 dict(title='Cartesian radwave 2-D (xy) P 3.3e-3 tau 1: E finest pair 1.31', text='Spatial (10-19 dt-only / space-only sweeps); rho is 1.98. 2-D only (x1 is 2.0).'),
]
keys = [kv.split('/', 1)[1].split('=', 1) for kv in open('/raven/ptmp/jinma/idort1020/tool/prodkeys_all.txt').read().split()]
conv_figs = [dict(file='conv_rw_cart.png', caption='Cartesian vet_sc radiation-modified acoustic wave, L1(rho) vs the exact damped eigenmode after one period, N 16-256 with nt = 2N (dt ~ dx), 6 regimes; PROD. Dashed: slope 2.'),
             dict(file='conv_rw_cart_E.png', caption='Same, L1(E).'),
             dict(file='conv_sp.png', caption='sp wedge (vet_col + vet_gd, production closure) linear-mode ladders, 4-5 levels: rsw/rswt (+ moving shell, + user-approved deck _F), fixed dt (X) and dt ~ dx (_C), lateral waves (16-256), grey atmosphere, and the rswt dt-only sweep. PROD.'),
             dict(file='conv_cs.png', caption='Cubed sphere and sp: radiation diffusion mode l = 2 (thick; atmosphere tau 20 with Marshak top; Eddington reference), space (N 8-64 panel cells, dt ~ dx) and time (N 16, dt/q), error classes all / seam / vertex on cs.'),
             dict(file='conv_lad_cart.png', caption='Cartesian suite ladders (rw_*, pulse_k*): X (space, fixed dt), C (dt ~ dx), T (dt refinement); PROD.')]
order = ['01_cb', '02_cbvar', '03_cbn4x', '04_cbg', '05_hsphere', '05_hsphere_sp', '06_nlte', '06_nlte_sp', '07_shd', '07_shd_sp', '08_shdg', '09_sdif', '09_sdif_sp', '10_ddf', '11_ddn', '13_eq1', '13_eq1_sp', '14_eq2', '14_eq2_sp',
         '15_eq3t', '18_shock', '19_cloud', '20_marsh', '20_marsh_sp', '21_pipe', '25_patm_sp', '26_patmT_sp', '27_pinc_sp', '28_pulse_lr_sp', '29_rsw_lr_sp', '30_g5',
         'ladder_sp_rsw', 'ladder_sp_rswt', 'ladder_sp_rswtv', 'ladder_sp_rlat', 'ladder_sp_rlat_om5', 'ladder_sp_rlatt', 'ladder_sp_patmT', 'ladder_sp_marsh_r1', 'ladder_sp_relax1_T', 'ladder_sp_exact']
import re
def clean(c):
    c = re.sub(r'\(FIN2 control[^)]*\)', '', c); c = re.sub(r'Caveat \(IN\).*$', '', c); c = c.replace('FIN3 final', 'PROD')
    c = re.sub(r'suite L1 \? and umbra \? \([^)]*\);', '', c); c = re.sub(r'Before \(FIN2 control[^.]*\.[^.]*\.', 'v10 (FIN3) 120 dirs: 0.004 / 0.005 / 0.005.', c)
    return c.strip()
other = [dict(title=k, file=F[k]['file'], caption=clean(F[k].get('caption', ''))) for k in order if k in F and F[k].get('file') and os.path.exists(D + '/fig/' + F[k]['file'])]

SEC = [
 dict(id='s411', nav='4.1.1 beams', title='4.1.1 Crossing beams', ref='cf. Ma+ Sect. 4.1.1, Fig. 2 (J21 4.2)', figs=[('01_cb.png', 'cf. Ma+ Fig. 2: E maps (exact / PROD / formal solution) and beam energy per row vs exact; ring 4, 64x256, c dt/dx 256.')],
      text='Two beams crossing in vacuum; the scheme must keep each beam narrow and its energy constant along the path.', tests=('Crossing beams',)),
 dict(id='s412', nav='4.1.2 sphere', title='4.1.2 Homogeneous sphere', ref='cf. Ma+ Sect. 4.1.2, Fig. 3 (J21 4.6)', figs=[('05_hsphere.png', 'cf. Ma+ Fig. 3: E, F and P_xx/E vs the exact ray solution (Cartesian homogeneous cylinder).')],
      text='An emitting optically thick body in vacuum: interior equilibrium, 1/r^2 fall-off outside.', tests=('Homogeneous sphere',)),
 dict(id='s413', nav='4.1.3 NLTE', title='4.1.3 Non-LTE scattering atmosphere', ref='cf. Ma+ Sect. 4.1.3 (J21 4.5)', figs=[('06_nlte.png', 'cf. Ma+ Sect. 4.1.3 figure: E(tau) for eps = 1e-1 ... 1e-8 vs the Eddington analytic solution.')],
      text='Scattering atmosphere with thermalisation parameter eps.', tests=('Non-LTE',)),
 dict(id='s414', nav='4.1.4 shadow', title='4.1.4 Static shadow', ref='cf. Ma+ Sect. 4.1.4, Fig. 5 (J21 4.3)', figs=[('07_shd.png', 'cf. Ma+ Fig. 5: E maps (exact, PROD, formal solution) and a cut at x = 0.4.')],
      text='An opaque clump shadowing a beam.', tests=('shadow',)),
 dict(id='s421', nav='4.2.1 diffusion', title='4.2.1 Dynamic diffusion', ref='cf. Ma+ Sect. 4.2.1 (J21 4.4)', figs=[('10_ddf.png', 'cf. Ma+ Sect. 4.2.1: advected diffusing pulse vs the analytic moving Gaussian, t = 16.')],
      text='A radiation pulse advected by moving gas while it diffuses.', tests=('Dynamic diffusion',)),
 dict(id='s422', nav='4.2.2 bubbles', title='4.2.2 Multi-source bubbles', ref='cf. Ma+ Sect. 4.2.2, Fig. 7', figs=[], text='EXCLUDED from the suite by user decision (2026-10-13); not rerun.'),
 dict(id='s431', nav='4.3.1 equilibrium', title='4.3.1 Thermal equilibrium', ref='cf. Ma+ Sect. 4.3.1, Fig. 8 (J21 4.1)', figs=[('13_eq1.png', 'cf. Ma+ Fig. 8, case 1: (T, E) = (1, 100) relaxation vs the ODE.'), ('14_eq2.png', 'cf. Ma+ Fig. 8, case 2: (T, E) = (100, 1).'), ('15_eq3t.png', 'Case 3: moving equilibrium v = 3, c = 10.')],
      text='Gas-radiation energy exchange in a uniform medium, at rest and moving.', tests=('Thermal equilibrium', 'Moving thermal')),
 dict(id='s432', nav='4.3.2 linear waves', title='4.3.2 Linear radiation-modified acoustic waves (second-order convergence)', ref='cf. Ma+ Sect. 4.3.2',
      figs=[('conv_rw_P001.png', 'cf. Ma+ Sect. 4.3.2, gas-dominated case: L1 error of rho (top) and E (bottom) after one period vs the exact damped eigenmode, N = 16-256 cells per wavelength, nt = 2N steps per period (dt ~ dx), tau_lambda = 0.01 / 1 / 100. PROD 1-D (Cartesian vet_sc) and PROD 2-D diagonal wave. Dashed: slope 2; p = order of the finest pair.'),
            ('conv_rw_P10.png', 'cf. Ma+ Sect. 4.3.2, radiation-dominated case: as above.'),
            ('conv_rw_ordermap.png', 'Extended grid: finest-pair order (N 128 -> 256) of rho and E for every (P_rad/P_gas, tau_lambda); green = second order.'),
            ('conv_rw_grid.png', 'Extended grid: L1 error vs N for each combination (rho circles, E triangles); dashed slope 2.')],
      text='Second order in rho in all six 1-D regimes; the 2-D diagonal wave E at P 3.3e-3, tau 1 is 1.31 (see the flags). Extended grid (10-21, production method, binary 4de8168d, x1, N 16-256, dt ~ dx): P_rad/P_gas 3.3e-3 - 3.33 x tau 0.01 - 1000. rho is second order except P 3.33 tau 10 (1.12). E LOSES ORDER IN THE THIN REGIME: finest-pair 0.66 / 0.66 / 0.76 at tau 0.01 and 1.55 / 1.58 at tau 0.1 (P <= 0.33); P 3.33 tau 10 E 1.37. P_rad/P_gas 33 is omitted: with c_light = 10 it has E/(rho c^2) = 1, outside the Newtonian regime of the test (the conservative2 write-back loop under time2_vstage diverges there).', tests=()),
 dict(id='s433', nav='4.3.3 shocks', title='4.3.3 Radiating shocks', ref='cf. Ma+ Sect. 4.3.3 (J21 4.8, JSD12 6.2)', figs=[('18_shock.png', 'cf. Ma+ Sect. 4.3.3: rho, T_gas, T_rad profiles of the M = 1.2 / 2 / 3 shocks vs the IDORT transport reference.')],
      text='Sub-critical radiating shocks with a radiative precursor.', tests=('Radiation shocks', 'rsx')),
 dict(id='s434', nav='4.3.4 cloud', title='4.3.4 Irradiated cloud', ref='cf. Ma+ Sect. 4.3.4, Fig. 14 (J21 Fig. 5)', figs=[('19_cloud.png', 'cf. Ma+ Fig. 14: density of the irradiated cloud at t = 0.06 / 0.10 / 0.16 with the IDORT beam.')],
      text='A dense cloud ablated by a strong beam (RHD).', tests=('Irradiated cloud',)),
 dict(id='sj', nav='Jiang tests', title='Jiang 2012 / 2021 tests not in the paper', ref='J21 4.4; JSD12 5.1, 5.2; G5 pulse', figs=[('09_sdif.png', 'J21 Sect. 4.4: static diffusion of a Gaussian pulse vs analytic.'), ('20_marsh.png', 'JSD12 Sect. 5.1: Marshak wave vs the S_N reference.'), ('21_pipe.png', 'JSD12 Sect. 5.2: crooked pipe, T_gas map and probe histories.'), ('30_g5.png', 'G5 2-D pulse: variance growth vs the exact moment solution.')],
      text='Cartesian tests from Jiang, Stone & Davis 2012 and Jiang 2021 that the IDORT paper does not repeat.', tests=('Static diffusion', 'Marshak wave', 'Crooked pipe', 'G5')),
 dict(id='sprod', nav='sp / cs', title='Production geometry: spherical wedge (vet_gd) and cubed sphere', ref='production closure vet_col + vet_gd, same PROD keys', figs=[('conv_sp.png', 'sp wedge: self-convergence of the linear radiation-acoustic modes (4-5 levels, dt ~ dx). Dashed: slope 2.'), ('conv_cs.png', 'Cubed sphere (all / seam / vertex cells) and sp wedge: l = 2 radiation diffusion mode, space (N 8-64, dt ~ dx) and time (N 16, dt/q). Dashed: slope 2.')],
      text='The spherical counterparts of the paper tests run on the production path; their numbers are in the summary table.', tests=('sp wedge',)),
]
DIAG = [['all suite run dirs (PROD arm p/; p_563/ = the 563ab474 attempt of the rerun rows)', '/raven/ptmp/jinma/idort1020/runs'], ['IDORT-row evaluations (eq, shocks, sp xord/tord/exact/lr)', '/raven/ptmp/jinma/idort1020/eval/out_p'],
        ['Cartesian radwave runs', '/raven/ptmp/jinma/idort1020/rw/runs'], ['cs/sp diffusion-mode runs', '/raven/ptmp/jinma/idort1020/cs/runs'],
        ['convergence data (json)', '/raven/ptmp/jinma/idort1020/page/{rwconv,spxconv,csconv,ladconv,conv_rows}.json'],
        ['time2_reject = retry non-admissible stage diagnostic (dbg_t2_admiss, ss/n32)', '/raven/ptmp/jinma/idort1020/diag/ss32']]
C = dict(sections=SEC, diag=DIAG, lede='Full IDORT suite (Ma, Pakmor, Justham & de Mink 2025, arXiv 2503.16627) plus the Jiang 2012/2021 additions, rerun with the updated production keys (RSG deck of 2026-10-20: vet_source = relaxed_esc, implicit_am_cs, time2_reject = retry with dt_retry) on the train-tip binary 563ab474: Cartesian vet_sc, spherical-wedge vet_gd and cubed-sphere variants, PROD arm only.',
         meta=dict(binary='athena_cpu_563ab474 (plain CPU, md5 0c9cbaca), /raven/ptmp/jinma/amdef_1019/bin; rerun rows: ' + T2BIN, commit='prodmerge-1019b train tip 563ab474 (5864f189 + aaa01b3f relaxed_esc + 2506e3b8 time2_reject retry + ee24be34 relaxed_esc T* fix)',
                   built='2026-10-20 12:18 CEST (gcc 13.3 + openmpi 5.0, Release, Kokkos OpenMP)',
                   runs='/raven/ptmp/jinma/idort1020 (runs/, rw/, cs/, eval/); Raven CPU, 2026-10-20', arm='p = deck <- args <- PROD keys (all in the deck, no command-line keys)'),
         conv_note='Order = log2 of successive L1 errors (finest pair decides PASS >= 1.8 / FLAG). Radwave errors are vs the exact eigenmode; sp, cs and ladder errors are self-convergence differences (restricted fine - coarse).',
         conv_figs=conv_figs, conv_table=[r for r in conv if r.get('arm') in ('PROD19', 'P', 'p19') and not r['case'].startswith(('rswv:', 'rswtv:', 'rswv_C', 'rswtv_C', 'marsh'))], flags=flags, summary=S, summary_note=summary_note,
         other_note='',
         other=other, regress=regress, keys=keys,
         keys_note='From /raven/ptmp/jinma/rsg_g25_1019/prod/R2/rsg_E_start1020r.athinput: every <rad_m1> key and <time> dt_retry / dt_retry_fac_t2 / dt_retry_grow / dt_retry_ngood, minus problem / unit / opacity / BC / grid / IC keys (tool/mkprod.py). Applied to every test deck after its own args.',
         drops_note='Dropped only where the test forbids them (logged per run in prod_drop.txt): implicit_marshak_face / implicit_face_weight (sp-only) on Cartesian decks; implicit_op_stencil on periodic-x1 decks; implicit_vimp / implicit_cs2_vimp* when the test has gas_feedback, coupling or dbg_gas_force off; vet_gd_gl_nmu / nphi when the test sets its own quadrature; <hydro> am_transport on meshes it does not support (Cartesian, 1-D sp, 2-D sp with sp_cart_polar_momentum); a test alias vet_col_source gives way to vet_source. Hydro fofc not applied (needs nghost 3).')
C.update(EXTRA)
# USER 10-21: show only the tests that appear in the IDORT paper (Ma+ Sect. 4; full-suite builder = mkcontent.py.1021full)
PAPER = ('Crossing beams', 'Homogeneous sphere/cylinder', 'Non-LTE scattering atmosphere', 'Static shadow', 'Dynamic diffusion',
         'Thermal equilibrium 1 / 2', 'Radiation shocks', 'Irradiated cloud')
C['summary'] = [r for r in S if r['test'].startswith(PAPER) and r['geom'] == 'Cartesian']
C['sections'] = [s for s in SEC if s['id'] not in ('sj', 'sprod')]
for s in C['sections']:
    if s['id'] == 's431':
        s['figs'] = [f for f in s['figs'] if f[0] != '15_eq3t.png']; s['tests'] = ('Thermal equilibrium 1',); s['text'] = 'Gas-radiation energy exchange in a uniform medium at rest.'
    if s['id'] == 's433': s['tests'] = ('Radiation shocks',)
C['conv_table'] = [r for r in C['conv_table'] if r['case'].startswith('radwave')]
C['flags'] = [f for f in flags if f['title'].startswith('Cartesian radwave')]
C['lede'] = ('IDORT paper tests only (Ma, Pakmor, Justham & de Mink 2025, arXiv 2503.16627, Sect. 4), rerun with the updated production keys '
             '(RSG deck of 2026-10-20: vet_source = relaxed_esc, implicit_am_cs, time2_reject = retry with dt_retry) on the train-tip binary 563ab474, '
             'Cartesian vet_sc, PROD arm only. Jiang-only tests, spherical / cubed-sphere counterparts and ladders are not shown.')
json.dump(C, open(D + '/content.json', 'w'), indent=1)
print('content: summary', len(S), 'conv', len(C['conv_table']), 'other', len(other))
