"""Tables (out/tables.md), plots (plots/*.png) and out/web.json for the Betelgeuse wind budget
on the Dent et al. 2026 SEM (arXiv:2608.19339v2 App. B).  usage: ../grains/venv/bin/python
analyse26.py   (needs out/gas_*.npz from gas.py, out/dust_*.pkl from dust.py, out/old24.npz
from dump_old24.py)"""
import os
import json
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import prof as P  # noqa: E402
import wind as W  # noqa: E402
import grainlib_bet as gl  # noqa: E402

B = '/orion/ptmp/jinma/rsg_wind_1008/betelgeuse26/'
O, PL = B + 'out/', B + 'plots/'
os.makedirs(PL, exist_ok=True)
XS = [1.15, 1.2, 1.3, 1.5, 2., 3., 5., 10., 20., 30.]
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
INK, INK2 = '#0b0b0b', '#52514e'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.grid': True,
                     'grid.color': '#e4e3df', 'grid.linewidth': 0.6, 'lines.linewidth': 2,
                     'legend.frameon': False, 'figure.facecolor': '#fcfcfb',
                     'axes.facecolor': '#fcfcfb'})
L = []
VTS = W.VTS26
TEXTS = ('harper', 'plaw', 'const')


def at(x, y, xs=XS):
    return [float(np.interp(np.log(v), np.log(x), y)) if x[0] <= v * 1.0001 and
            v <= x[-1] * 1.0001 else np.nan for v in xs]


def f(v, n=3):
    return '-' if v is None or not np.isfinite(v) else f'{v:.{n}g}'


def row(*c):
    L.append('| ' + ' | '.join(str(x) for x in c) + ' |')


def hdr(n):
    L.append('|' + '---|' * n)


def topaxis(ax, R, lab='r / R*  (R* = 812 Rsun)'):
    t = ax.secondary_xaxis('top', functions=(lambda r: r / R, lambda x: x * R))
    t.set_xlabel(lab)
    t.set_xticks([1, 1.5, 2, 3, 5, 10, 20, 30])
    t.set_xticklabels(['1', '1.5', '2', '3', '5', '10', '20', '30'])


st = P.stellar('C', 17.5)
R = st['R']
x = W.xgrid(1.15, 30., 600)
WS = {}
for t in TEXTS:
    P.TEXT = t
    WS[t] = W.steady26(st, x=x, Tmode=t)
P.TEXT = 'harper'
w = WS['harper']
g0 = P.G * st['M'] / R**2

# ------------------------------------------------------------------ parameters
L.append('# Betelgeuse wind budget on the Dent et al. 2026 SEM (arXiv:2608.19339v2, App. B)\n')
L.append('## Parameters\n')
row('quantity', 'value', 'source / note')
hdr(3)
for q, v, s in (
        ('d', '172 pc', 'Dent+2026 App. B (MacLeod+2025)'),
        ('R* (Rosseland tau = 2/3)', f"812 Rsun = {R:.4g} cm", 'App. B (42.61 mas K-band / 1.030)'),
        ('Teff', '3650 K', 'App. B (their MARCS model)'),
        ('L = 4 pi R^2 sigma Teff^4', f"{st['L']/P.LSUN:.4g} Lsun = {st['L']:.3e} erg/s",
         'derived (paper quotes no L)'),
        ('M', '17.5 Msun', 'App. B'),
        ('g(R*) = GM/R*^2', f"{st['g']:.3f} cm/s^2", 'paper: 0.73'),
        ('k_Edd = 4 pi c G M / L', f"{st['kE']:.3f} cm^2/g", ''),
        ('v_esc(R*)', f"{st['vesc']/1e5:.1f} km/s", ''),
        ('rho = n_H x 2.27e-24 g', 'mass per H 1.357 m_H (X = 0.737)',
         'brief; He-only A(He)=0.083 would give 1.332 m_H'),
        ('mu for P_gas', 'FastChem CE (lowP table), iterated on (T, rho)', 'wind.mu_of'),
        ('SEM wind term Mdot', f"{W.mdot26(st):.3g} Msun/yr (v_inf 9 km/s)",
         'paper 3.8e-6 (with 1.332 m_H: 3.79e-6)'),
        ('T_e(R)', 'Fig. 5 digitized (out/dent26_fig5_T.txt, +-50 K), 1.13-5 R*',
         'dip set to column minimum 1776 K at 1.19-1.20 R* (text: min ~1700 K)'),
        ('T beyond 5 R*', 'harper (baseline): Harper+2001 shape scaled to 1412 K; plaw: '
         'T ~ x^-0.67; const: 1412 K', 'T(10/20/30 R*) = ' + '; '.join(
             f"{t}: " + '/'.join(f"{v:.0f}" for v in P.T26(np.array([10., 20., 30.]), t))
             for t in TEXTS)),
        ('n_e/n_H', '2.5e-5 (1.144 R*) -> 5e-4 (2 R*), const beyond', 'not used by the force'
         ' calculation (LTE ck tables)')):
    row(q, v, s)

# scale height / turbulence in the hydrostatic chromosphere
Hobs = 0.06 * R
xm = P.XMIN
gm = g0 / xm**2
im = np.argmin(abs(w['x'] - 1.15))
a2 = P.KB * w['T'][im] / (w['mu'][im] * P.MH)
L.append('\n### Chromospheric scale height (H = 0.06 R* at R_min = 1.144 R*)\n')
L.append('Eq. (2): H = v_turb^2/(2 g*) (R_min/R*)^2 with v_turb = 19.5 km/s (most probable '
         'speed; P_turb = rho v_mp^2/2 = rho v_t^2 with v_t = v_mp/sqrt2 the 1-D rms).\n')
row('case', 'H / R*')
hdr(2)
for vt in (14e5, 19.5e5):
    for gs, lab in ((0.73, 'g = 0.73 (unreduced)'), (0.32, 'g* = 0.32 (reduced)')):
        row(f'Eq. 2, v_turb = {vt/1e5:g}, {lab}', f(vt**2 / (2 * gs) * xm**2 / R, 3))
L.append(f'\nEq. (2) reproduces 0.06 R* only with the UNREDUCED g = 0.73 (0.0599 R*); with '
         f'g* = 0.32 it gives 0.137 R*. The reduction 0.32/0.73 corresponds to an effective '
         f'Gamma = 1 - 0.32/0.73 = {1-0.32/0.73:.2f}; it is applied to the MARCS photosphere '
         f'only, not to the SEM chromosphere (whose H implies Gamma ~ 0 at v_mp = 19.5).\n')
row('required isotropic 1-D v_t for H_rho = 0.06 R* at 1.15 R* (thermal a^2 = kT/mu m_H '
    'included)', 'Gamma = 0', 'Gamma = 0.56')
hdr(3)
row('isothermal: v_t^2 = g(1-Gamma) H - a^2 [km/s]',
    *[f(np.sqrt(max(gm * (1 - G) * Hobs - a2, 0)) / 1e5, 3) for G in (0., 1 - 0.32 / 0.73)])
row('same as v_mp = sqrt2 v_t [km/s]',
    *[f(np.sqrt(2 * max(gm * (1 - G) * Hobs - a2, 0)) / 1e5, 3) for G in (0., 1 - 0.32 / 0.73)])
vt0 = W.vt_equiv(w, 0.)
L.append(f'\nThermal speed at 1.15 R*: a = {np.sqrt(a2)/1e5:.2f} km/s (T {w["T"][im]:.0f} K, '
         f'mu {w["mu"][im]:.2f}). Full calculation (T gradient, v dv/dr, local H_rho) '
         f'v_t,req(Gamma=0) at 1.15/1.2/1.3/1.4 R* = '
         + '/'.join(f(v / 1e5, 3) for v in at(x, vt0, [1.15, 1.2, 1.3, 1.4])) + ' km/s.\n')

# ------------------------------------------------------------------ Gamma_req
L.append('## Gamma_req (steady wind through the SEM: v(r) = 9 ((0.998 - x^-0.45)/0.998)^1.5 '
         'km/s, rho = total SEM rho, P = P_gas + rho v_t^2)\n')
L.append('Gamma_req = 1 + r^2/(GM) [v dv/dr + (1/rho) dP_gas/dr + v_t^2 dln rho/dr]. Inside '
         '~1.5 R* the exponential chromosphere carries most of the mass (f_chrom row) and the '
         'continuity velocity Mdot/(4 pi r^2 rho_total) falls below the SEM wind-term v: the '
         'wind velocity is ill-defined there, Gamma_req is the hydrostatic+turbulent balance '
         '(v dv/dr / g < 2e-4).\n')
row('T ext', 'v_t [km/s]', *[f'{v:g}' for v in XS], 'Gamma_req = 1 at r/R*')
hdr(3 + len(XS))
CROSS = {}
for t in TEXTS:
    for vt in VTS:
        y = WS[t]['Greq'][vt]
        cr = W.crossings(x, y)
        CROSS[(t, vt)] = cr
        row(t, f'{vt/1e5:g}', *[f(v, 3) for v in at(x, y)],
            ', '.join(f'{c:.3g}' for c in cr) or 'none')
L.append('')
row('quantity (T ext harper)', *[f'{v:g}' for v in XS])
hdr(1 + len(XS))
for lab, y, n in (('r [cm]', w['r'], 3), ('rho [g/cm^3]', w['rho'], 3), ('T [K]', w['T'], 4),
                  ('mu', w['mu'], 3), ('f_chrom = n_exp/n_H', w['f_chrom'], 2),
                  ('v SEM [km/s]', w['v'] / 1e5, 3), ('v continuity [km/s]', w['vc'] / 1e5, 3),
                  ('v dv/dr / g', w['term_acc'], 2), ('(1/rho)dP_gas/dr / g', w['term_gas'], 2),
                  ('H_rho / R*', w['Hrho'] / R, 3), ('v_t,req (Gamma_avail = 0) [km/s]',
                                                    vt0 / 1e5, 3)):
    row(lab, *[f(v, n) for v in at(x, y)])
md_c = 4 * np.pi * w['r']**2 * w['rho'] * w['v'] * P.YR / P.MSUN
L.append(f'\nContinuity check: 4 pi r^2 rho_total v_SEM = {md_c[np.argmin(abs(x-1.5))]:.3g} '
         f'(1.5 R*), {md_c[np.argmin(abs(x-2))]:.3g} (2), {md_c[np.argmin(abs(x-5))]:.3g} (5), '
         f'{md_c[-1]:.3g} (30 R*) Msun/yr vs the wind-term {W.mdot26(st):.3g}: they agree to '
         f'<0.1% beyond 2 R*; at 1.15 R* the total is {md_c[0]/W.mdot26(st):.1f}x (chromosphere).'
         f' v reaches only {w["v"][-1]/1e5:.2f} km/s at 30 R* (9 km/s asymptotically).\n')

# ------------------------------------------------------------------ gas force
GAS = {}
for tc, t in (('obs', 'harper'), ('obs', 'plaw'), ('obs', 'const'), ('re', 'harper'),
              ('300', 'harper'), ('600', 'harper'), ('900', 'harper'), ('1500', 'harper'),
              ('2000', 'harper')):
    GAS[(tc, t)] = np.load(O + f'gas_A_C17.5_dent26_{tc}_{t}.npz')
L.append('## Gamma_gas (k_F / k_Edd; LTE; DACE table A; M = 17.5)\n')
L.append('static = ck column (formal solution; grey-Lucy HSE photosphere below 1.027 R*, SEM '
         'density above, RE temperature below x_min = 1.144 R*, T fixed above: obs = Fig. 5 T, '
         're = LTE radiative equilibrium everywhere); thin = unattenuated stellar flux; Sob w = '
         'Sobolev with dv/dr of the SEM wind; Sob tN = Sobolev with |dv/dr| = v_t/H_rho '
         '(v_t = N km/s); A = unattenuated, B = window-attenuated flux. Fixed-T brackets '
         '(300-2000 K) apply for r >= 1.2 R*.\n')
row('T', 'T ext', 'estimate', *[f'{v:g}' for v in XS[:-1]], 'max 1.15-30 R* (at r/R*)')
hdr(4 + len(XS) - 1)
EST = [('thin', 'Gthin'), ('Sob w A', 'Gsob_w_A'), ('Sob w B', 'Gsob_w_B'),
       ('Sob t5 A', 'Gsob_t5_A'), ('Sob t10 A', 'Gsob_t10_A'), ('Sob t14 A', 'Gsob_t14_A'),
       ('Sob t19.5 A', 'Gsob_t19.5_A'), ('Sob t19.5 B', 'Gsob_t19.5_B')]
for (tc, t), d in GAS.items():
    ests = [('static', d['cr'], d['cGF'])] if tc in ('obs', 're') else []
    if tc != 're':
        ests += [(n, d['sx'], d[k]) for n, k in EST]
    for name, xx, gg in ests:
        m = (xx >= 1.15) & (xx <= 30)
        i = np.argmax(np.where(m, gg, -1))
        row(tc, t if tc == 'obs' else '-', name, *[f(v, 2) for v in at(xx, gg, XS[:-1])],
            f'{gg[i]:.3f} ({xx[i]:.1f})')
L.append('\nF/F* at the column top: ' + ', '.join(f'{k[0]}/{k[1]}: {v["cFratio"][-1]:.2f}'
                                                  for k, v in GAS.items()) +
         '. Static columns with T fixed ABOVE radiative equilibrium (the const extension, '
         '1412 K to 30 R*) are self-luminous in their bands (static 0.43 at 20 R* vs thin '
         '0.11): an energy-violating artefact, excluded from the bound below, as in the 2024 '
         'run.\n')


def gi(d, key, xg='sx'):
    return np.interp(np.log(x), np.log(d[xg] if xg else d['cr']), d[key])


# most generous gas force: max over static RE, static obs (harper / plaw), thin obs T,
# Sobolev t19.5 A, thin 1500 K (cool-phase bracket)
comp = {'static RE': np.interp(np.log(x), np.log(GAS[('re', 'harper')]['cr']),
                               GAS[('re', 'harper')]['cGF']),
        'static obs': np.maximum(*[np.interp(np.log(x), np.log(GAS[('obs', t)]['cr']),
                                             GAS[('obs', t)]['cGF']) for t in ('harper', 'plaw')]),
        'thin obs T': np.maximum.reduce([gi(GAS[('obs', t)], 'Gthin') for t in TEXTS]),
        'Sob t19.5 A': np.maximum.reduce([gi(GAS[('obs', t)], 'Gsob_t19.5_A') for t in TEXTS]),
        'thin 1500 K': gi(GAS[('1500', 'harper')], 'Gthin')}
Gmax = np.maximum.reduce(list(comp.values()))
Gbest = np.maximum(gi(GAS[('obs', 'harper')], 'Gsob_t19.5_A'),
                   np.interp(np.log(x), np.log(GAS[('obs', 'harper')]['cr']),
                             GAS[('obs', 'harper')]['cGF']))

# ------------------------------------------------------------------ dust
RS = {t: pickle.load(open(O + f'dust_steady_{t}.pkl', 'rb')) for t in TEXTS}
RB = pickle.load(open(O + 'dust_ball_harper.pkl', 'rb'))
L.append('## Dust (a): steady SEM-wind trajectory from 1.15 R* (kinematic; T_gas = Fig. 5 + '
         'extension; forsterite; Td radiative equilibrium)\n')
row('T ext', 'optics', 'alpha', 'n_seed/n_H', 'seeds grow at r/R* (cm)',
    'a = 0.1 um at r/R*', 't(1.15 R* -> r_0.1) [yr]', 'a at 2/3/5/10/40 R* [um]',
    'Gamma_dust at 2/3/5/10/20/30 R*', 'f_cond(40 R*)')
hdr(10)
SEED = {}
for t in TEXTS:
    for r in RS[t]:
        g_ = np.where(r['a'] > 1.05e-7)[0]
        i1 = np.where(r['a'] > 1e-5)[0]
        rs = r['r'][g_[0]] if len(g_) else np.nan
        SEED[(t, r['var'], r['alpha'], r['xseed'])] = (rs, r['r'][i1[0]] if len(i1) else np.nan)
        row(t, r['var'], f"{r['alpha']:g}", f"{r['xseed']:.0e}",
            f'{rs:.3g} ({rs*R:.3g})' if np.isfinite(rs) else '-',
            f(r['r'][i1[0]]) if len(i1) else '-',
            f(r['t'][i1[0]] / gl.YR, 2) if len(i1) else '-',
            '/'.join(f(v * 1e4, 2) for v in at(r['r'], r['a'], [2, 3, 5, 10, 40])),
            '/'.join(f(v, 2) for v in at(r['r'], r['gd'], [2, 3, 5, 10, 20, 30])),
            f(r['f'][-1], 2))
L.append('\nThe T extension beyond 5 R* changes nothing above (seeding is set by the grain '
         'temperature, which is radiative equilibrium with the stellar field; T_gas only enters '
         'the sticking speed and the sqrt(Td/Tg) factor).\n')

L.append('## Dust (b): ballistic parcels launched from R* at v0 (20-100 km/s, step 5), parcel '
         'density = SEM rho(r), T_gas = Fig. 5 (harper ext), launch flux Mdot(>v0) = 4 pi R^2 '
         'rho_ph v_con exp(-v0/v_con) with the FT eq. 7 fit to the SEM over 1.2-5 R*; outcome '
         'to 40 R*\n')
gs = gl.Star('bet26')
ft = P.ft_fit(st, 'dent26')
L.append(f"FT fit: v_esc/v_con = {ft['q']:.2f}, v_con = {ft['vcon']/1e5:.2f} km/s, rho_ph = "
         f"{ft['rho_ph']:.3g}, Mdot(>v_esc) = {ft['mdot_esc']:.3g} Msun/yr, rms "
         f"{ft['rms_dex']:.2f} dex (2024: q 10.35, v_con 7.96, 3.97e-6).\n")
row('optics', 'alpha', 'n_seed/n_H', 'Gamma_gas', 'dust-free escape v0 / Mdot',
    'min v0 for a wind with dust', 'Mdot(dust) [Msun/yr]', 'v_inf mass-wt [km/s]',
    'seeds grow at', 'a=0.1um at', 'max Gamma_d', 'timeouts')
hdr(12)
for var in ('pure', 'Fe3e-4', 'Fe1e-3'):
    for al in (0.1, 1.0):
        for xs in (1e-15, 1e-13):
            for gg in (0., 0.1):
                rr = sorted([r for r in RB if r['var'] == var and r['alpha'] == al and
                             r['xseed'] == xs and r['Ggas'] == gg], key=lambda r: r['v0'])
                v0s = np.array([r['v0'] for r in rr])
                wind_ = np.array([r['outcome'] == 'wind' for r in rr])
                nto = sum(r['outcome'] == 'timeout' for r in rr)
                vfree = gs.vesc * np.sqrt(1 - gg) * np.sqrt(1 - 1 / 40.) / 1e5
                mfree = gs.mdot_gt(vfree * 1e5) * gl.YR / gl.Msun
                dust_w = wind_ & (v0s < vfree)
                if dust_w.any():
                    k = np.where(dust_w)[0][0]
                    md = gs.mdot_gt(v0s[k] * 1e5) * gl.YR / gl.Msun
                    edges = np.concatenate([v0s - 2.5, [v0s[-1] + 2.5]]) * 1e5
                    dm = gs.mdot_gt(edges[:-1]) - gs.mdot_gt(edges[1:])
                    vi = np.array([r.get('vinf', np.nan) for r in rr])
                    ok = wind_ & np.isfinite(vi)
                    vinf = (dm[ok] * vi[ok]).sum() / dm[ok].sum() / 1e5
                    rs, r01 = rr[k]['r_seed'], rr[k]['r_01']
                    gdm = max(r.get('gdmax', 0) for r in rr if r['outcome'] == 'wind')
                else:
                    k, md, vinf, rs, r01 = None, 0., np.nan, np.nan, np.nan
                    gdm = max(r.get('gdmax', 0) or 0 for r in rr)
                row(var, f'{al:g}', f'{xs:.0e}', f'{gg:g}', f'{vfree:.1f} / {mfree:.2e}',
                    f(v0s[k]) if k is not None else 'none', f'{md:.2e}' if md else '0',
                    f(vinf, 3), f(rs), f(r01), f(gdm, 2), nto)

# ------------------------------------------------------------------ support budget
L.append('\n## Required non-radiative support (T ext harper unless stated)\n')
L.append('Gamma_gas,max = max{static RE, static obs (harper/plaw), thin at observed T (all '
         'ext), Sobolev v_t=19.5/H_rho flux A (all ext), thin at 1500 K} -- the most generous '
         'gas bound; v_t,eq = sqrt[(Gamma_req(v_t=0) - Gamma_avail) g H_rho]; P_w = rho '
         'v_t,eq^2.\n')
row('r/R*', *[f'{v:g}' for v in XS])
hdr(1 + len(XS))
row('Gamma_req (v_t = 0)', *[f(v, 3) for v in at(x, w['Greq'][0.])])
for k, v in comp.items():
    row(f'Gamma_gas {k}', *[f(q, 2) for q in at(x, v)])
row('Gamma_gas,max', *[f(q, 2) for q in at(x, Gmax)])
dpk = [q for q in RS['harper'] if q['var'] == 'pure' and q['alpha'] == 1. and q['xseed'] == 1e-13][0]
dfe = [q for q in RS['harper'] if q['var'] == 'Fe1e-3' and q['alpha'] == 0.1 and
       q['xseed'] == 1e-15][0]
dmid = [q for q in RS['harper'] if q['var'] == 'Fe3e-4' and q['alpha'] == 1. and
        q['xseed'] == 1e-15][0]
row('Gamma_dust steady, pure a=1 1e-13', *[f(q, 2) for q in at(dpk['r'], dpk['gd'])])
row('Gamma_dust steady, Fe3e-4 a=1 1e-15', *[f(q, 2) for q in at(dmid['r'], dmid['gd'])])
row('Gamma_dust steady, Fe1e-3 a=0.1 1e-15', *[f(q, 2) for q in at(dfe['r'], dfe['gd'])])
row('deficit Gamma_req - Gamma_gas,max', *[f(q, 2) for q in at(x, w['Greq'][0.] - Gmax)])
VTEQ = {t: W.vt_equiv(WS[t], Gmax) for t in TEXTS}
for t in TEXTS:
    row(f'v_t,eq [km/s] (gas only, ext {t})', *[f(q / 1e5, 3) for q in at(x, VTEQ[t])])
row('v_t,eq [km/s] (Gamma_avail = 0)', *[f(q / 1e5, 3) for q in at(x, vt0)])
vteq = VTEQ['harper']
Pw = w['rho'] * vteq**2
row('P_w [dyn/cm^2]', *[f(q, 2) for q in at(x, Pw)])
row('P_w / P_gas', *[f(q, 2) for q in at(x, Pw / w['P'])])
row('v_t,eq / v_mp-equivalent sqrt2 v_t,eq [km/s]', *[f(np.sqrt(2) * q / 1e5, 3)
                                                      for q in at(x, vteq)])

# ------------------------------------------------------------------ comparison with 2024
o = np.load(O + 'old24.npz')
Ro = float(o['R'])
L.append('\n## Side by side: Dent 2024 run (set A: R* 1014 Rsun, 222 pc, M 18, Mdot 2e-6, v_inf '
         'cap 12 km/s) vs Dent 2026 SEM (R* 812 Rsun, 172 pc, M 17.5, SEM v(r), Mdot 3.86e-6)\n')
L.append('Gamma_gas max for 2024 = max{static RE, Sobolev t10 A, thin 1500 K} (its original '
         'definition); for 2026 = the wider bound above. Dust seed radius = first growth, '
         'pure / alpha 1 / 1e-13 (range over all 12 optics/alpha/seed cases in brackets).\n')
sd24 = o['seeds'][:, 0]
sd26 = np.array([v[0] for k, v in SEED.items() if k[0] == 'harper'])


def blk(title, xs_new, xs_old, unit):
    L.append(f'\n### {title}\n')
    row(unit + ' (2026 / 2024)', *[f'{a:.3g} / {b:.3g}' for a, b in zip(xs_new, xs_old)])
    hdr(1 + len(xs_new))
    pairs = (('rho [g/cm^3]', w['rho'], o['rho'], 3), ('T [K]', w['T'], o['T'], 4),
             ('Gamma_req v_t=0', w['Greq'][0.], o['G0'], 3),
             ('Gamma_req v_t=10', w['Greq'][10e5], o['G10'], 3),
             ('v_t,req (Gamma_avail=0) [km/s]', vt0 / 1e5, o['vt0'] / 1e5, 3),
             ('v_t,eq after Gamma_gas,max [km/s]', vteq / 1e5, o['vteq'] / 1e5, 3),
             ('Gamma_gas,max', Gmax, o['Gmax'], 2), ('v [km/s]', w['v'] / 1e5, o['v'] / 1e5, 2))
    for lab, yn, yo, n in pairs:
        row(lab, *[f'{f(a, n)} / {f(b, n)}' for a, b in
                   zip(at(x, yn, xs_new), at(o['x'], yo, xs_old))])


blk('same r/R*', XS, XS, 'r/R*')
rc = np.array([1.2, 1.5, 2, 3, 5, 10, 20, 30]) * Ro       # 2024 radii in cm
blk('same r [cm] (2024 r/R* grid; 2026 x = r/(812 Rsun))', list(rc / R), list(rc / Ro),
    'x_2026 / x_2024')
L.append('\nr [cm] of the second block: ' + ', '.join(f'{v:.3g}' for v in rc) + '.')
L.append(f'\nDust seed radius (pure/a1/1e-13): 2026 {SEED[("harper", "pure", 1.0, 1e-13)][0]:.3g}'
         f' R* = {SEED[("harper", "pure", 1.0, 1e-13)][0]*R:.3g} cm; 2024 '
         f'{sd24[3]:.3g} R* = {sd24[3]*Ro:.3g} cm. All cases: 2026 {np.nanmin(sd26):.3g}-'
         f'{np.nanmax(sd26):.3g} R* ({np.nanmin(sd26)*R:.3g}-{np.nanmax(sd26)*R:.3g} cm); 2024 '
         f'{np.nanmin(sd24):.3g}-{np.nanmax(sd24):.3g} R* ({np.nanmin(sd24)*Ro:.3g}-'
         f'{np.nanmax(sd24)*Ro:.3g} cm). Labels 2024: ' + ', '.join(
             f'{a}={b:.3g}' for a, b in zip(o['seedlab'], sd24)) + '.\n')
L.append('Gamma_req crossings of 1 (2026, harper): ' + '; '.join(
    f'v_t {vt/1e5:g}: ' + (', '.join(f'{c:.3g}' for c in CROSS[('harper', vt)]) or 'none')
    for vt in VTS) + '. 2024 (A/18 dent 2e-6): v_t 0 at 19.3 R*, none for v_t >= 5.\n')

open(O + 'tables.md', 'w').write('\n'.join(L) + '\n')

# ------------------------------------------------------------------ plots
xx = np.exp(np.linspace(np.log(1.0), np.log(30.), 500))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
ax[0].plot(xx * R, P.rho(xx, 'dent26'), color=C[0], label='Dent 2026 SEM (total)')
ax[0].plot(xx[xx >= P.XMIN] * R, P.nH26_wind(xx[xx >= P.XMIN]) * P.MPH, color=C[0], lw=1.2, ls='--',
           label='  wind term only')
ax[0].plot(o['x'] * Ro, o['rho'], color=C[1], label='Dent 2024 (set A, 1014 Rsun)')
ax[0].axvline(P.XMIN * R, color=INK2, lw=0.8, ls=':')
ax[0].set(xscale='log', yscale='log', xlabel='r [cm]', ylabel='rho [g/cm^3]',
          xlim=(1.0 * R, 30 * R), ylim=(5e-18, 1e-10))
topaxis(ax[0], R)
ax[0].legend(fontsize=8)
for t, c in zip(TEXTS, (C[0], C[2], C[3])):
    ax[1].plot(xx * R, P.T26(xx, t), color=c, lw=2 if t == 'harper' else 1.2,
               label=f'2026 Fig. 5, ext {t}')
ax[1].plot(o['x'] * Ro, o['T'], color=C[1], label='2024 run T')
ax[1].set(xscale='log', yscale='log', xlabel='r [cm]', ylabel='T [K]', xlim=(1.0 * R, 30 * R),
          ylim=(150, 5000))
topaxis(ax[1], R)
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(PL + 'rho_T.png', dpi=150)
plt.close(fig)

fig, ax = plt.subplots(figsize=(8.5, 5))
r = w['r']
for vt, ls in zip(VTS, ('-', '--', '-.', ':', (0, (1, 3)))):
    ax.plot(r, np.maximum(w['Greq'][vt], 1e-4), color=INK, lw=1.6 if vt == 0 else 1.1, ls=ls,
            label=f'Gamma_req, v_t = {vt/1e5:g}')
ax.plot(o['r'], o['G0'], color=C[1], lw=1.2, label='Gamma_req v_t=0, Dent 2024 run')
ax.plot(r, Gmax, color=C[0], label='Gamma_gas,max (generous bound)')
ax.plot(r, comp['static RE'], color=C[6], lw=1.2, label='gas: static RE column')
ax.plot(r, gi(GAS[('obs', 'harper')], 'Gsob_t19.5_A'), color=C[3], lw=1.2,
        label='gas: Sobolev 19.5 km/s / H_rho')
ax.plot(dpk['r'] * R, np.maximum(dpk['gd'], 1e-4), color=C[5], label='dust: pure, a 1, 1e-13')
ax.plot(dmid['r'] * R, np.maximum(dmid['gd'], 1e-4), color=C[4],
        label='dust: Fe3e-4, a 1, 1e-15')
ax.plot(dfe['r'] * R, np.maximum(dfe['gd'], 1e-4), color=C[7], label='dust: Fe1e-3, a 0.1, 1e-15')
ax.axhline(1, color=INK2, lw=0.8)
ax.set(xscale='log', yscale='log', xlabel='r [cm]', ylabel='Gamma (force / gravity)',
       ylim=(1e-3, 10), xlim=(1.15 * R, 30 * R))
topaxis(ax, R)
ax.legend(fontsize=7.2, loc='lower right', ncol=2)
fig.suptitle('Betelgeuse, Dent 2026 SEM: required vs available Gamma (T ext harper)',
             fontsize=10, color=INK)
fig.tight_layout()
fig.savefig(PL + 'gamma.png', dpi=150)
plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 4.5))
for t, c in zip(TEXTS, (C[0], C[2], C[3])):
    ax.plot(r, VTEQ[t] / 1e5, color=c, lw=2 if t == 'harper' else 1.2,
            label=f'v_t,eq after Gamma_gas,max (ext {t})')
ax.plot(r, vt0 / 1e5, color=INK, lw=1.2, ls='--', label='v_t,req with no radiative force')
ax.plot(o['r'], o['vteq'] / 1e5, color=C[1], lw=1.2, label='Dent 2024 run, v_t,eq')
ax.axhspan(14, 19.5, color=INK2, alpha=0.12, lw=0)
ax.text(1.5 * R, 19.9, 'paper v_turb 14-19.5 km/s (most probable)', fontsize=8, color=INK2)
ax.axhspan(14 / np.sqrt(2), 19.5 / np.sqrt(2), color=C[0], alpha=0.08, lw=0)
ax.text(1.5 * R, 9.0, '  same as 1-D rms (/sqrt2)', fontsize=8, color=INK2)
ax.set(xscale='log', xlabel='r [cm]', ylabel='isotropic 1-D v_t [km/s]', xlim=(1.15 * R, 30 * R),
       ylim=(0, 30))
topaxis(ax, R)
ax.legend(fontsize=8, loc='upper right')
fig.tight_layout()
fig.savefig(PL + 'vt_req.png', dpi=150)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7.5, 4.5))
for i, var in enumerate(('pure', 'Fe3e-4', 'Fe1e-3')):
    for al, xs, ls in ((1., 1e-13, '-'), (0.1, 1e-15, '--')):
        q = [q for q in RS['harper'] if q['var'] == var and q['alpha'] == al and
             q['xseed'] == xs][0]
        ax.plot(q['r'] * R, q['a'] * 1e4, ls=ls, color=C[i],
                label=f'{var}, alpha {al:g}, seeds {xs:.0e}')
ax.axhline(0.1, color=INK2, lw=1)
for x0, x1, lab in ((1.4, 1.6, 'Haubois+19 ~1.5 R*'), (12., 14., 'MATISSE 13 R*'),
                    (24., 40., 'VISIR shell')):
    ax.axvspan(x0 * R, x1 * R, color=INK2, alpha=0.10, lw=0)
    ax.text(x0 * R * 1.02, 0.25, lab, fontsize=7.5, color=INK2, rotation=90)
ax.set(xscale='log', yscale='log', xlabel='r [cm]', ylabel='grain radius a [um]',
       ylim=(8e-4, 1), xlim=(1.15 * R, 40 * R))
topaxis(ax, R)
ax.legend(fontsize=7.5, loc='lower right')
fig.suptitle('Steady SEM wind (v -> 9 km/s), Dent 2026 rho and T', fontsize=10, color=INK)
fig.tight_layout()
fig.savefig(PL + 'grain_a.png', dpi=150)
plt.close(fig)

# ------------------------------------------------------------------ web.json


def rl(a, n=4):
    return [None if not np.isfinite(v) else float(f'{v:.{n}g}') for v in np.asarray(a, float)]


sub = slice(None, None, 3)
web = dict(
    meta=dict(model='Dent et al. 2026 SEM (arXiv:2608.19339v2 App. B)', R_star_cm=R,
              R_star_Rsun=812, d_pc=172, M_Msun=17.5, Teff_K=3650,
              L_Lsun=round(st['L'] / P.LSUN), mdot_Msun_yr=round(W.mdot26(st), 9),
              v_inf_kms=9, T_ext='harper (Harper+2001 shape beyond 5 R*)',
              units=dict(r_cm='cm', rho='g/cm^3', T='K', v='km/s', vt='km/s', a='um')),
    r_cm=rl(w['r'][sub]), x=rl(x[sub]), rho=rl(w['rho'][sub]), T=rl(w['T'][sub]),
    T_plaw=rl(WS['plaw']['T'][sub]), T_const=rl(WS['const']['T'][sub]),
    v_kms=rl(w['v'][sub] / 1e5), v_continuity_kms=rl(w['vc'][sub] / 1e5),
    Gamma_req={f'vt{vt/1e5:g}': rl(w['Greq'][vt][sub]) for vt in VTS},
    vt_req_noforce_kms=rl(vt0[sub] / 1e5), vt_req_kms=rl(vteq[sub] / 1e5),
    Gamma_gas_max=rl(Gmax[sub]), Gamma_gas_best=rl(Gbest[sub]),
    old2024=dict(r_cm=rl(o['r'][::3]), rho=rl(o['rho'][::3]), T=rl(o['T'][::3]),
                 Gamma_req_vt0=rl(o['G0'][::3]), vt_req_kms=rl(o['vteq'][::3])),
    dust={})
for lab, q in (('pure_a1_1e-13', dpk), ('Fe3e-4_a1_1e-15', dmid), ('Fe1e-3_a0.1_1e-15', dfe)):
    web['dust'][lab] = dict(r_cm=rl(q['r'][::2] * R), a_um=rl(q['a'][::2] * 1e4),
                            Gamma_dust=rl(q['gd'][::2]))
json.dump(web, open(O + 'web.json', 'w'), separators=(',', ':'))
print('\n'.join(L))
