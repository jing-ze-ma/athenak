"""Tables (out/tables.md) and plots (plots/*.png) for the Betelgeuse observed-density wind
analysis.  usage: ../grains/venv/bin/python analyse.py"""
import os
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import prof as P  # noqa: E402
import wind as W  # noqa: E402
import grainlib_bet as gl  # noqa: E402

B = '/orion/ptmp/jinma/rsg_wind_1008/betelgeuse/'
O, PL = B + 'out/', B + 'plots/'
os.makedirs(PL, exist_ok=True)
XS = [1.2, 1.5, 2., 3., 5., 10., 20., 30.]
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
INK, INK2 = '#0b0b0b', '#52514e'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.grid': True,
                     'grid.color': '#e4e3df', 'grid.linewidth': 0.6, 'lines.linewidth': 2,
                     'legend.frameon': False, 'figure.facecolor': '#fcfcfb',
                     'axes.facecolor': '#fcfcfb'})
L = []


def at(x, y, xs=XS):
    return [float(np.interp(np.log(v), np.log(x), y)) if x[0] <= v <= x[-1] else np.nan
            for v in xs]


def f(v, n=3):
    return '-' if v is None or not np.isfinite(v) else f'{v:.{n}g}'


def row(*c):
    L.append('| ' + ' | '.join(str(x) for x in c) + ' |')


def topaxis(ax, R):
    t = ax.secondary_xaxis('top', functions=(lambda r: r/R, lambda x: x*R))
    t.set_xlabel('r / R*  (set A, R* = 1014 Rsun)')
    t.set_xticks([1, 1.5, 2, 3, 5, 10, 20, 30])
    t.set_xticklabels(['1', '1.5', '2', '3', '5', '10', '20', '30'])


# ------------------------------------------------------------------ parameters
L.append('# Betelgeuse: observed extended-atmosphere density -> wind force budget\n')
L.append('## Parameters\n')
row('set', 'd [pc]', 'R [Rsun]', 'Teff [K]', 'L [Lsun]', 'M [Msun]', 'k_Edd [cm^2/g]',
    'v_esc [km/s]', 'n_H scaling', 'FT eq.7 fit to Dent 1.2-5 R*: v_esc/v_con, v_con, '
    'rho_ph, Mdot(>v_esc), rms')
L.insert(len(L), '|' + '---|'*10)
for s in 'AB':
    for M in P.SETS[s]['M']:
        st = P.stellar(s, M)
        ft = P.ft_fit(st)
        row(s, f(st['d']), f(P.SETS[s]['R'], 4), f(st['Teff'], 4), f(st['L']/P.LSUN),
            f(M), f(st['kE']), f(st['vesc']/1e5), 'x1' if s == 'A' else 'x(222/168)^0.5=1.15',
            f"{ft['q']:.2f}, {ft['vcon']/1e5:.2f} km/s, {ft['rho_ph']:.2e}, "
            f"{ft['mdot_esc']:.2e} Msun/yr, {ft['rms_dex']:.2f} dex")
L.append('\nSources: set A Dent et al. 2024 (R* 1014 Rsun at 222 pc, 42.5 mas); set B Joyce et '
         'al. 2020 (168 pc, 764 Rsun = same angular radius, so r/R* is unchanged and n_H ~ '
         'd^-1/2); n_H(r): Dent 2024 Fig. 3 digitized (1.08-5 R*), r^-2.85 beyond 5 R* '
         f'(spline slope at 5 R*: {P._slope5:.2f}); Harper+2001 Table 7 rescaled (x = 1.32 x_H, '
         'n_H x (131/d)^0.5), used at >= 1.5 R*; T: Dent Fig. 3 read off at 1.2-5 R*, Harper '
         'shape beyond (scaled x%.3f to match at 5 R*); rho = 2.27e-24 n_H; mu from FastChem '
         '(1.25 atomic-H warm gas, 2.3 once H2 forms beyond ~6-7 R*).\n' % P._fT)

# ------------------------------------------------------------------ Gamma_req
L.append('## Gamma_req (steady wind through the observed rho; v_inf cap 12 km/s)\n')
row('set/M', 'profile', 'Mdot', 'v_t [km/s]', *[f'{x:g} R*' for x in XS],
    'Gamma_req = 1 at r/R*', 'v at 2/5/10/20/30 R* [km/s]')
L.append('|' + '---|'*(5 + len(XS)))
WIND = {}
for s, M in (('A', 18.), ('A', 20.), ('B', 16.5), ('B', 19.)):
    st = P.stellar(s, M)
    for prof in ('dent', 'harper'):
        if prof == 'harper' and M != 18.:
            continue
        for md in (1e-6, 2e-6, 4e-6):
            w = W.steady(st, prof, md, x=W.xgrid(1.2 if prof == 'dent' else 1.5, 30., 600))
            WIND[(s, M, prof, md)] = w
            for vt in (0., 5e5, 10e5):
                g = w['Greq'][vt]
                cr = W.crossings(w['x'], g)
                row(f'{s}/{M:g}', prof, f'{md:.0e}', f'{vt/1e5:g}', *[f(v, 3) for v in
                                                                      at(w['x'], g)],
                    ', '.join(f'{c:.1f}' for c in cr) or 'none',
                    '/'.join(f(v, 2) for v in at(w['x'], w['v']/1e5, [2, 5, 10, 20, 30]))
                    if vt == 0 else '')
w = WIND[('A', 18., 'dent', 2e-6)]
L.append('\nTerms (A/18, Dent, 2e-6): (1/rho)dP_gas/dr / g = ' +
         ', '.join(f'{x:g}: {v:.3f}' for x, v in zip(XS, at(w['x'], w['term_gas']))) +
         '; v dv/dr / g = ' + ', '.join(f'{x:g}: {v:.3f}' for x, v in
                                         zip(XS, at(w['x'], w['term_acc']))) +
         '; H_rho/R* = ' + ', '.join(f'{x:g}: {v:.3f}' for x, v in
                                    zip(XS, at(w['x'], w['Hrho']/w['r']*w['x']))) +
         f"; v reaches 12 km/s at {', '.join(f'{k[3]:.0e}: {f(v["x_vinf"])}' for k, v in WIND.items() if k[:3] == ('A', 18., 'dent'))} R*.\n")


# ------------------------------------------------------------------ gas force
def gload(tab, sm, prof, tc):
    return np.load(O + f'gas_{tab}_{sm}_{prof}_{tc}.npz')


L.append('## Gamma_gas (k_F / k_Edd; LTE; M = 18 for set A, 16.5 for B -- scale by 18/20 '
         'resp. 16.5/19 for the upper masses)\n')
L.append('static = ck column (formal solution, T fixed: obs = observed T, re = LTE radiative '
         'equilibrium); thin = unattenuated stellar flux; Sob w = Sobolev with dv/dr of the '
         'steady wind (Mdot 2e-6), flux A; Sob t5/t10 = Sobolev with |dv/dr| = v_t/H_rho, flux A;'
         ' B = window-attenuated flux.\n')
row('table', 'set', 'profile', 'T', 'estimate', *[f'{x:g}' for x in XS[:-1]],
    'max 1.2-30 R* (at r/R*)')
L.append('|' + '---|'*(6 + len(XS) - 1))
GAS = {}
for tab in ('A', 'old'):
    for sm, prof, tcs in (('A18', 'dent', ['obs', 're', '300', '600', '900', '1500', '2000']),
                          ('A18', 'harper', ['obs', 're']), ('B16.5', 'dent', ['obs'])):
        for tc in tcs:
            d = gload(tab, sm, prof, tc)
            GAS[(tab, sm, prof, tc)] = d
            ests = [('static', d['cr'], d['cGF'])] if tc in ('obs', 're') else []
            if tc != 're':
                ests += [('thin', d['sx'], d['Gthin']), ('Sob w2e-6 A', d['sx'],
                                                         d['Gsob_w2e-06_A']),
                         ('Sob w2e-6 B', d['sx'], d['Gsob_w2e-06_B']),
                         ('Sob t5 A', d['sx'], d['Gsob_t5_A']),
                         ('Sob t10 A', d['sx'], d['Gsob_t10_A']),
                         ('Sob t10 B', d['sx'], d['Gsob_t10_B'])]
            for name, x, g in ests:
                m = (x >= 1.2) & (x <= 30)
                i = np.argmax(np.where(m, g, -1))
                row(tab, sm, prof, tc, name, *[f(v, 2) for v in at(x, g, XS[:-1])],
                    f'{g[i]:.3f} ({x[i]:.1f})')
L.append('\nStatic column checks: F/F* at the top ' + ', '.join(
    f'{k[0]}/{k[1]}/{k[2]}/{k[3]}: {v["cFratio"][-1]:.2f}' for k, v in GAS.items()
    if k[3] in ('obs', 're')) + ' (the grey-Lucy photosphere is not in ck RE: Gamma_F is '
    'a flux-weighted mean and does not depend on this normalisation). Fixed-T cool brackets '
    'are NOT given as static columns: isothermal gas above its RE temperature becomes '
    'self-luminous in its own strong bands and the static k_F exceeds the thin value by 4-6x '
    '(e.g. 1500 K: 0.95-1.5 at 20 R*), an energy-violating artefact.\n')

# ------------------------------------------------------------------ dust
RS = pickle.load(open(O + 'dust_steady.pkl', 'rb'))
RB = pickle.load(open(O + 'dust_ball.pkl', 'rb'))
L.append('## Dust (a): steady-wind trajectory through the observed rho (kinematic; T_gas = '
         'observed; v capped at 12 km/s)\n')
row('star', 'Mdot', 'optics', 'alpha', 'n_seed/n_H', 'seeds grow at r/R*', 'a = 0.1 um at r/R*',
    't(1.2 R* -> r_0.1) [yr]', 'a at 5/10/40 R* [um]', 'Gamma_dust at 2/3/5/10/20 R*',
    'f_cond(40 R*)')
L.append('|' + '---|'*11)
for r in RS:
    g = np.where(r['a'] > 1.05e-7)[0]
    i1 = np.where(r['a'] > 1e-5)[0]
    row(r['star'], f"{r['mdot']:.0e}", r['var'], f"{r['alpha']:g}", f"{r['xseed']:.0e}",
        f(r['r'][g[0]]) if len(g) else '-', f(r['r'][i1[0]]) if len(i1) else '-',
        f(r['t'][i1[0]]/gl.YR, 2) if len(i1) else '-',
        '/'.join(f(v*1e4, 2) for v in at(r['r'], r['a'], [5, 10, 40])),
        '/'.join(f(v, 2) for v in at(r['r'], r['gd'], [2, 3, 5, 10, 20])), f(r['f'][-1], 2))

L.append('\n## Dust (b): ballistic parcels launched from R* at v0 (20-100 km/s), parcel '
         'density = observed rho(r), T_gas = observed, launch flux Mdot(>v0) = 4 pi R^2 rho_ph '
         'v_con exp(-v0/v_con) with (rho_ph, v_con) of the FT fit above; outcome to 40 R*\n')
row('star', 'optics', 'alpha', 'n_seed/n_H', 'Gamma_gas', 'dust-free escape v0 [km/s] / '
    'Mdot', 'min v0 for a wind with dust', 'Mdot(dust) [Msun/yr]',
    'v_inf mass-wt [km/s]', 'seeds grow at (min v0)', 'a=0.1um at', 'max Gamma_d')
L.append('|' + '---|'*12)
BALL = {}
for star in ('betA18', 'betB16.5'):
    st = gl.Star(star)
    for var in ('pure', 'Fe3e-4', 'Fe1e-3'):
        for al in (0.1, 1.0):
            for xs in (1e-15, 1e-13):
                for gg in (0., 0.1):
                    rr = sorted([r for r in RB if r['star'] == star and r['var'] == var and
                                 r['alpha'] == al and r['xseed'] == xs and r['Ggas'] == gg],
                                key=lambda r: r['v0'])
                    v0s = np.array([r['v0'] for r in rr])
                    wind_ = np.array([r['outcome'] == 'wind' for r in rr])
                    vfree = st.vesc*np.sqrt(1 - gg)*np.sqrt(1 - 1/40.)/1e5
                    mfree = st.mdot_gt(vfree*1e5)*gl.YR/gl.Msun
                    dust_w = wind_ & (v0s < vfree)
                    if dust_w.any():
                        k = np.where(dust_w)[0][0]
                        md = st.mdot_gt(v0s[k]*1e5)*gl.YR/gl.Msun
                        # mass-weighted v_inf over all escaping v0 bins (dust or free)
                        edges = np.concatenate([v0s - 1.25, [v0s[-1] + 1.25]])*1e5
                        dm = st.mdot_gt(edges[:-1]) - st.mdot_gt(edges[1:])
                        vi = np.array([r.get('vinf', np.nan) for r in rr])
                        ok = wind_ & np.isfinite(vi)
                        vinf = (dm[ok]*vi[ok]).sum()/dm[ok].sum()/1e5
                        rs, r01 = rr[k]['r_seed'], rr[k]['r_01']
                        gdm = max(r.get('gdmax', 0) for r in rr if r['outcome'] == 'wind')
                    else:
                        k, md, vinf, rs, r01 = None, 0., np.nan, np.nan, np.nan
                        gdm = max(r.get('gdmax', 0) for r in rr)
                    BALL[(star, var, al, xs, gg)] = (md, vinf, v0s[k] if k is not None else
                                                     np.nan)
                    row(star, var, f'{al:g}', f'{xs:.0e}', f'{gg:g}',
                        f'{vfree:.1f} / {mfree:.2e}', f(v0s[k]) if k is not None else 'none',
                        f'{md:.2e}' if md else '0', f(vinf, 3), f(rs), f(r01), f(gdm, 2))

# ------------------------------------------------------------------ comparison
L.append('\n## Comparison (A/18, Dent, Mdot 2e-6): extra non-radiative support needed\n')
L.append('Gamma_avail,gas = max over {static RE (table A), Sobolev t10 flux A (table A), '
         'thin 1500 K} -- a generous upper bound; v_t,eq = sqrt[(Gamma_req(v_t=0) - '
         'Gamma_avail) g H_rho] is the isotropic turbulent/wave velocity whose pressure '
         'rho v_t^2 closes the deficit; P_w = rho v_t,eq^2.\n')
w = WIND[('A', 18., 'dent', 2e-6)]
x = w['x']
gA = GAS[('A', 'A18', 'dent', 're')]
gS = GAS[('A', 'A18', 'dent', 'obs')]
g15 = GAS[('A', 'A18', 'dent', '1500')]
gst = np.interp(np.log(x), np.log(gA['cr']), gA['cGF'])
gsob = np.interp(np.log(x), np.log(gS['sx']), gS['Gsob_t10_A'])
gsobw = np.interp(np.log(x), np.log(gS['sx']), gS['Gsob_w2e-06_A'])
gth15 = np.interp(np.log(x), np.log(g15['sx']), g15['Gthin'])
gav = np.maximum.reduce([gst, gsob, gth15])
dpure = [r for r in RS if r['star'] == 'betA18' and r['mdot'] == 2e-6 and r['var'] == 'pure'
         and r['alpha'] == 1. and r['xseed'] == 1e-13][0]
dfe = [r for r in RS if r['star'] == 'betA18' and r['mdot'] == 2e-6 and r['var'] == 'Fe1e-3'
       and r['alpha'] == 0.1 and r['xseed'] == 1e-15][0]
vte = W.vt_equiv(w, gav)
row('r/R*', *[f'{v:g}' for v in XS])
L.append('|' + '---|'*(len(XS) + 1))
row('Gamma_req (v_t = 0)', *[f(v, 3) for v in at(x, w['Greq'][0.])])
row('Gamma_gas static RE', *[f(v, 2) for v in at(x, gst)])
row('Gamma_gas Sobolev wind', *[f(v, 2) for v in at(x, gsobw)])
row('Gamma_gas Sobolev t10', *[f(v, 2) for v in at(x, gsob)])
row('Gamma_gas thin, 1500 K', *[f(v, 2) for v in at(x, gth15)])
row('Gamma_dust steady, pure a=1 1e-13', *[f(v, 2) for v in at(dpure['r'], dpure['gd'])])
row('Gamma_dust steady, Fe1e-3 a=0.1 1e-15', *[f(v, 2) for v in at(dfe['r'], dfe['gd'])])
row('deficit Gamma_req - Gamma_gas,max', *[f(v, 2) for v in at(x, w['Greq'][0.] - gav)])
row('v_t,eq [km/s] (gas only)', *[f(v/1e5, 3) for v in at(x, vte)])
row('P_w [dyn/cm^2]', *[f(v, 2) for v in at(x, w['rho']*np.nan_to_num(vte)**2)])
row('P_w / P_gas', *[f(v, 2) for v in at(x, w['rho']*np.nan_to_num(vte)**2/w['P'])])
for s, M in (('A', 20.), ('B', 16.5), ('B', 19.)):
    ww = WIND[(s, M, 'dent', 2e-6)]
    gs = gav*(18./M if s == 'A' else 1.)*(1. if s == 'A' else 16.5/M)
    if s == 'B':
        gB = GAS[('A', 'B16.5', 'dent', 'obs')]
        gs = np.maximum(np.interp(np.log(ww['x']), np.log(gB['cr']), gB['cGF']),
                        np.interp(np.log(ww['x']), np.log(gB['sx']), gB['Gsob_t10_A']))
        gs = np.maximum(gs, gth15*P.stellar('A', 18.)['kE']/P.stellar(s, 16.5)['kE'])
        gs = gs*16.5/M
    row(f'v_t,eq [km/s] {s}/{M:g}', *[f(v/1e5, 3) for v in at(ww['x'], W.vt_equiv(ww, gs))])
open(O + 'tables.md', 'w').write('\n'.join(L) + '\n')

# ------------------------------------------------------------------ plots
RA = P.stellar('A', 18.)['R']
RBs = P.stellar('B', 16.5)['R']
xx = np.exp(np.linspace(np.log(1.08), np.log(40.), 400))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
xd, lnd = P._xd, P._lnd
ax[0].plot(xd*RA, np.exp(lnd)*P.MPH, 'o', ms=3, color=INK2, label='Dent+2024 Fig. 3 (digitized)')
m5 = xx <= 5
ax[0].plot(xx[m5]*RA, P.rho(xx[m5]), color=C[0], label='Dent, smoothed (set A)')
ax[0].plot(xx[~m5]*RA, P.rho(xx[~m5]), '--', color=C[0], label='r^-2.85 extrapolation')
xh = xx[xx >= 1.5]
ax[0].plot(xh*RA, P.rho(xh, 'harper'), color=C[1], label='Harper+2001, rescaled (set A)')
ax[0].plot(xx*RBs, P.rho(xx, 'dent', 168.), color=C[2], label='Dent, set B (168 pc, x1.15)')
st = P.stellar('A', 18.)
ft = P.ft_fit(st)
xf = np.maximum(xx, 1.)
ax[0].plot(xx*RA, ft['rho_ph']*xf**-2*np.exp(-ft['q']*np.sqrt(1 - 1/xf)), ':', color=C[4],
           label=f"FT eq. 7 fit (v_esc/v_con = {ft['q']:.1f})")
for md, ls in ((1e-6, (0, (1, 2))), (4e-6, (0, (4, 2)))):
    ax[0].plot(xx*RA, md*P.MSUN/P.YR/(4*np.pi*(xx*RA)**2*12e5), ls=ls, color=INK2, lw=1,
               label=f'Mdot {md:.0e}, v = 12 km/s')
ax[0].set(xscale='log', yscale='log', xlabel='r [cm]', ylabel='rho [g/cm^3]',
          ylim=(1e-19, 1e-10))
topaxis(ax[0], RA)
ax[0].legend(fontsize=7.5)
ax[1].plot(xx*RA, P.Tobs(xx), color=C[0], label='Dent points + Harper shape (adopted)')
ax[1].plot(P.DENT_T[:, 0]*RA, P.DENT_T[:, 1], 'o', color=C[0], ms=5)
ax[1].plot(xh*RA, P.Tobs(xh, 'harper'), color=C[1], label='Harper+2001 (rescaled r)')
ax[1].axhspan(1500, 2000, color=C[4], alpha=0.12, lw=0)
ax[1].axhspan(200, 900, color=C[6], alpha=0.10, lw=0)
ax[1].text(1.15*RA, 1750, 'MOLsphere bracket', fontsize=8, color=INK2)
ax[1].text(1.15*RA, 500, 'non-LTE cool roots', fontsize=8, color=INK2)
ax[1].set(xscale='log', xlabel='r [cm]', ylabel='T [K]')
topaxis(ax[1], RA)
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(PL + 'rho_T.png', dpi=150)
plt.close(fig)

fig, ax = plt.subplots(figsize=(6.5, 4.3))
for i, md in enumerate((1e-6, 2e-6, 4e-6)):
    for s, M, ls in (('A', 18., '-'), ('B', 16.5, '--')):
        ww = WIND[(s, M, 'dent', md)]
        ax.plot(ww['r'], ww['v']/1e5, ls=ls, color=C[i],
                label=f'Mdot {md:.0e} ({s})' if True else None)
ax.axhspan(10, 15, color=INK2, alpha=0.10, lw=0)
ax.text(1.1*RA, 12.3, 'observed v_inf 10-15 km/s', fontsize=8, color=INK2)
ax.plot([1.2*RA, 5*RA], [5, 5], color=INK2, lw=1)
ax.text(1.25*RA, 5.4, 'observed mean outflow < 5 km/s inside 5 R*', fontsize=8, color=INK2)
ax.set(xscale='log', yscale='log', xlabel='r [cm]', ylabel='v [km/s]', ylim=(3e-3, 30))
topaxis(ax, RA)
ax.legend(fontsize=8, ncol=2, loc='lower right')
fig.tight_layout()
fig.savefig(PL + 'v.png', dpi=150)
plt.close(fig)

fig, ax = plt.subplots(figsize=(8.5, 5))
r = w['r']
ax.plot(r, w['Greq'][0.], color=INK, label='Gamma_req, v_t = 0')
ax.plot(r, np.maximum(w['Greq'][5e5], 1e-4), color=INK, lw=1.2, ls='--', label='Gamma_req, v_t = 5')
ax.plot(r, np.maximum(w['Greq'][10e5], 1e-4), color=INK, lw=1.2, ls=':', label='Gamma_req, v_t = 10')
ax.plot(gA['cr']*RA, gA['cGF'], color=C[0], label='gas: static ck column (LTE RE)')
ax.plot(gS['sx']*RA, gS['Gthin'], color=C[2], label='gas: thin (observed T)')
ax.plot(gS['sx']*RA, gS['Gsob_w2e-06_A'], color=C[1], label='gas: Sobolev, steady-wind dv/dr')
ax.plot(gS['sx']*RA, gS['Gsob_t10_A'], color=C[3], label='gas: Sobolev, 10 km/s / H_rho')
ax.plot(g15['sx']*RA, g15['Gthin'], color=C[4], label='gas: thin, T = 1500 K bracket')
ax.plot(dpure['r']*RA, np.maximum(dpure['gd'], 1e-4), color=C[5],
        label='dust (steady): pure, alpha 1, 1e-13')
ax.plot(dfe['r']*RA, np.maximum(dfe['gd'], 1e-4), color=C[7],
        label='dust (steady): Fe1e-3, alpha 0.1, 1e-15')
ax.set(xscale='log', yscale='log', xlabel='r [cm]', ylabel='Gamma (force / gravity)',
       ylim=(1e-3, 20), xlim=(1.2*RA, 30*RA))
topaxis(ax, RA)
ax.legend(fontsize=7.5, loc='lower right', ncol=2)
fig.suptitle('Set A, M = 18 Msun, Mdot = 2e-6 Msun/yr, Dent profile', fontsize=10, color=INK)
fig.tight_layout()
fig.savefig(PL + 'gamma.png', dpi=150)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7.5, 4.5))
for i, var in enumerate(('pure', 'Fe3e-4', 'Fe1e-3')):
    for al, xs, ls in ((1., 1e-13, '-'), (0.1, 1e-15, '--')):
        rr = [q for q in RS if q['star'] == 'betA18' and q['mdot'] == 2e-6 and q['var'] == var
              and q['alpha'] == al and q['xseed'] == xs][0]
        ax.plot(rr['r']*RA, rr['a']*1e4, ls=ls, color=C[i],
                label=f'{var}, alpha {al:g}, seeds {xs:.0e}')
ax.axhline(0.1, color=INK2, lw=1)
for x0, x1, lab in ((1.5, 2.0, 'Haubois+19'), (12., 14., 'MATISSE 13 R*'),
                    (24., 40., 'VISIR shell 0.5-1"')):
    ax.axvspan(x0*RA, x1*RA, color=INK2, alpha=0.10, lw=0)
    ax.text(x0*RA*1.02, 0.25, lab, fontsize=7.5, color=INK2, rotation=90)
ax.set(xscale='log', yscale='log', xlabel='r [cm]', ylabel='grain radius a [um]',
       ylim=(8e-4, 1), xlim=(1.2*RA, 40*RA))
topaxis(ax, RA)
ax.legend(fontsize=7.5, loc='lower right')
fig.suptitle('Steady wind, Mdot 2e-6, set A, observed rho and T', fontsize=10, color=INK)
fig.tight_layout()
fig.savefig(PL + 'grain_a.png', dpi=150)
plt.close(fig)
print('\n'.join(L))
