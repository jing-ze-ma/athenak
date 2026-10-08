"""Tables + plots from out/grid.pkl and the optics/chem tables.
usage: venv/bin/python report.py > out/report_tables.md"""
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa
import grainlib as gl  # noqa
from run_grid import V0, ALPHA, XSEED, VAR, TGAS, GGAS  # noqa

D = gl.D
PL = '/orion/ptmp/jinma/rsg_wind_1008/grains/plots/'
res = pickle.load(open(D + 'grid.pkl', 'rb'))
for _r in res:   # dust-driven wind: escapes, but the dust-free orbit is bound inside 20 R
    _r['ballistic'] = not (_r['rmax0'] < 20)
    _r['dw'] = (_r['outcome'] == 'wind') and not _r['ballistic'] and _r['vinf'] > 0
ch = gl.Chem()
STARS = list(gl.STARS)


def sel(**kw):
    return [r for r in res if all(r[k] == v for k, v in kw.items())]


# ---------------- optics + Td plot ----------------
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
cols = dict(zip(VAR, ['C0', 'C2', 'C1', 'C3']))
aa = np.geomspace(1e-7, 3e-4, 80)
for var in VAR:
    m = gl.Model('golden16', var, ch)
    ax[0].loglog(aa * 1e4, [m.kpr_dust(a) for a in aa], c=cols[var], label=var)
    for k, sn in enumerate(STARS):
        mm = m if sn == 'golden16' else gl.Model(sn, var, ch)
        for rr, ls in zip([3, 4, 5], ['-', '--', ':']):
            ax[1 + k].semilogx(aa * 1e4, [mm.Td_at(a, rr * mm.st.R) for a in aa], c=cols[var],
                               ls=ls, label='%s %dR' % (var, rr) if var == 'pure' else None)
for k, sn in enumerate(STARS):
    st = gl.Star(sn)
    for rr, ls in zip([3, 4, 5], ['-', '--', ':']):
        nH = st.rho_ft(rr * st.R) / ch.mH
        T = np.arange(700, 1500, 1.)
        Tc = T[np.argmin(np.abs([ch.lnS(t, nH, 0) for t in T]))]
        ax[1 + k].axhline(Tc, c='k', ls=ls, lw=0.8)
    ax[1 + k].set_title('%s: T_d(a) at 3/4/5 R (black: forsterite S=1, FT rho)' % sn,
                        fontsize=9)
    ax[1 + k].set_xlabel('a [um]')
    ax[1 + k].set_ylabel('T_d [K]')
    ax[1 + k].legend(fontsize=7)
ax[0].axhspan(300, 900, color='0.85', label='needed f_cond kappa (note sec. 4)')
ax[0].set_xlabel('a [um]')
ax[0].set_ylabel('kappa_pr [cm^2/g dust], Planck 4106 K flux mean')
ax[0].legend(fontsize=8)
fig.tight_layout()
fig.savefig(PL + 'optics_Td.png', dpi=110)

print('## kappa_pr per gram dust (golden16 Teff), cm^2/g; Q_pr\n')
print('| a [um] | ' + ' | '.join(VAR) + ' | Q_pr (pure) |\n' + '|---' * (len(VAR) + 2) + '|')
mp = {v: gl.Model('golden16', v, ch) for v in VAR}
for a in [1e-7, 1e-6, 1e-5, 3e-5, 5e-5, 1e-4, 3e-4]:
    print('| %.3g | ' % (a * 1e4) + ' | '.join('%.0f' % mp[v].kpr_dust(a) for v in VAR)
          + ' | %.3f |' % np.interp(np.log(a), mp['pure'].op.la, mp['pure'].qpr))

print('\n## T_d [K] (a = 1 nm / 0.1 um / 0.3 um) and forsterite T(S=1) at the FT density\n')
print('| star | r/R | T(S=1) | ' + ' | '.join(VAR) + ' |\n|---|---|---' + '|---' * len(VAR) + '|')
for sn in STARS:
    ms = {v: gl.Model(sn, v, ch) for v in VAR}
    st = ms['pure'].st
    for rr in [2, 3, 4, 5]:
        nH = st.rho_ft(rr * st.R) / ch.mH
        T = np.arange(700, 1500, 1.)
        Tc = T[np.argmin(np.abs([ch.lnS(t, nH, 0) for t in T]))]
        print('| %s | %d | %.0f | ' % (sn, rr, Tc) + ' | '.join(
            '/'.join('%.0f' % ms[v].Td_at(a, rr * st.R) for a in [1e-7, 1e-5, 3e-5])
            for v in VAR) + ' |')

# ---------------- outcome map ----------------
print('\n## Outcome map (dens = FT time-averaged): number of DUST-DRIVEN wind cases (dust-free orbit bound) out of '
      '|Tg| x |Ggas| = 12 per cell; rows v0 [km/s], columns n_seed/n_H\n')
for sn in STARS:
    for var in VAR:
        for al in ALPHA:
            rows = []
            for v0 in V0:
                c = []
                for xs in XSEED:
                    s = sel(star=sn, var=var, alpha=al, v0=v0, xseed=xs, dens='ft')
                    c.append(sum(r['dw'] for r in s))
                rows.append(c)
            rows = np.array(rows)
            if rows.sum() == 0:
                print('- %s %s alpha=%g: no wind anywhere' % (sn, var, al))
                continue
            print('\n%s %s alpha=%g (xseed %s)' % (sn, var, al, XSEED))
            for v0, c in zip(V0, rows):
                if c.sum():
                    print('  v0=%d: %s' % (v0, c))

fig, axs = plt.subplots(len(STARS) * 2, len(VAR) * len(ALPHA) // 1, figsize=(20, 8.5),
                        squeeze=False)
for i, sn in enumerate(STARS):
    for jd, dens in enumerate(['ft', 'mc']):
        for j, (var, al) in enumerate([(v, a) for v in VAR for a in ALPHA]):
            M = np.zeros((len(V0), len(XSEED)))
            for iv, v0 in enumerate(V0):
                for ix, xs in enumerate(XSEED):
                    s = sel(star=sn, var=var, alpha=al, v0=v0, xseed=xs, dens=dens)
                    M[iv, ix] = np.mean([r['dw'] for r in s])
            a_ = axs[2 * i + jd, j]
            a_.imshow(M, origin='lower', vmin=0, vmax=1, cmap='viridis', aspect='auto')
            a_.set_xticks(range(len(XSEED)))
            a_.set_xticklabels(['-16', '-15', '-14', '-13'], fontsize=6)
            a_.set_yticks(range(len(V0)))
            a_.set_yticklabels(V0, fontsize=6)
            if 2 * i + jd == 0:
                a_.set_title('%s a=%g' % (var, al), fontsize=7)
            if j == 0:
                a_.set_ylabel('%s %s\nv0 [km/s]' % (sn, dens), fontsize=7)
fig.suptitle('dust-driven wind fraction over (Tg x Gamma_gas) = 12 cases; x: log n_seed/n_H', fontsize=10)
fig.tight_layout()
fig.savefig(PL + 'outcome_map.png', dpi=110)

# ---------------- summary stats ----------------
print('\n## Wind counts by parameter (dens=ft / mc)\n')
for sn in STARS:
    for key, vals in [('var', VAR), ('alpha', ALPHA), ('xseed', XSEED), ('Tg', TGAS),
                      ('Ggas', GGAS)]:
        line = []
        for v in vals:
            nf = [sum(r['dw'] for r in sel(star=sn, dens=d, **{key: v}))
                  for d in ['ft', 'mc']]
            line.append('%s=%s: %d/%d' % (key, v, nf[0], nf[1]))
        print('- %s: %s' % (sn, '; '.join(line)))

print('\n## Minimum wind v0 and implied Mdot (dens=ft, Ggas=0.1, Tg=800)\n')
print('| star | var | alpha | xseed | v0_min | Mdot(>v0_min) [Msun/yr] | vinf [km/s] | '
      'amax [um] | fmax | Gd max | vd max [km/s] | t_0.1um [s] | t(r>3R) [s] |')
print('|---' * 13 + '|')
for sn in STARS:
    for var in VAR:
        for al in ALPHA:
            for xs in XSEED:
                s = sorted(sel(star=sn, var=var, alpha=al, xseed=xs, dens='ft', Ggas=0.1,
                               Tg=800.), key=lambda r: r['v0'])
                w = [r for r in s if r['dw']]
                if not w:
                    continue
                r = w[0]
                print('| %s | %s | %g | %g | %d | %.1e | %.0f | %.3f | %.2f | %.2f | %.0f | '
                      '%.1e | %.1e |' % (sn, var, al, xs, r['v0'], r['mdot'] / gl.Msun * gl.YR,
                                         r['vinf'] / 1e5, r['amax'] * 1e4, r['fmax'], r['gdmax'],
                                         r['vdmax'] / 1e5, r['t_01um'], r['t3']))

print('\n## Dust-free orbits (r_max/R, time at r>3R for dust-free = Gamma_gas only)\n')
print('| star | Ggas | ' + ' | '.join('v0=%d' % v for v in V0) + ' |\n|---|---' +
      '|---' * len(V0) + '|')
for sn in STARS:
    for gg in GGAS:
        cells = []
        for v0 in V0:
            s = sel(star=sn, var='Fe1e-3', alpha=0.01, xseed=1e-16, Tg=800., Ggas=gg, v0=v0,
                    dens='ft')[0]
            cells.append('%.2f / %.1e' % (s['rmax0'], s['t3']) if np.isfinite(s['rmax0'])
                         else 'unbound')
        print('| %s | %g | ' % (sn, gg) + ' | '.join(cells) + ' |')

# ---------------- representative trajectories ----------------
cases = [('golden16', 'pure', 85, 1.0, 1e-14, 800., 0.1),
         ('golden16', 'pure', 85, 0.1, 1e-14, 800., 0.1),
         ('golden16', 'pure', 80, 1.0, 1e-13, 800., 0.1),
         ('golden16', 'Fe3e-4', 85, 1.0, 1e-14, 800., 0.1),
         ('m20lgl5.5', 'pure', 60, 1.0, 1e-15, 800., 0.1),
         ('m20lgl5.5', 'pure', 60, 0.1, 1e-15, 800., 0.1),
         ('m20lgl5.5', 'Fe1e-3', 70, 1.0, 1e-14, 800., 0.1),
         ('m20lgl5.5', 'pure', 50, 1.0, 1e-14, 400., 0.3)]
fig, ax = plt.subplots(2, 2, figsize=(12, 8))
for k, (sn, var, v0, al, xs, tg, gg) in enumerate(cases):
    m = gl.Model(sn, var, ch)
    r = m.run(v0 * 1e5, al, xs, tg, gg)
    lab = '%s %s v0=%d a=%g x=%g -> %s' % (sn[:6], var, v0, al, xs, r['outcome'])
    yr = r['t'] / gl.YR
    ax[0, 0].plot(yr, r['r'], label=lab)
    ax[0, 1].semilogy(yr, r['a'] * 1e4)
    ax[1, 0].plot(yr, r['gd'] + gg)
    ax[1, 1].plot(yr, r['vd'] / 1e5)
ax[0, 0].set_ylabel('r/R')
ax[0, 0].set_ylim(1, 20)
ax[0, 0].legend(fontsize=6)
ax[0, 1].set_ylabel('a [um]')
ax[1, 0].set_ylabel('Gamma_gas + Gamma_dust')
ax[1, 0].axhline(1, c='k', lw=0.5)
ax[1, 0].set_ylim(0, 5)
ax[1, 1].set_ylabel('drift v_d [km/s] (diagnostic)')
for a_ in ax.ravel():
    a_.set_xlabel('t [yr]')
    a_.set_xlim(0, 60)
fig.tight_layout()
fig.savefig(PL + 'trajectories.png', dpi=110)
