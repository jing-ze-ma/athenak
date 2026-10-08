"""Figures + summary table for the column runs.  usage: python plots.py"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import rsglib as L  # noqa: E402

cases = ['marcs', 'hse1', 'hse3', 'hse10', 'w6v10', 'w6v30', 'w5v10', 'w5v30']
col = dict(zip(cases, plt.cm.tab10(np.arange(8))))


def ld(star, c, m, tag=''):
    return np.load(f'col_{star}_{c}_{m}{tag}.npz')


rows = []
for star in ('golden16', 'betelgeuse'):
    fig, ax = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for c in cases:
        for m, ls in (('clamp', '-'), ('extrap', '--')):
            d = ld(star, c, m)
            lab = c if m == 'clamp' else None
            ax[0, 0].plot(d['r'], d['GF'], ls, color=col[c], label=lab, lw=1.2)
            ax[0, 1].plot(np.maximum(d['tauR'], 1e-12), d['kF']/d['kR'], ls, color=col[c],
                          label=lab, lw=1.2)
        d = ld(star, c, 'clamp')
        ax[0, 0].plot(d['r'], d['kFthin']/d['kE'], ':', color=col[c], lw=0.8)
        ax[1, 0].plot(d['r'], d['T'], '-', color=col[c], label=c, lw=1.2)
        if c != 'marcs':
            ax[1, 0].plot(d['r'], d['T0'], ':', color=col[c], lw=0.8)
        for m in ('clamp', 'extrap'):
            dd = ld(star, c, m)
            ok = dd['tauR'] < 1
            i = np.argmax(np.where(ok, dd['GF'], -1))
            gtag = ''
            try:
                g0 = ld(star, c, m, '_grey')
                gtag = f'{np.max(np.where(g0["tauR"] < 1, g0["GF"], 0)):.3f}'
            except FileNotFoundError:
                pass
            rows.append((star, c, m, dd['GF'][i], dd['r'][i], dd['T'][i], dd['rho'][i],
                         dd['tauR'][i], dd['kF'][i]/dd['kR'][i], dd['kFthin'][i]/dd['kE'],
                         dd['Fratio'][-1], gtag,
                         dd['hist'][-1] if len(dd['hist']) else np.nan))
    ax[0, 0].axhline(1, color='k', lw=1.5)
    ax[0, 0].set(xscale='log', yscale='log', xlabel='r / R*', ylim=(1e-4, 3),
                 ylabel='Gamma_F = kappa_F / kappa_Edd',
                 title=f'{star}: solid clamp, dashed extrap (P<1e-8 bar), dotted thin limit')
    ax[0, 0].legend(fontsize=8, ncol=2)
    ax[0, 1].set(xscale='log', yscale='log', xlabel='tau_R (from top)',
                 ylabel='kappa_F / kappa_R', title='flux mean over Rosseland')
    ax[0, 1].invert_xaxis()
    ax[1, 0].set(xscale='log', xlabel='r / R*', ylabel='T [K]', ylim=(800, 6000),
                 title='T: solid after lambda iterations (MARCS: its own RE T), dotted grey')
    ax[1, 0].legend(fontsize=8, ncol=2)
    wc = L.band_centers()
    d0 = ld(star, 'hse1', 'clamp')
    ax[1, 1].step(range(L.NB), d0['fb_bb'], where='mid', color='k', label='B_nu(Teff)')
    for c in ('marcs', 'hse1', 'hse10', 'w5v10'):
        d = ld(star, c, 'clamp')
        ax[1, 1].step(range(L.NB), d['Fb_top'], where='mid', color=col[c], label=c)
    ax[1, 1].set_xticks(range(L.NB))
    ax[1, 1].set_xticklabels([f'{w:.2g}' for w in wc], rotation=45)
    ax[1, 1].set(xlabel='band centre [um] (descending)', ylabel='fraction of flux at top',
                 title='emergent band fluxes')
    ax[1, 1].legend(fontsize=8)
    fig.savefig(f'col_{star}.png', dpi=110)

with open('col_table.md', 'w') as f:
    f.write('| star | case | P<1e-8 bar | max Gamma_F (tau_R<1) | r/R* | T [K] | rho | tau_R |'
            ' kF/kR | thin Gamma_F there | F_top/F* | max Gamma_F grey T | last LI dT/T |\n')
    f.write('|---' * 13 + '|\n')
    for r in rows:
        f.write(f'| {r[0]} | {r[1]} | {r[2]} | {r[3]:.3f} | {r[4]:.3f} | {r[5]:.0f} | '
                f'{r[6]:.1e} | {r[7]:.1e} | {r[8]:.3g} | {r[9]:.3f} | {r[10]:.2f} | {r[11]} |'
                f' {r[12]:.2f} |\n')
print(open('col_table.md').read())

# band decomposition of the thin flux mean
with open('band_thin.md', 'w') as f:
    wc = L.band_centers()
    f.write('| T, rho | ' + ' | '.join(f'{w:.2g}um' for w in wc) + ' | total Gamma_F |\n')
    f.write('|---' * (L.NB + 2) + '|\n')
    st = L.star('golden16')
    fb = L.fband(st['Teff'])[0]
    for T, rho in ((1500, 1e-16), (2000, 1e-15), (2500, 1e-14), (3000, 1e-10),
                   (3500, 1e-9)):
        s = L.state(T, rho)
        c = fb*(L.GW*s['K']).sum(1)/st['kE']
        f.write(f'| {T}, {rho:.0e} | ' + ' | '.join(f'{x:.3f}' for x in c)
                + f' | {c.sum():.3f} |\n')
print(open('band_thin.md').read())
