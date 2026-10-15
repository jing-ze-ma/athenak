"""Diagnostics PNG of a standard IC column (tools/ic_star).  One y-scale per panel."""
import numpy as np

C1, C2, C3, C4, C5 = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#6250d6'
INK2, GRID = '#52514e', '#e4e3df'


def png(fn, col, pr, info, mesa=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    r, R = col['r'], col['rph']
    x = r/R
    m = pr['below'] & (r >= col['r'][0])
    fig, ax = plt.subplots(3, 2, figsize=(12, 11), constrained_layout=True)
    for a in ax.flat:
        a.grid(True, color=GRID, lw=0.6)
        a.tick_params(colors=INK2, labelsize=8)
        for sp in a.spines.values():
            sp.set_color(GRID)
        a.axvline(info['x1min']/R, color=INK2, lw=0.8, ls='--')
        if info.get('r_base'):
            a.axvspan(info['r_base']/R, info['r_top']/R, color='#f0efec', zorder=0)
        a.set_xlabel('r / R_ph', color=INK2, fontsize=9)
    xm = min(1.5, info['x1max']/R)
    a = ax[0, 0]
    a.semilogy(x, col['T'], color=C1, lw=1.6, label='T [K] (column)')
    if mesa is not None:
        a.semilogy(mesa['rc']/R, mesa['T'], color=C1, lw=1, ls=':', label='T MESA')
    a.set_xlim(x[0], xm)
    a.set_title('temperature (dashed line: x1min; shaded: Fe CZ)', fontsize=10)
    a.legend(fontsize=8, frameon=False)
    a = ax[0, 1]
    a.semilogy(x, col['rho'], color=C2, lw=1.6, label='rho [g/cc] (column)')
    if mesa is not None:
        a.semilogy(mesa['rc']/R, mesa['rho'], color=C2, lw=1, ls=':', label='rho MESA')
    a.axhline(info['dmax'], color=C4, lw=1, label='sponge dmax')
    a.axhline(info['rho_floor'], color=C5, lw=1, ls='--', label='column floor')
    a.set_xlim(x[0], info['x1max']/R*1.02)
    a.set_title('density, floor layer and sponge threshold', fontsize=10)
    a.legend(fontsize=8, frameon=False)
    a = ax[1, 0]
    a.plot(x[m], pr['nab'][m], color=C1, lw=1.6, label='nabla (column, MLT++)')
    a.plot(x[m], pr['nad'][m], color=C2, lw=1.2, label='nabla_ad')
    a.plot(x[m], pr['nrad'][m], color=C3, lw=1.2, label='nabla_rad')
    a.set_ylim(0, min(2.0, np.nanmax(pr['nrad'][m])*1.05))
    a.set_xlim(x[0], 1.0)
    a.set_title('gradients', fontsize=10)
    a.legend(fontsize=8, frameon=False)
    a = ax[1, 1]
    a.plot(x[m], col['fmlt'][m], color=C1, lw=1.6, label='F_MLT/F')
    a.plot(x[m], pr['gam_L'][m], color=C2, lw=1.2, label='Gamma_L = kappa L/(4 pi c G M)')
    a.plot(x[m], pr['gam_rad'][m], color=C3, lw=1.2, label='Gamma_rad = kappa F_r/(c g)')
    a.axhline(1.0, color=INK2, lw=0.6)
    a.set_xlim(x[0], 1.0)
    a.set_title('convective flux fraction and Eddington factors', fontsize=10)
    a.legend(fontsize=8, frameon=False)
    a = ax[2, 0]
    a.semilogy(x[m], pr['kap'][m], color=C4, lw=1.6, label='kappa_R [cm^2/g]')
    if info.get('r_fe'):
        a.axvline(info['r_fe']/R, color=C2, lw=1, label='Fe opacity peak')
    a.set_xlim(x[0], 1.0)
    a.set_title('opacity', fontsize=10)
    a.legend(fontsize=8, frameon=False)
    a = ax[2, 1]
    a.plot(x[m], pr['Lsum'][m]/col['L'], color=C1, lw=1.6,
           label='4 pi r^2 (F_r + F_MLT) / L')
    dd = pr['Fd']/col['Fr']
    mm = m & (col['fmlt'] < 1e-6) & (r < col['rj'])
    a.plot(x[mm], dd[mm], '.', ms=2, color=C2, label='F_diff(column T grad) / F_r (radiative rows)')
    if mesa is not None:
        a.plot(mesa['r']/R, mesa['L']/col['L'], color=C3, lw=1, ls=':', label='L_MESA(r) / L')
    a.set_ylim(0.9, 1.1)
    a.set_xlim(x[0], 1.0)
    a.set_title('luminosity consistency', fontsize=10)
    a.legend(fontsize=8, frameon=False)
    fig.suptitle('%s: standard IC column (tools/ic_star)' % info['name'], fontsize=11)
    fig.savefig(fn, dpi=110)
    plt.close(fig)
