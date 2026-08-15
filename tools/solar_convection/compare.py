#!/usr/bin/env python3
"""Compare the solar_convection ideal-gas run with the tabulated general-EOS run.

    python3 compare.py [snapshot_index|-1] [ideal_run] [table_run]

e.g.  python3 compare.py            -> last snapshot of runs 'ideal' and 'table'
      python3 compare.py -1 ideal192 table192

Reads soleos/<run>/bin/sun.hydro_w.*.bin, re-evaluates each run's own EOS offline
(tools/eoslib.py) to recover T, mu, Gamma_1 and grad_ad, and writes

    soleos/plots/profiles<tag>.png      horizontally averaged vertical structure
    soleos/plots/granulation<tag>.png   photospheric vz and T, both runs
    soleos/plots/history<tag>.png       mass and total energy vs time
    soleos/summary<tag>.txt             the numbers

where <tag> is empty for the default ideal/table pair and '_<ideal>_<table>' otherwise,
so a second resolution does not overwrite the first. Truncated dumps (see tools/dumps.py)
are skipped, and a negative index resolves to the last INTACT snapshot.

The two runs use IDENTICAL grids, boundary conditions, opacity and radiative transfer.
The only difference is the equation of state, so every difference below is the EOS.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(_HERE, os.pardir, os.pardir, 'vis', 'python'))
sys.path.insert(0, _HERE)
import bin_convert                                    # noqa: E402
import dumps                                          # noqa: E402
import eoslib                                         # noqa: E402

ROOT = os.environ.get('SOLEOS_ROOT', '/orion/u/jinma/ATHENAK/soleos')
GRID = os.path.join(ROOT, 'tools', 'eos_grid.bin')
PLOTS = os.path.join(ROOT, 'plots')
GRAV = 2.74e4          # |g|, cm/s^2 (constant, as in the pgen)
SIGMA = 5.6704e-5
TEFF = 5778.0          # the Teff the IC is built to radiate


# ---------------------------------------------------------------------------------------
def load(run, idx):
    """Assemble one snapshot onto the global (nx3, nx2, nx1) grid."""
    files = dumps.dump_files(os.path.join(ROOT, run))
    path = files[dumps.resolve(files, idx)]
    fd = bin_convert.read_binary(path)
    nx1, nx2, nx3 = fd['Nx1'], fd['Nx2'], fd['Nx3']
    m1, m2, m3 = fd['nx1_mb'], fd['nx2_mb'], fd['nx3_mb']
    out = {}
    for var in ['dens', 'velx', 'vely', 'velz', 'eint']:
        g = np.empty((nx3, nx2, nx1), dtype=np.float64)
        blocks = np.asarray(fd['mb_data'][var])
        for b in range(fd['n_mbs']):
            i, j, k = fd['mb_logical'][b][:3]
            g[k*m3:(k+1)*m3, j*m2:(j+1)*m2, i*m1:(i+1)*m1] = blocks[b]
        out[var] = g
    dx1 = (fd['x1max'] - fd['x1min'])/nx1
    out['z'] = fd['x1min'] + dx1*(np.arange(nx1) + 0.5)
    out['dz'] = dx1
    out['time'] = fd['time']
    out['file'] = os.path.basename(path)
    return out


def analyse(run, eos, idx):
    d = load(run, idx)
    rho = d['dens']
    espec = d['eint']/rho                       # eint is the internal energy DENSITY
    T = eos.temperature(rho, espec)
    p = eos.pressure(rho, espec)
    mu = eos.mu(rho, espec)
    g1 = eos.gamma1(rho, espec)
    gad = eos.grad_ad(rho, espec)
    cs = np.sqrt(g1*p/rho)
    vsq = d['velx']**2 + d['vely']**2 + d['velz']**2

    ax = (0, 1)                                  # average over the two horizontal axes
    prof = {
        'z': d['z'], 'dz': d['dz'], 'time': d['time'], 'file': d['file'],
        'rho': rho.mean(ax), 'T': T.mean(ax), 'p': p.mean(ax),
        'mu': mu.mean(ax), 'gamma1': g1.mean(ax), 'grad_ad': gad.mean(ax),
        'vx_rms': np.sqrt((d['velx']**2).mean(ax)),
        'v_rms': np.sqrt(vsq.mean(ax)),
        'mach': np.sqrt((vsq/cs**2).mean(ax)),
        'cs': cs.mean(ax),
    }
    # actual stratification, from the horizontal means (what Schwarzschild compares)
    lnp, lnT = np.log(prof['p']), np.log(prof['T'])
    prof['grad'] = np.gradient(lnT, lnp)
    # optical depth, integrated DOWN from the top with the pgen's own opacity
    kap = eoslib.get_kapr(rho, T).mean(ax)
    prof['kappa'] = kap
    dtau = kap*prof['rho']*d['dz']
    prof['tau'] = np.cumsum(dtau[::-1])[::-1]
    # photosphere: where tau crosses 2/3
    prof['z_ph'], prof['T_ph'] = photosphere(prof['z'], prof['tau'], prof['T'])
    # 3D fields kept for the granulation slices
    prof['_T3'] = T
    prof['_vx3'] = d['velx']
    prof['_tau3'] = np.cumsum((eoslib.get_kapr(rho, T)*rho*d['dz'])[:, :, ::-1],
                              axis=2)[:, :, ::-1]
    if isinstance(eos, eoslib.TableEOS):
        xh2, xhii = eos.ion_frac(rho, espec)
        prof['xh2'] = xh2.mean(ax)
        prof['xhii'] = xhii.mean(ax)
    return prof


def photosphere(z, tau, T):
    i = np.argmax(tau <= 2.0/3.0)
    if i == 0 or tau[0] < 2.0/3.0:
        return np.nan, np.nan
    w = (np.log(tau[i-1]) - np.log(2.0/3.0))/(np.log(tau[i-1]) - np.log(tau[i]))
    return z[i-1] + w*(z[i] - z[i-1]), T[i-1] + w*(T[i] - T[i-1])


def read_hst(run):
    f = os.path.join(ROOT, run, 'sun.hydro.hst')
    a = np.loadtxt(f)
    return {'t': a[:, 0], 'mass': a[:, 2], 'etot': a[:, 6],
            'ke': a[:, 7] + a[:, 8] + a[:, 9]}


# ---------------------------------------------------------------------------------------
def main():
    idx = int(sys.argv[1]) if len(sys.argv) > 1 else -1
    iname = sys.argv[2] if len(sys.argv) > 2 else 'ideal'
    tname = sys.argv[3] if len(sys.argv) > 3 else 'table'
    tag = '' if (iname, tname) == ('ideal', 'table') else '_%s_%s' % (iname, tname)
    os.makedirs(PLOTS, exist_ok=True)
    # EOS per run from its own input file, NOT by position: comparing two tabulated runs
    # (a control against a variant) would otherwise silently analyse the first as ideal.
    ideal = analyse(iname, eoslib.eos_for_run(os.path.join(ROOT, iname), GRID), idx)
    table = analyse(tname, eoslib.eos_for_run(os.path.join(ROOT, tname), GRID), idx)
    runs = [(iname, ideal, 'tab:blue'), (tname, table, 'tab:red')]

    # ---- profiles ---------------------------------------------------------------------
    fig, axs = plt.subplots(3, 3, figsize=(16, 12))
    zmax = ideal['z'][-1]

    def setup(ax, xlabel, logx=False):
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$z\ /\ z_{\rm max}$')
        if logx:
            ax.set_xscale('log')
        ax.grid(alpha=0.3)

    for name, r, c in runs:
        zz = r['z']/zmax
        axs[0, 0].plot(r['T'], zz, c, label=name)
        axs[0, 1].plot(r['rho'], zz, c, label=name)
        axs[0, 2].plot(r['mu'], zz, c, label=name)
        axs[1, 0].plot(r['gamma1'], zz, c, label=name)
        axs[1, 1].plot(r['grad_ad'], zz, c, ls='--', label=r'$\nabla_{\rm ad}$ ' + name)
        axs[1, 1].plot(r['grad'], zz, c, ls='-', label=r'$\nabla$ ' + name)
        axs[1, 2].plot(r['grad'] - r['grad_ad'], zz, c, label=name)
        axs[2, 0].plot(r['vx_rms']/1e5, zz, c, label=name)
        axs[2, 1].plot(r['mach'], zz, c, label=name)
        axs[2, 2].plot(r['tau'], zz, c, label=name)
        for ax in axs.flat:
            if not np.isnan(r['z_ph']):
                ax.axhline(r['z_ph']/zmax, color=c, lw=0.7, ls=':')

    setup(axs[0, 0], 'T [K]')
    setup(axs[0, 1], r'$\rho$ [g cm$^{-3}$]', logx=True)
    setup(axs[0, 2], r'$\mu$')
    setup(axs[1, 0], r'$\Gamma_1$')
    setup(axs[1, 1], r'$\nabla$ (solid), $\nabla_{\rm ad}$ (dashed)')
    setup(axs[1, 2], r'$\nabla-\nabla_{\rm ad}$  (>0 unstable)')
    setup(axs[2, 0], r'$v_{z,\rm rms}$ [km s$^{-1}$]')
    setup(axs[2, 1], 'Mach')
    setup(axs[2, 2], r'$\tau$', logx=True)
    axs[1, 2].axvline(0.0, color='k', lw=0.8)
    axs[2, 2].axvline(2.0/3.0, color='k', lw=0.8)
    if 'xhii' in table:
        tw = axs[0, 2].twiny()
        tw.plot(table['xhii'], table['z']/zmax, 'k-', lw=1, alpha=0.6)
        tw.plot(table['xh2'], table['z']/zmax, 'k--', lw=1, alpha=0.6)
        tw.set_xlabel(r'table: $x_{\rm H^+}$ (solid), $x_{\rm H_2}$ (dashed)', fontsize=9)
    for ax in axs.flat:
        ax.legend(fontsize=8)
    fig.suptitle('solar_convection: ideal gas vs tabulated general EOS (no radiation)   '
                 't = %.0f / %.0f s' % (ideal['time'], table['time']))
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, 'profiles%s.png' % tag), dpi=110)
    plt.close(fig)

    # ---- granulation ------------------------------------------------------------------
    fig, axs = plt.subplots(2, 2, figsize=(11, 10))
    for col, (name, r, _) in enumerate(runs):
        k = int(np.argmin(np.abs(r['tau'] - 2.0/3.0)))
        vz = r['_vx3'][:, :, k]/1e5
        tt = r['_T3'][:, :, k]
        lim = np.percentile(np.abs(vz), 99.0)
        im = axs[0, col].imshow(vz, origin='lower', cmap='seismic_r',
                                vmin=-lim, vmax=lim)
        axs[0, col].set_title('%s: $v_z$ [km/s] at $\\tau=2/3$ (z=%.2f)'
                              % (name, r['z'][k]/zmax))
        plt.colorbar(im, ax=axs[0, col], fraction=0.046)
        im = axs[1, col].imshow(tt, origin='lower', cmap='gist_heat')
        axs[1, col].set_title('%s: T [K] at $\\tau=2/3$' % name)
        plt.colorbar(im, ax=axs[1, col], fraction=0.046)
    fig.suptitle('photospheric granulation')
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, 'granulation%s.png' % tag), dpi=110)
    plt.close(fig)

    # ---- history ----------------------------------------------------------------------
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    for name, _, c in runs:
        h = read_hst(name)
        axs[0].plot(h['t'], 100*(h['mass']/h['mass'][0] - 1), c, label=name)
        axs[1].plot(h['t'], 100*(h['etot']/h['etot'][0] - 1), c, label=name)
        axs[2].plot(h['t'], h['ke'], c, label=name)
    for ax, lab in zip(axs, ['mass drift [%]', 'total energy drift [%]',
                             'kinetic energy [erg]']):
        ax.set_xlabel('t [s]')
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
        ax.legend()
    axs[2].set_yscale('log')
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, 'history%s.png' % tag), dpi=110)
    plt.close(fig)

    # ---- text summary -----------------------------------------------------------------
    lines = []

    def w(s=''):
        lines.append(s)
        print(s)

    w('solar_convection: ideal gas vs tabulated general EOS '
      '(H2 + H/He Saha, no radiation)')
    w('=' * 88)
    w('snapshot: %s %s (t=%.0f s), %s %s (t=%.0f s)'
      % (iname, ideal['file'], ideal['time'], tname, table['file'], table['time']))
    w()
    w('%-34s %14s %14s' % ('', iname, tname))
    w('-' * 88)

    def row(lab, a, b, fmt='%14.4g'):
        w('%-34s ' % lab + fmt % a + ' ' + fmt % b)

    for lab, key, i in [('base (z=0)  T [K]', 'T', 0),
                        ('base        rho [g/cm^3]', 'rho', 0),
                        ('base        p [erg/cm^3]', 'p', 0),
                        ('base        mu', 'mu', 0),
                        ('base        Gamma_1', 'gamma1', 0),
                        ('base        grad_ad', 'grad_ad', 0)]:
        row(lab, ideal[key][i], table[key][i])
    w()
    row('photosphere z/zmax', ideal['z_ph']/zmax, table['z_ph']/zmax)
    row('photosphere T [K]  (Teff=5778)', ideal['T_ph'], table['T_ph'])
    kph_i = int(np.argmin(np.abs(ideal['tau'] - 2.0/3.0)))
    kph_t = int(np.argmin(np.abs(table['tau'] - 2.0/3.0)))
    row('photosphere rho [g/cm^3]', ideal['rho'][kph_i], table['rho'][kph_t])
    row('photosphere mu', ideal['mu'][kph_i], table['mu'][kph_t])
    row('photosphere Gamma_1', ideal['gamma1'][kph_i], table['gamma1'][kph_t])
    row('photosphere grad_ad', ideal['grad_ad'][kph_i], table['grad_ad'][kph_t])
    w()
    row('min grad_ad in box', ideal['grad_ad'].min(), table['grad_ad'].min())
    row('min Gamma_1 in box', ideal['gamma1'].min(), table['gamma1'].min())
    row('max mu in box', ideal['mu'].max(), table['mu'].max())
    w()
    sup_i = int(np.sum(ideal['grad'] - ideal['grad_ad'] > 0))
    sup_t = int(np.sum(table['grad'] - table['grad_ad'] > 0))
    row('superadiabatic cells (of %d)' % len(ideal['z']), sup_i, sup_t, '%14d')
    row('CZ top z/zmax (last unstable)',
        ideal['z'][np.max(np.where(ideal['grad'] > ideal['grad_ad'])[0])]/zmax
        if sup_i else np.nan,
        table['z'][np.max(np.where(table['grad'] > table['grad_ad'])[0])]/zmax
        if sup_t else np.nan)
    w()
    row('v_z,rms at base [km/s]', ideal['vx_rms'][0]/1e5, table['vx_rms'][0]/1e5)
    row('v_z,rms at photosphere [km/s]',
        ideal['vx_rms'][kph_i]/1e5, table['vx_rms'][kph_t]/1e5)
    row('v_z,rms at top [km/s]', ideal['vx_rms'][-1]/1e5, table['vx_rms'][-1]/1e5)
    row('peak Mach', ideal['mach'].max(), table['mach'].max())
    w()
    for name, r, _ in runs:
        h = read_hst(name)
        w('%-9s history: mass %+.2f%%, totE %+.2f%%, final KE %.3e, t_end %.0f'
          % (name, 100*(h['mass'][-1]/h['mass'][0] - 1),
             100*(h['etot'][-1]/h['etot'][0] - 1), h['ke'][-1], h['t'][-1]))
    w()
    w('plots in %s (suffix %r)' % (PLOTS, tag))
    with open(os.path.join(ROOT, 'summary%s.txt' % tag), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
