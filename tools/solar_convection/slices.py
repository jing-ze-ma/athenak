#!/usr/bin/env python3
"""Vertical (x3-x1) 2D slices through the mid-x2 plane: vz, T, dtau, Mach.

    python3 slices.py <run> [first] [last]        e.g.  python3 slices.py ideal

Adapted from run/plotsun.py, with three changes:
  * reads soleos/<run>/bin directly (no .athdf round trip written next to the data),
  * T, p and c_s come from the run's OWN eos via eoslib, so this is valid for
    general_eos=table too (the original hardcoded p = 0.66667*eint, T = p/(rho Rgas),
    c_s from gamma=5/3),
  * the opacity is eoslib.get_kapr, which mirrors get_kapr() in solar_convection.cpp
    (C=20, composite H-/Kramers/e-scattering, 1e-2 floor). plotsun.py's tau panel used
    a STALE opacity (C=1e2, no floor, no Kramers) and was not trustworthy.

Colour limits are kept at plotsun.py's values so these are directly comparable with the
earlier solar-convection figures. Override with VZLIM / TMAX / MACHMAX env vars.
"""

import os
import sys

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt                       # noqa: E402

_HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(_HERE, os.pardir, os.pardir, 'vis', 'python'))
sys.path.insert(0, _HERE)
import bin_convert                                    # noqa: E402
import dumps                                          # noqa: E402
import eoslib                                         # noqa: E402

ROOT = os.environ.get('SOLEOS_ROOT', '/orion/u/jinma/ATHENAK/soleos')
GRID = os.path.join(ROOT, 'tools', 'eos_grid.bin')
STYLE = os.environ.get('MPLSTYLE',
                       os.path.join(_HERE, os.pardir, os.pardir, 'run', 'mesa2.mplstyle'))

VZLIM = float(os.environ.get('VZLIM', 2.0))
TMAX = float(os.environ.get('TMAX', 6000.0))
MACHMAX = float(os.environ.get('MACHMAX', 1.5))


def load(run, idx, files=None):
    """One snapshot on the global (nx3, nx2, nx1) grid, plus face coordinates."""
    if files is None:
        files = dumps.dump_files(os.path.join(ROOT, run))
    fd = bin_convert.read_binary(files[idx])
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
    out['x1f'] = np.linspace(fd['x1min'], fd['x1max'], nx1 + 1)
    out['x3f'] = np.linspace(fd['x3min'], fd['x3max'], nx3 + 1)
    out['dz'] = (fd['x1max'] - fd['x1min'])/nx1
    out['time'] = fd['time']
    out['n'] = len(files)
    return out


def fields(d, eos):
    """The four plotted quantities, all through the run's own EOS."""
    rho, e = d['dens'], d['eint']
    espec = e/rho
    T = eos.temperature(rho, espec)
    p = eos.pressure(rho, espec)
    cs = np.sqrt(eos.gamma1(rho, espec)*p/rho)
    v = np.sqrt(d['velx']**2 + d['vely']**2 + d['velz']**2)
    dtau = eoslib.get_kapr(rho, T)*rho*d['dz']
    return {'vz': d['velx']/1e5, 'tmp': T, 'tau': np.log10(dtau), 'mach': v/cs}


PANELS = [
    ('vz',   'seismic', -VZLIM, VZLIM, 'Vertical velocity$\\rm\\; [km\\; s^{-1}]$'),
    ('tmp',  'jet',     0.0,    TMAX,  '$\\rm Temperature\\; [K]$'),
    ('tau',  'jet',    -3.0,    3.0,   '$\\rm \\log_{10}\\, \\Delta\\tau$'),
    ('mach', 'inferno', 0.0,    MACHMAX, '$\\rm Mach\\; number$'),
]


def main():
    run = sys.argv[1] if len(sys.argv) > 1 else 'ideal'
    eos = eoslib.eos_for_run(os.path.join(ROOT, run), GRID)
    if os.path.exists(STYLE):
        plt.style.use(STYLE)
    color, bg = 'white', 'black'
    for k in ['text.color', 'axes.labelcolor', 'axes.edgecolor',
              'xtick.color', 'ytick.color']:
        mpl.rcParams[k] = color

    files = dumps.dump_files(os.path.join(ROOT, run))
    ntot = len(files)
    first = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    last = int(sys.argv[3]) if len(sys.argv) > 3 else ntot
    for op, _, _, _, _ in PANELS:
        os.makedirs(os.path.join(ROOT, 'plots', 'slices', run, op), exist_ok=True)

    for i in range(first, min(last, ntot)):
        if files[i] is None:
            print('%05d  SKIPPED: truncated dump' % i)
            continue
        d = load(run, i, files)
        f = fields(d, eos)
        x3c = 0.5*(d['x3f'][1:] + d['x3f'][:-1])
        x1c = 0.5*(d['x1f'][1:] + d['x1f'][:-1])
        px, py = np.meshgrid(x3c, x1c, indexing='ij')
        jmid = d['dens'].shape[1]//2
        for op, cmap, vmin, vmax, label in PANELS:
            fig = plt.figure(figsize=(8, 5), facecolor=bg)
            ax = plt.subplot()
            ax.set_facecolor(bg)
            im = plt.pcolormesh(px, py, f[op][:, jmid], vmin=vmin, vmax=vmax,
                                cmap=cmap, shading='nearest')
            cbar = plt.colorbar(im, pad=0.03)
            cbar.set_label(label, fontsize=15)
            ax.set_xlabel('y', fontsize=15)
            ax.set_ylabel('z', fontsize=15)
            ax.tick_params(axis='both', which='both', direction='out')
            plt.title('$\\rm t$ = %.1f' % d['time'], y=0.9, color=color)
            fig.savefig(os.path.join(ROOT, 'plots', 'slices', run, op,
                                     '%05d.png' % i),
                        facecolor=bg, transparent=False, bbox_inches='tight')
            plt.close()
        print('%05d  t=%8.0f  T[%.0f,%.0f]K  vz[%+.2f,%+.2f]km/s  Mach max %.2f'
              % (i, d['time'], f['tmp'].min(), f['tmp'].max(),
                 f['vz'].min(), f['vz'].max(), f['mach'].max()))
    print('saved slices %d..%d in %s'
          % (first, min(last, ntot)-1, os.path.join(ROOT, 'plots', 'slices', run)))


if __name__ == '__main__':
    main()
