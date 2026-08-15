#!/usr/bin/env python3
"""Sample T and vertical velocity on the CORRUGATED tau=2/3 surface (the real optical
photosphere) instead of a flat geometric height -> the true granulation appearance.

    python3 tau_surface.py <run> [first] [last]        e.g.  python3 tau_surface.py ideal

Adapted from run/plot_tau_surface.py, with three changes:
  * reads soleos/<run>/bin directly (no .athdf round trip written next to the data),
  * gets T from the run's OWN eos via eoslib, so it is correct for general_eos=table too
    (the original hardcoded the ideal-gas T = eint*(gamma-1)/(rho Rgas)),
  * writes to soleos/plots/tau_surface/<run>/.

The tau=2/3 crossing is linearly interpolated in tau between the bracketing cells, which
removes the stair-step artifact of nearest-cell sampling. Opacity is eoslib.get_kapr,
which mirrors get_kapr() in src/pgen/solar_convection.cpp exactly.
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

# colour limits, kept at the values run/plot_tau_surface.py used so the ideal maps look
# the same as before. Override on the command line with VMIN/VMAX env vars if a run
# lands somewhere else.
TLIM = (float(os.environ.get('TMIN', 5200.0)), float(os.environ.get('TMAX', 5700.0)))
VLIM = (float(os.environ.get('VMIN', -6.0)), float(os.environ.get('VMAX', 6.0)))


def load(run, idx, files=None):
    """One snapshot on the global (nx3, nx2, nx1) grid, plus face coordinates."""
    if files is None:
        files = dumps.dump_files(os.path.join(ROOT, run))
    fd = bin_convert.read_binary(files[idx])
    nx1, nx2, nx3 = fd['Nx1'], fd['Nx2'], fd['Nx3']
    m1, m2, m3 = fd['nx1_mb'], fd['nx2_mb'], fd['nx3_mb']
    out = {}
    for var in ['dens', 'velx', 'eint']:
        g = np.empty((nx3, nx2, nx1), dtype=np.float64)
        blocks = np.asarray(fd['mb_data'][var])
        for b in range(fd['n_mbs']):
            i, j, k = fd['mb_logical'][b][:3]
            g[k*m3:(k+1)*m3, j*m2:(j+1)*m2, i*m1:(i+1)*m1] = blocks[b]
        out[var] = g
    out['dz'] = (fd['x1max'] - fd['x1min'])/nx1
    out['x2f'] = np.linspace(fd['x2min'], fd['x2max'], nx2 + 1)
    out['x3f'] = np.linspace(fd['x3min'], fd['x3max'], nx3 + 1)
    out['time'] = fd['time']
    out['n'] = len(files)
    return out


def tau_surface(d, eos):
    """T and v_z interpolated onto tau = 2/3, counting tau down from the top."""
    rho, e, vz = d['dens'], d['eint'], d['velx']/1e5
    T = eos.temperature(rho, e/rho)
    dtau = eoslib.get_kapr(rho, T)*rho*d['dz']
    nx1 = rho.shape[2]
    tau = np.cumsum(dtau[:, :, ::-1], axis=2)[:, :, ::-1]   # from the top downward
    target = 2.0/3.0
    idx = np.argmax(tau[:, :, ::-1] >= target, axis=2)      # index counted from the top
    iph = np.clip((nx1 - 1) - idx, 0, nx1 - 2)              # keep iph+1 in range
    K, J = np.meshgrid(np.arange(rho.shape[0]), np.arange(rho.shape[1]), indexing='ij')
    tb = tau[K, J, iph]                                     # >= 2/3, below the crossing
    ta = tau[K, J, iph+1]                                   # <  2/3, above the crossing
    f = np.clip((target - ta)/np.maximum(tb - ta, 1e-30), 0.0, 1.0)
    Tsurf = (1.0-f)*T[K, J, iph+1] + f*T[K, J, iph]
    Vsurf = (1.0-f)*vz[K, J, iph+1] + f*vz[K, J, iph]
    # flag columns where tau never reaches 2/3 (optically thin all the way down)
    thin = tau[:, :, 0] < target
    return Tsurf, Vsurf, thin


def main():
    run = sys.argv[1] if len(sys.argv) > 1 else 'ideal'
    eos = eoslib.TableEOS(GRID) if run.startswith('table') else eoslib.IdealEOS()
    # TAG appends to the output directory name, so the same run can be rendered twice at
    # different colour limits without either set overwriting the other.
    outdir = os.path.join(ROOT, 'plots', 'tau_surface', run + os.environ.get('TAG', ''))
    os.makedirs(outdir, exist_ok=True)
    files = dumps.dump_files(os.path.join(ROOT, run))
    ntot = len(files)
    first = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    last = int(sys.argv[3]) if len(sys.argv) > 3 else ntot
    for i in range(first, min(last, ntot)):
        if files[i] is None:
            print('%05d  SKIPPED: truncated dump' % i)
            continue
        d = load(run, i, files)
        Tsurf, Vsurf, thin = tau_surface(d, eos)
        for op, fld, lim, cmap, lbl in [
                ('vz', Vsurf, VLIM, 'seismic_r', 'Vertical velocity [km s$^{-1}$]'),
                ('tmp', Tsurf, TLIM, 'gist_heat', 'Temperature [K]')]:
            fig = plt.figure(figsize=(6.2, 5.2))
            im = plt.pcolormesh(d['x2f'], d['x3f'], fld, vmin=lim[0], vmax=lim[1],
                                cmap=cmap, shading='flat')
            plt.colorbar(im, label=lbl, pad=0.02)
            plt.xlabel('y [cm]')
            plt.ylabel('z [cm]')
            plt.gca().set_aspect('equal')
            plt.title('%s   t = %.0f s   ($\\tau=2/3$ surface)' % (run, d['time']))
            plt.tight_layout()
            fig.savefig(os.path.join(outdir, '%s_%05d.png' % (op, i)), dpi=130,
                        bbox_inches='tight')
            plt.close()
        print('%05d  t=%8.0f  T[%.0f, %.0f] K  vz[%+.2f, %+.2f] km/s  thin columns %d'
              % (i, d['time'], Tsurf.min(), Tsurf.max(), Vsurf.min(), Vsurf.max(),
                 thin.sum()))
    print('saved tau=2/3 maps %d..%d in %s' % (first, min(last, ntot)-1, outdir))


if __name__ == '__main__':
    main()
