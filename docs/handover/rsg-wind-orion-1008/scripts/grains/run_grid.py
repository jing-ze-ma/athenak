"""Full parameter grid of parcel trajectories -> out/grid.pkl (summaries) and
out/traj_rep.pkl (representative trajectories). usage: venv/bin/python run_grid.py [nproc]"""
import sys
import pickle
import itertools
import numpy as np
from multiprocessing import Pool
import grainlib as gl

V0 = [20, 30, 40, 50, 60, 70, 80, 85, 90]          # km/s (20-50 briefed; 60-90 added)
ALPHA = [0.01, 0.1, 1.0]
XSEED = [1e-16, 1e-15, 1e-14, 1e-13]
VAR = ['lowk', 'pure', 'Fe3e-4', 'Fe1e-3']
TGAS = [400., 800., 1200., 1700.]
GGAS = [0.0, 0.1, 0.3]
DENS = ['ft', 'mc']
_models = {}


def model(star, var):
    k = (star, var)
    if k not in _models:
        _models[k] = gl.Model(star, var)
    return _models[k]


def one(p):
    star, v0, al, xs, var, tg, gg, dens = p
    m = model(star, var)
    r = m.run(v0 * 1e5, al, xs, tg, gg, dens=dens)
    grow = np.where(r['a'] > 1.05 * gl.A_SEED)[0]
    s = dict(star=star, v0=v0, alpha=al, xseed=xs, var=var, Tg=tg, Ggas=gg, dens=dens,
             outcome=r['outcome'], rmax=r['rmax'], t3=r['t3'], amax=r['amax'],
             fmax=r['fmax'], gdmax=r['gdmax'], vinf=r.get('vinf', np.nan),
             mdot=r.get('mdot', np.nan), vdmax=r['vd'].max() if len(grow) else 0.,
             r_seed=r['r'][grow[0]] if len(grow) else np.nan,
             rmax0=gl.rmax_dustfree(m.st, v0 * 1e5, gg), tend=r['t'][-1])
    i1 = np.where(r['a'] > 1e-5)[0]
    s['t_01um'] = (r['t'][i1[0]] - r['t'][grow[0]]) if len(i1) else np.nan
    i2 = np.where(r['gd'] + gg > 1)[0]
    s['r_G1'] = r['r'][i2[0]] if len(i2) else np.nan
    return s


if __name__ == '__main__':
    npr = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    P = list(itertools.product(list(gl.STARS), V0, ALPHA, XSEED, VAR, TGAS, GGAS, DENS))
    with Pool(npr) as pool:
        res = pool.map(one, P, chunksize=64)
    pickle.dump(res, open(gl.D + 'grid.pkl', 'wb'))
    print('done', len(res))
