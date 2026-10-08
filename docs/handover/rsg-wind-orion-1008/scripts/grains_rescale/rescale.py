"""Rescaled FT chromosphere: (i) f_rho on rho(r) and the launch flux, (ii) v_esc/v_con at
fixed rho_ph. Reuses ../grainlib.py unchanged. usage: ../venv/bin/python rescale.py [nproc]"""
import sys
import signal
import warnings
import pickle
import itertools
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, '/orion/ptmp/jinma/rsg_wind_1008/grains')
import grainlib as gl  # noqa

O = '/orion/ptmp/jinma/rsg_wind_1008/grains/rescale/out/'
V0 = list(np.arange(20., 100.01, 2.5))               # km/s
SCAL = [('f', 1.0), ('f', 0.3), ('f', 0.1), ('f', 0.03), ('f', 0.01),
        ('q', 12.5), ('q', 15.), ('q', 20.)]           # ('q', MESA) == ('f', 1)
OPT = ['lowk', 'pure', 'Fe3e-4', 'Fe1e-3']
ALPHA = [0.1, 1.0]
XSEED = [1e-15, 1e-13]
GGAS = [0.0, 0.1]
TGAS = [400., 1200.]
_m = {}


def model(star, var, sc):
    k = (star, var, sc)
    if k not in _m:
        m = gl.Model(star, var)
        st = m.st
        if sc[0] == 'f':
            st.rho_ph *= sc[1]          # scales rho(r) and mdot0 together
        else:
            st.v_con = st.vesc / sc[1]
        st.mdot0 = 4 * np.pi * st.R**2 * st.rho_ph * st.v_con
        _m[k] = m
    return _m[k]


class _TO(Exception):
    pass


def _alarm(sig, frm):
    raise _TO()


def one(p):
    warnings.simplefilter("ignore")
    star, sc, v0, var, al, xs, gg, tg, dens = p
    m = model(star, var, sc)
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(30)                     # per-run guard (one LSODA run hung)
    try:
        r = m.run(v0 * 1e5, al, xs, tg, gg, dens=dens)
    except _TO:
        return dict(star=star, sc=sc, v0=v0, var=var, alpha=al, xseed=xs, Ggas=gg, Tg=tg,
                    dens=dens, outcome='timeout', vinf=np.nan,
                    rmax0=gl.rmax_dustfree(m.st, v0 * 1e5, gg), vesc=m.st.vesc,
                    v_con=m.st.v_con, mdot0=m.st.mdot0)
    finally:
        signal.alarm(0)
    grow = np.where(r['a'] > 1.05 * gl.A_SEED)[0]
    i1 = np.where(r['a'] > 1e-5)[0]
    i2 = np.where(r['gd'] + gg > 1)[0]
    return dict(star=star, sc=sc, v0=v0, var=var, alpha=al, xseed=xs, Ggas=gg, Tg=tg,
                dens=dens, outcome=r['outcome'], vinf=r.get('vinf', np.nan),
                rmax=r['rmax'], rmax0=gl.rmax_dustfree(m.st, v0 * 1e5, gg),
                amax=r['amax'], gdmax=r['gdmax'], vdmax=r['vd'].max() if len(grow) else 0,
                r_seed=r['r'][grow[0]] if len(grow) else np.nan,
                r_01=r['r'][i1[0]] if len(i1) else np.nan,
                vr_01=r['v'][i1[0]] if len(i1) else np.nan,
                r_G1=r['r'][i2[0]] if len(i2) else np.nan,
                v_G1=r['v'][i2[0]] if len(i2) else np.nan,
                vesc=m.st.vesc, v_con=m.st.v_con, mdot0=m.st.mdot0)


if __name__ == '__main__':
    npr = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    P = [p + ('ft',) for p in itertools.product(list(gl.STARS), SCAL, V0, OPT, ALPHA, XSEED,
                                                 GGAS, TGAS)]
    P += [p + ('mc',) for p in itertools.product(list(gl.STARS), SCAL, V0, ['pure'], ALPHA,
                                                  XSEED, GGAS, TGAS)]
    with Pool(npr) as pool:
        res = list(pool.imap_unordered(one, P, chunksize=8))
    pickle.dump(res, open(O + 'grid.pkl', 'wb'))
    print('done', len(res))
