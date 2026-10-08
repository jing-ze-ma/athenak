"""Dust along (a) the steady-wind trajectory through the observed density and (b) FT-type
ballistic parcels whose density is the observed profile.  usage: python dust.py [nproc]"""
import sys
import signal
import pickle
import warnings
import itertools
import numpy as np
from multiprocessing import Pool
import grainlib_bet as gl
import prof as P

O = '/orion/ptmp/jinma/rsg_wind_1008/betelgeuse/out/'
STARS = ['betA18', 'betB16.5']
OPT = ['pure', 'Fe3e-4', 'Fe1e-3']
ALPHA = [0.1, 1.0]
XSEED = [1e-15, 1e-13]
MDOT = [1e-6, 2e-6, 4e-6]
V0 = list(np.arange(20., 100.01, 2.5))
GGAS = [0.0, 0.1]
VINF = 12e5
_m = {}


def model(star, var):
    if (star, var) not in _m:
        _m[(star, var)] = gl.Model(star, var)
    return _m[(star, var)]


class _TO(Exception):
    pass


def _alarm(sig, frm):
    raise _TO()


def steady_one(p):
    warnings.simplefilter('ignore')
    star, md, var, al, xs = p
    m = model(star, var)
    st = m.st
    Md = md * gl.Msun / gl.YR

    def vfun(r):
        return min(Md / (4 * np.pi * r**2 * st.rho_obs(r)), VINF)
    res = m.run_steady(vfun, al, xs, st.T_obs)
    res.update(star=star, mdot=md, var=var, alpha=al, xseed=xs)
    return res


def ball_one(p):
    warnings.simplefilter('ignore')
    star, v0, var, al, xs, gg = p
    m = model(star, var)
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(60)
    base = dict(star=star, v0=v0, var=var, alpha=al, xseed=xs, Ggas=gg,
                rmax0=gl.rmax_dustfree(m.st, v0 * 1e5, gg), vesc=m.st.vesc,
                mdot_gt=m.st.mdot_gt(v0 * 1e5) * gl.YR / gl.Msun)
    try:
        r = m.run(v0 * 1e5, al, xs, m.st.T_obs, gg, dens='obs', rout=40.)
    except _TO:
        base.update(outcome='timeout')
        return base
    finally:
        signal.alarm(0)
    grow = np.where(r['a'] > 1.05 * gl.A_SEED)[0]
    i1 = np.where(r['a'] > 1e-5)[0]
    base.update(outcome=r['outcome'], vinf=r.get('vinf', np.nan), rmax=r['rmax'],
                amax=r['amax'], gdmax=r['gdmax'],
                r_seed=r['r'][grow[0]] if len(grow) else np.nan,
                r_01=r['r'][i1[0]] if len(i1) else np.nan, t_end=r['t'][-1])
    return base


if __name__ == '__main__':
    npr = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    PS = list(itertools.product(STARS, MDOT, OPT, ALPHA, XSEED))
    PB = list(itertools.product(STARS, V0, OPT, ALPHA, XSEED, GGAS))
    with Pool(npr) as pool:
        RS = pool.map(steady_one, PS, chunksize=1)
        pickle.dump(RS, open(O + 'dust_steady.pkl', 'wb'))
        print('steady done', len(RS), flush=True)
        RB = pool.map(ball_one, PB, chunksize=1)
    pickle.dump(RB, open(O + 'dust_ball.pkl', 'wb'))
    print('ballistic done', len(RB), flush=True)
    for s in STARS:
        st = gl.Star(s)
        print(s, 'FT fit rho_ph %.3g v_con %.2f km/s vesc %.1f' % (st.rho_ph, st.v_con/1e5,
                                                                  st.vesc/1e5),
              P.ft_fit(P.stellar(s[4], float(s[5:]))))
