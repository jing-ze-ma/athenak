"""Recompute the Dent 2024 (set A, M 18, Mdot 2e-6, v_inf cap 12 km/s) Gamma_req, v_t,eq and
dust seed radii with the ORIGINAL betelgeuse/ modules (read-only import, no bytecode written)
-> out/old24.npz.  Gamma_avail,gas = max{static RE, Sobolev t10 A, thin 1500 K} as in the
original analyse.py."""
import sys
import pickle
sys.dont_write_bytecode = True
B = '/orion/ptmp/jinma/rsg_wind_1008/betelgeuse/'
sys.path.insert(0, B)
import numpy as np  # noqa: E402
import prof as P  # noqa: E402
import wind as W  # noqa: E402
assert P.__file__.startswith(B)
st = P.stellar('A', 18.)
w = W.steady(st, 'dent', 2e-6)
g = {k: np.load(B + f'out/gas_A_A18_dent_{k}.npz') for k in ('re', '1500', 'obs')}
x = w['x']


def at(xx, yy):
    return np.interp(np.log(x), np.log(xx), yy, left=np.nan, right=np.nan)


Gst = at(g['re']['cr'], g['re']['cGF'])
Gso = at(g['obs']['cr'], g['obs']['cGF'])
Gs10 = at(g['obs']['sx'], g['obs']['Gsob_t10_A'])
Gth15 = at(g['1500']['sx'], g['1500']['Gthin'])
Gth = at(g['obs']['sx'], g['obs']['Gthin'])
Gmax = np.nanmax(np.vstack([Gst, Gs10, Gth15]), axis=0)
vteq = W.vt_equiv(w, Gmax)
vt0 = W.vt_equiv(w, 0.)
RS = pickle.load(open(B + 'out/dust_steady.pkl', 'rb'))
seeds = []
for q in RS:
    if q['star'] != 'betA18' or q['mdot'] != 2e-6:
        continue
    gr = np.where(q['a'] > 1.05e-7)[0]
    i1 = np.where(q['a'] > 1e-5)[0]
    seeds.append((q['var'], q['alpha'], q['xseed'], q['r'][gr[0]] if len(gr) else np.nan,
                  q['r'][i1[0]] if len(i1) else np.nan))
np.savez('/orion/ptmp/jinma/rsg_wind_1008/betelgeuse26/out/old24.npz', x=x, r=w['r'],
         rho=w['rho'], T=w['T'], v=w['v'], G0=w['Greq'][0.], G5=w['Greq'][5e5],
         G10=w['Greq'][10e5], Gmax=Gmax, Gst=Gst, Gso=Gso, Gs10=Gs10, Gth15=Gth15, Gth=Gth,
         vteq=vteq, vt0=vt0, R=st['R'], M=st['M'],
         seeds=np.array([(s[3], s[4]) for s in seeds]),
         seedlab=np.array([f'{s[0]}/a{s[1]:g}/{s[2]:.0e}' for s in seeds]))
print('done', st['R'], len(seeds))
