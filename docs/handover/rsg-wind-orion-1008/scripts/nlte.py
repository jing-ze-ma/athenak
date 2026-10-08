"""Non-LTE (band-dependent thermalisation) correction to the optically thin radiative
equilibrium (RE) of RSG wind gas, after twophase.py Part 1.

Physics (two-level line ensemble per band b, thin gas in the diluted stellar field):
  net heating / gram / (4 pi)
    H(T) = sum_b [eps_b(T) kl_b + kc_b] [W B_b(Teff) - B_b(T)]
  kl_b = sum_g w_g k_line,bg  (Exo-FMS premixed table: molecular + atomic LINES only)
  kc_b = CIA + H- bf/ff (ck_lib.continuum minus Rayleigh) -> eps = 1 (true absorption)
  Rayleigh = pure scattering, removed (as in twophase.py).
  eps = C/(C + A), C = n_coll q;  n_coll = n_H2 + n_H + n_He (H2 is dissociated in the
  warm phase, 1e-7 of particles at 2500 K / 1e-10 bar; the n_H2-only variant is a flag).
  bands 0-7 (324-0.85 um): vib-rot, A = A_IR;  bands 8-10 (0.85-0.26 um): electronic
  molecular bands + atomic lines, A = A_el.
  Electronic-band cascade (fluorescence through vibrational levels):
    'a'   eps_el = C/(C+A_el)                     (pure scattering limit)
    'b'   eps_el = C/(C+A_el) + f eps_IR          (fraction f of the absorbed energy left
                                                    as vibrational excitation, thermalised
                                                    like an IR band, rest re-emitted)
    'bp'  eps_el = 1 - f                          (optimistic: all but f heats)
    'lte' eps = 1 everywhere (= twophase.py Part 1)
  The emission term carries the same eps (detailed balance: H = 0 at J = B).
Thermal time t_th = c_p / |4 pi dH/dT|_P (isobaric; c_v and isochoric for fixed rho),
  c_p per particle 5/2 k (atoms) or 7/2 k (H2), no dissociation latent heat (-> lower bound).
Run (one process per table):  TAB=lowP|old python nlte.py   -> res_<TAB>.pkl
"""
import os
import pickle
import sys
import numpy as np
sys.dont_write_bytecode = True
TAB = os.environ.get('TAB', 'lowP')
os.environ['CK_DATA'] = '/orion/ptmp/jinma/rsg_wind_1008/lowP/data/'
if TAB == 'lowP':
    os.environ['CK_KTABLE'] = 'Premixed_1x_g8_11_hiT2_lowP.txt'
    os.environ['CK_CETABLE'] = 'FastChem_ck_1x_int_hiT2_lowP.txt'
else:
    os.environ['CK_KTABLE'] = 'Premixed_1x_g8_11_hiT2.txt'
    os.environ['CK_CETABLE'] = 'FastChem_ck_1x_int_hiT2.txt'
import rsglib as L   # noqa: E402
ck = L.ck
OUT = os.path.dirname(os.path.abspath(__file__)) + '/'

TG = np.geomspace(300., 4000., 321)
PBAR = np.array([1e-11, 1e-10, 1e-9])
RHO = np.array([1e-16, 1e-15, 1e-14])
XR = np.array([2., 3., 4., 5.])
STARS = ('golden16', 'm20lgl5.5')
IEL = np.array([8, 9, 10])            # 0.85-0.26 um
QS = (1e-12, 1e-11, 1e-10)
AIR = (10., 100.)
AEL = (1e6, 1e8)
FC = (0.1, 0.5)


def dilution(x):
    return 0.5*(1 - np.sqrt(np.maximum(1 - 1/np.asarray(x, float)**2, 0)))


def Bb(T):
    T = np.atleast_1d(T)
    return L.SIG*T[:, None]**4*L.fband(T)/np.pi


def point(T, P=None, rho=None):
    """line kl (nb), absorptive continuum kc (nb), mu, rho, n_coll, n_H2, c_p, c_v (per g)."""
    mu = 2.3
    for _ in range(60):
        if rho is not None:
            Pd = rho*L.KB*T/(mu*L.MH)
        else:
            Pd = P*1e6
        v = ck.comp(L.CE, T, Pd*1e-6)
        if abs(v[0] - mu) < 1e-8*mu:
            break
        mu = 0.5*(mu + v[0])
    Pd = rho*L.KB*T/(mu*L.MH) if rho is not None else P*1e6
    rr = Pd*mu*L.MH/(L.KB*T)
    kc, rho_c, v = ck.continuum(L.CE, L.CIA, L.RAY, L.WL, T, Pd*1e-6)
    ntot = Pd/(L.KB*T)
    kray = sum(L.RAY[s]*v[s+1] for s in range(4))*ntot/rho_c
    kc = (kc - kray)*rho_c/rr
    kl = (L.GW*L.linek(T, Pd*1e-6, 'clamp')).sum(-1)
    xh2 = v[1]
    cpp = 2.5 + xh2           # per particle /k
    return dict(kl=kl, kc=kc, mu=v[0], rho=rr, P=Pd, ncoll=ntot*(v[1] + v[2] + v[3]),
                nh2=ntot*v[1], cp=cpp*L.KB/(v[0]*L.MH), cv=(cpp - 1)*L.KB/(v[0]*L.MH))


def tables():
    fn = OUT + f'grid_{TAB}.npz'
    if os.path.exists(fn):
        return dict(np.load(fn))
    d = {}
    for key, vals in (('P', PBAR), ('rho', RHO)):
        arr = {k: [] for k in ('kl', 'kc', 'mu', 'rho', 'P', 'ncoll', 'nh2', 'cp', 'cv')}
        for val in vals:
            row = [point(T, P=val) if key == 'P' else point(T, rho=val) for T in TG]
            for k in arr:
                arr[k].append([r[k] for r in row])
        for k in arr:
            d[key + '_' + k] = np.array(arr[k])
    np.savez(fn, **d)
    return d


def eps_bands(n, q, air, ael, casc, f):
    """(nT, nb) thermalisation probability."""
    C = n*q
    eir = C/(C + air)
    eel = C/(C + ael)
    if casc == 'b':
        eel = np.minimum(eel + f*eir, 1.)
    elif casc == 'bp':
        eel = np.full_like(eel, 1 - f)
    e = np.repeat(eir[:, None], L.NB, 1)
    e[:, IEL] = eel[:, None]
    return e


def heat(kl, kc, eps, W, Teff):
    """per gram / (4 pi), (nT,); also the separate line and continuum parts."""
    S = W*Bb(Teff) - Bb(TG)
    hl = (eps*kl*S).sum(-1)
    hc = (kc*S).sum(-1)
    return hl + hc, hl, hc


def roots(H):
    out = []
    lT = np.log(TG)
    dH = np.gradient(H, TG)
    for i in np.where(np.sign(H[:-1]) != np.sign(H[1:]))[0]:
        a = H[i]/(H[i] - H[i+1])
        T = float(np.exp(lT[i] + a*(lT[i+1] - lT[i])))
        out.append((T, bool(H[i+1] < H[i]), float((1-a)*dH[i] + a*dH[i+1]), i, a))
    return out


def scenarios():
    yield ('lte', None, None, None, 'lte', None)
    for q in QS:
        for air in AIR:
            for ael in AEL:
                yield (f'a q{q:.0e} AIR{air:.0f} Ael{ael:.0e}', q, air, ael, 'a', None)
                for f in FC:
                    yield (f'b{f} q{q:.0e} AIR{air:.0f} Ael{ael:.0e}', q, air, ael, 'b', f)
                    yield (f'bp{f} q{q:.0e} AIR{air:.0f} Ael{ael:.0e}', q, air, ael,
                           'bp', f)
    # n_H2-only collider flag (q = 1e-10, A_IR = 10, A_el = 1e6, scattering electronic)
    yield ('H2only a q1e-10 AIR10 Ael1e+06', 1e-10, 10., 1e6, 'a_h2', None)


def run():
    d = tables()
    res = []
    for sname in STARS:
        st = L.star(sname)
        for x in XR:
            W = dilution(x)
            for key, vals in (('P', PBAR), ('rho', RHO)):
                for iv, val in enumerate(vals):
                    kl, kc = d[key + '_kl'][iv], d[key + '_kc'][iv]
                    n = d[key + '_ncoll'][iv]
                    cap = d[key + '_cp'][iv] if key == 'P' else d[key + '_cv'][iv]
                    for name, q, air, ael, casc, f in scenarios():
                        if casc == 'lte':
                            eps = np.ones((len(TG), L.NB))
                        elif casc == 'a_h2':
                            eps = eps_bands(d[key + '_nh2'][iv], q, air, ael, 'a', f)
                        else:
                            eps = eps_bands(n, q, air, ael, casc, f)
                        H, hl, hc = heat(kl, kc, eps, W, st['Teff'])
                        rl = []
                        for T, stab, dHdT, i, a in roots(H):
                            c = (1-a)*cap[i] + a*cap[i+1]
                            tth = c/abs(4*np.pi*dHdT) if dHdT != 0 else np.inf
                            # fraction of the absorbed power that is continuum, at the root
                            S = W*Bb(st['Teff'])[0]
                            ab_l = ((1-a)*(eps[i]*kl[i]*S).sum() + a*(eps[i+1]*kl[i+1]*S).sum())
                            ab_c = (1-a)*(kc[i]*S).sum() + a*(kc[i+1]*S).sum()
                            rl.append(dict(T=T, stable=stab, tth=tth,
                                           fcont=ab_c/(ab_c + ab_l),
                                           rho=float(np.exp((1-a)*np.log(d[key + '_rho'][iv][i])
                                                            + a*np.log(d[key + '_rho'][iv][i+1]))),
                                           eir=float(eps[i, 0]), eel=float(eps[i, 9])))
                        res.append(dict(star=sname, x=x, key=key, val=val, scen=name,
                                        casc=casc, q=q, air=air, ael=ael, f=f, roots=rl,
                                        hl_edge=(float(H[0]), float(H[-1]))))
    with open(OUT + f'res_{TAB}.pkl', 'wb') as fh:
        pickle.dump(res, fh)
    return res


if __name__ == '__main__':
    run()
