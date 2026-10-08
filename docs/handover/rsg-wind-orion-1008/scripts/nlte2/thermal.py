"""Stage 2: net heating H(T) at fixed P = atoms (stage 1, atoms_res.pkl) + molecules
(NK93 cooling with a consistent two-level radiative-pumping factor) + IR vib-rot bands
0.85-4.4 um (k table x eps_vib) + electronic bands 0.26-0.85 um (k table x f_heat TiO)
+ continuum (H-/CIA, eps = 1, ck_lib.continuum via nlte.point); roots, stability,
t_th.     python thermal.py  -> thermal_res.pkl, tables printed by report.py

Two-level net collisional heating with a diluted field (occupation nbar) and Sobolev
escape beta (core penetration = beta):
   H = -Lambda_NK * (1 - phi) * (1 + nbar)/(1 + nbar s),  phi = nbar e^{dE/kT}/(1 + nbar),
   s = A beta/(A beta + C_ul) ~ 1/(1 + n_c/n_half);  Lambda_NK = n_c n_X L(T, n_c, Ntilde).
Exact for one transition; applied with a representative dE per coolant (bracketed).
"""
import os
import sys
import pickle
import numpy as np
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
sys.path.insert(0, HERE)
os.environ['TAB'] = 'lowP'
sys.path.insert(0, '/orion/ptmp/jinma/rsg_wind_1008/nlte')
import nlte as N0    # noqa: E402  (sets CK_* env, loads rsglib)
import mol as M      # noqa: E402
import comp as CP    # noqa: E402
L = N0.L
KB = 1.380649e-16
HCK = 1.438777
KMS = 1e5
BVIB = [3, 4, 5, 6, 7]       # 4.4-0.85 um
BEL = [8, 9, 10]             # 0.85-0.26 um
BFIR = [0, 1, 2]             # 324-4.4 um: replaced by NK93 + SiO
TIO = {'lo': (0.02, 50., 1e-12), 'mid': (0.07, 15., 1e-11), 'hi': (0.3, 5., 1e-10)}
A_SIO, E_SIO, KQ_SIO = 5.4, 1786., 1e-13      # assumed (see report)


def nbar(W, dE, Teff):
    return W/np.expm1(dE/Teff)


def two_level(Lam, W, dE, T, Teff, s):
    nb = nbar(W, dE, Teff)
    phi = nb*np.exp(dE/T)/(1 + nb)
    return -Lam*(1 - phi)*(1 + nb)/(1 + nb*s)


def molecules(T, P, W, Teff, vr, kfac=10., defac=1.):
    """per-volume net heating of each molecular coolant (erg cm^-3 s^-1)."""
    c = CP.comp(T, P, ['H', 'H2', 'He', 'H2O1', 'C1O1', 'O1Si1'])
    nH, nH2 = c['H'], c['H2']
    out = {}
    grad = 3*vr/KMS                              # (dv/dr + 2 v/r) in km/s per cm
    for sp, key in (('H2O', 'H2O1'), ('CO', 'C1O1')):
        nX = c[key]
        Nt = nX/grad
        nc = nH2 + nH                            # rotation: H ~ H2 (assumed)
        Lr, L0, LL, nh = M.rot(sp, T, nc, Nt)
        dE = defac*(min(0.3*T, 600.) if sp == 'H2O' else 2*5.53*np.sqrt(T/5.54))
        out[sp + ' rot'] = two_level(nc*nX*Lr, W, dE, T, Teff, 1/(1 + nc/nh))
        ncv = nH2 + kfac*nH                      # vibration: H quenches kfac x faster
        Lv, L0v, LLv, bet, Aeff = M.vib(sp, T, ncv, Nt)
        out[sp + ' vib'] = two_level(ncv*nX*Lv, W, M.HNU[sp]/KB, T, Teff,
                                     1/(1 + ncv*L0v/LLv))
    # H2 itself (X = H2)
    nc = nH2 + nH
    Lr, L0, LL, nh = M.rot('H2', T, nc, 0.)
    out['H2 rot'] = two_level(nc*nH2*Lr, W, defac*float(np.clip(0.5*T, 510., 2000.)), T, Teff,
                              1/(1 + nc/nh))
    ncv = nH2 + kfac*nH
    Lv, L0v, LLv, _, _ = M.vib('H2', T, ncv, 0.)
    out['H2 vib'] = two_level(ncv*nH2*Lv, W, 5987., T, Teff, 1/(1 + ncv*L0v/LLv))
    # SiO: rotation ~ CO tables with Ntilde x (3.1/0.11)^2; vibration two-level (assumed)
    nX = c['O1Si1']
    Nt = nX/grad*(3.1/0.11)**2
    Lr, L0, LL, nh = M.rot('CO', T, nH2 + nH, Nt)
    out['SiO rot'] = two_level((nH2 + nH)*nX*Lr, W, defac*2*1.04*np.sqrt(T/1.04), T, Teff,
                               1/(1 + (nH2 + nH)/nh))
    kd = KQ_SIO*(T/1000.)**0.5*(nH2 + kfac*nH)
    tau = (8.06e-4)**3*A_SIO/(8*np.pi)*nX*0.05/(vr)         # ~5% in one v=1-0 line
    bet = -np.expm1(-tau)/tau if tau > 1e-6 else 1.
    ku = kd*np.exp(-E_SIO/T)
    Lam = nX*ku*KB*E_SIO*A_SIO*bet/(A_SIO*bet + kd)
    out['SiO vib'] = two_level(Lam, W, E_SIO, T, Teff, A_SIO*bet/(A_SIO*bet + kd))
    # vib eps for the 0.85-4.4 um k-table bands: H2O nu2 quench rate and band A, beta
    Lv, L0v, LLv, bet, Aeff = M.vib('H2O', T, nH2 + kfac*nH, c['H2O1']/grad)
    kdn = L0v*np.exp(M.EV_K['H2O']/T)/M.HNU['H2O']
    Cq = (nH2 + kfac*nH)*kdn
    return out, (Cq, Aeff), c


def Bb(T):
    return L.SIG*T**4*L.fband(np.atleast_1d(T))[0]/np.pi


def grid_k():
    fn = HERE + 'kgrid.npz'
    if os.path.exists(fn):
        return dict(np.load(fn))
    import nlte2 as S1
    d = {k: [] for k in ('kl', 'kc', 'rho', 'cp', 'mu', 'Kg')}
    for P in S1.PBAR:
        rows = [N0.point(T, P=P) for T in S1.TG]
        for r, T in zip(rows, S1.TG):
            r['Kg'] = L.linek(T, P, 'clamp')
        for k in d:
            d[k].append([r[k] for r in rows])
    d = {k: np.array(v) for k, v in d.items()}
    np.savez(fn, **d)
    return d


def roots(TG, H):
    out = []
    lT = np.log(TG)
    for i in np.where(np.sign(H[:-1]) != np.sign(H[1:]))[0]:
        a = H[i]/(H[i] - H[i+1])
        out.append((float(np.exp(lT[i] + a*(lT[i+1] - lT[i]))), bool(H[i+1] < H[i]), i, a))
    return out


def run():
    import nlte2 as S1
    A1 = pickle.load(open(HERE + 'atoms_res.pkl', 'rb'))
    TG = A1['TG']
    kg = grid_k()
    res = []
    for r in A1['res']:
        iP = list(S1.PBAR).index(r['P'])
        W, Teff = r['W'], r['Teff']
        vr = r['dvdr']                       # = v/r
        kl, kc, rho, cp = kg['kl'][iP], kg['kc'][iP], kg['rho'][iP], kg['cp'][iP]
        Kg = kg['Kg'][iP]                                      # (nT, nb, ng)
        vD = np.sqrt(2*KB*TG/(20*1.6726e-24) + 2e5**2)          # as sobolev/sobolev.py
        tS = Kg*(rho*vD/vr)[:, None, None]
        betg = np.where(tS > 1e-6, -np.expm1(-tS)/np.maximum(tS, 1e-300), 1.)
        Hat = np.array([sum(p['hv'].values()) for p in r['pts']])
        Hsp = {k: np.array([p['hv'][k] for p in r['pts']]) for k in r['pts'][0]['hv']}
        for kfac in (1., 10., 100.):
            for tio in ('lo', 'mid', 'hi'):
                for defac in ((0.5, 1., 2.) if (kfac == 10. and tio == 'mid') else (1.,)):
                    comps = {k: np.zeros(len(TG)) for k in ('atoms', 'cont', 'vibband', 'elband')}
                    molc = {}
                    epsv = np.zeros(len(TG))
                    for it, T in enumerate(TG):
                        S = W*Bb(Teff) - Bb(T)
                        mo, (Cq, Aeff), c = molecules(T, r['P'], W, Teff, vr, kfac, defac)
                        bg = betg[it]
                        ev = (L.GW*Kg[it]*bg*Cq/(Cq + bg*Aeff)).sum(-1)/kl[it]   # (nb,)
                        kel = (L.GW*Kg[it]*bg).sum(-1)                            # (nb,)
                        for k, v in mo.items():
                            molc.setdefault(k, np.zeros(len(TG)))[it] = v
                        epsv[it] = ev[4]
                        fv, Av, qv = TIO[tio]
                        Ct = (c['H'] + c['H2'])*qv
                        fh = fv*Ct/(Ct + Av)
                        comps['cont'][it] = 4*np.pi*rho[it]*(kc[it]*S).sum()
                        comps['vibband'][it] = 4*np.pi*rho[it]*(ev[BVIB]*kl[it, BVIB]*S[BVIB]).sum()
                        comps['elband'][it] = 4*np.pi*rho[it]*fh*(kel[BEL]*S[BEL]).sum()
                    comps['atoms'] = Hat
                    Hm = sum(molc.values())
                    Htot = Hat + Hm + comps['cont'] + comps['vibband'] + comps['elband']
                    rl = []
                    dH = np.gradient(Htot, TG)
                    for T, stab, i, a in roots(TG, Htot):
                        dd = (1-a)*dH[i] + a*dH[i+1]
                        cpv = ((1-a)*cp[i]*rho[i] + a*cp[i+1]*rho[i+1])
                        allc = dict(comps, **molc, **{'at:' + k: v for k, v in Hsp.items()})
                        val = {k: float((1-a)*v[i] + a*v[i+1]) for k, v in allc.items()}
                        heat = {k: v for k, v in val.items() if v > 0 and k != 'atoms'}
                        cool = {k: v for k, v in val.items() if v < 0 and k != 'atoms'}
                        # effective eps: atomic+electronic vs IR (old-model definition)
                        Sx = [W*Bb(Teff) - Bb(T)]
                        kli = (1-a)*kl[i] + a*kl[i+1]
                        rr = (1-a)*rho[i] + a*rho[i+1]
                        den_el = 4*np.pi*rr*(kli[BEL]*Sx[0][BEL]).sum()
                        den_ir = 4*np.pi*rr*(kli[:8]*Sx[0][:8]).sum()
                        num_el = val['atoms'] + val['elband']
                        num_ir = sum(val[k] for k in molc) + val['vibband']
                        rl.append(dict(T=T, stable=stab, tth=cpv/abs(dd) if dd else np.inf,
                                       heat=heat, cool=cool,
                                       eps_el=num_el/den_el, eps_ir=num_ir/den_ir,
                                       epsv=float((1-a)*epsv[i] + a*epsv[i+1])))
                    res.append(dict(star=r['star'], x=r['x'], P=r['P'], ic=r['ic'],
                                    cfg=r['cfg'], kfac=kfac, tio=tio, defac=defac,
                                    roots=rl, H=Htot, Hat=Hat, comps=comps, molc=molc,
                                    Hsp=Hsp, rho=rho, W=W, R=r['R'], dvdr=vr))
    with open(HERE + 'thermal_res.pkl', 'wb') as fh:
        pickle.dump(dict(TG=TG, res=res), fh)


if __name__ == '__main__':
    run()
