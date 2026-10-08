"""Two-phase (warm + cool clump) gas + dust radiation-force test for RSG winds.

Part 1  optically thin radiative equilibrium (RE) of the gas in the diluted stellar field,
        sum_b sum_g w_g k^abs_bg(T, rho) [W(r) B_b(Teff) - B_b(T)] = 0
        (Rayleigh scattering removed from k^abs; continuum CIA + H- kept), all roots in
        [800, 4000] K at fixed rho and at fixed P (isobaric = Field criterion), stability
        from the sign of dH/dT along the constraint.
Part 2  two-phase layer on the MESA columns: pressure balance + mass conservation, gas
        Sobolev force per phase (sobolev.gamma_sob, flux A) and the thin value.
Part 3  dust in the cool phase: two condensation criteria, delta f_cond kappa_d scan,
        clump porosity (1 - e^-tau)/tau.
Run:  python twophase.py   (writes part1_*.npz/png, part2_*.npz, tables_*.md, plots)
"""
import itertools
import os
import sys
import numpy as np
sys.dont_write_bytecode = True
os.environ.setdefault('CK_DATA',
                      '/orion/u/jinma/ATHENAK/athenak/docs/handover/rsg-ck-1008/data/')
import rsglib as L      # noqa: E402
import sobolev as S     # noqa: E402
ck = L.ck

OUT = os.path.dirname(os.path.abspath(__file__)) + '/'
MESA_D = '/orion/ptmp/jinma/rsg_wind_1008/mesa/'
TG = np.geomspace(300., 4000., 321)   # extended below 800 K: cold branch
RHOG = np.logspace(-17, -12, 26)
PG = np.logspace(-7, -1, 31)              # dyn/cm^2
XR = np.array([1.5, 2., 3., 4., 5.])


def dilution(x):
    return 0.5*(1 - np.sqrt(np.maximum(1 - 1/np.asarray(x, float)**2, 0)))


def kray(T, P, mu):
    """Rayleigh part of the per-gram continuum (nb), as continuum() builds it."""
    v = ck.comp(L.CE, T, P*1e-6)
    ntot = P/(L.KB*T)
    rho = P*mu*L.MH/(L.KB*T)
    return sum(L.RAY[s]*v[s+1] for s in range(4))*ntot/rho


def state_abs(T, rho=None, P=None, mode='clamp'):
    """L.state at (T, rho) or (T, P); adds Kabs = K - Rayleigh."""
    if rho is None:
        mu = 2.3
        for _ in range(40):
            rho = P*mu*L.MH/(L.KB*T)
            mu1 = ck.comp(L.CE, T, P*1e-6)[0]
            if abs(mu1 - mu) < 1e-7*mu:
                break
            mu = mu1
        rho = P*mu1*L.MH/(L.KB*T)
    s = L.state(T, rho, mode)
    s['rho'] = rho
    s['Kabs'] = s['K'] - kray(T, s['P'], s['mu'])[:, None]
    return s


def Bb(T):
    """band-integrated Planck intensity B_b(T), (n, nb)."""
    T = np.atleast_1d(T)
    return L.SIG*T[:, None]**4*L.fband(T)/np.pi


def heating(Kabs, T, W, Teff):
    """net heating per gram / (4 pi): (n,) for Kabs (n, nb, ng) at gas T (n,)."""
    kb = (L.GW*Kabs).sum(-1)
    return (kb*(W*Bb(Teff) - Bb(T))).sum(-1)


def roots(T, H):
    """all sign changes of H(T) (log-T linear interpolation); stable if dH/dT < 0."""
    out = []
    lT = np.log(T)
    for i in np.where(np.sign(H[:-1]) != np.sign(H[1:]))[0]:
        a = H[i]/(H[i] - H[i+1])
        out.append((float(np.exp(lT[i] + a*(lT[i+1] - lT[i]))), bool(H[i+1] < H[i])))
    return out


def part1_tables(mode):
    fn = OUT + f'part1_K_{mode}.npz'
    if os.path.exists(fn):
        return dict(np.load(fn))
    nb, ng = L.NB, L.NG
    Kr = np.zeros((len(RHOG), len(TG), nb, ng))
    Kp = np.zeros((len(PG), len(TG), nb, ng))
    mur = np.zeros((len(RHOG), len(TG)))
    mup = np.zeros((len(PG), len(TG)))
    rhop = np.zeros((len(PG), len(TG)))
    for (i, rho), (j, T) in itertools.product(enumerate(RHOG), enumerate(TG)):
        s = state_abs(T, rho=rho, mode=mode)
        Kr[i, j], mur[i, j] = s['Kabs'], s['mu']
    for (i, P), (j, T) in itertools.product(enumerate(PG), enumerate(TG)):
        s = state_abs(T, P=P, mode=mode)
        Kp[i, j], mup[i, j], rhop[i, j] = s['Kabs'], s['mu'], s['rho']
    d = dict(Kr=Kr, Kp=Kp, mur=mur, mup=mup, rhop=rhop)
    np.savez(fn, **d)
    return d


def part1(star, mode):
    st = L.star(star)
    tab = part1_tables(mode)
    res = {'rho': {}, 'P': {}}
    for x in XR:
        W = dilution(x)
        for key, K in (('rho', tab['Kr']), ('P', tab['Kp'])):
            res[key][x] = []
            for i in range(K.shape[0]):
                H = heating(K[i], TG, W, st['Teff'])
                res[key][x].append(roots(TG, H))
    return st, tab, res


def two_stable(rl):
    st = [t for t, s in rl if s]
    return (max(st), min(st)) if len(st) >= 2 else None


# ---------------- Part 2 / 3 ----------------
SPECIES = {  # T_cond range [K], delta_max (dust/gas mass), p for the grain T
    'Al2O3': dict(Tc=(1400., 1600.), delta=1.0e-4),
    'silicate': dict(Tc=(1000., 1200.), delta=4.0e-3),
}
FF = (0.01, 0.1, 0.3)
FCOND = (0.1, 0.3, 1.0)
KD = (300., 1000., 3000.)
LCL = (0., 0.01, 0.1)       # clump size / r (0 = no porosity)
TCSCAN = (1000., 1200., 1400.)
MODELS = ('m15lgl5.1', 'm20lgl5.3', 'm20lgl5.5')


def tgrain(Teff, W, p):
    return Teff*W**(1./(4 + p))


def porosity(tau):
    tau = np.asarray(tau, float)
    return np.where(tau > 1e-8, -np.expm1(-tau)/np.maximum(tau, 1e-300), 1 - 0.5*tau)


def phase_split(rho_bar, Tw, Tc, f, mode):
    """pressure balance + mass conservation, FastChem mu iterated."""
    chi = Tw/Tc
    for _ in range(30):
        rw = rho_bar/(f*chi + 1 - f)
        rc = chi*rw
        sw, sc = L.state(Tw, rw, mode), L.state(Tc, rc, mode)
        chi1 = (Tw*sc['mu'])/(Tc*sw['mu'])
        if abs(chi1 - chi) < 1e-6*chi:
            break
        chi = chi1
    return rw, rc, sw, sc


def isobaric_pair(st, tab, x, rho_bar, f, mode):
    """(Tw, Tc) the two stable isobaric roots at r, with P chosen so that
    f rho_c + (1-f) rho_w = rho_bar.  None if no P on the grid has two stable roots."""
    W = dilution(x)
    best = []
    for i in range(len(PG)):
        H = heating(tab['Kp'][i], TG, W, st['Teff'])
        pr = two_stable(roots(TG, H))
        if pr is None:
            continue
        Tw, Tc = pr
        rw = PG[i]*np.interp(Tw, TG, tab['mup'][i])*L.MH/(L.KB*Tw)
        rc = PG[i]*np.interp(Tc, TG, tab['mup'][i])*L.MH/(L.KB*Tc)
        best.append((np.log(f*rc + (1-f)*rw), Tw, Tc, PG[i]))
    if not best:
        return None
    best = np.array(best)
    lr = np.log(rho_bar)
    if lr < best[:, 0].min() or lr > best[:, 0].max():
        return None   # rho_bar outside the two-phase P band (no extrapolation)
    o = np.argsort(best[:, 0])
    return (float(np.interp(lr, best[o, 0], best[o, 1])),
            float(np.interp(lr, best[o, 0], best[o, 2])),
            float(np.exp(np.interp(lr, best[o, 0], np.log(best[o, 3])))))


def load_col(m, case, mode):
    d = np.load(f'{MESA_D}col_{m}_{case}_{mode}.npz')
    st = L.star(m)
    st = dict(st, Teff=float(d['Teff']), kE=float(d['kE']))
    r = d['r']*st['R']
    rho, T = d['rho'], d['T']
    H = np.abs(1./np.gradient(np.log(rho), r))
    x = d['r']
    out = dict(x=XR, rho=np.exp(np.interp(XR, x, np.log(rho))),
               T=np.interp(XR, x, T), H=np.exp(np.interp(XR, x, np.log(H))),
               vc=float(d['v']), xmax=float(x.max()))
    return st, out


def gas_force(st, x, rho, T, dvdr, mode):
    r = np.array([x*st['R']])
    K = S.opacities(np.array([T]), np.array([rho]), mode)[0]
    G, _ = S.gamma_sob(st, r, np.array([rho]), np.array([T]), dvdr, mode, 2e5, 20., 'A',
                       K=K)
    return float(G[0]), float(S.gamma_thin(st, K)[0])


def part2(tabs):
    rows = []
    for m, case, mode in itertools.product(MODELS, ('mft', 'mft12'), ('clamp', 'extrap')):
        st, c = load_col(m, case, mode)
        tab = tabs[mode]
        for ix, x in enumerate(XR):
            if x > c['xmax']:
                continue
            rb, Tb, H, vc = c['rho'][ix], c['T'][ix], c['H'][ix], c['vc']
            W = dilution(x)
            for f in FF:
                pairs = []
                ip = isobaric_pair(st, tab, x, rb, f, mode)
                if ip is not None:
                    pairs.append(('RE', ip[0], ip[1]))
                pairs += [('scan', Tb, tc) for tc in TCSCAN if tc < Tb]
                for src, Tw, Tc in pairs:
                    rw, rc, sw, sc = phase_split(rb, Tw, Tc, f, mode)
                    fm = f*rc/rb
                    for gname, dv in (('v/H', vc/H), ('v/r', vc/(x*st['R']))):
                        Gw, Gwt = gas_force(st, x, rw, Tw, dv, mode)
                        Gc, Gct = gas_force(st, x, rc, Tc, dv, mode)
                        rows.append(dict(model=m, case=case, mode=mode, x=x, f=f, src=src,
                                         Tw=Tw, Tc=Tc, Tbar=Tb, rhobar=rb, rw=rw, rc=rc,
                                         fm=fm, grad=gname, Gw=Gw, Gwthin=Gwt, Gcgas=Gc,
                                         Gcthin=Gct, kE=st['kE'], Teff=st['Teff'],
                                         R=st['R'], W=W))
    return rows


def dust_rows(rows):
    """expand each Part-2 row over species x criterion x f_cond x kappa_d x l."""
    out = []
    for rw in rows:
        r = rw['x']*rw['R']
        for sp, spd in SPECIES.items():
            for crit in ('gas', 'grain_p1', 'grain_p0', 'grain_pm1'):
                if crit == 'gas':
                    Tcmp = rw['Tc']
                else:
                    p = {'grain_p1': 1., 'grain_p0': 0., 'grain_pm1': -1.}[crit]
                    Tcmp = tgrain(rw['Teff'], rw['W'], p)
                for Tcond, tag in ((spd['Tc'][1], 'opt'), (spd['Tc'][0], 'nom')):
                    ok = Tcmp < Tcond
                    for fc, kd, lc in itertools.product(FCOND, KD, LCL):
                        kdust = spd['delta']*fc*kd if ok else 0.
                        kgas = rw['Gcgas']*rw['kE']
                        tau = (kgas + kdust)*rw['rc']*lc*r
                        por = porosity(tau)
                        Gc = (kgas + kdust)*por/rw['kE']
                        Gbar = rw['fm']*Gc + (1 - rw['fm'])*rw['Gw']
                        out.append(dict(rw, sp=sp, crit=crit, Tcond=tag, cond=ok, fcond=fc,
                                        kd=kd, l=lc, tau=float(tau), por=float(por),
                                        Gdust=kdust*por/rw['kE'], Gc=float(Gc),
                                        Gbar=float(Gbar)))
    return out


def x_required(rw, lc):
    """minimum delta f_cond kappa_d [cm^2/g gas] for Gamma_c > 1 (inf if impossible)."""
    kE, kg = rw['kE'], rw['Gcgas']*rw['kE']
    s = rw['rc']*lc*rw['x']*rw['R']
    if lc == 0:
        return max(kE - kg, 0.)
    if 1./s <= kE:                      # max effective opacity 1/(rho_c l)
        return np.inf
    from scipy.optimize import brentq
    g = (lambda X: (kg + X)*porosity((kg + X)*s) - kE)
    if g(0.) >= 0:
        return 0.
    hi = kE
    while g(hi) < 0:
        hi *= 2
    return brentq(lambda X: g(X), 0., hi, rtol=1e-6) if hi < 1e30 else np.inf
