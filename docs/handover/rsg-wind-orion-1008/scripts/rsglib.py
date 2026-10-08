"""RSG molecular radiation-force feasibility: opacity lookups on the Exo-FMS ck tables.

All cgs.  Uses ck_lib.py (copy of data/exo_fms_ck/tools/ck_lib.py) for the table readers,
continuum (CIA + Rayleigh + H- John 1988) and band Planck fractions.
P below the table edge (1e-8 bar) has two modes:
  'clamp'  : line k and composition at 1e-8 bar (upper bound on molecular abundance)
  'extrap' : band-mean line k extrapolated in log P from the two lowest P nodes, slope
             clipped >= 0 (opacity may only FALL toward lower P: dissociation); the
             g-shape is kept from 1e-8 bar.  Lower bound.
"""
import os
import sys
import numpy as np
sys.dont_write_bytecode = True
os.environ.setdefault('CK_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data') + '/')
import ck_lib as ck  # noqa: E402

G = 6.674e-8
C = 2.99792458e10
MSUN = 1.989e33
LSUN = 3.828e33
RSUN = 6.957e10
SIG = 5.670374419e-5
KB = 1.380649e-16
MH = 1.6726e-24
YR = 3.156e7

KT = ck.read_ktable(ck.DATA + 'ck/Premixed_1x_g8_11_hiT2.txt')
CE = ck.read_ce(ck.DATA + 'CE_tables/FastChem_ck_1x_int_hiT2.txt')
CIA = ck.read_cia()
RAY = ck.read_ray()
WL = KT['wl']
GW = KT['gw']
NB, NG = len(WL) - 1, len(GW)
LK = np.log10(np.maximum(KT['K'], 1e-99))
LTg = np.log10(KT['T'])
LPg = np.log10(KT['P'])          # bar
# band-mean (over g) line k, for the extrap mode
LKBAR = np.log10(np.maximum((KT['K'] * GW).sum(-1), 1e-99))   # (nT, nP, nb)

# band Planck fractions on a fine T grid (fast lookup)
_TF = np.linspace(300., 12000., 2341)
_FB = ck.planck_frac(_TF, WL)    # (nT, nb)


def fband(T):
    T = np.atleast_1d(T)
    return np.array([np.interp(T, _TF, _FB[:, b]) for b in range(NB)]).T  # (n, nb)


def band_centers():
    return 2. / (1. / WL[:-1] + 1. / WL[1:])


def linek(T, Pbar, mode='clamp'):
    """(nb, ng) line k [cm^2/g] at scalar T, P[bar]."""
    lT = np.log10(T)
    lP = np.log10(Pbar)
    i, a = ck.interp_lin(lT, LTg, None)
    j, b = ck.interp_lin(lP, LPg, None)
    lk = (1-a)*((1-b)*LK[i, j]+b*LK[i, j+1]) + a*((1-b)*LK[i+1, j]+b*LK[i+1, j+1])
    if mode == 'extrap' and lP < LPg[0]:
        kb0 = (1-a)*LKBAR[i, 0] + a*LKBAR[i+1, 0]
        kb1 = (1-a)*LKBAR[i, 1] + a*LKBAR[i+1, 1]
        slope = np.maximum((kb1 - kb0) / (LPg[1] - LPg[0]), 0.)
        lk = lk + (slope * (lP - LPg[0]))[:, None]
    return 10**lk


def state(T, rho, mode='clamp'):
    """Return dict: P [dyn], mu, K (nb,ng) total per-gram, kc (nb), kline (nb,ng)."""
    mu = 2.3
    for _ in range(30):
        P = rho * KB * T / (mu * MH)
        v = ck.comp(CE, T, P * 1e-6)
        if abs(v[0] - mu) < 1e-6 * mu:
            break
        mu = v[0]
    P = rho * KB * T / (mu * MH)
    kc, rho_c, v = ck.continuum(CE, CIA, RAY, WL, T, P * 1e-6)
    # continuum() derives rho from P with its own mu; rescale per-gram to our rho
    kc = kc * rho_c / rho
    kl = linek(T, P * 1e-6, mode)
    return dict(P=P, mu=mu, K=kl + kc[:, None], kc=kc, kline=kl)


def means(T, K):
    """Rosseland (gas T) and Planck (gas T) means of K (nb, ng)."""
    fp = fband(np.array([1.01*T, 0.99*T]))
    wb = SIG*(1.01*T)**4*fp[0] - SIG*(0.99*T)**4*fp[1]
    inv = (GW / K).sum(1)
    ok = wb > 0
    kR = wb[ok].sum() / (wb[ok]*inv[ok]).sum()
    kP = (fband(T)[0] * (GW*K).sum(1)).sum()
    return kR, kP


def kappa_thin_flux(K, Trad):
    """Flux mean in the optically thin limit for a B_nu(Trad) stellar spectrum."""
    return (fband(Trad)[0] * (GW*K).sum(1)).sum()


# ---------------- stellar Rosseland table used by the red_giant runs ----------------
_OP = os.path.join(ck.DATA, 'rosseland_gs98_x0.7_z0.014.txt')  # PATCH orion: honour CK_DATA


def _read_op():
    hdr = [ln for ln in open(_OP) if ln.startswith('#')]
    g = None
    for i, ln in enumerate(hdr):
        if 'grid:' in ln:
            g = hdr[i+1].strip('# \n').split()
    nT, nD = int(g[0]), int(g[1])
    lT0, dlT, lD0, dlD = map(float, g[2:])
    dat = np.loadtxt(_OP, comments='#').ravel()[:nT*nD].reshape(nT, nD)
    return lT0 + dlT*np.arange(nT), lD0 + dlD*np.arange(nD), dat


OPT, OPD, OPK = _read_op()


def kappa_aesopus(T, rho):
    from scipy.interpolate import RegularGridInterpolator
    f = RegularGridInterpolator((OPT, OPD), OPK, bounds_error=False, fill_value=None)
    lR = np.log10(rho) - 3*np.log10(T) + 18
    return 10**f([[np.log10(T), np.log10(rho)]])[0], lR


# ---------------- MESA profiles (read-only) ----------------
MESA_DIR = '/orion/ptmp/jinma/Arepo/ICs_golden/ic/mesa-profile/'
MESA = {k: MESA_DIR + 'RSGg/' + k + '.data' for k in (
    'm10lgl4.4a2', 'm10lgl4.6', 'm10lgl4.8', 'm10lgl5.0', 'm15lgl5.1', 'm20lgl5.3',
    'm20lgl5.4', 'm20lgl5.5')}
MESA['m26_p455'] = MESA_DIR + 'golden/m26_p455.data'


def read_mesa(name):
    """Header dict + columns, arrays SURFACE FIRST (MESA order), cgs."""
    ln = open(MESA[name]).read().splitlines()
    h = dict(zip(ln[1].split(), ln[2].split()))
    cols = ln[5].split()
    d = np.loadtxt(MESA[name], skiprows=6)
    c = dict(zip(cols, d.T))
    out = dict(M=float(h['star_mass'])*MSUN, L=float(h['photosphere_L'])*LSUN,
               R=float(h['photosphere_r'])*RSUN, Teff=float(h['Teff']),
               r=10**c['logR']*RSUN, rho=10**c['logRho'], T=10**c['logT'], P=10**c['logP'],
               prad=c['prad'], kap=c['opacity'], cv=c['conv_vel'])
    # tau from the surface (MESA surface = photosphere, tau = 2/3)
    dr = -np.diff(out['r'])
    kr = out['kap']*out['rho']
    out['tau'] = 2/3. + np.concatenate([[0.], np.cumsum(0.5*(kr[1:] + kr[:-1])*dr)])
    return out


def ft_vcon(m):
    """Fuller & Tsuna (2024): conv_vel where tau = P_rad c / (P v_con), outermost
    crossing scanning inward from the surface."""
    cv = m['cv']
    crit = np.where(cv > 0, m['prad']*C/(m['P']*np.maximum(cv, 1e-99)), np.inf)
    f = np.log(m['tau']) - np.log(crit)          # < 0 near the surface
    k = np.where(f >= 0)[0]
    info = {}
    if len(k) == 0:
        j = np.where(cv > 0)[0][0]
        info.update(ok=False, note='no crossing; nearest convective point')
        return cv[j], m['r'][j], m['tau'][j], info
    k = k[0]
    if k == 0 or not np.isfinite(f[k-1]):
        info.update(ok=(k > 0), note='crossing adjacent to cv=0 zone; zone value')
        return cv[k], m['r'][k], m['tau'][k], info
    a = f[k-1]/(f[k-1] - f[k])
    info.update(ok=True, note='')
    return ((1-a)*cv[k-1] + a*cv[k], (1-a)*m['r'][k-1] + a*m['r'][k],
            np.exp((1-a)*np.log(m['tau'][k-1]) + a*np.log(m['tau'][k])), info)


# ---------------- stars ----------------
def star(name):
    if name in MESA:
        m = read_mesa(name)
        M, L_, R = m['M'], m['L'], m['R']
        return dict(name=name, M=M, L=L_, R=R, Teff=(L_/(4*np.pi*SIG*R**2))**0.25,
                    Teff_mesa=m['Teff'], kE=4*np.pi*C*G*M/L_, g=G*M/R**2)
    if name == 'golden16':
        M = 15.2596 * MSUN
        kE = 1.73439
        L = 4*np.pi*C*G*M/kE
        R = 669.875 * RSUN
    elif name == 'betelgeuse':
        M = 18. * MSUN
        L = 10**5.1 * LSUN
        Teff = 3600.
        R = np.sqrt(L/(4*np.pi*SIG*Teff**4))
        kE = 4*np.pi*C*G*M/L
    Teff = (L/(4*np.pi*SIG*R**2))**0.25
    return dict(name=name, M=M, L=L, R=R, Teff=Teff, kE=kE, g=G*M/R**2)


def read_marcs(fn):
    lines = open(fn).read().splitlines()
    i1 = [i for i, ln in enumerate(lines) if ln.startswith(' k lgTauR  lgTau5')][0]
    i2 = [i for i, ln in enumerate(lines) if ln.startswith(' k lgTauR    KappaRoss')][0]
    n = i2 - i1 - 1
    a = np.array([ln.split() for ln in lines[i1+1:i1+1+n]], float)
    b = np.array([ln.split() for ln in lines[i2+1:i2+1+n]], float)
    hdr = {}
    for ln in lines[:10]:
        if 'Radius' in ln:
            hdr['R'] = float(ln.split()[0])
        if 'Luminosity' in ln:
            hdr['L'] = float(ln.split()[0]) * LSUN
        if 'Teff' in ln:
            hdr['Teff'] = float(ln.split()[0])
        if 'Surface gravity' in ln:
            hdr['g'] = float(ln.split()[0])
        if 'Mass' in ln:
            hdr['M'] = float(ln.split()[0]) * MSUN
    return dict(hdr=hdr, lgtauR=a[:, 1], depth=a[:, 3], T=a[:, 4], Pg=a[:, 6],
                kR=b[:, 2], rho=b[:, 3], mu=b[:, 4])
