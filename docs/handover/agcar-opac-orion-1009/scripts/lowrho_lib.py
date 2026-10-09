"""LTE Planck / Rosseland means at very low density (AG Car thin atmosphere).
Composition, Saha-Boltzmann ionisation, partition functions, continuum and line opacity.
All cgs.  See NOTE-2026-10-09-orion-lowrho-planck.md for the physics choices."""
import os
import numpy as np
import numba as nb

B = '/orion/ptmp/jinma/agcar_opac_1009/'
h, kB, c, me, qe = 6.62607015e-27, 1.380649e-16, 2.99792458e10, 9.1093837e-28, 4.80320471e-10
amu, eV, sigT = 1.66053907e-24, 1.602176634e-12, 6.6524587e-25
sigSB = 5.670374419e-5
PIE2MC = np.pi*qe**2/(me*c)               # 0.02654 cm^2 Hz
RYD = 13.605693                           # eV
SAHA = (2*np.pi*me*kB/h**2)**1.5          # x T^1.5 cm^-3
X_H, Z_MET = 0.36, 0.02
XI_TURB = float(os.environ.get('XI_TURB', '2.0e5'))   # microturbulence, cm/s
NMAX_H = int(os.environ.get('NMAX_H', '30'))   # hydrogenic level cutoff (H I, He II)

# GS98 number fractions within the metals (LANL grevsau1, = data/stellar_opac/planck_tools)
GS98 = {6: 2.45187e-01, 7: 6.15882e-02, 8: 5.00608e-01, 10: 8.90220e-02, 11: 1.58306e-03,
        12: 2.81512e-02, 13: 2.18523e-03, 14: 2.62723e-02, 15: 2.08688e-04, 16: 1.58306e-02,
        17: 2.34152e-04, 18: 1.85993e-03, 19: 9.76107e-05, 20: 1.69628e-03, 21: 1.09521e-06,
        22: 7.75349e-05, 23: 7.40453e-06, 24: 3.46336e-04, 25: 1.81760e-04, 26: 2.34152e-02,
        27: 6.15882e-05, 28: 1.31673e-03, 29: 1.20087e-05, 30: 2.94780e-05}
AW = {1: 1.008, 2: 4.0026, 6: 12.011, 7: 14.007, 8: 15.999, 10: 20.180, 11: 22.990,
      12: 24.305, 13: 26.982, 14: 28.085, 15: 30.974, 16: 32.06, 17: 35.45, 18: 39.948,
      19: 39.098, 20: 40.078, 21: 44.956, 22: 47.867, 23: 50.942, 24: 51.996, 25: 54.938,
      26: 55.845, 27: 58.933, 28: 58.693, 29: 63.546, 30: 65.38}
ZL = [1, 2] + sorted(GS98)
NST = 10


def composition():
    """number of nuclei of each element per gram of gas"""
    mbar = sum(GS98[z]*AW[z] for z in GS98)
    N = {1: X_H/(AW[1]*amu), 2: (1 - X_H - Z_MET)/(AW[2]*amu)}
    for z in GS98:
        N[z] = Z_MET*GS98[z]/mbar/amu
    return N


NPG = composition()


class Atoms:
    def __init__(self, f=B + 'work/atoms_%s.npz' % os.environ.get('LINES', 'cd1')):
        d = np.load(f)
        self.d = {k: d[k] for k in d.files}
        self.ions = [(z, s) for z in ZL for s in range(0, min(z, NST) + 1)]

    def levels(self, z, s):
        """(g, E [eV], ip [eV], ground-config flag); hydrogenic n-levels for H I, He II"""
        k = f'{z}_{s}'
        ip = self.d[k + '_ip'][0]
        if (z == 1 and s == 0) or (z == 2 and s == 1):
            n = np.arange(1, NMAX_H + 1)
            Zc = s + 1
            return 2.0*n**2, ip*(1 - 1/n**2), ip, n == 1, n
        return self.d[k + '_g'], self.d[k + '_E'], ip, self.d[k + '_gc'], None

    def U(self, z, s, T):
        g, E, ip, _, _ = self.levels(z, s)
        return np.sum(g*np.exp(-E*eV/(kB*T)))


def saha(at, T, rho, UT=None, npg=None):
    """LTE ionisation.  npg: nuclei per gram per element (default NPG; pass a depleted
    copy for atoms locked in molecules).  returns n_e, dict (z,s)->n [cm^-3], n_Hminus"""
    if npg is None:
        npg = NPG
    if UT is None:
        UT = {k: at.U(*k, T) for k in at.ions}
    kT = kB*T/eV
    lnS = np.log(SAHA*T**1.5)
    # log of the stage ratios times n_e : ln(n_{s+1} n_e / n_s)
    lr = {}
    for z in ZL:
        smax = min(z, NST)
        lr[z] = np.array([np.log(2*UT[(z, s+1)]/UT[(z, s)]) + lnS
                          - at.d[f'{z}_{s}_ip'][0]/kT for s in range(smax)])
    ntot = {z: npg[z]*rho for z in ZL}

    def fracs(lne):
        out = {}
        q = 0.0
        for z in ZL:
            c = np.concatenate([[0.0], np.cumsum(lr[z] - lne)])
            c -= c.max()
            f = np.exp(c)
            f /= f.sum()
            out[z] = f
            q += ntot[z]*np.sum(np.arange(len(f))*f)
        return out, q

    def hminus(nH, ne):
        return nH*ne*np.exp(0.754/kT - lnS)/(2*UT[(1, 0)])   # Saha, g(H-) = 1

    lo, hi = np.log(1e-30), np.log(2*sum(ntot[z]*z for z in ZL))
    for _ in range(200):
        mid = 0.5*(lo + hi)
        f, q = fracs(mid)
        q -= hminus(ntot[1]*f[1][0], np.exp(mid))
        if q > np.exp(mid):
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-10:
            break
    ne = np.exp(0.5*(lo + hi))
    f, _ = fracs(np.log(ne))
    n = {(z, s): ntot[z]*f[z][s] for z in ZL for s in range(len(f[z]))}
    return ne, n, hminus(n[(1, 0)], ne)


# ------------------------------------------------------------------ continuum
def planck_nu(nu, T):
    x = h*nu/(kB*T)
    return 2*h*nu**3/c**2/np.expm1(x)


def dBdT_nu(nu, T):
    x = h*nu/(kB*T)
    ex = np.exp(-x)
    return 2*h*nu**3/c**2*x/T*ex/(1 - ex)**2


def bf_ion(at, z, s, nu, T, U):
    """bound-free cross section per particle of ion (z,s) [cm^2], Boltzmann-averaged over
    its levels, NOT yet corrected for stimulated emission.  Ground configuration: Verner
    et al. (1996) fit shifted by the level energy; all other levels: hydrogenic Kramers
    with Z_eff = s+1, n_eff from the binding energy, g_bf = 1."""
    g, E, ip, gc, nl = at.levels(z, s)
    hnu = h*nu/eV
    w = g*np.exp(-E*eV/(kB*T))/U
    out = np.zeros_like(nu)
    k = f'{z}_{s}'
    Zc = s + 1
    kr = []
    for j in range(len(g)):
        if w[j] < 1e-30:
            continue
        eb = ip - E[j]
        if eb <= 0:
            continue
        if gc[j] and (k + '_vE') in at.d and nl is None:
            vE, vS = at.d[k + '_vE'], at.d[k + '_vS']
            out += w[j]*np.interp(hnu + E[j], vE, vS, left=0.0, right=0.0)
        elif nl is not None and nl[j] == 1:
            vE, vS = at.d[k + '_vE'], at.d[k + '_vS']
            out += w[j]*np.interp(hnu, vE, vS, left=0.0, right=0.0)
        else:
            kr.append((eb, w[j]*Zc*np.sqrt(RYD/eb)*eb**3))
    if kr:      # Kramers levels: sigma = 7.907e-18/Zc^2/hnu^3 * sum_{eb<=hnu} w neff eb^3
        kr = np.array(sorted(kr))
        cs = np.concatenate([[0.0], np.cumsum(kr[:, 1])])
        out += 7.907e-18/Zc**2/hnu**3*cs[np.searchsorted(kr[:, 0], hnu, side='right')]
    return out


def hminus_bf(nu):
    """Gray (2005) H- bound-free per H- ion [cm^2] (2250-16419 A; flat below 2250 A)"""
    lam = np.clip(c/nu*1e8, 2250.0, None)
    a = [1.99654, -1.18267e-5, 2.64243e-6, -4.40524e-10, 3.23992e-14, -1.39568e-18,
         2.78701e-23]
    s = sum(a[i]*lam**i for i in range(7))*1e-18
    return np.where(c/nu*1e8 < 16419.0, np.clip(s, 0, None), 0.0)


def hminus_ff(nu, T):
    """Gray (2005) H- free-free per H atom per unit electron pressure [cm^4/dyn],
    lambda clipped to 1823 A .. 15 um (lambda^2 scaling beyond)"""
    lam0 = c/nu*1e8
    lam = np.clip(lam0, 1823.0, 1.5e5)
    L = np.log10(lam)
    th = np.log10(5040.0/T)
    f0 = -2.2763 - 1.6850*L + 0.76661*L**2 - 0.053346*L**3
    f1 = 15.2827 - 9.2846*L + 1.99381*L**2 - 0.142631*L**3
    f2 = -197.789 + 190.266*L - 67.9775*L**2 + 10.6913*L**3 - 0.625151*L**4
    k = 10**(-26 + f0 + f1*th + f2*th**2)
    return k*np.where(lam0 > 1.5e5, (lam0/1.5e5)**2, 1.0)


def ff_hydrogenic(nu, T, ne, n_ion_by_charge):
    """free-free (Kramers, g_ff = 1) per cm, incl. stimulated emission.
    n_ion_by_charge: sum_s s^2 n_s"""
    return 3.692e8*n_ion_by_charge*ne/np.sqrt(T)/nu**3*(-np.expm1(-h*nu/(kB*T)))


def rayleigh(nu):
    """(sigma_H, sigma_He, sigma_H2) [cm^2]: Kurucz/ATLAS H (Dalgarno), He, Dalgarno &
    Williams H2; lambda clipped at 1300 A (H, H2) and 600 A (He) on the blue side"""
    lam = c/nu*1e8
    lh = np.clip(lam, 1300.0, None)
    sH = 5.799e-13/lh**4 + 1.422e-6/lh**6 + 2.784/lh**8
    sH2 = 8.14e-13/lh**4 + 1.28e-6/lh**6 + 1.61/lh**8
    le = np.clip(lam, 600.0, None)
    sHe = 5.484e-14/le**4*(1 + 2.44e5/le**2 + 5.94e10/(le**2 - 2.90e5)**2)**2
    return sH, sHe, sH2


# ------------------------------------------------------------------ lines
@nb.njit(cache=True)
def _humlicek(x, y):
    """Humlicek (1982) W4 complex probability function, real part = Voigt H(a=y, v=x)"""
    t = complex(y, -x)
    s = abs(x) + y
    if s >= 15.0:
        w = t*0.5641896/(0.5 + t*t)
    elif s >= 5.5:
        u = t*t
        w = t*(1.410474 + u*0.5641896)/(0.75 + u*(3.0 + u))
    elif y >= 0.195*abs(x) - 0.176:
        w = (16.4955 + t*(20.20933 + t*(11.96482 + t*(3.778987 + t*0.5642236))))/(
            16.4955 + t*(38.82363 + t*(39.27121 + t*(21.69274 + t*(6.699398 + t)))))
    else:
        u = t*t
        w = np.exp(u) - t*(36183.31 - u*(3321.9905 - u*(1540.787 - u*(219.0313 - u*(
            35.76683 - u*(1.320522 - u*0.56419))))))/(32066.6 - u*(24322.84 - u*(
                9022.228 - u*(2186.181 - u*(364.2191 - u*(61.57037 - u*(1.841439 - u)))))))
    return w.real


@nb.njit(cache=True)
def add_lines(sig, lnu0, dln, nu, S, dnuD, a, nwin):
    """sig[i] += S phi(nu_i), Voigt profile (Doppler width dnuD, damping a), sampled at
    the grid points (opacity sampling), window +- nwin Doppler widths."""
    n = sig.shape[0]
    for l in range(nu.shape[0]):
        x0 = (np.log(nu[l]) - lnu0)/dln
        wpts = nwin[l]*dnuD[l]/nu[l]/dln
        i0 = max(0, int(x0 - wpts))
        i1 = min(n - 1, int(x0 + wpts) + 1)
        norm = S[l]/(1.7724538509*dnuD[l])
        for i in range(i0, i1 + 1):
            nui = np.exp(lnu0 + i*dln)
            v = (nui - nu[l])/dnuD[l]
            sig[i] += norm*_humlicek(v, a[l])


def line_strengths(at, z, s, T, U):
    """per particle of ion (z,s): nu, S = (pi e^2/mc) gf exp(-El/kT)(1-exp(-h nu/kT))/U,
    Doppler width, damping parameter a (radiative only; Gamma from gfall, else classical)"""
    k = f'{z}_{s}_lines'
    if k not in at.d:
        return None
    nu, gf, El, gR = at.d[k]
    S = PIE2MC*gf*np.exp(-El*eV/(kB*T))*(-np.expm1(-h*nu/(kB*T)))/U
    m = AW[z]*amu
    b = np.sqrt(2*kB*T/m + XI_TURB**2)
    dnuD = nu*b/c
    gam = np.where(gR > 0, gR, 0.2223e16/(c/nu*1e8)**2)
    a = gam/(4*np.pi*dnuD)
    return nu, S, dnuD, a
