"""Sobolev (velocity-gradient) desaturation of the RSG correlated-k radiation force.

Each (band b, g-point g) is a line ensemble covering the fraction w_g of band b with total
opacity k_bg [cm^2/g] (lines + continuum at that g).  Sobolev optical depth over the
velocity-coherence length L_S = v_D/|dv/dr|:
    tau_S,bg = k_bg rho v_D / |dv/dr|,  v_D = sqrt(2 k T/(mu_abs m_H) + xi^2)
Radiative acceleration in Eddington units (CAK escape factor, radial streaming):
    Gamma_Sob(r) = (1/kappa_Edd) sum_b F_b/F sum_g w_g k_bg (1 - exp(-tau_S))/tau_S
Band flux fractions:
    'A' : unattenuated stellar f_b(Teff)  (upper bound; Gamma_Sob <= Gamma_thin exactly)
    'B' : f_b exp(-tau_win,b(r)), tau_win,b = int_{r_ph}^{r} k_b,g0 rho dr' (g0 = the band's
          smallest-k g-point; r_ph = tau_R=2/3 radius), renormalised to sum 1.
Usage (library):  G, contrib = gamma_sob(st, r, rho, T, dvdr, mode, xi, mu_abs, flux)
  st from rsglib.star(); r [cm]; rho [g/cm3]; T [K]; dvdr [1/s]; xi [cm/s].
  contrib (n, nb) sums over b to G.  Pass K=(n, nb, ng) to skip the opacity lookups and
  r_ph [cm] for flux 'B' (default r[0]).
"""
import numpy as np
import rsglib as L

SIGREF = 0.2   # cm^2/g, CAK reference opacity for the force multiplier M = Gamma_Sob/Gamma_e


def opacities(T, rho, mode):
    K = np.zeros((len(T), L.NB, L.NG))
    kR = np.zeros(len(T))
    for i in range(len(T)):
        K[i] = L.state(T[i], rho[i], mode)['K']
        kR[i] = L.means(T[i], K[i])[0]
    return K, kR


def vdopp(T, mu_abs, xi):
    return np.sqrt(2*L.KB*T/(mu_abs*L.MH) + xi**2)


def esc(t):
    t = np.asarray(t, float)
    out = np.ones_like(t)
    big = t > 1e-6
    out[big] = -np.expm1(-t[big])/t[big]
    out[~big] = 1 - 0.5*t[~big]
    return out


def flux_frac(st, r, rho, K, flux, r_ph=None):
    n = len(r)
    fb = np.tile(L.fband(st['Teff'])[0], (n, 1))
    if flux == 'A':
        return fb
    r_ph = r[0] if r_ph is None else r_ph
    chi = K[:, :, 0]*rho[:, None]                     # (n, nb) window opacity per cm
    dt = np.zeros((n, L.NB))
    dt[1:] = 0.5*(chi[1:] + chi[:-1])*np.diff(r)[:, None]
    above = (r > r_ph)
    dt[~above] = 0.
    # the interval straddling r_ph: keep only the part above r_ph
    i0 = np.argmax(above) if above.any() else n
    if 0 < i0 < n:
        frac = (r[i0] - r_ph)/(r[i0] - r[i0-1])
        dt[i0] *= frac
    tau = np.cumsum(dt, axis=0)
    f = fb*np.exp(-tau)
    return f/f.sum(1, keepdims=True)


def gamma_sob(st, r, rho, T, dvdr, mode='clamp', xi=2e5, mu_abs=20., flux='A', K=None,
              r_ph=None):
    r, rho, T = np.asarray(r, float), np.asarray(rho, float), np.asarray(T, float)
    dvdr = np.abs(np.broadcast_to(np.asarray(dvdr, float), r.shape))
    if K is None:
        K = opacities(T, rho, mode)[0]
    vD = vdopp(T, mu_abs, xi)
    with np.errstate(divide='ignore', invalid='ignore'):
        ls = np.where(dvdr > 0, vD/dvdr, np.inf)      # Sobolev length [cm]
    tauS = K*(rho*ls)[:, None, None]
    e = np.where(np.isfinite(tauS), esc(np.where(np.isfinite(tauS), tauS, 1.)), 0.)
    kb = (L.GW*K*e).sum(2)                            # (n, nb)
    Fb = flux_frac(st, r, rho, K, flux, r_ph)
    contrib = Fb*kb/st['kE']
    return contrib.sum(1), contrib


def gamma_thin(st, K):
    return (L.fband(st['Teff'])[0]*(L.GW*K).sum(2)).sum(1)/st['kE']


def force_multiplier(st, G, rho, T, dvdr, mu_abs=20., xi=2e5):
    """CAK M(t) = Gamma_Sob/Gamma_e, t = SIGREF rho v_D/|dv/dr|, Gamma_e = SIGREF/kappa_Edd."""
    t = SIGREF*rho*vdopp(T, mu_abs, xi)/np.abs(dvdr)
    return G*st['kE']/SIGREF, t


def lucy_tau(st, r, rho, kR):
    import rt
    return rt.lucy_T(st, r, rho, kR)


def r_photo(r, tauR):
    """tau_R = 2/3 radius (tauR measured from the top, decreasing in r); r[0] if thin."""
    if tauR[0] < 2/3.:
        return r[0]
    i = np.where(tauR >= 2/3.)[0].max()
    if i >= len(r) - 1:
        return r[-1]
    lt = np.log(tauR[i:i+2] + 1e-300)
    return np.interp(np.log(2/3.), lt[::-1], r[i:i+2][::-1])


def beta_wind(st, mode, beta, vinf, mdot, v0=3e5, n=300, x0=1.001, npass=3):
    """beta-law wind 1.001..10 R; grey Lucy T with the ck kappa_R (3 passes)."""
    R = st['R']
    r = R*np.logspace(np.log10(x0), 1., n)
    x = 1 - R/r
    v = v0 + (vinf - v0)*x**beta
    dvdr = (vinf - v0)*beta*x**(beta - 1)*R/r**2
    rho = mdot/(4*np.pi*r**2*v)
    W = 0.5*(1 - np.sqrt(np.maximum(1 - (R/r)**2, 0)))
    T = st['Teff']*W**0.25
    for _ in range(npass):
        K, kR = opacities(T, rho, mode)
        T, tau = lucy_tau(st, r, rho, kR)
    K, kR = opacities(T, rho, mode)
    dtr = np.zeros(n)
    dtr[:-1] = 0.5*(kR[1:]*rho[1:] + kR[:-1]*rho[:-1])*np.diff(r)
    tauR = np.cumsum(dtr[::-1])[::-1]
    return dict(r=r, v=v, dvdr=dvdr, rho=rho, T=T, K=K, kR=kR, tauR=tauR)
