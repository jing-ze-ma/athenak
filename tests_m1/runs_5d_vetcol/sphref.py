"""Independent 1-D spherical grey formal solution (long characteristics) for T-S6/T-S4.

Continuum problem: chi(r) = C / r^2 (C = rho0 kappa_t r_in^2), S(r) given as a function,
inner sphere r_in with the diffusion intensity I = E_b + 3 F_b mu_f / c, vacuum at r_top.
Each (r, mu) ray is integrated exactly in geometry: tau(s) = (C/p) atan(s/p) along the
line of impact parameter p, the source on a grid uniform in tau (Simpson with the
exponential kernel), started from the thermalised intensity S + dS/dtau where the ray is
more than TMAX deep.  Angles: Gauss-Legendre in mu on [mu_c, 1] (core rays) and [0, mu_c]
(rays with p > r_in) outgoing, [0, 1] incoming.  No p-ray grid, no short characteristics:
nothing is shared with rad_m1_vetcol.cpp but the problem definition.
"""
import numpy as np

TMAX = 60.0


def _tau(s, p, C):
    if p < 1e-14:
        return C*(-1.0/s)  # p = 0 is never used (Gauss nodes)
    return (C/p)*np.arctan(s/p)


def _s_of_tau(t, p, C):
    return p*np.tan(p*t/C)


def ray_intensity(sa, sb, p, C, Sfun, Ia, nt=2001):
    """I at s = sb for a ray entering at s = sa (< sb) with intensity Ia."""
    ta, tb = _tau(sa, p, C), _tau(sb, p, C)
    ttot = tb - ta
    therm = False
    if C == 0.0 or ttot < 1e-12:
        return Ia
    if ttot > TMAX:
        ta = tb - TMAX
        therm = True
    t = np.linspace(ta, tb, nt)
    s = _s_of_tau(t, p, C)
    r = np.sqrt(p*p + s*s)
    S = Sfun(r)
    ker = np.exp(-(tb - t))
    h = t[1] - t[0]
    w = np.ones(nt)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    integ = h/3.0*np.sum(w*S*ker)
    if therm:
        Ia = S[0] + (S[1] - S[0])/h
    return Ia*np.exp(-(tb - ta)) + integ


def moments(r, rin, rtop, C, Sfun, Eb, Fb, c, ng=48, nt=2001):
    """J, H, K (E units) at radius r."""
    xg, wg = np.polynomial.legendre.leggauss(ng)
    xg = 0.5*(xg + 1.0)
    wg = 0.5*wg
    muc = np.sqrt(max(1.0 - (rin/r)**2, 0.0))
    J = H = K = 0.0
    ztop = lambda p: np.sqrt(rtop*rtop - p*p)
    # outgoing core rays mu in [muc, 1]
    for x, w in zip(xg, wg):
        mu = muc + (1.0 - muc)*x
        ww = (1.0 - muc)*w
        p = r*np.sqrt(1.0 - mu*mu)
        zf = np.sqrt(max(rin*rin - p*p, 0.0))
        muf = zf/rin
        I0 = max(Eb + 3.0*Fb*muf/c, 0.0)
        I = ray_intensity(zf, r*mu, p, C, Sfun, I0, nt)
        J += 0.5*ww*I
        H += 0.5*ww*mu*I
        K += 0.5*ww*mu*mu*I
    # outgoing rays that pass their tangent point, mu in [0, muc]
    for x, w in zip(xg, wg):
        mu = muc*x
        ww = muc*w
        p = r*np.sqrt(1.0 - mu*mu)
        I = ray_intensity(-ztop(p), r*mu, p, C, Sfun, 0.0, nt)
        J += 0.5*ww*I
        H += 0.5*ww*mu*I
        K += 0.5*ww*mu*mu*I
    # incoming rays mu in [0, 1]
    for x, w in zip(xg, wg):
        mu = x
        p = r*np.sqrt(1.0 - mu*mu)
        I = ray_intensity(-ztop(p), -r*mu, p, C, Sfun, 0.0, nt)
        J += 0.5*w*I
        H -= 0.5*w*mu*I
        K += 0.5*w*mu*mu*I
    return J, H, K


def eddington_E(r, rho0, kt, fin, rin, rout, c, q, n=2.0):
    """the pgen's atm_init = eddington E(r) (rad_m1_tests2.cpp, sph_atm)"""
    return (3.0*rho0*kt*fin*rin**(n + 2.0)/(c*(n + 1.0))
            * (r**(-(n + 1.0)) - rout**(-(n + 1.0))) + fin*rin*rin/(rout*rout*c*q))
