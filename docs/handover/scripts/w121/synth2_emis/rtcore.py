"""rtcore.py: shared constants, Planck helpers, geometry, and the angle-dependent formal solution.

Radiative transfer (all arms): per column, along the LOCAL VERTICAL (plane-parallel slant
path, no curvature of the ray), absorption only (no scattering), LTE source B(T), source
linear in tau within each cell (face values interpolated linearly in tau between cell
centres), I_in = S at the start face (optically thick; the start face is the deepest face
where every chain has vertical tau_from_top > TAUSTART).  Each layer's emission is weighted by
r^2 of that layer (the transported quantity is J = r^2 I), so that the column luminosity
behaves like the GCM's spherical two-stream (ck_spherical); J/r_top^2 is the intensity at
the top face.  Disk flux x d^2 = sum_col J_top(mu) mu dOmega.
"""
import numpy as np

h, c, kB, sig = 6.62607015e-27, 2.99792458e10, 1.380649e-16, 5.670374419e-5
RJ, RSUN = 7.1492e9, 6.957e10
RP = 1.742*RJ              # Sing+2024 (the input's R_p), as synth.py
RS = 1.461*RSUN            # Sing+2024
TS = 6628.0                # Sing+2024 T_eff
LOGG_S, FEH_S = 4.24, 0.13
AOR = 3.7844
TAUSTART = 40.0
SD = '/viper/ptmp2/jinma/w121prod_0929/synth2_emis/'


def planck_lam(lam_cm, T):
    """B_lambda [erg/s/cm2/sr/cm]; lam (..., nl) broadcast with T."""
    x = h*c/(lam_cm*kB*T)
    return 2*h*c*c/lam_cm**5/np.expm1(np.minimum(x, 700.0))


class BandPlanck:
    """int_band B_lambda dlambda as a function of T (table in log T, fine grid in lambda)."""

    def __init__(self, edges_um, Tmin=50.0, Tmax=2.0e4, nT=3000, nl=2000):
        self.lT = np.linspace(np.log10(Tmin), np.log10(Tmax), nT)
        T = 10**self.lT
        e = np.sort(np.asarray(edges_um, float))
        tab = np.zeros((nT, len(e) - 1))
        for b in range(len(e) - 1):
            lam = np.geomspace(e[b], e[b + 1], nl)*1e-4
            B = planck_lam(lam[None, :], T[:, None])
            tab[:, b] = np.trapezoid(B, lam, axis=1)
        self.tab = tab            # ascending-wavelength band order
        self.edges = e

    def __call__(self, T):
        x = (np.log10(np.clip(T, 10**self.lT[0], 10**self.lT[-1])) - self.lT[0])/(
            self.lT[1] - self.lT[0])
        i = np.clip(x.astype(int), 0, len(self.lT) - 2)
        f = (x - i)[..., None]
        return self.tab[i]*(1 - f) + self.tab[i + 1]*f


def gauss01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5*(x + 1), 0.5*w


def phases(nph=72):
    phi = np.arange(nph)/nph
    lo = np.pi - 2*np.pi*phi                     # sub-observer longitude (as synth.py)
    obs = np.stack([-np.cos(lo), -np.sin(lo), 0*lo], -1)
    return phi, obs


def face_source(S, dtau):
    """S, dtau: (nlay, ...) cell values -> (nlay+1, ...) face values (linear in tau)."""
    shp = np.broadcast_shapes(S.shape[1:], dtau.shape[1:])
    Sf = np.empty((S.shape[0] + 1,) + shp, np.result_type(S, dtau))
    Sf[0] = S[0]
    Sf[-1] = S[-1]
    a, b = dtau[:-1].astype(np.float64), dtau[1:].astype(np.float64)
    ab = a + b
    Sf[1:-1] = np.where(ab > 0, (S[:-1]*b + S[1:]*a)/np.maximum(ab, 1e-300),
                        0.5*(S[:-1] + S[1:]))
    return Sf


def formal(dtau, Sf, mu, taustart=TAUSTART):
    """Emergent J at the top face.
    dtau (nlay, *ch) vertical cell optical depths; Sf (nlay+1, *chS) face sources broadcastable
    to (nlay+1, *ch); mu (nmu,) -> J (*ch, nmu)."""
    nlay = dtau.shape[0]
    ttop = np.cumsum(dtau[::-1], axis=0)[::-1]          # tau from the top to the lower face
    deep = ttop.reshape(nlay, -1).min(axis=1) > taustart
    s = int(np.nonzero(deep)[0].max()) if deep.any() else 0
    shp = np.broadcast_shapes(dtau.shape[1:], Sf.shape[1:])
    J = np.broadcast_to(Sf[s][..., None], shp + (len(mu),)).astype(np.float64).copy()
    imu = 1.0/np.asarray(mu, float)
    for l in range(s, nlay):
        t = dtau[l][..., None]*imu
        e = np.exp(-t)
        a = np.where(t > 1e-4, -np.expm1(-t)/np.maximum(t, 1e-30), 1 - 0.5*t)
        J *= e
        J += Sf[l + 1][..., None]*(1 - a) + Sf[l][..., None]*(a - e)
    return J


def formal32(dtau, Sf, mu, taustart=TAUSTART):
    """formal() in float32 with in-place arithmetic (same scheme; 2nd-order Taylor below t=1e-2)."""
    nlay = dtau.shape[0]
    ttop = np.cumsum(dtau[::-1].astype(np.float64), axis=0)[::-1]
    deep = ttop.reshape(nlay, -1).min(axis=1) > taustart
    s = int(np.nonzero(deep)[0].max()) if deep.any() else 0
    shp = np.broadcast_shapes(dtau.shape[1:], Sf.shape[1:]) + (len(mu),)
    Sf = Sf.astype(np.float32)
    J = np.empty(shp, np.float32)
    J[...] = Sf[s][..., None]
    imu = (1.0/np.asarray(mu, float)).astype(np.float32)
    t = np.empty(shp, np.float32)
    e = np.empty(shp, np.float32)
    a = np.empty(shp, np.float32)
    tmp = np.empty(shp, np.float32)
    for l in range(s, nlay):
        np.multiply(dtau[l].astype(np.float32)[..., None], imu, out=t)
        np.negative(t, out=e)
        np.exp(e, out=e)
        # a = (1 - e)/t, Taylor 1 - t/2 + t^2/6 for t < 1e-2
        np.subtract(1.0, e, out=a)
        np.maximum(t, 1e-30, out=tmp)
        np.divide(a, tmp, out=a)
        small = t < 1e-2
        if small.any():
            a[small] = 1 - 0.5*t[small] + t[small]**2/6
        J *= e
        s0 = Sf[l][..., None]
        s1 = Sf[l + 1][..., None]
        np.subtract(1.0, a, out=tmp)
        tmp *= s1
        J += tmp
        np.subtract(a, e, out=tmp)
        tmp *= s0
        J += tmp
    return J
