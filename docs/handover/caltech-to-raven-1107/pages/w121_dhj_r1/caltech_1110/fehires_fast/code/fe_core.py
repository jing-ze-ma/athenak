"""fe_core.py: Fe I line list, partition functions, Fe ionisation, continuum opacities.

All cgs.  Line list: Kurucz gfall gf2600.all (Fe I, lab gf where available; per-line strengths
match pRT 3.4 lbl Fe to 0.3 %, see val_lines.py).  Partition functions: direct level sums over
every level in gf2600.pos
(Fe I) and gf2601.pos (Fe II).  Electron density, H/H2 split and mu: the run's own composition
model (EOSCompositionModel, comptab.cpp -> data/comptab_*.bin.gz; WASP-121b 3x: X 0.7188
Y 0.242 [M/H] 0.4771).  Portable version: all inputs are in DATA (data/ next to this file);
the raw Kurucz files are only needed to REBUILD the cached line list / partition functions
(env FEHIRES_KURUCZ); call configure(cfg) with a planet.load() dict before use.
"""
import gzip
import os
import numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data') + '/'
KUR = os.environ.get('FEHIRES_KURUCZ', 'kurucz') + '/'
PRT_FE = os.environ.get('FEHIRES_PRT_FE', '')   # pRT lbl Fe .h5, validation scripts only
KB = 1.380649e-16
HP = 6.62607015e-27
ME = 9.1093837e-28
MU_ = 1.66053907e-24
CL = 2.99792458e10
EV = 1.602176634e-12
HCK = 1.4387769              # hc/k [cm K]
SIGF = 0.0265400193          # pi e^2/(m_e c) [cm^2 Hz]
M_FE = 55.845*MU_
CHI_FE = 7.9024*EV           # Fe I ionisation energy (NIST ASD value)
XH, YHE = 0.7188, 0.2420     # run composition (eos_xh, eos_yhe)
A_FE_SUN = 3.2e-5            # Fe/H by number, solar: eos_composition.hpp metal_donor table
COMPTAB = DATA + 'comptab_3x.bin.gz'
LAM0, LAM1 = 380.0e-7, 790.0e-7   # ESPRESSO range used [cm]
DLN = 1.0e-6                 # ln-lambda pixel (R = 1e6)
DV = CL*DLN                  # 0.2998 km/s


def configure(cfg):
    """set the run composition from a planet.load() dict (defaults = WASP-121b 3x)."""
    global XH, YHE, A_FE_SUN, COMPTAB
    XH, YHE, A_FE_SUN, COMPTAB = cfg['eos_xh'], cfg['eos_yhe'], cfg['a_fe'], cfg['comptab']


def air2vac(lam_nm):
    """Kurucz gives air wavelengths above 200 nm; IAU (Morton 2000) refractivity, inverted."""
    lam = np.asarray(lam_nm, float)*10.0          # A, air
    vac = lam.copy()
    for _ in range(4):
        s2 = (1e4/vac)**2
        n = 1 + 8.34254e-5 + 2.406147e-2/(130 - s2) + 1.5998e-4/(38.9 - s2)
        vac = lam*n
    return np.where(lam_nm > 200, vac, lam)/10.0


def _parse(fn):
    lam, lgf, e1, j1, e2, j2 = [], [], [], [], [], []
    for ln in open(fn):
        try:
            lam.append(float(ln[0:11]))
            lgf.append(float(ln[11:18]))
            e1.append(abs(float(ln[24:36])))
            j1.append(float(ln[36:41]))
            e2.append(abs(float(ln[52:64])))
            j2.append(float(ln[64:69]))
        except ValueError:
            for q in (lam, lgf, e1, j1, e2, j2):
                if len(q) > len(j2):
                    q.pop()
    return [np.array(q) for q in (lam, lgf, e1, j1, e2, j2)]


def _levels(fn):
    _, _, e1, j1, e2, j2 = _parse(fn)
    E = np.concatenate([e1, e2])
    J = np.concatenate([j1, j2])
    key = np.round(E, 2)*100 + J/100.0
    _, ix = np.unique(key, return_index=True)
    return E[ix], 2*J[ix] + 1


_PF = {}


def partition(T, ion):
    """U(T) by direct level sum (ion 0 = Fe I from gf2600.pos, 1 = Fe II from gf2601.pos)."""
    if ion not in _PF:
        cache = DATA + 'pf_fe%d.npz' % ion
        if os.path.exists(cache):
            d = np.load(cache)
            _PF[ion] = (d['lt'], d['lu'])
        else:
            E, g = _levels(KUR + ('gf2600.pos' if ion == 0 else 'gf2601.pos'))
            lt = np.linspace(np.log10(300.0), np.log10(2e4), 600)
            lu = np.log10([(g*np.exp(-HCK*E/10**t)).sum() for t in lt])
            np.savez(cache, lt=lt, lu=lu, nlev=len(E))
            _PF[ion] = (lt, lu)
    lt, lu = _PF[ion]
    return 10**np.interp(np.log10(T), lt, lu)


def lines(l0=LAM0, l1=LAM1, pad=1e-4):
    """Fe I lines with vacuum wavelengths in [l0, l1] (cm, padded); returns lam [cm], gf,
    E_low [eV]."""
    cache = DATA + 'fe1_lines.npz'
    if os.path.exists(cache):
        d = np.load(cache)
        return d['lam'], d['gf'], d['El']
    lam, lgf, e1, j1, e2, j2 = _parse(KUR + 'gfall_2600.all')
    lv = air2vac(lam)*1e-7
    m = (lv > l0*(1 - pad)) & (lv < l1*(1 + pad))
    El = np.minimum(e1, e2)[m]*HP*CL/EV
    gr, gw = damping(KUR + 'gfall_2600.all')
    gr, gw = gr[m], gw[m]
    # missing radiative damping: classical 2.223e15/lambda_A^2; missing vdW: Unsold-like
    # default log gamma_W = -7.8 (the gfall median for lines that have one)
    gr = np.where(gr > 0, 10**gr, 2.223e15/(lv[m]*1e8)**2)
    gw = np.where(gw < 0, 10**gw, 10**-7.8)
    out = dict(lam=lv[m], gf=10**lgf[m], El=El, gr=gr, gw=gw)
    np.savez(cache, **out)
    return out['lam'], out['gf'], out['El']


def damping(fn):
    """log10 gamma_rad [1/s], log10 gamma_vdW per perturber at 1e4 K [cm^3/s] (gfall cols)."""
    gr, gw = [], []
    for ln in open(fn):
        try:
            float(ln[0:11])
            float(ln[11:18])
            float(ln[24:36])
            float(ln[36:41])
            float(ln[52:64])
            float(ln[64:69])
        except ValueError:
            continue
        try:
            gr.append(float(ln[80:86]))
        except ValueError:
            gr.append(0.0)
        try:
            gw.append(float(ln[92:98]))
        except ValueError:
            gw.append(0.0)
    return np.array(gr), np.array(gw)


def line_damping():
    d = np.load(DATA + 'fe1_lines.npz')
    return d['gr'], d['gw']


class Comp:
    """bilinear (log rho, log T) lookup of the run composition model (comptab.bin)."""
    def __init__(self, fn=None):
        fn = fn or COMPTAB
        with (gzip.open(fn, 'rb') if fn.endswith('.gz') else open(fn, 'rb')) as fh:
            buf = fh.read()
        nx, ny = (int(v) for v in np.frombuffer(buf, np.int32, 2))
        self.x0, self.dx, self.y0, self.dy = np.frombuffer(buf, np.float64, 4, 8)
        d = np.frombuffer(buf, np.float64, offset=40).reshape(ny, nx, 5)
        self.nx, self.ny, self.d = nx, ny, d

    def __call__(self, rho, T):
        x = np.clip((np.log10(rho) - self.x0)/self.dx, 0, self.nx - 1.000001)
        y = np.clip((np.log10(T) - self.y0)/self.dy, 0, self.ny - 1.000001)
        ix, iy = x.astype(int), y.astype(int)
        fx, fy = (x - ix)[..., None], (y - iy)[..., None]
        d = self.d
        v = ((d[iy, ix]*(1 - fx) + d[iy, ix + 1]*fx)*(1 - fy)
             + (d[iy + 1, ix]*(1 - fx) + d[iy + 1, ix + 1]*fx)*fy)
        return {'xe': 10**v[..., 0], 'xh2': v[..., 1], 'xhii': v[..., 2], 'mu': v[..., 3],
                'pspec': 10**v[..., 4]}


def saha_fe(T, ne):
    """n(Fe II)/n(Fe I)."""
    return (2*partition(T, 1)/partition(T, 0)*(2*np.pi*ME*KB*T/HP**2)**1.5
            * np.exp(-CHI_FE/(KB*T))/ne)


def species(rho, T, comp, afe=A_FE_SUN, tcond=0.0):
    """number densities: FeI, Fe(total), e, H, H2, He; plus x(FeI)/U(FeI)."""
    c = comp(rho, T)
    ntot = rho/(c['mu']*MU_)
    ne = c['xe']*ntot
    nht = XH*rho/(1.00794*MU_)
    nH = nht*np.clip(1 - c['xh2'] - c['xhii'], 0, 1)
    nH2 = 0.5*c['xh2']*nht
    nHe = YHE*rho/(4.002602*MU_)
    nfe = afe*nht*(T >= tcond)
    r = saha_fe(T, np.maximum(ne, 1e-30))
    nfe1 = nfe/(1 + r)
    return dict(nfe1=nfe1, nfe=nfe, ne=ne, nH=nH, nH2=nH2, nHe=nHe, ntot=ntot,
                fe1_over_U=nfe1/partition(T, 0))


# --- continuum (absorption + scattering out of the beam), per cm, at wavelengths lam [cm]
_AHM = [1.99654, -1.18267e-5, 2.64243e-6, -4.40524e-10, 3.23992e-14, -1.39568e-18, 2.78701e-23]


def hminus_bf_xsec(lam):
    """Gray (Photospheres) H- bf cross section per H- [cm^2], as pRT."""
    la = lam*1e8
    s = sum(a*la**i for i, a in enumerate(_AHM))*1e-18
    return np.where(la <= 1.64e4, np.maximum(s, 0), 0.0)


def hminus_ff_xsec(lam, T, pe):
    """Gray H- ff, cm^2 per neutral H atom, as pRT (T clamped at >= 2500 K, where pRT sets 0)."""
    la = np.log10(lam*1e8)
    th = np.log10(5040.0/np.maximum(T, 2500.0))
    f0 = -2.2763 - 1.6850*la + 0.76661*la**2 - 0.053346*la**3
    f1 = 15.2827 - 9.2846*la + 1.99381*la**2 - 0.142631*la**3
    f2 = -197.789 + 190.266*la - 67.9775*la**2 + 10.6913*la**3 - 0.625151*la**4
    return 1e-26*pe*10**(f0 + f1*th + f2*th**2)


def rayleigh_xsec(lam):
    """H2 and He from pRT's own Rayleigh routine (per gram -> per particle); atomic H
    (not in pRT) from the Kurucz ATLAS HRAYOP form; cm^2 per particle.
    H2/He at the 9 continuum nodes are cached in data/rayleigh_nodes.npz (pRT 3.4 values), so
    pRT is not needed at run time."""
    la = np.atleast_1d(lam*1e8).astype(float)
    c = np.load(DATA + 'rayleigh_nodes.npz')
    if c['lam'].shape == np.shape(lam) and np.array_equal(c['lam'], lam):
        return c['h2'], c['he'], 5.799e-13/la**4 + 1.422e-6/la**6 + 2.784/la**8
    from petitRADTRANS.fortran_radtrans_core import fortran_radtrans_core as fcore
    one = np.ones(1)

    def prt(sp, m):
        return fcore.compute_rayleigh_scattering_opacities(
            sp, one, la, 2.33*one, 3000*one, 1e-3*one)[:, 0]*m*MU_
    h2 = prt('H2', 2.01588)
    he = prt('He', 4.002602)
    h = 5.799e-13/la**4 + 1.422e-6/la**6 + 2.784/la**8
    return h2, he, h


def alpha_cont(lam, T, sp):
    """continuum extinction coefficient [1/cm], shape (..., nlam)."""
    T = np.asarray(T)[..., None]
    nH, ne, nH2, nHe = (np.asarray(sp[k])[..., None] for k in ('nH', 'ne', 'nH2', 'nHe'))
    nhm = nH*ne*0.25*(HP**2/(2*np.pi*ME*KB*T))**1.5*np.exp(0.754*EV/(KB*T))
    h2, he, h = rayleigh_xsec(lam)
    return (nhm*hminus_bf_xsec(lam) + nH*hminus_ff_xsec(lam, T, ne*KB*T)
            + nH2*h2 + nHe*he + nH*h)


def hipass(y, w=1001):
    """continuum removal as for normalised data: subtract a smoothed running minimum over
    w pixels (w = 1001 px = 300 km/s)."""
    from scipy.ndimage import minimum_filter1d, uniform_filter1d
    base = uniform_filter1d(minimum_filter1d(y, w, mode='nearest'), w, mode='nearest')
    return y - base


def template_hp(fwhm_kms):
    """CCF template: pRT 3.4 lbl Fe I transit depth (isothermal 3000 K, solar Fe, all neutral;
    template.py -> data/template_raw.npz 'raw'), continuum-removed (hipass) and broadened by a
    Gaussian of FWHM fwhm_kms.  For 4.28 km/s this is bit-identical to the viper template.npz
    'tbhp' used on 10-07."""
    raw = hipass(np.load(DATA + 'template_raw.npz')['raw'])
    return np.convolve(raw, lsf_kernel(fwhm_kms), mode='same')


def lsf_kernel(fwhm_kms):
    sig = fwhm_kms*1e5/2.3548/DV
    k = np.arange(-int(5*sig), int(5*sig) + 1)
    ker = np.exp(-0.5*(k/sig)**2)
    return ker/ker.sum()
