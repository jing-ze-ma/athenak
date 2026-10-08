"""[betelgeuse26 COPY: adds set C + profile dent26 (Dent+2026 SEM), see the end]
Betelgeuse observed extended-atmosphere profiles and stellar parameter sets.

Density: Dent et al. 2024 (arXiv:2404.06501) Fig. 3 (updated Harper model), digitized
  n_H(r/R*), R* = photospheric (42.5 mas) = 1014 Rsun at 222 pc, 1.08-5 R*; smoothed by a
  cubic spline in (ln x, ln n); beyond 5 R* extrapolated as n ~ x^-2.85 (the observed
  slope beyond ~2 R*, section 7 of the note).
  Harper, Brown & Lim 2001 Table 7 (second profile): their R* = 56 mas = 1.32 x the
  photospheric R*, d = 131 pc -> x_phot = 1.32 x_H, n_H x (131/d)^0.5; used for x >= 1.5
  (inside that the table is their 11 um 'photosphere'); runs to 52.8 R*.
Temperature: Dent Fig. 3 read off at 1.2/1.5/2/3/4/5 R* (2570/3040/3800/3290/2660/2260 K;
  section 7 of the note), held at 2570 K inside 1.2 R*; beyond 5 R* Harper's T(x/1.32)
  scaled to 2260 K at 5 R*. Harper profile uses Harper's own T(x/1.32).
rho = n_H x 2.27e-24 g (brief).
Set B (Joyce et al. 2020, 168 pc, 764 Rsun): same angular radius (764 Rsun at 168 pc =
  42.3 mas), so x is unchanged and n_H scales as d^-1/2 -> x (222/168)^0.5.
"""
import os
import numpy as np
from scipy.interpolate import UnivariateSpline, PchipInterpolator

G, C, SIG, KB, MH = 6.674e-8, 2.99792458e10, 5.670374419e-5, 1.380649e-16, 1.6726e-24
MSUN, LSUN, RSUN, YR = 1.989e33, 3.828e33, 6.957e10, 3.156e7
MPH = 2.27e-24          # g per H
DS = '/orion/ptmp/jinma/rsg_wind_1008/dens_survey/'
SLOPE = -2.85

# Harper, Brown & Lim 2001 Table 7 (r/R*_H, T_e [K], log n_H) -- transcribed from h01_tb7.png
H01 = np.array([
    [1.005, 2973, 14.82], [1.01, 2720, 14.44], [1.02, 2637, 14.23], [1.03, 2563, 14.02],
    [1.04, 2499, 13.81], [1.05, 2448, 13.61], [1.06, 2414, 13.40], [1.08, 2400, 12.99],
    [1.09, 2418, 12.78], [1.10, 2450, 12.58], [1.15, 2744, 11.59], [1.20, 3043, 10.78],
    [1.25, 3307, 10.38], [1.30, 3529, 10.22], [1.35, 3681, 10.12], [1.40, 3768, 10.04],
    [1.45, 3811, 9.97], [1.50, 3805, 9.90], [1.55, 3788, 9.84], [1.60, 3771, 9.78],
    [1.65, 3750, 9.73], [1.70, 3727, 9.68], [1.80, 3688, 9.58], [1.90, 3584, 9.50],
    [2.00, 3447, 9.42], [2.20, 3221, 9.29], [2.40, 3017, 9.17], [2.70, 2764, 9.02],
    [3.00, 2569, 8.89], [3.50, 2295, 8.71], [4.00, 2072, 8.56], [4.50, 1884, 8.43],
    [5.00, 1732, 8.32], [5.50, 1600, 8.22], [6.00, 1453, 8.12], [6.50, 1330, 8.04],
    [7.00, 1225, 7.96], [8.00, 1056, 7.83], [9.00, 926, 7.71], [10.0, 824, 7.61],
    [15.0, 525, 7.21], [20.0, 382, 6.94], [25.0, 298, 6.73], [30.0, 243, 6.56],
    [35.0, 205, 6.42], [40.0, 177, 6.29]])
FH = 1.32
DENT_T = np.array([[1.2, 2570.], [1.5, 3040.], [2., 3800.], [3., 3290.], [4., 2660.],
                   [5., 2260.]])

SETS = {
    # name: d [pc], R [Rsun], Teff, L [Lsun] (None -> 4 pi R^2 sigma Teff^4), masses
    'A': dict(d=222., R=1014., Teff=3650., L=None, M=(18., 20.),
              src='d 222 pc, R 1014 Rsun (Dent+2024 convention), Teff 3650 K, '
                  'L = 4 pi R^2 sigma Teff^4, M 18-20'),
    'C': dict(d=172., R=812., Teff=3650., L=None, M=(17.5,),
              src='Dent+2026 SEM: d 172 pc, R 812 Rsun (Rosseland tau=2/3), Teff 3650 K '
                  '(their MARCS model), L = 4 pi R^2 sigma Teff^4, M 17.5'),
    'B': dict(d=168., R=764., Teff=None, L=1e5, M=(16.5, 19.),
              src='Joyce+2020: d 168 pc, R 764 Rsun, L 1e5 Lsun, M 16.5-19'),
}


def stellar(name, M):
    s = SETS[name]
    R = s['R'] * RSUN
    if s['L'] is None:
        L = 4*np.pi*R**2*SIG*s['Teff']**4
        Te = s['Teff']
    else:
        L = s['L'] * LSUN
        Te = (L/(4*np.pi*R**2*SIG))**0.25
    Mg = M*MSUN
    return dict(name=f'{name}{M:g}', set=name, M=Mg, L=L, R=R, Teff=Te, d=s['d'],
                kE=4*np.pi*C*G*Mg/L, g=G*Mg/R**2, vesc=np.sqrt(2*G*Mg/R))


def _dent_raw():
    d = np.loadtxt(DS + 'dent24_fig3_digitized.txt')
    o = np.argsort(d[:, 0])
    x, ln = d[o, 0], d[o, 1]*np.log(10.)
    keep = np.concatenate([[True], np.diff(x) > 1e-6])
    return x[keep], ln[keep]


_xd, _lnd = _dent_raw()
_spl = UnivariateSpline(np.log(_xd), _lnd, k=3, s=len(_xd)*(0.02*np.log(10))**2)
X5 = 5.0
_ln5 = float(_spl(np.log(X5)))
_slope5 = float(_spl.derivative()(np.log(X5)))

# temperature: Dent points + Harper beyond 5 R* scaled to match at 5 R*
_TH = PchipInterpolator(np.log(FH*H01[:, 0]), np.log(H01[:, 1]))
_fT = DENT_T[-1, 1]/np.exp(_TH(np.log(X5)))
_xo = FH*H01[:, 0]
_sel = _xo > X5*1.04
_Tx = np.concatenate([DENT_T[:, 0], _xo[_sel]])
_Ty = np.concatenate([DENT_T[:, 1], H01[_sel, 1]*_fT])
_TD = PchipInterpolator(np.log(_Tx), np.log(_Ty))
_HN = PchipInterpolator(np.log(FH*H01[9:, 0]), H01[9:, 2]*np.log(10.))   # from r_H=1.10


def nH(x, prof='dent', d=222.):
    """n_H [cm^-3] at x = r/R*(photospheric)."""
    x = np.atleast_1d(np.asarray(x, float))
    if prof == 'dent':
        lx = np.log(np.clip(x, _xd[0], None))
        ln = np.where(x <= X5, _spl(np.minimum(lx, np.log(X5))),
                      _ln5 + SLOPE*(lx - np.log(X5)))
        return np.exp(ln)*np.sqrt(222./d)
    if prof == 'dent26':
        return nH26(x)*NSCALE
    if prof == 'harper':
        lx = np.log(np.clip(x, 1.5, FH*H01[-1, 0]))
        return np.exp(_HN(lx))*np.sqrt(131./d)
    raise ValueError(prof)


def Tobs(x, prof='dent'):
    x = np.atleast_1d(np.asarray(x, float))
    if prof == 'dent26':
        return T26(x)
    if prof == 'dent':
        return np.exp(_TD(np.log(np.clip(x, 1.2, FH*H01[-1, 0]))))
    return np.exp(_TH(np.log(np.clip(x, FH*1.10, FH*H01[-1, 0]))))


def rho(x, prof='dent', d=222.):
    return nH(x, prof, d)*MPH


def dlog(f, x, h=1e-3):
    """d ln f / d ln x numerically (f callable of x)."""
    return (np.log(f(x*np.exp(h))) - np.log(f(x*np.exp(-h))))/(2*h)


def ft_fit(st, prof='dent', x0=1.2, x1=5.0):
    """Fit FT eq. 7  n = n_ph x^-2 exp(-q sqrt(1-1/x)) to the observed profile on [x0, x1];
    returns q (= v_esc/v_con), v_con, rho_ph, Mdot(>v_esc) = 4 pi R^2 rho_ph v_con e^-q."""
    x = np.exp(np.linspace(np.log(x0), np.log(x1), 200))
    y = np.log(nH(x, prof, st['d'])) + 2*np.log(x)
    A = np.vstack([np.ones_like(x), -np.sqrt(1 - 1/x)]).T
    (lnph, q), res, *_ = np.linalg.lstsq(A, y, rcond=None)
    rms = np.sqrt(np.mean((A @ np.array([lnph, q]) - y)**2))
    rho_ph = np.exp(lnph)*MPH
    vcon = st['vesc']/q
    md = 4*np.pi*st['R']**2*rho_ph*vcon*np.exp(-q)
    return dict(q=q, vcon=vcon, rho_ph=rho_ph, mdot_esc=md*YR/MSUN, rms_dex=rms/np.log(10))


# ---------------------------------------------------------------- Dent et al. 2026 SEM
# arXiv:2608.19339v2 App. B eq. (1): n_H = n_chrom exp[-(x - x_min)/h] + n_* x^-2
#   (0.998 - x^-gamma)^-delta, x = R/R*, R* = 812 Rsun; inside x_min the exponential is
#   continued inward (only used where the grey hydrostatic photosphere is joined, gas.py).
# T_e: Fig. 5 digitized (digitize_fig5.py -> out/dent26_fig5_T.txt, +-50 K), 1.13-5 R*;
#   the V-shaped minimum is set to the column minimum (1776 K at 1.19-1.20; text: ~1700 K).
#   Beyond 5 R* (env BET26_TEXT): 'harper' (default) = Harper+2001 T(x/1.32) shape scaled to
#   1412 K at 5 R*; 'plaw' = T ~ x^-0.67 (log slope of Fig. 5 over 4.6-5 R*); 'const' = 1412 K.
NCH, XMIN, HCH, NST, GAM, DEL = 1.34e12, 1.144, 0.06, 2.96e9, 0.45, 1.5
NSCALE = float(os.environ.get('BET26_NSCALE', 1.0))
TEXT = os.environ.get('BET26_TEXT', 'harper')
_t26 = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'out/dent26_fig5_T.txt'))
_t26x, _t26T = _t26[:, 0], _t26[:, 1].copy()
_dip = (_t26x > 1.185) & (_t26x < 1.205)
_t26T[_dip] = _t26[_dip, 2].min()
T5_26 = float(_t26T[-1])
PL26 = -0.67


def nH26(x):
    x = np.atleast_1d(np.asarray(x, float))
    return NCH*np.exp(-(x - XMIN)/HCH) + nH26_wind(x)


def nH26_wind(x):
    x = np.atleast_1d(np.asarray(x, float))
    xx = np.maximum(x, XMIN)        # the wind term is singular at 1.0045 R*; frozen inside x_min
    return NST*xx**-2*(0.998 - xx**-GAM)**-DEL


def v26(x, vinf=9e5):
    """velocity implied by the SEM wind term for a steady wind with v -> vinf."""
    x = np.atleast_1d(np.asarray(x, float))
    return vinf*(np.clip(0.998 - x**-GAM, 0, None)/0.998)**DEL


def T26(x, mode=None):
    mode = TEXT if mode is None else mode
    x = np.atleast_1d(np.asarray(x, float))
    Ti = np.interp(x, _t26x, _t26T)
    if mode == 'harper':
        Te = T5_26*np.exp(_TH(np.log(np.clip(x, 5., FH*H01[-1, 0]))) - _TH(np.log(5.)))
    elif mode == 'plaw':
        Te = T5_26*(np.maximum(x, 5.)/5.)**PL26
    elif mode == 'const':
        Te = np.full_like(x, T5_26)
    else:
        raise ValueError(mode)
    return np.where(x <= 5., Ti, Te)
