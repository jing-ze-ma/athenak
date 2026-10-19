"""Microphysics of the standard IC column (tools/ic_star).

Numerically IDENTICAL (same expressions, same operation order) to the AG Car builders:
  * TableEOS  = make_ic_mlt_star_ge.TableEOS (Viper lbv_1008/agcar/geos/scripts; the code's
    gas-only general-EOS table dump, bicubic RectBivariateSpline of log10 e/rho and log10 p/rho,
    (log T, log rho) CLAMPED to the dumped grid, clamp counter);
  * read_ross = the he_star_m1 table reader + the pgen's he_opac_logd_min extension
    (log kappa continued below the grid with the first interval's slope clamped to [0, 1]);
    bilinear (kx = ky = 1) in (log T, log rho), as the pgen lookup;
  * Star.thermo / Star.mlt_gradient = make_ic_mlt_star_ge.Star (gas EOS + radiation, local
    Bohm-Vitense MLT, MLT++ of Paxton+2013 sect. 7.2 / MESA r15140 adjust_gradT_excess).

MLT++ forms (identical for pp_alpha = 1 up to the last bit of the arithmetic order):
  'mesa' (standard): af = f2 + (1 - f2)(1 - pp_alpha); nabla <- af nabla + (1 - af) nabla_ad
          where nabla - nabla_ad > f1 (make_ic_agcar_ge.py, AG Car A: pp_alpha 0.970755);
  'f2'   (make_ic_agcar_mltpp.py StarPP, AG Car B repro and the hg145 ic2 build):
          nabla <- nabla_ad + f2 (nabla - nabla_ad) where nabla_rad > nabla_ad and
          nabla - nabla_ad > f1  (pp_alpha must be 1).
"""
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import brentq

CL = 2.99792458e10
SIG = 5.670374419e-5
A_RAD_ = 7.5657332503e-15
KB = 1.3806488e-16          # = Units::k_boltzmann_cgs
MU_ = 1.660538921e-24       # = Units::atomic_mass_unit_cgs
G_CODE = 6.67408e-8         # = Units::grav_constant_cgs
RSUN = 6.957e10
LSUN = 3.828e33
GMSUN = 1.3271244e26        # the AG Car builders' G Msun
MSUN = 1.98847e33


class TableEOS:
    """the run's gas-only table dump; (log T, log rho) clamped to the grid (derivative 0
    across a clamped axis); counts calls and clamped calls (coverage report)."""
    def __init__(self, fn, verbose=True, clamp=True):
        with open(fn) as fh:
            lines = [fh.readline() for _ in range(6)]
        nx, ny, xmin, dx, ymin, dy = [float(v) for v in lines[2].split()[1:]]
        nx, ny = int(nx), int(ny)
        d = np.loadtxt(fn, comments='#')
        x = xmin + dx*np.arange(nx)
        y = ymin + dy*np.arange(ny)
        self.xr, self.yr = (x[0], x[-1]), (y[0], y[-1])
        self.sle = RectBivariateSpline(y, x, d[:, 0].reshape(ny, nx), kx=3, ky=3)
        self.slp = RectBivariateSpline(y, x, d[:, 1].reshape(ny, nx), kx=3, ky=3)
        self.desc = 'table ' + fn
        self.fn = fn
        self.ncall = 0
        self.nclamp = 0
        self.verbose = verbose
        self.clamp = clamp      # False: plain spline (extrapolates; make_ic_mlt_star.py)
        self.comp = lines[4].split()[1:] if lines[3].startswith('# comp') else None

    def _ev(self, sp, lt, lr, dx, dy):
        self.ncall += 1
        if self.verbose and self.ncall % 200000 == 0:
            print('    TableEOS: %d calls, %d clamped; last log T %s log rho %s'
                  % (self.ncall, self.nclamp, np.ravel(lt)[:1], np.ravel(lr)[:1]), flush=True)
            if self.ncall > 2.0e8:
                raise RuntimeError('TableEOS: more than 2e8 calls')
        lt2 = np.clip(lt, *self.yr)
        lr2 = np.clip(lr, *self.xr)
        cl = np.any(lt2 != lt) or np.any(lr2 != lr)
        if cl:
            self.nclamp += 1
        if not self.clamp:
            return sp.ev(lt, lr, dx=dx, dy=dy)
        v = sp.ev(lt2, lr2, dx=dx, dy=dy)
        if (dx and np.any(lt2 != lt)) or (dy and np.any(lr2 != lr)):
            v = np.where((lt2 != lt) if dx else (lr2 != lr), 0.0, v)
        return v

    def lp(self, lt, lr, dx=0, dy=0):       # log10 p/rho; dx: d/dlog T, dy: d/dlog rho
        return self._ev(self.slp, lt, lr, dx, dy)

    def le(self, lt, lr, dx=0, dy=0):       # log10 e/rho (specific)
        return self._ev(self.sle, lt, lr, dx, dy)


def read_ross_raw(fn):
    """(nT, nD, lT0, dlT, lD0, dlD, values[nT, nD]) of a he_star_m1 table"""
    with open(fn) as fh:
        for ln in fh:
            if ln.startswith('# ') and len(ln.split()) == 7 and ln.split()[1].isdigit():
                nT, nD, lT0, dlT, lD0, dlD = [float(v) for v in ln.split()[1:]]
    nT, nD = int(nT), int(nD)
    v = np.loadtxt(fn, comments='#').reshape(nT, nD)
    return nT, nD, lT0, dlT, lD0, dlD, v


def read_ross(fn, ld_ext=None):
    """he_star_m1 table format; ld_ext: HsReadOpacityTable's extension to low rho
    (log10 kappa continued in log rho with the first interval's slope clamped to [0, 1])"""
    nT, nD, lT0, dlT, lD0, dlD, v = read_ross_raw(fn)
    if ld_ext is not None and ld_ext < lD0 - 1e-9*dlD:
        nadd = int(np.ceil((lD0 - ld_ext)/dlD - 1e-9))
        sl = np.clip((v[:, 1] - v[:, 0])/dlD, 0.0, 1.0)
        ext = v[:, :1] - sl[:, None]*(nadd - np.arange(nadd))[None, :]*dlD
        v = np.concatenate([ext, v], axis=1)
        nD += nadd
        lD0 -= nadd*dlD
    return lT0 + dlT*np.arange(nT), lD0 + dlD*np.arange(nD), v


def q_hopf(tau):
    """Hopf function, Mihalas' fit (|err| < 1e-3)"""
    return 0.7104 - 0.1331*np.exp(-3.4488*tau)


class Star:
    """GM, L, R_ph, EOS, kappa_R and the MLT/MLT++ controls of one star"""
    def __init__(self, gm, lum, rstar, eos_file, ross, ld_ext, alpha, pp_alpha, f1, f2,
                 pp_form='mesa', verbose=True, eos_clamp=True, mlt_option='KW',
                 alt_hp=False):
        self.GM, self.LUM, self.RSTAR = gm, lum, rstar
        self.eos = TableEOS(eos_file, verbose, eos_clamp)
        self.ross_lims = None
        ty, tx, kv = read_ross(ross, ld_ext)
        self.ross_lims = (ty[0], ty[-1], tx[0], tx[-1])
        self.sk = RectBivariateSpline(ty, tx, kv, kx=1, ky=1)
        self.ALPHA = alpha
        self.PP_ALPHA = pp_alpha
        self.PP_F1 = f1
        self.PP_F2 = f2
        self.PP_MAXLOGT = 8.0       # MESA gradT_excess_max_logT default (controls.defaults:2399)
        if pp_form not in ('mesa', 'f2'):
            raise ValueError('pp_form must be mesa or f2')
        if pp_form == 'f2' and pp_alpha not in (0.0, 1.0):
            raise ValueError("pp_form 'f2' is the pp_alpha = 1 rule (or 0 = off)")
        self.PP_FORM = pp_form
        if mlt_option not in ('KW', 'Cox'):
            raise ValueError("mlt_option must be 'Cox' (MESA default) or 'KW' (legacy)")
        self.MLT_OPTION = mlt_option
        self.ALT_HP = alt_hp

    def kap(self, T, rho):
        return 10.0**self.sk.ev(np.log10(T), np.log10(rho))

    def thermo(self, T, rho):
        """gas EOS + separate radiation (total P = Pg + aT^4/3, e = eg + aT^4/rho).
        Returns P, Pg, chi_rho, chi_T (of the TOTAL), c_P, delta, nabla_ad (of the
        total)."""
        eos = self.eos
        lt, lr = np.log10(T), np.log10(rho)
        pg = rho*10.0**eos.lp(lt, lr)
        dlp_dlr, dlp_dlt = eos.lp(lt, lr, dy=1), eos.lp(lt, lr, dx=1)
        eg = 10.0**eos.le(lt, lr)
        deg_dT = eg*eos.le(lt, lr, dx=1)/T
        pr = A_RAD_*T**4/3.0
        P = pg + pr
        chir = pg*(1.0 + dlp_dlr)/P
        chit = (pg*dlp_dlt + 4.0*pr)/P
        cv = deg_dT + 4.0*A_RAD_*T**3/rho
        cp = cv + P*chit**2/(rho*T*chir)
        dlt = chit/chir
        nad = P*dlt/(rho*T*cp)
        return P, pg, chir, chit, cp, dlt, nad

    def mlt_gradient(self, T, rho, r, Ffrac=1.0, pp=True):
        """local MLT (+ MLT++). Returns nabla, F_rad/F, F_MLT/F, v_mlt (textbook,
        unboosted), kappa, P, Hp, nabla_ad, nabla_rad, Pg.  pp=False: textbook MLT.
        MLT_OPTION 'KW' (legacy AG Car builders): Kippenhahn-Weigert cubic in xi with U;
        'Cox' (MESA r15140 default, mlt_get_results.f90:1088-1385): Cox & Giuli 14.78-14.82
        with A = K_conv/K_rad, a0 = 9/4, Zeta = Gamma^3/B^3, nabla = (1 - Zeta) nabla_rad +
        Zeta nabla_ad.  ALT_HP (MESA alt_scale_height_flag, default true): H_p = min(P/(rho g),
        sqrt(P/G)/rho) for the mixing length (and the returned Hp)."""
        P, pg, chir, chit, cp, dlt, nad = self.thermo(T, rho)
        lt, lr = np.log10(T), np.log10(rho)
        kap = 10.0**self.sk.ev(lt, lr)
        g = self.GM/r**2
        F = Ffrac*self.LUM/(4*np.pi*r**2)
        Hp = P/(rho*g)
        if self.ALT_HP:
            Hp = min(Hp, np.sqrt(P/G_CODE)/rho)
        nrad = 3.0*kap*P*F/(4.0*A_RAD_*CL*g*T**4)
        W = nrad - nad
        if W <= 0.0:
            return nrad, 1.0, 0.0, 0.0, kap, P, Hp, nad, nrad, pg
        lm = self.ALPHA*Hp
        if self.MLT_OPTION == 'KW':
            U = 3.0*A_RAD_*CL*T**3/(cp*rho**2*kap*lm**2)*np.sqrt(8.0*Hp/(g*dlt))

            def f(x):
                return x**3 + 8.0*U/9.0*x**2 + 16.0*U**2/9.0*x - 8.0*U/9.0*W
            hi = max(1e-12, min(np.sqrt(W), (8.0*U/9.0*W)**(1.0/3.0)) * 1.0001 + 1e-30)
            while f(hi) < 0.0:
                hi *= 2.0
            xi = brentq(f, 0.0, hi, xtol=1e-14, rtol=1e-13)
            nab = nad + 2.0*U*xi + xi**2
            v = lm*np.sqrt(g*dlt/(8.0*Hp))*xi    # (the textbook bubble speed, unboosted)
        else:                                   # MESA 'Cox' (C&G 1968 ch. 14)
            x = dlt*rho/(2*P)                   # Q = chiT/chiRho = delta
            kc = cp*g*lm*lm*rho*np.sqrt(x)/9.0
            kr = (4*A_RAD_*CL/3)*T**3/(kap*rho)
            a0 = 9.0/4.0
            A = kc/kr
            B3 = (A*A/a0)*W
            delta = a0*B3
            f_ = -2 + 9*a0 + 27*a0*a0*delta
            f0 = f_ if f_ > 1e100 else f_*f_ + 4*(-1 + 3*a0)**3
            f0 = f_ if f0 <= 0 else (np.sqrt(f0) if f_ <= 1e100 else f0)
            f1 = -2 + 9*a0 + 27*a0*a0*delta + f0
            if f1 < 0:
                return nrad, 1.0, 0.0, 0.0, kap, P, Hp, nad, nrad, pg
            f1 = f1**(1.0/3.0)
            f2 = 2*2.0**(1.0/3.0)*(1 - 3*a0)/f1
            Gam = (4.0**(1.0/3.0)*f1 + f2 - 2)/(6*a0)
            zeta = min(max(Gam**3/B3, 0.0), 1.0)
            nab = (1.0 - zeta)*nrad + zeta*nad
            v = self.ALPHA*np.sqrt(dlt*P/(8*rho))*Gam/A if A > 0 else 0.0   # MESA conv_vel (:1317-1319)
        if pp and self.PP_ALPHA > 0.0 and nab - nad > self.PP_F1 and lt <= self.PP_MAXLOGT:
            if self.PP_FORM == 'mesa':
                af = self.PP_F2 + (1.0 - self.PP_F2)*(1.0 - self.PP_ALPHA)
                nab = af*nab + (1.0 - af)*nad
            else:
                nab = nad + self.PP_F2*(nab - nad)
        return nab, nab/nrad, 1.0 - nab/nrad, v, kap, P, Hp, nad, nrad, pg


PP_CTRL = dict(lambda1=1.0, beta1=0.35, lambda2=0.5, beta2=0.25, dlambda=0.1, dbeta=0.1,
               xa1_center=0.0)


def gradT_excess_alpha(beta, lam, c=PP_CTRL):
    """MESA r15140 star/private/mlt_info.f90:1427-1460 (+ the snapping at 1480-1481): the
    fraction of the full MLT++ boost in the (beta, lambda) plane; beta = min P_gas/P x
    (1 + xa(1, centre)), lambda = min(1, max L_rad/L_Edd from the surface to the beta-minimum
    cell).  No time smoothing (a static IC)."""
    l1, b1, l2, b2 = c['lambda1'], c['beta1'], c['lambda2'], c['beta2']
    dl, db = c['dlambda'], c['dbeta']
    if l1 < 0:
        al = 1.0
    elif lam >= l1:
        al = 1.0 if beta <= b1 else ((b1 + db - beta)/db if beta < b1 + db else 0.0)
    elif lam >= l2:
        bl = b2 + (lam - l2)*(b1 - b2)/(l1 - l2)
        al = 1.0 if beta <= bl else ((bl + db - beta)/db if beta < bl + db else 0.0)
    elif lam > l2 - dl:
        al = 1.0 if beta <= b2 else ((lam - (l2 - dl))/dl if beta < b2 + db else 0.0)
    else:
        al = 0.0
    if al < 1e-4:
        al = 0.0
    if al > 0.9999:
        al = 1.0
    return al
