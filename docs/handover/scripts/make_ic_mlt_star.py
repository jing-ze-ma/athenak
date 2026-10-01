#!/usr/bin/env python3
"""IC column for pgen he_star_m1, GENERAL star: the 'mlt' mode of make_ic_he_presn_m1.py
with every star- and microphysics-specific input a parameter.

Structure (same equations as make_ic_he_presn_m1.py mode mlt, which stays as it is):
  * photosphere at R_ph (tau = 2/3; R_ph given, or from L and Teff), T from the grey Hopf
    relation T^4 = (3/4) Teff^4 (tau + q(tau));
  * ABOVE R_ph, --atm:
      hopf  (the He star): Hopf T(tau) continued upward, gas hydrostatic
            dPg/dr = -rho (g - kappa_R F/c), F = L/(4 pi r^2); rho_R SHOT so that tau ->
            TAU_TOP (1e-4) where rho reaches RHO_TOP (1e-14); integrated 5e9 cm past it;
      iso   (Ma, Bildsten & Jiang 2026, sect. 2.2.2): isothermal at T_ph = T(tau = 2/3),
            hydrostatic with the local g (--atm_rad 1 [default]: g - kappa_R F/c, the
            balance the pgen's force_reference = wb_arad reference holds; 0: g alone, the
            paper's exp(-(r - r_ph) g rho_ph/P_ph)), down to the density floor
            --rho_floor,
            constant rho_floor above; rho_R SHOT so that the atmosphere's tau(R_ph) = 2/3;
  * BELOW R_ph: Hopf + gas hydrostatic balance with the radiation force of F down to
    tau = TAU_J, then (T, rho) inward with dPtot/dr = -rho g, dlnT/dlnPtot = nabla (local
    MLT, alpha --alpha, Bohm-Vitense bubble model), F_rad = F nabla/nabla_rad,
    F_MLT = F - F_rad, F = L/(4 pi r^2): the radiation force comes from F_rad only (the
    gas balance is dPg/dr = -rho g + rho kappa_R F_rad/c), down to --rin.
EOS: --eos table:<dump file> (the run's gas-only table dump, log10 e/rho and log10 p/rho
on (log rho, log T)) or --eos ideal:<mu>:<gamma>; radiation (a T^4/3, a T^4) is added to
both.
Opacity: --ross <table> (the he_star_m1 table format), bilinear in (log T, log rho),
optionally extended to low rho as the pgen's problem/he_opac_logd_min does (--ld_ext).
Output: <out>/<name>: r rho eint F_r(radiative part) E=aT^4 T F_MLT/F  (he_ic_cols = 5,
mlt_flux_frozen = true reads column 7), plus <out>/mlt_struct.npz.

Presets: --preset he reproduces make_ic_he_presn_m1.py mode mlt bit for bit (the 4.0 Msun
He star); --preset bsg is the 20 Msun BSG of Ma, Bildsten & Jiang 2026.
"""
import argparse
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import brentq

CL = 2.99792458e10
SIG = 5.670374419e-5
A_RAD_ = 7.5657332503e-15
KB = 1.3806488e-16          # = Units::k_boltzmann_cgs
MU_ = 1.660538921e-24       # = Units::atomic_mass_unit_cgs
RSUN = 6.957e10
LSUN = 3.828e33
GMSUN = 1.3271244e26

PRESETS = {
    'he': dict(gm=4.18143e26, lum=2.3066e38, rstar=2.3717e11, teff=None,
               eos='table:/viper/ptmp2/jinma/hepresn_0929/'
                   'eos_table_he_x0_y0.98_z0.02_logd-14.txt',
               ross='/viper/ptmp2/jinma/caltech_handover_0926/athenak_data/he_box/'
                    'rosseland_he_x0.0_z0.02.txt',
               atm='hopf', rin=0.5*2.3717e11 - 9.0e9, drn=2.0e6, tauj=3.0,
               name='ic_he_presn_m1_mlt.txt', ld_ext=None),
    'bsg': dict(gm=20.0*GMSUN, lum=LSUN*10.0**5.161, rstar=None, teff=15835.0,
                eos='ideal:0.62:1.6666666666666667',
                ross='/viper/ptmp2/jinma/bsg_1001/opac/rosseland_tops_x0.7_z0.008.txt',
                atm='iso', rin=17.0*RSUN, drn=2.0e8, tauj=3.0,
                name='ic_bsg_mlt.txt', ld_ext=-18.0),
}


class TableEOS:
    """the run's gas-only table dump (make_ic_he_presn_m1.read_dump)"""
    def __init__(self, fn):
        with open(fn) as fh:
            lines = [fh.readline() for _ in range(6)]
        nx, ny, xmin, dx, ymin, dy = [float(v) for v in lines[2].split()[1:]]
        nx, ny = int(nx), int(ny)
        d = np.loadtxt(fn, comments='#')
        x = xmin + dx*np.arange(nx)
        y = ymin + dy*np.arange(ny)
        self.sle = RectBivariateSpline(y, x, d[:, 0].reshape(ny, nx), kx=3, ky=3)
        self.slp = RectBivariateSpline(y, x, d[:, 1].reshape(ny, nx), kx=3, ky=3)
        self.desc = 'table ' + fn

    def lp(self, lt, lr, dx=0, dy=0):       # log10 p/rho; dx: d/dlog T, dy: d/dlog rho
        return self.slp.ev(lt, lr, dx=dx, dy=dy)

    def le(self, lt, lr, dx=0, dy=0):       # log10 e/rho (specific)
        return self.sle.ev(lt, lr, dx=dx, dy=dy)


class IdealEOS:
    """ideal gas, p = rho k T/(mu m_u), e = p/(gamma - 1)"""
    def __init__(self, mu, gam):
        self.c = np.log10(KB/(mu*MU_))
        self.cg = np.log10(gam - 1.0)
        self.desc = 'ideal mu %.6g gamma %.10g' % (mu, gam)

    def lp(self, lt, lr, dx=0, dy=0):
        if dy:
            return 0.0*lt
        if dx:
            return 1.0 + 0.0*lt
        return self.c + lt

    def le(self, lt, lr, dx=0, dy=0):
        if dy:
            return 0.0*lt
        if dx:
            return 1.0 + 0.0*lt
        return self.c - self.cg + lt


def read_ross(fn, ld_ext=None):
    """he_star_m1 table format; ld_ext: HsReadOpacityTable's extension to low rho
    (log10 kappa continued in log rho with the first interval's slope clamped to
    [0, 1])"""
    with open(fn) as fh:
        for ln in fh:
            if ln.startswith('# ') and len(ln.split()) == 7 and ln.split()[1].isdigit():
                nT, nD, lT0, dlT, lD0, dlD = [float(v) for v in ln.split()[1:]]
    nT, nD = int(nT), int(nD)
    v = np.loadtxt(fn, comments='#').reshape(nT, nD)
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
    def __init__(self, a):
        self.GM, self.LUM = a.gm, a.lum
        if a.rstar is not None:
            self.RSTAR = a.rstar
        else:
            self.RSTAR = np.sqrt(self.LUM/(4*np.pi*SIG*a.teff**4))
        if a.eos.startswith('table:'):
            self.eos = TableEOS(a.eos[6:])
        else:
            _, mu, gam = a.eos.split(':')
            self.eos = IdealEOS(float(mu), float(gam))
        ty, tx, kv = read_ross(a.ross, a.ld_ext)
        self.sk = RectBivariateSpline(ty, tx, kv, kx=1, ky=1)
        self.ALPHA = a.alpha

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

    def mlt_gradient(self, T, rho, r, Ffrac=1.0):
        """local MLT (make_ic_he_presn_m1.mlt_gradient). Returns nabla, F_rad/F, F_MLT/F,
        v_mlt, kappa, P, Hp, nabla_ad, nabla_rad, Pg."""
        P, pg, chir, chit, cp, dlt, nad = self.thermo(T, rho)
        lt, lr = np.log10(T), np.log10(rho)
        kap = 10.0**self.sk.ev(lt, lr)
        g = self.GM/r**2
        F = Ffrac*self.LUM/(4*np.pi*r**2)
        Hp = P/(rho*g)
        nrad = 3.0*kap*P*F/(4.0*A_RAD_*CL*g*T**4)
        W = nrad - nad
        if W <= 0.0:
            return nrad, 1.0, 0.0, 0.0, kap, P, Hp, nad, nrad, pg
        lm = self.ALPHA*Hp
        U = 3.0*A_RAD_*CL*T**3/(cp*rho**2*kap*lm**2)*np.sqrt(8.0*Hp/(g*dlt))

        def f(x):
            return x**3 + 8.0*U/9.0*x**2 + 16.0*U**2/9.0*x - 8.0*U/9.0*W
        hi = max(1e-12, min(np.sqrt(W), (8.0*U/9.0*W)**(1.0/3.0)) * 1.0001 + 1e-30)
        while f(hi) < 0.0:
            hi *= 2.0
        xi = brentq(f, 0.0, hi, xtol=1e-14, rtol=1e-13)
        nab = nad + 2.0*U*xi + xi**2
        v = lm*np.sqrt(g*dlt/(8.0*Hp))*xi
        return nab, nab/nrad, 1.0 - nab/nrad, v, kap, P, Hp, nad, nrad, pg


def build(a):
    st = Star(a)
    GM, LUM, RSTAR, eos, sk = st.GM, st.LUM, st.RSTAR, st.eos, st.sk
    Teff = (LUM/(4*np.pi*RSTAR**2*SIG))**0.25
    TAU_J = a.tauj
    TAU_TOP = 1.0e-4
    RHO_TOP = 1.0e-14
    print('  star: GM %.6e L %.6e R_ph %.6e (%.4f Rsun) Teff %.1f; EOS %s; kappa_R %s'
          % (GM, LUM, RSTAR, RSTAR/RSUN, Teff, eos.desc, a.ross))

    def kap_of(T_, rho_):
        return 10.0**sk.ev(np.log10(T_), np.log10(rho_))

    def rhs1(rr, y):                    # Hopf atmosphere, y = (ln rho, tau)
        lnr, tau = y
        tau = max(tau, 1.0e-8)
        rho_ = np.exp(lnr)
        F = LUM/(4*np.pi*rr**2)
        q = q_hopf(tau)
        T_ = (3.0*F/(4.0*SIG)*(tau + q))**0.25
        kap_ = kap_of(T_, rho_)
        dtau = -kap_*rho_
        dlnT = 0.25*(-2.0/rr + (1.0 + 0.1331*3.4488*np.exp(-3.4488*tau))*dtau/(tau + q))
        lt_, lr_ = np.log10(T_), np.log10(rho_)
        p_ = rho_*10.0**eos.lp(lt_, lr_)
        geff = GM/rr**2 - kap_*F/CL
        dlnp = -rho_*geff/p_
        dlnr = (dlnp - eos.lp(lt_, lr_, dx=1)*dlnT)/(1.0 + eos.lp(lt_, lr_, dy=1))
        return [dlnr, dtau]

    Tph = (3.0*(LUM/(4*np.pi*RSTAR**2))/(4.0*SIG)*(2.0/3.0 + q_hopf(2.0/3.0)))**0.25
    rho_fl = a.rho_floor

    def rhs_iso(rr, y):                 # isothermal atmosphere at T_ph, y = (ln rho, tau)
        rho_ = np.exp(y[0])
        F = LUM/(4*np.pi*rr**2)
        kap_ = kap_of(Tph, rho_)
        lt_, lr_ = np.log10(Tph), np.log10(rho_)
        p_ = rho_*10.0**eos.lp(lt_, lr_)
        geff = GM/rr**2 - (kap_*F/CL if a.atm_rad else 0.0)
        return [-rho_*geff/p_/(1.0 + eos.lp(lt_, lr_, dy=1)), -kap_*rho_]

    def rhs_re(rr, y):                  # grey RE atmosphere, y = (ln rho, tau)
        lnr, tau = y
        rho_ = np.exp(lnr)
        T_, dlnT = t_re(rr, tau, rho_)
        F = LUM/(4*np.pi*rr**2)
        kap_ = kap_of(T_, rho_)
        lt_, lr_ = np.log10(T_), np.log10(rho_)
        p_ = rho_*10.0**eos.lp(lt_, lr_)
        geff = GM/rr**2 - kap_*F/CL
        dlnp = -rho_*geff/p_
        dlnr = (dlnp - eos.lp(lt_, lr_, dx=1)*dlnT)/(1.0 + eos.lp(lt_, lr_, dy=1))
        return [dlnr, -kap_*rho_]

    def t_re(rr, tau, rho_):
        # T^4 = (3/4) T_eff^4 (tau + q(tau)) 2 W(r): Hopf at r = R_ph, diluted by the
        # geometric factor W = (1 - sqrt(1 - (R/r)^2))/2 above (W = 1/2 below R)
        tq = max(tau, 0.0)
        x = min(RSTAR/rr, 1.0)
        w2 = 1.0 - np.sqrt(max(1.0 - x*x, 0.0))
        T_ = (0.75*Teff**4*(tq + q_hopf(tq))*w2)**0.25
        dq = 0.1331*3.4488*np.exp(-3.4488*tq)
        dtau = -kap_of(T_, rho_)*rho_ if tau > 0.0 else 0.0
        dw2 = (-x*x/rr/np.sqrt(max(1.0 - x*x, 1e-30))) if x < 1.0 else 0.0
        dlnT = 0.25*((1.0 + dq)*dtau/(tq + q_hopf(tq)) + dw2/w2)
        return T_, dlnT

    if a.atm == 're':
        def shoot(lnrho_R):
            s = solve_ivp(rhs_re, [RSTAR, a.rtop], [lnrho_R, 2.0/3.0],
                          method='LSODA', rtol=1e-10, atol=1e-13)
            return s

        def resid(lnrho_R):             # tau left at rtop: 0 wanted
            s = shoot(lnrho_R)
            return s.y[1][-1] if s.success else -1.0
    elif a.atm == 'hopf':
        def top_event(rr, y):
            return y[0] - np.log(RHO_TOP)
        top_event.terminal = True

        def shoot(lnrho_R):
            s = solve_ivp(rhs1, [RSTAR, RSTAR + 3.0e11], [lnrho_R, 2.0/3.0],
                          events=top_event, method='LSODA', rtol=1e-10, atol=1e-13)
            return s

        def resid(lnrho_R):
            s = shoot(lnrho_R)
            if s.t_events[0].size == 0:
                return -1.0             # never reaches the floor (Gamma > 1): too dense
            return s.y_events[0][0][1] - TAU_TOP
    else:
        def top_event(rr, y):
            return y[0] - np.log(rho_fl)
        top_event.terminal = True

        def blow_event(rr, y):          # Gamma > 1: rho grows outward, give up
            return y[0] - np.log(1.0e-6)
        blow_event.terminal = True

        def shoot(lnrho_R):
            s = solve_ivp(rhs_iso, [RSTAR, RSTAR + 2.0*RSTAR], [lnrho_R, 2.0/3.0],
                          events=[top_event, blow_event], method='LSODA', rtol=1e-10,
                          atol=1e-13)
            return s

        def resid(lnrho_R):             # tau left at the floor: 0 wanted
            s = shoot(lnrho_R)
            if s.t_events[0].size == 0:
                return -1.0
            return s.y_events[0][0][1]
    a_, b_ = np.log(a.rho_lo), np.log(a.rho_hi)
    fa, fb = resid(a_), resid(b_)
    assert fa*fb < 0, (fa, fb)
    lnrho_R = brentq(resid, a_, b_, xtol=1e-12)
    rho_R = np.exp(lnrho_R)
    print('  mlt: Teff %.1f K, photosphere rho_R = %.5e (atm %s)' % (Teff, rho_R, a.atm))

    DRN = a.drn
    su = shoot(lnrho_R)
    if a.atm == 're':
        rx = np.arange(RSTAR, a.rtop, DRN)
        su = solve_ivp(rhs_re, [RSTAR, rx[-1]], [lnrho_R, 2.0/3.0], t_eval=rx,
                       method='LSODA', rtol=1e-10, atol=1e-13)
        assert su.success, su.message
        r_u = su.t[1:]
        rho_u = np.exp(su.y[0][1:])
        T_u = np.array([t_re(rr, tt, dd)[0]
                        for rr, tt, dd in zip(r_u, su.y[1][1:], rho_u)])
        print('  RE atmosphere to %.4e (%.3f Rsun): rho %.3e, T %.1f K, tau %.2e there'
              % (r_u[-1], r_u[-1]/RSUN, rho_u[-1], T_u[-1], su.y[1][-1]))
    elif a.atm == 'hopf':
        rtop_atm = su.t_events[0][0]
        rx = np.arange(RSTAR, rtop_atm + 5.0e9, DRN)
        su = solve_ivp(rhs1, [RSTAR, rx[-1]], [lnrho_R, 2.0/3.0], t_eval=rx,
                       method='LSODA', rtol=1e-10, atol=1e-13)
        assert su.success, su.message
        tau_u = np.maximum(su.y[1][1:], 1e-8)
        r_u = su.t[1:]
        T_u = (3.0*(LUM/(4*np.pi*r_u**2))/(4.0*SIG)*(tau_u + q_hopf(tau_u)))**0.25
        rho_u = np.exp(su.y[0][1:])
    else:
        rtop_atm = su.t_events[0][0]
        rx = np.arange(RSTAR, rtop_atm, DRN)
        su = solve_ivp(rhs_iso, [RSTAR, rtop_atm], [lnrho_R, 2.0/3.0], t_eval=rx,
                       method='LSODA', rtol=1e-10, atol=1e-13)
        assert su.success, su.message
        r_f = np.arange(rx[-1] + DRN, a.rtop, DRN*a.fl_step)
        r_u = np.concatenate([su.t[1:], r_f])
        rho_u = np.concatenate([np.exp(su.y[0][1:]), np.full(r_f.size, rho_fl)])
        T_u = np.full(r_u.size, Tph)
        Hph = 10.0**eos.lp(np.log10(Tph), np.log10(rho_R))/(GM/RSTAR**2)
        print('  iso atmosphere: T_ph %.1f K, H(R_ph; g only) %.4e cm (%.4f Rsun), floor '
              '%.1e at r = %.6e (%.4f Rsun), constant above to %.4e'
              % (Tph, Hph, Hph/RSUN, rho_fl, rtop_atm, rtop_atm/RSUN, a.rtop))

    def evj(rr, y):
        return y[1] - TAU_J
    evj.terminal = True
    rd = np.arange(RSTAR, 0.9*RSTAR, -DRN)
    s1 = solve_ivp(rhs1, [RSTAR, 0.9*RSTAR], [lnrho_R, 2.0/3.0], t_eval=rd, events=evj,
                   method='LSODA', rtol=1e-10, atol=1e-13)
    rj = s1.t_events[0][0]
    yj = s1.y_events[0][0]
    Tj = (3.0*(LUM/(4*np.pi*rj**2))/(4.0*SIG)*(yj[1] + q_hopf(yj[1])))**0.25
    rhoj = np.exp(yj[0])
    print('  mlt: Hopf -> MLT join at tau = %.2f, r = %.6e (r/R %.4f), T %.4e rho %.4e'
          % (TAU_J, rj, rj/RSTAR, Tj, rhoj))

    def rhs2(rr, y):                    # y = (ln T, ln rho, tau); MLT interior
        T_, rho_ = np.exp(y[0]), np.exp(y[1])
        nab, ffr, fmlt, v, kap_, P, Hp, nad, nrad, pg = st.mlt_gradient(T_, rho_, rr)
        g = GM/rr**2
        dlnP = -rho_*g/P
        dlnT = nab*dlnP
        P_, pg_, chir, chit, cp, dlt, nad_ = st.thermo(T_, rho_)
        dlnr = (dlnP - chit*dlnT)/chir
        return [dlnT, dlnr, -kap_*rho_]

    rin = a.rin
    rn = np.arange(rj, rin, -DRN)
    s2 = solve_ivp(rhs2, [rj, rin], [np.log(Tj), np.log(rhoj), yj[1]], t_eval=rn,
                   method='LSODA', rtol=1e-9, atol=1e-12)
    assert s2.success, s2.message
    m1 = s1.t > rj
    r_h = s1.t[m1]
    tau_h = np.maximum(s1.y[1][m1], 1e-8)
    T_h = (3.0*(LUM/(4*np.pi*r_h**2))/(4.0*SIG)*(tau_h + q_hopf(tau_h)))**0.25
    rho_h = np.exp(s1.y[0][m1])
    r = np.concatenate([s2.t[::-1], r_h[::-1], r_u])
    T = np.concatenate([np.exp(s2.y[0])[::-1], T_h[::-1], T_u])
    rho = np.concatenate([np.exp(s2.y[1])[::-1], rho_h[::-1], rho_u])
    nn = len(r)
    ffr = np.ones(nn)
    fmlt = np.zeros(nn)
    vmlt = np.zeros(nn)
    nab_ = np.zeros(nn)
    nrad_ = np.zeros(nn)
    nad_ = np.zeros(nn)
    Hp_ = np.zeros(nn)
    isc = r < rj
    for i in np.where(isc)[0]:
        o = st.mlt_gradient(T[i], rho[i], r[i])
        nab_[i], ffr[i], fmlt[i], vmlt[i], _, _, Hp_[i], nad_[i], nrad_[i], _ = o
    Fr = ffr*LUM/(4*np.pi*r**2)
    lt, lr = np.log10(T), np.log10(rho)
    eint = rho*10.0**eos.le(lt, lr)
    pg = rho*10.0**eos.lp(lt, lr)
    kap = 10.0**sk.ev(lt, lr)
    g = GM/r**2
    tau = np.zeros(nn)
    tau[-1] = 1.0e-4 if a.atm == 'hopf' else 0.0
    for i in range(nn - 2, -1, -1):
        tau[i] = tau[i+1] + 0.5*(r[i+1] - r[i])*(rho[i]*kap[i] + rho[i+1]*kap[i+1])
    gam = kap*Fr/CL/g
    os.makedirs(a.out, exist_ok=True)
    fn = os.path.join(a.out, a.name)
    E = A_RAD_*T**4
    with open(fn, 'w') as fh:
        fh.write('# he_star_m1 IC (make_ic_mlt_star.py, he_ic_cols = 5): atm %s, '
                 'inward MLT '
                 '(alpha %.2f); GM = %.6e L = %.6e R_ph = %.6e; EOS %s\n'
                 '# r[cm] rho eint F_r(radiative part) E=aT^4  | T[K] F_MLT/F '
                 '(informational, cols 6,7)\n'
                 % (a.atm, st.ALPHA, GM, LUM, RSTAR, eos.desc))
        for row in zip(r, rho, eint, Fr, E, T, fmlt):
            fh.write('%.10e %.10e %.10e %.10e %.10e %.10e %.6e\n' % row)
    np.savez(os.path.join(a.out, 'mlt_struct.npz'), r=r, rho=rho, T=T, Fr=Fr, fmlt=fmlt,
             tau=tau, nab=nab_, nrad=nrad_, nad=nad_, Hp=Hp_, kap=kap, vmlt=vmlt, pg=pg,
             gam=gam, rstar=RSTAR, gm=GM, lum=LUM, rj=rj)
    print('wrote', fn)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--preset', choices=sorted(PRESETS), default='he')
    p.add_argument('--out', default='.')
    for k in ('gm', 'lum', 'rstar', 'teff', 'rin', 'drn', 'tauj', 'ld_ext'):
        p.add_argument('--' + k, type=float)
    p.add_argument('--eos')
    p.add_argument('--ross')
    p.add_argument('--atm', choices=['hopf', 'iso', 're'])
    p.add_argument('--name')
    p.add_argument('--alpha', type=float, default=1.5)
    p.add_argument('--atm_rad', type=int, default=1)
    p.add_argument('--rho_floor', type=float, default=1.0e-16)
    p.add_argument('--rtop', type=float, default=85.0*RSUN)
    p.add_argument('--fl_step', type=float, default=100.0)
    p.add_argument('--rho_lo', type=float, default=1.0e-11)
    p.add_argument('--rho_hi', type=float, default=3.0e-9)
    a = p.parse_args()
    for k, v in PRESETS[a.preset].items():
        if getattr(a, k, None) is None:
            setattr(a, k, v)
    build(a)


if __name__ == '__main__':
    main()
