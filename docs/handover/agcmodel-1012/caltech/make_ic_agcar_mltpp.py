#!/usr/bin/env python3
"""AG Car (general EOS) MLT++ IC column for pgen he_star_m1: a Raven reconstruction of
lbv_1008/agcar/geos/scripts/make_ic_agcar_ge.py (Viper only), built on make_ic_mlt_star.py
(same directory: its table EOS, Rosseland reader, Hopf q(tau) and local MLT are used).

icB2-1011 (worker K).  Two modes:

  --mode inward  (reproduce ic_agcar_B_ge.txt): R_ph from L, Teff; grey Hopf T(tau) with
      the gas HSE dPg/dr = -rho (g - kappa_R F/c), FULL F, from R_ph (tau 2/3) inward to
      tau = --tauj (3); below it (T, rho) inward with dPtot/dr = -rho g, dlnT/dlnPtot =
      nabla_MLT++ on a log grid (dln r = --dlnr from r_join down to --rin);
      above R_ph an isothermal (T_ph = T_Hopf(2/3)) gravity-only atmosphere (rho_R shot so
      that its tau(R_ph) = 2/3) down to --rho_floor, constant floor above (step
      --fl_step DRN) to --rtop; DRN = --dlnr R_ph.
  --mode outward (the new column): the rows of --base (the old file) with r < r_join are
      copied VERBATIM (identical text), and from (r_join, T_j, rho_j, tau = tauj) the same
      MLT++ equations are integrated OUTWARD (log grid continued) to the new R_ph'; above
      it an isothermal atmosphere at T_ph' = T(R_ph'), starting at rho(R_ph'), floor and
      E as the old file.  MLT++ (and F_MLT) thus extends to tau_R = 2/3 instead of 3.
      --stop internal (default): R_ph' where the old column's tau scale (tau(r_join) =
          tauj, the Hopf layer's) reaches 2/3; the cap is gravity-only (--cap_rad 0, as the
          old file) and has its own tau 0.19, so the file's true tau_R at R_ph' is 0.19.
      --stop atm: R_ph' where the cap on top has tau 2/3 (true tau_R 2/3, as the old
          shoot of rho_R); needs --cap_rad 1 (cap with g - kappa F/c: a gravity-only cap
          has tau < 0.48 even on r_join).

MLT++ as recovered from the old column (Paxton+2013 sect. 7.2, f1 1e-4, f2 1e-3, alpha 1):
the local MLT gradient (alpha_MLT --alpha, Bohm-Vitense) is replaced by
    nabla = nabla_ad + f2 (nabla_MLT - nabla_ad)   where nabla_MLT - nabla_ad > f1,
and left alone otherwise.  (Fit to ic_agcar_B_ge.txt: (nabla - nabla_ad)/(nabla_MLT -
nabla_ad) = 1.000e-3 +- 1e-5 in every convective row with excess > 1e-4, = 1 below; the
switch sits between excesses 9.68e-5 and 1.004e-4.)
F_r = F nabla/nabla_rad, F_MLT/F = 1 - nabla/nabla_rad, F = L/(4 pi r^2).

E column: a T^4 at and below R_ph; above it E = (F/c) (1 + (c a T_ph^4/F_ph - 1)(R_ph/r)^2)
(recovered from the old file: E/(F/c) - 1 ~ r^-2 above R_ph).  The run replaces E there by
the sph1d solution (he_ic_rad = sph1d), so this column only matters with he_ic_rad = file.

Output columns as ic_agcar_B_ge.txt: r rho eint F_r(radiative) E T F_MLT/F.
"""
import argparse
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_ic_mlt_star as M   # noqa: E402

RSUN = M.RSUN


class StarPP(M.Star):
    """make_ic_mlt_star.Star with the MLT++ excess reduction"""
    def __init__(self, a):
        super().__init__(a)
        self.f1, self.f2, self.pp = a.f1, a.f2, a.mltpp

    def mlt_gradient(self, T, rho, r, Ffrac=1.0):
        o = list(super().mlt_gradient(T, rho, r, Ffrac))
        nab, nad, nrad = o[0], o[7], o[8]
        if self.pp and nrad > nad and nab - nad > self.f1:
            nab = nad + self.f2*(nab - nad)
            o[0], o[1], o[2] = nab, nab/nrad, 1.0 - nab/nrad
        return tuple(o)


def build(a):
    st = StarPP(a)
    GM, LUM, eos, sk = st.GM, st.LUM, st.eos, st.sk
    CL, SIG, AR = M.CL, M.SIG, M.A_RAD_
    RSTAR = st.RSTAR
    Teff = (LUM/(4*np.pi*RSTAR**2*SIG))**0.25
    DRN = a.dlnr*RSTAR
    TAU_J = a.tauj
    q_hopf = M.q_hopf

    def kap_of(T_, rho_):
        return 10.0**sk.ev(np.log10(T_), np.log10(rho_))

    def rhs1(rr, y):                    # Hopf layer, y = (ln rho, tau); full F
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

    def rhs2(rr, y):                    # MLT++ interior, y = (ln T, ln rho, tau)
        T_, rho_ = np.exp(y[0]), np.exp(y[1])
        nab, ffr, fmlt, v, kap_, P, Hp, nad, nrad, pg = st.mlt_gradient(T_, rho_, rr)
        g = GM/rr**2
        dlnP = -rho_*g/P
        dlnT = nab*dlnP
        P_, pg_, chir, chit, cp, dlt, nad_ = st.thermo(T_, rho_)
        dlnr = (dlnP - chit*dlnT)/chir
        return [dlnT, dlnr, -kap_*rho_]

    def iso_atm(Tph, r0, rho0, rad=0):
        """isothermal gravity-only (--atm_rad 0) or g - kappa F/c atmosphere at Tph from
        r0; returns the shot solution object pieces"""
        rho_fl = a.rho_floor

        def rhs_iso(rr, y):
            rho_ = np.exp(y[0])
            F = LUM/(4*np.pi*rr**2)
            kap_ = kap_of(Tph, rho_)
            lt_, lr_ = np.log10(Tph), np.log10(rho_)
            p_ = rho_*10.0**eos.lp(lt_, lr_)
            geff = GM/rr**2 - (kap_*F/CL if rad else 0.0)
            return [-rho_*geff/p_/(1.0 + eos.lp(lt_, lr_, dy=1)), -kap_*rho_]

        def top_event(rr, y):
            return y[0] - np.log(rho_fl)
        top_event.terminal = True

        def blow_event(rr, y):
            return y[0] - np.log(1.0e-6)
        blow_event.terminal = True

        def shoot(lnrho_R):
            return solve_ivp(rhs_iso, [r0, r0 + 2.0*r0], [lnrho_R, 2.0/3.0],
                             events=[top_event, blow_event], method='LSODA', rtol=1e-10,
                             atol=1e-13)
        if rho0 is None:                # shoot rho_R: the atmosphere's tau(R_ph) = 2/3
            def resid(lnrho_R):
                s = shoot(lnrho_R)
                if s.t_events[0].size == 0:
                    return -1.0
                return s.y_events[0][0][1]
            lo, hi = np.log(a.rho_lo), np.log(a.rho_hi)
            assert resid(lo)*resid(hi) < 0
            lnrho_R = brentq(resid, lo, hi, xtol=1e-12)
        else:
            lnrho_R = np.log(rho0)
        su = shoot(lnrho_R)
        rtop_atm = su.t_events[0][0]
        rx = np.arange(r0, rtop_atm, DRN)
        su = solve_ivp(rhs_iso, [r0, rtop_atm], [lnrho_R, 2.0/3.0], t_eval=rx,
                       method='LSODA', rtol=1e-10, atol=1e-13)
        assert su.success, su.message
        r_f = np.arange(rx[-1] + DRN, a.rtop, DRN*a.fl_step)
        r_u = np.concatenate([su.t[1:], r_f])
        rho_u = np.concatenate([np.exp(su.y[0][1:]), np.full(r_f.size, rho_fl)])
        tau_atm = 2.0/3.0 - su.y[1][-1]     # the atmosphere's own optical depth
        return np.exp(lnrho_R), r_u, rho_u, np.full(r_u.size, Tph), rtop_atm, tau_atm

    def iso_tau(Tph, r0, rho0, rad=0):
        """optical depth of the isothermal atmosphere (as iso_atm) on top of (r0, rho0)"""
        def rhs_iso(rr, y):
            rho_ = np.exp(y[0])
            F = LUM/(4*np.pi*rr**2)
            kap_ = kap_of(Tph, rho_)
            lt_, lr_ = np.log10(Tph), np.log10(rho_)
            p_ = rho_*10.0**eos.lp(lt_, lr_)
            geff = GM/rr**2 - (kap_*F/CL if rad else 0.0)
            return [-rho_*geff/p_/(1.0 + eos.lp(lt_, lr_, dy=1)), -kap_*rho_]

        def top_event(rr, y):
            return y[0] - np.log(a.rho_floor)
        top_event.terminal = True
        s = solve_ivp(rhs_iso, [r0, 3.0*r0], [np.log(rho0), 0.0], events=top_event,
                      method='LSODA', rtol=1e-10, atol=1e-13)
        assert s.t_events[0].size == 1
        return -s.y_events[0][0][1]

    # ---- the Hopf layer from R_ph inward to tau = TAU_J (needs rho_R of the atmosphere)
    Tph = (3.0*(LUM/(4*np.pi*RSTAR**2))/(4.0*SIG)*(2.0/3.0 + q_hopf(2.0/3.0)))**0.25
    rho_R, r_u, rho_u, T_u, rtop_atm, tau_atm = iso_atm(Tph, RSTAR, None, a.atm_rad)
    print('  star: GM %.8e L %.8e R_ph %.8e (%.4f Rsun) Teff %.2f T_ph %.2f rho_R %.6e'
          % (GM, LUM, RSTAR, RSTAR/RSUN, Teff, Tph, rho_R))

    def evj(rr, y):
        return y[1] - TAU_J
    evj.terminal = True
    rd = np.arange(RSTAR, 0.9*RSTAR, -DRN)
    s1 = solve_ivp(rhs1, [RSTAR, 0.9*RSTAR], [np.log(rho_R), 2.0/3.0], t_eval=rd,
                   events=evj, method='LSODA', rtol=1e-10, atol=1e-13, dense_output=True)
    rj = s1.t_events[0][0]
    yj = s1.y_events[0][0]
    Tj = (3.0*(LUM/(4*np.pi*rj**2))/(4.0*SIG)*(yj[1] + q_hopf(yj[1])))**0.25
    rhoj = np.exp(yj[0])
    print('  join at tau %.3f: r_j %.10e (r/R %.5f) T %.6e rho %.6e'
          % (TAU_J, rj, rj/RSTAR, Tj, rhoj))

    def mlt_cols(r, T, rho, isc):
        nn = len(r)
        ffr, fmlt = np.ones(nn), np.zeros(nn)
        for i in np.where(isc)[0]:
            o = st.mlt_gradient(T[i], rho[i], r[i])
            ffr[i], fmlt[i] = o[1], o[2]
        return ffr, fmlt

    def ecol(r, T, rph, tph):
        F = LUM/(4*np.pi*r**2)
        E = AR*T**4
        up = r > rph
        Fph = LUM/(4*np.pi*rph**2)
        E[up] = F[up]/CL*(1.0 + (CL*AR*tph**4/Fph - 1.0)*(rph/r[up])**2)
        return E

    hdr0 = ''
    if a.mode == 'inward':
        if a.grid_from:                 # the radii of an existing column below r_j
            rb = np.loadtxt(a.grid_from, usecols=0)
            rn = np.concatenate([[rj], rb[rb < rj*(1.0 - 1e-12)][::-1]])
        else:                           # log grid anchored at r_j
            nk = int(np.floor(np.log(rj/a.rin)/a.dlnr + 1e-9))
            rn = rj*np.exp(-a.dlnr*np.arange(nk + 1))
            rn = rn[rn >= a.rin]
        s2 = solve_ivp(rhs2, [rj, rn[-1]], [np.log(Tj), np.log(rhoj), yj[1]], t_eval=rn,
                       method='LSODA', rtol=a.rtol2, atol=a.atol2)
        assert s2.success, s2.message
        # Hopf rows on r_j + k DRN (r_j < r < R_ph) and the R_ph row itself (the old file)
        r_h = rj + DRN*np.arange(1, int(np.ceil((RSTAR - rj)/DRN)) + 1)
        r_h = np.concatenate([r_h[r_h < RSTAR*(1.0 - 1e-12)], [RSTAR]])
        yh = s1.sol(r_h)
        yh[:, -1] = [np.log(rho_R), 2.0/3.0]
        tau_h = np.maximum(yh[1], 1e-8)
        T_h = (3.0*(LUM/(4*np.pi*r_h**2))/(4.0*SIG)*(tau_h + q_hopf(tau_h)))**0.25
        rho_h = np.exp(yh[0])
        r = np.concatenate([s2.t[::-1], r_h, r_u])
        T = np.concatenate([np.exp(s2.y[0])[::-1], T_h, T_u])
        rho = np.concatenate([np.exp(s2.y[1])[::-1], rho_h, rho_u])
        ffr, fmlt = mlt_cols(r, T, rho, r < rj)
        rph, tph = RSTAR, Tph
        rows_keep = []
    else:
        # the old file below r_join, verbatim
        rows_keep, rk = [], []
        with open(a.base) as fh:
            for ln in fh:
                if ln.startswith('#'):
                    continue
                v = float(ln.split()[0])
                if v < rj*(1.0 - 1e-12):
                    rows_keep.append(ln)
                    rk.append([float(x) for x in ln.split()])
        rk = np.array(rk)
        print('  kept %d rows of %s below r_j (last r %.10e, tau-join row dropped)'
              % (len(rows_keep), a.base, rk[-1, 0]))

        tau_end = 2.0/3.0 if a.stop == 'internal' else -a.tauj   # (atm: search beyond)

        def evph(rr, y):
            return y[2] - tau_end
        evph.terminal = True
        nk = int(np.ceil(np.log(1.2*RSTAR/rj)/a.dlnr))
        ro = rj*np.exp(a.dlnr*np.arange(nk + 1))
        so = solve_ivp(rhs2, [rj, ro[-1]], [np.log(Tj), np.log(rhoj), yj[1]], t_eval=ro,
                       events=evph, method='LSODA', rtol=a.rtol2, atol=a.atol2,
                       dense_output=True)
        assert so.status == 1 or a.stop == 'atm', 'tau 2/3 not reached'
        r_ev = so.t_events[0][0] if so.status == 1 else so.t[-1]
        if a.stop == 'internal':        # the old column's tau scale: tau(r_j) = tauj
            rph = r_ev
        else:                           # atm: TRUE tau_R = 2/3 at R_ph', i.e. the iso
            def hres(rr):               # atmosphere on top has tau 2/3 (as the old shoot)
                y = so.sol(rr)
                return iso_tau(np.exp(y[0]), rr, np.exp(y[1]), a.cap_rad) - 2.0/3.0
            r_ev = so.t[-1] if so.status == 0 else r_ev
            assert hres(rj) > 0.0, hres(rj)
            k = 1
            while so.t[k] < r_ev and hres(so.t[k]) > 0.0:
                k += max(1, int(0.002/a.dlnr))
            k = min(k, len(so.t) - 1)
            assert hres(so.t[k]) < 0.0, 'atm stop not bracketed'
            rph = brentq(hres, so.t[max(k - int(0.002/a.dlnr), 0)], so.t[k], xtol=1e-6*DRN)
        yph = so.sol(rph)
        tph, rhoph = np.exp(yph[0]), np.exp(yph[1])
        mo = so.t < rph*(1.0 - 0.05*a.dlnr)     # (no sliver cell at R_ph')
        ro_ = np.concatenate([so.t[mo], [rph]])
        T_o = np.concatenate([np.exp(so.y[0][mo]), [tph]])
        rho_o = np.concatenate([np.exp(so.y[1][mo]), [rhoph]])
        tau_int = so.y[2][mo]           # the old column's tau scale (tau(r_j) = tauj)
        _, r_u, rho_u, T_u, rtop_atm, tau_atm = iso_atm(tph, rph, rhoph, a.cap_rad)
        tau_rj = tau_atm + (yj[1] - tau_int[-1]) + (tau_int[-1] - float(so.sol(rph)[2]))
        print('  outward MLT++ (stop %s): R_ph\' %.10e (%.5f R_ph0, %.4f Rsun) T_ph\' %.2f '
              'rho %.6e; Teff\' %.2f; iso atmosphere tau above R_ph\' %.4f, floor at %.6e;'
              ' true tau at r_j %.4f' % (a.stop, rph, rph/RSTAR, rph/RSUN, tph, rhoph,
                                         (LUM/(4*np.pi*rph**2*SIG))**0.25, tau_atm,
                                         rtop_atm, tau_rj))
        r = np.concatenate([ro_, r_u])
        T = np.concatenate([T_o, T_u])
        rho = np.concatenate([rho_o, rho_u])
        ffr, fmlt = mlt_cols(r, T, rho, r <= rph*(1.0 + 1e-12))
        hdr0 = ('# icB2-1011 (make_ic_agcar_mltpp.py --mode outward): rows r < r_join %.10e '
                '(tau 3 of the old column) copied verbatim from %s; MLT++ continued OUTWARD '
                'from the join to tau_R 2/3 (stop %s): R_ph\' %.10e (old %.10e), T_ph\' %.4f;'
                ' iso atmosphere T_ph\' (%s) from rho(R_ph\') %.6e (its own tau '
                '%.4f), floor %.1e, true tau at r_join %.4f; E above R_ph\' as the old file\n'
                % (rj, a.base, a.stop, rph, RSTAR, tph,
                   'g - kappa F/c' if a.cap_rad else 'gravity only', rhoph, tau_atm,
                   a.rho_floor, tau_rj))
    lt, lr = np.log10(T), np.log10(rho)
    eint = rho*10.0**eos.le(lt, lr)
    Fr = ffr*LUM/(4*np.pi*r**2)
    E = ecol(r, T, rph, tph)
    fn = os.path.join(a.out, a.name)
    os.makedirs(a.out, exist_ok=True)
    with open(fn, 'w') as fh:
        fh.write('# he_star_m1 IC (lbv_1008/agcar/geos/scripts/make_ic_agcar_ge.py, he_ic_cols '
                 '= 5): AG Car log L %.4f Teff %.1f M %.2f Msun, INWARD MLT++ column from '
                 'R_ph(tau 2/3) %.8e cm; alpha_MLT %.2f, MLT++ (Paxton+2013 sect 7.2, MESA '
                 'r15140 adjust_gradT_excess) alpha 1 f1 %g f2 %g; GM %.8e L %.8e; base '
                 '%.6e cm; EOS table %s; kappa_R %s; iso atmosphere T_ph %.2f gravity only, '
                 'E above R_ph diluted\n'
                 % (np.log10(LUM/M.LSUN), a.teff, GM/M.GMSUN, RSTAR, st.ALPHA, a.f1, a.f2,
                    GM, LUM, a.base_r, a.eos[6:], a.ross, Tph))
        if hdr0:
            fh.write(hdr0)
        fh.write('# r[cm] rho eint F_r(radiative part) E  | T[K] F_MLT/F (cols 6,7)\n')
        for ln in rows_keep:
            fh.write(ln)
        for row in zip(r, rho, eint, Fr, E, T, fmlt):
            fh.write('%.10e %.10e %.10e %.10e %.10e %.10e %.6e\n' % row)
    # diagnostics
    kap = 10.0**sk.ev(lt, lr)
    np.savez(fn + '.npz', r=r, rho=rho, T=T, Fr=Fr, fmlt=fmlt, kap=kap, rph=rph, rj=rj,
             rstar=RSTAR, tph=tph)
    print('wrote', fn)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--mode', choices=['inward', 'outward'], default='inward')
    p.add_argument('--base', default='/raven/ptmp/jinma/agcar_1008/files/ic_agcar_B_ge.txt')
    p.add_argument('--out', default='.')
    p.add_argument('--name', default='ic_agcar_B_ge_repro.txt')
    p.add_argument('--gm', type=float, default=55.0*M.GMSUN)
    p.add_argument('--lum', type=float, default=M.LSUN*10.0**6.17)
    p.add_argument('--teff', type=float, default=20000.0)
    p.add_argument('--rstar', type=float, default=None)
    p.add_argument('--eos', default='table:/raven/ptmp/jinma/agcar_1008/shakeB_ana/eos_dump.txt')
    p.add_argument('--ross',
                   default='/raven/ptmp/jinma/agcar_1008/files/rosseland_tops_gs98_x0.36_z0.02.txt')
    p.add_argument('--ld_ext', type=float, default=None)
    p.add_argument('--alpha', type=float, default=1.5)
    p.add_argument('--mltpp', type=int, default=1)
    p.add_argument('--f1', type=float, default=1.0e-4)
    p.add_argument('--f2', type=float, default=1.0e-3)
    p.add_argument('--tauj', type=float, default=3.0)
    p.add_argument('--dlnr', type=float, default=2.0e-4)
    p.add_argument('--base_r', type=float, default=2.490606e12)
    p.add_argument('--rin', type=float, default=0.95*2.490606e12)
    p.add_argument('--rtop', type=float, default=350.0*RSUN)
    p.add_argument('--fl_step', type=float, default=50.0)
    p.add_argument('--atm_rad', type=int, default=0)
    p.add_argument('--cap_rad', type=int, default=0,
                   help='outward: the iso atmosphere above R_ph\' with g - kappa F/c (1) or g')
    p.add_argument('--rho_floor', type=float, default=1.0e-16)
    p.add_argument('--rho_lo', type=float, default=1.0e-11)
    p.add_argument('--rho_hi', type=float, default=3.0e-9)
    p.add_argument('--stop', choices=['atm', 'internal'], default='internal',
                   help='outward: end MLT++ where the iso atmosphere on top has tau 2/3 '
                        '(atm, true tau_R 2/3) or where the old tau scale reaches 2/3')
    p.add_argument('--grid_from', default=None,
                   help='inward: take the radii below r_join from this column file')
    p.add_argument('--rtol2', type=float, default=1.0e-9)
    p.add_argument('--atol2', type=float, default=1.0e-12)
    a = p.parse_args()
    build(a)


if __name__ == '__main__':
    main()
