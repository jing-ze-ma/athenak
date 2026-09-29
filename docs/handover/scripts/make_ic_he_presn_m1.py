#!/usr/bin/env python3
"""IC for pgen he_star_m1: 4.0 Msun presupernova He star (Woosley 2019) for the GAS-ONLY
tabulated EOS with implicit M1 radiation (eos_radiation = false).

Input : column_he4_presn_sph.txt (bench/hestar_presn; r p T rho kappa nabla nabla_ad
        nabla_rad tau v_mlt Hp N2 Prad/Pgas), the EOS table dump of the run's own table
        (log10 e/rho, log10 p/rho on (log rho, log T)), the Rosseland table.
Method: T(r) is the column's.  The radiative flux is the DIFFUSIVE part of the column,
        F_r = L/(4 pi r^2) * min(1, nabla/nabla_rad)   (= L/(4 pi r^2) - F_MLT),
        and rho(r) is integrated inward from the top node (rho = the column's rho there)
        from gas-pressure hydrostatic balance with the radiation force of F_r,
          dp_gas(rho,T)/dr = -rho (GM/r^2 - kappa_R(rho,T) F_r/c),
        which is what the module's force_reference = wb_arad reference removes.
Output: r[cm] rho eint F_r[erg/cm^2/s]  (eint = rho 10^le(rho,T), gas only); the pgen sets
        E = a T(rho,eint)^4 itself.
        MODE rad (the STATIC radiative-equilibrium star, F_r = L/(4 pi r^2) everywhere):
        the column's T(r) carries only L - F_MLT in the FeCZ, so its F_r is not divergence
        free and a static gas heats/cools at ~1e-3/s (the missing convection).  Below the
        radius where the column's tau = JOIN_TAU the (rho, T) pair is instead integrated
        inward with the radiative-equilibrium gradient of the FULL luminosity,
          dT/dr = -3 kappa_R(rho,T) rho L/(16 pi a c r^2 T^3),
        plus the same gas hydrostatic balance (Gamma = kappa L/(4 pi r^2 c g) may exceed 1:
        density inversion).  The star is then a genuine steady state of the radiative
        transfer + hydrostatic balance (convectively unstable, which needs a seed).
        MODE rz: as rad but only the radiative zone below the FeCZ is re-integrated (the
        FeCZ keeps the column T(r), F_r = L - F_MLT; it is not stationary without
        convection).
Usage : make_ic_he_presn_m1.py [out_dir] [col|rad|rz|hopf|relax|state PREV_IC RUN_DIR K]   (file ic_he_presn_m1_<mode>.txt)
"""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline, CubicSpline

BENCH = '/viper/u2/jinma/ATHENAK/bench'
COL = BENCH + '/hestar_presn/column_he4_presn_sph.txt'
# the run's own table (<hydro>/eos_table_dump of the he_star_m1 input, logd -14..-3)
DUMP = '/viper/ptmp2/jinma/hepresn_0929/eos_table_he_x0_y0.98_z0.02_logd-14.txt'
ROSS = ('/viper/ptmp2/jinma/caltech_handover_0926/athenak_data/he_box/'
        'rosseland_he_x0.0_z0.02.txt')
GM = 4.18143e26           # g_surf R^2 of the column header (M = 6.2651e33 g)
LUM = 2.3066e38
CL = 2.99792458e10
RSTAR = 2.3717e11
HERE = os.path.dirname(os.path.abspath(__file__))
JOIN_TAU = 10.0


def read_dump(fn):
    with open(fn) as fh:
        lines = [fh.readline() for _ in range(6)]
    nx, ny, xmin, dx, ymin, dy = [float(v) for v in lines[2].split()[1:]]
    nx, ny = int(nx), int(ny)
    d = np.loadtxt(fn, comments='#')
    x = xmin + dx*np.arange(nx)
    y = ymin + dy*np.arange(ny)
    return x, y, d[:, 0].reshape(ny, nx), d[:, 1].reshape(ny, nx)


def read_ross(fn):
    with open(fn) as fh:
        for ln in fh:
            if ln.startswith('# grid'):
                pass
            if ln.startswith('# ') and len(ln.split()) == 7 and ln.split()[1].isdigit():
                nT, nD, lT0, dlT, lD0, dlD = [float(v) for v in ln.split()[1:]]
    nT, nD = int(nT), int(nD)
    v = np.loadtxt(fn, comments='#').reshape(nT, nD)
    return lT0 + dlT*np.arange(nT), lD0 + dlD*np.arange(nD), v


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    mode = sys.argv[2] if len(sys.argv) > 2 else 'col'
    x, y, le, lp = read_dump(DUMP)
    sle = RectBivariateSpline(y, x, le, kx=3, ky=3)
    slp = RectBivariateSpline(y, x, lp, kx=3, ky=3)
    ty, tx, kv = read_ross(ROSS)
    sk = RectBivariateSpline(ty, tx, kv, kx=1, ky=1)
    global sk_
    sk_ = sk
    if mode == 'state':
        return main_state(out, sys.argv[3], sys.argv[4], int(sys.argv[5]))
    if mode == 'relax':
        return main_relax(out, sys.argv[3], sys.argv[4], int(sys.argv[5]), sle, slp)

    d = np.loadtxt(COL, comments='#')
    d = d[np.argsort(d[:, 0])]
    r, T, rho_c = d[:, 0], d[:, 2], d[:, 3]
    nab, nrad = d[:, 5], d[:, 7]
    Fr = LUM/(4*np.pi*r**2)*np.clip(nab/nrad, 0.0, 1.0)
    if mode == 'hopf':
        return main_hopf(out, d, sle, slp, sk)
    lTs = CubicSpline(r, np.log(T))
    Fs = CubicSpline(r, Fr*r**2)         # smooth r^2 F

    def rhs(rr, lnrho):
        rho = np.exp(lnrho[0])
        lt = np.log10(np.exp(lTs(rr)))
        lr = np.log10(rho)
        lpv = slp.ev(lt, lr)
        dlp_dlr = slp.ev(lt, lr, dy=1)     # x = log rho is the SECOND axis
        dlp_dlt = slp.ev(lt, lr, dx=1)
        p = rho*10.0**lpv
        dlnp_dlnr = 1.0 + dlp_dlr
        dlnp_dlnT = dlp_dlt
        kap = 10.0**sk.ev(lt, lr)
        geff = GM/rr**2 - kap*Fs(rr)/rr**2/CL
        dlnT_dr = lTs(rr, 1)
        dlnp_dr = (-rho*geff/p) if p > 0 else 0.0
        drho = (dlnp_dr - dlnp_dlnT*dlnT_dr)/dlnp_dlnr
        return [drho]

    sol = solve_ivp(rhs, [r[-1], r[0]], [np.log(rho_c[-1])], t_eval=r[::-1],
                    method='LSODA', rtol=1e-9, atol=1e-12)
    assert sol.success, sol.message
    rho = np.exp(sol.y[0][::-1])
    # headroom above the column: isothermal (T = T_top) upward integration from the top
    # node, F_r = its last value * (r_top/r)^2, out to rho ~ 1e-13 (the table's floor)
    rx = r[-1] + np.arange(1, 1201)*4.0e6
    lt_top = np.log(T[-1])
    fr2_top = Fr[-1]*r[-1]**2

    def rhs_x(rr, lnrho):
        rho_ = np.exp(lnrho[0])
        lt_ = np.log10(np.exp(lt_top))
        lr_ = np.log10(rho_)
        p_ = rho_*10.0**slp.ev(lt_, lr_)
        kap_ = 10.0**sk.ev(lt_, lr_)
        geff_ = GM/rr**2 - kap_*fr2_top/rr**2/CL
        return [(-rho_*geff_/p_)/(1.0 + slp.ev(lt_, lr_, dy=1))]

    sx = solve_ivp(rhs_x, [r[-1], rx[-1]], [np.log(rho_c[-1])], t_eval=rx,
                   method='LSODA', rtol=1e-10, atol=1e-12)
    assert sx.success, sx.message
    rho_x = np.exp(sx.y[0])
    for lim in (3e-14, 1e-14):
        j = np.searchsorted(-rho_x, -lim)
        print('  headroom: rho = %.1e at r = %.6e (r/R %.4f), +%.3e cm above the column'
              % (lim, rx[min(j, len(rx) - 1)], rx[min(j, len(rx) - 1)]/RSTAR,
                 rx[min(j, len(rx) - 1)] - r[-1]))
    r_col_n = len(r)
    r = np.concatenate([r, rx])
    T = np.concatenate([T, np.full(len(rx), T[-1])])
    rho_c = np.concatenate([rho_c, np.full(len(rx), np.nan)])
    rho = np.concatenate([rho, rho_x])
    Fr = np.concatenate([Fr, fr2_top/rx**2])
    if mode in ('rad', 'rz'):
        A_RAD = 7.5657332503e-15
        tau_c = np.concatenate([d[np.argsort(d[:, 0])][:, 8], np.zeros(len(rx))])
        if mode == 'rad':
            j = int(np.argmin(abs(tau_c[:r_col_n] - JOIN_TAU)))
        else:
            # rz: only the RADIATIVE zone below the FeCZ is re-integrated: join at the
            # first node above 0.5 R where the column's diffusive fraction nab/nrad < 0.9995
            dd = d[np.argsort(d[:, 0])]
            ratio = dd[:, 5]/dd[:, 7]
            j = int(np.where((dd[:, 0] > 0.5*RSTAR) & (ratio < 0.9995))[0][0]) - 1
        Fr_col = Fr.copy()
        print('  rad mode: join at r = %.6e (r/R %.4f), column tau %.3g' %
              (r[j], r[j]/RSTAR, tau_c[j]))
        Fr = LUM/(4*np.pi*r**2) if mode == 'rad' else Fr_col

        def rhs2(rr, y):
            lnT, lnr = y
            T_ = np.exp(lnT)
            rho_ = np.exp(lnr)
            lt_ = np.log10(T_)
            lr_ = np.log10(rho_)
            kap_ = 10.0**sk.ev(lt_, lr_)
            p_ = rho_*10.0**slp.ev(lt_, lr_)
            dlnT = -3.0*kap_*rho_*LUM/(16*np.pi*A_RAD*CL*rr**2*T_**4)
            geff_ = GM/rr**2 - kap_*LUM/(4*np.pi*rr**2)/CL
            dlnp = -rho_*geff_/p_
            dlp_dlr_ = slp.ev(lt_, lr_, dy=1)
            dlp_dlt_ = slp.ev(lt_, lr_, dx=1)
            dlnr = (dlnp - dlp_dlt_*dlnT)/(1.0 + dlp_dlr_)
            return [dlnT, dlnr]

        s2 = solve_ivp(rhs2, [r[j], r[0]], [np.log(T[j]), np.log(rho[j])],
                       t_eval=r[:j+1][::-1], method='LSODA', rtol=1e-9, atol=1e-12)
        assert s2.success, s2.message
        T = T.copy()
        T[:j+1] = np.exp(s2.y[0][::-1])
        rho[:j+1] = np.exp(s2.y[1][::-1])
        print('  rad mode: T/T_col at r_in=0.5R %.4f, rho/rho_col %.4f; T(0.35R) %.4e'
              % (T[np.argmin(abs(r - 0.5*RSTAR))]/d[np.argsort(d[:, 0])][:, 2][
                  np.argmin(abs(r - 0.5*RSTAR))], (rho/rho_c)[np.argmin(
                      abs(r - 0.5*RSTAR))], T[0]))
    lt = np.log10(T)
    lr = np.log10(rho)
    eint = rho*10.0**sle.ev(lt, lr)
    pg = rho*10.0**slp.ev(lt, lr)
    kap = 10.0**sk.ev(lt, lr)
    g = GM/r**2
    gam = kap*Fr/CL/g
    print('nodes %d  r %.5e .. %.5e' % (len(r), r[0], r[-1]))
    print('rho/rho_col: min %.4f max %.4f (top %.4f, r_in=0.5R %.4f)'
          % (np.nanmin(rho/rho_c), np.nanmax(rho/rho_c), rho[r_col_n-1]/rho_c[r_col_n-1],
             (rho/rho_c)[np.argmin(abs(r - 0.5*RSTAR))]))
    print('Gamma = kappa F_r/(c g): max %.4f at r/R = %.4f' %
          (gam.max(), r[np.argmax(gam)]/RSTAR))
    # optical depth from the top of the file with the table kappa
    tau = np.zeros_like(r)
    for i in range(len(r) - 2, -1, -1):
        tau[i] = tau[i+1] + 0.5*(r[i+1] - r[i])*(rho[i]*kap[i] + rho[i+1]*kap[i+1])
    for t0 in (1e-2, 2.0/3.0, 1.0, 100.0):
        i = np.argmin(abs(tau - t0))
        print('  tau(table kappa) = %.3g at r = %.6e (r/R %.4f) rho %.3e T %.4e'
              % (t0, r[i], r[i]/RSTAR, rho[i], T[i]))
    hp = pg/(rho*g)
    m = (r > 0.5*RSTAR) & (r < 2.41e11)
    print('gas-only Hp = p_gas/(rho g): min %.3e at r/R %.4f (mesh range)'
          % (hp[m].min(), r[m][np.argmin(hp[m])]/RSTAR))
    os.makedirs(out, exist_ok=True)
    fn = os.path.join(out, 'ic_he_presn_m1_%s.txt' % mode)
    with open(fn, 'w') as fh:
        fh.write('# he_star_m1 IC: 4.0 Msun He star (Woosley 2019), GAS-ONLY table EOS\n'
                 '# T(r) of column_he4_presn_sph.txt; rho from gas hydrostatic balance\n'
                 '# with\n'
                 '# the radiation force of F_r (make_ic_he_presn_m1.py); GM = %.6e\n'
                 '# r[cm]  rho[g/cm^3]  eint[erg/cm^3]  F_r[erg/cm^2/s] (diffusive)\n'
                 % GM)
        for a, b, c, e, tt in zip(r, rho, eint, Fr, T):
            fh.write('%.10e %.10e %.10e %.10e %.10e\n' % (a, b, c, e, tt))
    print('wrote', fn)


SIG = 5.670374419e-5
A_RAD_ = 7.5657332503e-15


def q_hopf(tau):
    """Hopf function, Mihalas' fit (|err| < 1e-3)"""
    return 0.7104 - 0.1331*np.exp(-3.4488*tau)


def main_hopf(out, d, sle, slp, sk):
    """STATIC radiative-equilibrium star with the exact grey (Hopf) atmosphere:
    F_r = L/(4 pi r^2) everywhere; for tau < TAU_J the plane-parallel LTE relation
    T^4 = (3 F/4 sigma) (tau + q(tau)), below it dE/dr = -3 kappa rho F/c (diffusion),
    both with the gas hydrostatic balance.  The pair (rho, tau) is pinned at the column's
    top node (rho = the column's, tau = the column's tau_top = 1e-4)."""
    d = d[np.argsort(d[:, 0])]
    rc, rho_c, tau_c = d[:, 0], d[:, 3], d[:, 8]
    r_pin, rho_pin, tau_pin = rc[-1], rho_c[-1], tau_c[-1]
    TAU_J = JOIN_TAU
    rx = r_pin + np.arange(1, 1201)*4.0e6

    def kap_of(T_, rho_):
        return 10.0**sk.ev(np.log10(T_), np.log10(rho_))

    def rhs1(rr, y):
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
        p_ = rho_*10.0**slp.ev(lt_, lr_)
        geff = GM/rr**2 - kap_*F/CL
        dlnp = -rho_*geff/p_
        dlnr = (dlnp - slp.ev(lt_, lr_, dx=1)*dlnT)/(1.0 + slp.ev(lt_, lr_, dy=1))
        return [dlnr, dtau]

    def rhs2(rr, y):
        lnT, lnr = y
        T_, rho_ = np.exp(lnT), np.exp(lnr)
        lt_, lr_ = np.log10(T_), np.log10(rho_)
        kap_ = kap_of(T_, rho_)
        p_ = rho_*10.0**slp.ev(lt_, lr_)
        dlnT = -3.0*kap_*rho_*LUM/(16*np.pi*A_RAD_*CL*rr**2*T_**4)
        geff = GM/rr**2 - kap_*LUM/(4*np.pi*rr**2)/CL
        dlnp = -rho_*geff/p_
        dlnr = (dlnp - slp.ev(lt_, lr_, dx=1)*dlnT)/(1.0 + slp.ev(lt_, lr_, dy=1))
        return [dlnT, dlnr]

    # upward from the pin (region 1)
    su = solve_ivp(rhs1, [r_pin, rx[-1]], [np.log(rho_pin), tau_pin], t_eval=rx,
                   method='LSODA', rtol=1e-9, atol=1e-13)
    assert su.success, su.message
    # downward from the pin to tau = TAU_J (event), then region 2
    def ev(rr, y):
        return y[1] - TAU_J
    ev.terminal = True
    nodes = rc[::-1]
    s1 = solve_ivp(rhs1, [r_pin, rc[0]], [np.log(rho_pin), tau_pin], t_eval=nodes,
                   method='LSODA', rtol=1e-9, atol=1e-13, events=ev)
    assert s1.success, s1.message
    rj = s1.t_events[0][0]
    yj = s1.y_events[0][0]
    Fj = LUM/(4*np.pi*rj**2)
    Tj = (3.0*Fj/(4.0*SIG)*(yj[1] + q_hopf(yj[1])))**0.25
    print('  hopf: join at tau = %.1f, r = %.6e (r/R %.4f), T %.4e rho %.4e'
          % (TAU_J, rj, rj/RSTAR, Tj, np.exp(yj[0])))
    n1 = s1.t.size
    r1 = s1.t
    tau1 = s1.y[1]
    rho1 = np.exp(s1.y[0])
    F1 = LUM/(4*np.pi*r1**2)
    T1 = (3.0*F1/(4.0*SIG)*(np.maximum(tau1, 1e-8) + q_hopf(np.maximum(tau1, 1e-8))))**0.25
    rest = nodes[nodes < r1[-1]]
    s2 = solve_ivp(rhs2, [rj, rc[0]], [np.log(Tj), yj[0]], t_eval=rest,
                   method='LSODA', rtol=1e-9, atol=1e-13)
    assert s2.success, s2.message
    # assemble ascending
    Fu = LUM/(4*np.pi*rx**2)
    tau_u = np.maximum(su.y[1], 1e-8)
    Tu = (3.0*Fu/(4.0*SIG)*(tau_u + q_hopf(tau_u)))**0.25
    r_all = np.concatenate([s2.t[::-1], r1[::-1], rx])
    T_all = np.concatenate([np.exp(s2.y[0])[::-1], T1[::-1], Tu])
    rho_all = np.concatenate([np.exp(s2.y[1])[::-1], rho1[::-1], np.exp(su.y[0])])
    o = np.argsort(r_all)
    r, T, rho = r_all[o], T_all[o], rho_all[o]
    Fr = LUM/(4*np.pi*r**2)
    lt, lr = np.log10(T), np.log10(rho)
    eint = rho*10.0**sle.ev(lt, lr)
    pg = rho*10.0**slp.ev(lt, lr)
    kap = 10.0**sk.ev(lt, lr)
    g = GM/r**2
    tau = np.zeros_like(r)
    for i in range(len(r) - 2, -1, -1):
        tau[i] = tau[i+1] + 0.5*(r[i+1] - r[i])*(rho[i]*kap[i] + rho[i+1]*kap[i+1])
    print('nodes %d  r %.5e .. %.5e' % (len(r), r[0], r[-1]))
    gam = kap*Fr/CL/g
    print('Gamma = kappa F/(c g): max %.4f at r/R = %.4f' % (gam.max(), r[np.argmax(gam)]/RSTAR))
    i5 = np.argmin(abs(r - 0.5*RSTAR))
    print('r_in = 0.5R: T %.4e rho %.4e ; T(0.35R) %.4e' % (T[i5], rho[i5], T[0]))
    for t0 in (1e-2, 2.0/3.0, 1.0, 100.0):
        i = np.argmin(abs(tau - t0))
        print('  tau(table kappa) = %.3g at r = %.6e (r/R %.4f) rho %.3e T %.4e'
              % (t0, r[i], r[i]/RSTAR, rho[i], T[i]))
    for lim in (3e-14, 1e-14):
        i = np.searchsorted(-rho, -lim)
        print('  rho = %.1e at r = %.6e (r/R %.4f)' % (lim, r[min(i, len(r)-1)],
                                                     r[min(i, len(r)-1)]/RSTAR))
    os.makedirs(out, exist_ok=True)
    fn = os.path.join(out, 'ic_he_presn_m1_hopf.txt')
    with open(fn, 'w') as fh:
        fh.write('# he_star_m1 IC (mode hopf): STATIC radiative-equilibrium star, F_r = L/(4 pi r^2)\n'
                 '# everywhere, grey (Hopf) LTE atmosphere above tau = %g, GAS-ONLY table EOS;\n'
                 '# see make_ic_he_presn_m1.py.  GM = %.6e\n'
                 '# r[cm]  rho[g/cm^3]  eint[erg/cm^3]  F_r[erg/cm^2/s]  T[K] (T is informational)\n' % (TAU_J, GM))
        for a, b, c, e, tt in zip(r, rho, eint, Fr, T):
            fh.write('%.10e %.10e %.10e %.10e %.10e\n' % (a, b, c, e, tt))
    print('wrote', fn)


def main_relax(out, prev_fn, run_dir, k, sle, slp):
    """One fixed-point step of the discrete M1 equilibrium: the radiation temperature
    T_E = (E/a)^(1/4) of the output number k of a column run started from the previous
    IC replaces the previous T(r) (the ratio is applied on the mesh, constant outside);
    rho is re-integrated from gas hydrostatic balance with T_E and F_r = L/(4 pi r^2)."""
    sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'vis', 'python'))
    from bin_convert import read_binary
    pv = np.loadtxt(prev_fn, comments='#')
    r, rho0, e0, F0, T0 = pv[:, 0], pv[:, 1], pv[:, 2], pv[:, 3], pv[:, 4]
    f = read_binary(os.path.join(run_dir, 'bin', 'hepresn.m1.%05d.bin' % k))
    n = f['Nx1']
    rc = f['x1min'] + (np.arange(n) + 0.5)*(f['x1max'] - f['x1min'])/n
    E = np.asarray(f['mb_data']['m1_e'])[0].mean(axis=(0, 1))
    Te = (E/A_RAD_)**0.25
    ratio = Te/np.interp(rc, r, T0)
    T = T0*np.interp(r, rc, ratio)          # clamped outside the mesh
    print('  relax: t = %.2f s, T_E/T_prev - 1 on the mesh: min %.3e max %.3e'
          % (f['time'], (ratio - 1).min(), (ratio - 1).max()))
    Fr = LUM/(4*np.pi*r**2)
    lTs = CubicSpline(r, np.log(T))
    Fs = CubicSpline(r, Fr*r**2)

    def rhs(rr, lnrho):
        rho = np.exp(lnrho[0])
        lt = np.log10(np.exp(lTs(rr)))
        lr = np.log10(rho)
        p = rho*10.0**slp.ev(lt, lr)
        kap = 10.0**sk_.ev(lt, lr)
        geff = GM/rr**2 - kap*Fs(rr)/rr**2/CL
        dlnp_dr = -rho*geff/p
        dlp_dlr = slp.ev(lt, lr, dy=1)
        dlp_dlt = slp.ev(lt, lr, dx=1)
        return [(dlnp_dr - dlp_dlt*lTs(rr, 1))/(1.0 + dlp_dlr)]

    ip = int(np.argmin(abs(r - 2.4325e11)))
    sd = solve_ivp(rhs, [r[ip], r[0]], [np.log(rho0[ip])], t_eval=r[:ip+1][::-1],
                   method='LSODA', rtol=1e-9, atol=1e-12)
    su = solve_ivp(rhs, [r[ip], r[-1]], [np.log(rho0[ip])], t_eval=r[ip:],
                   method='LSODA', rtol=1e-9, atol=1e-12)
    assert sd.success and su.success
    rho = np.concatenate([np.exp(sd.y[0])[::-1], np.exp(su.y[0])[1:]])
    eint = rho*10.0**sle.ev(np.log10(T), np.log10(rho))
    print('  relax: rho/rho_prev - 1: min %.3e max %.3e' % ((rho/rho0 - 1).min(),
                                                          (rho/rho0 - 1).max()))
    fn = os.path.join(out, os.path.basename(prev_fn).replace('.txt', '') + '_r%d.txt' % k)
    with open(fn, 'w') as fh:
        fh.write('# he_star_m1 IC, relaxed one step from %s using %s output %d\n'
                 '# r[cm]  rho[g/cm^3]  eint[erg/cm^3]  F_r[erg/cm^2/s]  T[K]\n'
                 % (os.path.basename(prev_fn), run_dir, k))
        for a, b, c, e, tt in zip(r, rho, eint, Fr, T):
            fh.write('%.10e %.10e %.10e %.10e %.10e\n' % (a, b, c, e, tt))
    print('wrote', fn)


def main_state(out, prev_fn, run_dir, k):
    """The state of output k of a (velocity-damped) column run as the IC: rho, eint (total
    energy minus kinetic minus rho Phi), F_r and E at the cell centres; outside the mesh
    the previous IC scaled by the edge ratio.  Written with 5 columns (he_ic_cols = 5).
    It is a fixed point of the DISCRETE scheme at that radial resolution (nx1), so it is a
    test IC for the thin-column gate, not a production one."""
    sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'vis', 'python'))
    from bin_convert import read_binary
    pv = np.loadtxt(prev_fn, comments='#')
    r, rho0, e0, F0, T0 = pv[:, 0], pv[:, 1], pv[:, 2], pv[:, 3], pv[:, 4]
    E0 = A_RAD_*T0**4
    fu = read_binary(os.path.join(run_dir, 'bin', 'hepresn.hydro_u.%05d.bin' % k))
    fm = read_binary(os.path.join(run_dir, 'bin', 'hepresn.m1.%05d.bin' % k))
    n = fu['Nx1']
    dx = (fu['x1max'] - fu['x1min'])/n
    rc = fu['x1min'] + (np.arange(n) + 0.5)*dx
    md = fu['mb_data']
    ave = lambda a: np.asarray(a)[0].mean(axis=(0, 1))
    rho, mom, ener = ave(md['dens']), ave(md['mom1']), ave(md['ener'])
    eint = ener - 0.5*mom**2/rho - rho*GM*(1.0/fu['x1min'] - 1.0/rc)
    E = ave(fm['mb_data']['m1_e'])
    F = ave(fm['mb_data']['m1_f1'])
    print('  state: t = %.1f s, max|v| %.3e cm/s' % (fu['time'], np.abs(mom/rho).max()))

    def ext(a_state, a_prev, lo):
        ref = np.interp(rc[0] if lo else rc[-1], r, a_prev)
        return a_prev*(a_state[0 if lo else -1]/ref)
    low, high = r < rc[0], r > rc[-1]
    nr = np.concatenate([r[low], rc, r[high]])
    out_cols = []
    for st, pr in ((rho, rho0), (eint, e0), (F, F0), (E, E0)):
        out_cols.append(np.concatenate([ext(st, pr, True)[low], st, ext(st, pr, False)[high]]))
    fn = os.path.join(out, os.path.basename(prev_fn).replace('.txt', '') + '_state%d.txt' % k)
    with open(fn, 'w') as fh:
        fh.write('# he_star_m1 IC (he_ic_cols = 5), the state of output %d of %s\n'
                 '# r[cm]  rho  eint  F_r  E\n' % (k, run_dir))
        for row in zip(nr, *out_cols):
            fh.write('%.10e %.10e %.10e %.10e %.10e\n' % row)
    print('wrote', fn)


if __name__ == '__main__':
    main()
