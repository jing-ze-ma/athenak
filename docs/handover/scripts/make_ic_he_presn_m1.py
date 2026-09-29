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
Usage : make_ic_he_presn_m1.py [out_dir]
"""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline, CubicSpline

BENCH = '/viper/u2/jinma/ATHENAK/bench'
COL = BENCH + '/hestar_presn/column_he4_presn_sph.txt'
DUMP = BENCH + '/analysis_0916/iso3d_w5/hestar/eos_table_box_w5.txt'
ROSS = ('/viper/ptmp2/jinma/caltech_handover_0926/athenak_data/he_box/'
        'rosseland_he_x0.0_z0.02.txt')
GM = 4.18143e26           # g_surf R^2 of the column header (M = 6.2651e33 g)
LUM = 2.3066e38
CL = 2.99792458e10
RSTAR = 2.3717e11


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
    x, y, le, lp = read_dump(DUMP)
    sle = RectBivariateSpline(y, x, le, kx=3, ky=3)
    slp = RectBivariateSpline(y, x, lp, kx=3, ky=3)
    ty, tx, kv = read_ross(ROSS)
    sk = RectBivariateSpline(ty, tx, kv, kx=1, ky=1)

    d = np.loadtxt(COL, comments='#')
    d = d[np.argsort(d[:, 0])]
    r, T, rho_c = d[:, 0], d[:, 2], d[:, 3]
    nab, nrad = d[:, 5], d[:, 7]
    Fr = LUM/(4*np.pi*r**2)*np.clip(nab/nrad, 0.0, 1.0)
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
                    method='LSODA', rtol=1e-10, atol=1e-12)
    assert sol.success, sol.message
    rho = np.exp(sol.y[0][::-1])
    # headroom above the column: isothermal (T = T_top) upward integration from the top
    # node, F_r = its last value * (r_top/r)^2, out to rho ~ 1e-13 (the table's floor)
    rx = r[-1] + np.arange(1, 601)*4.0e6
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
    for lim in (3e-13, 1e-13):
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
        print('  tau(table kappa, from file top) = %.3g at r = %.6e (r/R %.4f) rho %.3e T %.4e'
              % (t0, r[i], r[i]/RSTAR, rho[i], T[i]))
    hp = pg/(rho*g)
    m = (r > 0.5*RSTAR) & (r < 2.41e11)
    print('gas-only Hp = p_gas/(rho g): min %.3e at r/R %.4f (mesh range)'
          % (hp[m].min(), r[m][np.argmin(hp[m])]/RSTAR))
    os.makedirs(out, exist_ok=True)
    fn = os.path.join(out, 'ic_he_presn_m1_gas.txt')
    with open(fn, 'w') as fh:
        fh.write('# he_star_m1 IC: 4.0 Msun He star (Woosley 2019), GAS-ONLY tabulated EOS\n'
                 '# T(r) of column_he4_presn_sph.txt; rho from gas hydrostatic balance with\n'
                 '# the radiation force of F_r (make_ic_he_presn_m1.py); GM = %.6e\n'
                 '# r[cm]  rho[g/cm^3]  eint[erg/cm^3]  F_r[erg/cm^2/s] (diffusive part)\n'
                 % GM)
        for a, b, c, e in zip(r, rho, eint, Fr):
            fh.write('%.10e %.10e %.10e %.10e\n' % (a, b, c, e))
    print('wrote', fn)


if __name__ == '__main__':
    main()
