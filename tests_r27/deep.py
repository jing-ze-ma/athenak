#!/usr/bin/env python3
"""tests_r27: size a DEEPER inner wall for the 3-D He4 presupernova wedge.

Reads the two 1-D IC columns (tests_r14 MLT model, tests_r23 radiative model), the
general-EOS table dump and the Rosseland table, and prints, for a list of candidate
r_in, everything needed to choose one: enclosed mass, thermal time of the added
layer, T, rho, beta, tau, c_s, and the transverse CFL at 96^2 / 192^2.

usage: deep.py [r_in/R ...]
"""
import sys

import numpy as np

B = '/viper/u2/jinma/ATHENAK/bench/wt_he4'
ICM = B + '/tests_r14/ic_he4_tall_neww.txt'
ICR = B + '/tests_r23/ic_rad_f100.txt'
EOSTBL = ('/viper/u2/jinma/ATHENAK/bench/analysis_0916/iso3d_w5/hestar/'
          'eos_table_box_w5.txt')
OPAC = '/viper/u2/jinma/ATHENAK/bench/hestar_fecz/rosseland_he_x0.0_z0.02.txt'

RS = 2.3717e11
LSTAR = 2.3066e38
MSUN = 1.989e33
TURN = 4705.0
ARAD = 7.5657332503e-15
RAD_LO, RAD_HI = np.log10(1.0e-11), np.log10(1.0e-10)


def load_tbl(fn):
    with open(fn) as fh:
        ln = [fh.readline() for _ in range(6)]
    nx, ny, xmin, dx, ymin, dy = [float(v) for v in ln[2].split()[1:]]
    nx, ny = int(nx), int(ny)
    d = np.loadtxt(fn, comments='#')
    return dict(nx=nx, ny=ny, xmin=xmin, dx=dx, ymin=ymin, dy=dy,
                E=d[:, 0].reshape(ny, nx), P=d[:, 1].reshape(ny, nx))


TB = load_tbl(EOSTBL)


def bil(F, x, y):
    gx = np.clip((x - TB['xmin'])/TB['dx'], 0, TB['nx']-1.001)
    gy = np.clip((y - TB['ymin'])/TB['dy'], 0, TB['ny']-1.001)
    ix = np.floor(gx).astype(np.int32)
    iy = np.floor(gy).astype(np.int32)
    u, v = gx-ix, gy-iy
    return ((1-u)*(1-v)*F[iy, ix] + u*(1-v)*F[iy, ix+1]
            + (1-u)*v*F[iy+1, ix] + u*v*F[iy+1, ix+1])


def weight(x):
    s = (x - RAD_LO)/(RAD_HI - RAD_LO)
    return np.where(s <= 0, 0.0, np.where(s >= 1, 1.0, s*s*(3-2*s)))


def temp_of(rho, eint):
    x = np.log10(np.maximum(rho, 1e-300))
    w = weight(x)
    lo = np.full(np.shape(rho), 2.5)
    hi = np.full(np.shape(rho), 7.0)
    for _ in range(60):
        mid = 0.5*(lo+hi)
        f = rho*10.0**bil(TB['E'], x, mid) + w*ARAD*(10.0**mid)**4 - eint
        lo = np.where(f < 0.0, mid, lo)
        hi = np.where(f < 0.0, hi, mid)
    return 10.0**(0.5*(lo+hi))


def thermo(rho, T):
    """gas + (weighted) radiation: p_gas, p_tot, e_tot, Gamma1, cs."""
    x = np.log10(rho)
    y = np.log10(T)
    w = weight(x)
    pg = rho*10.0**bil(TB['P'], x, y)
    eg = rho*10.0**bil(TB['E'], x, y)
    pr = w*ARAD*T**4/3.0
    er = w*ARAD*T**4
    p, e = pg+pr, eg+er
    h = 0.01
    # chi_rho, chi_T on the GAS part from the table; radiation analytic
    dpg_dlr = (rho*10.0**bil(TB['P'], x+h, y) - rho*10.0**bil(TB['P'], x-h, y))/(2*h)
    dpg_dlt = (rho*10.0**bil(TB['P'], x, y+h) - rho*10.0**bil(TB['P'], x, y-h))/(2*h)
    deg_dlt = (rho*10.0**bil(TB['E'], x, y+h) - rho*10.0**bil(TB['E'], x, y-h))/(2*h)
    # note: bil(P) is log10(p/rho) so rho*10**.. at rho fixed; the rho-derivative
    # must also carry the explicit rho prefactor -> add p itself
    chi_rho = (dpg_dlr + pg)/p
    chi_T = (dpg_dlt + 4.0*pr)/p
    cv = (deg_dlt + 4.0*er)/T          # de/dT at constant rho (per volume)
    g1 = chi_rho + chi_T**2*p/(T*cv)
    return pg, p, e, g1, np.sqrt(g1*p/rho)


def load_opac(fn):
    with open(fn) as fh:
        ln = [fh.readline() for _ in range(8)]
    nT, nD, lTmin, dlT, lDmin, dlD = [float(v) for v in ln[6].split()[1:]]
    nT, nD = int(nT), int(nD)
    K = np.loadtxt(fn, comments='#').reshape(nT, nD)
    return dict(nT=nT, nD=nD, lTmin=lTmin, dlT=dlT, lDmin=lDmin, dlD=dlD, K=K)


OP = load_opac(OPAC)


def kappa(rho, T):
    gx = np.clip((np.log10(T) - OP['lTmin'])/OP['dlT'], 0, OP['nT']-1.001)
    gy = np.clip((np.log10(rho) - OP['lDmin'])/OP['dlD'], 0, OP['nD']-1.001)
    ix = np.floor(gx).astype(np.int32)
    iy = np.floor(gy).astype(np.int32)
    u, v = gx-ix, gy-iy
    K = OP['K']
    return 10.0**((1-u)*(1-v)*K[ix, iy] + u*(1-v)*K[ix+1, iy]
                  + (1-u)*v*K[ix, iy+1] + u*v*K[ix+1, iy+1])


def col(fn):
    d = np.loadtxt(fn, comments='#')
    r, rho, e = d[:, 0], d[:, 1], d[:, 2]
    T = temp_of(rho, e)
    pg, p, et, g1, cs = thermo(rho, T)
    ka = kappa(rho, T)
    # tau from the top down
    dr = np.diff(r)
    dtau = 0.5*(ka[1:]*rho[1:] + ka[:-1]*rho[:-1])*dr
    tau = np.concatenate((np.cumsum(dtau[::-1])[::-1], [0.0]))
    # enclosed quantities, integrated outward
    dm = 4*np.pi*(0.5*(r[1:]+r[:-1]))**2*0.5*(rho[1:]+rho[:-1])*dr
    de = 4*np.pi*(0.5*(r[1:]+r[:-1]))**2*0.5*(e[1:]+e[:-1])*dr
    return dict(r=r, rho=rho, e=e, T=T, pg=pg, p=p, g1=g1, cs=cs, ka=ka, tau=tau,
                cm=np.concatenate(([0.0], np.cumsum(dm))),
                ce=np.concatenate(([0.0], np.cumsum(de))))


def at(c, r):
    return {k: np.interp(r, c['r'], v) for k, v in c.items()
            if k != 'r' and len(v) == len(c['r'])}


def between(c, key, ra, rb):
    return np.interp(rb, c['r'], c[key]) - np.interp(ra, c['r'], c[key])


if __name__ == '__main__':
    fr = [float(v) for v in sys.argv[1:]] or [0.50, 0.45, 0.42, 0.40, 0.38, 0.365]
    for nm, fn in (('MLT  (r14)', ICM), ('RAD  (r23)', ICR)):
        c = col(fn)
        print('=== %s : file spans %.5e .. %.5e cm (%.4f .. %.4f R)'
              % (nm, c['r'][0], c['r'][-1], c['r'][0]/RS, c['r'][-1]/RS))
        m_env = between(c, 'cm', 0.5*RS, 1.25*RS)
        e_env = between(c, 'ce', 0.5*RS, 1.25*RS)
        print('    present envelope 0.50-1.25 R: M = %.4e g = %.4f Msun, '
              'E_int = %.4e erg, E/L = %.4f turnover'
              % (m_env, m_env/MSUN, e_env, e_env/LSTAR/TURN))
        print('  %6s %10s %10s %10s %9s %8s %8s %9s %9s %9s %8s %8s %8s %8s'
              % ('r/R', 'r[cm]', 'T[K]', 'rho', 'logrho', 'logT', 'beta', 'tau',
                 'cs[cm/s]', 'kappa', 'dM[Ms]', 'dM/Menv', 'dE/L[to]', 'dt96'))
        for f in fr:
            r = f*RS
            a = at(c, r)
            beta = a['pg']/a['p']
            dm = between(c, 'cm', r, 0.5*RS)
            de = between(c, 'ce', r, 0.5*RS)
            dth96 = (np.pi/2)/96
            dth192 = (np.pi/2)/192
            dt96 = 0.3*r*dth96/a['cs']
            dt192 = 0.3*r*dth192/a['cs']
            print('  %6.3f %10.4e %10.4e %10.3e %9.3f %8.3f %8.4f %9.3e %9.3e '
                  '%9.3e %8.4f %8.2f %8.3f %8.2f/%.2f'
                  % (f, r, a['T'], a['rho'], np.log10(a['rho']), np.log10(a['T']),
                     beta, a['tau'], a['cs'], a['ka'], dm/MSUN, dm/m_env,
                     de/LSTAR/TURN, dt96, dt192))
        print('    EOS table covers log rho [%.1f, %.1f], log T [%.2f, %.2f]'
              % (TB['xmin'], TB['xmin']+(TB['nx']-1)*TB['dx'],
                 TB['ymin'], TB['ymin']+(TB['ny']-1)*TB['dy']))
        print()
