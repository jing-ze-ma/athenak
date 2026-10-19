#!/usr/bin/env python3
"""S0 of the accretor RHD port (docs/dev/accretor_rhd_design.md section 1.1): the gainer's
radiative envelope as a 1-D column in the Roche equipotential coordinate
    psi = Phi_s - Phi_wb   [(km/s)^2],  Phi_s = Phi(R_acc, phi = 90 deg),
built with the run's OWN general-EOS table (the <hydro>/eos_table_dump of the run's
<hydro> block) and a Rosseland table.  Read by ry_per_accretor (problem/thermo = general,
problem/env_ic_file): columns 1-3 (psi, T, rho) are used; the rest are diagnostics.

Physics (grey, LTE, plane-parallel atmosphere joined to diffusion, one pass):
  * geometry along the phi = 90 deg cut of the Roche potential of the pgen (spin 1,
    Phi_wb = Phi below the envelope top): r_eq(psi), g = dPhi/dr;
  * constant luminosity L = 4 pi R_acc^2 sigma T_eff^4 (R_acc = the phi = 90 radius of the
    photospheric equipotential), F = L/(4 pi r_eq^2);
  * hydrostatic balance with gas + radiation pressure, dP_tot/dpsi = rho;
  * radiative transfer dP_rad/dpsi = kappa_R rho F/(c g) (exact diffusion; in the
    plane-parallel top it IS the Eddington grey atmosphere T^4 = 3/4 T_eff^4 (tau + 2/3),
    started at tau0 with P_rad(tau0) = a T^4/3 of that law);
  * psi_start of the top is shot so that T = T_eff (tau = 2/3) falls on psi = 0;
  * above the start point: isothermal at T(tau0), gas-pressure hydrostatic.
Convection is NOT included: the Schwarzschild check (nabla_rad vs nabla_ad with gas +
radiation) is only reported.

Usage: make_ic_accretor_column.py EOS_DUMP ROSS_TABLE OUT [key=value ...]
  keys: m_acc m_don a_sep r_acc t_eff tau0 r_deep psi_top_cph2 cs_ph
Output: OUT (ascii): psi T rho P_gas P_rad kappa_R r_eq F Gamma_edd nabla nabla_ad
"""
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import brentq

# the pgen's constants (ry_per_accretor.cpp)
kG, kMsun, kRsun, kVel = 6.674e-8, 1.989e33, 6.957e10, 1.0e5
A_RAD, C_L, SIG = 7.5657e-15, 2.99792458e10, 5.6704e-5
LSUN = 3.828e33


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
    hdr = None
    with open(fn) as fh:
        for ln in fh:
            if not ln.startswith('#'):
                break
            w = ln.split()
            if len(w) == 7 and w[1].isdigit():
                hdr = [float(v) for v in w[1:]]
    nT, nD, lT0, dlT, lD0, dlD = hdr
    nT, nD = int(nT), int(nD)
    v = np.loadtxt(fn, comments='#').reshape(nT, nD)
    return lT0 + dlT*np.arange(nT), lD0 + dlD*np.arange(nD), v


def main():
    eosf, rossf, out = sys.argv[1:4]
    kv = dict(a.split('=') for a in sys.argv[4:])
    ma = float(kv.get('m_acc', 16.0))
    md = float(kv.get('m_don', 18.0))
    asep = float(kv.get('a_sep', 32.5575))
    racc = float(kv.get('r_acc', 9.00129260997162))
    teff = float(kv.get('t_eff', 28300.0))
    tau0 = float(kv.get('tau0', 1.0e-5))
    rdeep = float(kv.get('r_deep', 5.9))
    cph = float(kv.get('cs_ph', 19.4132))
    psitop = -float(kv.get('psi_top_cph2', 20.0))*cph*cph

    gu = kG*kMsun/(kRsun*kVel*kVel)
    gma, gmd = gu*ma, gu*md
    om = np.sqrt(gu*(ma + md)/asep**3)
    xcm = asep*md/(ma + md)

    def phi90(r):
        return -gma/r - gmd/np.sqrt(r*r + asep*asep) - 0.5*om*om*(r*r + xcm*xcm)

    def dphi90(r):           # (km/s)^2 / Rsun
        return gma/r**2 + gmd*r/(r*r + asep*asep)**1.5 - om*om*r

    phis = phi90(racc)

    def psi_of_r(r):
        return phis - phi90(r)

    def r_of_psi(psi):
        return brentq(lambda r: psi_of_r(r) - psi, 0.3*racc, 1.45*racc, xtol=1e-14)

    lum = 4.0*np.pi*(racc*kRsun)**2*SIG*teff**4

    # EOS and opacity
    xe, ye, le, lp = read_dump(eosf)
    sp_p = RectBivariateSpline(ye, xe, lp, kx=3, ky=3)    # log10 p/rho (logT, logrho)
    sp_e = RectBivariateSpline(ye, xe, le, kx=3, ky=3)
    lTo, lDo, ko = read_ross(rossf)
    sp_k = RectBivariateSpline(lTo, lDo, ko, kx=1, ky=1)

    def kap(rho, T):
        lt = np.clip(np.log10(T), lTo[0], lTo[-1])
        ld = np.clip(np.log10(rho), lDo[0], lDo[-1])
        return 10.0**sp_k(lt, ld, grid=False)

    def pgas(rho, T):
        return rho*10.0**sp_p(np.log10(T), np.log10(rho), grid=False)

    def rho_of(pg, T):
        lt = np.log10(T)
        f = lambda ld: ld + sp_p(lt, ld, grid=False) - np.log10(pg)
        return 10.0**brentq(f, xe[0], xe[-1], xtol=1e-13)

    def rhs(psi, y):
        pg, pr = np.exp(y)
        T = (3.0*pr/A_RAD)**0.25
        rho = rho_of(pg, T)
        r = r_of_psi(psi)
        g = dphi90(r)*kVel*kVel/kRsun                       # cm/s^2
        F = lum/(4.0*np.pi*(r*kRsun)**2)
        dprad = kap(rho, T)*rho*F/(C_L*g)*kVel*kVel          # per unit psi [(km/s)^2]
        dptot = rho*kVel*kVel
        return [(dptot - dprad)/pg, dprad/pr]

    def start_state(psi0):
        r = r_of_psi(psi0)
        g = dphi90(r)*kVel*kVel/kRsun
        T0 = (0.75*teff**4*(tau0 + 2.0/3.0))**0.25
        F = lum/(4.0*np.pi*(r*kRsun)**2)
        pr = A_RAD*T0**4/3.0
        pg = 1.0e2
        for _ in range(50):                    # tau0 = kappa P_gas/(g (1 - Gamma))
            rho = rho_of(pg, T0)
            k = kap(rho, T0)
            gam = k*F/(C_L*g)
            pg = tau0*g*(1.0 - gam)/k
        return [np.log(pg), np.log(pr)], T0

    def teff_event(psi, y):
        return (3.0*np.exp(y[1])/A_RAD)**0.25 - teff
    teff_event.terminal = True

    def shoot(psi0):
        y0, _ = start_state(psi0)
        s = solve_ivp(rhs, (psi0, 2.0e3), y0, method='LSODA', rtol=1e-10, atol=1e-12,
                      events=teff_event)
        if len(s.t_events[0]) == 0:      # photosphere deeper than the window
            return 2.0e3 + (teff - (3.0*np.exp(s.y[1][-1])/A_RAD)**0.25)
        return s.t_events[0][0]

    # psi_start such that T = T_eff at psi = 0 (the photospheric equipotential)
    a, b = -8000.0, -20.0
    psi0 = brentq(shoot, a, b, xtol=1e-6)
    psi_deep = psi_of_r(rdeep)
    y0, T0 = start_state(psi0)
    # output grid: fine through the atmosphere/photosphere, geometric below
    g1 = np.arange(np.ceil(psi0), 4000.0, 1.0)
    g2 = np.geomspace(4000.0, psi_deep, 12000)[1:]
    grid = np.concatenate([[psi0], g1[g1 > psi0], g2])
    s = solve_ivp(rhs, (psi0, psi_deep), y0, method='LSODA', rtol=1e-10, atol=1e-12,
                  t_eval=grid, dense_output=False)
    if not s.success:
        raise SystemExit('integration failed: ' + s.message)
    psi = s.t
    pg, pr = np.exp(s.y)
    T = (3.0*pr/A_RAD)**0.25
    rho = np.array([rho_of(a_, b_) for a_, b_ in zip(pg, T)])
    # top: isothermal at T0, gas-only hydrostatic, psi_top .. psi0
    pt = np.arange(psitop, psi0, 2.0)
    # (ideal-gas isothermal: the gas there is ionized; floored at 1e-17 g/cc, below the
    # hot ambient's pressure anyway, so only T of this part is used by the pgen)
    c2 = pg[0]/rho[0]
    rhot = np.maximum(rho[0]*np.exp((pt - psi0)*kVel*kVel/c2), 1.0e-17)
    pgt = rhot*c2
    psi = np.concatenate([pt, psi])
    T = np.concatenate([np.full(len(pt), T0), T])
    rho = np.concatenate([rhot, rho])
    pg = np.concatenate([pgt, pg])
    pr = np.concatenate([np.full(len(pt), pr[0]), pr])
    r = np.array([r_of_psi(q) for q in psi])
    F = lum/(4.0*np.pi*(r*kRsun)**2)
    kR = kap(rho, T)
    g = dphi90(r)*kVel*kVel/kRsun
    gam = kR*F/(C_L*g)
    # gradients and the adiabatic gradient (gas + radiation), numerically from the table
    ptot = pg + pr
    with np.errstate(divide='ignore', invalid='ignore'):
        nab = np.gradient(np.log(T), np.log(ptot))
    nab = np.where(np.isfinite(nab), nab, 0.0)     # the isothermal top: P const at the floor

    def nabla_ad(rh, t):
        h = 1e-4

        def u(rr, tt):       # specific energy, gas + radiation
            return 10.0**sp_e(np.log10(tt), np.log10(rr), grid=False) + A_RAD*tt**4/rr

        def P(rr, tt):
            return pgas(rr, tt) + A_RAD*tt**4/3.0
        u_r = (u(rh*(1+h), t) - u(rh*(1-h), t))/(2*h*rh)
        u_t = (u(rh, t*(1+h)) - u(rh, t*(1-h)))/(2*h*t)
        p_r = (P(rh*(1+h), t) - P(rh*(1-h), t))/(2*h*rh)
        p_t = (P(rh, t*(1+h)) - P(rh, t*(1-h)))/(2*h*t)
        p0 = P(rh, t)
        dtdr = (p0/rh**2 - u_r)/u_t
        return (dtdr/t)/((p_r + p_t*dtdr)/p0)
    nad = np.array([nabla_ad(a_, b_) for a_, b_ in zip(rho, T)])

    # checks
    i0 = np.argmin(np.abs(psi))
    hres = np.gradient(ptot, psi*kVel*kVel)/rho - 1.0
    sel = (psi > psi0 + 10) & (psi < psi_deep*0.999)
    conv = sel & (nab > nad) & (psi > 0)
    print('# accretor column: M_a %.3f M_d %.3f a %.4f Rsun, R_acc %.6f (phi 90), T_eff %.0f K'
          % (ma, md, asep, racc, teff))
    print('# L = %.4e erg/s = %.4e Lsun (log %.3f)' % (lum, lum/LSUN, np.log10(lum/LSUN)))
    print('# psi_start %.3f (tau0 %.1e, T0 %.1f K), psi_deep %.1f (r %.3f)'
          % (psi0, tau0, T0, psi_deep, rdeep))
    print('# photosphere (psi 0): r %.6f T %.1f K rho %.4e g/cc (env13: 2.82e-9), Gamma_edd %.4f'
          % (r[i0], T[i0], rho[i0], gam[i0]))
    for rr in (9.0, 8.5, 8.0, 7.5, 7.0, 6.5, 6.3, 6.0):
        j = np.argmin(np.abs(r - rr))
        print('#   r %.3f psi %.5g T %.4e rho %.4e Prad/Pgas %.4f kappa %.4f Gamma %.3f '
              'nabla %.4f nabla_ad %.4f' % (r[j], psi[j], T[j], rho[j], pr[j]/pg[j], kR[j],
                                           gam[j], nab[j], nad[j]))
    print('# HSE residual |dP_tot/dpsi/rho - 1| max (finite differences on the output grid) '
          '%.3e (median %.3e)' % (np.max(np.abs(hres[sel])), np.median(np.abs(hres[sel]))))
    print('# smoothness: max |d ln T| between rows %.4f, max |d ln rho| %.4f'
          % (np.max(np.abs(np.diff(np.log(T)))), np.max(np.abs(np.diff(np.log(rho))))))
    if conv.any():
        segs = []
        idx = np.where(conv)[0]
        start = idx[0]
        for a_, b_ in zip(idx[:-1], idx[1:]):
            if b_ != a_ + 1:
                segs.append((start, a_))
                start = b_
        segs.append((start, idx[-1]))
        for a_, b_ in segs:
            print('# SCHWARZSCHILD UNSTABLE: r %.4f .. %.4f, T %.3g .. %.3g K, max nabla-nabla_ad'
                  ' %.4f' % (r[b_], r[a_], T[a_], T[b_], np.max(nab[a_:b_+1] - nad[a_:b_+1])))
    else:
        print('# Schwarzschild: stable everywhere below the photosphere')
    hdr = ('accretor envelope column (make_ic_accretor_column.py): EOS %s, opacity %s\n'
           'M_a %.4g M_d %.4g a %.6g R_acc %.10g T_eff %.1f L %.6e erg/s tau0 %.1e\n'
           'psi[(km/s)^2] T[K] rho[g/cc] P_gas P_rad[dyn/cm^2] kappa_R[cm^2/g] r_eq[Rsun] '
           'F[erg/cm^2/s] Gamma_edd nabla nabla_ad' % (eosf, rossf, ma, md, asep, racc, teff,
                                                    lum, tau0))
    np.savetxt(out, np.column_stack([psi, T, rho, pg, pr, kR, r, F, gam, nab, nad]),
               header=hdr, fmt='%.12e')


if __name__ == '__main__':
    main()
