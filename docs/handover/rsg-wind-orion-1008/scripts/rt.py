"""(b) 1-D static spherical RT per (band, g-point) through RSG atmospheres / winds.

Formal solution: short characteristics (linear source, log-mean chi) along p-z rays
(core rays + one tangent ray per shell), LTE pure absorption S = B_band(T) (Rayleigh and
e- scattering counted as absorption).  Inner boundary: diffusion I+ = B + mu dB/dtau at
tau_R ~ 1e3.  88 problems (11 bands x 8 g) vectorised.
Temperature: grey Lucy T^4 = Teff^4 (W(r) + 3/4 tau'), tau' = int kR rho (R/r)^2 dr,
then N_LI lambda iterations to radiative equilibrium with the ck opacities
(sum w k (J - B(T)) = 0) in shells with tau_R < TAU_RE; rho(r) held fixed.
usage: python rt.py STAR CASE [MODE]   e.g. python rt.py golden16 hse1 clamp
CASE: hse1 hse3 hse10 | w6v10 w6v30 w5v10 w5v30 | marcs
"""
import sys
import numpy as np
import rsglib as L

import os
N_LI = int(os.environ.get('N_LI', 14))
TAG = os.environ.get('TAG', '')
TAU_RE = 3.0


def hse(st, s, mode, ptop=1e-7, rmax=30.):
    """Grey Lucy hydrostatic atmosphere with scale height x s.  Returns ascending-r arrays."""
    R, M, Te = st['R'], st['M'], st['Teff']
    lP = np.arange(np.log10(ptop), 6.0, 0.025)
    P = 10**lP
    rtop = R * 1.01
    for it in range(4):
        n = len(P)
        r = np.zeros(n)
        T = np.zeros(n)
        rho = np.zeros(n)
        tau = np.zeros(n)
        kR = np.zeros(n)

        def tstate(Pk, rk, tk):
            W = 0.5*(1 - np.sqrt(max(1 - (R/rk)**2, 0.))) if rk > R else 0.5
            Tk = Te*(W + 0.75*tk)**0.25
            mu = L.ck.comp(L.CE, Tk, Pk*1e-6)[0]
            rk_ = Pk*mu*L.MH/(L.KB*Tk)
            s_ = L.state(Tk, rk_, mode)
            k_ = L.means(Tk, s_['K'])[0]
            return Tk, rk_, k_
        r[0] = rtop
        T[0], rho[0], kR[0] = tstate(P[0], r[0], 0.)
        tau[0] = kR[0]*P[0]*s/(L.G*M/r[0]**2)
        T[0], rho[0], kR[0] = tstate(P[0], r[0], tau[0])
        nend = n
        for k in range(1, n):
            dP = P[k] - P[k-1]
            rk, tk, rhok, kk = r[k-1], tau[k-1], rho[k-1], kR[k-1]
            for _ in range(3):
                rm = 0.5*(r[k-1] + rk)
                g = L.G*M/rm**2/s
                rhom = np.sqrt(rho[k-1]*rhok)
                dr = dP/(rhom*g)
                rk = r[k-1] - dr
                km = np.sqrt(kR[k-1]*kk)
                tk = tau[k-1] + km*rhom*dr*(R/rm)**2
                Tk, rhok, kk = tstate(P[k], rk, tk)
            r[k], tau[k], T[k], rho[k], kR[k] = rk, tk, Tk, rhok, kk
            if tk > 100. or Tk > 9800.:
                nend = k + 1
                break
        r, T, rho, tau, kR, P_ = r[:nend], T[:nend], rho[:nend], tau[:nend], kR[:nend], \
            P[:nend]
        r23 = np.interp(np.log(2/3.), np.log(np.maximum(tau, 1e-30)), r)
        shift = R - r23
        rtop = min(rtop + shift, rmax*R)
    # drop layers above rmax (only for very extended cases) -- keep monotone
    o = np.argsort(r)
    return dict(r=r[o], T=T[o], rho=rho[o], tau=tau[o], kR=kR[o], P=P_[o])


def mesa_interior(st, mode, dlP=0.025):
    """MESA (r, rho, T) resampled in log P (dlogP = 0.025, as hse) from the surface
    inward to ck tau_R > 100 or T > 9800 K.  Returns ascending-r arrays."""
    m = L.read_mesa(st['name'])
    lP = np.log10(m['P'])                       # increasing inward
    lPg = np.arange(lP[0], lP[-1], dlP)
    lr = np.interp(lPg, lP, np.log(m['r']))
    lrho = np.interp(lPg, lP, np.log(m['rho']))
    lT = np.interp(lPg, lP, np.log(m['T']))
    r, rho, T = np.exp(lr), np.exp(lrho), np.exp(lT)
    tau = 0.
    kprev = None
    nend = len(r)
    for k in range(len(r)):
        kk = L.means(T[k], L.state(T[k], rho[k], mode)['K'])[0]
        if kprev is not None:
            tau += np.sqrt(kk*kprev)*np.sqrt(rho[k]*rho[k-1])*(r[k-1] - r[k])
        kprev = kk
        if tau > 100. or T[k] > 9800.:
            nend = k + 1
            break
    r, rho, T = r[:nend][::-1], rho[:nend][::-1], T[:nend][::-1]
    return dict(r=r, rho=rho, T=T, m=m, tau_in=tau)


def lucy_T(st, r, rho, kR):
    """Recompute grey Lucy T on an arbitrary ascending-r structure."""
    R, Te = st['R'], st['Teff']
    dt = np.zeros_like(r)
    chi = kR*rho*(R/r)**2
    dt[:-1] = 0.5*(chi[1:] + chi[:-1])*np.diff(r)
    tau = np.cumsum(dt[::-1])[::-1] - dt[-1]
    W = np.where(r > R, 0.5*(1 - np.sqrt(np.maximum(1 - (R/r)**2, 0))), 0.5)
    return Te*(W + 0.75*tau)**0.25, tau


def build(star, case, mode):
    st = L.star(star)
    R = st['R']
    if case == 'marcs':
        fn = {'golden16': 'dl_marcs/s4000_g+0.0_m15._t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o'
              '+0.00_r+0.00_s+0.00.mod',
              'betelgeuse': 'dl_marcs/s3600_g-0.5_m15._t05_st_z+0.00_a+0.00_c+0.00_n+0.00'
              '_o+0.00_r+0.00_s+0.00.mod'}[star]
        m = L.read_marcs(fn)
        r = R - m['depth']
        o = np.argsort(r)
        return st, dict(r=r[o], T=m['T'][o], rho=m['rho'][o], fixT=True,
                        src=fn.split('/')[-1][:22])
    if case.startswith('hse'):
        s = float(case[3:])
        a = hse(st, s, mode)
        return st, dict(r=a['r'], T=a['T'], rho=a['rho'], fixT=False)
    if case.startswith('mft'):
        # MESA interior (r, rho, T) + Fuller & Tsuna (2024) eq. 7 chromosphere to R_d = 5 R,
        # dust wind beyond, to 10 R.  case 'mft' = v_con the FT way, 'mft12' = v_esc/v_con
        # fixed 12.5.  MESA T kept for r <= R; Lucy grey T only above R.
        a = mesa_interior(st, mode)
        m = a['m']
        Rd = 5.*R
        vesc = np.sqrt(2*L.G*st['M']/R)
        vft, rft, tft, info = L.ft_vcon(m)
        vcon = vft if case == 'mft' else vesc/12.5
        ratio = vesc/vcon
        vinf = 30e5
        rho_ph = m['rho'][0]
        mdot = 4*np.pi*R**2*rho_ph*vcon*np.exp(-ratio*np.sqrt(1 - R/Rd))
        rw = R * np.logspace(np.log10(1.0), 1., 140)
        rc = R * (1. + np.logspace(-4., np.log10(0.6), 160))
        r = np.unique(np.concatenate([a['r'], rw[1:], rc]))
        x = np.maximum(r/R, 1.)
        rho_ft = rho_ph*x**-2*np.exp(-ratio*np.sqrt(1 - 1/x))
        vout = np.sqrt(vcon**2 + (vinf**2 - vcon**2)*np.maximum(1 - Rd/r, 0.))
        rho_dw = mdot/(4*np.pi*r**2*vout)
        inside = r <= R*(1 + 1e-12)
        rho = np.where(inside, np.exp(np.interp(r, a['r'], np.log(a['rho']))),
                       np.where(r <= Rd, rho_ft, rho_dw))
        keep = r <= 10.*R*1.0001
        r, rho, inside = r[keep], rho[keep], inside[keep]
        Tm = np.exp(np.interp(r, a['r'], np.log(a['T'])))
        T = Tm.copy()
        for _ in range(3):
            kR = np.array([L.means(T[i], L.state(T[i], rho[i], mode)['K'])[0]
                           for i in range(len(r))])
            Tl, tau = lucy_T(st, r, rho, kR)
            T = np.where(inside, Tm, Tl)
        return st, dict(r=r, T=T, rho=rho, fixT=False, mdot=mdot, v=vcon, rho_ph=rho_ph,
                        vesc=vesc, v_ft=vft, r_ft=rft/R, tau_ft=tft, ft_ok=info['ok'],
                        cv_max_outer=m['cv'][m['r'] > 0.9*R].max(), ratio=ratio)
    if case.startswith('ft'):
        # Fuller & Tsuna (2024) eq. 7 extended chromosphere on top of hse1, dust wind
        # beyond R_d = 5 R (eq. 15 Mdot), truncated at 10 R
        a = hse(st, 1., mode)
        ratio = {'ft8': 8., 'ft10': 10., 'ft12': 12.5, 'ft15': 15.}[case]
        Rd = 5.*R
        vesc = np.sqrt(2*L.G*st['M']/R)
        vcon = vesc/ratio
        vinf = 30e5
        rho_ph = np.exp(np.interp(R, a['r'], np.log(a['rho'])))
        mdot = 4*np.pi*R**2*rho_ph*vcon*np.exp(-ratio*np.sqrt(1 - R/Rd))
        rw = R * np.logspace(np.log10(1.0), 1., 140)
        rc = R * (1. + np.logspace(-4., np.log10(0.6), 160))
        r = np.unique(np.concatenate([a['r'], rw, rc]))
        lrh = np.interp(r, a['r'], np.log(a['rho']), right=-200.)
        rho_h = np.exp(lrh)
        x = np.maximum(r/R, 1.)
        rho_ft = rho_ph*x**-2*np.exp(-ratio*np.sqrt(1 - 1/x))
        vout = np.sqrt(vcon**2 + (vinf**2 - vcon**2)*np.maximum(1 - Rd/r, 0.))
        rho_dw = mdot/(4*np.pi*r**2*vout)
        rho = np.where(r < R, rho_h,
                       np.where(r <= Rd, np.maximum(rho_h, rho_ft), np.maximum(rho_h, rho_dw)))
        keep = r <= 10.*R*1.0001
        r, rho = r[keep], rho[keep]
        T = np.interp(r, a['r'], a['T'])
        for _ in range(3):
            kR = np.array([L.means(T[i], L.state(T[i], rho[i], mode)['K'])[0]
                           for i in range(len(r))])
            T, tau = lucy_T(st, r, rho, kR)
        return st, dict(r=r, T=T, rho=rho, fixT=False, mdot=mdot, v=vcon, rho_ph=rho_ph,
                        vesc=vesc)
    # wind on top of hse1
    a = hse(st, 1., mode)
    mdot = {'6': 1e-6, '5': 1e-5}[case[1]] * L.MSUN / L.YR
    v = float(case.split('v')[1]) * 1e5
    rw = R * np.logspace(np.log10(1.0), 1., 140)
    r = np.unique(np.concatenate([a['r'], rw]))
    lrh = np.interp(r, a['r'], np.log(a['rho']), right=-200.)
    rho = np.maximum(np.exp(lrh), mdot/(4*np.pi*r**2*v))
    r = r[r <= 10.*R*1.0001]
    rho = rho[:len(r)]
    # grey T on the combined structure, kR from the hse T first, then re-evaluated
    T = np.interp(r, a['r'], a['T'])
    for _ in range(3):
        kR = np.array([L.means(T[i], L.state(T[i], rho[i], mode)['K'])[0]
                       for i in range(len(r))])
        T, tau = lucy_T(st, r, rho, kR)
    return st, dict(r=r, T=T, rho=rho, fixT=False, mdot=mdot, v=v)


def opac(T, rho, mode):
    K = np.zeros((len(T), L.NB, L.NG))
    kR = np.zeros(len(T))
    kP = np.zeros(len(T))
    for i in range(len(T)):
        K[i] = L.state(T[i], rho[i], mode)['K']
        kR[i], kP[i] = L.means(T[i], K[i])
    return K, kR, kP


def Bband(T):
    return L.SIG*T[:, None]**4*L.fband(T)/np.pi    # (N, nb), band-integrated intensity


def formal(r, chi, S, nc=16):
    """chi, S: (N, nq).  Returns J, H (N, nq) band-integrated (S units)."""
    N, nq = chi.shape
    mu0 = np.cos(np.linspace(0., 0.5*np.pi*0.97, nc))
    pc = r[0]*np.sqrt(1 - mu0**2)
    p = np.concatenate([pc, r])            # rays
    kray = np.concatenate([np.full(nc, -1), np.arange(N)])  # tangent shell, -1 = core
    nr = len(p)
    z = np.sqrt(np.maximum(r[:, None]**2 - p[None, :]**2, 0.))   # (N, nr)
    valid = r[:, None] >= p[None, :]*(1 - 1e-12)

    def step(Iold, i0, i1, sel):
        dz = np.abs(z[i1, sel] - z[i0, sel])[None, :]
        c0 = chi[i0][:, None]
        c1 = chi[i1][:, None]
        lm = np.where(np.abs(c0 - c1) > 1e-8*c0, (c0 - c1)/np.log(c0/c1), c0)
        dt = np.maximum(dz*lm, 1e-30)
        e0 = -np.expm1(-dt)
        e1 = dt - e0
        wn = np.where(dt > 1e-4, e1/dt, 0.5*dt)
        wo = np.where(dt > 1e-4, e0 - e1/dt, 0.5*dt)
        return Iold*np.exp(-dt) + wo*S[i0][:, None] + wn*S[i1][:, None]

    Im = np.zeros((N, nq, nr))
    Ip = np.zeros((N, nq, nr))
    I = np.zeros((nq, nr))
    for i in range(N-2, -1, -1):
        sel = valid[i]
        I[:, sel] = step(I[:, sel], i+1, i, sel)
        Im[i] = I
    # outward
    I = np.zeros((nq, nr))
    core = kray < 0
    dtr = np.maximum(chi[0]*(r[1] - r[0]), 1e-30)
    dBdt = (S[0] - S[1])/dtr
    I[:, core] = S[0][:, None] + mu0[None, :]*dBdt[:, None]
    t0 = kray == 0
    I[:, t0] = Im[0][:, t0]
    Ip[0] = I
    for i in range(1, N):
        tg = kray == i
        I[:, tg] = Im[i][:, tg]
        sel = valid[i] & ~tg
        I[:, sel] = step(I[:, sel], i-1, i, sel)
        Ip[i] = I
    J = np.zeros((N, nq))
    H = np.zeros((N, nq))
    for i in range(N):
        sel = valid[i]
        mu = z[i, sel]/r[i]
        o = np.argsort(mu)
        mu = mu[o]
        a = Ip[i][:, sel][:, o]
        b = Im[i][:, sel][:, o]
        J[i] = 0.5*np.trapezoid(a + b, mu, axis=1)
        H[i] = 0.5*np.trapezoid((a - b)*mu, mu, axis=1)
    return J, H


def run(star, case, mode='clamp'):
    st, a = build(star, case, mode)
    r, T, rho = a['r'], a['T'].copy(), a['rho']
    N = len(r)
    T0 = T.copy()
    hist = []
    nli = 0 if a['fixT'] else N_LI
    for it in range(nli + 1):
        K, kR, kP = opac(T, rho, mode)
        chi = (K*rho[:, None, None]).reshape(N, -1)
        S = np.repeat(Bband(T), L.NG, axis=1)
        J, H = formal(r, chi, S)
        # tau_R from the top
        dtr = np.zeros(N)
        dtr[:-1] = 0.5*(kR[1:]*rho[1:] + kR[:-1]*rho[:-1])*np.diff(r)
        tauR = np.cumsum(dtr[::-1])[::-1]
        if it == nli:
            break
        w = np.tile(L.GW, L.NB)
        kq = K.reshape(N, -1)
        Tn = T.copy()
        for i in np.where(tauR < TAU_RE)[0]:
            rhs = (w*kq[i]*J[i]).sum()
            kb = (w*kq[i]).reshape(L.NB, L.NG).sum(1)
            lo, hi = 300., 12000.
            for _ in range(40):
                tm = 0.5*(lo + hi)
                f = (kb*Bband(np.array([tm]))[0]).sum() - rhs
                lo, hi = (lo, tm) if f > 0 else (tm, hi)
            Tn[i] = 0.5*(lo + hi)
        Tn = np.sqrt(Tn*T)   # geometric damping
        dT = np.max(np.abs(Tn - T)[tauR < TAU_RE]/T[tauR < TAU_RE])
        hist.append(dT)
        T = Tn
    w = np.tile(L.GW, L.NB)
    kq = K.reshape(N, -1)
    Fq = 4*np.pi*H
    F = (w*Fq).sum(1)
    grad = (w*kq*Fq).sum(1)/L.C
    kF = (w*kq*Fq).sum(1)/F
    Fb = (w*Fq).reshape(N, L.NB, L.NG).sum(2)
    Fstar = st['L']/(4*np.pi*r**2)
    res = dict(r=r/st['R'], rho=rho, T=T, T0=T0, tauR=tauR, kR=kR, kP=kP, kF=kF,
               GF=kF/st['kE'], Gloc=grad/(L.G*st['M']/r**2), Fratio=F/Fstar,
               Fb_top=Fb[-1]/F[-1], fb_bb=L.fband(st['Teff'])[0], hist=np.array(hist),
               kFthin=np.array([L.kappa_thin_flux(K[i], st['Teff']) for i in range(N)]),
               kE=st['kE'], Teff=st['Teff'],
               **{k: a[k] for k in ('mdot', 'v', 'rho_ph', 'vesc', 'v_ft', 'r_ft', 'tau_ft',
                     'ft_ok', 'cv_max_outer', 'ratio') if k in a})
    np.savez(f'col_{star}_{case}_{mode}{TAG}.npz', **res)
    return res


if __name__ == '__main__':
    star, case = sys.argv[1], sys.argv[2]
    mode = sys.argv[3] if len(sys.argv) > 3 else 'clamp'
    res = run(star, case, mode)
    r = res
    i = np.argmax(np.where(r['tauR'] < 1., r['GF'], 0))
    print(f'{star} {case} {mode}: N={len(r["r"])} LI dT/T hist {np.round(r["hist"], 3)}')
    print(f'  F/F* at tauR~1: {np.interp(0., -np.log(r["tauR"]+1e-30), r["Fratio"]):.3f}'
          f'  top: {r["Fratio"][-1]:.3f}')
    print(f'  max Gamma_F (tauR<1) = {r["GF"][i]:.3f} at r/R*={r["r"][i]:.4f} '
          f'T={r["T"][i]:.0f} rho={r["rho"][i]:.2e} tauR={r["tauR"][i]:.2e} '
          f'kF/kR={r["kF"][i]/r["kR"][i]:.3g}  thin GF there {r["kFthin"][i]/r["kE"]:.3f}')
