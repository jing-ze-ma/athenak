#!/usr/bin/env python3
"""1-D spherical radiation-HYDRO toy: Athena++ implicit DO radiation (as toy.py, now with the
flow velocity: explicit advection flux adv_vel*I, HLL with (mu c - adv), comoving source with
tran_coef, lab-frame billing of E and F to the gas exactly as GetHydroSourceTerms/AddSourceTerms)
coupled to a vl2 + PLM + HLLC 1-D spherical hydro with a point mass (Athena++ structure:
stage 1 donor/0.5 dt + rad solve (0.5 dt from I^n); stage 2 PLM/dt + rad solve (dt from I^n)).
Inner BC (hydro): 'his' = bd_data rho,P ghosts + HLLC zero-contact ghost v (iterated on the PLM
face state), or 'reflect'.  Inner BC (rad): 'bath' (his I+-, + ComToLab with ghost v unless
--noboost), 'bathP1', 'mirror' (bc_swap 3).  Outer: his outflow / vacuum.
Every switch of the mechanism hunt is a command-line option (see main()).
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, argparse, json, time
from toy import (opac, taufac, RSUN, RHO_U, T_U, CRAT, PRAT, T_UNIT, DAY, F_U, GAMMA, PROF, ARAD)

GM_CODE = 2.84560241e+03
G1 = GAMMA - 1.0


def vanleer(dl, dr):
    p = dl * dr
    return np.where(p > 0, 2 * p / np.where(p > 0, dl + dr, 1.0), 0.0)


def hllc(rl, ul, pl, rr, ur, pr):
    """Athena++-style HLLC (pressure-based wave speeds). returns mass, mom, energy flux, S*"""
    cl = np.sqrt(GAMMA * pl / rl); cr = np.sqrt(GAMMA * pr / rr)
    rhoa = 0.5 * (rl + rr); ca = 0.5 * (cl + cr)
    pmid = np.maximum(0.5 * (pl + pr + (ul - ur) * rhoa * ca), 0.0)
    ql = np.where(pmid <= pl, 1.0, np.sqrt(1 + (GAMMA + 1) / (2 * GAMMA) * (pmid / pl - 1)))
    qr = np.where(pmid <= pr, 1.0, np.sqrt(1 + (GAMMA + 1) / (2 * GAMMA) * (pmid / pr - 1)))
    sl = ul - cl * ql; sr = ur + cr * qr
    sm = (pr - pl + rl * ul * (sl - ul) - rr * ur * (sr - ur)) / (rl * (sl - ul) - rr * (sr - ur))
    El = pl / G1 + 0.5 * rl * ul**2; Er = pr / G1 + 0.5 * rr * ur**2
    Fl = np.array([rl * ul, rl * ul**2 + pl, (El + pl) * ul]); Fr = np.array([rr * ur, rr * ur**2 + pr, (Er + pr) * ur])
    Ul = np.array([rl, rl * ul, El]); Ur = np.array([rr, rr * ur, Er])
    def star(U, s, u, r, p):
        fac = r * (s - u) / (s - sm)
        E = U[2]
        return np.array([fac, fac * sm, fac * (E / r + (sm - u) * (sm + p / (r * (s - u))))])
    Usl = star(Ul, sl, ul, rl, pl); Usr = star(Ur, sr, ur, rr, pr)
    F = np.where(sl >= 0, Fl, np.where(sm >= 0, Fl + sl * (Usl - Ul), np.where(sr > 0, Fr + sr * (Usr - Ur), Fr)))
    return F, sm


class RHD:
    def __init__(s, a):
        s.a = a
        N = a.nx; s.N = N
        s.rf = np.linspace(a.rmin, a.rmax, N + 1); s.dx = np.diff(s.rf)
        r0, r1 = s.rf[:-1], s.rf[1:]
        s.rv = 0.75 * (r1**4 - r0**4) / (r1**3 - r0**3)
        s.V = (r1**3 - r0**3) / 3.0; s.A = s.rf**2
        s.csrc = 1.5 * (r1**2 - r0**2) / (r1**3 - r0**3)        # coord_src1_i
        dx0 = s.dx[0]
        gf = a.rmin - dx0 * np.arange(3)                            # ghost faces
        s.rg = np.array([0.75 * (gf[i]**4 - gf[i + 1]**4) / (gf[i]**3 - gf[i + 1]**3) for i in range(2)])
        x, w = np.polynomial.legendre.leggauss(2 * a.nmu)
        s.mu = x; s.w = 0.5 * w; s.NA = len(x); s.pos = s.mu > 0
        s.cw = s.w[s.pos].sum(); s.cmw = (s.mu * s.w)[s.pos].sum()
        pr = PROF[:, 0] / RSUN; s.prmax = pr[-1]
        Tp = PROF[:, 1] / T_U; Rp = PROF[:, 2] / RHO_U; Fp = PROF[:, 5] / F_U
        lin = lambda x, y: np.interp(x, pr, y, left=None)
        def interp_ext(x, y):      # bsg_star.cpp extrapolates linearly below the first point; profile starts 16.5
            return np.interp(x, pr, y)
        T = interp_ext(s.rv, Tp); rho = interp_ext(s.rv, Rp); F = interp_ext(s.rv, Fp)
        out = s.rv > s.prmax
        if out.any():
            G = 6.674e-8 * RHO_U * T_UNIT**2
            Mmax = PROF[-1, 3] / (RHO_U * RSUN**3)
            gmax = G * Mmax / s.prmax**2
            gl = gmax * (s.rv[out] / s.prmax)**2
            T[out] = Tp[-1]
            rho[out] = np.maximum(Rp[-1] * np.exp(-gl * (s.rv[out] - s.prmax) / Tp[-1]), a.dfloor)
            F[out] = 0.0
        s.Tg = interp_ext(s.rg, Tp); s.rhog = interp_ext(s.rg, Rp); s.Fg = interp_ext(s.rg, Fp)
        s.Lref = s.Fg[0] * s.rg[0]**2
        s.rho = rho.copy(); P = rho * T
        if a.noise > 0:
            rng = np.random.default_rng(a.seed)
            amp = np.where(rho > 0.01, a.noise, 0.0)
            s.rho = rho * (1 + amp * (rng.random(N) - 0.5))
        s.mom = np.zeros(N)
        s.Eg = P / G1
        s.I = s.ic_intensity(T**4, F)
        s.vg = np.zeros(2)        # inner hydro ghost velocities
        s.t = 0.0; s.cycle = 0
        s.E0 = s.etot()
        s.led = dict(bin=0.0, bout=0.0, hin=0.0, hout=0.0, reset=0.0)

    # ---------------------------------------------------------------- helpers
    def ic_intensity(s, E, F):
        E = np.atleast_1d(E); F = np.atleast_1d(F)
        I = np.empty((len(E), s.NA))
        I[:, s.pos] = (0.5 * (E / s.cw + F / s.cmw))[:, None]
        I[:, ~s.pos] = (0.5 * (E / s.cw - F / s.cmw))[:, None]
        return I

    def prim(s, rho, mom, Eg):
        v = mom / rho
        P = np.maximum((Eg - 0.5 * rho * v**2) * G1, s.a.pfloor)
        return rho, v, P

    def Tgas(s):
        r, v, P = s.prim(s.rho, s.mom, s.Eg)
        return np.maximum(P / r, s.a.tfloor)

    def etot(s):
        return np.sum(s.V * (s.Eg + PRAT * (s.I @ s.w) - s.rho * GM_CODE / s.rv))

    # ---------------------------------------------------------------- hydro
    def ghosts(s, rho, v, P):
        """returns extended arrays with 2 ghosts each side"""
        a = s.a
        if a.hbc == 'reflect':
            rgh = rho[[1, 0]]; vgh = -v[[1, 0]]; pgh = P[[1, 0]]
        else:
            rgh = s.rhog[::-1].copy(); pgh = (s.rhog * s.Tg)[::-1].copy(); vgh = np.array([0.0, s.vg[0]])
        # outer (his Outflow_X1)
        x1 = s.rv[-1]
        ro = np.empty(2); vo = np.empty(2); po = np.empty(2)
        rprev = rho[-1]; rprev_r = x1
        for i in range(2):
            x1g = s.rv[-1] + (i + 1) * s.dx[-1]
            if v[-1] < 0:
                ro[i] = rho[-1]; vo[i] = 0.0
            else:
                ro[i] = rprev * rprev_r**2 / x1g**2; vo[i] = v[-1]
            rprev = ro[i]; rprev_r = x1g
            po[i] = P[-1]
        R = np.concatenate([rgh, rho, ro]); Vv = np.concatenate([vgh, v, vo]); Pp = np.concatenate([pgh, P, po])
        return R, Vv, Pp

    def recon(s, q, order):
        if order == 1:
            return q[1:-2], q[2:-1]              # faces 0..N: left state from cell k-1, right from k
        dq = vanleer(q[1:-1] - q[:-2], q[2:] - q[1:-1])   # slopes for cells -1..N
        qc = q[1:-1]
        ql = qc[:-1] + 0.5 * dq[:-1]            # left of face k (cell k-1)
        qr = qc[1:] - 0.5 * dq[1:]
        return ql, qr

    def hydro_rhs(s, rho, mom, Eg, order):
        r, v, P = s.prim(rho, mom, Eg)
        if s.a.hbc == 'his':
            s.zero_contact(r, v, P, order)
        R, Vv, Pp = s.ghosts(r, v, P)
        rl, rr = s.recon(R, order); ul, ur = s.recon(Vv, order); pl, pr = s.recon(Pp, order)
        rl = np.maximum(rl, s.a.dfloor); rr = np.maximum(rr, s.a.dfloor)
        pl = np.maximum(pl, s.a.pfloor); pr = np.maximum(pr, s.a.pfloor)
        F, sm = hllc(rl, ul, pl, rr, ur, pr)
        s.mflux = F[0]
        dF = (s.A[1:] * F[:, 1:] - s.A[:-1] * F[:, :-1]) / s.V
        drho = -dF[0]; dmom = -dF[1] + P * 2 * s.csrc; dE = -dF[2]
        g = GM_CODE * s.csrc / s.rv
        dmom -= r * g
        dE -= 0.5 * (F[0, :-1] + F[0, 1:]) * g
        s.hface = (F[2, 0] * s.A[0], F[2, -1] * s.A[-1])
        return drho, dmom, dE, (r, v, P)

    def zero_contact(s, r, v, P, order):
        """bsg_star.cpp: iterate the innermost ghost v so the HLLC contact at face 0 is zero"""
        for it in range(100):
            R, Vv, Pp = s.ghosts(r, v, P)
            sl = slice(0, 4)
            rl, rr = s.recon(R[:5], order); ul, ur = s.recon(Vv[:5], order); pl, pr = s.recon(Pp[:5], order)
            rl, rr, ul, ur, pl, pr = rl[0], rr[0], ul[0], ur[0], pl[0], pr[0]
            cl = np.sqrt(GAMMA * pl / rl); cr = np.sqrt(GAMMA * pr / rr)
            rhoa = 0.5 * (rl + rr); ca = 0.5 * (cl + cr)
            pmid = 0.5 * (pl + pr + (ul - ur) * rhoa * ca)
            ql = 1.0 if pmid <= pl else np.sqrt(1 + (GAMMA + 1) / (2 * GAMMA) * (pmid / pl - 1))
            qr = 1.0 if pmid <= pr else np.sqrt(1 + (GAMMA + 1) / (2 * GAMMA) * (pmid / pr - 1))
            des = (pr - pl - rr * cr * qr * ur) / (rl * cl * ql)
            num = pl - pr + rl * cl * ql * ul + rr * cr * qr * ur
            if abs(num) <= 1e-14 * max(pl, pr):
                break
            s.vg[0] += des - ul

    # ---------------------------------------------------------------- radiation
    def rad_solve(s, wght, I_ini, I_cur, rho, v, Tg0, order):
        a = s.a; N, NA = s.N, s.NA
        ss, sa, spp = opac(rho, Tg0)
        if a.exch == 'ross':
            spp = sa.copy()
        chi = ss + sa
        if a.hbc == 'reflect':
            rgh, Tgh = rho[0], Tg0[0]
        else:
            rgh, Tgh = s.rhog[0], s.Tg[0]
        gs, ga, _ = opac(np.array([rgh]), np.array([Tgh])); chig = (gs + ga)[0]
        chiL = np.concatenate([[chig], chi]); chiR = np.concatenate([chi, [chi[-1]]])
        dxL = np.concatenate([[s.dx[0]], s.dx]); dxR = np.concatenate([s.dx, [s.dx[-1]]])
        tau = a.taucell * (dxL * chiL + dxR * chiR)
        fr, fl = taufac(tau)
        # advection velocity at faces (GetTaufactorAdv), ghost v included
        if a.hbc == 'reflect':
            vgi = -v[0]
        else:
            vgi = s.vg[0]
        vo = v[-1] if v[-1] >= 0 else 0.0
        vext = np.concatenate([[vgi], v, [vo]])
        t2 = tau * tau
        fadv = np.where(t2 > 1e-3, 1 - np.exp(-t2), t2 - 0.5 * t2 * t2)
        if a.noadv:
            fadv[:] = 0.0
        if a.wall_noadv:
            fadv[0] = 0.0
        adv = fadv * 0.5 * (vext[:-1] + vext[1:])
        # inner ghost intensities
        if a.bc == 'bath' or a.bc == 'bathP1':
            if a.bc == 'bath':
                Ig = s.ic_intensity(s.Tg[0]**4, s.Fg[0])[0]
            else:
                Ig = s.Tg[0]**4 + 3 * s.mu * s.Fg[0]
            if not a.noboost:      # ComToLab with the ghost velocity
                vv = vgi
                lor = 1 / np.sqrt(1 - (vv / CRAT)**2)
                tc = lor * (1 - s.mu * vv / CRAT)
                Ig = Ig / tc**4
            Igc = Ig; mirror = False
        else:
            mirror = True
            Igc = 6 * s.mu * s.Fg[0]
        # explicit advection flux of I_cur
        Iext = np.vstack([(I_cur[0, ::-1] + Igc)[None] if mirror else Ig[None], I_cur,
                          np.where(s.pos, I_cur[-1], 0.0)[None]])
        if order == 1:
            IL, IR = Iext[:-1], Iext[1:]
        else:
            Ie2 = np.vstack([Iext[:1], Iext, Iext[-1:]])
            dq = vanleer(Ie2[1:-1] - Ie2[:-2], Ie2[2:] - Ie2[1:-1])
            IL = Iext[:-1] + 0.5 * dq[:-1]; IR = Iext[1:] - 0.5 * dq[1:]
        afl = np.where(adv[:, None] > 0, adv[:, None] * IL, adv[:, None] * IR)
        advdiv = -wght * (s.A[1:, None] * afl[1:] - s.A[:-1, None] * afl[:-1]) / s.V[:, None]
        # implicit HLL coefficients
        vel = CRAT * s.mu[None, :]
        frk, flk = fr[:, None], fl[:, None]
        smax = np.where(vel > 0, frk * vel, -flk * vel); smin = np.where(vel > 0, -flk * vel, frk * vel)
        if a.flux == 'upwind':
            smax = np.where(vel > 0, vel, 0.0 * vel); smin = np.where(vel > 0, 0.0 * vel, vel)
        u = vel - adv[:, None]
        dd = smax - smin
        Ac = smax * (u - smin) / dd; Bc = smin * (smax - u) / dd
        if a.flux == 'central':
            Ac = 0.5 * u * np.ones_like(Ac); Bc = 0.5 * u * np.ones_like(Bc)
        if a.flux == 'compact':          # HLL + extra compact diffusion boosting eps
            pass
        # Doppler factors
        lor = 1 / np.sqrt(1 - (v / CRAT)**2)
        tc = lor[:, None] * (1 - s.mu[None, :] * v[:, None] / CRAT)       # N x NA
        if a.nodoppler:
            tc = np.ones_like(tc)
        if a.nodop_wall > 0:
            tc[:a.nodop_wall] = 1.0
        if a.nodop_from > 0:
            tc[a.nodop_from:] = 1.0
        wcm = s.w[None, :] / tc**2; wcm = wcm / wcm.sum(1, keepdims=True)
        cdt = CRAT * wght
        V = s.V; Af = s.A; dt = wght
        nu = NA + 1; ntot = N * nu
        ii = np.arange(N)
        iI = ii[:, None] * nu + np.arange(NA)[None, :]; iT = ii * nu + NA
        # transport matrix (per unit "V*I")
        rows = []; cols = []; vals = []
        for n in range(NA):
            diag = dt * (Af[1:] * Ac[1:, n] - Af[:-1] * Bc[:-1, n])
            if s.mu[n] > 0:
                diag[-1] += dt * Af[-1] * Bc[-1, n]
            rows += [iI[:, n], iI[:-1, n], iI[1:, n]]
            cols += [iI[:, n], iI[1:, n], iI[:-1, n]]
            vals += [diag, dt * Af[1:-1] * Bc[1:-1, n], -dt * Af[1:-1] * Ac[1:-1, n]]
            if mirror:
                rows.append(np.array([iI[0, n]])); cols.append(np.array([iI[0, NA - 1 - n]]))
                vals.append(np.array([-dt * Af[0] * Ac[0, n]]))
        Ltr = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(ntot, ntot))
        bI = np.zeros((N, NA)); bI[0] = -dt * Af[0] * Ac[0] * Igc
        I = I_cur.copy(); T = Tg0.copy()
        clipmask = None
        for it in range(a.newton):
            X = np.zeros(ntot); X[iI] = I; X[iT] = T
            LX = (Ltr @ X)[iI] + bI                         # V * transport(I)
            Jcm = np.sum(wcm * tc**4 * I, 1)
            T4 = T**4
            src = (chi[:, None] - spp[:, None]) * Jcm[:, None] + spp[:, None] * T4[:, None]
            RI = V[:, None] * ((I - I_ini - advdiv) + cdt * (chi[:, None] * tc * I - src / tc**3)) + LX
            if a.clip:
                # Athena++ clips I* = I_ini + divflx at TINY before the source (per angle)
                Istar = I_ini + advdiv - LX / V[:, None] + 0.0
                # NOTE: the split own-cell part makes I* differ from this; see --clip doc in REPORT
                clipmask = Istar < 1e-20
            RT = V * (rho * (T - Tg0) / G1 + PRAT * cdt * spp * (T4 - Jcm))
            R = np.zeros(ntot); R[iI] = RI; R[iT] = RT
            # local Jacobian
            dJ = (wcm * tc**4)                                # dJcm/dI_m  (N x NA)
            blk = np.zeros((N, nu, nu))
            for n in range(NA):
                blk[:, n, :NA] = -(V * cdt)[:, None] * ((chi - spp)[:, None] / tc[:, n:n + 1]**3) * dJ
                blk[:, n, n] += V * (1 + cdt * chi * tc[:, n])
                blk[:, n, NA] = -V * cdt * spp * 4 * T**3 / tc[:, n]**3
            blk[:, NA, :NA] = -(V * PRAT * cdt * spp)[:, None] * dJ
            blk[:, NA, NA] = V * (rho / G1 + PRAT * cdt * spp * 4 * T**3)
            br = (ii[:, None, None] * nu + np.arange(nu)[None, :, None]) * np.ones((1, 1, nu), int)
            bc = (ii[:, None, None] * nu + np.arange(nu)[None, None, :]) * np.ones((1, nu, 1), int)
            Jl = sp.csr_matrix((blk.ravel(), (br.ravel(), bc.ravel())), shape=(ntot, ntot))
            Jm = (Ltr + Jl).tocsr()
            d = np.abs(Jm.diagonal()); Dr = sp.diags(1.0 / d)
            dX = spl.spsolve((Dr @ Jm).tocsc(), -(Dr @ R))
            dI = dX[iI]; dT = dX[iT]
            lam = 1.0
            bad = (T + dT) < 0.5 * T
            if bad.any():
                lam = min(1.0, np.min(0.5 * T[bad] / np.abs(dT[bad])))
            I = I + lam * dI; T = T + lam * dT
            if lam == 1.0 and np.max(np.abs(dT) / T) < a.tol and np.max(np.abs(dI).sum(1) / np.abs(I).sum(1)) < a.tol:
                break
        s.nit = it + 1
        # billing (GetHydroSourceTerms): delta = I_new - (I_ini + divflx - const*I_new), lab weights
        LX = (Ltr @ np.concatenate([np.zeros(0)]) if False else None)
        X = np.zeros(ntot); X[iI] = I
        trans = ((Ltr @ X)[iI] + bI) / V[:, None]
        dsrc = I - (I_ini + advdiv - trans)
        dE = dsrc @ s.w; dF = dsrc @ (s.w * s.mu)
        # boundary fluxes (implicit + explicit advective), energy per unit time*area
        Igf = (I[0, ::-1] + Igc) if mirror else Igc
        F0 = (Ac[0] * Igf + Bc[0] * I[0] + afl[0]) @ s.w
        Iout = np.where(s.pos, I[-1], 0.0)
        FN = (Ac[-1] * I[-1] + Bc[-1] * Iout + afl[-1]) @ s.w
        s.rad_f0 = F0; s.rad_fN = FN
        s.Fface = (Ac * np.vstack([Igf[None], I]) + Bc * np.vstack([I, Iout[None]])) @ s.w
        s.Fadv = afl @ s.w
        s.chi = chi
        if a.debug:
            Jc = np.sum(wcm * tc**4 * I, 1)
            s.dbg = dict(I_ini=I_ini.copy(), advdiv=advdiv.copy(), trans=trans.copy(), I=I.copy(), tc=tc.copy(),
                         chi=chi.copy(), spp=spp.copy(), Jcm=Jc, T=T.copy(), cdt=cdt, v=v.copy(), Ac=Ac.copy(),
                         Bc=Bc.copy(), adv=adv.copy(), Ig=np.array(Igc), dE=dE.copy(), dF=dF.copy(), rho=rho.copy(), Tg0=Tg0.copy())
        return I, dE, dF

    def vrad(s, v_stage, rho, I_ini, wght, Tg0):
        """velocity handed to the radiation solve (Athena++: the stage velocity = 'stage')"""
        m = s.a.vsrc
        if m == 'stage':
            return v_stage
        if m == 'old':
            return s.v_n
        if m == 'pred':   # stage v + the radiation-force impulse over the solve, from the start-of-step comoving flux
            ss, sa, _ = opac(rho, Tg0)
            E = I_ini @ s.w; F = I_ini @ (s.w * s.mu)
            Fcm = F - 4.0 / 3.0 * (s.v_n / CRAT) * E
            return v_stage + wght * PRAT * (ss + sa) * Fcm / rho
        raise ValueError(m)

    def add_rad_source(s, rho, mom, Eg, dE, dF):
        mom = mom - PRAT * dF / CRAT
        es = -PRAT * dE
        ek = 0.5 * mom**2 / rho
        eint = Eg + es - ek
        bad = eint < 0
        Eg = Eg + es
        if bad.any():
            Tn = s._Tnew_fallback
            Eg = np.where(bad, Tn * rho / G1 + ek, Eg)
        return mom, Eg

    # ---------------------------------------------------------------- full step (vl2)
    def step(s):
        a = s.a
        r, v, P = s.prim(s.rho, s.mom, s.Eg)
        cs = np.sqrt(GAMMA * P / r)
        dt = a.cfl * np.min(s.dx / (np.abs(v) + cs))
        if a.dtfix > 0:
            dt = a.dtfix / T_UNIT
        I_n = s.I.copy(); U_n = (s.rho.copy(), s.mom.copy(), s.Eg.copy())
        s.v_n = v.copy()
        E_before = s.etot()
        # stage 1
        drho, dmom, dEg, _ = s.hydro_rhs(*U_n, order=1)
        rho1 = np.maximum(U_n[0] + 0.5 * dt * drho, a.dfloor); mom1 = U_n[1] + 0.5 * dt * dmom; Eg1 = U_n[2] + 0.5 * dt * dEg
        if not a.frozen_rad:
            r1, v1, P1 = s.prim(rho1, mom1, Eg1)
            T1 = np.maximum(P1 / r1, a.tfloor)
            s._Tnew_fallback = T1
            I1, dE1, dF1 = s.rad_solve(0.5 * dt, I_n, I_n, r1, s.vrad(v1, r1, I_n, 0.5 * dt, T1), T1, order=1)
            mom1, Eg1 = s.add_rad_source(rho1, mom1, Eg1, dE1, dF1)
        else:
            I1 = I_n
        # stage 2
        hf1 = s.hface
        drho, dmom, dEg, _ = s.hydro_rhs(rho1, mom1, Eg1, order=a.xorder)
        hf = s.hface
        rho2 = np.maximum(U_n[0] + dt * drho, a.dfloor); mom2 = U_n[1] + dt * dmom; Eg2 = U_n[2] + dt * dEg
        if a.frozen_hydro:
            rho2, mom2 = U_n[0].copy(), U_n[1].copy()
            Eg2 = U_n[2].copy()
        r2, v2, P2 = s.prim(rho2, mom2, Eg2)
        T2 = np.maximum(P2 / r2, a.tfloor)
        s._Tnew_fallback = T2
        I2, dE2, dF2 = s.rad_solve(dt, I_n, I1, r2, s.vrad(v2, r2, I_n, dt, T2), T2, order=a.rad_xorder)
        mom2, Eg2 = s.add_rad_source(rho2, mom2, Eg2, dE2, dF2)
        s.rho, s.mom, s.Eg, s.I = rho2, mom2, Eg2, I2
        s.t += dt; s.cycle += 1; s.dt = dt
        s.led['bin'] += dt * PRAT * s.rad_f0 * s.A[0]
        s.led['bout'] += dt * PRAT * s.rad_fN * s.A[-1]
        s.led['hin'] += dt * hf[0]
        s.led['hout'] += dt * hf[1]

    def diag(s):
        r, v, P = s.prim(s.rho, s.mom, s.Eg)
        T = P / r
        J = s.I @ s.w; F = s.I @ (s.w * s.mu)
        L = lambda k: (s.Fface[k] / CRAT + s.Fadv[k] / CRAT) * s.rf[k]**2 / s.Lref
        # r(rho = 1e-10 g/cc) = 0.1 code
        k = np.where(r * RHO_U > 1e-10)[0]
        rs = s.rv[k[-1]] if len(k) else np.nan
        return dict(t=s.t / DAY, dt=s.dt * T_UNIT, T0=T[0] * T_U, T1=T[1] * T_U, T2=T[2] * T_U, T3=T[3] * T_U,
                    rho0=r[0], rho1=r[1], rho2=r[2], v0=v[0] * np.sqrt(8.314462618e7 / 0.62 * 1e5) / 1e5,
                    Lw=L(0), L1=L(1), L2=L(2), L5=L(5), L50=L(50), Lout=L(s.N),
                    Fc0=F[0] * s.rv[0]**2 / s.Lref, Fc1=F[1] * s.rv[1]**2 / s.Lref, Fc2=F[2] * s.rv[2]**2 / s.Lref,
                    dE=s.etot() - s.E0, bin=s.led['bin'], bout=s.led['bout'], hin=s.led['hin'], hout=s.led['hout'],
                    rsurf=rs, vg=s.vg[0] * 3.6618e6 / 1e5, nit=s.nit)


ECGS = 4 * np.pi * RSUN**3 * ARAD * T_U**4 / PRAT   # code energy -> erg (full sphere)


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('--nx', type=int, default=384); p.add_argument('--rmin', type=float, default=20.0)
    p.add_argument('--rmax', type=float, default=80.0); p.add_argument('--nmu', type=int, default=4)
    p.add_argument('--taucell', type=float, default=4.0); p.add_argument('--cfl', type=float, default=0.3)
    p.add_argument('--dtfix', type=float, default=0.0, help='fixed dt [s] (0 = cfl)')
    p.add_argument('--days', type=float, default=3.0)
    p.add_argument('--bc', default='bath', help='bath | bathP1 | mirror')
    p.add_argument('--hbc', default='his', help='his | reflect')
    p.add_argument('--noboost', action='store_true'); p.add_argument('--nodoppler', action='store_true')
    p.add_argument('--nodop_wall', type=int, default=0, help='Doppler (tran_coef) off in the first K cells')
    p.add_argument('--nodop_from', type=int, default=0, help='Doppler off from cell K outward')
    p.add_argument('--noadv', action='store_true'); p.add_argument('--wall_noadv', action='store_true')
    p.add_argument('--flux', default='hll', help='hll | upwind | central')
    p.add_argument('--exch', default='planck'); p.add_argument('--clip', action='store_true')
    p.add_argument('--noise', type=float, default=0.05); p.add_argument('--seed', type=int, default=0)
    p.add_argument('--dfloor', type=float, default=1e-11); p.add_argument('--pfloor', type=float, default=1e-15)
    p.add_argument('--tfloor', type=float, default=0.05)
    p.add_argument('--xorder', type=int, default=2); p.add_argument('--rad_xorder', type=int, default=2)
    p.add_argument('--frozen_hydro', action='store_true'); p.add_argument('--frozen_rad', action='store_true')
    p.add_argument('--newton', type=int, default=12); p.add_argument('--tol', type=float, default=3e-7)
    p.add_argument('--out', default='run'); p.add_argument('--every', type=float, default=0.05)
    p.add_argument('--prof', type=float, default=0.25)
    p.add_argument('--debug', action='store_true')
    p.add_argument('--vsrc', default='stage', help='velocity in the rad solve: stage (Athena++) | old (v^n) | pred')
    p.add_argument('--save', default=''); p.add_argument('--load', default='')
    a = p.parse_args(argv)
    s = RHD(a)
    if a.load:
        z = np.load(a.load)
        s.rho, s.mom, s.Eg, s.I, s.vg, s.t = z['rho'], z['mom'], z['Eg'], z['I'], z['vg'], float(z['t'])
        s.E0 = float(z['E0'])
        for k in s.led: s.led[k] = float(z['led_' + k])
    hist = []; profs = []
    t0 = time.time()
    def dump():
        r, v, P = s.prim(s.rho, s.mom, s.Eg)
        profs.append(dict(t=s.t / DAY, rho=r[:80].tolist(), v=v[:80].tolist(), T=(P / r)[:80].tolist(),
                          J=(s.I @ s.w)[:80].tolist(), F=(s.I @ (s.w * s.mu))[:80].tolist(),
                          Ff=(s.Fface[:81] / CRAT).tolist() if hasattr(s, 'Fface') else [],
                          Fa=(s.Fadv[:81] / CRAT).tolist() if hasattr(s, 'Fadv') else [],
                          chi=s.chi[:80].tolist() if hasattr(s, 'chi') else [],
                          rhoall=r.tolist(), I=s.I[:12].tolist()))
    dump()
    nexth = a.every; nextp = a.prof
    tend = s.t / DAY + a.days if a.load else a.days
    while s.t / DAY < tend:
        try:
            s.step()
        except Exception as e:
            print('EXCEPTION', repr(e)[:200]); print('DEAD at', s.t / DAY); break
        ok = np.all(np.isfinite(s.Eg)) and np.all(np.isfinite(s.I)) and s.dt > 1e-6 * 212 / T_UNIT
        if s.t / DAY >= nexth or not ok:
            nexth += a.every
            d = s.diag(); hist.append(d)
            print(f"t {d['t']:.3f} dt {d['dt']:.0f} T0 {d['T0']:.3e} T1 {d['T1']:.3e} T2 {d['T2']:.3e} "
                  f"rho0-2 {d['rho0']:.3g} {d['rho1']:.3g} {d['rho2']:.3g} v0 {d['v0']:+.3f} vg {d['vg']:+.3f} "
                  f"Lw {d['Lw']:+.1f} L1 {d['L1']:+.1f} L5 {d['L5']:+.1f} L50 {d['L50']:+.2f} Lout {d['Lout']:.2f} "
                  f"dE {d['dE']*ECGS:+.3e} bin {d['bin']*ECGS:+.3e} hin {d['hin']*ECGS:+.2e} "
                  f"rs {d['rsurf']:.2f} nit {d['nit']} ({time.time()-t0:.0f}s)", flush=True)
        if s.t / DAY >= nextp:
            nextp += a.prof; dump()
        if not ok:
            print('DEAD at', s.t / DAY); break
    dump()
    if a.save:
        np.savez(a.save, rho=s.rho, mom=s.mom, Eg=s.Eg, I=s.I, vg=s.vg, t=s.t, E0=s.E0,
                 **{'led_' + k: v for k, v in s.led.items()})
    json.dump(dict(args=vars(a), rv=s.rv.tolist(), rf=s.rf.tolist(), hist=hist, profs=profs, Lref=s.Lref,
                   ECGS=ECGS), open(a.out + '.json', 'w'))


if __name__ == '__main__':
    main()
