#!/usr/bin/env python3
"""1-D spherical implicit discrete-ordinates radiation + gas energy toy that copies the
Athena++ (e47fabb7/053bc73) implicit radiation discretisation used for Linhao Ma's BSG:

  * per-angle HLL face flux with the tau-factor f_r/f_l (GetTaufactor, tau_scheme 1),
    tau = taucell * (dx sigma_L + dx sigma_R)   (implicit_first_order.cpp)
  * backward Euler for I over dt (stage-2 of vl2 restarts from ir1 with wght = dt)
  * cell source exactly as AbsorptionScattering (v = 0):
      (I - I0) + dt*div = c dt [ -(sa+ss) I + (ss+sa-sp) J + sp T^4 ]
      rho (T - T0)/(g-1) = prat c dt sp (J - T^4)
  * 1-D angle set: 2*nmu Gauss-Legendre on [-1,1], w = wGL/2 (angulargrid.cpp, angle_flag 0)
  * his inner bath ghost I+- = 0.5 (E/Sum w+ +- F/Sum mu w+), E = T_g^4, F = profile (no boost)
  * outer: copy outgoing, zero incoming.
  * opacities lagged (evaluated at the start of the step), as Athena++.
The solve is a full Newton on all (I, T) unknowns (direct = the converged fixed point).
Optional: 1-D hydro (frozen by default), the per-angle RHS clip, alternative fluxes/BCs.
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, argparse, json, sys, time

# ------------------------------------------------------------------ units (Athena++ radiation.cpp)
RSUN = 6.957e10; RHO_U = 1e-9; T_U = 1e5; MU = 0.62
R_IDEAL = 8.314462618e7 / MU
CS = np.sqrt(R_IDEAL * T_U)
CRAT = 2.99792458e10 / CS
ARAD = 7.565733250033928e-15
PRAT = ARAD * T_U**3 / (RHO_U * R_IDEAL)
T_UNIT = RSUN / CS                      # s
DAY = 86400.0 / T_UNIT                  # code time per day
F_U = ARAD * T_U**4 * 2.99792458e10     # cgs flux per code Fr
GAMMA = 5.0 / 3.0
KES = 0.34

D = np.load('/resnick/groups/carnegie_poc/jingze/athpp_toy_1104/inputs.npz')
PROF = D['prof']; KR = D['kr']; KP = D['kp']; LT0, DLT, LD0, DLD = D['tab']


def opac(rho_code, T_code):
    """bilinear in (log T, log rho) of the padded tables; returns sigma_s, sigma_a, sigma_p (code)."""
    rho = rho_code * RHO_U; T = T_code * T_U
    x = (np.log10(T) - LT0) / DLT; y = (np.log10(rho) - LD0) / DLD
    x = np.clip(x, 0, KR.shape[0] - 1.000001); y = np.clip(y, 0, KR.shape[1] - 1.000001)
    i = x.astype(int); j = y.astype(int); fx = x - i; fy = y - j
    def bil(K):
        return 10**((1-fx)*(1-fy)*K[i, j] + fx*(1-fy)*K[i+1, j] + (1-fx)*fy*K[i, j+1] + fx*fy*K[i+1, j+1])
    kr = bil(KR); kp = bil(KP)
    ss = KES * rho * RSUN
    sa = (kr - KES) * rho * RSUN
    spl_ = kp * rho * RSUN
    return ss, sa, spl_


def taufac(tau):
    """GetTaufactor tau_scheme 1: returns (f_r same direction, f_l opposite)."""
    t2 = tau * tau
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        fr = np.where(t2 > 1e-3, (1 - np.exp(-t2)) / np.maximum(t2, 1e-300), 1 - 0.5 * t2)
        fl = np.where(t2 > 1e-3, (1 - np.exp(-t2 * t2)) / np.maximum(t2, 1e-300), t2)
    return np.sqrt(fr), np.sqrt(fl)


class Toy:
    def __init__(s, a):
        s.a = a
        N = a.nx; s.N = N
        s.rf = np.linspace(a.rmin, a.rmax, N + 1)
        s.dx = np.diff(s.rf)
        r0, r1 = s.rf[:-1], s.rf[1:]
        s.rv = 0.75 * (r1**4 - r0**4) / (r1**3 - r0**3)          # Athena++ spherical x1v
        s.V = (r1**3 - r0**3) / 3.0
        s.A = s.rf**2
        # ghost centre (inner), one ghost needed for the first-order implicit flux
        g0, g1 = a.rmin - s.dx[0], a.rmin
        s.rg = 0.75 * (g1**4 - g0**4) / (g1**3 - g0**3)
        # angles
        x, w = np.polynomial.legendre.leggauss(2 * a.nmu)
        s.mu = x; s.w = 0.5 * w; s.NA = len(x)
        s.pos = s.mu > 0
        s.cw = s.w[s.pos].sum(); s.cmw = (s.mu * s.w)[s.pos].sum()
        # IC (profile interpolated linearly in r, as bsg_star.cpp)
        pr = PROF[:, 0] / RSUN; s.prmax = pr[-1]
        Tp = PROF[:, 1] / T_U; Rp = PROF[:, 2] / RHO_U; Fp = PROF[:, 5] / F_U
        s.prof = (pr, Tp, Rp, Fp)
        T = np.interp(s.rv, pr, Tp); rho = np.interp(s.rv, pr, Rp); F = np.interp(s.rv, pr, Fp)
        out = s.rv > s.prmax
        if out.any():   # his atmosphere: isothermal, g ~ (r/rmax)^2, rho floor
            G = 6.674e-8 / (1.0 / (RHO_U * T_UNIT**2))   # grav const in code (rho,l,t units)
            Mmax = PROF[-1, 3] / (RHO_U * RSUN**3)
            gmax = G * Mmax / s.prmax**2
            gl = gmax * (s.rv[out] / s.prmax)**2
            T[out] = Tp[-1]
            rho[out] = np.maximum(Rp[-1] * np.exp(-gl * (s.rv[out] - s.prmax) / Tp[-1]), a.dfloor)
            F[out] = 0.0
        s.Tg = np.interp(s.rg, pr, Tp); s.rhog = np.interp(s.rg, pr, Rp); s.Fg = np.interp(s.rg, pr, Fp)
        s.Lref = s.Fg * s.rg**2            # code luminosity / 4pi at the ghost (= profile L)
        s.T = T.copy()
        s.rho = rho.copy()
        if a.noise > 0:     # bsg_star.cpp: rho*(1+0.05*(U-0.5)) where rho>0.01, P unperturbed
            rng = np.random.default_rng(a.seed)
            m = rho > 0.01
            amp = np.where(m, a.noise, 0.0)
            fac = 1 + amp * (rng.random(N) - 0.5)
            s.rho = rho * fac
            s.T = T / fac                      # eint = tem*rho_unpert/(g-1) -> T = P/rho
        s.I = s.ic_intensity(T**4, F)        # radiation from the unperturbed T
        s.vel = np.zeros(N)
        s.t = 0.0
        s.E0tot = s.etot()
        s.bflux_int = 0.0     # time-integrated inner-face energy flux (into domain >0), code*area
        s.oflux_int = 0.0

    def ic_intensity(s, E, F):
        I = np.empty((len(E), s.NA))
        I[:, s.pos] = (0.5 * (E / s.cw + F / s.cmw))[:, None]
        I[:, ~s.pos] = (0.5 * (E / s.cw - F / s.cmw))[:, None]
        return I

    def ghost_I(s, Iwall):
        a = s.a
        if a.bc == 'bath':
            return s.ic_intensity(np.array([s.Tg**4]), np.array([s.Fg]))[0]
        if a.bc == 'bathP1':     # same E,F but P1 shape E+3 mu F
            return s.Tg**4 + 3 * s.mu * s.Fg
        if a.bc == 'mirror':     # bc_swap 3: I_g(mu) = I_0(-mu) + 6 mu F_in  (handled implicitly)
            return None
        raise ValueError(a.bc)

    def etot(s):
        return np.sum(s.V * (s.rho * s.T / (GAMMA - 1) + PRAT * (s.I @ s.w)))

    # -------------------------------------------------------------- one backward-Euler step
    def step(s, dt):
        a = s.a; N, NA = s.N, s.NA
        ss, sa, spp = opac(s.rho, s.T)
        if a.exch == 'ross':
            spp = sa.copy()
        sg_s, sg_a, _ = opac(np.array([s.rhog]), np.array([s.Tg]))
        chi = ss + sa
        chig = (sg_s + sg_a)[0]
        # face coefficients F_k = A_k I_{k-1} + B_k I_k, k = 0..N
        chiL = np.concatenate([[chig], chi]); chiR = np.concatenate([chi, [chi[-1]]])
        dxL = np.concatenate([[s.dx[0]], s.dx]); dxR = np.concatenate([s.dx, [s.dx[-1]]])
        tau = a.taucell * (dxL * chiL + dxR * chiR)
        fr, fl = taufac(tau)
        if a.flux == 'upwind':
            fr[:] = 1.0; fl[:] = 0.0
        elif a.flux == 'central':
            fr[:] = 1.0; fl[:] = 1.0      # gives central with full diffusion? see below
        v = CRAT * s.mu[None, :]
        av = np.abs(v)
        frk, flk = fr[:, None], fl[:, None]
        sm = frk + flk
        Apos = v * frk * (1 + flk) / sm; Bpos = v * flk * (1 - frk) / sm
        Aneg = -av * flk * (1 - frk) / sm; Bneg = -av * frk * (1 + flk) / sm
        Ac = np.where(v > 0, Apos, Aneg); Bc = np.where(v > 0, Bpos, Bneg)
        if a.flux == 'central':
            Ac = 0.5 * v * np.ones((N + 1, 1)); Bc = 0.5 * v * np.ones((N + 1, 1))
        if a.flux == 'hllfix':     # HLL tau-factor, but with an extra compact diffusion of strength a.kdiff
            pass
        Ig = s.ghost_I(s.I[0])
        nu = NA + 1
        n_tot = N * nu
        I0 = s.I.copy(); T0 = s.T.copy()
        I = s.I.copy(); T = s.T.copy()
        cdt = CRAT * dt
        V = s.V; Af = s.A
        idx = lambda i, n: i * nu + n
        # static linear part of the transport for the Jacobian
        rows = []; cols = []; vals = []
        ii = np.arange(N)
        for n in range(NA):
            # cell i: + dt*(A_{i+1} F_{i+1} - A_i F_i)
            # F_{i+1} = Ac[i+1] I_i + Bc[i+1] I_{i+1};  F_i = Ac[i] I_{i-1} + Bc[i] I_i
            diag = dt * (Af[1:] * Ac[1:, n] - Af[:-1] * Bc[:-1, n])
            up = dt * Af[1:-1] * Bc[1:-1, n]          # coef of I_{i+1} in row i (i<N-1)
            lo = -dt * Af[1:-1] * Ac[1:-1, n]         # coef of I_{i-1} in row i (i>0)
            # outer ghost: I_N = I_{N-1} for mu>0 (copy), 0 for mu<0
            if s.mu[n] > 0:
                diag[-1] += dt * Af[-1] * Bc[-1, n]
            rows += [idx(ii, n), idx(ii[:-1], n), idx(ii[1:], n)]
            cols += [idx(ii, n), idx(ii[1:], n), idx(ii[:-1], n)]
            vals += [diag, up, lo]
        if a.bc == 'mirror':
            # I_g(mu_n) = I_0(-mu_n) + 6 mu F_in ; partner of n is NA-1-n
            for n in range(NA):
                rows.append(np.array([idx(0, n)])); cols.append(np.array([idx(0, NA - 1 - n)]))
                vals.append(np.array([-dt * Af[0] * Ac[0, n]]))
        rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
        Ltr = sp.csr_matrix((vals, (rows, cols)), shape=(n_tot, n_tot))
        # constant transport RHS from the inner ghost
        bvec = np.zeros(n_tot)
        if a.bc == 'mirror':
            Fin = s.Lref / s.rg**2
            Ig_const = 6 * s.mu * Fin
        else:
            Ig_const = Ig
        bvec[idx(0, np.arange(NA))] = -dt * Af[0] * Ac[0, :] * Ig_const
        cell_rows = np.repeat(np.arange(N), NA)
        for it in range(a.newton):
            Ivec = np.zeros(n_tot)
            Ivec[idx(ii[:, None], np.arange(NA)[None, :]).ravel()] = I.ravel()
            Ivec[idx(ii, NA)] = T
            J = I @ s.w
            T4 = T**4
            # residuals
            R = Ltr @ Ivec + bvec
            RI = R[idx(ii[:, None], np.arange(NA)[None, :])]
            RI = RI + V[:, None] * ((I - I0) + cdt * (chi[:, None] * I - (chi - spp)[:, None] * J[:, None] - spp[:, None] * T4[:, None]))
            RT = V * (s.rho * (T - T0) / (GAMMA - 1) + PRAT * cdt * spp * (T4 - J))
            Rv = np.zeros(n_tot)
            Rv[idx(ii[:, None], np.arange(NA)[None, :]).ravel()] = RI.ravel()
            Rv[idx(ii, NA)] = RT
            # local Jacobian blocks
            r_ = []; c_ = []; v_ = []
            for n in range(NA):
                for m in range(NA):
                    val = V * (cdt * (-(chi - spp) * s.w[m]))
                    if m == n:
                        val = val + V * (1 + cdt * chi)
                    r_.append(idx(ii, n)); c_.append(idx(ii, m)); v_.append(val)
                r_.append(idx(ii, n)); c_.append(idx(ii, NA)); v_.append(-V * cdt * spp * 4 * T**3)
                r_.append(idx(ii, NA)); c_.append(idx(ii, n)); v_.append(-V * PRAT * cdt * spp * s.w[n])
            r_.append(idx(ii, NA)); c_.append(idx(ii, NA))
            v_.append(V * (s.rho / (GAMMA - 1) + PRAT * cdt * spp * 4 * T**3))
            Jl = sp.csr_matrix((np.concatenate(v_), (np.concatenate(r_), np.concatenate(c_))), shape=(n_tot, n_tot))
            Jm = (Ltr + Jl).tocsc()
            dX = spl.spsolve(Jm, -Rv)
            dI = dX[idx(ii[:, None], np.arange(NA)[None, :])]
            dT = dX[idx(ii, NA)]
            # damp T to stay positive
            lam = 1.0
            neg = (T + dT) < 0.5 * T
            if neg.any():
                lam = min(1.0, np.min(0.5 * T[neg] / np.abs(dT[neg])))
            I = I + lam * dI; T = T + lam * dT
            err = np.max(np.abs(dT) / T)
            if err < a.tol and lam == 1.0:
                break
        s.newton_it = it + 1
        # boundary energy flux bookkeeping (inner face, outward positive), code units * r^2
        if a.bc == 'mirror':
            Igf = I[0, ::-1] + Ig_const
        else:
            Igf = Ig
        Fface0 = Ac[0] * Igf + Bc[0] * I[0]
        s.Lin_face = (Fface0 @ s.w) * Af[0] / CRAT       # code Fr-units * r^2 (compare s.Lref)
        Iout = np.where(s.pos, I[-1], 0.0)
        FfN = Ac[-1] * I[-1] + Bc[-1] * Iout
        s.Lout_face = (FfN @ s.w) * Af[-1] / CRAT
        s.bflux_int += dt * PRAT * (Fface0 @ s.w) * Af[0]
        s.oflux_int += dt * PRAT * (FfN @ s.w) * Af[-1]
        s.I = I; s.T = T; s.t += dt
        # face fluxes everywhere (for diagnostics)
        Iext_l = np.vstack([Igf[None, :], I]); Iext_r = np.vstack([I, Iout[None, :]])
        s.Fface = ((Ac * Iext_l + Bc * Iext_r) @ s.w) / CRAT   # Fr code units at faces 0..N
        s.chi = chi; s.spp = spp; s.sa = sa

    def diag(s):
        J = s.I @ s.w; F = s.I @ (s.w * s.mu)
        return dict(t=s.t / DAY, T0=s.T[0] * T_U, T1=s.T[1] * T_U, T2=s.T[2] * T_U, T3=s.T[3] * T_U,
                    Er0=J[0] ** 0.25 * T_U, Er1=J[1] ** 0.25 * T_U, Er2=J[2] ** 0.25 * T_U,
                    Lin=s.Lin_face / s.Lref, Lout=s.Lout_face / s.Lref,
                    Fc0=F[0] * s.rv[0]**2 / s.Lref, Fc1=F[1] * s.rv[1]**2 / s.Lref,
                    L5=s.Fface[5] * s.rf[5]**2 / s.Lref, L50=s.Fface[50] * s.rf[50]**2 / s.Lref,
                    L150=s.Fface[150] * s.rf[150]**2 / s.Lref,
                    dE=(s.etot() - s.E0tot), bint=s.bflux_int, oint=s.oflux_int, nit=s.newton_it)


def ecgs(x):
    """code energy (sum V*(...)) -> erg, full sphere"""
    return x * 4 * np.pi * RSUN**3 * ARAD * T_U**4 / PRAT


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('--nx', type=int, default=384); p.add_argument('--rmin', type=float, default=20.0)
    p.add_argument('--rmax', type=float, default=80.0); p.add_argument('--nmu', type=int, default=4)
    p.add_argument('--taucell', type=float, default=4.0); p.add_argument('--dt', type=float, default=212.0)
    p.add_argument('--days', type=float, default=3.0); p.add_argument('--bc', default='bath')
    p.add_argument('--flux', default='hll'); p.add_argument('--exch', default='planck')
    p.add_argument('--noise', type=float, default=0.0); p.add_argument('--seed', type=int, default=0)
    p.add_argument('--dfloor', type=float, default=1e-11)
    p.add_argument('--newton', type=int, default=20); p.add_argument('--tol', type=float, default=1e-10)
    p.add_argument('--out', default='run'); p.add_argument('--every', type=float, default=0.05)
    p.add_argument('--prof', type=float, default=0.5, help='profile dump cadence (d)')
    a = p.parse_args(argv)
    s = Toy(a)
    dt = a.dt / T_UNIT
    nstep = int(round(a.days * DAY / dt))
    hist = []; profs = []
    ev = max(1, int(round(a.every * DAY / dt))); pv = max(1, int(round(a.prof * DAY / dt)))
    t0 = time.time()
    def dump():
        J = s.I @ s.w; F = s.I @ (s.w * s.mu)
        profs.append(dict(t=s.t / DAY, T=s.T.tolist(), J=J.tolist(), F=F.tolist(),
                          Ff=(s.Fface if hasattr(s, 'Fface') else np.zeros(s.N + 1)).tolist(),
                          I=s.I[:12].tolist()))
    dump()
    for k in range(1, nstep + 1):
        s.step(dt)
        if not np.all(np.isfinite(s.T)):
            print('NaN at', s.t / DAY); break
        if k % ev == 0:
            d = s.diag(); hist.append(d)
            print(f"t {d['t']:.3f} T0 {d['T0']:.4e} T1 {d['T1']:.4e} T2 {d['T2']:.4e} Lin {d['Lin']:+.2f} "
                  f"Fc0 {d['Fc0']:+.2f} L5 {d['L5']:+.2f} L50 {d['L50']:+.2f} Lout {d['Lout']:.2f} "
                  f"dE {ecgs(d['dE']):+.3e} bint {ecgs(d['bint']):+.3e} nit {d['nit']} ({time.time()-t0:.0f}s)",
                  flush=True)
        if k % pv == 0:
            dump()
    json.dump(dict(args=vars(a), rv=s.rv.tolist(), rf=s.rf.tolist(), hist=hist, profs=profs,
                   Lref=s.Lref, rho=s.rho.tolist(), unitsF=F_U), open(a.out + '.json', 'w'))


if __name__ == '__main__':
    main()
