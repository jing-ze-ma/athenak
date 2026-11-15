"""scheme.py: 1-D cell-centred implicit two-moment (E, F) grey radiation, SDIRK2 in time,
closure (D) and streaming weight (s_ray) from a steady formal solution on the current state.

  dE/dt + dF/dx = c ka (B - E)          dF/dt + c^2 d(D E)/dx = -c chi F,  chi = ka + ks
Face f between cells L, R:
  F_E,f = (F_L+F_R)/2 - alpha_f (a_f/2) [E]_f + w_rc RC_f
  P_f   = c^2 (D_L E_L + D_R E_R)/2 - alpha_f (a_f/2) [F]_f
  [u]_f = u_R^- - u_L^+  (1st order: u_R - u_L; 2nd: limited linear reconstruction)
  RC_f  = theta_f h c^2 [ (g_L+g_R)/2 - ((DE)_R-(DE)_L)/dx ], g_i central gradient of DE,
  theta_f = 1/(1 + h c chi_f), h = the stage step (gamma dt), a_f = c sqrt(max(D_L, D_R)).
Optional gas: rho_cv dT/dt = -c ka (B(T) - E), B linearised about T0 (B = B0 + Bp (T-T0)).
Each SDIRK stage = npic Picard passes; each pass: rays + limiter selection from the
current iterate (lagged), assemble the linear system, one sparse solve."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numba import njit

GAM = 1.0 - 1.0/np.sqrt(2.0)
_mu, _w = np.polynomial.legendre.leggauss(16)       # 8 per hemisphere
MU, WT = _mu, 0.5*_w                                # sum WT = 1: J = sum WT I (E units)


@njit(cache=True)
def _cell(Ic, S0, S1, Sm, t):
    """Exact integration across one cell, source linear in tau from S0 (entry) to S1 (exit),
    Sm = cell mean.  Returns (cell-averaged I, exit I)."""
    if t < 1e-4:
        e, g, q = 1.0 - t + 0.5*t*t, 1.0 - 0.5*t + t*t/6.0, 0.5 - t/6.0
    else:
        e = np.exp(-t)
        g = -np.expm1(-t)/t
        q = (1.0 - g)/t
    dS = S1 - S0
    return Sm + (Ic - S0)*g - dS*q, S1 + (Ic - S0)*e - dS*g


@njit(cache=True)
def _sweep2(S, Sl, Sr, chi, dx, mus, periodic, Iin):
    """Exact integration of every ray across the cells.  Source per cell linear between the
    face values Sl[i], Sr[i] (left/right face of cell i; Sl = Sr = S gives piecewise-constant
    S).  Returns cell-averaged Iavg[m, i] and face Ifc[m, f].
    Periodic: exact wrap (pass 1 gives the transmission A and I_end, I_0 = I_end/(1-A))."""
    nm, n = mus.size, S.size
    Iavg = np.zeros((nm, n))
    Ifc = np.zeros((nm, n + 1))
    for m in range(nm):
        mu = mus[m]
        amu = abs(mu)
        I0 = Iin[m]
        if periodic:
            Ie, lnA = 0.0, 0.0
            for jj in range(n):
                i = jj if mu > 0 else n - 1 - jj
                S0, S1 = (Sl[i], Sr[i]) if mu > 0 else (Sr[i], Sl[i])
                t = chi[i]*dx/amu
                Ie = _cell(Ie, S0, S1, S[i], t)[1]
                lnA -= t
            A = np.exp(lnA)
            I0 = Ie/(1.0 - A) if 1.0 - A > 1e-12 else 0.0
        Ic = I0
        if mu > 0:
            Ifc[m, 0] = I0
        else:
            Ifc[m, n] = I0
        for jj in range(n):
            i = jj if mu > 0 else n - 1 - jj
            S0, S1 = (Sl[i], Sr[i]) if mu > 0 else (Sr[i], Sl[i])
            Iavg[m, i], Ic = _cell(Ic, S0, S1, S[i], chi[i]*dx/amu)
            if mu > 0:
                Ifc[m, i + 1] = Ic
            else:
                Ifc[m, i] = Ic
    return Iavg, Ifc


SLIN = False    # True: piecewise-linear source (face values = mean of the two cells)


def smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x*x*(3.0 - 2.0*x)


class Problem:
    """Grid, opacities, boundary data.  bc = 'periodic' or 'inflow' (left ghost = (Eg, Fg, Dg),
    right = outflow copy; rays: left incoming beam Ibeam in the largest mu, vacuum at right)."""

    def __init__(self, n, L=1.0, c=1.0, ka=0.0, ks=0.0, B=0.0, bc='periodic',
                 ghost=None, gas=None):
        self.n, self.dx, self.c = n, L/n, c
        self.x = (np.arange(n) + 0.5)*self.dx
        self.ka = np.broadcast_to(np.asarray(ka, float), (n,)).copy()
        self.ks = np.broadcast_to(np.asarray(ks, float), (n,)).copy()
        self.chi = self.ka + self.ks
        self.B = np.broadcast_to(np.asarray(B, float), (n,)).copy()
        self.bc = bc
        if ghost is None:      # inflow beam consistent with the rays' beam (largest GL mu)
            mm = MU.max()
            ghost = (1.0, c*mm, mm*mm)
        self.ghost = ghost
        self.gas = gas      # dict(rho_cv, T0, B0, Bp) or None
        self.nv = 3 if gas else 2


# ---------------------------------------------------------------- closure from the rays
RAYX = {}       # extras of the last rays() call: kinetic face flux, ray-described mask


def rays(pb, E, F, T=None):
    """Formal solution on the state: S = (ka B + ks E)/chi.  Returns D, C1 per cell.
    Cells the rays leave empty (J < 1e-3 E) take neighbour values (8 passes); beyond that
    D = max(1/3, f^2) (beam -> 1, isotropic -> 1/3), C1 = 1."""
    B = pb.B if T is None else pb.gas['B0'] + pb.gas['Bp']*(T - pb.gas['T0'])
    chi = pb.chi
    S = np.where(chi > 0, (pb.ka*B + pb.ks*E)/np.where(chi > 0, chi, 1.0), 0.0)
    Iin = np.zeros(MU.size)
    if pb.bc == 'inflow':
        Iin[np.argmax(MU)] = pb.ghost[0]/WT[np.argmax(MU)]
    if SLIN:
        if pb.bc == 'periodic':
            Sf = 0.5*(S + np.roll(S, 1))
            Sf = np.r_[Sf, Sf[0]]
        else:
            Sf = np.r_[S[0], 0.5*(S[1:] + S[:-1]), S[-1]]
        # cells with chi = 0 have S = 0 and do not interact; keep their faces neutral
        Sl, Sr = Sf[:-1].copy(), Sf[1:].copy()
    if not SLIN:
        Sl, Sr = S, S
    Ia, Ifc = _sweep2(S, Sl, Sr, chi, pb.dx, MU, pb.bc == 'periodic', Iin)
    RAYX['Fray'] = pb.c*((WT*MU) @ Ifc)         # rays' net flux through faces 0..n
    J = WT @ Ia
    K = (WT*MU**2) @ Ia
    C1n = WT @ np.abs(Ia - S[None, :])
    Emax = max(np.abs(E).max(), 1e-300)
    ok = (J > 1e-3*np.maximum(E, 0.0)) & (J > 1e-14*Emax)
    D = np.where(ok, K/np.where(ok, J, 1.0), 0.0)
    C1 = np.where(ok, C1n/np.where(ok, J, 1.0), 0.0)
    per = pb.bc == 'periodic'
    for _ in range(8):
        if ok.all():
            break
        Dn, Cn, cnt = np.zeros_like(D), np.zeros_like(D), np.zeros_like(D)
        for sh in (1, -1):
            if per:
                o, d, cc = np.roll(ok, sh), np.roll(D, sh), np.roll(C1, sh)
            else:
                o, d, cc = np.zeros_like(ok), np.zeros_like(D), np.zeros_like(D)
                if sh == 1:
                    o[1:], d[1:], cc[1:] = ok[:-1], D[:-1], C1[:-1]
                else:
                    o[:-1], d[:-1], cc[:-1] = ok[1:], D[1:], C1[1:]
            Dn += o*d
            Cn += o*cc
            cnt += o
        new = (~ok) & (cnt > 0)
        D[new], C1[new] = Dn[new]/cnt[new], Cn[new]/cnt[new]
        ok = ok | new
    RAYX['ok'] = ok.copy()
    f = np.where(E > 1e-14*Emax, np.abs(F)/(pb.c*np.where(E > 0, E, 1.0)), 0.0)
    D = np.where(ok, D, np.maximum(1.0/3.0, np.minimum(f, 1.0)**2))
    C1 = np.where(ok, C1, 1.0)
    return np.clip(D, 0.0, 1.0), C1      # clip only matters for non-realizable states


# ---------------------------------------------------------------- spatial operator
def _ext(pb, u, gval):
    """State on extended cells e = -2..n+1 (index e+2)."""
    n = pb.n
    e = np.arange(-2, n + 2)
    if pb.bc == 'periodic':
        return u[e % n]
    out = u[np.clip(e, 0, n - 1)].astype(float)
    out[:2] = gval
    return out


def _select(ue, lim):
    """Lagged limiter selection: slope_e = pl*(u_e - u_{e-1}) + pr*(u_{e+1} - u_e) for the
    extended cells e = -1..n (array index e+1).  lim: 'minmod', 'mc' or 'vl' (van Leer,
    slope 2 dl dr/(dl+dr) = [dr/(dl+dr)] dl + [dl/(dl+dr)] dr, weights lagged)."""
    dl = ue[1:-1] - ue[:-2]
    dr = ue[2:] - ue[1:-1]
    same = dl*dr > 0
    pl, pr = np.zeros_like(dl), np.zeros_like(dl)
    if lim == 'vl':
        den = np.where(same, dl + dr, 1.0)
        pl[same], pr[same] = (dr/den)[same], (dl/den)[same]
    elif lim == 'minmod':
        useL = np.abs(dl) <= np.abs(dr)
        pl[same & useL], pr[same & ~useL] = 1.0, 1.0
    else:
        cand = np.stack([np.abs(2*dl), np.abs(2*dr), np.abs(0.5*(dl + dr))])
        k = np.argmin(cand, axis=0)
        pl[same & (k == 0)] = 2.0
        pr[same & (k == 1)] = 2.0
        pl[same & (k == 2)], pr[same & (k == 2)] = 0.5, 0.5
    return pl, pr


def jameson(ue, eps=1e-12):
    """Smoothness sensor |u+ - 2u + u-|/(|u+ - u| + |u - u-| + eps) on extended cells
    e = -1..n (index e+1): ~dx in smooth fields, ~1 at jumps and 2dx modes."""
    d2 = np.abs(ue[2:] - 2.0*ue[1:-1] + ue[:-2])
    d1 = np.abs(ue[2:] - ue[1:-1]) + np.abs(ue[1:-1] - ue[:-2])
    return d2/(d1 + eps*(np.abs(ue).max() + 1e-300))


def alpha_weight(cand, tau, theta, s, psi=None):
    A1 = np.where(tau > 1e-8, -np.expm1(-tau)/np.where(tau > 1e-8, tau, 1.0), 1.0 - 0.5*tau)
    if cand == 'A1':
        return A1
    if cand == 'A2':                       # 1 - (1-A1)(1-theta)
        return 1.0 - (1.0 - A1)*(1.0 - theta)
    if cand == 'A2max':
        return np.maximum(A1, theta)
    if cand == 'A3':                       # A1 with the ray streaming weight
        return np.maximum(A1, s)
    if cand == 'A3b':                      # A2 with the ray streaming weight
        return np.maximum(1.0 - (1.0 - A1)*(1.0 - theta), s)
    if cand == 'rayflux':                  # computed in assemble, passed as psi
        return psi
    if cand == 'S1':                       # linear: A1 * smoothness sensor, max ray weight
        return np.maximum(A1*psi, s)
    if cand == 'A1sq':
        return A1**2
    if cand == 'one':
        return np.ones_like(tau)
    if cand == 'zero':
        return np.zeros_like(tau)
    raise ValueError(cand)


def assemble(pb, U, hs, opt):
    """Linear system dU/dt = M U + b at the iterate U (closure, limiter lagged).
    opt: cand, rc ('oma'|'omt'|'one'|'off'), recon (None|'minmod'|'mc'), Dfix."""
    n, c, dx = pb.n, pb.c, pb.dx
    E, F = U[:n], U[n:2*n]
    T = U[2*n:] if pb.gas else None
    if opt.get('Dfix') is not None:
        D, C1 = np.full(n, opt['Dfix']), np.zeros(n)
    else:
        D, C1 = rays(pb, E, F, T)
    if opt.get('Dguard'):           # realizability floor D >= f^2 (state flux; NOT rays-only)
        f = np.abs(F)/(c*np.maximum(E, 1e-300))
        D = np.maximum(D, np.minimum(f, 1.0)**2)
    Eg, Fg, Dg = pb.ghost
    # characteristic inflow ghost: incoming u+ = Eg + Fg/ag given, outgoing u- from cell 0:
    #   E_g = Eg + (E0 - F0/ag)/2,  F_g = Fg - ag (E0 - F0/ag)/2   (cst, coef E0, coef F0)
    ag = c*np.sqrt(Dg)
    gmapE, gmapF = (Eg, 0.5, -0.5/ag), (Fg, -0.5*ag, 0.5)
    gE = gmapE[0] + gmapE[1]*E[0] + gmapE[2]*F[0]
    gF = gmapF[0] + gmapF[1]*E[0] + gmapF[2]*F[0]
    De, Ce, chie = _ext(pb, D, Dg), _ext(pb, C1, 1.0), _ext(pb, pb.chi, pb.chi[0])
    per = pb.bc == 'periodic'
    nf = n if per else n + 1
    f = np.arange(nf)                      # face f: L = f-1, R = f (extended cell index)
    iL, iR = f - 1 + 2, f + 2              # positions in extended arrays
    chif = 0.5*(chie[iL] + chie[iR])
    tau = chif*dx
    theta = 1.0/(1.0 + hs*c*chif)
    s = smoothstep((np.maximum(Ce[iL], Ce[iR]) - 0.1)/0.4)
    psi = None
    if opt['cand'] == 'rayflux':
        # alpha_f = clip((Fc - Fray)/((a/2)(E_R - E_L)), 0, A1); faces with no ray information
        # on either side (rays left both cells empty) take A1.  Uses the iterate (stage start
        # with npic = 1) and its rays; fixed during the stage.
        Ee, Fe = _ext(pb, E, gE), _ext(pb, F, gF)
        okE = _ext(pb, RAYX['ok'].astype(float), 1.0) if opt.get('Dfix') is None \
            else np.ones(n + 4)
        Fray = RAYX['Fray'][:nf]
        dE = Ee[iR] - Ee[iL]
        Fc = 0.5*(Fe[iL] + Fe[iR])
        aa = c*np.sqrt(np.maximum(De[iL], De[iR]))
        big = np.abs(dE) > 1e-12*(np.abs(Ee[iL]) + np.abs(Ee[iR]) + 1e-300)
        A1 = alpha_weight('A1', tau, theta, s)
        raw = np.where(big, (Fc - Fray)/np.where(big, 0.5*aa*dE, 1.0), 0.0)
        psi = np.where((okE[iL] + okE[iR]) > 0, np.clip(raw, 0.0, A1), A1)
    if opt['cand'] == 'S1':
        pc = np.maximum(jameson(_ext(pb, E, gE)), jameson(_ext(pb, F/c, gF/c)))
        psi = np.maximum(pc[f], pc[f + 1])          # cells L = f-1, R = f
    al = alpha_weight(opt['cand'], tau, theta, s, psi)
    RAYX['alpha'] = al
    asp = opt.get('aspeed', 'sqrtD')
    a = c*np.sqrt(np.maximum(De[iL], De[iR])) if asp == 'sqrtD' \
        else np.full(nf, c*(2.0 if asp == '2c' else 1.0))
    wrc = {'oma': 1.0 - al, 'omt': 1.0 - theta, 'one': np.ones(nf),
           'off': np.zeros(nf)}[opt.get('rc', 'oma')]
    # jump coefficients on offsets k = -2,-1,0,+1 (cells f-2 .. f+1)
    def jump(ue):
        J = np.zeros((nf, 4))
        J[:, 1], J[:, 2] = -1.0, 1.0
        if opt.get('recon'):
            pl, pr = _select(ue, opt['recon'])   # index e+1
            plR, prR = pl[f + 1], pr[f + 1]      # cell R = f
            plL, prL = pl[f], pr[f]              # cell L = f-1
            J[:, 1] += 0.5*plR
            J[:, 2] += -0.5*(plR - prR)
            J[:, 3] += -0.5*prR
            J[:, 0] += 0.5*plL
            J[:, 1] += -0.5*(plL - prL)
            J[:, 2] += -0.5*prL
        return J
    JE = jump(_ext(pb, E, gE))
    JF = jump(_ext(pb, F, gF))
    lo = opt.get('lo')
    if lo is not None and lo.any():
        # low-order faces (positivity limiter): 1st-order states, alpha = 1, speed lo_speed,
        # no Rhie-Chow; D and C1 are unchanged
        J1 = np.zeros(4)
        J1[1], J1[2] = -1.0, 1.0
        JE[lo], JF[lo] = J1, J1
        al = np.where(lo, 1.0, al)
        a = np.where(lo, c*opt.get('lo_speed', 1.0), a)
        wrc = np.where(lo, 0.0, wrc)
    q = np.array([-0.25, -0.25 + 1.0, 0.25 - 1.0, 0.25])/dx     # RC on DE
    offs = np.arange(-2, 2)
    eidx = f[:, None] + offs[None, :]                            # extended cell e
    Dk = De[eidx + 2]
    cEE = -(al*a/2)[:, None]*JE + (wrc*theta*hs*c*c)[:, None]*q[None, :]*Dk
    cEF = np.zeros((nf, 4))
    cEF[:, 1] = cEF[:, 2] = 0.5
    cFE = np.zeros((nf, 4))
    cFE[:, 1], cFE[:, 2] = 0.5*c*c*Dk[:, 1], 0.5*c*c*Dk[:, 2]
    if lo is not None and lo.any() and opt.get('lo_mode') == 'char':
        # characteristic upwinding of E +- F/c with speeds +-c on low-order faces: the face
        # pressure average uses the free-streaming factor 1 (cell D and C1 unchanged)
        cFE[lo, 1], cFE[lo, 2] = 0.5*c*c, 0.5*c*c
    cFF = -(al*a/2)[:, None]*JF
    # map extended cells to unknowns or ghost constants
    if per:
        col = eidx % n
        const = np.zeros_like(eidx, dtype=bool)
    else:
        col = np.clip(eidx, 0, n - 1)
        const = eidx < 0
    rows, cols, vals = [], [], []
    b = np.zeros(pb.nv*n)
    for rowoff, cE, cF, gq in ((0, cEE, cEF, None), (n, cFE, cFF, None)):
        for cell, sg in ((f, 1.0/dx), (f - 1, -1.0/dx)):
            valid = (cell >= 0) & (cell < n) if not per else np.ones(nf, bool)
            rc_ = (cell % n) if per else cell
            for coef, var, gv in ((cE, 0, gmapE), (cF, n, gmapF)):
                R = np.broadcast_to((rowoff + rc_)[:, None], coef.shape)
                m = valid[:, None] & ~const
                rows.append(R[m]); cols.append(var + col[m]); vals.append(sg*coef[m])
                mc = valid[:, None] & const
                if mc.any():
                    np.add.at(b, R[mc], sg*coef[mc]*gv[0])
                    for gc, gcol in ((gv[1], 0), (gv[2], n)):
                        rows.append(R[mc]); cols.append(np.full(R[mc].shape, gcol))
                        vals.append(sg*coef[mc]*gc)
    i = np.arange(n)
    rows += [i, n + i]
    cols += [i, n + i]
    vals += [-c*pb.ka, -c*pb.chi]
    if pb.gas:
        g = pb.gas
        rows += [i, 2*n + i, 2*n + i]
        cols += [2*n + i, i, 2*n + i]
        vals += [c*pb.ka*g['Bp'], c*pb.ka/g['rho_cv'], -c*pb.ka*g['Bp']/g['rho_cv']]
        b0 = g['B0'] - g['Bp']*g['T0']
        b[:n] += c*pb.ka*b0
        b[2*n:] -= c*pb.ka*b0/g['rho_cv']
    else:
        b[:n] += c*pb.ka*pb.B
    M = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(pb.nv*n, pb.nv*n))
    return M, b


HIST = None      # set to a list to record Picard convergence: (pass, rel change, n switches)


def stage_solve(pb, rhs, guess, hs, opt):
    """Solve Y = rhs + hs (M(Y) Y + b(Y)) by npic Picard passes (closure+limiter lagged).
    npic = 'conv': iterate until the relative change < 1e-12 (max 60 passes)."""
    Y = guess
    I = sp.identity(pb.nv*pb.n, format='csr')
    npic = opt.get('npic', 3)
    nmax = 60 if npic == 'conv' else npic
    hist = []
    for k in range(nmax):
        M, b = assemble(pb, Y, hs, opt)
        Yn = spla.spsolve((I - hs*M).tocsc(), rhs + hs*b)
        ch = np.abs(Yn - Y).max()/max(np.abs(Yn).max(), 1e-300)
        hist.append(ch)
        Y = Yn
        if npic == 'conv' and ch < 1e-12:
            break
    if HIST is not None:
        HIST.append(hist)
    return Y


def realizable(pb, V, tol=1e-12):
    n = pb.n
    e, f = V[:n], V[n:2*n]
    ref = max(e.max(), 1e-300)
    return e.min() >= -tol*ref and (np.abs(f) - pb.c*e).max() <= tol*pb.c*ref


STATS = {'steps': 0, 'stages': 0, 'marked': 0, 'resolves': 0, 'be': 0, 'befail': 0,
         'allfaces': 0}


def bad_cells(pb, V, tol=1e-12):
    n = pb.n
    e, f = V[:n], V[n:2*n]
    ref = max(e.max(), 1e-300)
    return (e < -tol*ref) | (np.abs(f) - pb.c*e > tol*pb.c*ref) | ~np.isfinite(e + f)


def mood_stage(pb, rhs, guess, hs, opt, lo):
    """Stage solve with the a posteriori positivity limiter: faces of non-realizable cells
    (plus opt['mood_nb'] neighbours) switch to the low-order flux, re-solve, up to 3 times;
    the marks in lo accumulate."""
    Y = stage_solve(pb, rhs, guess, hs, dict(opt, lo=lo))
    STATS['stages'] += 1
    per = pb.bc == 'periodic'
    nf = pb.n if per else pb.n + 1
    for _ in range(opt.get('mood_rounds', 3)):
        bad = np.nonzero(bad_cells(pb, Y))[0]
        if bad.size == 0:
            break
        nb = opt.get('mood_nb', 0)
        fa = np.concatenate([bad + k for k in range(-nb, 2 + nb)])   # faces i-nb .. i+1+nb
        fa = fa % nf if per else fa[(fa >= 0) & (fa < nf)]
        new = ~lo[fa]
        STATS['marked'] += int(new.sum())
        lo[fa] = True
        STATS['resolves'] += 1
        Y = stage_solve(pb, rhs, Y, hs, dict(opt, lo=lo))
    return Y


def step(pb, U, h, opt):
    """One SDIRK2 step (Alexander, L-stable, stiffly accurate); opt['ti']='be' for BE.
    opt['fallback']: if a stage is not realizable, redo the step with BE (1st order, same
    alpha), then BE with alpha = 1, a = c.  Returns (U_new, stage-1 state or None)."""
    if opt.get('mood'):
        STATS['steps'] += 1
        nf = pb.n if pb.bc == 'periodic' else pb.n + 1
        lo = np.zeros(nf, bool)
        if opt.get('ti') != 'be':
            Y1 = mood_stage(pb, U, U, GAM*h, opt, lo)
            k1 = (Y1 - U)/(GAM*h)
            Y2 = mood_stage(pb, U + (1.0 - GAM)*h*k1, Y1, GAM*h, opt, lo)
            if not (bad_cells(pb, Y1).any() or bad_cells(pb, Y2).any()):
                return Y2, Y1
            STATS['be'] += 1
        Y = mood_stage(pb, U, U, h, opt, np.zeros(nf, bool))
        if bad_cells(pb, Y).any() and opt.get('mood_all'):
            # last resort: the whole step low order (BE, every face)
            STATS['allfaces'] = STATS.get('allfaces', 0) + 1
            Y = stage_solve(pb, U, U, h, dict(opt, lo=np.ones(nf, bool)))
        if bad_cells(pb, Y).any():
            STATS['befail'] += 1
        return Y, None
    if opt.get('ti') == 'be':
        return stage_solve(pb, U, U, h, opt), None
    Y1 = stage_solve(pb, U, U, GAM*h, opt)
    k1 = (Y1 - U)/(GAM*h)
    Y2 = stage_solve(pb, U + (1.0 - GAM)*h*k1, Y1, GAM*h, opt)
    if opt.get('fallback') and not (realizable(pb, Y1) and realizable(pb, Y2)):
        o = dict(opt, recon=None)
        Y2 = stage_solve(pb, U, U, h, o)
        if not realizable(pb, Y2):
            Y2 = stage_solve(pb, U, U, h, dict(o, cand='one', aspeed='c'))
        return Y2, None
    return Y2, Y1


def run(pb, E, F, tend, h, opt, T=None, check=None):
    """Advance to tend with steps <= h; optional check(U, stageY) after each step."""
    U = np.r_[E, F] if T is None else np.r_[E, F, T]
    t = 0.0
    nst = int(np.ceil(tend/h - 1e-9))
    hh = tend/nst
    for _ in range(nst):
        U, Y1 = step(pb, U, hh, opt)
        t += hh
        if check is not None:
            check(U, Y1)
    return U
