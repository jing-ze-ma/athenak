"""stag.py: staggered two-moment scheme, E at cells, F at faces, F equation solved at every face.

  E_i:  dE/dt = -(Phi_{i+1/2} - Phi_{i-1/2})/dx + c ka (B - E)
        Phi_f = F_f - alpha_f (a_f/2) [E]_f            ([E]_f: jump of reconstructed E at face f)
  F_f:  dF/dt = -(Psi_R - Psi_L)/dx - c chi_f F_f
        Psi_i = c^2 D_i E_i - alpha_i (a_i/2) [F]_i    ([F]_i: jump of reconstructed F across cell i)
alpha = max(A(tau), s_ray) at faces (tau_f, max of the two cells' s) and at cells (tau_i, s_i);
a_f = c sqrt(max(D_L, D_R)), a_i = c sqrt(D_i).  No Rhie-Chow: alpha -> 0 gives the compact
staggered diffusion.  D, C1 from the rays on (E, F interpolated to cells).
Periodic: n faces (face f = left face of cell f).  Inflow: n+1 faces; left ghost cell from the
characteristic split (incoming beam given, outgoing E - F/ag from cell 0 / face 0), right copy.
Low-order cells/faces (positivity limiter): alpha = 1, first-order jumps, speed c; lo_mode 'char'
also uses factor 1 in the pressure c^2 E of a low-order cell (upwinding of E +- F/c).
The operator is affine in (E, F, T) once closure, limiter weights, alpha and the marks are frozen
from the iterate; the sparse matrix is recovered from a few colored probe evaluations."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scheme as S

GAM = S.GAM
STATS = {'steps': 0, 'stages': 0, 'marked': 0, 'resolves': 0, 'be': 0, 'allfaces': 0, 'befail': 0}


def nfaces(pb):
    return pb.n if pb.bc == 'periodic' else pb.n + 1


def split(pb, U):
    n, nf = pb.n, nfaces(pb)
    return U[:n], U[n:n + nf], (U[n + nf:] if pb.gas else None)


def cellF(pb, F):
    """F interpolated to cells."""
    return 0.5*(F + np.roll(F, -1)) if pb.bc == 'periodic' else 0.5*(F[:-1] + F[1:])


def ghostE(pb, E, F):
    """Left ghost-cell E of the inflow boundary (characteristic split)."""
    Eg, Fg, Dg = pb.ghost
    ag = pb.c*np.sqrt(Dg)
    return Eg + 0.5*(E[0] - F[0]/ag)


def faceE(pb, E, F=None):
    """E interpolated to faces (inflow: face 0 uses the ghost cell, face n the edge cell)."""
    if pb.bc == 'periodic':
        return 0.5*(E + np.roll(E, 1))
    e0 = E[0] if F is None else 0.5*(E[0] + ghostE(pb, E, F))
    return np.r_[e0, 0.5*(E[1:] + E[:-1]), E[-1]]


def ext(pb, E, F):
    """Extended cells e = -2..n+1 (index e+2) and faces f = -2..n+2 (index f+2); linear in E, F."""
    n = pb.n
    if pb.bc == 'periodic':
        return E[np.arange(-2, n + 2) % n], F[np.arange(-2, n + 3) % n]
    Eg, Fg, Dg = pb.ghost
    ag = pb.c*np.sqrt(Dg)
    um = E[0] - F[0]/ag                        # outgoing characteristic at the boundary
    Ee = np.r_[[Eg + 0.5*um]*2, E, [E[-1]]*2]
    Fe = np.r_[[Fg - 0.5*ag*um]*2, F, [F[-1]]*2]
    return Ee, Fe


def jumps(ue, w, lo):
    """Jump of the reconstructed field at the interfaces between consecutive entries of ue.
    w = (pl, pr) frozen limiter weights per entry (index matching ue[1:-1]) or None (1st order);
    returns jump for interfaces between ue[k] and ue[k+1], k = 1..len-3."""
    d = ue[2:-1] - ue[1:-2]                    # u_R - u_L
    if w is None:
        return d
    pl, pr = w
    s = pl*(ue[1:-1] - ue[:-2]) + pr*(ue[2:] - ue[1:-1])      # slopes of ue[1:-1]
    jr = d - 0.5*s[1:] - 0.5*s[:-1]
    return np.where(lo, d, jr)


def params(pb, U, hs, opt):
    """Frozen quantities at the iterate: D, s_ray, alpha and speeds at faces and cells, limiter
    weights, low-order masks."""
    n, c, dx = pb.n, pb.c, pb.dx
    E, F, T = split(pb, U)
    if opt.get('Dfix') is not None:
        D, C1 = np.full(n, opt['Dfix']), np.zeros(n)
    else:
        D, C1 = S.rays(pb, E, cellF(pb, F), T)
    s = S.smoothstep((C1 - 0.1)/0.4)
    per = pb.bc == 'periodic'
    nf = nfaces(pb)
    if per:
        iL, iR = np.arange(nf) - 1, np.arange(nf) % n
    else:
        iL, iR = np.clip(np.arange(nf) - 1, 0, n - 1), np.clip(np.arange(nf), 0, n - 1)
    chif = 0.5*(pb.chi[iL] + pb.chi[iR])
    A = lambda tau: np.where(tau > 1e-8, -np.expm1(-tau)/np.maximum(tau, 1e-300), 1.0)
    lo_f = opt.get('lo_f', np.zeros(nf, bool))
    lo_c = opt.get('lo_c', np.zeros(n, bool))
    p = dict(D=D,
             al_f=np.where(lo_f, 1.0, np.maximum(A(chif*dx), np.maximum(s[iL], s[iR]))),
             a_f=np.where(lo_f, c, c*np.sqrt(np.maximum(D[iL], D[iR]))),
             al_c=np.where(lo_c, 1.0, np.maximum(A(pb.chi*dx), s)),
             a_c=np.where(lo_c, c, c*np.sqrt(D)),
             Dp=np.where(lo_c & (opt.get('lo_mode') == 'char'), 1.0, D),
             chif=chif, lo_f=lo_f, lo_c=lo_c)
    Ee, Fe = ext(pb, E, F)
    if opt.get('recon'):
        p['wE'] = S._select(Ee, opt['recon'])            # entries Ee[1:-1]
        p['wF'] = S._select(Fe, opt['recon'])
    return p


def operator(pb, U, p):
    """dU/dt for state U with the frozen parameters p (affine in U)."""
    n, c, dx = pb.n, pb.c, pb.dx
    nf = nfaces(pb)
    per = pb.bc == 'periodic'
    E, F, T = split(pb, U)
    Ee, Fe = ext(pb, E, F)
    # E flux at faces 0..nf-1: cells L = f-1, R = f -> Ee indices f+1, f+2
    wE = p.get('wE')
    jE = jumps(Ee, wE, np.zeros(n + 1, bool))[:nf] if wE is not None else (Ee[2:-1] - Ee[1:-2])[:nf]
    if wE is not None:
        jE = np.where(p['lo_f'], (Ee[2:-1] - Ee[1:-2])[:nf], jE)
    Phi = F - 0.5*p['al_f']*p['a_f']*jE
    PhiR = np.roll(Phi, -1) if per else Phi[1:]
    PhiL = Phi if per else Phi[:-1]
    dE = -(PhiR - PhiL)/dx
    # F 'flux' at cells -1..n: faces i (Fe index i+2) and i+1 (index i+3)
    wF = p.get('wF')
    dFj = Fe[2:-1] - Fe[1:-2]                  # cells -1..n  (F_{i+1} - F_i)
    jF = jumps(Fe, wF, np.zeros(dFj.size, bool)) if wF is not None else dFj
    cells = np.arange(-1, n + 1)
    ci = cells % n if per else np.clip(cells, 0, n - 1)
    Dp, alc, ac = p['Dp'][ci], p['al_c'][ci], p['a_c'][ci]
    if wF is not None:
        jF = np.where(p['lo_c'][ci], dFj, jF)
    Ec = Ee[1:-1]                               # cells -1..n
    if not per:
        Dp = Dp.copy()
        Dp[0] = pb.ghost[2]
    Psi = c*c*Dp*Ec - 0.5*alc*ac*jF             # index i+1 for cell i
    # face f between cells f-1 (Psi index f) and f (index f+1)
    dF = -(Psi[1:nf + 1] - Psi[:nf])/dx - c*p['chif']*F
    if pb.gas:
        g = pb.gas
        B = g['B0'] + g['Bp']*(T - g['T0'])
        ex = c*pb.ka*(B - E)
        return np.r_[dE + ex, dF, -ex/g['rho_cv']]
    return np.r_[dE + c*pb.ka*(pb.B - E), dF]


def assemble(pb, U, p):
    """Sparse M, b with operator(U) = M U + b, by colored probes (stencil radius < 4)."""
    n, nf = pb.n, nfaces(pb)
    per = pb.bc == 'periodic'
    N = U.size
    b = operator(pb, np.zeros(N), p)
    pos = np.r_[np.arange(n), np.arange(nf) - 0.5]
    if pb.gas:
        pos = np.r_[pos, np.arange(n)]
    blocks = [(0, n, 0.0), (n, nf, -0.5)] + ([(n + nf, n, 0.0)] if pb.gas else [])
    K = 8 if n % 8 == 0 else (10 if n % 10 == 0 else None)
    rows, cols, vals = [], [], []
    for off, size, shift in blocks:
        Kb = K if K is not None else size
        for col in range(Kb):
            v = np.zeros(N)
            v[off + col:off + size:Kb] = 1.0
            r = operator(pb, v, p) - b
            nz = np.nonzero(r)[0]
            # column of this color nearest to the row's position
            q = pos[nz] - shift - col
            j = col + Kb*np.round(q/Kb).astype(int)
            if per:
                j = j % size
            else:
                j = np.where(j < 0, j + Kb, np.where(j >= size, j - Kb, j))
            rows.append(nz); cols.append(off + j); vals.append(r[nz])
    M = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(N, N))
    return M, b


def stage_solve(pb, rhs, guess, hs, opt):
    Y = guess
    I = sp.identity(guess.size, format='csr')
    for _ in range(opt.get('npic', 3)):
        p = params(pb, Y, hs, opt)
        M, b = assemble(pb, Y, p)
        Y = spla.spsolve((I - hs*M).tocsc(), rhs + hs*b)
    return Y


def upwindE(pb, E, F):
    """E at faces from the upwind cell (with respect to the sign of F_f)."""
    n = pb.n
    if pb.bc == 'periodic':
        return np.where(F >= 0, np.roll(E, 1), E)
    EL, ER = np.r_[ghostE(pb, E, F), E], np.r_[E, E[-1]]
    return np.where(F >= 0, EL, ER)


FACECHECK = 'avg'     # face test |F_f| <= c E_f with E_f = mean of the cells ('avg') or upwind ('upw')


def bad(pb, V, tol=1e-12):
    """Non-realizable cells and faces with co-located values: (E_i, F at cells), (E at faces, F_f)."""
    E, F, _ = split(pb, V)
    ref = max(E.max(), 1e-300)
    c = pb.c
    bc = (E < -tol*ref) | (np.abs(cellF(pb, F)) - c*E > tol*c*ref) | ~np.isfinite(E)
    Ef = faceE(pb, E, F) if FACECHECK == 'avg' else upwindE(pb, E, F)
    bf = (np.abs(F) - c*Ef > tol*c*ref) | ~np.isfinite(F)
    return bc, bf


def mood_stage(pb, rhs, guess, hs, opt, lo_f, lo_c):
    """Stage solve + a posteriori marking: a bad cell marks itself and its two faces, a bad face
    marks itself and its two cells; re-solve up to 3 times, marks accumulate."""
    o = dict(opt, lo_f=lo_f, lo_c=lo_c)
    Y = stage_solve(pb, rhs, guess, hs, o)
    STATS['stages'] += 1
    n, nf = pb.n, nfaces(pb)
    per = pb.bc == 'periodic'
    for _ in range(3):
        bc, bf = bad(pb, Y)
        if not (bc.any() or bf.any()):
            break
        ic, jf = np.nonzero(bc)[0], np.nonzero(bf)[0]
        mf = np.r_[ic, ic + 1, jf]
        mc = np.r_[ic, jf - 1, jf]
        mf = mf % nf if per else mf[(mf >= 0) & (mf < nf)]
        mc = mc % n if per else mc[(mc >= 0) & (mc < n)]
        STATS['marked'] += int((~lo_f[mf]).sum() + (~lo_c[mc]).sum())
        lo_f[mf], lo_c[mc] = True, True
        STATS['resolves'] += 1
        Y = stage_solve(pb, rhs, Y, hs, o)
    return Y


def step(pb, U, h, opt):
    """SDIRK2 step; with opt['mood'] the positivity limiter and the BE / all-low-order chain."""
    n, nf = pb.n, nfaces(pb)
    z = lambda: (np.zeros(nf, bool), np.zeros(n, bool))
    if not opt.get('mood'):
        Y1 = stage_solve(pb, U, U, GAM*h, opt)
        k1 = (Y1 - U)/(GAM*h)
        return stage_solve(pb, U + (1 - GAM)*h*k1, Y1, GAM*h, opt), Y1
    STATS['steps'] += 1
    lo_f, lo_c = z()
    Y1 = mood_stage(pb, U, U, GAM*h, opt, lo_f, lo_c)
    k1 = (Y1 - U)/(GAM*h)
    Y2 = mood_stage(pb, U + (1 - GAM)*h*k1, Y1, GAM*h, opt, lo_f, lo_c)
    if not (any(x.any() for x in bad(pb, Y1)) or any(x.any() for x in bad(pb, Y2))):
        return Y2, Y1
    STATS['be'] += 1
    lo_f, lo_c = z()
    Y = mood_stage(pb, U, U, h, opt, lo_f, lo_c)
    if any(x.any() for x in bad(pb, Y)):
        STATS['allfaces'] += 1
        Y = stage_solve(pb, U, U, h, dict(opt, lo_f=np.ones(nf, bool), lo_c=np.ones(n, bool)))
    if any(x.any() for x in bad(pb, Y)):
        STATS['befail'] += 1
    return Y, None


def run(pb, E, F, tend, h, opt, T=None, check=None):
    U = np.r_[E, F] if T is None else np.r_[E, F, T]
    nst = int(np.ceil(tend/h - 1e-9))
    for _ in range(nst):
        U, Y1 = step(pb, U, tend/nst, opt)
        if check is not None:
            check(U, Y1)
    return U
