"""Offline preconditioner study on dumped implicit M1 systems (runs_5m_precond).

usage: python study.py '<dump glob>' [candidate ...]
Right-preconditioned BiCGStab from x0 with the code's test max|r| < tol max|b|.
"""
import sys
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from m1op import Op


def thomas(a, b, c, r):
    """tridiagonal solve along the last axis (a: lower, c: upper; a[...,0], c[...,-1]
    ignored)"""
    n = r.shape[-1]
    cp = np.empty_like(r)
    dp = np.empty_like(r)
    cp[..., 0] = c[..., 0]/b[..., 0]
    dp[..., 0] = r[..., 0]/b[..., 0]
    for i in range(1, n):
        den = b[..., i] - a[..., i]*cp[..., i-1]
        cp[..., i] = c[..., i]/den
        dp[..., i] = (r[..., i] - a[..., i]*dp[..., i-1])/den
    x = np.empty_like(r)
    x[..., -1] = dp[..., -1]
    for i in range(n-2, -1, -1):
        x[..., i] = dp[..., i] - cp[..., i]*x[..., i+1]
    return x


class Precs:
    def __init__(self, op):
        self.op = op
        f = op.f
        self.TA, self.TB, self.TC = f['TA'], f['TB'], f['TC']
        N3, N2, N1 = op.shape
        K, J = np.meshgrid(np.arange(N3), np.arange(N2), indexing='ij')
        self.kl = K % op.n3
        self.jl = J % op.n2
        self.col = (self.kl + self.jl) & 1          # block-local colour
        self.colg = (K + J) & 1                      # global colour
        self.C = {'jm': f['CJM'], 'jp': f['CJP'], 'km': f['CKM'], 'kp': f['CKP']}

    def line(self, r):
        return thomas(self.TA, self.TB, self.TC, r.reshape(self.op.shape)).ravel()

    def trans(self, z, local):
        """sum_nb C_nb z_nb (the transverse 5-point coupling), block-local or global"""
        C = self.C
        out = np.zeros_like(z)
        sh = {'jm': (-1, 1), 'jp': (1, 1), 'km': (-1, 0), 'kp': (1, 0)}
        for key, (d, ax) in sh.items():
            zn = np.roll(z, -d, axis=ax)             # z at (j+d) or (k+d)
            t = C[key]*zn
            if local:
                lo = self.jl if ax == 1 else self.kl
                n = self.op.n2 if ax == 1 else self.op.n3
                edge = (lo == 0) if d < 0 else (lo == n - 1)
                t = np.where(edge[..., None], 0.0, t)
            out += t
        return out

    def rbgs(self, r, local=True, sym=False, nsweep=1):
        r3 = r.reshape(self.op.shape)
        colr = self.col if local else self.colg
        z = np.zeros_like(r3)
        order = [0, 1] * nsweep + ([0] if sym else [])
        for c in order:
            rhs = r3 - self.trans(z, local)
            zc = thomas(self.TA, self.TB, self.TC, rhs)
            m = (colr == c)[..., None]
            z = np.where(m, zc, z)
        return z.ravel()


class Level:
    """one level of a (j,k) semicoarsening multigrid with x1-line red-black GS"""
    def __init__(self, A, shape, local=None):
        self.A = A.tocsr()
        self.shape = shape
        N3, N2, N1 = shape
        n = N1*N2*N3
        K, J, I = np.meshgrid(np.arange(N3), np.arange(N2), np.arange(N1), indexing='ij')
        colid = (K*N2 + J).ravel()
        Ac = self.A.tocoo()
        same = colid[Ac.row] == colid[Ac.col]
        D = sp.csr_matrix((Ac.data[same], (Ac.row[same], Ac.col[same])), shape=(n, n))
        self.O = sp.csr_matrix((Ac.data[~same], (Ac.row[~same], Ac.col[~same])),
                               shape=(n, n))
        if local is not None:     # local = (n2, n3): drop couplings across blocks from
            b2, b3 = local        # the smoother (they stay in A): the code's rbgs
            blk = ((K//b3)*1000 + J//b2).ravel()
            oc = self.O.tocoo()
            keep = blk[oc.row] == blk[oc.col]
            self.O = sp.csr_matrix((oc.data[keep], (oc.row[keep], oc.col[keep])),
                                   shape=(n, n))
            col = (((K % b3) + (J % b2)) & 1).ravel()
        else:
            col = ((K + J) & 1).ravel()
        self.lu = spla.splu(D.tocsc())
        self.red = col == 0
        self.nx = shape

    def smooth(self, r, z, order=(0, 1)):
        for c in order:
            zc = self.lu.solve(r - self.O @ z)
            m = self.red if c == 0 else ~self.red
            z = np.where(m, zc, z)
        return z


def agg_P(shape, c2, c3):
    N3, N2, N1 = shape
    M3, M2 = -(-N3//c3), -(-N2//c2)
    K, J, I = np.meshgrid(np.arange(N3)//c3, np.arange(N2)//c2, np.arange(N1),
                          indexing='ij')
    agg = ((K*M2 + J)*N1 + I).ravel()
    n = N1*N2*N3
    P = sp.csr_matrix((np.ones(n), (np.arange(n), agg)), shape=(n, N1*M2*M3))
    return P, (M3, M2, N1)


class MG:
    def __init__(self, op, A, nmin=4, local=True, post=True, cf=1.0, smoothP=0.0,
                 nlev=99):
        self.levels, self.P = [], []
        shape = op.shape
        lv = Level(A, shape, (op.n2, op.n3) if local else None)
        self.levels.append(lv)
        while min(shape[0], shape[1]) > nmin and len(self.levels) < nlev:
            P, cshape = agg_P(shape, 2 if shape[1] > 1 else 1, 2 if shape[0] > 1 else 1)
            if smoothP > 0.0:     # smoothed aggregation: P <- (I - w D^-1 A) P
                Dinv = sp.diags(1.0/lv.A.diagonal())
                P = (P - smoothP*(Dinv @ (lv.A @ P))).tocsr()
            Ac = (P.T @ lv.A @ P).tocsr()
            self.P.append(P)
            shape = cshape
            lv = Level(Ac, shape)
            self.levels.append(lv)
        self.clu = spla.splu(self.levels[-1].A.tocsc())
        self.post, self.cf = post, cf
        print("  MG levels:", [l.shape for l in self.levels], flush=True)

    def vcycle(self, r, l=0):
        if l == len(self.levels) - 1:
            return self.clu.solve(r)
        lv = self.levels[l]
        z = lv.smooth(r, np.zeros_like(r))
        rr = r - lv.A @ z
        z = z + self.cf*(self.P[l] @ self.vcycle(self.P[l].T @ rr, l + 1))
        if self.post:
            z = lv.smooth(r, z, order=(1, 0))
        return z


def blk_ids(shape, b2, b3):
    N3, N2, N1 = shape
    K, J, I = np.meshgrid(np.arange(N3), np.arange(N2), np.arange(N1), indexing='ij')
    return ((K//b3)*100000 + J//b2).ravel()


def drop_cross(A, bid):
    c = A.tocoo()
    keep = bid[c.row] == bid[c.col]
    return sp.csr_matrix((c.data[keep], (c.row[keep], c.col[keep])), shape=A.shape)


class MGB:
    """block-local (j,k) semicoarsening MG: every level lives inside the MeshBlock (the
    aggregates never cross a block face).  Level 0 smooths block-locally (rbgs_fwd) and
    forms its residual with the full A (one halo).  Levels >= 1 use the Galerkin
    operator with the cross-block couplings DROPPED (no communication), except the
    bottom (one cell per block and x1 index) which, with glob=True, keeps them and is
    solved exactly over all blocks (the only global step)."""
    def __init__(self, op, A, glob=True, post=False, nsm=1, fres_local=False,
                 maxlev=99, bot='exact'):
        shape = op.shape
        b2, b3 = op.n2, op.n3
        self.post, self.nsm = post, nsm
        self.levels, self.P, self.Ares = [], [], []
        lv = Level(A, shape, (b2, b3))
        self.levels.append(lv)
        Afull = A
        self.Ares.append(drop_cross(A, blk_ids(shape, b2, b3)) if fres_local else A)
        while b2 > 1 or b3 > 1:
            N3, N2, N1 = shape
            c2, c3 = -(-b2//2), -(-b3//2)
            K, J, I = np.meshgrid(np.arange(N3), np.arange(N2), np.arange(N1),
                                  indexing='ij')
            Jc = (J//b2)*c2 + (J % b2)//2
            Kc = (K//b3)*c3 + (K % b3)//2
            M2, M3 = (N2//b2)*c2, (N3//b3)*c3
            agg = ((Kc*M2 + Jc)*N1 + I).ravel()
            n = N1*N2*N3
            P = sp.csr_matrix((np.ones(n), (np.arange(n), agg)), shape=(n, N1*M2*M3))
            Afull = (P.T @ Afull @ P).tocsr()
            shape, b2, b3 = (M3, M2, N1), c2, c3
            Aloc = drop_cross(Afull, blk_ids(shape, b2, b3))
            self.P.append(P)
            if len(self.levels) >= maxlev and not (b2 == 1 and b3 == 1):
                lvb = Level(Aloc, shape, (b2, b3))
                if bot == 'line':
                    self.clu = type('L', (), {'solve': lambda s_, r: lvb.lu.solve(r)})()
                else:     # one rbgs_fwd sweep
                    self.clu = type('L', (), {'solve': lambda s_, r:
                                              lvb.smooth(r, np.zeros_like(r))})()
                self.levels.append(None)
                break
            if b2 == 1 and b3 == 1:
                Abot = Afull if glob else Aloc
                self.clu = spla.splu(Abot.tocsc())
                self.levels.append(None)
            else:
                self.levels.append(Level(Aloc, shape, (b2, b3)))
                self.Ares.append(Aloc)
        print("  MGB levels:", len(self.levels), "bottom", shape, flush=True)

    def vcycle(self, r, l=0):
        if self.levels[l] is None:
            return self.clu.solve(r)
        lv = self.levels[l]
        z = np.zeros_like(r)
        for _ in range(self.nsm):
            z = lv.smooth(r, z)
        rr = r - self.Ares[l] @ z
        z = z + self.P[l] @ self.vcycle(self.P[l].T @ rr, l + 1)
        if self.post:
            z = lv.smooth(r, z, order=(1, 0))
        return z


def bicgstab(A, b, x0, prec, tol=1e-10, maxit=400):
    bs = np.max(np.abs(b))
    x = x0.copy()
    r = b - A @ x
    rh = r.copy()
    rho = alpha = omega = 1.0
    v = np.zeros_like(b)
    p = np.zeros_like(b)
    hist = [np.max(np.abs(r))/bs]
    nit = 0
    while hist[-1] > tol and nit < maxit:
        nit += 1
        rhon = rh @ r
        beta = (rhon/rho)*(alpha/omega)
        p = r + beta*(p - omega*v)
        y = prec(p)
        v = A @ y
        alpha = rhon/(rh @ v)
        s = r - alpha*v
        zz = prec(s)
        t = A @ zz
        omega = (t @ s)/(t @ t)
        x += alpha*y + omega*zz
        r = s - omega*t
        rho = rhon
        hist.append(np.max(np.abs(r))/bs)
    return x, nit, hist


def coarse_setup(op, A, a1, a2, a3):
    """piecewise-constant aggregation P (aggregates a1 x a2 x a3 cells), Galerkin A_c"""
    N3, N2, N1 = op.shape
    K, J, I = np.meshgrid(np.arange(N3)//a3, np.arange(N2)//a2, np.arange(N1)//a1,
                          indexing='ij')
    c3, c2, c1 = -(-N3//a3), -(-N2//a2), -(-N1//a1)
    agg = ((K*c2 + J)*c1 + I).ravel()
    P = sp.csr_matrix((np.ones(op.N), (np.arange(op.N), agg)), shape=(op.N, c1*c2*c3))
    Ac = (P.T @ A @ P).tocsc()
    lu = spla.splu(Ac)
    return P, lu


def main():
    op = Op(sys.argv[1])
    A = op.A
    pc = Precs(op)
    print(f"shape {op.shape} nst {op.nst} blocks {op.nb2}x{op.nb3}")
    s1 = op.st[1][..., 0]
    s2 = op.st[2][..., -1]
    print(f"x1 boundary coefs (must be 0): {np.abs(s1).max():.3e} {np.abs(s2).max():.3e}")
    # consistency of the preconditioner rows with the stencil
    for nm, s in (('TA', 1), ('TB', 0), ('TC', 2), ('CJM', 3), ('CJP', 4), ('CKM', 5),
                  ('CKP', 6)):
        sgn = 1.0
        d = np.abs(op.f[nm] - sgn*op.st[s]).max()/np.abs(op.st[s]).max()
        print(f"  {nm} vs {sgn:+.0f} st[{s}]: rel diff {d:.2e}")
    cands = sys.argv[2:] or ['line', 'rbgs_fwd']
    if cands[0] == 'diag':
        diagnose(op, A, pc)
        return
    for cn in cands:
        prec = make_prec(cn, op, A, pc)
        t0 = time.time()
        x, nit, h = bicgstab(A, op.b, op.x0, prec)
        print(f"{cn:28s} it {nit:4d}  r0 {h[0]:.2e}  ({time.time()-t0:.1f} s)", flush=True)


def make_prec(cn, op, A, pc):
    if cn == 'line':
        return pc.line
    if cn == 'rbgs_fwd':
        return lambda r: pc.rbgs(r)
    if cn == 'rbgs':
        return lambda r: pc.rbgs(r, sym=True)
    if cn == 'rbgs_fwd_glob':
        return lambda r: pc.rbgs(r, local=False)
    if cn.startswith('rbgs_n'):       # rbgs_n<k>: k forward sweeps, block-local
        k = int(cn[6:])
        return lambda r: pc.rbgs(r, nsweep=k)
    if cn.startswith('rbgs_n') or cn.startswith('rbglob_n'):
        k = int(cn[8:])
        return lambda r: pc.rbgs(r, local=False, nsweep=k)
    if cn.startswith('tl_'):          # tl_<a1>_<a2>_<a3>: rbgs_fwd + coarse (mult.)
        a1, a2, a3 = (int(v) for v in cn[3:].split('_'))
        P, lu = coarse_setup(op, A, a1, a2, a3)

        def f(r):
            z = pc.rbgs(r)
            rr = r - A @ z
            return z + P @ lu.solve(P.T @ rr)
        return f
    if cn.startswith('mgb'):   # mgb[_loc][_post][_ns<k>][_fl]
        mg = MGB(op, A, glob='loc' not in cn, post='post' in cn,
                 nsm=int(cn.split('_ns')[1][0]) if '_ns' in cn else 1,
                 fres_local='_fl' in cn,
                 maxlev=int(cn.split('_L')[1][0]) if '_L' in cn else 99,
                 bot='line' if '_bl' in cn else 'rb')
        return mg.vcycle
    if cn.startswith('mg'):   # mg[_nopost][_glob][_cf<x>][_sp<w>][_nl<k>]
        kw = dict(post='nopost' not in cn, local='glob' not in cn)
        for t in cn.split('_'):
            if t.startswith('cf'):
                kw['cf'] = float(t[2:])
            if t.startswith('sp'):
                kw['smoothP'] = float(t[2:])
            if t.startswith('nl'):
                kw['nlev'] = int(t[2:])
        mg = MG(op, A, **kw)
        return mg.vcycle
    raise ValueError(cn)


def diagnose(op, A, pc, nits=(5, 10, 20)):
    """where the slow error lives: row excess and coupling ratios by depth, then the
    error left after n rbgs_fwd BiCGStab iterations, by depth and sideways wavenumber"""
    st = op.st
    offs = np.abs(st[1:]).sum(axis=0)
    sig = (st[0] - offs)/st[0]
    x1c = (np.abs(st[1]) + np.abs(st[2]))/st[0]
    trc = np.abs(st[3:7]).sum(axis=0)/st[0]
    N3, N2, N1 = op.shape
    print("depth profile (i: mean over j,k): row excess sigma, x1 coupling, sideways "
          "coupling (fractions of the diagonal), |b|")
    for i in list(range(0, N1, max(1, N1//12))) + [N1-1]:
        print(f"  i={i:3d} sigma {sig[..., i].mean():.3e} (min {sig[..., i].min():.2e})"
              f"  x1 {x1c[..., i].mean():.3e}  side {trc[..., i].mean():.3e}"
              f"  |b| {np.abs(op.f['b'][..., i]).mean():.3e}")
    xs, nit, _ = bicgstab(A, op.b, op.x0, lambda r: pc.rbgs(r), tol=1e-14)
    print(f"reference solve to 1e-14: {nit} it")
    for n in nits:
        x, _, h = bicgstab(A, op.b, op.x0, lambda r: pc.rbgs(r), tol=0.0, maxit=n)
        e = (xs - x).reshape(op.shape)
        r = (op.b - A @ x).reshape(op.shape)
        ei = np.sqrt((e**2).mean(axis=(0, 1)))
        ri = np.abs(r).max(axis=(0, 1))
        imx = int(np.argmax(ri))
        print(f"after {n} it: max|r|/max|b| {h[-1]:.2e} at depth i={imx};"
              f" rel err max {np.abs(e).max()/np.abs(xs).max():.2e}")
        print("   |r|max by depth (every 8th): "
              + " ".join(f"{v/np.abs(op.b).max():.1e}" for v in ri[::8]))
        # sideways spectrum of the error at the depth of the largest residual
        ek = np.abs(np.fft.fft2(e[..., imx]))**2
        k3 = np.fft.fftfreq(N3)*N3
        k2 = np.fft.fftfreq(N2)*N2
        kk = np.sqrt(k3[:, None]**2 + k2[None, :]**2)
        tot = ek.sum()
        for lo, hi in ((0, 0.5), (0.5, 1.5), (1.5, 2.5), (2.5, 4.5), (4.5, 8.5),
                       (8.5, 16.5), (16.5, 1e9)):
            msk = (kk >= lo) & (kk < hi)
            print(f"     |k| in [{lo},{hi}): {ek[msk].sum()/tot:.3f}")
        # error vs block-local position (block edges vs interior) at that depth
        jl = np.arange(N2) % op.n2
        edge = (jl == 0) | (jl == op.n2 - 1)
        print(f"     rms err at block-edge columns / interior: "
              f"{np.sqrt((e[:, edge, :]**2).mean()):.2e} / "
              f"{np.sqrt((e[:, ~edge, :]**2).mean()):.2e}")


if __name__ == '__main__':
    main()
