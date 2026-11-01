"""Hodge-split deep horizontal filter of the radial-gauge edge integrals (hclean_1031).

Per radial cell c the horizontal face fluxes are the 1-form g(c) = L(c+1) - L(c) on the
edges of the shell grid (gauge invariant; both copies of a panel-seam edge carry the same
value x orientation sign).  On the closed cubed-sphere shell g = d0 chi + d1^T psi exactly
(no harmonic part): d0 chi is curl free (horizontally solenoidal flux: leaves B_r alone),
d1^T psi is the minimum-norm part with curl = the B_r flux increment dP1(c).
The filter:
  - chi(c) -> chi + wc(c) (M_v - 1) chi       (vertex tensor 1-2-1, M_v)
  - P1(n)  -> P1 + wf(n) (M_c - 1) P1         (cell tensor 1-2-1 on the radial fluxes)
  - psi(c) from the new dP1(c);  g_new = d0 chi_new + d1^T psi_new
  - L_new(0) = L(0) + d1^T (d1 d1^T)^-1 dP1(0);  L_new(n+1) = L_new(n) + g_new(n)
div B = 0 by construction (curl of edge integrals); every field on and above the first
face with wf = 0 and every cell with wc = 0 (and wf = 0 on both faces) is unchanged.
M = (F1 F2 + (F1 F2)^T)/2, Fk = 1-2-1 along the k-th grid direction (opposite neighbour
pairs on the global graph; panel corners (valence 3) are left unfiltered); M is symmetric,
so it conserves the total radial flux of every shell."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.spatial import cKDTree

FRAME = (((0, 1, 0), (0, 0, 1), (1, 0, 0)), ((-1, 0, 0), (0, 0, 1), (0, 1, 0)),
         ((0, -1, 0), (0, 0, 1), (-1, 0, 0)), ((0, 1, 0), (-1, 0, 0), (0, 0, 1)),
         ((1, 0, 0), (0, 0, 1), (0, -1, 0)), ((0, 1, 0), (1, 0, 0), (0, 0, -1)))


def pmap(p, xi, eta):
    a, b, n = (np.array(v, float) for v in FRAME[p])
    x, y = np.tan(xi), np.tan(eta)
    return (a*x[..., None] + b*y[..., None] + n)/np.sqrt(1 + x*x + y*y)[..., None]


def _dedupe(X, tol=1e-9):
    """global ids of points X (n, 3): first occurrence is canonical"""
    t = cKDTree(X)
    gid = -np.ones(len(X), int)
    ng = 0
    for i in range(len(X)):
        if gid[i] >= 0:
            continue
        for j in t.query_ball_point(X[i], tol):
            gid[j] = ng
        ng += 1
    return gid, ng


class Shell:
    def __init__(self, N):
        self.N = N
        fq = 0.25*np.pi*(-1.0 + 2.0*np.arange(N + 1)/N)
        cq = 0.25*np.pi*(-1.0 + 2.0*(np.arange(N) + 0.5)/N)
        P = np.arange(6)
        hh = 1e-7
        # vertices (p, K, J) at (xi = fq[J], eta = fq[K])
        XV = np.stack([pmap(p, fq[None, :] + 0*fq[:, None], fq[:, None] + 0*fq[None, :])
                       for p in P])                              # (6, N+1, N+1, 3)
        self.vid, self.nv = _dedupe(XV.reshape(-1, 3))
        self.vid = self.vid.reshape(6, N + 1, N + 1)
        # edges: L2 (p, K, j) along xi at (cq[j], fq[K]); L3 (p, k, J) along eta
        X2 = np.stack([pmap(p, cq[None, :] + 0*fq[:, None], fq[:, None] + 0*cq[None, :])
                       for p in P])                              # (6, N+1, N, 3)
        D2 = np.stack([pmap(p, cq[None, :] + hh + 0*fq[:, None], fq[:, None] + 0*cq[None, :])
                       - pmap(p, cq[None, :] - hh + 0*fq[:, None], fq[:, None]
                              + 0*cq[None, :]) for p in P])
        X3 = np.stack([pmap(p, fq[None, :] + 0*cq[:, None], cq[:, None] + 0*fq[None, :])
                       for p in P])                              # (6, N, N+1, 3)
        D3 = np.stack([pmap(p, fq[None, :] + 0*cq[:, None], cq[:, None] + hh + 0*fq[None, :])
                       - pmap(p, fq[None, :] + 0*cq[:, None], cq[:, None] - hh
                              + 0*fq[None, :]) for p in P])
        XE = np.concatenate((X2.reshape(-1, 3), X3.reshape(-1, 3)))
        DE = np.concatenate((D2.reshape(-1, 3), D3.reshape(-1, 3)))
        DE /= np.linalg.norm(DE, axis=-1)[:, None]
        eid, self.ne = _dedupe(XE)
        canon = np.zeros(self.ne, int)
        canon[eid[::-1]] = np.arange(len(eid))[::-1]             # first occurrence
        sg = np.sign((DE*DE[canon[eid]]).sum(-1))
        n2 = X2[..., 0].size
        self.e2, self.s2 = eid[:n2].reshape(6, N + 1, N), sg[:n2].reshape(6, N + 1, N)
        self.e3, self.s3 = eid[n2:].reshape(6, N, N + 1), sg[n2:].reshape(6, N, N + 1)
        # d0 (E x V) from the canonical copies
        rows, cols, vals = [], [], []
        for p in P:
            for (e, s, va, vb) in ((self.e2[p], self.s2[p], self.vid[p][:, :-1],
                                    self.vid[p][:, 1:]),
                                   (self.e3[p], self.s3[p], self.vid[p][:-1, :],
                                    self.vid[p][1:, :])):
                rows += [e.ravel(), e.ravel()]
                cols += [vb.ravel(), va.ravel()]
                vals += [s.ravel(), -s.ravel()]
        d0 = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows),
                                                   np.concatenate(cols))),
                           shape=(self.ne, self.nv)).tocsr()
        d0.data = np.sign(d0.data)                     # duplicates (seam copies) -> +-1
        self.d0 = d0
        # d1 (C x E): P1 = L3[k, j+1] - L3[k, j] - L2[k+1, j] + L2[k, j]
        self.nc = 6*N*N
        cid = np.arange(self.nc).reshape(6, N, N)
        rows, cols, vals = [], [], []
        for (e, s, sl, f) in ((self.e3, self.s3, np.s_[:, :, 1:], 1),
                              (self.e3, self.s3, np.s_[:, :, :-1], -1),
                              (self.e2, self.s2, np.s_[:, 1:, :], -1),
                              (self.e2, self.s2, np.s_[:, :-1, :], 1)):
            rows.append(cid.ravel()); cols.append(e[sl].ravel())
            vals.append(f*s[sl].ravel())
        self.d1 = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows),
                                                        np.concatenate(cols))),
                                shape=(self.nc, self.ne)).tocsr()
        self.cid = cid
        # Poisson operators (pinned at node 0; right-hand sides sum to zero)
        Lc = (self.d1 @ self.d1.T).tocsc()[1:, 1:]
        Lv = (self.d0.T @ self.d0).tocsc()[1:, 1:]
        self.solc, self.solv = spl.splu(Lc), spl.splu(Lv)
        # smoothers
        XC = np.stack([pmap(p, cq[None, :] + 0*cq[:, None], cq[:, None] + 0*cq[None, :])
                       for p in P])
        exc = np.stack([pmap(p, cq[None, :] + hh + 0*cq[:, None], cq[:, None]
                             + 0*cq[None, :]) for p in P]) - XC
        exv = np.stack([pmap(p, fq[None, :] + hh + 0*fq[:, None], fq[:, None]
                             + 0*fq[None, :]) for p in P]) - XV
        Ac = abs(self.d1) @ abs(self.d1).T          # cell adjacency via shared edges
        Av = abs(self.d0).T @ abs(self.d0)
        self.Mc, self.ncorner_c = self._smoother(Ac, XC.reshape(-1, 3), exc.reshape(-1, 3),
                                                 np.arange(self.nc))
        # vertex positions / local xi direction from the first copy
        first = np.zeros(self.nv, int)
        flat = self.vid.ravel()
        first[flat[::-1]] = np.arange(len(flat))[::-1]
        self.Mv, self.ncorner_v = self._smoother(Av, XV.reshape(-1, 3)[first],
                                                 exv.reshape(-1, 3)[first],
                                                 np.arange(self.nv))

    @staticmethod
    def _smoother(A, X, ex, ids):
        A = sp.csr_matrix(A)
        A.setdiag(0); A.eliminate_zeros()
        n = A.shape[0]
        r1, c1, r2, c2 = [], [], [], []
        ncorner = 0
        for i in range(n):
            nb = A.indices[A.indptr[i]:A.indptr[i + 1]]
            if len(nb) != 4:
                ncorner += 1
                continue
            d = X[nb] - X[i]
            d -= (d @ X[i])[:, None]*X[i]
            d /= np.linalg.norm(d, axis=-1)[:, None]
            e = ex[i] - (ex[i] @ X[i])*X[i]
            e /= np.linalg.norm(e)
            al = np.abs(d @ e)
            o = np.argsort(-al)
            pa, pb = nb[o[:2]], nb[o[2:]]
            if (d[o[0]] @ d[o[1]]) > -0.5 or (d[o[2]] @ d[o[3]]) > -0.5:
                raise SystemExit('smoother: neighbour pairing failed at %d' % i)
            r1 += [i, i]; c1 += list(pa)
            r2 += [i, i]; c2 += list(pb)
        I = sp.identity(n, format='csr')

        def F(r, c):
            # symmetric 1-2-1: pair weights symmetrised, so every F is doubly stochastic
            # (conserves the shell sum; constants kept); differs from the plain stencil
            # only where the two nodes disagree on the direction (panel seams)
            N_ = sp.coo_matrix((np.full(len(r), 0.25), (r, c)), shape=(n, n)).tocsr()
            N_ = 0.5*(N_ + N_.T)
            dg = np.asarray(N_.sum(1)).ravel()
            return I - sp.diags(dg) + N_
        P = F(r1, c1) @ F(r2, c2)
        return (0.5*(P + P.T)).tocsr(), ncorner

    # panel arrays <-> global edge vectors
    def to_edges(self, A2, A3):
        """A2 (6, N+1, N, R), A3 (6, N, N+1, R) -> (E, R) from the canonical copies,
        and the max copy mismatch"""
        R = A2.shape[-1]
        g = np.zeros((self.ne, R))
        cnt = np.zeros(self.ne)
        v2 = (A2*self.s2[..., None]).reshape(-1, R)
        v3 = (A3*self.s3[..., None]).reshape(-1, R)
        e = np.concatenate((self.e2.ravel(), self.e3.ravel()))
        v = np.concatenate((v2, v3))
        np.add.at(g, e, v); np.add.at(cnt, e, 1)
        g /= cnt[:, None]
        mis = np.abs(v - g[e]).max()
        return g, mis

    def to_panels(self, g):
        R = g.shape[-1]
        return (g[self.e2.ravel()].reshape(self.e2.shape + (R,))*self.s2[..., None],
                g[self.e3.ravel()].reshape(self.e3.shape + (R,))*self.s3[..., None])

    def solve_c(self, rhs):
        x = np.zeros((self.nc,) + rhs.shape[1:])
        x[1:] = self.solc.solve(rhs[1:])
        return x

    def solve_v(self, rhs):
        x = np.zeros((self.nv,) + rhs.shape[1:])
        x[1:] = self.solv.solve(rhs[1:])
        return x


def taper(n, nfull, nzero):
    i = np.arange(n)
    return np.where(i <= nfull, 1.0, np.where(
        i >= nzero, 0.0, np.cos(0.5*np.pi*(i - nfull)/max(nzero - nfull, 1))**2))


def hodge_filter(L2, L3, N, cfull, czero, passes=1, filter_br=True):
    """L2 (6, N+1, N, I+1), L3 (6, N, N+1, I+1) radial-gauge edge integrals; cells
    c <= cfull full, cos^2 taper to 0 at cell czero; radial faces n <= cfull + 1 full,
    0 from face czero.  Returns new L2, L3 and diagnostics."""
    sh = Shell(N)
    nf = L2.shape[-1]
    wc = taper(nf - 1, cfull, czero)
    wf = taper(nf, cfull + 1, czero)
    if not filter_br:
        wf = 0*wf
    G, mis_g = sh.to_edges(L2[..., 1:] - L2[..., :-1], L3[..., 1:] - L3[..., :-1])
    # radial fluxes per face from the panel L (curl is gauge-free per panel)
    P1 = (L3[:, :, 1:] - L3[:, :, :-1]) - (L2[:, 1:] - L2[:, :-1])   # (6, N, N, nf)
    P1 = P1.reshape(sh.nc, nf)
    chk_d1d0 = abs(sh.d1 @ sh.d0).max()
    dP1 = sh.d1 @ G                                       # (C, nf-1)
    chk_curl = np.abs(dP1 - (P1[:, 1:] - P1[:, :-1])).max()/np.abs(P1).max()
    psi = sh.solve_c(dP1)
    Gp = sh.d1.T @ psi
    chi = sh.solve_v(sh.d0.T @ (G - Gp))
    chk_hodge = np.abs(sh.d0 @ chi + Gp - G).max()/np.abs(G).max()
    Mv, Mc = sh.Mv, sh.Mc
    chin, P1n = chi.copy(), P1.copy()
    for _ in range(passes):
        chin = chin + wc*(Mv @ chin - chin)
        P1n = P1n + wf*(Mc @ P1n - P1n)
    dP1n = P1n[:, 1:] - P1n[:, :-1]
    Gn = sh.d0 @ chin + sh.d1.T @ sh.solve_c(dP1n)
    keep = (wc == 0) & (wf[:-1] == 0) & (wf[1:] == 0)
    Gn[:, keep] = G[:, keep]
    dL0 = sh.d1.T @ sh.solve_c((P1n[:, 0] - P1[:, 0])[:, None])     # (E, 1)
    dL = np.concatenate((dL0, dL0 + np.cumsum(Gn - G, axis=1)), axis=1)  # (E, nf)
    d2, d3 = sh.to_panels(dL)
    L2n, L3n = L2 + d2, L3 + d3
    P1c = ((L3n[:, :, 1:] - L3n[:, :, :-1]) - (L2n[:, 1:] - L2n[:, :-1])).reshape(sh.nc, nf)
    diag = dict(wc=wc, wf=wf, mis_g=mis_g/np.abs(G).max(), d1d0=chk_d1d0,
                curl=chk_curl, hodge=chk_hodge, ncorner=(sh.ncorner_c, sh.ncorner_v),
                p1_target=np.abs(P1c - P1n).max()/np.abs(P1).max(),
                monopole=np.abs(P1n.sum(0)).max()/np.abs(P1).max(),
                nv=sh.nv, ne=sh.ne, nc=sh.nc)
    return L2n, L3n, diag
