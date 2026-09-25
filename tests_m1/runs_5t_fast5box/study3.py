"""Offline: Fourier sideways coarse space per x1 layer (runs_5t_fast5box).

usage: python3 study3.py '<dump glob>' cand ...
cand: fK[r][_mgbL3] : modes |k2|,|k3| <= K (tensor, 'r' = |k2|+|k3| <= K), coarse
first, then rbgs_fwd or mgbL3; also any study2.py candidate.
"""
import os
import sys
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "runs_5m_precond"))
sys.path.insert(0, os.path.join(HERE, "..", "runs_5p_coarse2"))
import study as S          # noqa: E402
import study2 as S2        # noqa: E402
from m1op import Op        # noqa: E402


def basis1(N, K):
    j = np.arange(N) + 0.5
    b = [np.ones(N)]
    for k in range(1, K + 1):
        b.append(np.cos(2*np.pi*k*j/N))
        b.append(np.sin(2*np.pi*k*j/N))
    return b, [0] + [k for k in range(1, K + 1) for _ in (0, 1)]


def four_P(op, K, radial=False, ret_modes=False):
    N3, N2, N1 = op.shape
    b2, k2 = basis1(N2, K)
    b3, k3 = basis1(N3, K) if N3 > 1 else ([np.ones(1)], [0])
    modes = [(a, c) for a in range(len(b2)) for c in range(len(b3))
             if not radial or k2[a] + k3[c] <= K]
    nm = len(modes)
    rows, cols, vals = [], [], []
    idx = np.arange(op.N).reshape(op.shape)
    for q, (a, c) in enumerate(modes):
        ph = np.outer(b3[c], b2[a])            # (N3, N2)
        for i in range(N1):
            rows.append(idx[:, :, i].ravel())
            cols.append(np.full(N3*N2, i*nm + q))
            vals.append(ph.ravel())
    P = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows),
                       np.concatenate(cols))), shape=(op.N, N1*nm))
    if ret_modes:
        return P, nm, [(k2[a], k3[c]) for (a, c) in modes]
    return P, nm


def decouple(Ac, nm, keep):
    """drop the coupling between modes q, q' unless keep(q, q')"""
    c = Ac.tocoo()
    ok = np.array([keep(r % nm, s % nm) for r, s in zip(c.row, c.col)])
    return sp.csr_matrix((c.data[ok], (c.row[ok], c.col[ok])), shape=Ac.shape)


def make(cn, op, A, pc):
    if cn.startswith('f'):
        parts = cn.split('_')
        t = parts[0][1:]
        radial = t.endswith('r')
        K = int(t.rstrip('r'))
        P, nm, md = four_P(op, K, radial, True)
        Ac = (P.T @ A @ P).tocsc()
        if 'mc' in parts:     # mode-diagonal, from the layer-mean stencil coefficients
            parts.remove('mc')
            opm = type('O', (), {})()
            opm.st = np.broadcast_to(op.st.mean(axis=(1, 2), keepdims=True),
                                     op.st.shape).copy()
            opm.shape, opm.N, opm.nst = op.shape, op.N, op.nst
            Am = Op.assemble(opm, range(op.nst))
            Ac = (P.T @ Am @ P).tocsc()
            parts.append('dg')
        if 'dg' in parts:     # mode-diagonal coarse operator
            Ac = decouple(Ac, nm, lambda q, p: q == p).tocsc()
            parts.remove('dg')
        if 'pr' in parts:     # couplings kept within the same (k2, k3)
            Ac = decouple(Ac, nm, lambda q, p: md[q] == md[p]).tocsc()
            parts.remove('pr')
        add = 'add' in parts
        if add:
            parts.remove('add')
        proj = 'proj' in parts
        if proj:
            parts.remove('proj')
            PtP = sp.diags(1.0/np.asarray((P.multiply(P)).sum(axis=0)).ravel())
        lu = spla.splu(Ac)
        rest = '_'.join(parts[1:])
        last = rest.endswith('last')
        if last:
            rest = rest[:-5].rstrip('_')
        inner = S2.make(rest, op, A, pc) if rest else (lambda r: pc.rbgs(r))
        print(f"  {cn}: {nm} modes/layer, coarse n {Ac.shape[0]}", flush=True)
        if add:
            def f(r):
                return inner(r) + P @ lu.solve(P.T @ r)
            return f
        if proj:
            def f(r):
                g = P.T @ r
                return P @ lu.solve(g) + inner(r - P @ (PtP @ g))
            return f
        if last:
            def f(r):
                z = inner(r)
                return z + P @ lu.solve(P.T @ (r - A @ z))
            return f

        def f(r):
            xg = P @ lu.solve(P.T @ r)
            return xg + inner(r - A @ xg)
        return f
    return S2.make(cn, op, A, pc)


def main():
    op = Op(sys.argv[1])
    A = op.A
    pc = S.Precs(op)
    print(f"shape {op.shape} nst {op.nst} blocks {op.nb2}x{op.nb3}", flush=True)
    for cn in sys.argv[2:]:
        prec = make(cn, op, A, pc)
        t0 = time.time()
        x0 = op.x0*float(os.environ.get('X0SCALE', '1'))
        x, nit, h = S.bicgstab(A, op.b, x0, prec)
        print(f"{cn:28s} it {nit:4d}  r0 {h[0]:.2e}  ({time.time()-t0:.1f} s)",
              flush=True)


if __name__ == '__main__':
    main()
