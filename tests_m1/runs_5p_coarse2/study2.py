"""Offline study of global coarse corrections on dumped implicit M1 systems
(runs_5p_coarse2; builds on runs_5m_precond/study.py).

usage: python3 study2.py '<dump glob>' cand [cand ...]
cands: rbgs_fwd | mgbL<k>[p] (block-local levels, k levels, p = post-smoothing) |
       gc[_mgbL<k>[p]] (global per-layer sideways mean, coarse first, then the smoother
       or the block-local levels) | mg_glob[_nopost] | mg[_nopost] (global levels,
       block-local smoother) | any study.py candidate
"""
import os
import sys
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                                "runs_5m_precond"))
import study as S          # noqa: E402
from m1op import Op        # noqa: E402


def mean_P(op, gb2=1, gb3=1):
    N3, N2, N1 = op.shape
    K, J, X1 = np.meshgrid(np.arange(N3), np.arange(N2), np.arange(N1), indexing='ij')
    b = ((K*gb3)//N3)*gb2 + (J*gb2)//N2
    agg = (X1*(gb2*gb3) + b).ravel()
    return sp.csr_matrix((np.ones(op.N), (np.arange(op.N), agg)),
                         shape=(op.N, N1*gb2*gb3))


def make(cn, op, A, pc):
    if cn.startswith('gc'):
        parts = cn.split('_')
        gb = [int(t[1:]) for t in parts if t.startswith('b') and t[1:].isdigit()]
        gb2 = gb[0] if gb else 1
        P = mean_P(op, gb2, gb2 if op.shape[0] > 1 else 1)
        lu = spla.splu((P.T @ A @ P).tocsc())
        rest = [t for t in parts[1:] if not (t.startswith('b') and t[1:].isdigit())]
        inner = make('_'.join(rest), op, A, pc) if rest else (lambda r: pc.rbgs(r))

        def f(r):
            xg = P @ lu.solve(P.T @ r)
            return xg + inner(r - A @ xg)
        return f
    if cn.startswith('mgbL'):
        k = int(cn[4])
        mg = S.MGB(op, A, glob=False, post=cn.endswith('p'), maxlev=k, bot='rb')
        return mg.vcycle
    return S.make_prec(cn, op, A, pc)


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
