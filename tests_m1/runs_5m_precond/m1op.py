"""Load an implicit_dump_op file set and build the global operator (runs_5m_precond).

A dump m1op.c<N>.r<rank>.bin holds, per MeshBlock, the stored stencil st(0..nst-1),
b, x0 and the preconditioner rows TA, TB, TC, CJM, CJP, CKM, CKP (active cells).
Assumes one MeshBlock in x1, a uniform level, and periodic x2 / x3 (both study cases).
"""
import glob
import numpy as np
import scipy.sparse as sp

# stencil slot -> (dk, dj, di)
OFF = {0: (0, 0, 0), 1: (0, 0, -1), 2: (0, 0, 1), 3: (0, -1, 0), 4: (0, 1, 0),
       5: (-1, 0, 0), 6: (1, 0, 0),
       7: (0, -1, -1), 8: (0, -1, 1), 9: (0, 1, -1), 10: (0, 1, 1),
       11: (-1, 0, -1), 12: (-1, 0, 1), 13: (1, 0, -1), 14: (1, 0, 1),
       15: (-1, -1, 0), 16: (-1, 1, 0), 17: (1, -1, 0), 18: (1, 1, 0),
       19: (0, 0, -2), 20: (0, 0, 2), 21: (0, -2, 0), 22: (0, 2, 0),
       23: (-2, 0, 0), 24: (2, 0, 0)}
PRE = ['b', 'x0', 'TA', 'TB', 'TC', 'CJM', 'CJP', 'CKM', 'CKP']


class Op:
    def __init__(self, pattern):
        files = sorted(glob.glob(pattern))
        blocks = []
        for fn in files:
            with open(fn, 'rb') as f:
                nmb, nst, n1, n2, n3 = np.fromfile(f, np.int32, 5)
                meta = np.fromfile(f, np.int32, 5*nmb).reshape(nmb, 5)
                data = np.fromfile(f, np.float64).reshape(nmb, nst + 9, n3, n2, n1)
            for m in range(nmb):
                blocks.append((meta[m], data[m]))
        self.nst, self.n1, self.n2, self.n3 = int(nst), int(n1), int(n2), int(n3)
        nb2 = max(b[0][2] for b in blocks) + 1
        nb3 = max(b[0][3] for b in blocks) + 1
        self.nb2, self.nb3 = nb2, nb3
        N1, N2, N3 = n1, n2*nb2, n3*nb3
        self.shape = (N3, N2, N1)
        self.st = np.zeros((nst,) + self.shape)
        self.f = {k: np.zeros(self.shape) for k in PRE}
        self.blk = np.zeros(self.shape, dtype=np.int32)   # block id per cell
        for bid, (mt, d) in enumerate(blocks):
            _, lx1, lx2, lx3, _ = mt
            sl = (slice(lx3*n3, (lx3 + 1)*n3), slice(lx2*n2, (lx2 + 1)*n2), slice(None))
            self.st[(slice(None),) + sl] = d[:nst]
            for q, k in enumerate(PRE):
                self.f[k][sl] = d[nst + q]
            self.blk[sl] = bid
        self.b = self.f['b'].ravel()
        self.x0 = self.f['x0'].ravel()
        self.N = self.b.size
        self.A = self.assemble(range(nst))

    def assemble(self, slots, mask=None):
        """sparse matrix of the given stencil slots; mask(dk,dj,di) -> bool keeps a
        coupling only where it is True (arrays over cells)"""
        N3, N2, N1 = self.shape
        idx = np.arange(self.N).reshape(self.shape)
        K, J, I = np.meshgrid(np.arange(N3), np.arange(N2), np.arange(N1), indexing='ij')
        rows, cols, vals = [], [], []
        for s in slots:
            dk, dj, di = OFF[s]
            c = self.st[s]
            ii = I + di
            ok = (ii >= 0) & (ii < N1) & (c != 0.0)
            if mask is not None:
                ok &= mask(dk, dj, di)
            kk = (K + dk) % N3
            jj = (J + dj) % N2
            rows.append(idx[ok])
            cols.append(idx[kk[ok], jj[ok], ii[ok]])
            vals.append(c[ok])
        r = np.concatenate(rows)
        return sp.csr_matrix((np.concatenate(vals), (r, np.concatenate(cols))),
                             shape=(self.N, self.N))
