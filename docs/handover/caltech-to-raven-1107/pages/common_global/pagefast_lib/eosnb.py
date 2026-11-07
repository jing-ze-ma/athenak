"""pagefast-1102: numba kernel for eos_gen.Table.temperature's row bisection (the numpy version did ~10 full-array gather
passes over the table per panel).  Same arithmetic per cell: row(j) = (1 - fx) * fe[j, i] + fx * fe[j, i + 1] in float64
(separate multiply/add, no fastmath, so no FMA contraction), the same bisection from (0, ny - 1) while jhi - jlo > 1 with
'row(jm) > ft' moving jhi, else jlo.  Returns jlo, jhi, row(jlo), row(jhi), row(0).  nogil: the per-panel threads run it."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pylib'))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'numba_cache'))
import numpy as np
import numba


@numba.njit(cache=True, nogil=True, fastmath=False)
def bisect(fe, i, fx, ft, ny):
    n = i.shape[0]
    jlo = np.zeros(n, np.int64); jhi = np.empty(n, np.int64)
    flo = np.empty(n); fhi = np.empty(n); f0 = np.empty(n)
    for c in range(n):
        ic = i[c]; a = 1 - fx[c]; b = fx[c]; t = ft[c]
        lo = 0; hi = ny - 1
        while hi - lo > 1:
            jm = (lo + hi) // 2
            if a * fe[jm, ic] + b * fe[jm, ic + 1] > t:
                hi = jm
            else:
                lo = jm
        jlo[c] = lo; jhi[c] = hi
        flo[c] = a * fe[lo, ic] + b * fe[lo, ic + 1]
        fhi[c] = a * fe[hi, ic] + b * fe[hi, ic + 1]
        f0[c] = a * fe[0, ic] + b * fe[0, ic + 1]
    return jlo, jhi, flo, fhi, f0
