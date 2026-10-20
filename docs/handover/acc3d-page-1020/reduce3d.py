#!/usr/bin/env python3
"""ACC3D page 10-20: per-dump reduction for the 3-D isosurface movie (READ ONLY on the run).
usage: reduce3d.py <out_dir>  -> out_dir/red_NNNNN.npz for every hydro_w dump (static + stream) not yet reduced.
Keeps r <= d_L1: log10 rho [code] and log10 T_gas [K] (EOS-table inversion, tgas.py; float16, (phi, theta, r)),
true cell centres, time."""
import os
import sys
import numpy as np
from multiprocessing import Pool
import common as C

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)


def work(a):
    t, f, ph = a
    idx = '%05d' % C.dnum(f)
    o = os.path.join(OUT, 'red_%s.npz' % idx)
    if os.path.exists(o) and not C.FORCE:
        return None
    try:
        h = C.load(f, ['dens', 'eint'])
    except Exception as e:  # noqa: BLE001
        return 'skip %s %s' % (idx, e)
    ir = int(np.searchsorted(h['x1v'], C.DL1))
    rho = h['dens'][:, :, :ir].astype(np.float64)
    lT = C.tgas().logT(rho, h['eint'][:, :, :ir].astype(np.float64))
    lr = np.log10(np.maximum(rho, 1e-30))
    np.savez(o + '.tmp.npz', t=h['Time'], r=h['x1v'][:ir], th=h['x2v'], ph=h['x3v'],
             lrho=lr.astype(np.float16), lT=lT.astype(np.float16))
    os.replace(o + '.tmp.npz', o)
    return 'reduced %s t %.4f' % (idx, h['Time'])


if __name__ == '__main__':
    C.tgas()
    L = C.dumps('hydro_w')
    n = 0
    with Pool(C.NP) as p:
        for m in p.imap_unordered(work, L):
            if m:
                n += 1
    print('reduce3d: %d dumps, %d new' % (len(L), n))
