#!/usr/bin/env python3
"""Code-in-the-loop relaxation step for the he_star_m1 column (recipe of
bench/hestar_fecz/relax1d_3msun/analysis/make_ic.py, adapted to the sp column + M1).

usage: relax_ic_he_presn_m1.py RUNDIR PREV_IC STRUCT_NPZ OUT_IC [K]
Takes output K (default: the last) of a thin-column run (hydro_u + m1 every cycle), and
writes the next IC (he_ic_cols = 5, 7 columns  r rho eint F_r E T_rad fmlt):
  * nodes at the cell centroids (= the code's x1v, RadialCentroid) inside the mesh:
    rho, eint = ener - ke - rho Phi, E of the run (the velocity is NOT carried: the next
    segment starts from rest, which removes the breathing/kappa mode between segments);
  * outside the mesh: PREV_IC scaled by the edge ratios (keeps the ghost padding's shape);
  * F_r and fmlt are the TARGET of the IC structure (F_r + F_MLT = L/(4 pi r^2)), not the
    run's flux: they define the reference force and the frozen MLT deposit.
The pgen (he_ic_balance = true) then re-balances rho discretely at the run's T."""
import os
import sys

import numpy as np

sys.path.insert(0, '/viper/ptmp2/jinma/wt_hepresn/tests_hepresn')
import analyze_gate_mlt as A  # noqa: E402

A_RAD = 7.5657332503e-15


def main():
    run, prev, npzf, out = sys.argv[1:5]
    k = int(sys.argv[5]) if len(sys.argv) > 5 else -1
    npz = np.load(npzf)
    par = A.read_mesh(run)
    hu = A.load(run, 'hydro_u', None)
    m1 = A.load(run, 'm1', None)
    t, cyc, u = hu[k]
    E = m1[k][2]['m1_e']
    nx = u['dens'].size
    re_ = A.edges(par, nx)
    rc = A.centroid(re_[:-1], re_[1:])
    phi = A.GM * (1.0 / re_[0] - 1.0 / rc)
    rho = u['dens']
    eint = u['ener'] - 0.5 * (u['mom1']**2 + u['mom2']**2 + u['mom3']**2) / rho - rho * phi
    p = np.loadtxt(prev, comments='#')
    rp = p[:, 0]
    lo, hi = rp < rc[0], rp > rc[-1]
    def edge(col, vals, i):
        return vals[i] / np.exp(np.interp(rc[i], rp, np.log(p[:, col])))
    r = np.concatenate([rp[lo], rc, rp[hi]])
    cols = []
    for c_, v in ((1, rho), (2, eint), (4, E)):
        cols.append(np.concatenate([p[lo, c_] * edge(c_, v, 0), v, p[hi, c_] * edge(c_, v, -1)]))
    rho_n, e_n, E_n = cols
    fm = np.interp(r, npz['r'], npz['fmlt'])
    Fr = A.LUM * (1.0 - fm) / (4 * np.pi * r**2)
    Tr = (E_n / A_RAD)**0.25
    with open(out, 'w') as fh:
        fh.write('# he_star_m1 IC, relaxed: %s output %d (t = %.3f s, cycle %d), v = 0; '
                 'he_ic_cols = 5\n# r rho eint F_r(target) E T_rad fmlt(target)\n'
                 % (os.path.abspath(run), k, t, cyc))
        for row in zip(r, rho_n, e_n, Fr, E_n, Tr, fm):
            fh.write('%.12e %.12e %.12e %.12e %.12e %.10e %.8e\n' % row)
    print('wrote %s from t = %.2f; max |rho/rho_prev - 1| %.2e, |E/E_prev - 1| %.2e (mesh)'
          % (out, t, np.abs(rho / np.exp(np.interp(rc, rp, np.log(p[:, 1]))) - 1).max(),
             np.abs(E / np.exp(np.interp(rc, rp, np.log(p[:, 4]))) - 1).max()))


if __name__ == '__main__':
    main()
