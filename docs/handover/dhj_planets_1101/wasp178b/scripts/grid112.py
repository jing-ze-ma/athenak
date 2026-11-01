"""grid112.py (dhj_planets_1101): radial grid of the WASP-121b hydro fresh deck class
(C32 x nx1 112, x1max 2.1e10, 8 coefficients + 2 sech^2 bumps) carried to another planet.

Reference = w121hfresh_1101/deck/hfresh.athinput on its own fresh-start column
(ic_w121_3x_f5.txt, ap 1.124210e10, grav 1172.7924, point mass; the 3x EOS dump, whose
  table starts at log rho -14: the column is tabulated to 1e-9 bar and continued
  isothermally above, as planet_setup.p_of_r does).
  1. cells per pressure scale height n_H(p) = H/dr of every W121 cell, H = (p/rho)/g(r);
  2. p_top = p(x1max = 2.1e10) on that column (isothermal continuation above the IC top);
  3. new planet: x1max = r(p_top) on ITS IC column (the same pressure, i.e. the same
     number of top scale heights above 1e-6 bar), target dr_t(r) = H(r)/n_H(p(r));
  4. minimax (as w121_mhd_1003/gridmm.py) of max(dr/dr_t) over the 8 coefficients + 2
     bumps for NX1 cells; the smallest NX1 with max ratio <= 1 hits every W121 n_H.
usage: [X1MAX=r] grid112.py PROFILE GRAV AP NX1 [NX1 ...]  (writes gridkeys_<N>.txt, faces_<N>.txt)
  X1MAX (dhj_planets_1101 addition): explicit x1max instead of r(W121 p_top); cells/H targets
  above the W121 top pressure are held at the W121 top value.
"""
import os
import sys
import numpy as np
from scipy.optimize import minimize
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import planet_setup as ps  # noqa: E402
import dhj_remap as dr_  # noqa: E402

EOSF = '/viper/ptmp2/jinma/w121prod_0929/met3/dump/eos_table.txt'   # W121 reference (3x)
# dhj_planets_1101: the planet column's EOS dump dir (env PLANET_EOSD; default the 3x dump)
EOSD_P = os.environ.get('PLANET_EOSD', '/viper/ptmp2/jinma/w121prod_0929/met3/dump')
EOSF_P = EOSD_P + '/eos_table.txt'
W_IC = '/viper/ptmp2/jinma/w121_s10fresh_1027/pkg/ic/ic_w121_3x_f5.txt'
W_AP, W_G, W_X1MAX, W_N = 1.124210e10, 1172.7924, 2.1e10, 112
W_C = [-1.509683, -1.314625, -13.749263, 79.204795, -274.939323, 523.658095, -465.600806,
       153.292986, -1.290044, 0.520389, 0.146721, -2.426785, 0.689686, 0.104620, 0, 0, 0, 0]
BAR = 1e6


def cells(col, gm, ap, e):
    rc = dr_.centroids(e)
    dr = np.diff(e)
    p = ps.p_of_r(col, gm, ap, rc)
    H = np.interp(np.log10(p*BAR), col.lp, col.pr)/(gm/rc**2)
    return rc, dr, p, H


def r_of_p_ext(col, gm, ap, pbar):
    """r(p) on the column; above its top the isothermal continuation (inverse of
    planet_setup.p_of_r there; planet_setup.r_of_p clamps at the column top)."""
    lp = np.log10(pbar*BAR)
    phi = np.where(lp < col.lp[0], col.phi[0] + col.pr[0]*np.log(10.0)*(col.lp[0] - lp),
                   np.interp(lp, col.lp, col.phi))
    return 1.0/(1.0/ap - phi/gm)


def table(p, nh, dr, label):
    out = ['%s: cells per H (min / median) and dr range per pressure band' % label]
    for lo, hi in ((1e-9, 1e-6), (1e-6, 1e-4), (1e-4, 1e-2), (1e-2, 1), (1, 10), (10, 300)):
        m = (p >= lo) & (p < hi)
        if m.any():
            out.append('  %7.0e..%-7.0e bar: n %3d  cells/H min %.2f med %.2f  dr %.3e..%.3e'
                       % (lo, hi, m.sum(), nh[m].min(), np.median(nh[m]), dr[m].min(),
                          dr[m].max()))
    m = p >= 1e-6
    out.append('  all p >= 1e-6 bar: cells/H min %.2f;  cells above 1e-6 bar: %d'
               % (nh[m].min(), (~m).sum()))
    return '\n'.join(out)


def main():
    prof, grav, ap = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    ns = [int(v) for v in sys.argv[4:]]
    eos = ps.EOSTable(EOSF)
    colw = ps.Column(W_IC, eos, ptop=1e-9)
    gmw = W_G*W_AP**2
    ew = dr_.edges(W_N, W_AP, W_X1MAX, np.array(W_C, float))
    rcw, drw, pw, Hw = cells(colw, gmw, W_AP, ew)
    nhw = Hw/drw
    ptop = ps.p_of_r(colw, gmw, W_AP, W_X1MAX)
    print(table(pw, nhw, drw, 'WASP-121b hfresh (nx1 112, x1max 2.1e10)'))
    print('W121 p(x1max) on the IC column = %.3e bar; r(1e-6 bar) = %.5e; check r(p(x1max)) = '
          '%.6e' % (ptop, ps.r_of_p(colw, gmw, W_AP, 1e-6), r_of_p_ext(colw, gmw, W_AP, ptop)))
    lpw = np.log10(pw)[::-1]
    nhw_r = nhw[::-1]

    def nh_t(p):     # W121 cells/H at pressure p (ends held)
        return np.interp(np.log10(p), lpw, nhw_r)

    eosp = ps.EOSTable(EOSF_P) if EOSF_P != EOSF else eos
    print('planet EOS dump: %s' % EOSF_P)
    col = ps.Column(prof, eosp, ptop=1e-9)
    gm = grav*ap**2
    x1max = float(r_of_p_ext(col, gm, ap, ptop))
    if os.environ.get('X1MAX'):    # dhj_planets_1101: explicit top (puffy planets: Hill cap)
        print('same-pressure rule would give x1max = %.6e; env X1MAX overrides' % x1max)
        x1max = float(os.environ['X1MAX'])
        ptop = float(ps.p_of_r(col, gm, ap, x1max))
    print('planet: ap %.6e grav %.4f  x1max = r(%.3e bar) = %.6e  (x1max/ap %.4f, '
          'W121 %.4f)' % (ap, grav, ptop, x1max, x1max/ap, W_X1MAX/W_AP))
    rr = np.linspace(ap, x1max, 4001)
    pr_ = ps.p_of_r(col, gm, ap, rr)
    Hr = np.interp(np.log10(pr_*BAR), col.lp, col.pr)/(gm/rr**2)
    drt_r = Hr/nh_t(pr_)

    def drt(r):
        return np.interp(r, rr, drt_r)

    # cumulative target cell count -> least-squares start (as grid.py)
    ncum = np.concatenate([[0], np.cumsum(np.diff(rr)/(0.5*(drt_r[1:] + drt_r[:-1])))])
    print('target cell count (exact W121 cells/H): %.1f' % ncum[-1])
    for N in ns:
        xi = np.arange(N + 1)/N
        rt = np.interp(xi*ncum[-1], ncum, rr)
        ut = (rt - ap)/(x1max - ap)
        A = np.array([xi**k*(1 - xi) for k in range(1, 9)]).T
        w = 1.0/np.gradient(ut)
        c, *_ = np.linalg.lstsq(A*w[:, None], (ut - xi)*w, rcond=None)

        def ev(p):
            cc = np.concatenate([p[:8], p[8:14], np.zeros(4)])
            e = dr_.edges(N, ap, x1max, cc)
            d = np.diff(e)
            if np.any(d <= 0):
                return None
            return e, d

        def obj(p, beta=60.0):
            r = ev(p)
            if r is None:
                return 1e3
            e, d = r
            lr = np.log(d/drt(0.5*(e[1:] + e[:-1])))
            m = lr.max()
            return m + np.log(np.exp(beta*(lr - m)).sum())/beta
        best = None
        starts = [np.array(W_C[:14], float)]       # the W121 hfresh keys themselves
        for b1x, b2x in ((0.52, 0.69), (0.4, 0.62), (0.6, 0.95), (0.55, 0.9), (0.65, 0.85)):
            starts.append(np.concatenate([c, [0.0, b1x, 0.03, 0.0, b2x, 0.03]]))
        for p0 in starts:
            res = minimize(obj, p0, method='Nelder-Mead',
                           options=dict(maxiter=40000, maxfev=40000, xatol=1e-7, fatol=1e-9))
            res = minimize(obj, res.x, method='Powell', options=dict(maxiter=20000))
            e, d = ev(res.x)
            rmax = (d/drt(0.5*(e[1:] + e[:-1]))).max()
            if best is None or rmax < best[0]:
                best = (rmax, res.x)
        rmax, p = best
        p = np.round(p, 6)
        p[[10, 13]] = np.abs(p[[10, 13]])
        e, d = ev(p)
        rmax = (d/drt(0.5*(e[1:] + e[:-1]))).max()
        rc, drn, pn, Hn = cells(col, gm, ap, e)
        print('\nN %d: max dr/target %.4f (6-decimal keys), dr min %.3e max %.3e'
              % (N, rmax, d.min(), d.max()))
        print(table(pn, Hn/drn, drn, 'planet nx1 %d' % N))
        names = ['f_stretch_r_c%d' % k for k in range(1, 9)] + \
            ['f_stretch_r_b1_amp', 'f_stretch_r_b1_x', 'f_stretch_r_b1_w',
             'f_stretch_r_b2_amp', 'f_stretch_r_b2_x', 'f_stretch_r_b2_w']
        keys = 'mesh/nx1=%d meshblock/nx1=%d mesh/x1min=%.6e mesh/x1max=%.6e ' % (
            N, N, ap, x1max) + ' '.join('mesh/%s=%.6f' % (n, v) for n, v in zip(names, p))
        print('KEYS ' + keys)
        open(os.path.join(os.getcwd(), 'gridkeys_%d.txt' % N), 'w').write(keys + '\n')
        np.savetxt(os.path.join(os.getcwd(), 'faces_%d.txt' % N), e, fmt='%.8e')
        # radial dt estimate (CFL 0.3, c_s + 1 km/s), both planets
        cs = np.interp(np.log10(pn*BAR), col.lp, col.cs)
        csw = np.interp(np.log10(pw*BAR), colw.lp, colw.cs)
        dtn = 0.3*np.min(drn/(cs + 1e5))
        dtw = 0.3*np.min(drw/(csw + 1e5))
        print('radial dt est (CFL 0.3, IC column, c_s + 1 km/s): planet %.2f s, W121 %.2f s,'
              ' ratio %.3f' % (dtn, dtw, dtn/dtw))


if __name__ == '__main__':
    main()
