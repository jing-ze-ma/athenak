#!/usr/bin/env python3
"""Turn a TOPS gray-mean dump (fetch_tops.py output: T[keV] rho kR kP, one
row per point) into the ASCII table AthenaK reads -- the same layout as
data/stellar_opac/rosseland_gs98_x0.7_z0.014.txt:

    # <free comment lines>
    # nT nD lTmin dlT lDmin dlD      <- FIRST numeric comment line, the grid
    # <more comment lines>
    nT*nD values of log10 kappa [cm^2/g], one per line, T SLOWEST

Axes are log10 T[K] (uniform, nT nodes from lTmin step dlT) and log10 rho
[g/cm^3] (uniform, nD nodes from lDmin step dlD); the run must then set
<hydro>/rad_kappa_src = table_rho.  Default grid = the repo's:
217 x 281, log T 2.6..8.0 step 0.025, log rho -14..0 step 0.05.

The TOPS source grid is irregular in log T (the database temperature list) and
regular in log rho, so the conversion is a bilinear interpolation in
(log10 T, log10 rho) of log10 kappa, with EDGE-FILL outside the source box
(exactly what merge_rosseland.py does for OPLIB).

usage: convert_tops.py IN.dat OUT.txt {planck|rosseland} ["header comment"]
"""
import sys
import numpy as np

KEV_K = 1.1604518e7      # 1 keV in kelvin

LT0, DLT, NT = 2.6, 0.025, 217
LD0, DLD, ND = -14.0, 0.05, 281


def read_tops(fn):
    d = np.loadtxt(fn)
    T = np.unique(d[:, 0])
    R = np.unique(d[:, 1])
    kR = np.full((len(T), len(R)), np.nan)
    kP = np.full((len(T), len(R)), np.nan)
    ti = {v: i for i, v in enumerate(T)}
    ri = {v: i for i, v in enumerate(R)}
    for t, r, a, b in d:
        kR[ti[t], ri[r]] = a
        kP[ti[t], ri[r]] = b
    assert np.isfinite(kR).all() and np.isfinite(kP).all(), 'ragged TOPS grid'
    return np.log10(T*KEV_K), np.log10(R), np.log10(kR), np.log10(kP)


def bilin(xs, ys, F, X, Y):
    """bilinear in (xs, ys) with edge-fill; xs, ys ascending, may be uneven."""
    X = np.clip(X, xs[0], xs[-1])
    Y = np.clip(Y, ys[0], ys[-1])
    i = np.clip(np.searchsorted(xs, X) - 1, 0, len(xs)-2)
    j = np.clip(np.searchsorted(ys, Y) - 1, 0, len(ys)-2)
    fx = (X - xs[i])/(xs[i+1] - xs[i])
    fy = (Y - ys[j])/(ys[j+1] - ys[j])
    return ((1-fx)*(1-fy)*F[i, j] + fx*(1-fy)*F[i+1, j]
            + (1-fx)*fy*F[i, j+1] + fx*fy*F[i+1, j+1])


def write_table(fn, K, comments):
    lT = LT0 + DLT*np.arange(NT)
    lD = LD0 + DLD*np.arange(ND)
    with open(fn, 'w') as f:
        for c in comments:
            f.write('# %s\n' % c)
        f.write('# grid: nT nD lTmin dlT lDmin dlD\n')
        f.write('# %d %d %g %g %g %g\n' % (NT, ND, LT0, DLT, LD0, DLD))
        f.write('# then nT*nD rows, T slowest: log10 kappa [cm^2/g]\n')
        for i in range(NT):
            for j in range(ND):
                f.write('%.5f\n' % K[i, j])
    return lT, lD


if __name__ == '__main__':
    src, out, which = sys.argv[1], sys.argv[2], sys.argv[3]
    extra = sys.argv[4] if len(sys.argv) > 4 else ''
    lT_s, lD_s, lkR, lkP = read_tops(src)
    F = lkP if which.startswith('p') else lkR
    lT = LT0 + DLT*np.arange(NT)
    lD = LD0 + DLD*np.arange(ND)
    TT, DD = np.meshgrid(lT, lD, indexing='ij')
    K = bilin(lT_s, lD_s, F, TT, DD)
    name = 'Planck' if which.startswith('p') else 'Rosseland'
    write_table(out, K, [
        'AthenaK stellar %s opacity table, log10 kappa_%s [cm^2/g] on '
        '(log10 T[K], log10 rho[g/cm3])' % (name, name[0]),
        'source: LANL TOPS/ATOMIC gray means, fetched with fetch_tops.py from %s'
        % src.split('/')[-1],
        extra,
        'VALID BOX: log T %.3f..%.3f, log rho %.2f..%.2f -- outside it the values '
        'are EDGE-FILLED and are NOT data.' % (lT_s[0], lT_s[-1], lD_s[0],
                                               lD_s[-1]),
        'The Planck mean EXCLUDES scattering; the Rosseland mean includes it.',
        'Read with <hydro>/rad_kappa_src = table_rho.'])
    bad = (TT < lT_s[0]) | (TT > lT_s[-1]) | (DD < lD_s[0]) | (DD > lD_s[-1])
    print('wrote %s ; %.1f %% of nodes edge-filled' % (out, 100*bad.mean()))
