#!/usr/bin/env python3
"""Build the `ic_profile` (r rho eint) for the 4 Msun presupernova He star from
column_he4_presn_sph.txt, in the TAPERED general EOS the box production runs.

This is the spherical twin of the box pipeline, which is two scripts:

  hestar_fecz/make_ic_he.py         rho := the column's rho, and eint is solved so that
                                    the RUN's own table EOS returns the column's TOTAL
                                    pressure:  p_table(rho, eint) = p_column.  What must
                                    be preserved is HYDROSTATIC BALANCE, i.e. p(r), not
                                    T(r): the column was integrated with a fixed
                                    mu = 1.3423 ideal gas + aT^4/3, which is not a state
                                    of the run's Saha/partition-function table.
  wt_he4_adi/tools/box_convection/mk_taper_ic.py
                                    re-expresses that profile in the THIN-REGION RADIATION
                                    TAPER at FIXED (rho, T):  e = e_gas(rho,T)
                                    + w(rho) a T^4, with w the smoothstep of
                                    src/utils/rad_taper.hpp between eos_rad_rho_lo and
                                    eos_rad_rho_hi.  Feeding the untapered e to the
                                    tapered EOS instead would throw the radiation energy
                                    into the gas and multiply T by ~7 at the top.

THE ONE DIFFERENCE, stated plainly: make_ic_he.py performed the p_table = p_column solve
THROUGH THE CODE (a zero-cycle probe run wrote <problem>/column_dump, and `refine` did one
secant step per node, 4-6 passes to |dp/p| < 1e-8).  Here the same solve is done
OFFLINE by bisecting on log10 T in the EOS TABLE DUMP, the same tabulated surface
the run interpolates -- but with BILINEAR interpolation, where eos_table.hpp uses a
bicubic Hermite with the tabulated derivatives.  The difference between the two
interpolants on this grid (d log rho = 0.02, d log T = 0.004) is the accuracy quoted
in eos_table.cpp's own dump comment ("bilinear interpolation on this grid is far more
accurate than the plot that consumes it"); it is measured below by re-evaluating p
with the bilinear surface and is
NOT the 1e-8 that a probe-driven secant reaches.  One `refine` pass against a real
zero-cycle dump of the spherical pgen will close that gap once the reader exists.

Usage:  python3 make_ic_sph.py            -> ic_he4_presn_sph.txt
"""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
COL = os.path.join(HERE, 'column_he4_presn_sph.txt')
DUMP = ('/viper/u2/jinma/ATHENAK/bench/analysis_0916/iso3d_w5/hestar/'
        'eos_table_box_w5.txt')
OUT = os.path.join(HERE, 'ic_he4_presn_sph.txt')

A_RAD = 7.5657332503e-15
RHO_HI = 1.93e-9        # <hydro>/eos_rad_rho_hi : w = 1 (full LTE aT^4) at and above
RHO_LO = 4.97e-10       # <hydro>/eos_rad_rho_lo : w = 0 (the two-stream owns Erad) below
# THE TEMPERATURE GATE on that taper, <hydro>/eos_rad_t_hi / eos_rad_t_lo.  w is raised to
# max(w_rho, w_T) so that a cell HOTTER than T_HI keeps the full LTE aT^4 in the EOS
# whatever its density does -- see src/utils/rad_taper.hpp::WeightGated and tests_r8.  The
# two temperatures are the run's OWN T (the table T solved below, not the column's
# ideal-gas T) at the SAME two optical depths the density window came from, tau 3 and 0.3.
# T_HI = 0 leaves the gate off and reproduces the pre-gate file exactly.
T_HI = 6.3281e4         # <hydro>/eos_rad_t_hi : T_table at tau = 3.00 (r = 2.33866e11)
T_LO = 4.5213e4         # <hydro>/eos_rad_t_lo : T_table at tau = 0.30 (r = 2.38258e11)
# the planned mesh (GAP_ANALYSIS.md, FIRST-RUN PRODUCTION GRID)
X1MIN, X1MAX = 9.4868e10, 2.4057e11


# ---- read_dump / bilin / weight are mk_taper_ic.py's, unchanged --------------------
def read_dump(fn):
    with open(fn) as fh:
        lines = [fh.readline() for _ in range(6)]
    nx, ny, xmin, dx, ymin, dy = [float(v) for v in lines[2].split()[1:]]
    nx, ny = int(nx), int(ny)
    d = np.loadtxt(fn, comments='#')
    le = d[:, 0].reshape(ny, nx)
    lp = d[:, 1].reshape(ny, nx)
    x = xmin + dx*np.arange(nx)
    y = ymin + dy*np.arange(ny)
    return x, y, le, lp


def bilin(x, y, f, xq, yq):
    ix = np.clip(np.searchsorted(x, xq) - 1, 0, len(x) - 2)
    iy = np.clip(np.searchsorted(y, yq) - 1, 0, len(y) - 2)
    tx = (xq - x[ix])/(x[ix+1] - x[ix])
    ty = (yq - y[iy])/(y[iy+1] - y[iy])
    return ((1-tx)*(1-ty)*f[iy, ix] + tx*(1-ty)*f[iy, ix+1]
            + (1-tx)*ty*f[iy+1, ix] + tx*ty*f[iy+1, ix+1])


def weight(lrho, llo, lhi):
    s = np.clip((lrho - llo)/(lhi - llo), 0.0, 1.0)
    return s*s*(3.0 - 2.0*s)


def main():
    x, y, le, lp = read_dump(DUMP)
    d = np.loadtxt(COL, comments='#')
    r, p, T, rho = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    tau = d[:, 8]
    o = np.argsort(r)
    r, p, T, rho, tau = r[o], p[o], T[o], rho[o], tau[o]
    lr = np.log10(rho)

    print('column %s: %d nodes, r = %.6e .. %.6e cm' % (os.path.basename(COL), len(r),
                                                        r[0], r[-1]))
    print('  rho %.4e .. %.4e   T %.4e .. %.4e' % (rho.min(), rho.max(), T.min(),
                                                   T.max()))
    print('  EOS table dump covers log10 rho %.2f .. %.2f, log10 T %.3f .. %.3f'
          % (x[0], x[-1], y[0], y[-1]))
    inside = ((lr > x[0]) & (lr < x[-1]) & (np.log10(T) > y[0])
              & (np.log10(T) < y[-1]))
    print('  column inside the table: %s (%d of %d nodes)'
          % (bool(inside.all()), inside.sum(), len(r)))
    print('  planned mesh x1 = %.5e .. %.5e ; margins below/above: %.4e / %.4e cm'
          % (X1MIN, X1MAX, X1MIN - r[0], r[-1] - X1MAX))

    # ---- solve p_table_untapered(rho, T) = p_column for T, by bisection in log10 T
    def ptab(lt):
        return rho*10.0**bilin(x, y, lp, lr, lt) + A_RAD*10.0**(4.0*lt)/3.0
    lo = np.full_like(lr, y[0] + 1e-9)
    hi = np.full_like(lr, y[-1] - 1e-9)
    assert np.all(ptab(lo) < p) and np.all(ptab(hi) > p), 'p_column outside the table'
    for _ in range(200):
        mid = 0.5*(lo + hi)
        f = ptab(mid) - p
        lo = np.where(f < 0.0, mid, lo)
        hi = np.where(f < 0.0, hi, mid)
    lt = 0.5*(lo + hi)
    Tt = 10.0**lt
    rel = np.abs(ptab(lt)/p - 1.0)
    print('\n  bisection residual |p_bilin(rho,T)/p_col - 1|: max %.3e  rms %.3e'
          % (rel.max(), np.sqrt((rel**2).mean())))
    # how much of a difference the interpolant itself makes: the same node values read
    # with a bicubic spline instead of bilinear.  This is an ESTIMATE of the gap to the
    # run's Hermite interpolation, which uses tabulated derivatives the dump omits.
    try:
        from scipy.interpolate import RectBivariateSpline
        sp = RectBivariateSpline(y, x, lp, kx=3, ky=3)
        pc = rho*10.0**sp.ev(lt, lr) + A_RAD*10.0**(4.0*lt)/3.0
        se = RectBivariateSpline(y, x, le, kx=3, ky=3)
        dg = np.abs(10.0**se.ev(lt, lr)/10.0**bilin(x, y, le, lr, lt) - 1.0)
        pdiff = np.abs(pc/p - 1.0).max()
        print('  bilinear vs bicubic on the SAME dump: |dp/p| max %.3e ; |de_gas/e_gas|'
              ' max %.3e  <- the real accuracy of this offline solve'
              % (pdiff, dg.max()))
    except ImportError:
        pdiff = float('nan')
        print('  (scipy unavailable: interpolant comparison skipped)')
    rt = Tt/T
    print('  T_table/T_column: min %.4f  max %.4f  (the mu = 1.3423 ideal gas of the'
          ' column vs the table\'s partition functions)' % (rt.min(), rt.max()))

    eg = rho*10.0**bilin(x, y, le, lr, lt)
    er = A_RAD*Tt**4
    wrho = weight(lr, np.log10(RHO_LO), np.log10(RHO_HI))
    if T_HI > 0.0:
        wt = weight(lt, np.log10(T_LO), np.log10(T_HI))
        w = np.maximum(wrho, wt)
    else:
        wt = np.zeros_like(wrho)
        w = wrho
    eint = eg + w*er
    # the two thresholds the gate WANTS, in the run's own temperature: T_table at the two
    # optical depths the density window came from.  Printed so that the input file's
    # eos_rad_t_hi / eos_rad_t_lo can be set to them (and checked against them).
    for tt in (3.0, 0.30):
        i = int(np.argmin(np.abs(tau - tt)))
        print('  tau = %4.2f at r = %.6e (r/R %.4f): rho = %.4e  T_col = %.5e  '
              'T_table = %.5e  w_rho = %.4f  w_T = %.4f'
              % (tt, r[i], r[i]/2.3717e11, rho[i], T[i], Tt[i], wrho[i], wt[i]))
    print('  gate: w != w_rho at %d of %d nodes; max |w - w_rho| = %.4f'
          % (int((w != wrho).sum()), len(w), float(np.abs(w - wrho).max())))
    print('  taper: w < 1 above r = %.6e cm (rho < %.3e); w = 0 above r = %.6e cm'
          % (r[w < 1.0 - 1e-12].min() if np.any(w < 1.0 - 1e-12) else np.nan, RHO_HI,
             r[w <= 0.0].min() if np.any(w <= 0.0) else np.nan))
    print('  e_rad/e_gas at the top node: %.3f ; eint(tapered)/eint(untapered) there:'
          ' %.4f' % (er[-1]/eg[-1], eint[-1]/(eg[-1] + er[-1])))
    print('  eint monotone inward: %s' % bool(np.all(np.diff(eint) < 0.0)))

    print('\n  %12s %8s %12s %12s %12s %8s %12s %12s'
          % ('r[cm]', 'r/R', 'rho', 'T_col', 'T_table', 'w', 'e_untapered', 'eint'))
    for i in list(range(0, len(r), max(1, len(r)//16))) + [len(r) - 1]:
        print('  %12.5e %8.4f %12.4e %12.4e %12.4e %8.4f %12.4e %12.4e'
              % (r[i], r[i]/2.3717e11, rho[i], T[i], Tt[i], w[i], eg[i] + er[i],
                 eint[i]))

    hdr = (
        '# hestar_presn: SPHERICAL ic_profile for the 4.0 Msun He star presupernova\n'
        '# (Woosley 2019), built from column_he4_presn_sph.txt by make_ic_sph.py.\n'
        '# THREE COLUMNS, ascending r:   r[cm]  rho[g/cm^3]  eint[erg/cm^3]\n'
        '# rho is the column unchanged.  eint is the PRIMITIVE internal energy density\n'
        '# w0(IEN) of AthenaK\'s general EOS with <hydro>/eos_radiation = true and the\n'
        '# thin-region radiation taper <hydro>/eos_rad_rho_hi = %.3g,\n'
        '# eos_rad_rho_lo = %.3g (src/utils/rad_taper.hpp):\n'
        '#     eint = e_gas(rho,T) + w a T^4,  w = max(w_rho, w_T) with the TEMPERATURE\n'
        '# GATE <hydro>/eos_rad_t_hi = %.5g, eos_rad_t_lo = %.5g K: w_rho is a\n'
        '# smoothstep in log10 rho, w_T one in log10 T, and the max keeps a HOT cell in\n'
        '# LTE whatever\n'
        '# its density does (rad_taper::WeightGated; tests_r8).\n'
        '# T at each node is NOT the column\'s T: it is solved so that the run\'s own\n'
        '# table returns the column\'s TOTAL pressure with the UNTAPERED radiation,\n'
        '#     p_gas_table(rho,T) + a T^4/3 = p_column(r),\n'
        '# i.e. hydrostatic balance is what is preserved, exactly as in\n'
        '# hestar_fecz/make_ic_he.py; the taper is then applied at fixed (rho,T), as in\n'
        '# tools/box_convection/mk_taper_ic.py.  The (1-w) a T^4/3 that leaves p at the\n'
        '# top is carried instead by the two-stream force <problem>/rt_rad_force=true.\n'
        '# EOS table dump: %s\n'
        '# Solve done OFFLINE with BILINEAR interpolation of that dump; the box line\n'
        '# refined through a zero-cycle probe run; the bisection residual on the\n'
        '# bilinear surface is %.1e, but bilinear vs bicubic on the same nodes differs\n'
        '# by up to %.1e in p, which is the honest accuracy of this file.\n'
        '# Star: M = 6.2651e33 g, L = 2.3066e38 erg/s, Teff = 48978 K, R = 2.3717e11cm,\n'
        '# X = 0, Y = 0.98, Z = 0.02.  Covers the planned mesh %.5e .. %.5e cm.\n'
        '# r[cm]  rho[g/cm^3]  eint[erg/cm^3]\n'
        % (RHO_HI, RHO_LO, T_HI, T_LO, DUMP, rel.max(), pdiff, X1MIN, X1MAX))
    with open(OUT, 'w') as fh:
        fh.write(hdr)
        for a, b, cc in zip(r, rho, eint):
            fh.write('%.10e %.10e %.10e\n' % (a, b, cc))
    print('\nwrote %s: %d nodes, r = %.6e .. %.6e cm' % (OUT, len(r), r[0], r[-1]))


if __name__ == '__main__':
    main()
