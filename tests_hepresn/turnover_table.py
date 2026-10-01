#!/usr/bin/env python3
"""Per-dump table of the he_star_m1 3-D wedge (porosity: steady or runaway?).
usage: turnover_table.py RUNDIR INPUT [--every N] [--first N] [--tt 4700] [--nums a b ..]
Uses the hydro_w + m1 bin pairs (every N-th, default 4 = 1 turnover at dt 1175 s) and the
user.hst + run*.log of RUNDIR (gates.py readers). Columns:
  t/tt; hst over the preceding window: L_top/L_in, M_tot/M0, Mdot_top (g/s, wedge), <dt>,
  counts of Picard NON-CONVERGED (NC) and NEWTON-FALLBACK (NFB); from the dump:
  dt proxy (min dx/(|v|+c_gas), s) and its r/R + direction, drho = rms(rho')/<rho> median
  in the zones rad (r/R < 0.66), FeCZ (0.66-0.91), top (0.91 - r_ph), porosity factor
  P = <F_1>/F_diff(<rho>,<T>) FeCZ median, L_rad/L and L_conv/L = (Lcg + Lcr)/L FeCZ
  median, v_r'/v_MLT FeCZ median, Gamma of the shell mean state
  (kappa(<rho>,<T>) L / (4 pi c G M)) max and the outermost r/R with Gamma > 1,
  photosphere r/R (tau = 1 of the shell-mean state from the outer edge), porous top
  (outermost r/R < r_ph with drho > 0.1), hk4 = rho' power fraction at <= 4 cells of
  this grid,
  hk64 = at |k| >= 16 (= <= 4 cells of the 64 x 64 grid), FeCZ medians."""
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, '/viper/ptmp2/jinma/wt_hepresn/vis/python')
sys.path.insert(0, '/viper/ptmp2/jinma/wt_hepresn/tests_hepresn')
from bin_convert import read_binary  # noqa: E402
import analyze_gate_mlt as A  # noqa: E402
import analyze_3d as B  # noqa: E402
import gates as G  # noqa: E402
import porosity as P  # noqa: E402

C = 2.99792458e10
L = A.LUM
MIC = '/viper/ptmp2/jinma/hepresn_0929/mlt/mlt_struct.npz'


def arg(name, default):
    return type(default)(sys.argv[sys.argv.index(name) + 1]) \
        if name in sys.argv else default


def hst_window(h, ev, t0, t1):
    s = (h['time'] > t0) & (h['time'] <= t1)
    if not s.any():
        return [np.nan] * 5 + [0, 0]
    e = np.array([q for q in ev if t0 < q[0] <= t1] or [(0, -1)])
    return [(h['L_top'][s] / h['L_in'][s]).mean(), h['M_tot'][s][-1] / h['M_tot'][0],
            h['Mdot_top'][s].mean(), h['dt'][s].mean(), 0,
            int((e[:, 1] == 1).sum()), int((e[:, 1] == 0).sum())]


def shell_fft(rho, rm, geo):
    nmb, n3, n2, n1 = rho.shape
    t0 = sorted(set(np.round(geo[:, 2], 8)))
    p0 = sorted(set(np.round(geo[:, 4], 8)))
    sh = np.empty((len(p0) * n3, len(t0) * n2, n1))
    for m in range(nmb):
        it = t0.index(np.round(geo[m, 2], 8))
        ip = p0.index(np.round(geo[m, 4], 8))
        sh[ip * n3:(ip + 1) * n3, it * n2:(it + 1) * n2, :] = rho[m] / rm
    N3, N2 = sh.shape[:2]
    pw = np.abs(np.fft.fft2(sh - sh.mean(axis=(0, 1)), axes=(0, 1)))**2
    k3 = np.abs(np.fft.fftfreq(N3) * N3)[:, None]
    k2 = np.abs(np.fft.fftfreq(N2) * N2)[None, :]
    tot = pw.sum(axis=(0, 1)) + 1e-300
    hi4 = np.maximum(k3 / N3, k2 / N2) >= 0.25
    hi64 = np.maximum(k3, k2) >= 16
    return ((pw * hi4[:, :, None]).sum(axis=(0, 1)) / tot,
            (pw * hi64[:, :, None]).sum(axis=(0, 1)) / tot)


def main():
    run, inp = sys.argv[1], sys.argv[2]
    tt = arg('--tt', 4700.0)
    every = arg('--every', 4)
    first = arg('--first', 0)
    par = B.read_mesh(inp)
    re_ = A.edges(par, int(par['nx1']))
    rc = A.centroid(re_[:-1], re_[1:])
    x = rc / A.RSTAR
    kap = P.kap_table()
    slp = B.eos_pg()
    npz = np.load(MIC)
    h = G.read_hst(os.path.join(run, 'hepresn.user.hst'))
    # restart links append: keep the last record of any time (later runs override)
    _, keep = np.unique(h['time'][::-1], return_index=True)
    keep = len(h['time']) - 1 - keep
    h = {k: v[keep] for k, v in h.items()}
    ev = []
    for fn in sorted(glob.glob(os.path.join(run, 'run*.log')), key=os.path.getmtime):
        ev += G.read_log(fn)
    fh = sorted(glob.glob(os.path.join(run, 'bin', '*.hydro_w.*.bin')))
    nums = [int(re.findall(r'\.(\d+)\.bin', f)[0]) for f in fh]
    if '--nums' in sys.argv:
        sel = [int(a) for a in sys.argv[sys.argv.index('--nums') + 1:]
               if not a.startswith('--')]
    else:
        sel = [n for n in nums if n >= first and (n - first) % every == 0]
    zr, zf = x < 0.66, (x > 0.66) & (x < 0.91)
    print('t/tt Ltop/L M/M0 Mdot_top <dt> NC NFB | dtprox@r/R(dir) | drho rad/FeCZ/top '
          '| P Lrad/L Lconv/L vr/vMLT | Gmax r(G>1) r_ph r_por | hk4 hk64')
    tprev = None
    for n in sel:
        fhn = os.path.join(run, 'bin', 'hepresn.hydro_w.%05d.bin' % n)
        fmn = os.path.join(run, 'bin', 'hepresn.m1.%05d.bin' % n)
        o = B.shells(fhn, fmn, re_, npz, slp)
        dh, dm = read_binary(fhn), read_binary(fmn)
        rho = np.asarray(dh['mb_data']['dens'])
        E = np.asarray(dm['mb_data']['m1_e'])
        F1 = np.asarray(dm['mb_data']['m1_f1'])
        geo = np.asarray(dh['mb_geometry'])
        nmb, n3, n2, n1 = rho.shape
        w = np.empty((nmb, n3, n2, 1))
        for m in range(nmb):
            t = np.linspace(geo[m, 2], geo[m, 3], n2 + 1)
            w[m] = (np.cos(t[:-1]) - np.cos(t[1:]))[None, :, None]
        ws = w.sum()

        def mean(q):
            return (q * w).sum(axis=(0, 1, 2)) / ws
        rm, Em = mean(rho), mean(E)
        Tm = (Em / P.AR)**0.25
        km = kap(Tm, rm)
        Fmean = -C / (3 * km * rm) * np.gradient(Em, rc)
        drho = np.sqrt(mean((rho - rm)**2)) / rm
        por = mean(F1) / Fmean
        gam = km * L / (4 * np.pi * C * A.GM)
        dtau = km * rm * np.diff(re_)
        tau = np.cumsum(dtau[::-1])[::-1]
        rph = x[np.where(tau >= 1.0)[0].max()] if (tau >= 1).any() else np.nan
        zt = (x > 0.91) & (x < rph)
        g1 = np.where(gam > 1)[0]
        rpo = np.where((drho > 0.1) & (x < rph))[0]
        hk4, hk64 = shell_fft(rho, rm, geo)
        t = float(dh['time'])
        t0 = tprev if tprev is not None else t - every * 1175.0
        tprev = t
        hw = hst_window(h, ev, t0, t)
        print('%5.2f %6.3f %6.4f %8.2e %6.2f %4d %5d | %6.2f@%.3f(%d) | %.3f/%.3f/%.3f '
              '| %5.2f %5.3f %+6.3f %.2e | %5.2f %.3f %.3f %.3f | %.3f %.3f'
              % (t / tt, hw[0], hw[1], hw[2], hw[3], hw[5], hw[6], o['dtmin'],
                 o['dtloc'][1], o['dtloc'][0], np.median(drho[zr]), np.median(drho[zf]),
                 np.median(drho[zt]) if zt.any() else np.nan, np.median(por[zf]),
                 np.median(o['Lrad'][zf]) / L,
                 np.median((o['Lcg'] + o['Lcr'])[zf]) / L,
                 np.median(o['vrms'][zf] / o['vmlt'][zf]), gam[zf].max(),
                 x[g1.max()] if len(g1) else np.nan, rph,
                 x[rpo.max()] if len(rpo) else np.nan, np.median(hk4[zf]),
                 np.median(hk64[zf])), flush=True)


if __name__ == '__main__':
    main()
