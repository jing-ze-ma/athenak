#!/usr/bin/env python3
"""Column gate analysis (he_star_m1 thin column, outputs every cycle, any radial grid).
usage: analyze_col.py RUNDIR STRUCT_NPZ [--brief]
Zones from the IC structure (mlt_struct.npz of make_ic mlt): rad = 0.5 R .. FeCZ bottom
(wall cell excluded), FeCZ = F_MLT/F > 1e-3, top = FeCZ top .. tau 1, thin = tau < 1,
wall = cell 0.  (i) f/(rho g) = d(rho v_r)/dt/(rho g): first step, max over per-step
values, value from the cumulative change at the end; (ii) shell L = 4 pi r^2 F_r (+ L_MLT
of the IC when mlt_flux_frozen); (iii) energy ledger; (iv) FeCZ heating measured vs the
prediction -div(A F_r)/V of the IC flux (non-frozen run: equals the missing div F_MLT)."""
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, '/viper/ptmp2/jinma/wt_hepresn/vis/python')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_gate_mlt as A  # noqa: E402

GM, LUM, RSTAR = A.GM, A.LUM, A.RSTAR


def zones_of(rc, npz):
    tau = np.interp(rc, npz['r'], npz['tau'])
    fm = np.interp(rc, npz['r'], npz['fmlt'])
    rr = npz['r'][npz['fmlt'] > 1e-3]
    b, t = rr.min(), rr.max()
    z = [('rad', (rc < b)), ('FeCZ', (rc >= b) & (rc <= t)),
         ('top-tau1', (rc > t) & (tau >= 1)), ('tau<1', tau < 1)]
    z[0][1][0] = False
    return z, b, t, tau, fm


def main():
    run, npz = sys.argv[1], np.load(sys.argv[2])
    brief = '--brief' in sys.argv
    par = A.read_mesh(run)
    frozen = False
    for ln in open(os.path.join(run, 'col.athinput')):
        if ln.split('#')[0].strip().replace(' ', '') == 'mlt_flux_frozen=true':
            frozen = True
    rl = open(os.path.join(run, 'run.log'), errors='replace').read()
    if 'mlt_flux_frozen=true' in rl.replace(' ', '') or 'mlt_flux_frozen: max' in rl:
        frozen = True
    hu = A.load(run, 'hydro_u', None)
    m1 = A.load(run, 'm1', None)
    nx = hu[0][2]['dens'].size
    re_ = A.edges(par, nx)
    rc = A.centroid(re_[:-1], re_[1:])
    g = GM / rc**2
    rho0 = hu[0][2]['dens']
    zones, rb, rt, tau, fm = zones_of(rc, npz)
    steps = [(a, b) for a, b in zip(hu[:-1], hu[1:]) if b[0] > a[0]]
    dts = np.array([b[0] - a[0] for a, b in steps])
    print('%s: nx1 %d, dx %.2e..%.2e, FeCZ %.4f-%.4f R, t_end %.1f, %d steps, '
          'dt first %.3f min %.3f max %.3f, frozen MLT %s'
          % (os.path.basename(run), nx, np.diff(re_).min(), np.diff(re_).max(),
             rb / RSTAR, rt / RSTAR, hu[-1][0], len(steps), dts[0], dts.min(),
             dts.max(), frozen))
    fs = [(b[2]['mom1'] - a[2]['mom1']) / (b[0] - a[0]) / (rho0 * g) for a, b in steps]
    fend = (hu[-1][2]['mom1'] - hu[0][2]['mom1']) / (hu[-1][0] - hu[0][0]) / (rho0 * g)
    print('(i) f/(rho g)   %-9s %10s %10s %10s  (r/R of max step value)'
          % ('zone', 'first', 'max step', 'cum end'))
    for n_, z in zones + [('wall', np.arange(nx) == 0)]:
        mx = np.array([np.abs(f[z]).max() for f in fs])
        k = int(mx.argmax())
        print('               %-9s %10.2e %10.2e %10.2e  (%.4f at t %.0f)'
              % (n_, np.abs(fs[0][z]).max(), mx.max(), np.abs(fend[z]).max(),
                 rc[z][np.abs(fs[k][z]).argmax()] / RSTAR, steps[k][1][0]))
    # (ii) shell luminosity
    Lm = LUM * np.interp(rc, npz['r'], npz['fmlt']) if frozen else 0.0
    L0 = 4 * np.pi * rc**2 * m1[0][2]['m1_f1']
    L = 4 * np.pi * rc**2 * m1[-1][2]['m1_f1']
    print('(ii) (L_rad%s)/L_star - 1 at the end: ' % (' + L_MLT' if frozen else '')
          + '  '.join('%s [%+.2e, %+.2e]' % (n_, ((L + Lm) / LUM - 1)[z].min(),
                                             ((L + Lm) / LUM - 1)[z].max())
                      for n_, z in zones)
          + ' ; max |L/L(0) - 1| %s' % '  '.join(
              '%s %.1e' % (n_, np.abs(L / L0 - 1)[z].max()) for n_, z in zones))
    hst = glob.glob(os.path.join(run, '*.user.hst'))
    if hst:
        with open(hst[0]) as fh:
            hdr = [ln for ln in fh if ln.startswith('#')]
        names = re.findall(r'\[\d+\]=(\S+)', hdr[-1])
        d = np.loadtxt(hst[0], ndmin=2)
        col = {n: d[:, i] for i, n in enumerate(names)}
        t = col['time']
        net = col['L_in'] - col['L_top']
        integ = np.concatenate([[0.0],
                                np.cumsum(0.5 * np.diff(t) * (net[1:] + net[:-1]))])
        de = col['Etot'] - col['Etot'][0]
        print('(iii) ledger (dEtot - int(L_in - L_top)dt)/(L_in t) at t %.1f: %.3e ; '
              'L_top/L_in %.6f ; M_tot change %.2e'
              % (t[-1], (de[-1] - integ[-1]) / (col['L_in'][-1] * t[-1]),
                 col['L_top'][-1] / col['L_in'][-1],
                 col['M_tot'][-1] / col['M_tot'][0] - 1))
    # (iv) FeCZ heating: e_int change rate vs -div(A F_r)/V of the IC's F_r
    vol = (re_[1:]**3 - re_[:-1]**3) / 3.0
    Fr_f = np.exp(np.interp(re_, npz['r'], np.log(npz['Fr'])))
    gpred = -(re_[1:]**2 * Fr_f[1:] - re_[:-1]**2 * Fr_f[:-1]) / vol
    if frozen:
        fm_f = np.interp(re_, npz['r'], npz['fmlt'])
        Fm_f = Fr_f * fm_f / (1 - fm_f)
        gpred = gpred - (re_[1:]**2 * Fm_f[1:]
                         - re_[:-1]**2 * Fm_f[:-1]) / vol
    u0 = hu[0][2]
    phi = GM * (1.0 / re_[0] - 1.0 / rc)
    ei0 = (u0['ener'] - 0.5 * u0['mom1']**2 / u0['dens'] - u0['dens'] * phi
           + m1[0][2]['m1_e'])
    z = zones[1][1]
    for k in sorted(set([1, len(hu) - 1])):
        t_, c_, u = hu[k]
        ei = (u['ener'] - 0.5 * u['mom1']**2 / u['dens'] - u['dens'] * phi
              + m1[k][2]['m1_e'])
        gm = (ei - ei0) / (t_ - hu[0][0])
        sel = (z & (np.abs(gpred) > 0.2 * np.abs(gpred[z]).max())
               if np.abs(gpred[z]).max() > 0 else z)
        print('(iv) t %7.1f FeCZ heating d(e_int + E)/dt: pred |max| %.3e, '
              'meas |max| %.3e erg/cm3/s; meas/pred median %.3f [%.3f, %.3f]; '
              'rel. rate max |de/e|/t %.2e /s'
              % (t_, np.abs(gpred[z]).max(), np.abs(gm[z]).max(),
                 np.median(gm[sel] / gpred[sel]) if not frozen else np.nan,
                 (gm[sel] / gpred[sel]).min() if not frozen else np.nan,
                 (gm[sel] / gpred[sel]).max() if not frozen else np.nan,
                 np.abs(gm / ei0)[z].max()))
    if not brief:
        m = re.findall(r'newton fallbacks=([0-9.e+-]+) \(([0-9.e+-]+) per', rl)
        print('(v) gas-Newton fallbacks', m[-1] if m else 'none', '; locations r/R:',
              sorted(set('%.4f' % (float(v) / RSTAR) for v in
                         re.findall(r'NEWTON-FALLBACK.* r=([0-9.e+-]+)', rl)))[:12])
        print('    FATAL/nan/nonconv lines:',
              len(re.findall(r'FATAL|nan|NaN|non-conv|nonconv', rl)))


if __name__ == '__main__':
    main()
