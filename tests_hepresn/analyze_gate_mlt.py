#!/usr/bin/env python3
"""Gate 1a analysis for the STRETCHED-grid thin column (mlt / hopf ICs).
usage: analyze_gate_mlt.py RUNDIR [STRUCT_NPZ]
  RUNDIR      run dir with col.athinput (stretch keys read from it), hepresn.*.hst and
              bin/hepresn.hydro_u.*.bin, hepresn.m1.*.bin (one output per cycle)
  STRUCT_NPZ  mlt_struct.npz of make_ic_he_presn_m1.py mlt (r, tau, fmlt, ...); used for the
              zone edges (tau = 1) and for the predicted FeCZ heating -div F_r of the IC
Prints: (i) net radial force f/(rho g) = d(rho v_r)/dt/(rho g) by zone: first step | max over
per-step values | cumulative value at the end; (ii) shell luminosity per zone; (iii) the energy
ledger; (iv) gas heating rate of the FeCZ against the prediction from the IC.
Zones (r/R): rad 0.5-0.635 (wall cell 0 excluded), FeCZ 0.635-0.968, 0.968 - tau 1, tau < 1,
wall cell 0 separately.
"""
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, '/viper/ptmp2/jinma/wt_hepresn/vis/python')
from bin_convert import read_binary  # noqa: E402

GM = 4.18143e26
LUM = 2.3066e38
RSTAR = 2.3717e11


def read_mesh(run):
    par = {}
    sec = None
    for ln in open(os.path.join(run, 'col.athinput')):
        ln = ln.split('#')[0].strip()
        if ln.startswith('<'):
            sec = ln.strip('<>')
        elif '=' in ln and sec == 'mesh':
            k, v = [x.strip() for x in ln.split('=', 1)]
            par[k] = v
    return par


def edges(par, nx):
    r0, r1 = float(par['x1min']), float(par['x1max'])
    xi = np.arange(nx + 1) / nx
    u = xi.copy()
    if par.get('use_grid_stretch_r_poly', 'false') == 'true':
        for k in range(1, 9):
            c = float(par.get('f_stretch_r_c%d' % k, 0.0))
            u += c * xi**k * (1 - xi)
        for b in (1, 2):   # local bumps (coordinates/grid_stretch.hpp)
            a = float(par.get('f_stretch_r_b%d_amp' % b, 0.0))
            if a != 0.0:
                xb = float(par['f_stretch_r_b%d_x' % b])
                w = float(par['f_stretch_r_b%d_w' % b])
                u += a * w * (np.tanh((xi - xb) / w) - (1 - xi) * np.tanh(-xb / w)
                              - xi * np.tanh((1 - xb) / w))
    return r0 + (r1 - r0) * u


def centroid(rl, rr):
    q = rl / rr
    return 0.25 * (q * q + 1.0) / ((1.0 / 3.0) * (q * q + q + 1.0)) * (rr + rl)


def load(run, tag, re_):
    fs = sorted(glob.glob(os.path.join(run, 'bin', 'hepresn.%s.*.bin' % tag)))
    out = []
    for f in fs:
        d = read_binary(f)
        v = {k: np.asarray(a)[0].mean(axis=(0, 1)) for k, a in d['mb_data'].items()}
        out.append((d['time'], d['cycle'], v))
    return out


def main():
    run = sys.argv[1]
    npz = np.load(sys.argv[2]) if len(sys.argv) > 2 else None
    par = mesh = read_mesh(run)
    hu = load(run, 'hydro_u', None)
    m1 = load(run, 'm1', None)
    nx = hu[0][2]['dens'].size
    re_ = edges(par, nx)
    rc = centroid(re_[:-1], re_[1:])
    dxc = np.diff(re_)
    t0, c0, u0 = hu[0]
    g = GM / rc**2
    rho0 = u0['dens']
    x = rc / RSTAR
    if npz is not None:
        tau = np.interp(rc, npz['r'], npz['tau'])
        fm = np.interp(rc, npz['r'], npz['fmlt'])
    else:
        tau = np.full(nx, 10.0)
        fm = np.zeros(nx)
    zones = [('rad 0.5-0.635 R', (x >= 0.5) & (x < 0.635)),
             ('FeCZ 0.635-0.968 R', (x >= 0.635) & (x < 0.968)),
             ('0.968 R - tau 1', (x >= 0.968) & (tau >= 1.0)),
             ('tau < 1', tau < 1.0)]
    zones[0][1][0] = False
    print('outputs: %d hydro_u, %d m1; nx1 = %d, r = %.5e .. %.5e; dx %.3e .. %.3e; '
          'final t = %.2f (cycle %d)' % (len(hu), len(m1), nx, rc[0], rc[-1], dxc.min(),
                                         dxc.max(), hu[-1][0], hu[-1][1]))
    steps = [(a, b) for a, b in zip(hu[:-1], hu[1:]) if b[0] > a[0]]
    print('\n(i) net radial force f/(rho g) = d(rho v_r)/dt/(rho g)')
    print('    %-22s %12s %12s %14s %10s' % ('zone', 'first step', 'max step', 'cum. at end',
                                          'max at r/R'))
    for name, z in zones:
        vals = [np.abs((b[2]['mom1'] - a[2]['mom1']) / (b[0] - a[0]) / (rho0 * g))[z]
                for a, b in steps]
        first = vals[0].max()
        mx = np.array([v.max() for v in vals])
        fend = (hu[-1][2]['mom1'] - u0['mom1']) / (hu[-1][0] - t0) / (rho0 * g)
        jj = np.argmax(np.abs(fend[z]))
        print('    %-22s %12.3e %12.3e %14.3e %10.4f'
              % (name, first, mx.max(), np.abs(fend[z]).max(), x[z][jj]))
    fw = [(b[2]['mom1'] - a[2]['mom1'])[0] / (b[0] - a[0]) / (rho0 * g)[0] for a, b in steps]
    print('    wall cell 0:           %12.3e %12.3e' % (abs(fw[0]), np.abs(fw).max()))
    print('    dt: first %.4g, min %.4g, max %.4g, median %.4g over %d steps'
          % (steps[0][1][0] - steps[0][0][0], min(b[0] - a[0] for a, b in steps),
             max(b[0] - a[0] for a, b in steps),
             np.median([b[0] - a[0] for a, b in steps]), len(steps)))
    print('\n(ii) shell luminosity 4 pi r^2 F_r (cell centred)')
    F0 = m1[0][2]['m1_f1']
    L0 = 4 * np.pi * rc**2 * F0
    for name, z in zones:
        print('    t=0 %-20s L/L_star min %.6f max %.6f' % (name, (L0 / LUM)[z].min(),
                                                            (L0 / LUM)[z].max()))
    for t, c, u in m1[1:][::max(1, len(m1) // 8)]:
        L = 4 * np.pi * rc**2 * u['m1_f1']
        print('    t %9.2f  max|L/L(0)-1|: ' % t + '  '.join(
            '%s %.2e' % (name.split()[0], np.abs(L / L0 - 1)[z].max()) for name, z in zones))
    L = 4 * np.pi * rc**2 * m1[-1][2]['m1_f1']
    print('    final: L/L_star at the end: ' + '  '.join(
        '%s [%.4f, %.4f]' % (name.split()[0], (L / LUM)[z].min(), (L / LUM)[z].max())
        for name, z in zones))
    hst = glob.glob(os.path.join(run, '*.user.hst'))
    if hst:
        with open(hst[0]) as fh:
            hdr = [ln for ln in fh if ln.startswith('#')]
        names = re.findall(r'\[\d+\]=(\S+)', hdr[-1])
        d = np.loadtxt(hst[0])
        col = {n: d[:, i] for i, n in enumerate(names)}
        t = col['time']
        de = col['Etot'] - col['Etot'][0]
        net = col['L_in'] - col['L_top']
        integ = np.concatenate([[0.0], np.cumsum(0.5 * np.diff(t) * (net[1:] + net[:-1]))])
        print('\n(iii) energy ledger: dEtot vs int (L_in - L_top) dt')
        i = -1
        print('    final t %.3f: dEtot %.4e int %.4e (dE - int)/(L_in t) = %.4e; L_top/L_in %.6f'
              % (t[i], de[i], integ[i], (de[i] - integ[i]) / (col['L_in'][i] * t[i]),
                 col['L_top'][i] / col['L_in'][i]))
        dt_ = np.diff(t)
        dE_ = np.diff(col['Etot'])
        trap = 0.5 * dt_ * (net[1:] + net[:-1])
        m = dt_ > 0
        r_ = np.abs((dE_ - trap) / (col['L_in'][1:] * dt_))[m]
        print('    per step |dE - trap|/(L_in dt): max %.3e (step %d), after step 3 %.3e'
              % (r_.max(), int(np.argmax(r_)), r_[3:].max()))
        print('    Mdot_top max %.4e Min_top min %.4e M_tot change %.3e'
              % (col['Mdot_top'].max(), col['Min_top'].min(),
                 col['M_tot'][-1] / col['M_tot'][0] - 1))
    if npz is not None:
        print('\n(iv) FeCZ gas heating: predicted -div F_r of the IC vs measured d e_int/dt')
        # predicted: G = -div F_r = -(L_r(r_r) - L_r(r_l))/V with L_r = L (1 - fmlt) (IC F_r)
        fml = np.interp(re_, npz['r'], npz['fmlt'])
        Lr = LUM * (1 - fml)
        vol = 4 * np.pi / 3 * (re_[1:]**3 - re_[:-1]**3)
        gpred = -(Lr[1:] - Lr[:-1]) / vol
        z = zones[1][1]
        for k in (1, len(hu) // 8, len(hu) // 2, len(hu) - 1):
            k = max(1, min(k, len(hu) - 1))
            t, c, u = hu[k]
            mom, dens, ener = u['mom1'], u['dens'], u['ener']
            e_int = ener - 0.5 * mom**2 / dens - dens * GM * (1.0 / re_[0] - 1.0 / rc)
            e0 = u0['ener'] - 0.5 * u0['mom1']**2 / u0['dens'] - u0['dens'] * GM * (
                1.0 / re_[0] - 1.0 / rc)
            gm = (e_int - e0) / (t - t0)
            sel = z & (np.abs(gpred) > 0.2 * np.abs(gpred[z]).max())
            print('    t = %8.3f: measured/predicted over |pred| > 0.2 max in the FeCZ: '
                  'median %.3f, min %.3f, max %.3f ; max pred %.3e, max meas %.3e '
                  'erg/cm3/s ; all-FeCZ rms(meas-pred)/rms(pred) %.3f'
                  % (t, np.median(gm[sel] / gpred[sel]), (gm[sel] / gpred[sel]).min(),
                     (gm[sel] / gpred[sel]).max(), np.abs(gpred[z]).max(),
                     np.abs(gm[z]).max(),
                     np.sqrt(np.mean((gm[z] - gpred[z])**2)) / np.sqrt(np.mean(gpred[z]**2))))
        # the radiative-source measure at the last output: -div F_r of the M1 field
        Fm = m1[min(1, len(m1) - 1)][2]['m1_f1']
        print('    pred FeCZ range of -div F_r: %.3e .. %.3e erg/cm3/s (heating where >0), '
              'rho e/t_heat at the rate max: e/G = %.1f s'
              % (gpred[z].min(), gpred[z].max(),
                 (u0['ener'] / np.maximum(np.abs(gpred), 1e-30))[z].min()))


if __name__ == '__main__':
    main()
