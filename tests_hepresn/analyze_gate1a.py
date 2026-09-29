#!/usr/bin/env python3
"""Gate 1a analysis of a he_star_m1 thin-column run (docs: tests_hepresn/README.md).
usage: analyze_gate1a.py RUNDIR R_INT [R_TAU_LUM]
  RUNDIR   run directory with hepresn.hst and bin/hepresn.hydro_u.*.bin, hepresn.m1.*.bin
  R_INT    radius of tau = 1 (the pgen prints 'r_int (tau = 1) = ...')
Prints: (i) net radial force f/(rho g) = d(rho v_r)/dt/(rho g) per cell, from the momentum
change since t = 0 divided by t, for every output; (ii) shell luminosity 4 pi r^2 F_r
against its t = 0 value and against L; (iii) the energy ledger from the history file.
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


def load(run, tag):
    fs = sorted(glob.glob(os.path.join(run, 'bin', 'hepresn.%s.*.bin' % tag)))
    out = []
    for f in fs:
        d = read_binary(f)
        v = {k: np.asarray(a).mean(axis=(0, 1, 2)) if False else
             np.asarray(a)[0].mean(axis=(0, 1)) for k, a in d['mb_data'].items()}
        r = 0.5*(d['x1min'] + d['x1max'])
        n = d['Nx1']
        dx = (d['x1max'] - d['x1min'])/n
        rc = d['x1min'] + (np.arange(n) + 0.5)*dx
        out.append((d['time'], d['cycle'], rc, v))
    return out


def main():
    run = sys.argv[1]
    rint = float(sys.argv[2])
    hu = load(run, 'hydro_u')
    m1 = load(run, 'm1')
    t0, c0, rc, u0 = hu[0]
    g = GM/rc**2
    rho0 = u0['dens']
    below = rc < rint
    print('outputs: %d hydro_u, %d m1; nx1 = %d, r = %.5e .. %.5e, cells below r_int: %d'
          % (len(hu), len(m1), len(rc), rc[0], rc[-1], below.sum()))
    print('\n(i) net radial force per cell f/(rho g) = d(rho v_r)/dt/(rho g), '
          'from the change since t=0')
    print('    %10s %6s %14s %14s %14s' % ('t[s]', 'cycle', 'max|.| tau>=1',
                                          'max|.| 4..N-4', 'max|.| all'))
    for t, c, r_, u in hu[1:]:
        f = (u['mom1'] - u0['mom1'])/(t - t0)/(rho0*g)
        print('    %10.3f %6d %14.4e %14.4e %14.4e'
              % (t, c, np.abs(f[below]).max(), np.abs(f[4:-4]).max(), np.abs(f).max()))
    # per-step force from consecutive outputs
    zones = [(0.0, 0.8), (0.8, 0.9), (0.9, 0.97), (0.97, 9.9)]
    print('    by zone (r/R ranges, tau>=1 only), max|f|/(rho g) since t=0:')
    for lo, hi in zones:
        z = below & (rc/RSTAR >= lo) & (rc/RSTAR < hi)
        row = []
        for t, c, r_, u in hu[1:]:
            f = (u['mom1'] - u0['mom1'])/(t - t0)/(rho0*g)
            row.append(np.abs(f[z]).max())
        print('      %.2f-%.2f R: ' % (lo, min(hi, 1.0)) +
              ' '.join('%.1e' % x for x in row))
    print('    consecutive:')
    for a, b in zip(hu[:-1], hu[1:]):
        f = (b[3]['mom1'] - a[3]['mom1'])/(b[0] - a[0])/(rho0*g)
        print('    t %10.3f -> %10.3f  dt %9.3f  max|.| tau>=1 %12.4e at r/R %.4f'
              % (a[0], b[0], b[0] - a[0], np.abs(f[below]).max(),
                 rc[below][np.argmax(np.abs(f[below]))]/RSTAR))
    print('\n(ii) shell luminosity 4 pi r^2 F_r (cell centred, lab)')
    F0 = m1[0][3]['m1_f1']
    L0 = 4*np.pi*rc**2*F0
    print('    t=0: L/L_star: min %.6f max %.6f (tau>=1)'
          % ((L0/LUM)[below].min(), (L0/LUM)[below].max()))
    for t, c, r_, u in m1[1:]:
        L = 4*np.pi*rc**2*u['m1_f1']
        print('    t %9.2f  max|L/L(0)-1| tau>=1: %.3e   max|L/L_star-1| tau>=1: %.3e '
              '(all radiative r<0.635R: %.3e)'
              % (t, np.abs(L/L0 - 1)[below].max(), np.abs(L/LUM - 1)[below].max(),
                 np.abs(L/LUM - 1)[rc < 0.635*RSTAR].max()))
    hst = glob.glob(os.path.join(run, '*.user.hst'))
    if hst:
        with open(hst[0]) as fh:
            hdr = [ln for ln in fh if ln.startswith('#')]
        names = re.findall(r'\[\d+\]=(\S+)', hdr[-1])
        d = np.loadtxt(hst[0])
        col = {n: d[:, i] for i, n in enumerate(names)}
        t = col['time']
        de = col['Etot'] - col['Etot'][0]
        integ = np.concatenate([[0.0], np.cumsum(0.5*np.diff(t)*(
            (col['L_in'] - col['L_top'])[1:] + (col['L_in'] - col['L_top'])[:-1]))])
        print('\n(iii) energy ledger: dEtot vs int (L_in - L_top) dt')
        print('    %10s %14s %14s %14s %14s' % ('t', 'dEtot', 'int(Lin-Ltop)dt',
                                                'diff/(L_in t)', 'L_top/L_in'))
        for i in range(0, len(t), max(1, len(t)//12)):
            if t[i] > 0:
                print('    %10.3f %14.5e %14.5e %14.4e %14.6f'
                      % (t[i], de[i], integ[i], (de[i] - integ[i])/(col['L_in'][i]*t[i]),
                         col['L_top'][i]/col['L_in'][i]))
        i = -1
        print('    final t %.3f: (dE - int)/(L_in t) = %.4e' %
              (t[i], (de[i] - integ[i])/(col['L_in'][i]*t[i])))
        dt_ = np.diff(t)
        dE_ = np.diff(col['Etot'])
        right = dt_*(col['L_in'][1:] - col['L_top'][1:])
        trap = 0.5*dt_*((col['L_in'] - col['L_top'])[1:] + (col['L_in'] - col['L_top'])[:-1])
        m = dt_ > 0
        print('    per step |dE - dt (L_in - L_top)|/(L_in dt): right-rect max %.3e, '
              'trapezoid max %.3e (steps %d); after step 3: %.3e / %.3e'
              % (np.abs((dE_ - right)/(col['L_in'][1:]*dt_))[m].max(),
                 np.abs((dE_ - trap)/(col['L_in'][1:]*dt_))[m].max(), m.sum(),
                 np.abs((dE_ - right)/(col['L_in'][1:]*dt_))[m][3:].max(),
                 np.abs((dE_ - trap)/(col['L_in'][1:]*dt_))[m][3:].max()))
        print('    Mdot_top max %.4e Min_top min %.4e M_tot change %.3e' %
              (col['Mdot_top'].max(), col['Min_top'].min(),
               col['M_tot'][-1]/col['M_tot'][0] - 1))


if __name__ == '__main__':
    main()
