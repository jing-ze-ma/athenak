#!/usr/bin/env python3
"""Porosity test on he_star_m1 3-D dumps (hydro_w + m1 bin pairs).
usage: porosity.py RUNDIR INPUT DUMPNUM [DUMPNUM ...]
Per shell (every 8th radial cell printed, FeCZ summary line):
  drho = rms(rho')/<rho>;  L1 = 4 pi r^2 <F_1> / L (M1 radial flux, the run's);
  Lloc = 4 pi r^2 <-c/(3 kappa rho) dE/dr> / L (local diffusion flux per cell);
  Lmean = 4 pi r^2 (-c/(3 kappa(<rho>,<T>) <rho>) d<E>/dr) / L (diffusion flux of the
  shell-mean state); P = L1/Lmean (porosity factor); corr(F_1', rho');
  hk = fraction of the horizontal rho' variance at wavelengths <= 4 cells (|k| >= n/4 of
  the 64 x 64 shell), i.e. unresolved lumpiness."""
import os
import sys

import numpy as np

sys.path.insert(0, '/viper/ptmp2/jinma/wt_hepresn/vis/python')
sys.path.insert(0, '/viper/ptmp2/jinma/wt_hepresn/tests_hepresn')
from bin_convert import read_binary  # noqa: E402
import analyze_gate_mlt as A  # noqa: E402
import analyze_3d as B  # noqa: E402

AR = 7.5657332503e-15
C = 2.99792458e10
L = 2.3066e38
KT = ('/viper/ptmp2/jinma/caltech_handover_0926/athenak_data/he_box/'
      'rosseland_he_x0.0_z0.02.txt')


def kap_table():
    hdr = [ln for ln in open(KT) if ln.startswith('#')]
    g = [float(v) for v in hdr[-2].split()[1:]]
    nT, nD, lT0, dlT, lD0, dlD = int(g[0]), int(g[1]), g[2], g[3], g[4], g[5]
    k = np.loadtxt(KT, comments='#').reshape(nT, nD)

    def kap(T, rho):
        x = np.clip((np.log10(T) - lT0) / dlT, 0, nT - 1.000001)
        y = np.clip((np.log10(rho) - lD0) / dlD, 0, nD - 1.000001)
        i, j = x.astype(int), y.astype(int)
        fx, fy = x - i, y - j
        v = ((1 - fx) * (1 - fy) * k[i, j] + fx * (1 - fy) * k[i + 1, j]
             + (1 - fx) * fy * k[i, j + 1] + fx * fy * k[i + 1, j + 1])
        return 10.0**v
    return kap


def main():
    run, inp = sys.argv[1], sys.argv[2]
    par = B.read_mesh(inp)
    re_ = A.edges(par, int(par['nx1']))
    rc = A.centroid(re_[:-1], re_[1:])
    kap = kap_table()
    for num in sys.argv[3:]:
        dh = read_binary(os.path.join(run, 'bin', 'hepresn.hydro_w.%05d.bin' % int(num)))
        dm = read_binary(os.path.join(run, 'bin', 'hepresn.m1.%05d.bin' % int(num)))
        rho = np.asarray(dh['mb_data']['dens'])
        E = np.asarray(dm['mb_data']['m1_e'])
        F1 = np.asarray(dm['mb_data']['m1_f1'])
        geo = np.asarray(dh['mb_geometry'])
        nmb, n3, n2, n1 = rho.shape
        w = np.empty((nmb, n3, n2, 1))
        for m in range(nmb):
            t = np.linspace(geo[m, 2], geo[m, 3], n2 + 1)
            w[m] = ((np.cos(t[:-1]) - np.cos(t[1:])))[None, :, None]
        ws = w.sum()

        def mean(q):
            return (q * w).sum(axis=(0, 1, 2)) / ws
        T = (np.maximum(E, 1e-30) / AR)**0.25
        kr = kap(T, rho) * rho
        dEdr = np.gradient(E, rc, axis=3)
        Floc = -C / (3 * kr) * dEdr
        rm, Em = mean(rho), mean(E)
        Tm = (Em / AR)**0.25
        Fmean = -C / (3 * kap(Tm, rm) * rm) * np.gradient(Em, rc)
        a = 4 * np.pi * rc**2 / L
        L1, Lloc, Lmn = a * mean(F1), a * mean(Floc), a * Fmean
        drho = np.sqrt(mean((rho - rm)**2)) / rm
        fp, rp = F1 - mean(F1), rho - rm
        corr = mean(fp * rp) / np.sqrt(mean(fp**2) * mean(rp**2) + 1e-300)
        # assemble the 64 x 64 shells: blocks ordered by (phi_min, theta_min)
        t0 = sorted(set(np.round(geo[:, 2], 8)))
        p0 = sorted(set(np.round(geo[:, 4], 8)))
        sh = np.empty((len(p0) * n3, len(t0) * n2, n1))
        for m in range(nmb):
            it = t0.index(np.round(geo[m, 2], 8))
            ip = p0.index(np.round(geo[m, 4], 8))
            sh[ip * n3:(ip + 1) * n3, it * n2:(it + 1) * n2, :] = rho[m] / rm
        N3, N2 = sh.shape[:2]
        P = np.abs(np.fft.fft2(sh - sh.mean(axis=(0, 1)), axes=(0, 1)))**2
        k3 = np.abs(np.fft.fftfreq(N3) * N3)[:, None]
        k2 = np.abs(np.fft.fftfreq(N2) * N2)[None, :]
        hi = (np.maximum(k3 / N3, k2 / N2) >= 0.25)
        hk = (P * hi[:, :, None]).sum(axis=(0, 1)) / (P.sum(axis=(0, 1)) + 1e-300)
        x = rc / A.RSTAR
        fz = (x > 0.66) & (x < 0.91)
        print('== dump %s t %.0f s: FeCZ median drho %.3f  L1 %.3f  Lloc %.3f  '
              'Lmean %.3f  '
              'P=L1/Lmean %.3f  corr(F1,rho) %+.2f  hk %.2f'
              % (num, dh['time'], np.median(drho[fz]), np.median(L1[fz]),
                 np.median(Lloc[fz]),
                 np.median(Lmn[fz]), np.median((L1 / Lmn)[fz]), np.median(corr[fz]),
                 np.median(hk[fz])))
        for i in range(0, n1, 12):
            print('  r/R %.3f drho %.3e L1 %.3f Lloc %.3f Lmean %.3f P %.3f '
                  'corr %+.2f hk %.2f'
                  % (x[i], drho[i], L1[i], Lloc[i], Lmn[i], L1[i] / Lmn[i], corr[i],
                     hk[i]))


if __name__ == '__main__':
    main()
