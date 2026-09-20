#!/usr/bin/env python3
"""tests_r19: POROSITY test for the wedge9d luminosity excess.  READ-ONLY on the run.

Does the horizontally inhomogeneous density at nearly uniform T make the shell-averaged
radiative flux exceed the flux of the shell-mean state?

usage: poros.py [dumpindex ...]
"""
import os
import sys

import numpy as np

sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert as bc   # noqa: E402

RUN = '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r11/wedge9d'
EOSTBL = ('/viper/u2/jinma/ATHENAK/bench/analysis_0916/iso3d_w5/hestar/'
          'eos_table_box_w5.txt')
OPAC = '/viper/u2/jinma/ATHENAK/bench/hestar_fecz/rosseland_he_x0.0_z0.02.txt'
RS = 2.3717e11
LSTAR = 2.3066e38
TURN = 4705.0
ARAD = 7.5657332503e-15
CLIGHT = 2.99792458e10
SIGSB = 5.670374419e-5
RAD_LO, RAD_HI = np.log10(1.0e-11), np.log10(1.0e-10)   # eos_rad_rho_lo/hi
CPOLY = [0.992525, -0.404821, -0.372255, -1.499759]

# ------------------------------------------------------------------ EOS table (dt3d.py)


def load_tbl(fn):
    with open(fn) as fh:
        ln = [fh.readline() for _ in range(6)]
    nx, ny, xmin, dx, ymin, dy = [float(v) for v in ln[2].split()[1:]]
    nx, ny = int(nx), int(ny)
    d = np.loadtxt(fn, comments='#')
    E = d[:, 0].reshape(ny, nx)
    return dict(nx=nx, ny=ny, xmin=xmin, dx=dx, ymin=ymin, dy=dy, E=E)


TB = load_tbl(EOSTBL)


def bil(F, ix, iy, u, v):
    return ((1-u)*(1-v)*F[iy, ix] + u*(1-v)*F[iy, ix+1]
            + (1-u)*v*F[iy+1, ix] + u*v*F[iy+1, ix+1])


def cellidx(x, y):
    gx = (x - TB['xmin'])/TB['dx']
    gy = (y - TB['ymin'])/TB['dy']
    ix = np.clip(np.floor(gx).astype(np.int32), 0, TB['nx']-2)
    iy = np.clip(np.floor(gy).astype(np.int32), 0, TB['ny']-2)
    return ix, iy, gx-ix, gy-iy


def weight(x):
    s = (x - RAD_LO)/(RAD_HI - RAD_LO)
    return np.where(s <= 0, 0.0, np.where(s >= 1, 1.0, s*s*(3-2*s)))


def temp_of(rho, eint):
    """solve e_gas(rho,T) + w(rho) a T^4 = eint, bisection in log10 T (dt3d.py)."""
    x = np.log10(np.maximum(rho, 1e-300))
    w = weight(x)
    lo = np.full(rho.shape, 2.5)
    hi = np.full(rho.shape, 7.0)
    for _ in range(50):
        mid = 0.5*(lo+hi)
        ix, iy, u, v = cellidx(x, mid)
        f = rho*10.0**bil(TB['E'], ix, iy, u, v) + w*ARAD*(10.0**mid)**4 - eint
        lo = np.where(f < 0.0, mid, lo)
        hi = np.where(f < 0.0, hi, mid)
    return 10.0**(0.5*(lo+hi))


# ------------------------------------------------------------------ opacity (relax_ic.py)
def load_opac(fn):
    hdr, vals = [], []
    for ln in open(fn):
        (hdr if ln.startswith('#') else vals).append(ln)
    g = [h for h in hdr if h.strip().startswith('# 2')]
    nT, nD, lTmin, dlT, lDmin, dlD = [float(v) for v in g[-1].strip('# \n').split()]
    k = np.array([float(v) for v in vals]).reshape(int(nT), int(nD))
    return k, lTmin, dlT, int(nT), lDmin, dlD, int(nD)


OP = load_opac(OPAC)


def kappa(T, rho):
    k, lTmin, dlT, nT, lDmin, dlD, nD = OP
    x = np.clip((np.log10(T) - lTmin)/dlT, 0, nT - 1.000001)
    y = np.clip((np.log10(rho) - lDmin)/dlD, 0, nD - 1.000001)
    i0, j0 = np.floor(x).astype(int), np.floor(y).astype(int)
    fx, fy = x - i0, y - j0
    return 10.0**((1-fx)*(1-fy)*k[i0, j0] + fx*(1-fy)*k[i0+1, j0]
                  + (1-fx)*fy*k[i0, j0+1] + fx*fy*k[i0+1, j0+1])


# ------------------------------------------------------------------ grid (dt3d.py)
def stretch(xi):
    u = xi.copy()
    xik = xi.copy()
    for kk in range(1, 5):
        u += CPOLY[kk-1]*xik*(1.0-xi)
        xik = xik*xi
    return u


def grid(h):
    n1, n2, n3 = h['Nx1'], h['Nx2'], h['Nx3']
    rf = h['x1min'] + (h['x1max']-h['x1min'])*stretch(np.linspace(0.0, 1.0, n1+1))
    rl, rr = rf[:-1], rf[1:]
    q = rl/rr
    rc = 0.25*(q*q+1.0)/((1.0/3.0)*(q*q+q+1.0))*(rr+rl)
    dr = rr-rl
    thf = np.linspace(h['x2min'], h['x2max'], n2+1)
    tl, tr = thf[:-1], thf[1:]
    thc = -((np.sin(tr)-tr*np.cos(tr)) - (np.sin(tl)-tl*np.cos(tl)))/(np.cos(tr)-np.cos(tl))
    wj = np.cos(tl) - np.cos(tr)          # solid-angle weight per theta row
    return rc, dr, thc, wj


def load(idx):
    d = bc.read_binary(os.path.join(RUN, 'bin', 'he4.hydro_w.%05d.bin' % idx))
    n1, n2, n3 = d['Nx1'], d['Nx2'], d['Nx3']
    out = {}
    for nm in ('dens', 'eint', 'velx'):
        g = np.empty((n3, n2, n1), dtype=np.float64)
        src = d['mb_data'][nm]
        nj, nk = d['nx2_mb'], d['nx3_mb']
        for m in range(d['n_mbs']):
            _, lx2, lx3, _ = d['mb_logical'][m]
            g[lx3*nk:(lx3+1)*nk, lx2*nj:(lx2+1)*nj, :] = src[m]
        out[nm] = g
    return d, out


# ------------------------------------------------------------------ analysis
def shmean(F, W):
    """area-weighted horizontal mean -> (n1,).  F (n3,n2,n1), W (n2,)."""
    return np.tensordot(W, F.sum(axis=0), axes=([0], [0]))/(W.sum()*F.shape[0])


def nyq_share(X):
    """share of variance at 2-cell (Nyquist) wavelength in theta and phi, and at
    wavelengths <= 4 cells, for a 2-D periodic shell map X (n3,n2)."""
    n3, n2 = X.shape
    A = np.fft.fft2(X - X.mean())
    P = (np.abs(A)**2)
    P[0, 0] = 0.0
    tot = P.sum()
    if tot <= 0:
        return 0.0, 0.0, 0.0
    k3 = np.fft.fftfreq(n3)*n3
    k2 = np.fft.fftfreq(n2)*n2
    K3, K2 = np.meshgrid(k3, k2, indexing='ij')
    nyq2 = P[:, n2//2].sum()/tot            # 2-cell in theta
    nyq3 = P[n3//2, :].sum()/tot            # 2-cell in phi
    small = P[(np.abs(K2) >= n2/4) | (np.abs(K3) >= n3/4)].sum()/tot   # <= 4 cells
    return nyq2, nyq3, small


def analyse(idx, rows, spec_rows):
    d, g = load(idx)
    rc, dr, thc, wj = grid(d)
    n1, n2, n3 = d['Nx1'], d['Nx2'], d['Nx3']
    rho, eint = g['dens'], g['eint']
    T = np.empty_like(rho)
    for a in range(0, n1, 16):
        b = min(a+16, n1)
        T[:, :, a:b] = temp_of(rho[:, :, a:b], eint[:, :, a:b])
    kap = kappa(T, rho)
    kr = kap*rho
    # per-column diffusion flux
    dTdr = np.gradient(T, rc, axis=2)
    Fcol = -(4.0*ARAD*CLIGHT/3.0)*T**3/kr*dTdr
    # shell means
    rm = shmean(rho, wj)
    Tm = shmean(T, wj)
    krm = shmean(kr, wj)
    ikrm = shmean(1.0/kr, wj)
    Fcm = shmean(Fcol, wj)
    T3m = shmean(T**3, wj)
    kapm = kappa(Tm, rm)
    dTmdr = np.gradient(Tm, rc)
    Fmean = -(4.0*ARAD*CLIGHT/3.0)*Tm**3/(kapm*rm)*dTmdr
    Fmix = -(4.0*ARAD*CLIGHT/3.0)*T3m*dTmdr*ikrm       # <1/(kr)> with shell-mean dT/dr
    Freq = LSTAR/(4.0*np.pi*rc**2)
    P1 = ikrm*krm
    P2 = ikrm*(kapm*rm)
    rrms = np.sqrt(shmean((rho/rm - 1.0)**2, wj))
    trms = np.sqrt(shmean((T/Tm - 1.0)**2, wj))
    # weighted correlation of Fcol with rho on each shell
    W = np.broadcast_to(wj[None, :, None], rho.shape)
    a1 = Fcol - Fcm
    a2 = rho - rm
    cov = shmean(a1*a2, wj)
    cc = cov/np.sqrt(shmean(a1*a1, wj)*shmean(a2*a2, wj))
    # limiter and optical depth per cell
    dtau = kr*dr[None, None, :]
    dtaum = shmean(dtau, wj)
    flim = np.abs(Fcol)/(4.0*SIGSB*T**4)
    flimm = shmean(flim, wj)
    flimmax = flim.max(axis=(0, 1))
    spec = {}
    for i in spec_rows:
        spec[i] = (nyq_share(rho[:, :, i]/rm[i]), nyq_share(T[:, :, i]/Tm[i]))
    return dict(t=d['time'], rc=rc, rm=rm, Tm=Tm, P1=P1, P2=P2, rrms=rrms, trms=trms,
                Fcm=Fcm, Fmean=Fmean, Fmix=Fmix, Freq=Freq, cc=cc, dtaum=dtaum,
                flimm=flimm, flimmax=flimmax, spec=spec, Tmax=T.max(axis=(0, 1)))


if __name__ == '__main__':
    ids = [int(a) for a in sys.argv[1:]] or [5, 14, 16, 20]
    r0 = None
    outs = []
    for i in ids:
        if r0 is None:
            dd = bc.read_binary(os.path.join(RUN, 'bin', 'he4.hydro_w.%05d.bin' % i))
            rc0, _, _, _ = grid(dd)
            del dd
            r0 = rc0/RS
            rows = [k for k in range(0, len(r0), 4) if 0.55 <= r0[k] <= 1.001]
            spec_rows = [k for k in rows if 0.58 <= r0[k] <= 0.99][::3]
        A = analyse(i, rows, spec_rows)
        outs.append((i, A))
        print('\n' + '='*100)
        print('dump %d   t = %.1f s = %.4f turnover' % (i, A['t'], A['t']/TURN))
        print('%4s %7s %10s %10s %8s %8s %8s %8s %9s %9s %9s %7s %10s %9s'
              % ('i', 'r/R', '<rho>', '<T>', 'P1', 'P2', 'rms_rho', 'rms_T',
                 '<Fcol>/Fr', 'Fmean/Fr', 'Fmix/Fr', 'corr', 'dtau_cell', 'Flim_f'))
        for k in rows:
            print('%4d %7.4f %10.3e %10.4e %8.3f %8.3f %8.4f %8.5f %9.4f %9.4f %9.4f '
                  '%7.3f %10.3e %9.2e'
                  % (k, r0[k], A['rm'][k], A['Tm'][k], A['P1'][k], A['P2'][k],
                     A['rrms'][k], A['trms'][k], A['Fcm'][k]/A['Freq'][k],
                     A['Fmean'][k]/A['Freq'][k], A['Fmix'][k]/A['Freq'][k],
                     A['cc'][k], A['dtaum'][k], A['flimm'][k]))
        print('  max cell limiter argument f/(4 sigma T^4) over rows 0.55-1.0 R: %.3e'
              % max(A['flimmax'][k] for k in rows))
        print('\n  --- 2-cell power share (Nyquist theta / Nyquist phi / <=4 cells) ---')
        print('%4s %7s   %-28s %-28s' % ('i', 'r/R', 'rho/<rho>', 'T/<T>'))
        for k in spec_rows:
            (a2, a3, asm), (b2, b3, bsm) = A['spec'][k]
            print('%4d %7.4f   %7.4f %7.4f %7.4f      %7.4f %7.4f %7.4f'
                  % (k, r0[k], a2, a3, asm, b2, b3, bsm))

    print('\n' + '='*100)
    print('RATIO <Fcol>/Fmean (the porosity enhancement of the shell-averaged flux)')
    hdr = '%4s %7s' % ('i', 'r/R')
    for i, A in outs:
        hdr += ' %9.3ft' % (A['t']/TURN)
    print(hdr)
    for k in rows:
        ln = '%4d %7.4f' % (k, r0[k])
        for i, A in outs:
            ln += ' %10.4f' % (A['Fcm'][k]/A['Fmean'][k])
        print(ln)
    print('\nRATIO Fmix/Fmean  ( = P2, shell-mean dT/dr with per-column 1/(kappa rho) )')
    print(hdr)
    for k in rows:
        ln = '%4d %7.4f' % (k, r0[k])
        for i, A in outs:
            ln += ' %10.4f' % (A['Fmix'][k]/A['Fmean'][k])
        print(ln)
