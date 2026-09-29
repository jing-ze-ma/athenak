#!/usr/bin/env python3
"""3-D wedge diagnostics for he_star_m1 (bin dumps of hydro_w and m1).
usage: analyze_3d.py RUNDIR INPUT STRUCT_NPZ [--all] [--every N]
Per dump (default: the last; --all = every dump, --every N = every N-th):
  * per shell (volume-weighted horizontal means, primes = deviation from the shell mean):
    v_r' rms / v_MLT(IC), corr(v_r', T'), convective luminosity of the fluctuations
    L_cg = 4 pi r^2 <v_r' (e+p)_gas'> and L_cr = 4 pi r^2 <v_r' (4/3) E'>, radiative
    L_rad = 4 pi r^2 <F_r>, kinetic L_k = 4 pi r^2 <v_r rho v^2/2>; T = (E/a)^(1/4)
    (T_gas = T_rad to the stiff exchange below the photosphere), p_gas from the run's EOS
    table at (rho, T);
  * the hydro dt proxy min dx_d/(|v_d| + c_gas) (c_gas^2 = 5/3 p_gas/rho) and where;
Writes RUNDIR/diag3d/<dump>.npz with the shell profiles."""
import glob
import os
import sys

import numpy as np
from scipy.interpolate import RectBivariateSpline

sys.path.insert(0, '/viper/ptmp2/jinma/wt_hepresn/vis/python')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bin_convert import read_binary  # noqa: E402
import analyze_gate_mlt as A  # noqa: E402

AR = 7.5657332503e-15
RSTAR = 2.3717e11
DUMP = '/viper/ptmp2/jinma/hepresn_0929/eos_table_he_x0_y0.98_z0.02_logd-14.txt'


def read_mesh(fn):
    par, sec = {}, None
    for ln in open(fn):
        ln = ln.split('#')[0].strip()
        if ln.startswith('<'):
            sec = ln.strip('<>')
        elif '=' in ln and sec == 'mesh':
            k, v = [x.strip() for x in ln.split('=', 1)]
            par[k] = v
    return par


def eos_pg():
    with open(DUMP) as fh:
        lines = [fh.readline() for _ in range(6)]
    nx, ny, xmin, dx, ymin, dy = [float(v) for v in lines[2].split()[1:]]
    nx, ny = int(nx), int(ny)
    d = np.loadtxt(DUMP, comments='#')
    x = xmin + dx * np.arange(nx)
    y = ymin + dy * np.arange(ny)
    # columns: log10 e (per volume?) and log10 p/rho: as make_ic_he_presn_m1.read_dump
    slp = RectBivariateSpline(y, x, d[:, 1].reshape(ny, nx), kx=1, ky=1)
    return slp


def blocks(d, name):
    return np.asarray(d['mb_data'][name])


def shells(fn_h, fn_m, re_, npz, slp):
    dh, dm = read_binary(fn_h), read_binary(fn_m)
    rho, v1, v2, v3 = [blocks(dh, k) for k in ('dens', 'velx', 'vely', 'velz')]
    eg = blocks(dh, 'eint')
    E = blocks(dm, 'm1_e')
    F1 = blocks(dm, 'm1_f1')
    nmb, n3, n2, n1 = rho.shape
    geo = np.asarray(dh['mb_geometry'])
    w = np.empty((nmb, n3, n2, 1))
    th = np.empty((nmb, n3, n2, 1))
    for m in range(nmb):
        t = np.linspace(geo[m, 2], geo[m, 3], n2 + 1)
        dp = (geo[m, 5] - geo[m, 4]) / n3
        w[m] = ((np.cos(t[:-1]) - np.cos(t[1:])) * dp)[None, :, None]
        th[m] = (0.5 * (t[:-1] + t[1:]))[None, :, None]
    ws = w.sum()

    def mean(q):
        return (q * w).sum(axis=(0, 1, 2)) / ws

    rc = A.centroid(re_[:-1], re_[1:])
    T = (np.maximum(E, 1e-30) / AR)**0.25
    pg = rho * 10.0**slp.ev(np.log10(T).ravel(), np.log10(rho).ravel()).reshape(rho.shape)
    h = eg + pg
    vr = v1 - mean(v1)
    Tp = T - mean(T)
    hp = h - mean(h)
    Ep = E - mean(E)
    area = 4 * np.pi * rc**2
    out = dict(t=dh['time'], cycle=dh['cycle'], r=rc,
               vrms=np.sqrt(mean(vr**2)), vhrms=np.sqrt(mean(v2**2 + v3**2)),
               corr=mean(vr * Tp) / np.sqrt(mean(vr**2) * mean(Tp**2) + 1e-300),
               Lcg=area * mean(vr * hp), Lcr=area * mean(vr * 4.0 / 3.0 * Ep),
               Lrad=area * mean(F1), Lk=area * mean(v1 * 0.5 * rho * (v1**2 + v2**2 + v3**2)),
               vmlt=np.interp(rc, npz['r'], npz['vmlt']), vr_mean=mean(v1))
    # dt proxy
    cs = np.sqrt(5.0 / 3.0 * pg / rho)
    dx1 = np.diff(re_)[None, None, None, :]
    dth = (geo[0, 3] - geo[0, 2]) / n2
    dph = (geo[0, 5] - geo[0, 4]) / n3
    rr = rc[None, None, None, :]
    dts = [dx1 / (np.abs(v1) + cs), rr * dth / (np.abs(v2) + cs),
           rr * np.sin(th) * dph / (np.abs(v3) + cs)]
    best = min(range(3), key=lambda q: dts[q].min())
    idx = np.unravel_index(dts[best].argmin(), dts[best].shape)
    out['dtmin'] = dts[best].min()
    out['dtloc'] = (best + 1, rc[idx[3]] / RSTAR, float(th[idx[0], idx[1], idx[2], 0]),
                    float(np.sqrt(v1[idx]**2 + v2[idx]**2 + v3[idx]**2)), float(cs[idx]))
    out['vmax'] = float(np.sqrt(v1**2 + v2**2 + v3**2).max())
    out['vmax_r'] = rc[np.unravel_index(np.sqrt(v1**2 + v2**2 + v3**2).argmax(), v1.shape)[3]] / RSTAR
    return out


def main():
    run, inp, npzf = sys.argv[1:4]
    every = 1
    if '--every' in sys.argv:
        every = int(sys.argv[sys.argv.index('--every') + 1])
    npz = np.load(npzf)
    slp = eos_pg()
    par = read_mesh(inp)
    nx = int(par['nx1'])
    re_ = A.edges(par, nx)
    fh = sorted(glob.glob(os.path.join(run, 'bin', '*.hydro_w.*.bin')))
    fm = sorted(glob.glob(os.path.join(run, 'bin', '*.m1.*.bin')))
    pairs = list(zip(fh, fm))
    if '--all' not in sys.argv and '--every' not in sys.argv:
        pairs = pairs[-1:]
    os.makedirs(os.path.join(run, 'diag3d'), exist_ok=True)
    L = A.LUM
    x = None
    for n, (a, b) in enumerate(pairs[::every]):
        o = shells(a, b, re_, npz, slp)
        np.savez(os.path.join(run, 'diag3d', os.path.basename(a).replace('.bin', '.npz')), **o)
        x = o['r'] / RSTAR
        fz = (x > 0.66) & (x < 0.91)
        print('t %9.1f cyc %6d | FeCZ: vrms/vMLT %.2e  corr %+.2f  (Lcg+Lcr)/L %+.3f  Lrad/L %.3f'
              '  Lk/L %+.1e | vmax %.2e at r/R %.3f | dt proxy %.2f s dir %d r/R %.3f'
              % (o['t'], o['cycle'], np.median(o['vrms'][fz] / o['vmlt'][fz]),
                 np.median(o['corr'][fz]), np.median((o['Lcg'] + o['Lcr'])[fz]) / L,
                 np.median(o['Lrad'][fz]) / L, np.median(o['Lk'][fz]) / L, o['vmax'],
                 o['vmax_r'], o['dtmin'], o['dtloc'][0], o['dtloc'][1]))
    if '--profile' in sys.argv and x is not None:
        for q in range(0, len(x), max(1, len(x) // 24)):
            print('  r/R %.3f vrms %.2e vMLT %.2e corr %+.2f Lcg/L %+.3f Lcr/L %+.3f Lrad/L %.3f'
                  % (x[q], o['vrms'][q], o['vmlt'][q], o['corr'][q], o['Lcg'][q] / L,
                     o['Lcr'][q] / L, o['Lrad'][q] / L))


if __name__ == '__main__':
    main()
