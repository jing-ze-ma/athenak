#!/usr/bin/env python3
"""tests_r17x: WHERE the 3-D wedge7 hydro dt is set.  READ-ONLY.

Per-cell dt = CFL * min_d dx_d/(|v_d| + c_s) on the spherical-polar wedge, with the
run's own EOS (general table + density-tapered radiation), reproducing
src/eos/eos_table.hpp EvalFromLogs and src/hydro/hydro_newdt.cpp.

usage: dt3d.py <dumpindex> [<dumpindex> ...]
"""
import os
import sys

import numpy as np

sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert as bc   # noqa: E402

RUN = '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r11/wedge7'
EOSTBL = ('/viper/u2/jinma/ATHENAK/bench/analysis_0916/iso3d_w5/hestar/'
          'eos_table_box_w5.txt')
RS = 2.3717e11
TURN = 4705.0
ARAD = 7.5657332503e-15
CFL = 0.15
DFLOOR = 1.0e-13
RAD_LO, RAD_HI = np.log10(1.0e-11), np.log10(1.0e-10)
CPOLY = [0.992525, -0.404821, -0.372255, -1.499759]

# ---------------------------------------------------------------- EOS table


def load_tbl(fn):
    with open(fn) as fh:
        ln = [fh.readline() for _ in range(6)]
    nx, ny, xmin, dx, ymin, dy = [float(v) for v in ln[2].split()[1:]]
    nx, ny = int(nx), int(ny)
    d = np.loadtxt(fn, comments='#')
    E = d[:, 0].reshape(ny, nx)      # log10(e_gas/rho)
    P = d[:, 1].reshape(ny, nx)      # log10(p_gas/rho)
    return dict(nx=nx, ny=ny, xmin=xmin, dx=dx, ymin=ymin, dy=dy, E=E, P=P,
                Ey=np.gradient(E, dy, axis=0), Px=np.gradient(P, dx, axis=1),
                Py=np.gradient(P, dy, axis=0))


TB = load_tbl(EOSTBL)


def bil(F, gx, gy, ix, iy, u, v):
    return ((1-u)*(1-v)*F[iy, ix] + u*(1-v)*F[iy, ix+1]
            + (1-u)*v*F[iy+1, ix] + u*v*F[iy+1, ix+1])


def cellidx(x, y):
    """bicell indices + LOCAL coords, unclamped -> bilinear extrapolation off-table
    (the code linearly continues at the boundary slope; this is the same to 1st order)"""
    gx = (x - TB['xmin'])/TB['dx']
    gy = (y - TB['ymin'])/TB['dy']
    ix = np.clip(np.floor(gx).astype(np.int32), 0, TB['nx']-2)
    iy = np.clip(np.floor(gy).astype(np.int32), 0, TB['ny']-2)
    return ix, iy, gx-ix, gy-iy


def weight(x):
    """rad_taper::Weight(log10 rho) and dw/dx"""
    s = (x - RAD_LO)/(RAD_HI - RAD_LO)
    w = np.where(s <= 0, 0.0, np.where(s >= 1, 1.0, s*s*(3-2*s)))
    dw = np.where((s <= 0) | (s >= 1), 0.0, 6*s*(1-s)/(RAD_HI-RAD_LO))
    return w, dw


def egas_of(x, y, rho):
    ix, iy, u, v = cellidx(x, y)
    return rho*10.0**bil(TB['E'], x, y, ix, iy, u, v)


def temp_of(rho, eint):
    """solve e_gas(rho,T) + w(rho) a T^4 = eint for T, bisection in log10 T."""
    x = np.log10(np.maximum(rho, 1e-300))
    w, _ = weight(x)
    lo = np.full(rho.shape, 2.5)
    hi = np.full(rho.shape, 7.0)
    for _ in range(60):
        mid = 0.5*(lo+hi)
        T = 10.0**mid
        f = egas_of(x, mid, rho) + w*ARAD*T**4 - eint
        lo = np.where(f < 0.0, mid, lo)
        hi = np.where(f < 0.0, hi, mid)
    return 10.0**(0.5*(lo+hi))


def eos_state(rho, eint):
    """T, p_tot, Gamma1, c_s exactly as EvalFromLogs + SoundSpeedFromP build them."""
    T = temp_of(rho, eint)
    x, y = np.log10(np.maximum(rho, 1e-300)), np.log10(T)
    ix, iy, u, v = cellidx(x, y)
    ev = bil(TB['E'], x, y, ix, iy, u, v)
    evy = bil(TB['Ey'], x, y, ix, iy, u, v)
    pv = bil(TB['P'], x, y, ix, iy, u, v)
    pvx = bil(TB['Px'], x, y, ix, iy, u, v)
    pvy = bil(TB['Py'], x, y, ix, iy, u, v)
    egas = rho*10.0**ev
    pgas = rho*10.0**pv
    cv_g = (egas/rho)*evy/T
    w, dwdx = weight(x)
    erad0 = ARAD*T**4
    erad = w*erad0
    prad = erad/3.0
    p = pgas + prad
    chir = pgas*(1.0+pvx)/p + (erad0/3.0)*dwdx*np.log10(np.e)/p
    chit = (pgas*pvy + 4.0*prad)/p
    cv = cv_g + 4.0*erad/(rho*T)
    g1 = chir + p*chit*chit/(rho*T*cv)
    cs = np.sqrt(np.maximum(g1*p/rho, 0.0))
    return T, p, g1, cs

# ---------------------------------------------------------------- grid


def stretch(xi):
    u = xi.copy()
    xik = xi.copy()
    for k in range(1, 5):
        u += CPOLY[k-1]*xik*(1.0-xi)
        xik = xik*xi
    return u


def grid(h):
    n1, n2, n3 = h['Nx1'], h['Nx2'], h['Nx3']
    r0, r1 = h['x1min'], h['x1max']
    rf = r0 + (r1-r0)*stretch(np.linspace(0.0, 1.0, n1+1))
    rl, rr = rf[:-1], rf[1:]
    q = rl/rr
    rc = 0.25*(q*q+1.0)/((1.0/3.0)*(q*q+q+1.0))*(rr+rl)      # Mignone centroid
    dr = rr-rl
    thf = np.linspace(h['x2min'], h['x2max'], n2+1)
    tl, tr = thf[:-1], thf[1:]
    thc = -((np.sin(tr)-tr*np.cos(tr)) - (np.sin(tl)-tl*np.cos(tl)))/(np.cos(tr)-np.cos(tl))
    dth = tr-tl
    phf = np.linspace(h['x3min'], h['x3max'], n3+1)
    phc = 0.5*(phf[1:]+phf[:-1])
    dph = phf[1]-phf[0]
    return rc, dr, thc, dth, phc, dph

# ---------------------------------------------------------------- assembly


def load(idx):
    fn = os.path.join(RUN, 'bin', 'he4.hydro_w.%05d.bin' % idx)
    d = bc.read_binary(fn)
    n1, n2, n3 = d['Nx1'], d['Nx2'], d['Nx3']
    out = {}
    for nm in ('dens', 'velx', 'vely', 'velz', 'eint'):
        g = np.full((n3, n2, n1), np.nan, dtype=np.float64)
        src = d['mb_data'][nm]
        nj, nk = d['nx2_mb'], d['nx3_mb']
        for m in range(d['n_mbs']):
            _, lx2, lx3, _ = d['mb_logical'][m]
            g[lx3*nk:(lx3+1)*nk, lx2*nj:(lx2+1)*nj, :] = src[m]
        assert not np.isnan(g).any()
        out[nm] = g
    return d, out


def analyse(idx, nworst=10):
    d, g = load(idx)
    t = d['time']
    rc, dr, thc, dth, phc, dph = grid(d)
    n1, n2, n3 = d['Nx1'], d['Nx2'], d['Nx3']
    rho, eint = g['dens'], g['eint']
    # EOS in radial slabs to keep memory modest
    T = np.empty_like(rho)
    cs = np.empty_like(rho)
    for a in range(0, n1, 24):
        b = min(a+24, n1)
        T[:, :, a:b], _, _, cs[:, :, a:b] = eos_state(rho[:, :, a:b], eint[:, :, a:b])
    dx1 = dr[None, None, :]
    dx2 = (rc[None, :]*dth[:, None]).T[None, :, :]                 # (1,n2,n1)
    dx2 = (dth[:, None]*rc[None, :])[None, :, :]
    sn = np.abs(np.sin(thc))
    dx3 = (sn[:, None]*rc[None, :]*dph)[None, :, :]
    s1 = dx1/(np.abs(g['velx'])+cs)
    s2 = dx2/(np.abs(g['vely'])+cs)
    s3 = dx3/(np.abs(g['velz'])+cs)
    sig = CFL*np.minimum(np.minimum(s1, s2), s3)
    lim = np.argmin(np.stack([s1*np.ones_like(sig), s2*np.ones_like(sig),
                              s3*np.ones_like(sig)]), axis=0)
    return dict(t=t, d=d, rho=rho, T=T, cs=cs, sig=sig, lim=lim,
                v=(g['velx'], g['vely'], g['velz']),
                rc=rc, dr=dr, thc=thc, dth=dth, phc=phc, dph=dph,
                n=(n1, n2, n3), nworst=nworst)


def report(A):
    n1, n2, n3 = A['n']
    rho, T, cs, sig = A['rho'], A['T'], A['cs'], A['sig']
    vx, vy, vz = A['v']
    rc = A['rc']
    rmean = rho.mean(axis=(0, 1))
    tmean = T.mean(axis=(0, 1))
    print('\n############ t = %.1f s = %.4f turnover   cycle %d'
          % (A['t'], A['t']/TURN, A['d']['cycle']))
    print('  dt_min = %.4f s   (CFL %.2f)' % (sig.min(), CFL))
    flat = sig.ravel()
    order = np.argsort(flat)[:A['nworst']]
    print('\n  --- %d worst cells ---' % A['nworst'])
    print('   dt[s]   r/R    th[d]  ph[d]   i   j   k  jb kb   rho        T'
          '        |v|      v_r      cs       dir  <rho>_sh   <T>_sh   rho/<rho> T/<T>')
    for f in order:
        k, j, i = np.unravel_index(f, sig.shape)
        vv = np.sqrt(vx[k, j, i]**2 + vy[k, j, i]**2 + vz[k, j, i]**2)
        print('  %7.4f %6.4f %6.2f %6.2f %4d %3d %3d %3d %3d %9.3e %8.3e %8.2e '
              '%8.2e %8.2e  %s %9.3e %8.3e %8.3f %7.3f'
              % (sig[k, j, i], rc[i]/RS, np.degrees(A['thc'][j]),
                 np.degrees(A['phc'][k]), i, j, k, j % 24, k % 24,
                 rho[k, j, i], T[k, j, i], vv, vx[k, j, i], cs[k, j, i],
                 'r th ph'.split()[A['lim'][k, j, i]], rmean[i], tmean[i],
                 rho[k, j, i]/rmean[i], T[k, j, i]/tmean[i]))
    return rmean, tmean


def shell_table(As, every=6):
    n1 = As[0]['n'][0]
    rc = As[0]['rc']
    print('\n=== dt_min per radial shell (s), every %d shells ===' % every)
    hdr = '   i   r/R     dr[cm]  '
    for A in As:
        hdr += '  dt(%.1ft)  dir  ' % (A['t']/TURN)
    print(hdr + '  <rho>@2.0t   <T>@2.0t   nfloor(2.0/2.5/3.0)')
    shmin = [A['sig'].min(axis=(0, 1)) for A in As]
    shdir = []
    for A in As:
        idx = A['sig'].reshape(-1, n1).argmin(axis=0)
        lm = A['lim'].reshape(-1, n1)
        shdir.append(np.array(['r', 'th', 'ph'])[lm[idx, np.arange(n1)]])
    nfl = [((A['rho'] <= 1.05*DFLOOR).sum(axis=(0, 1))) for A in As]
    rmean = As[0]['rho'].mean(axis=(0, 1))
    tmean = As[0]['T'].mean(axis=(0, 1))
    rows = sorted(set(list(range(0, n1, every)) + [n1-1]
                      + [int(np.argmin(s)) for s in shmin]))
    for i in rows:
        mark = '*' if any(i == int(np.argmin(s)) for s in shmin) else ' '
        ln = '%s%3d %6.4f %9.3e ' % (mark, i, rc[i]/RS, As[0]['dr'][i])
        for a in range(len(As)):
            ln += ' %8.4f %-4s' % (shmin[a][i], shdir[a][i])
        ln += ' %9.3e %9.3e  %s' % (rmean[i], tmean[i],
                                    '/'.join(str(f[i]) for f in nfl))
        print(ln)
    return shmin


def shell_pdf(A, shells):
    print('\n--- shell PDFs at t = %.3f turnover ---' % (A['t']/TURN))
    print('   i   r/R    q      min       1%%       median      99%%       max')
    for i in shells:
        for nm, F in (('rho/<rho>', A['rho'][:, :, i]), ('T/<T>', A['T'][:, :, i])):
            f = (F/F.mean()).ravel()
            p = np.percentile(f, [0, 1, 50, 99, 100])
            print('  %3d %6.4f %-9s %9.3e %9.3e %9.3e %9.3e %9.3e'
                  % (i, A['rc'][i]/RS, nm, p[0], p[1], p[2], p[3], p[4]))
        F = A['cs'][:, :, i]
        p = np.percentile(F.ravel(), [0, 1, 50, 99, 100])
        print('  %3d %6.4f %-9s %9.3e %9.3e %9.3e %9.3e %9.3e'
              % (i, A['rc'][i]/RS, 'c_s', p[0], p[1], p[2], p[3], p[4]))
        for nm, F in (('|v_th|', np.abs(A['v'][1][:, :, i])),
                      ('|v_ph|', np.abs(A['v'][2][:, :, i]))):
            p = np.percentile(F.ravel(), [0, 1, 50, 99, 100])
            print('  %3d %6.4f %-9s %9.3e %9.3e %9.3e %9.3e %9.3e'
                  % (i, A['rc'][i]/RS, nm, p[0], p[1], p[2], p[3], p[4]))


if __name__ == '__main__':
    ids = [int(a) for a in sys.argv[1:]] or [4, 5, 6]
    As = [analyse(i) for i in ids]
    for A in As:
        report(A)
    shell_table(As)
    lim_shells = sorted(set(int(np.argmin(A['sig'].min(axis=(0, 1)))) for A in As))
    for A in As:
        shell_pdf(A, lim_shells)
