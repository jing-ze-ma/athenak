#!/usr/bin/env python3
"""RSGPAGE 10-19: per-dump reduction of the RSG g25 cubed-sphere production run (READ ONLY on the run).
usage: reduce.py <run_dir> <out_dir>   writes out_dir/red_NNNNN.npz for every hydro_w+m1 dump pair not yet reduced.
Cubed sphere (gnomonic equiangular, AthenaK src/coordinates/cubed_sphere.hpp): MeshBlocks are written panel-major (gid order,
nmb/6 per panel), lx2/lx3 = block position inside the panel, x2,x3 in [-1,1] -> xi,eta = pi/4 x2,x3; direction on panel p is
(a tan xi + b tan eta + n)/delta.  Velocity: velx = v_r, (vely, velz) on the unit tangents e1, e2 with e1.e2 = -sin xi sin eta.
Contents: shell means, luminosity L(r), photosphere R_ph/T_eff per column (lat-lon map), equatorial + meridional slices
(nearest cell), lat-lon maps of v_r at a few radii."""
import os
import sys
import glob
import time
import numpy as np
sys.path.insert(0, os.environ.get('VISPY', '/raven/u/jinma/ATHENAK/wt-accbig/vis/python'))
import bin_convert as bc  # noqa: E402

RUN, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
A = 7.565733e-15
SIG = 5.670374e-5
RS = 6.957e10
RPH0 = 5.80481e13          # IC photosphere (ic/pass3b)
RMAX_SL = 1.5              # slices kept to r < 1.5 R_ph(IC)
NANG = 720                 # slice angular samples
NLAT, NLON = 180, 360      # lat-lon map
MAPR = (0.07, 0.1, 0.3, 0.6, 0.9, 0.98)   # v_r map radii / R_ph(IC)

FR = {0: ((0, 1, 0), (0, 0, 1), (1, 0, 0)), 1: ((-1, 0, 0), (0, 0, 1), (0, 1, 0)), 2: ((0, -1, 0), (0, 0, 1), (-1, 0, 0)),
      3: ((0, 1, 0), (-1, 0, 0), (0, 0, 1)), 4: ((1, 0, 0), (0, 0, 1), (0, -1, 0)), 5: ((0, 1, 0), (1, 0, 0), (0, 0, -1))}
FA = np.array([FR[p][0] for p in range(6)], float)
FB = np.array([FR[p][1] for p in range(6)], float)
FN = np.array([FR[p][2] for p in range(6)], float)


def assemble(d, var):
    """-> array (6, N, N, nr) indexed [panel, k(eta), j(xi), i]"""
    nmb = len(d['mb_logical'])
    per = nmb // 6
    b0 = d['mb_data'][var][0]
    n3, n2, n1 = b0.shape
    nb = int(round(np.sqrt(per)))
    N = n2 * nb
    P = np.empty((6, N, N, n1), np.float32)
    for m in range(nmb):
        p = m // per
        lx1, lx2, lx3, lev = d['mb_logical'][m]
        P[p, lx3 * n3:(lx3 + 1) * n3, lx2 * n2:(lx2 + 1) * n2, :] = d['mb_data'][var][m]
    return P


def lookup(q, N):
    """unit vectors q (..., 3) -> (panel, k, j) nearest cell indices"""
    dots = q @ FN.T
    p = np.argmax(dots, axis=-1)
    qa = np.einsum('...c,...c->...', q, FA[p]); qb = np.einsum('...c,...c->...', q, FB[p]); qn = np.einsum('...c,...c->...', q, FN[p])
    xi = np.arctan(qa / qn); eta = np.arctan(qb / qn)
    j = np.clip(((xi / (np.pi / 4) + 1) * 0.5 * N).astype(int), 0, N - 1)
    k = np.clip(((eta / (np.pi / 4) + 1) * 0.5 * N).astype(int), 0, N - 1)
    return p, k, j


fh = sorted(glob.glob(os.path.join(os.environ.get('BINDIR') or os.path.join(RUN, 'bin'), '*.hydro_w.*.bin')))
if os.environ.get('ONLY'):
    fh = [x for x in fh if x.split('.')[-2] in os.environ['ONLY'].split(',')]
geo = None
for f in fh:
    idx = f.split('.')[-2]
    fm = f.replace('hydro_w', 'm1')
    o = os.path.join(OUT, 'red_%s.npz' % idx)
    if os.path.exists(o) or not os.path.exists(fm):
        continue
    if time.time() - max(os.path.getmtime(f), os.path.getmtime(fm)) < 180:   # may still be being written
        continue
    t0 = time.time()
    try:
        h = bc.read_binary(f)
        m = bc.read_binary(fm)
    except Exception as e:
        print('skip', idx, e)
        continue
    t = float(h['time'])
    r = np.asarray(h['mb_x1v'][0], float) if 'mb_x1v' in h else None
    rf = np.asarray(h['mb_x1f'][0], float) if 'mb_x1f' in h else None
    if r is None:
        g = np.array([ln.split() for ln in open(os.path.join(RUN, 'rsg.x1grid.txt')) if not ln.startswith('#')], float)
        r = g[3:-3, 2]; rf = np.append(g[3:-3, 1], g[-3, 1])
    rho = assemble(h, 'dens'); vr = assemble(h, 'velx'); v2 = assemble(h, 'vely'); v3 = assemble(h, 'velz')
    ei = assemble(h, 'eint')
    del h
    E = assemble(m, 'm1_e'); F1 = assemble(m, 'm1_f1')
    del m
    N = rho.shape[1]
    if geo is None:
        x2c = -1 + (np.arange(N) + 0.5) * 2.0 / N
        XI, ETA = np.meshgrid(x2c * np.pi / 4, x2c * np.pi / 4, indexing='xy')   # [k, j]: XI varies with j
        X, Y = np.tan(XI), np.tan(ETA)
        dl = np.sqrt(1 + X ** 2 + Y ** 2)
        dOm = (1 + X ** 2) * (1 + Y ** 2) / dl ** 3 * (np.pi / 2 / N) ** 2       # exact-ish cell solid angle
        W = np.broadcast_to(dOm, (6, N, N)).copy()
        cosg = -np.sin(XI) * np.sin(ETA)
        Q = (FA[:, None, None, :] * X[None, :, :, None] + FB[:, None, None, :] * Y[None, :, :, None] + FN[:, None, None, :]) / dl[None, :, :, None]
        geo = 1
        print('cs N %d, sum dOmega/4pi = %.6f' % (N, W.sum() / 4 / np.pi), flush=True)
    Wt = W.sum()
    sm = lambda q: np.tensordot(W, q, axes=([0, 1, 2], [0, 1, 2])) / Wt  # noqa: E731
    rhom = sm(rho)
    mvr = sm(rho * vr) / rhom
    vrrms = np.sqrt(np.maximum(sm(rho * vr ** 2) / rhom - mvr ** 2, 0))
    vt2 = v2 ** 2 + v3 ** 2 + 2 * v2 * v3 * cosg[None, :, :, None]
    vtrms = np.sqrt(sm(rho * vt2) / rhom)
    Em = sm(E); Fr = sm(F1)
    Lr = 4 * np.pi * r ** 2 * Fr
    Trad = (np.maximum(E, 1e-30) / A) ** 0.25
    Tm = sm(Trad)
    eim = sm(ei)
    hent = (5. / 3.) * ei / rho
    Fenth = sm(rho * vr * hent) - sm(rho * vr) * sm(hent)
    Lconv = 4 * np.pi * r ** 2 * Fenth
    lrho = np.log10(rho)
    drho_rms = np.sqrt(sm((lrho - sm(lrho)[None, None, None, :]) ** 2))   # rms of log rho fluctuation per shell
    # photosphere per column (grey Eddington; flux measured at 1.2-1.4 R_ph)
    sel = (r > 1.2 * RPH0) & (r < 1.4 * RPH0)
    Lloc = (F1[..., sel] * r[sel] ** 2).mean(axis=-1)
    i0 = np.searchsorted(r, 0.7 * RPH0)
    G = SIG * Trad[..., i0:] ** 4 * r[i0:] ** 2 - Lloc[..., None]
    neg = G < 0
    n = np.argmax(neg, axis=-1)
    ok = neg.any(axis=-1) & (n > 0)
    n = np.maximum(n, 1)
    ga = np.take_along_axis(G, (n - 1)[..., None], -1)[..., 0]; gb = np.take_along_axis(G, n[..., None], -1)[..., 0]
    ra, rb = r[i0 + n - 1], r[i0 + n]
    Rph = np.where(ok, ra + (rb - ra) * ga / (ga - gb), np.nan)
    Teff = (np.maximum(Lloc, 0) / Rph ** 2 / SIG) ** 0.25
    del G, neg
    # slices
    ns = np.searchsorted(r, RMAX_SL * RPH0)
    ang = (np.arange(NANG) + 0.5) * 2 * np.pi / NANG
    qe = np.stack([np.cos(ang), np.sin(ang), 0 * ang], -1)          # equator (z = 0), angle = phi
    qm = np.stack([np.sin(ang), 0 * ang, np.cos(ang)], -1)          # meridian x-z plane, angle from +z
    out = dict(t=t, r=r, rf=rf, ns=ns, N=N, ang=ang,
               rho=rhom, mvr=mvr, vrrms=vrrms, vtrms=vtrms, E=Em, Trad=Tm, eint=eim, Fr=Fr, Lr=Lr, Lconv=Lconv,
               dlrho_rms=drho_rms,
               Rph_med=np.nanmedian(Rph), Rph_p2=np.nanpercentile(Rph, 2), Rph_p98=np.nanpercentile(Rph, 98),
               Teff_med=np.nanmedian(Teff), Rph_nan=int((~np.isfinite(Rph)).sum()),
               nonfinite=int((~np.isfinite(rho)).sum() + (~np.isfinite(vr)).sum() + (~np.isfinite(E)).sum()),
               vmax=float(np.nanmax(np.abs(vr))), rhomin=float(np.nanmin(rho)))
    lrho0 = None
    for key, qq in (('eq', qe), ('me', qm)):
        p, k, j = lookup(qq, N)
        out[key + '_lrho'] = lrho[p, k, j, :ns].astype(np.float16)
        out[key + '_T'] = (Trad[p, k, j, :ns] / 1e3).astype(np.float16)
        out[key + '_vr'] = (vr[p, k, j, :ns] / 1e5).astype(np.float16)
        out[key + '_vt'] = (np.sqrt(np.maximum(vt2[p, k, j, :ns], 0)) / 1e5).astype(np.float16)
        out[key + '_lE'] = np.log10(np.maximum(E[p, k, j, :ns], 1e-30)).astype(np.float16)
        out[key + '_drho'] = (lrho[p, k, j, :ns] - np.log10(rhom[:ns])[None, :]).astype(np.float16)
    lat = np.linspace(-90, 90, NLAT + 1); lat = 0.5 * (lat[1:] + lat[:-1])
    lon = np.linspace(-180, 180, NLON + 1); lon = 0.5 * (lon[1:] + lon[:-1])
    LA, LO = np.meshgrid(np.radians(lat), np.radians(lon), indexing='ij')
    ql = np.stack([np.cos(LA) * np.cos(LO), np.cos(LA) * np.sin(LO), np.sin(LA)], -1)
    p, k, j = lookup(ql, N)
    out['map_Rph'] = (Rph[p, k, j] / RPH0).astype(np.float32)
    out['map_Teff'] = Teff[p, k, j].astype(np.float32)
    out['map_r'] = np.array(MAPR)
    for x in MAPR:
        i = int(np.searchsorted(r, x * RPH0))
        out['map_vr_%03d' % int(round(100 * x))] = (vr[p, k, j, i] / 1e5).astype(np.float16)
        out['map_drho_%03d' % int(round(100 * x))] = (lrho[p, k, j, i] - np.log10(rhom[i])).astype(np.float16)
    np.savez(o + '.tmp.npz', **out)
    os.replace(o + '.tmp.npz', o)
    print('reduced', idx, 't %.3f d' % (t / 86400), 'R_ph med %.4f R_ph0' % (out['Rph_med'] / RPH0),
          'Teff med %.0f' % out['Teff_med'], 'vr rms max %.3f km/s' % (vrrms.max() / 1e5), 'nonfinite', out['nonfinite'],
          '%.0f s' % (time.time() - t0), flush=True)
