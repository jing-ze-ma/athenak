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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eos_gen  # noqa: E402
TBL = eos_gen.Table()          # rsgE deck <hydro>: X 0.653848 Y 0.326395 a_metal 16, H2 off, ion on, logd -21..-1, logT 3.3..6.8, dlog 0.02
TFLOOR = float(os.environ.get('TFLOOR', '2000'))
VER = 2                        # reduction format (2: T_gas, native rings + shells at float32)

RUN, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
A = 7.565733e-15
SIG = 5.670374e-5
RS = 6.957e10
RPH0 = float(os.environ.get('RPH0', '5.80019e13'))   # IC photosphere (rsgE: ic_rsg_g25_deep.txt header; Raven pass3b 5.80481e13)
AGE = float(os.environ.get('AGE', '120'))           # skip dumps younger than this [s]
RMAX_SL = 1.5              # slices kept to r < 1.5 R_ph(IC)
NANG = 720                 # slice angular samples
NLAT, NLON = 180, 360      # lat-lon map
MAPR = (0.07, 0.11, 0.3, 0.6, 0.9, 0.98, 1.085)   # native-shell radii / R_ph(IC)

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


# BINDIR may be a ':'-separated list (rsgE: dumps move /work/hdd -> /projects after link 1); first dir holding an index wins
_seen = {}
for _bd in (os.environ.get('BINDIR') or os.path.join(RUN, 'bin')).split(':'):
    for _f in sorted(glob.glob(os.path.join(_bd, '*.hydro_w.*.bin'))):
        _seen.setdefault(_f.split('.')[-2], _f)
fh = [_seen[k] for k in sorted(_seen)]
if os.environ.get('ONLY'):
    fh = [x for x in fh if x.split('.')[-2] in os.environ['ONLY'].split(',')]
geo = None
for f in fh:
    idx = f.split('.')[-2]
    fm = f.replace('hydro_w', 'm1')
    o = os.path.join(OUT, 'red_%s.npz' % idx)
    if os.path.exists(o):
        try:
            if int(np.load(o).get('ver', 1)) >= VER:
                continue
        except Exception:
            pass
    if not os.path.exists(fm):
        continue
    try:
        if time.time() - max(os.path.getmtime(f), os.path.getmtime(fm)) < AGE:   # may still be being written
            continue
    except OSError:                                   # dump moved / deleted meanwhile: retry next time
        continue
    t0 = time.time()
    try:
        h = bc.read_binary(f)
        m = bc.read_binary(fm)
    except Exception as e:                            # incl. a dump deleted mid-read: skip, retried next run
        print('skip', idx, e, flush=True)
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
    Tg = np.empty_like(rho); Tlo = np.zeros(rho.shape, bool)
    for pp in range(rho.shape[0]):                    # T_gas(rho, e_int) with the run's own general EOS (eos_gen.py)
        Tg[pp], Tlo[pp] = TBL.temperature(rho[pp], ei[pp])
    Tfl = Tlo | (Tg <= TFLOOR * 1.002)                # at / below the 2000 K floor (or under the table)
    del Tlo
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
    Tgm = sm(Tg)
    Tfloor_frac = sm(Tfl.astype(np.float32))
    # slices: native rings only (no resampling), float32, r < RMAX_SL R_ph(IC)
    ns = int(np.searchsorted(r, RMAX_SL * RPH0))
    out = dict(ver=VER, t=t, r=r, rf=rf, ns=ns, N=N,
               rho=rhom, mvr=mvr, vrrms=vrrms, vtrms=vtrms, E=Em, Trad=Tm, Tgas=Tgm, Tfloor_frac=Tfloor_frac,
               eint=eim, Fr=Fr, Lr=Lr, Lconv=Lconv, dlrho_rms=drho_rms,
               Rph_med=np.nanmedian(Rph), Rph_p2=np.nanpercentile(Rph, 2), Rph_p98=np.nanpercentile(Rph, 98),
               Teff_med=np.nanmedian(Teff), Rph_nan=int((~np.isfinite(Rph)).sum()),
               nonfinite=int((~np.isfinite(rho)).sum() + (~np.isfinite(vr)).sum() + (~np.isfinite(E)).sum()
                             + (~np.isfinite(Tg)).sum()),
               nfloorT=int(Tfl.sum()), vmax=float(np.nanmax(np.abs(vr))), rhomin=float(np.nanmin(rho)))
    # equatorial ring = panels 4,0,1,2 at eta = 0 (mean of the two cell rows next to x3 = 0); b = +z there, so eta = 0 is
    # exactly z = 0 and phi = phi0(panel) + xi.  Meridional ring (y = 0) = panels 5,0,3,2 at xi = 0 (mean of the two
    # columns next to x2 = 0); angle psi = atan2(z, x) = -90+eta, eta, 90+eta, 180-eta (panel 2 reversed).
    h2 = (N // 2 - 1, N // 2)
    xie = (-1 + np.arange(N + 1) * 2.0 / N) * np.pi / 4
    edges = lambda ph0s: np.concatenate([ph0 + xie[(0 if c == 0 else 1):] for c, ph0 in enumerate(ph0s)])  # noqa: E731

    def ring_eq(q):
        return np.concatenate([0.5 * (q[pp, h2[0], :, :ns] + q[pp, h2[1], :, :ns]) for pp in (4, 0, 1, 2)], axis=0)

    def ring_me(q):
        parts = [0.5 * (q[pp, :, h2[0], :ns] + q[pp, :, h2[1], :ns]) for pp in (5, 0, 3, 2)]
        parts[3] = parts[3][::-1]
        return np.concatenate(parts, axis=0)

    out['eqa_phi'] = edges((-np.pi / 2, 0, np.pi / 2, np.pi))
    out['mea_psi'] = edges((-np.pi / 2, 0, np.pi / 2, np.pi))
    for key, ring in (('eqa', ring_eq), ('mea', ring_me)):
        out[key + '_lrho'] = np.log10(ring(rho)).astype(np.float32)
        out[key + '_vr'] = (ring(vr) / 1e5).astype(np.float32)
        out[key + '_Tg'] = ring(Tg).astype(np.float32)
        out[key + '_Tfl'] = (ring(Tfl.astype(np.float32)) > 0.25)
    # native shells (6, N, N) [panel, k(eta), j(xi)] for the Mollweide maps
    out['map_r'] = np.array(MAPR)
    out['sh_Rph'] = (Rph / RPH0).astype(np.float32)
    out['sh_Teff'] = Teff.astype(np.float32)
    for x in MAPR:
        i = int(np.argmin(abs(r - x * RPH0)))
        tag = '%03d' % int(round(100 * x)) if x < 1 else '%04d' % int(round(1000 * x))
        out['sh_vr_' + tag] = (vr[..., i] / 1e5).astype(np.float32)
        out['sh_drho_' + tag] = (lrho[..., i] - np.log10(rhom[i])).astype(np.float32)
        out['sh_dT_' + tag] = (Tg[..., i] / Tgm[i] - 1).astype(np.float32)
        out['sh_Tfl_' + tag] = Tfl[..., i]
    del Tfl
    np.savez(o + '.tmp.npz', **out)
    os.replace(o + '.tmp.npz', o)
    print('reduced', idx, 't %.3f d' % (t / 86400), 'R_ph med %.4f R_ph0' % (out['Rph_med'] / RPH0),
          'Teff med %.0f' % out['Teff_med'], 'T floor cells %d' % out['nfloorT'], 'vr rms max %.3f km/s' % (vrrms.max() / 1e5), 'nonfinite', out['nonfinite'],
          '%.0f s' % (time.time() - t0), flush=True)
