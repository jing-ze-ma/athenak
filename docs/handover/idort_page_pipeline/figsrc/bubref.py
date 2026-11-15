"""idort_page Q2: time-dependent (retarded) and steady free-streaming references for the IDORT multi-source bubbles.

20 opaque discs (IDORT create.py, suiteF_1012/ana/bub_sources.json: x, y, r, T) held at T from t = 0 in vacuum, c = 10.
Exact intensity at (P, t) along the arrival direction u: I = T_i^4 of the FIRST disc hit going backward from P along -u,
if the emission left the disc surface at t_e = t - s/c >= 0 (s = path length to the hit), else 0 (the disc is opaque, so a
disc behind it never contributes).  2-D (in-plane rays, IDORT/our ring set): s = s_inplane.  3-D (cylinders, z-directions,
Gauss sets): s = s_inplane / sin(theta).  Inside a disc E = T^4 (isotropic), F = 0.
E = <I> (a = 1), F = c <I u> (F in units with |F| <= c E).
Boundaries:  'per'  = periodic unit box (images of every disc, as in the suite runs, 5x5 tiles);
             'open' = no images (IDORT: outflow on all boundaries, paper Sect. 4.2.2, Config REFLECTIVE_X/Y = 2).
Times: 0.005 / 0.025 / 0.1 (= 0.05 / 0.25 / 1 t_c, IDORT Fig. 7) and 'inf' (steady).
usage: python bubref.py [n_every]   -> cache/bubref_<every>.npz (sampled on the 512^2 cell centres every <every> cells)
"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

C = 10.0
SRC = np.array(json.load(open('/raven/ptmp/jinma/suiteF_1012/ana/bub_sources.json')))
TS = (0.005, 0.025, 0.1, np.inf)
NPHI, NMU = 2048, 16
PHI = (np.arange(NPHI) + 0.5)*2*np.pi/NPHI
UX, UY = np.cos(PHI), np.sin(PHI)
MU, WMU = np.polynomial.legendre.leggauss(NMU)
ST = np.sqrt(1 - MU**2)


def first_hit(px, py, tiles):
    """px, py: (n,) points. returns s (n, NPHI) in-plane distance to the first disc hit backward along -u, B = T^4"""
    n = px.size
    best = np.full((n, NPHI), np.inf)
    B = np.zeros((n, NPHI))
    for sx in tiles:
        for sy in tiles:
            for (x, y, r, T) in SRC:
                cx = (x + sx - px)[:, None]
                cy = (y + sy - py)[:, None]
                # backward ray P - s u hits |Q - Cc| = r : s^2 + 2 s (c.u) + |c|^2 - r^2 = 0 with c = Cc - P
                # -> near root s = -(c.u) - sqrt((c.u)^2 - |c|^2 + r^2), the disc lying toward -u (c.u < 0)
                cu = cx*UX[None, :] + cy*UY[None, :]
                cc = cx*cx + cy*cy - r*r
                disc = cu*cu - cc
                ok = (disc > 0) & (cu < 0)
                s = np.where(ok, -cu - np.sqrt(np.where(ok, disc, 0.0)), np.inf)
                s = np.where(s >= 0, s, np.inf)
                hit = s < best
                best = np.where(hit, s, best)
                B = np.where(hit, T**4, B)
    return best, B


def chunk(args):
    px, py, bc = args
    tiles = (-2, -1, 0, 1, 2) if bc == 'per' else (0,)
    s, B = first_hit(px, py, tiles)
    nt = len(TS)
    out = {k: np.zeros((nt, px.size)) for k in ('E2', 'F2x', 'F2y', 'E3', 'F3x', 'F3y')}
    inside = np.zeros(px.size, bool)
    Tin = np.zeros(px.size)
    for (x, y, r, T) in SRC:
        dx = px - x
        dy = py - y
        if bc == 'per':
            dx = (dx + 0.5) % 1 - 0.5
            dy = (dy + 0.5) % 1 - 0.5
        m = dx*dx + dy*dy < r*r
        inside |= m
        Tin = np.where(m, T, Tin)
    for k, t in enumerate(TS):
        lim = C*t
        m2 = np.where(s <= lim, B, 0.0)
        out['E2'][k] = m2.mean(1)
        out['F2x'][k] = C*(m2*UX).mean(1)
        out['F2y'][k] = C*(m2*UY).mean(1)
        e3 = np.zeros(px.size)
        f3x = np.zeros(px.size)
        f3y = np.zeros(px.size)
        for a in range(NMU):
            m3 = np.where(s <= lim*ST[a], B, 0.0)
            e3 += WMU[a]/2*m3.mean(1)
            f3x += WMU[a]/2*ST[a]*C*(m3*UX).mean(1)
            f3y += WMU[a]/2*ST[a]*C*(m3*UY).mean(1)
        out['E3'][k], out['F3x'][k], out['F3y'][k] = e3, f3x, f3y
    for key in out:
        if key.startswith('E'):
            out[key][:, inside] = Tin[inside]**4
        else:
            out[key][:, inside] = 0.0
    return out


def make(every=4, nproc=16):
    n = 512
    xs = -0.5 + (np.arange(n) + 0.5)/n
    sub = xs[every//2::every]
    X, Y = np.meshgrid(sub, sub)          # [j, i] = (y, x)
    px, py = X.ravel(), Y.ravel()
    res = {}
    with Pool(nproc) as pool:
        for bc in ('per', 'open'):
            idx = np.array_split(np.arange(px.size), 4*nproc)
            parts = pool.map(chunk, [(px[i], py[i], bc) for i in idx])
            for key in parts[0]:
                res['%s_%s' % (key, bc)] = np.concatenate([p[key] for p in parts], axis=1).reshape(len(TS), sub.size, sub.size)
            print(bc, 'done', flush=True)
    fn = '/raven/ptmp/jinma/idort_page/cache/bubref_%d.npz' % every
    np.savez(fn, x=sub, ts=np.array(TS), **res)
    print('wrote', fn)


if __name__ == '__main__':
    make(int(sys.argv[1]) if len(sys.argv) > 1 else 4, int(sys.argv[2]) if len(sys.argv) > 2 else 16)
