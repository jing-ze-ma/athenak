#!/usr/bin/env python3
"""he_remap_rst.py: remap an AthenaK restart file onto a finer angular grid (x2, x3 by an
integer factor: 1 = identity, 2 = 64^2 -> 128^2) and, optionally, onto a new radial
StretchRPoly grid (same x1min, x1max, ng); writes a restart file the code reads with -r.
Written for the he_star_m1 spherical-polar wedge (single level, x2 AND x3 periodic, one
MeshBlock along x1, general EOS, implicit M1).  The parser follows src/outputs/restart.cpp
(rt-integration 69e8f149) and refuses blocks it does not know.

  he_remap_rst.py IN.rst OUT.rst --factor 2 --mbx2 32 --mbx3 32
                  [--grid NEW.npy --nx1 N]  [--keep-caches] [--ghosts wrap|source]
                  [--dt-factor 0.5] [--force-math] [--no-slopes]

RADIAL (--grid: 18 StretchRPoly parameters, written into the embedded input as
f_stretch_r_*): conservative overlap remap of a minmod-limited linear reconstruction in
V = r^3/3; see radial_remap().  The first link must use time/restart_refill_ghosts=true.
ANGULAR, per coarse cell (exact conservation to round-off):
  hydro: rho linear (minmod, volume centroids); v_i and h = E/rho - v^2/2 linear about the
      MASS centroid; a uniform shift of h restores E (it takes the sub-cell kinetic-energy
      excess out of the internal energy).  A coarse cell whose child would get e_int <
      0.5 x the parent's specific e_int is injected (children = parent).
  rad_m1 E, F: conservative linear; injected where a child would get |F| > c E.
  f0x1: linear in angle (face areas); f0x2/f0x3: fine faces on coarse faces copy them,
      mid faces take the mean of the two, then linear along the other angle.
  pgen block (he_star_m1 adaptive-MLT F_sub, radial): copied (radially: r^2 F linear).
  ctr_mem (thin-relax closure memory): injected.
DROPPED unless --keep-caches: the general-EOS caches (wtemp, wder, eint), the M1 step
predictor and the rad_signal_speed inputs.  The reader accepts files without them and
rebuilds them from u0.  Injecting them would be wrong: with both EOS caches present the
first ConsToPrim is frozen (Hydro::c2p_freeze_derived) and every child would keep its
parent's pressure.  --keep-caches exists for the identity gate only.
GHOSTS: --ghosts wrap fills the angular ghosts by the periodic wrap of the new global
array; --ghosts source (identity only) copies the source ghosts.
HEADER: nx1/nx2/nx3 of <mesh> and <meshblock>, RegionSize, RegionIndcs, nmb_total,
root_level and the Z-ordered logical locations are rewritten; time and ncycle are kept.
A restart's first cycle takes the header dt as is (driver.cpp) and dt = 0 is fatal, so
the header dt is --dt-factor x the source dt; the next cycle may at most double it.
"""
import argparse
import re
import struct
import sys

import numpy as np

MAGIC = {b'PGENST01': 'pgen', b'RTWARM01': 'warm', b'M1PRED01': 'pred',
         b'M1TIME2A': 't2', b'M1ONEP01': 'onep', b'M1MRWIN1': 'mr', b'EINTRST1': 'eint',
         b'CKRST001': 'ck', b'M1CTRLX1': 'ctr', b'M1RSSIG1': 'rss'}
NIND = 19   # ints in RegionIndcs (mesh.hpp)


# ------------------------------------------------------------------------------- reading
class Rst:
    pass


def read_rst(fn):
    r = Rst()
    raw = np.memmap(fn, dtype=np.uint8, mode='r')
    h = bytes(raw[:400000])
    loc = h.find(b'<par_end>')
    if loc < 0:
        sys.exit('no <par_end>')
    o = loc + 10
    r.text = h[:o].decode()
    r.nmb, r.root = struct.unpack_from('<ii', h, o)
    o += 8
    r.msize = list(struct.unpack_from('<9d', h, o))
    o += 72
    r.mind = list(struct.unpack_from('<%di' % NIND, h, o))
    o += 4 * NIND
    r.bind = list(struct.unpack_from('<%di' % NIND, h, o))
    o += 4 * NIND
    r.time, r.dt = struct.unpack_from('<2d', h, o)
    o += 16
    r.ncycle, = struct.unpack_from('<i', h, o)
    o += 4
    r.lloc = np.array(struct.unpack_from('<%di' % (5 * r.nmb), h, o)).reshape(r.nmb, 5)
    o += 20 * r.nmb
    r.cost = np.array(struct.unpack_from('<%df' % r.nmb, h, o), dtype=np.float32)
    o += 4 * r.nmb
    r.blocks = []           # (magic bytes, payload bytes), in file order
    while h[o:o + 8] in MAGIC:
        nb, = struct.unpack_from('<Q', h, o + 8)
        r.blocks.append((h[o:o + 8], h[o + 16:o + 16 + nb]))
        o += 16 + nb
    if h[o:o + 8].isalnum():
        sys.exit('unknown marked block %r' % h[o:o + 8])
    r.dsize, = struct.unpack_from('<Q', h, o)
    o += 8
    if o + r.nmb * r.dsize != raw.size:
        sys.exit(
            'size mismatch: header %d + %d x %d != %d' %
            (o, r.nmb, r.dsize, raw.size))
    r.data0 = o
    ng, n1, n2, n3 = r.bind[0], r.bind[1], r.bind[2], r.bind[3]
    r.ng = ng
    r.o1, r.o2, r.o3 = n1 + 2 * ng, n2 + 2 * ng, n3 + 2 * ng
    if n2 <= 1 or n3 <= 1:
        sys.exit('needs a 3-D mesh')
    hd = {MAGIC[m]: p for m, p in r.blocks}
    for k in ('warm', 't2', 'onep', 'mr', 'ck'):
        if k in hd:
            sys.exit('block %s not supported' % k)
    c = r.o1 * r.o2 * r.o3
    r.npred = int(np.frombuffer(hd['pred'][:8], '<i4')[1]) if 'pred' in hd else 0
    r.neint = int(np.frombuffer(hd['eint'][:8], '<i4')[0]) if 'eint' in hd else 0
    r.nctr = int(np.frombuffer(hd['ctr'][:8], '<i4')[1]) if 'ctr' in hd else 0
    r.nrss = int(np.frombuffer(hd['rss'][:8], '<i4')[1]) if 'rss' in hd else 0
    f1 = (r.o1 + 1) * r.o2 * r.o3
    f2 = r.o1 * (r.o2 + 1) * r.o3
    f3 = r.o1 * r.o2 * (r.o3 + 1)
    nr = r.dsize // 8
    # layout (restart.cpp): hydro | m1 u0 | f0x1 | f0x2 | f0x3 | wtemp | wder(2) | pred |
    # eint | ctr | rss   (general EOS, implicit transverse M1; no mhd/rad/turb/z4c here)
    ns = get_param(r.text, 'hydro', 'nscalars')
    r.nhyd = 5 + (int(ns) if ns else 0)
    rest = nr - f1 - f2 - f3 - c * (r.nhyd + 4 + r.npred + r.neint + r.nctr + r.nrss)
    if rest not in (0, 3 * c):
        sys.exit('cannot place the record layout (rest %d, c %d)' % (rest, c))
    r.has_wt = rest == 3 * c
    r.sec = [('hyd', r.nhyd * c, (r.nhyd, r.o3, r.o2, r.o1)),
             ('m1', 4 * c, (4, r.o3, r.o2, r.o1)),
             ('f1', f1, (r.o3, r.o2, r.o1 + 1)),
             ('f2', f2, (r.o3, r.o2 + 1, r.o1)),
             ('f3', f3, (r.o3 + 1, r.o2, r.o1))]
    if r.has_wt:
        r.sec += [('wt', c, (1, r.o3, r.o2, r.o1)), ('wd', 2 * c, (2, r.o3, r.o2, r.o1))]
    for nm, n in (('pred', r.npred), ('eint', r.neint), ('ctr', r.nctr), ('rss', r.nrss)):
        if n:
            r.sec.append((nm, n * c, (n, r.o3, r.o2, r.o1)))
    assert sum(s[1] for s in r.sec) == nr
    r.raw = raw
    return r


def mb_data(r, m):
    """dict of the MeshBlock record m's sections (float64 views)"""
    a = np.frombuffer(
        r.raw,
        dtype='<f8',
        count=r.dsize //
        8,
        offset=r.data0 +
        m *
        r.dsize)
    d, o = {}, 0
    for nm, n, shp in r.sec:
        d[nm] = a[o:o + n].reshape(shp)
        o += n
    return d


def gather(r):
    """global ACTIVE-angle arrays (radial ghosts included): cells (n, N3, N2, o1),
    x1 faces (N3, N2, o1+1), x2 faces (lower faces, periodic) (N3, N2, o1), x3 faces
    (N3, N2, o1)"""
    ng = r.ng
    b2, b3 = r.bind[2], r.bind[3]
    N2, N3 = r.mind[2], r.mind[3]
    G = {}
    for nm, n, shp in r.sec:
        if nm in ('f1',):
            G[nm] = np.empty((N3, N2, r.o1 + 1))
        elif nm in ('f2', 'f3'):
            G[nm] = np.empty((N3, N2, r.o1))
        else:
            G[nm] = np.empty((shp[0], N3, N2, r.o1))
    for m in range(r.nmb):
        lx2, lx3 = r.lloc[m, 1], r.lloc[m, 2]
        d = mb_data(r, m)
        J, K = slice(lx2 * b2, (lx2 + 1) * b2), slice(lx3 * b3, (lx3 + 1) * b3)
        a2, a3 = slice(ng, ng + b2), slice(ng, ng + b3)
        for nm in G:
            if nm == 'f1':
                G[nm][K, J] = d[nm][a3, a2]
            elif nm in ('f2', 'f3'):
                G[nm][K, J] = d[nm][a3, a2]       # lower faces of the active cells
            else:
                G[nm][:, K, J] = d[nm][:, a3, a2]
    return G


# ------------------------------------------------------------------------------- remap
def minmod(a, b):
    return np.where(a * b > 0.0, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)


def slopes(q, ax):
    """minmod slope per coarse cell along axis ax (periodic)"""
    return minmod(np.roll(q, -1, axis=ax) - q, q - np.roll(q, 1, axis=ax))


def child_geom(x2min, x2max, N2, f):
    """theta children of each coarse cell: volume fractions W (N2, f) and offsets xi
    (N2, f) of the child volume centroids from the coarse volume centroid, in units of the
    coarse dtheta.  f = 1 gives W = 1, xi = 0 exactly."""
    e = np.linspace(x2min, x2max, N2 * f + 1)
    lo, hi = e[:-1], e[1:]
    w = np.cos(lo) - np.cos(hi)
    cen = (np.sin(hi) - hi * np.cos(hi) - np.sin(lo) + lo * np.cos(lo)) / w
    w = w.reshape(N2, f)
    cen = cen.reshape(N2, f)
    W = w / w.sum(axis=1, keepdims=True)
    cc = (W * cen).sum(axis=1, keepdims=True)
    xi = (cen - cc) / ((x2max - x2min) / N2)
    if f == 1:
        W[:] = 1.0
        xi[:] = 0.0
    return W, xi


def prolong_cells(q, s2, s3, xi2, xi3, f):
    """q, s2, s3: (N3, N2, R) coarse value and slopes; xi2 (N2, f), xi3 (f,): offsets.
    Returns (N3*f, N2*f, R) with child (K, a, J, b) = q + s2 xi2[J, b] + s3 xi3[a]."""
    N3, N2, R = q.shape
    out = (q[:, None, :, None, :] + s2[:, None, :, None, :] * xi2[None, None, :, :, None]
           + s3[:, None, :, None, :] * xi3[None, :, None, None, None])
    return out.reshape(N3 * f, N2 * f, R)


def inject(q, f):
    N3, N2, R = q.shape
    return np.broadcast_to(q[:, None, :, None, :],
                           (N3, f, N2, f, R)).reshape(N3 * f, N2 * f, R)


def coarse_sum(x, W2, f):
    """volume-weighted sum over the f x f children (fine (N3f, N2f, R) -> (N3, N2, R));
    W2 (N2, f) theta weights, phi weights 1/f"""
    N3f, N2f, R = x.shape
    N3, N2 = N3f // f, N2f // f
    y = x.reshape(N3, f, N2, f, R) * W2[None, None, :, :, None] / f
    return y.sum(axis=(1, 3))


def remap(r, G, f, clight, use_slopes=True, force_math=False):
    N2, N3 = r.mind[2], r.mind[3]
    W2, xi2 = child_geom(r.msize[1], r.msize[4], N2, f)
    xi3 = (np.arange(f) + 0.5) / f - 0.5
    stats = {}
    S = 1.0 if use_slopes else 0.0
    out = {}
    # ---- hydro
    u = G['hyd']
    rho = u[0]
    v = [u[1 + n] / rho for n in range(3)]
    h = u[4] / rho - 0.5 * (v[0]**2 + v[1]**2 + v[2]**2)
    rho_f = prolong_cells(rho, S * slopes(rho, 1), S * slopes(rho, 0), xi2, xi3, f)
    # mass centroid offsets per coarse cell
    m_c = coarse_sum(rho_f, W2, f)                                # = rho to round-off
    xi2_f = xi2.reshape(1, N2 * f, 1)
    xi3_f = np.tile(xi3, N3).reshape(N3 * f, 1, 1)
    xm2 = coarse_sum(rho_f * xi2_f, W2, f) / m_c
    xm3 = coarse_sum(rho_f * xi3_f, W2, f) / m_c

    def spec(q):
        s2, s3 = S * slopes(q, 1), S * slopes(q, 0)
        base = inject(q - s2 * xm2 - s3 * xm3, f)
        return base + inject(s2, f) * xi2_f + inject(s3, f) * xi3_f

    v_f = [spec(x) for x in v]
    h_f = spec(h)
    ke_f = 0.5 * (v_f[0]**2 + v_f[1]**2 + v_f[2]**2)
    E_c = u[4]
    dh = (E_c - coarse_sum(rho_f * (h_f + ke_f), W2, f)) / m_c
    h_f = h_f + inject(dh, f)
    # e_int guard: Phi from the file's e_int cache when present
    bad = np.zeros_like(rho, dtype=bool)
    if 'eint' in G:
        eint_c = G['eint'][0]
        phi = h - eint_c / rho
        se_f = h_f - inject(phi, f)
        sec = inject(eint_c / rho, f)
        badf = (se_f < 0.5 * sec) | (rho_f <= 0.0)
        stats['phi_check'] = phi
    else:
        badf = rho_f <= 0.0
    bad |= badf.reshape(N3, f, N2, f, -1).any(axis=(1, 3))
    if f == 1 and not force_math:
        bad[:] = True          # identity: plain copy (bitwise)
    stats['hyd_fallback'] = int(bad.sum())
    stats['hyd_fb_i'] = bad.sum(axis=(0, 1))
    B = inject(bad, f)
    hyd = np.empty((u.shape[0], N3 * f, N2 * f, u.shape[3]))
    hyd[0] = np.where(B, inject(rho, f), rho_f)
    for n in range(3):
        hyd[1 + n] = np.where(B, inject(u[1 + n], f), rho_f * v_f[n])
    hyd[4] = np.where(B, inject(u[4], f), rho_f * (h_f + ke_f))
    for n in range(5, u.shape[0]):      # passive scalars: specific, mass centroid
        hyd[n] = np.where(B, inject(u[n], f), rho_f * spec(u[n] / rho))
    out['hyd'] = hyd
    # ---- rad_m1 u0
    w = G['m1']
    m1 = np.empty((4, N3 * f, N2 * f, w.shape[3]))
    for n in range(4):
        m1[n] = prolong_cells(w[n], S * slopes(w[n], 1), S * slopes(w[n], 0), xi2, xi3, f)
    fmag = np.sqrt(m1[1]**2 + m1[2]**2 + m1[3]**2)
    badf = (m1[0] <= 0.0) | (fmag > clight * m1[0] * (1.0 + 1e-12))
    badm = badf.reshape(N3, f, N2, f, -1).any(axis=(1, 3))
    if f == 1 and not force_math:
        badm[:] = True
    stats['m1_fallback'] = int(badm.sum())
    stats['m1_fb_i'] = badm.sum(axis=(0, 1))
    stats['m1_parent_fc'] = float((np.sqrt(w[1]**2 + w[2]**2 + w[3]**2) /
                                   (clight * w[0])).max())
    Bm = inject(badm, f)
    for n in range(4):
        m1[n] = np.where(Bm, inject(w[n], f), m1[n])
    out['m1'] = m1
    # ---- face fluxes
    q = G['f1']
    out['f1'] = prolong_cells(q, S * slopes(q, 1), S * slopes(q, 0), xi2, xi3, f)
    q = G['f2']                              # (N3, N2 faces, R)
    if f == 1:
        out['f2'] = q.copy()
    else:
        qf = np.empty((N3, N2 * f, q.shape[2]))
        for b in range(f):
            qf[:, b::f] = q + (b / f) * (np.roll(q, -1, axis=1) - q)
        s3 = S * slopes(qf, 0)
        out['f2'] = (qf[:, None] + s3[:, None] * xi3[None, :, None, None]).reshape(
            N3 * f, N2 * f, q.shape[2])
    q = G['f3']                              # (N3 faces, N2, R)
    if f == 1:
        out['f3'] = q.copy()
    else:
        qf = np.empty((N3 * f, N2, q.shape[2]))
        for a in range(f):
            qf[a::f] = q + (a / f) * (np.roll(q, -1, axis=0) - q)
        s2 = S * slopes(qf, 1)
        # x3 faces: area element r dr dtheta, uniform
        xe = (np.arange(f) + 0.5) / f - 0.5
        out['f3'] = (qf[:, :, None, :] + s2[:, :, None, :] *
                     xe[None, None, :, None]).reshape(N3 * f, N2 * f, q.shape[2])
    # ---- injected (ctr memory; caches only with --keep-caches)
    for nm in ('ctr', 'wt', 'wd', 'pred', 'eint', 'rss'):
        if nm in G:
            out[nm] = np.stack([inject(G[nm][n], f) for n in range(G[nm].shape[0])])
    return out, stats


# ------------------------------------------------------------------------------- radial
# the code's StretchRPoly (mesh: use_grid_stretch_r_poly; 8 polynomial coefficients, two
# tanh bumps and the lch plateau); ghost faces = the same map at xi outside [0, 1]
SKEYS = (['f_stretch_r_c%d' % k for k in range(1,
                                               9)] + ['f_stretch_r_b1_amp',
                                                      'f_stretch_r_b1_x',
                                                      'f_stretch_r_b1_w',
                                                      'f_stretch_r_b2_amp',
                                                      'f_stretch_r_b2_x',
                                                      'f_stretch_r_b2_w',
                                                      'f_stretch_r_p_amp',
                                                      'f_stretch_r_p_xa',
                                                      'f_stretch_r_p_xb',
                                                      'f_stretch_r_p_w'])


def _lch(z):
    a = np.abs(z)
    return a + np.log1p(np.exp(-2 * a)) - np.log(2.0)


def _gpl(xa, xb, w, xi):
    return 0.5 * w * (_lch((xi - xa) / w) - _lch(-xa / w) -
                      _lch((xi - xb) / w) + _lch(-xb / w))


def ufun(p, xi):
    u = xi.copy()
    xik = xi.copy()
    for k in range(8):
        u = u + p[k] * xik * (1 - xi)
        xik = xik * xi
    for b in range(2):
        a, xb, w = p[8 + 3 * b:11 + 3 * b]
        if a != 0:
            u = u + a * w * (np.tanh((xi - xb) / w) - (1 - xi) *
                             np.tanh(-xb / w) - xi * np.tanh((1 - xb) / w))
    a, xa, xb, w = p[14:18]
    if a != 0:
        u = u + a * (_gpl(xa, xb, w, xi) - xi * _gpl(xa, xb, w, 1.0))
    return u


def get_param(text, block, par):
    m = re.search(r'<%s>\n(?:(?!<).*\n)*?%s\s*=\s*(\S+)' %
                  (re.escape(block), re.escape(par)), text)
    return None if m is None else m.group(1)


def faces_of(p, n1, ng, x1min, x1max):
    xi = np.arange(-ng, n1 + ng + 1) / float(n1)
    return x1min + (x1max - x1min) * ufun(np.asarray(p, dtype=float), xi)


def overlap_mats(fo, fn):
    """A, B (n_new x n_old) in V = r^3/3: the content of new cell n of a field that
    is q_o + s_o (V - Vc_o) in old cell o is sum_o q_o A + s_o B"""
    Vo, Vn = fo**3 / 3.0, fn**3 / 3.0
    no, nn = Vo.size - 1, Vn.size - 1
    A = np.zeros((nn, no))
    B = np.zeros((nn, no))
    Vc = 0.5 * (Vo[1:] + Vo[:-1])
    for n in range(nn):
        a, b = Vn[n], Vn[n + 1]
        lo = max(np.searchsorted(Vo, a, side='right') - 1, 0)
        hi = min(np.searchsorted(Vo, b, side='left'), no)
        for o in range(lo, hi):
            x0, x1 = max(a, Vo[o]), min(b, Vo[o + 1])
            if x1 > x0:
                A[n, o] = x1 - x0
                B[n, o] = 0.5 * ((x1 - Vc[o])**2 - (x0 - Vc[o])**2)
    return A, B, Vo, Vn


def radial_remap(r, G, pnew, n1new, use_slopes=True):
    """conservative radial remap of the gathered state onto a new StretchRPoly grid (same
    x1min, x1max, ng).  Active cells: overlap integration of a minmod-limited linear
    reconstruction in V = r^3/3 (exact conservation per column).  An OLD cell's slope is
    dropped (all hydro fields, or all M1 fields) wherever a new cell it overlaps would get
    e_int < 0.5 x the piecewise-constant remap of e_int, or |F| > c E: exact conservation,
    and positivity by Cauchy-Schwarz.  Conserved: rho, m_r, r m_theta, r m_phi (L_z), U =
    E - rho Phi (thermal + kinetic; E rebuilt with the NEW cells' Phi = a + b/x1v fit to
    the file), E_rad, F_r, r F_theta, r F_phi.  x1-face fluxes f0x1: r^2 f linear in r;
    x2/x3 face fluxes: linear in r at the new centroids; ctr_mem: nearest cell.  Radial
    GHOSTS: the old ghost values in order -- the first link must refill them
    (time/restart_refill_ghosts=true).  Returns the new G; updates r in place."""
    ng = r.ng
    n1 = r.bind[1]
    x1min, x1max = r.msize[0], r.msize[3]
    pold = [float(get_param(r.text, 'mesh', k) or 0.0) for k in SKEYS]
    fo = faces_of(pold, n1, ng, x1min, x1max)
    pnew = [float('%.12e' % v) for v in pnew]    # exactly what the input file will say
    fn = faces_of(pnew, n1new, ng, x1min, x1max)
    if np.any(np.diff(fn) <= 0):
        sys.exit('new grid not monotone')
    xo = 0.75 * (fo[1:]**4 - fo[:-1]**4) / (fo[1:]**3 - fo[:-1]**3)
    xn = 0.75 * (fn[1:]**4 - fn[:-1]**4) / (fn[1:]**3 - fn[:-1]**3)
    ao, an = slice(ng, ng + n1), slice(ng, ng + n1new)
    A, B, Vo, Vn = overlap_mats(fo[ng:ng + n1 + 1], fn[ng:ng + n1new + 1])
    dVo = np.diff(Vo)
    dVn = np.diff(Vn)
    Vc = 0.5 * (Vo[1:] + Vo[:-1])
    S = 1.0 if use_slopes else 0.0
    st = {}
    Aov = (A > 0.0)

    def slope_of(q):
        qa = q[..., ao]
        dq = np.diff(qa, axis=-1) / np.diff(Vc)
        s = np.zeros_like(qa)
        s[..., 1:-1] = minmod(dq[..., 1:], dq[..., :-1])
        # one-sided at the inner end (the fallback keeps it safe); none at the outer end
        # (floor region: no extrapolation below the floor)
        s[..., 0] = dq[..., 0]
        return S * s

    def remap1(q, s=None, zero=None):
        """q (..., o1 old) -> (..., o1 new): active by overlap, ghosts copied; s = slopes
        (computed when None), zero = OLD active cells whose slope is dropped"""
        qa = q[..., ao]
        if s is None:
            s = slope_of(q)
        if zero is not None:
            s = np.where(zero, 0.0, s)
        cont = qa @ A.T + s @ B.T
        out = np.empty(q.shape[:-1] + (n1new + 2 * ng,))
        out[..., an] = cont / dVn
        out[..., :ng] = q[..., :ng]
        out[..., ng + n1new:] = q[..., ng + n1:]
        return out

    def old_of(badn):
        return (badn.astype(float) @ Aov.astype(float)) > 0.0

    def tot(q, dv, sl):
        return (q[..., sl] * dv).sum(axis=-1)
    # Phi = a + b/x1v from the file (the etotgrav energy carries rho Phi)
    u = G['hyd']
    rho = u[0]
    ke = 0.5 * (u[1]**2 + u[2]**2 + u[3]**2) / rho
    phi = (u[4] - ke - G['eint'][0]) / rho
    pm = np.median(phi, axis=(0, 1))
    Am = np.vstack([np.ones_like(xo), 1.0 / xo]).T
    cf = np.linalg.lstsq(Am[ao], pm[ao], rcond=None)[0]
    st['phi_fit'] = (cf, float(np.abs((Am @ cf - pm) / pm)[ao].max()))
    phin = cf[0] + cf[1] / xn
    U = u[4] - rho * phi
    qs = [rho, u[1], u[2] * xo, u[3] * xo, U] + [u[n] for n in range(5, u.shape[0])]
    ss = [slope_of(q) for q in qs]
    eint_r = remap1(G['eint'][0], s=np.zeros_like(G['eint'][0][..., ao]))
    zero = np.zeros(rho[..., ao].shape, dtype=bool)
    for it in range(8):
        res = [remap1(q, s, zero) for q, s in zip(qs, ss)]
        rn, mr, mt, mp, Un = res[:5]
        mt, mp = mt / xn, mp / xn
        ken = 0.5 * (mr**2 + mt**2 + mp**2) / rn
        badn = ((Un - ken) < 0.5 * eint_r)[..., an] | (rn[..., an] <= 0.0)
        if not badn.any():
            break
        zero |= old_of(badn)
    st['rad_hyd_slope_dropped_old_cells'] = int(zero.sum())
    st['rad_hyd_bad_after'] = int(badn.sum())
    hyd = np.empty(u.shape[:3] + (n1new + 2 * ng,))
    hyd[0], hyd[1], hyd[2], hyd[3] = rn, mr, mt, mp
    for n in range(5, u.shape[0]):
        hyd[n] = res[n]
    eint_n = Un - ken
    st['rad_eint_min_ratio'] = float((eint_n / eint_r)[..., an].min())
    hyd[4] = Un + hyd[0] * phin
    hyd[4][..., :ng] = u[4][..., :ng]
    hyd[4][..., ng + n1new:] = u[4][..., ng + n1:]
    st['rad_rho_min'] = float(hyd[0][..., an].min())
    for nm, qo, qn in (('mass', rho, hyd[0]), ('U', U, Un), ('m_r', u[1], hyd[1])):
        to, tn = tot(qo, dVo, ao), tot(qn, dVn, an)
        st['rad_cons_' + nm] = float((np.abs(tn - to) / np.abs(to).max()).max())
    Gn = {'hyd': hyd, 'eint': eint_n[None]}
    w = G['m1']
    m = re.search(r'<rad_m1>\n(?:(?!<).*\n)*?c_light\s*=\s*(\S+)', r.text)
    cl = float(m.group(1))
    qs = [w[0], w[1], w[2] * xo, w[3] * xo]
    ss = [slope_of(q) for q in qs]
    zero = np.zeros(w[0][..., ao].shape, dtype=bool)
    for it in range(6):
        e0, f1, f2, f3 = [remap1(q, s, zero) for q, s in zip(qs, ss)]
        f2, f3 = f2 / xn, f3 / xn
        badn = ((np.sqrt(f1**2 + f2**2 + f3**2) > cl *
                e0 * (1.0 + 1e-12)) | (e0 <= 0.0))[..., an]
        if not badn.any():
            break
        zero |= old_of(badn)
    st['rad_m1_slope_dropped_old_cells'] = int(zero.sum())
    m1 = np.stack([e0, f1, f2, f3])
    to, tn = tot(w[0], dVo, ao), tot(m1[0], dVn, an)
    st['rad_cons_Erad'] = float((np.abs(tn - to) / to).max())
    fm = np.sqrt(m1[1]**2 + m1[2]**2 + m1[3]**2)
    over = fm > cl * m1[0]
    st['rad_F_max_ratio_before_clip'] = float((fm / (cl * m1[0]))[..., an].max())
    st['rad_F_clipped'] = int(over[..., an].sum())
    fac = np.where(over, cl * m1[0] / np.maximum(fm, 1e-300), 1.0)
    for n in (1, 2, 3):
        m1[n] = m1[n] * fac
    Gn['m1'] = m1
    # face fluxes
    g = G['f1'] * fo[None, None, :]**2
    Gn['f1'] = np.stack([np.stack([np.interp(fn, fo, g[k, j]) for j in range(g.shape[1])])
                         for k in range(g.shape[0])]) / fn[None, None, :]**2
    for nm in ('f2', 'f3'):
        q = G[nm]
        Gn[nm] = np.apply_along_axis(lambda c: np.interp(xn, xo, c), -1, q)
    near = np.abs(xn[:, None] - xo[None, :]).argmin(axis=1)
    for nm in ('ctr', 'wt', 'wd', 'pred', 'rss'):
        if nm in G:
            Gn[nm] = G[nm][..., near]
    # header / text / pgen state
    r.text = set_param(r.text, 'mesh', 'nx1', n1new)
    r.text = set_param(r.text, 'meshblock', 'nx1', n1new)
    for k, v in zip(SKEYS, pnew):
        r.text = set_param(r.text, 'mesh', k, '%.12e' % v)
    r.msize[6] = (x1max - x1min) / n1new
    for ind in (r.mind, r.bind):
        ind[1] = n1new
        ind[5] = ng + n1new - 1
        if ind[10]:
            ind[10] = n1new // 2
            ind[14] = ng + n1new // 2 - 1
    r.o1 = n1new + 2 * ng
    blocks = []
    for mg, pl in r.blocks:
        if MAGIC[mg] == 'pgen':
            hdr = np.frombuffer(pl[:8], '<i4')
            nfc = int(hdr[0])
            if nfc != fo.size:
                sys.exit('pgen state: %d faces, grid has %d' % (nfc, fo.size))
            v = np.frombuffer(pl[8:], '<f8')
            tl, fs = v[0], v[1:]
            fsn = np.interp(fn, fo, fs * fo**2) / fn**2
            pl = (np.array([fn.size, int(hdr[1])], '<i4').tobytes()
                  + np.array([tl], '<f8').tobytes() + fsn.astype('<f8').tobytes())
            st['pgen_fsub_max'] = float(np.abs(fs).max())
        blocks.append((mg, pl))
    r.blocks = blocks
    r.fo, r.fn = fo, fn
    return Gn, st


# ------------------------------------------------------------------------------- sector
def extract_sector(r, G, x2s, x2e, x3s, x3e):
    """cut the angular sub-sector [x2s, x2e] x [x3s, x3e] (on source cell edges) out of
    the gathered state, WITH a one-cell halo of true source neighbours (periodic wrap of
    the source), so that the angular prolongation's limited slopes at the sector edges use
    the real neighbours; strip_sector() removes the halo after remap().  Updates r: the
    geometry seen by remap() is the halo'd box; r.sector holds the final one."""
    N2, N3 = r.mind[2], r.mind[3]
    d2 = (r.msize[4] - r.msize[1])/N2
    d3 = (r.msize[5] - r.msize[2])/N3
    j0 = (x2s - r.msize[1])/d2
    j1 = (x2e - r.msize[1])/d2
    k0 = (x3s - r.msize[2])/d3
    k1 = (x3e - r.msize[2])/d3
    idx = []
    for v in (j0, j1, k0, k1):
        if abs(v - round(v)) > 1e-8:
            sys.exit('sector bounds are not on source cell edges (%.9f)' % v)
        idx.append(int(round(v)))
    j0, j1, k0, k1 = idx
    if not (0 <= j0 < j1 <= N2 and 0 <= k0 < k1 <= N3):
        sys.exit('sector outside the source mesh')
    J = np.arange(j0 - 1, j1 + 1) % N2
    K = np.arange(k0 - 1, k1 + 1) % N3
    Gs = {}
    for nm, a in G.items():
        if a.ndim == 3:
            Gs[nm] = a[np.ix_(K, J)]
        else:
            Gs[nm] = a[:, K][:, :, J]
    r.sector = (x2s, x2e, x3s, x3e, j1 - j0, k1 - k0)
    r.msize[1], r.msize[4] = x2s - d2, x2e + d2
    r.msize[2], r.msize[5] = x3s - d3, x3e + d3
    r.mind[2], r.mind[3] = j1 - j0 + 2, k1 - k0 + 2
    return Gs


def strip_sector(r, out, f):
    """remove the prolonged halo (f fine cells per side) and set the sector geometry"""
    x2s, x2e, x3s, x3e, n2, n3 = r.sector
    o = {}
    for nm, a in out.items():
        if a.ndim == 3:
            o[nm] = np.ascontiguousarray(a[f:f + n3*f, f:f + n2*f])
        else:
            o[nm] = np.ascontiguousarray(a[:, f:f + n3*f, f:f + n2*f])
    r.msize[1], r.msize[4], r.msize[2], r.msize[5] = x2s, x2e, x3s, x3e
    r.mind[2], r.mind[3] = n2, n3
    r.text = set_param(r.text, 'mesh', 'x2min', '%.16g' % x2s)
    r.text = set_param(r.text, 'mesh', 'x2max', '%.16g' % x2e)
    r.text = set_param(r.text, 'mesh', 'x3min', '%.16g' % x3s)
    r.text = set_param(r.text, 'mesh', 'x3max', '%.16g' % x3e)
    return o


# ------------------------------------------------------------------------------- writing
def zorder(nb2, nb3):
    """Z-ordered logical locations (lx1=0) of an nb2 x nb3 root grid, as the tree orders
    them (x1 bit lowest, then x2, then x3)"""
    def key(l2, l3):
        k = 0
        for b in range(16):
            k |= ((l2 >> b) & 1) << (3 * b + 1)
            k |= ((l3 >> b) & 1) << (3 * b + 2)
        return k
    ll = [(l2, l3) for l3 in range(nb3) for l2 in range(nb2)]
    return sorted(ll, key=lambda t: key(*t))


def set_param(text, block, par, val):
    pat = re.compile(r'(<%s>\n(?:(?!<).*\n)*?)(%s\s*=\s*)(\S+)' % (re.escape(block),
                                                                   re.escape(par)))
    new, n = pat.subn(lambda m: m.group(1) + m.group(2) + str(val), text, count=1)
    if n != 1:
        sys.exit('parameter %s/%s not found' % (block, par))
    return new


def write_rst(r, out, fn, f, mbx2, mbx3, keep, ghosts, keep_dt, src_mb=None, dtfac=0.5):
    ng = r.ng
    N2, N3 = r.mind[2] * f, r.mind[3] * f
    if N2 % mbx2 or N3 % mbx3:
        sys.exit('meshblock size does not divide the mesh')
    nb2, nb3 = N2 // mbx2, N3 // mbx3
    nmb = nb2 * nb3
    import math
    root = max(1, math.ceil(math.log2(max(nb2, nb3, 1)))) if max(nb2, nb3) > 1 else 0
    text = r.text
    if f != 1 or mbx2 != r.bind[2] or mbx3 != r.bind[3]:
        text = set_param(text, 'mesh', 'nx2', N2)
        text = set_param(text, 'mesh', 'nx3', N3)
        text = set_param(text, 'meshblock', 'nx2', mbx2)
        text = set_param(text, 'meshblock', 'nx3', mbx3)
    msize = list(r.msize)
    msize[7] = (msize[4] - msize[1]) / N2
    msize[8] = (msize[5] - msize[2]) / N3
    if f == 1:
        msize[7], msize[8] = r.msize[7], r.msize[8]
    mind = list(r.mind)
    mind[2], mind[3] = N2, N3
    mind[6], mind[7] = ng, ng + N2 - 1
    mind[8], mind[9] = ng, ng + N3 - 1
    bind = list(r.bind)
    bind[2], bind[3] = mbx2, mbx3
    bind[6], bind[7] = ng, ng + mbx2 - 1
    bind[8], bind[9] = ng, ng + mbx3 - 1
    if bind[10] or bind[11]:      # coarse indices are set on the MeshBlock: keep the rule
        bind[11], bind[12] = mbx2 // 2, mbx3 // 2
        bind[15], bind[16] = ng, ng + mbx2 // 2 - 1
        bind[17], bind[18] = ng, ng + mbx3 // 2 - 1
    o1, o2, o3 = r.o1, mbx2 + 2 * ng, mbx3 + 2 * ng
    keepset = {'hyd', 'm1', 'f1', 'f2', 'f3', 'ctr'}
    if keep:
        keepset |= {'wt', 'wd', 'pred', 'eint', 'rss'}
    shapes = {'hyd': (r.nhyd, o3, o2, o1), 'm1': (4, o3, o2, o1), 'f1': (o3, o2, o1 + 1),
              'f2': (o3, o2 + 1, o1), 'f3': (o3 + 1, o2, o1), 'wt': (1, o3, o2, o1),
              'wd': (2, o3, o2, o1), 'pred': (r.npred, o3, o2, o1),
              'eint': (r.neint, o3, o2, o1), 'ctr': (r.nctr, o3, o2, o1),
              'rss': (r.nrss, o3, o2, o1)}
    order = [s[0] for s in r.sec if s[0] in keepset]
    dsize = sum(int(np.prod(shapes[nm])) for nm in order) * 8
    blocks = []
    dropped = {'pred': 'pred', 'eint': 'eint', 'rss': 'rss'}
    for mg, pl in r.blocks:
        nm = MAGIC[mg]
        if nm in dropped and not keep:
            continue
        blocks.append((mg, pl))
    ll = zorder(nb2, nb3)
    with open(fn, 'wb') as fh:
        fh.write(text.encode())
        fh.write(struct.pack('<ii', nmb, root))
        fh.write(struct.pack('<9d', *msize))
        fh.write(struct.pack('<%di' % NIND, *mind))
        fh.write(struct.pack('<%di' % NIND, *bind))
        fh.write(struct.pack('<2d', r.time, r.dt if keep_dt else dtfac * r.dt))
        fh.write(struct.pack('<i', r.ncycle))
        for l2, l3 in ll:
            fh.write(struct.pack('<5i', 0, l2, l3, root, int(r.lloc[0, 4])))
        fh.write(np.full(nmb, r.cost[0], dtype='<f4').tobytes())
        for mg, pl in blocks:
            fh.write(mg)
            fh.write(struct.pack('<Q', len(pl)))
            fh.write(pl)
        fh.write(struct.pack('<Q', dsize))
        jj = np.arange(o2) - ng
        kk = np.arange(o3) - ng
        jf = np.arange(o2 + 1) - ng
        kf = np.arange(o3 + 1) - ng
        for m, (l2, l3) in enumerate(ll):
            J = (l2 * mbx2 + jj) % N2
            K = (l3 * mbx3 + kk) % N3
            JF = (l2 * mbx2 + jf) % N2
            KF = (l3 * mbx3 + kf) % N3
            if ghosts == 'source':
                sd = mb_data(r, src_mb[(l2, l3)])
            for nm in order:
                if nm == 'f1':
                    a = out[nm][np.ix_(K, J)]
                elif nm == 'f2':
                    a = out[nm][np.ix_(K, JF)]
                elif nm == 'f3':
                    a = out[nm][np.ix_(KF, J)]
                else:
                    a = out[nm][:, K][:, :, J]
                a = np.ascontiguousarray(a, dtype='<f8')
                if ghosts == 'source':
                    s = np.array(sd[nm])
                    act = np.zeros(s.shape, dtype=bool)
                    if nm in ('f1',):
                        act[ng:ng + mbx3, ng:ng + mbx2] = True
                    elif nm == 'f2':
                        act[ng:ng + mbx3, ng:ng + mbx2] = True
                    elif nm == 'f3':
                        act[ng:ng + mbx3, ng:ng + mbx2] = True
                    else:
                        act[:, ng:ng + mbx3, ng:ng + mbx2] = True
                    a = np.where(act, a, s)
                assert a.shape == shapes[nm], (nm, a.shape, shapes[nm])
                fh.write(a.tobytes())
    return nmb, dsize


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('inp')
    ap.add_argument('out')
    ap.add_argument('--factor', type=int, default=2)
    ap.add_argument('--mbx2', type=int, default=None)
    ap.add_argument('--mbx3', type=int, default=None)
    ap.add_argument('--keep-caches', action='store_true')
    ap.add_argument('--ghosts', default='wrap', choices=['wrap', 'source'])
    ap.add_argument('--keep-dt', action='store_true')
    ap.add_argument(
        '--dt-factor',
        type=float,
        default=0.5,
        help='header dt = this x the source dt (the first cycle takes it as is)')
    ap.add_argument('--no-slopes', action='store_true')
    ap.add_argument('--force-math', action='store_true')
    ap.add_argument('--stats', default=None, help='npz with per-cell diagnostics')
    ap.add_argument(
        '--grid',
        default=None,
        help='.npy: 18 StretchRPoly parameters of the NEW '
        'radial grid (same x1min/x1max/ng); with --nx1')
    ap.add_argument('--nx1', type=int, default=None)
    ap.add_argument('--sector', type=float, nargs=4, default=None,
                    metavar=('X2MIN', 'X2MAX', 'X3MIN', 'X3MAX'),
                    help='cut this angular sector (source cell edges) before the '
                    'angular prolongation; the new mesh is periodic on it as before')
    a = ap.parse_args()
    r = read_rst(a.inp)
    m = re.search(r'<rad_m1>\n(?:(?!<).*\n)*?c_light\s*=\s*(\S+)', r.text)
    clight = float(m.group(1))
    print('read %s: t=%.6f cycle %d dt %.6g, %d MeshBlocks %dx%dx%d, nhydro %d, pred %d, '
          'eint %d, ctr %d, rss %d, blocks %s' % (
              a.inp, r.time, r.ncycle, r.dt, r.nmb, r.bind[1], r.bind[2], r.bind[3],
              r.nhyd, r.npred, r.neint, r.nctr, r.nrss, [MAGIC[b[0]] for b in r.blocks]))
    G = gather(r)
    if a.grid:
        G, sr = radial_remap(r, G, np.load(a.grid), a.nx1, use_slopes=not a.no_slopes)
        print('radial: %d -> %d cells; ' % (len(r.fo) - 1 - 2 * r.ng, a.nx1)
              + ', '.join('%s %s' % (k, v) for k, v in sr.items()))
    if a.sector:
        G = extract_sector(r, G, *a.sector)
    out, st = remap(r, G, a.factor, clight, use_slopes=not a.no_slopes,
                    force_math=a.force_math)
    if a.sector:
        out = strip_sector(r, out, a.factor)
    print(
        'fallbacks (coarse cells -> injection): hydro %d, m1 %d of %d; max parent |F|/cE '
        '%.4g' %
        (st['hyd_fallback'],
         st['m1_fallback'],
            G['hyd'][0].size,
            st['m1_parent_fc']))
    if a.grid:
        rf = r.fn[r.ng:r.ng + a.nx1 + 1]
        xv = 0.5 * (rf[1:] + rf[:-1]) / (r.msize[0] / 0.4)
        a_ = slice(r.ng, r.ng + a.nx1)
        ncol = G['hyd'][0].shape[0] * G['hyd'][0].shape[1]
        for lo, hi in ((0.40, 0.55), (0.55, 0.70),
                       (0.70, 0.85), (0.85, 1.05), (1.05, 3.0)):
            m = (xv >= lo) & (xv < hi)
            print('   r/R %.2f-%.2f: angular fallback fraction hydro %.4f, m1 %.4f' % (
                lo, hi, st['hyd_fb_i'][a_][m].sum() / (ncol * m.sum()),
                st['m1_fb_i'][a_][m].sum() / (ncol * m.sum())))
    mbx2 = a.mbx2 or r.bind[2]
    mbx3 = a.mbx3 or r.bind[3]
    src_mb = None
    if a.ghosts == 'source':
        if a.factor != 1 or mbx2 != r.bind[2] or mbx3 != r.bind[3]:
            sys.exit('--ghosts source needs the identity remap')
        src_mb = {(int(r.lloc[m, 1]), int(r.lloc[m, 2])): m for m in range(r.nmb)}
    nmb, ds = write_rst(r, out, a.out, a.factor, mbx2, mbx3, a.keep_caches, a.ghosts,
                        a.keep_dt, src_mb, a.dt_factor)
    print('wrote %s: %d MeshBlocks %dx%dx%d, record %d bytes, caches %s, dt %s' % (
        a.out, nmb, r.bind[1], mbx2, mbx3, ds, 'kept' if a.keep_caches else 'dropped',
        'kept' if a.keep_dt else '%.3g x source' % a.dt_factor))
    if a.stats:
        np.savez(a.stats, phi=st.get('phi_check'))


if __name__ == '__main__':
    main()
