#!/usr/bin/env python3
"""Radial remap of a cubed-sphere deep_hot_jupiter_rt HYDRO restart onto a new radial
grid.

usage:
  dhj_remap.py RST OUTDIR [--grid "mesh/nx1=256 meshblock/nx1=256 mesh/x1min=... c1..c8"]
               [--method avg|point] [--interp-identity] [--check-only]

Writes OUTDIR/remap.dat (the file read by problem/remap_file, src/pgen/
deep_hot_jupiter_rt.cpp) and OUTDIR/remap.athinput: the input embedded in RST with the
grid keys of --grid, time/start_time = t(RST) and problem/remap_file set, and nothing else
changed (output counters and last_time carry on).  Run it as a FRESH start:
  athena -i OUTDIR/remap.athinput -d RUNDIR
--grid takes the planet_setup.py override string (G_PROD in grid_w121*.env); omitted =
the restart's own grid (identity remap: a straight copy of the active cells, bitwise).
x1min/x1max must be those of the restart.  HORIZONTAL (dhj-remap-h): --grid may also
raise mesh/nx2 = mesh/nx3 by an integer factor f (C32 -> C256: f = 8) and set any
meshblock/nx2, nx3 that divides them (one MeshBlock in r stays required).  Each old cell
is prolonged conservatively: MC-limited piecewise linear in (xi, eta) (the equiangular
chart), slopes from the face neighbours incl. the restart's ghosts (at a panel seam the
code's resampled cross-panel ghosts, vectors in this panel's basis), child offsets
centred on the exact-solid-angle mean so sum_children q dOmega = q dOmega_parent to
round-off for rho, e and the three stored (covariant) momentum components (these via
u = m/rho, MC-limited, plus one shift per parent).  E_tot = e + KE + rho Phi is NOT
conserved by default: with momentum conserved, any sub-cell velocity profile raises KE
(Cauchy-Schwarz; +3e-3 of KE at rot 300), and rho Phi at the child centres differs from
the parent's by the quadrature of the centrifugal term (-7e-4 of sum rho Phi).
--henergy total takes both from e (e x f per parent, f capped at 1 +- FTOL): at rot 300
that cools shear zones at 1e-4..1e-3 bar by 10-18 % (KE term), so it is not the default.
dhj-remap-3x: e and p are prolonged through e/rho, p/rho plus one shift per parent
(--hthermal specific; T stays bounded by the neighbours), and L_z is restored to
round-off per old radial shell by one uniform rotation increment (--hlz shell; C32 ->
C128 at rot 300: max 1.1 m/s, |dOmega| < 1e-8 rad/s).
RAISED TOP (dhj-remap-3x): x1max may grow (x1min fixed).  Above the old top-cell
centroid every column is filled hydrostatically and isothermally with the code's
potential (rot_potential incl.) at T = max(T_top, --tfill-min 1000 K), p continuous,
e/p and u_h of the top cell, u_r = 0 (--eos-table gives T_top, linearly extrapolated
below eos_logd_min as the code does).  Mass, e and the column momenta (m1; r m2, r m3,
i.e. L_z) are rescaled over the old domain only (the straddling cell split by volume);
the fill above is reported separately.  All 18 StretchRPoly terms (8 coefficients,
bumps b1/b2, plateau) are read; unknown f_stretch_r_* keys are fatal.
The
radial remap below then runs on every fine column.  Records are written in the code's
Z order of the new layout (panel, Morton at the root level); the pgen checks every
LogicalLocation.  Streaming: one new MeshBlock at a time.  --selftest = gate (c).

Method (per column; every column shares the radial grid), default 'avg':
  the restart holds cell AVERAGES.  ln rho and ln e (e = internal energy density, the
  EINTRST1 slab: the code's own w0(IEN)) are converted to centroid values
  (- dr^2/24 (L'' + L'^2)), interpolated at the new centroids with a monotone cubic
  (PCHIP, through the old centroids, ghosts included so nothing is extrapolated) and
  converted back to averages of the new cells.  u = m/rho (the stored components,
  covariant on the cubed sphere; the radial remap does not change the basis) is
  interpolated the same way, linearly.  Then per column: rho x f_m (column mass exact),
  e x f_e (column internal energy exact), u + du (column momentum exact); a common
  factor keeps T and the relative HSE residual.  KE then changes by the sub-cell velocity
  profile (~1e-3 of KE, ~1e-7 of e + KE) and the discrete gravitational energy by the
  finer placement of the mass (a few % of sum rho Phi, ~1e-3 of E_tot): neither is put
  back into e, which would heat the gas spuriously.
  'point' = the radab_0924 prototype (no average <-> centroid conversion).
  A fully conservative remap of the cumulative mass / e (PCHIP of the cumulative) was
  tried and dropped: it doubles the HSE residual in the deep and is ill-posed in the
  first and last old cells (remap_0929/RESULTS.md).
"""
import argparse
import os
import re
import struct
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

NSTR = 8
NBUMP = 2
AVG_ITMAX = 200   # fixed-point iterations of the average -> centroid inversion
FTOL = 0.1     # horizontal --henergy total: |f - 1| cap (the thin top, see main)
BAR = 1.0e6


# ------------------------------------------------------------------ input text helpers
def parse_input(txt):
    """-> {block: {key: value-string}} from athinput text (comments stripped)"""
    d, blk = {}, None
    for ln in txt.splitlines():
        s = ln.split('#', 1)[0].strip()
        if not s:
            continue
        m = re.match(r'<(.+)>', s)
        if m:
            blk = m.group(1).strip()
            d.setdefault(blk, {})
            continue
        if '=' in s and blk is not None:
            k, v = s.split('=', 1)
            d[blk][k.strip()] = v.strip()
    return d


def set_key(txt, blk, key, val, note='remap'):
    """set <blk>/key = val in athinput text (replace the line or add it after <blk>)"""
    lines = txt.splitlines()
    out, inblk, done = [], False, False
    for ln in lines:
        m = re.match(r'\s*<(.+)>', ln)
        if m:
            if inblk and not done:
                out.append('%s = %s   # %s' % (key, val, note))
                done = True
            inblk = (m.group(1).strip() == blk)
        elif inblk and not done:
            s = ln.split('#', 1)[0]
            if '=' in s and s.split('=', 1)[0].strip() == key:
                out.append('%-20s = %s   # %s' % (key, val, note))
                done = True
                continue
        out.append(ln)
    if not done:
        if inblk:
            out.append('%s = %s   # %s' % (key, val, note))
        else:
            out += ['<%s>' % blk, '%s = %s   # %s' % (key, val, note)]
    return '\n'.join(out) + '\n'


def grid_of(par):
    """(nx1, x1min, x1max, c): c = the 18 StretchRPoly terms of grid_stretch.hpp, read as
    Mesh does: f_stretch_r_c1..c8, then bumps b1, b2 (amp, x, w; only when _amp is
    given), then the plateau (p_amp, p_xa, p_xb, p_w; only when p_amp is given).  A
    term left out builds a different grid silently, so every term is read here."""
    m = par['mesh']
    c = [float(m.get('f_stretch_r_c%d' % k, 0.0)) for k in range(1, NSTR+1)]
    for b in range(1, NBUMP+1):
        kb = 'f_stretch_r_b%d' % b
        if kb + '_amp' in m:
            c += [float(m[kb + '_amp']), float(m[kb + '_x']), float(m[kb + '_w'])]
        else:
            c += [0.0, 0.0, 0.0]
    if 'f_stretch_r_p_amp' in m:
        c += [float(m['f_stretch_r_p_' + k]) for k in ('amp', 'xa', 'xb', 'w')]
    else:
        c += [0.0, 0.0, 0.0, 0.0]
    known = (['f_stretch_r_c%d' % k for k in range(1, NSTR+1)] +
             ['f_stretch_r_b%d_%s' % (b, s) for b in range(1, NBUMP+1)
              for s in ('amp', 'x', 'w')] +
             ['f_stretch_r_p_' + k for k in ('amp', 'xa', 'xb', 'w')])
    unk = [k for k in m if k.startswith('f_stretch_r_') and k not in known]
    if unk:
        raise SystemExit('unknown radial-stretch keys (grid would be wrong): %s' % unk)
    poly = m.get('use_grid_stretch_r_poly', 'false').lower() in ('true', '1')
    if not poly or m.get('use_grid_stretch_r', 'false').lower() in ('true', '1'):
        raise SystemExit('only the polynomial radial stretch is supported')
    return (int(m['nx1']), float(m['x1min']), float(m['x1max']), np.array(c))


def _logcosh(x):
    ax = np.abs(x)
    return ax + np.log1p(np.exp(-2.0*ax)) - 0.69314718055994530942


def _plat_g(xa, xb, w, xi):
    return 0.5*w*(_logcosh((xi - xa)/w) - _logcosh(-xa/w)
                  - _logcosh((xi - xb)/w) + _logcosh(-xb/w))


def edges(nx1, r0, r1, c, ng=0):
    """cell edges as Coordinates builds them (LeftEdgeX + StretchRPoly incl. the bumps
    and the plateau, same operation order), ghosts incl."""
    i = np.arange(-ng, nx1+ng+1)
    t = i/float(nx1)                           # LeftEdgeX (cell_locations.hpp), exactly
    x = (t*r1 - t*r0) - (0.5*r1 - 0.5*r0) + (0.5*r0 + 0.5*r1)
    xi = (x - r0)/(r1 - r0)
    u = xi.copy()
    xik = xi.copy()
    for k in range(NSTR):
        u += c[k]*xik*(1.0 - xi)
        xik = xik*xi
    for b in range(NBUMP):
        a, xb, w = c[NSTR + 3*b:NSTR + 3*b + 3]
        if a == 0.0:
            continue
        u += a*w*(np.tanh((xi - xb)/w) - (1.0 - xi)*np.tanh(-xb/w)
                  - xi*np.tanh((1.0 - xb)/w))
    pa, pxa, pxb, pw = c[NSTR + 3*NBUMP:NSTR + 3*NBUMP + 4]
    if pa != 0.0:
        u += pa*(_plat_g(pxa, pxb, pw, xi) - xi*_plat_g(pxa, pxb, pw, 1.0))
    return r0 + (r1 - r0)*u


def centroids(e):
    q = e[:-1]/e[1:]
    return 0.25*(q*q + 1.0)/((q*q + q + 1.0)/3.0)*(e[1:] + e[:-1])


# ------------------------------------------------------------------ restart reader
MAGICS = (b'PGENST01', b'RTWARM01', b'EINTRST1', b'CKRST001', b'M1PRED01', b'M1TIME2A',
          b'M1ONEP01', b'M1MRWIN1', b'M1CTRLX1')


def read_rst(fn):
    with open(fn, 'rb') as f:
        b = f.read()
    i = b.find(b'<par_end>')
    i = b.find(b'\n', i) + 1
    txt = b[:i].decode(errors='replace')
    nmb, root = struct.unpack('<ii', b[i:i+8])
    o = i + 8 + 9*8
    mbi = struct.unpack('<19i', b[o+76:o+152])
    o += 152
    t, dt = struct.unpack('<dd', b[o:o+16])
    ncyc = struct.unpack('<i', b[o+16:o+20])[0]
    o += 20
    lloc = np.frombuffer(b[o:o+20*nmb], '<i4').reshape(nmb, 5).copy()
    par = parse_input(txt)
    ng, nx1, nx2, nx3 = mbi[0], mbi[1], mbi[2], mbi[3]
    O1, O2, O3 = nx1+2*ng, nx2+2*ng, nx3+2*ng
    NC = O1*O2*O3
    # which optional sections exist (they sit in the step-3 header, before the data)
    hdr_end = o + 24*nmb + 4096
    sec = [mg.decode() for mg in MAGICS if b.find(mg, o, hdr_end) > 0]
    for bad in ('RTWARM01', 'M1PRED01', 'M1TIME2A', 'M1CTRLX1'):
        if bad in sec:
            raise SystemExit('restart carries section %s: slab order not handled' % bad)
    if 'mhd' in par:
        raise SystemExit('MHD restart: only hydro restarts are supported')
    for nv in range(8, 64):
        ds = 8*NC*nv
        st = len(b) - nmb*ds
        if st > i and int(np.frombuffer(b[st-8:st], '<u8')[0]) == ds:
            break
    else:
        raise SystemExit('restart layout not recognised: ' + fn)
    a = np.frombuffer(b[st:], '<f8').reshape(nmb, nv, O3, O2, O1)
    h = dict(text=txt, par=par, nmb=nmb, t=t, dt=dt, ncycle=ncyc, lloc=lloc, ng=ng,
             nx1=nx1, nx2=nx2, nx3=nx3, nv=nv, sections=sec)
    h['u'] = a[:, 0:5]
    h['p'] = a[:, 6]
    h['eint'] = a[:, 8] if 'EINTRST1' in sec else None
    return h


# ------------------------------------------------------------------ geometry
def cs_columns(h):
    """cos_cell (tangent-basis angle) and solid angle of every column (nmb, nx3, nx2)"""
    nmb, nx2, nx3 = h['nmb'], h['nx2'], h['nx3']
    mesh = h['par']['mesh']
    n2, n3 = int(mesh['nx2']), int(mesh['nx3'])     # per panel
    cosc = np.empty((nmb, nx3, nx2))
    dom = np.empty_like(cosc)
    for m in range(nmb):
        lx2, lx3 = h['lloc'][m, 1], h['lloc'][m, 2]
        x2 = -1.0 + 2.0*(lx2*nx2 + np.arange(nx2) + 0.5)/n2
        x3 = -1.0 + 2.0*(lx3*nx3 + np.arange(nx3) + 0.5)/n3
        X, Y = np.meshgrid(np.tan(np.pi/4*x2), np.tan(np.pi/4*x3))   # (nx3, nx2)
        C, D = np.sqrt(1 + X*X), np.sqrt(1 + Y*Y)
        cosc[m] = -X*Y/(C*D)
        dxi = np.pi/4*2.0/n2
        deta = np.pi/4*2.0/n3
        dom[m] = dxi*deta*(1 + X*X)*(1 + Y*Y)/(1 + X*X + Y*Y)**1.5
    return cosc, dom


def kinetic(rho, m1, m2, m3, c):
    """0.5 m_i v^i with the gnomonic metric on the angular pair (c broadcast on r)"""
    det = 1.0 - c*c
    v1 = m1/rho
    v2 = (m2 - c*m3)/(rho*det)
    v3 = (m3 - c*m2)/(rho*det)
    return 0.5*(m1*v1 + m2*v2 + m3*v3)


# ------------------------------------------------------------------ remap kernels
def avg_corr(L, x, dx, cap=0.02):
    """ln(cell average) - ln(value at the centroid) of a profile exp(L(x)) over cells of
    width dx: dx^2/24 (L'' + L'^2) (second order; the spherical r^2 weight adds a term
    ~ H/r smaller).  Capped at +-cap: the floor-noise top has factor-10 jumps per cell."""
    d1 = np.gradient(L, x, axis=-1, edge_order=2)
    d2 = np.gradient(d1, x, axis=-1, edge_order=2)
    return np.clip(dx**2/24.0*(d2 + d1*d1), -cap, cap)


# -------------------------------------------------------- horizontal (dhj-remap-h)
# z components of PanelFrame (a, b, n) in src/coordinates/cubed_sphere.hpp, per panel
PANEL_Z = ((0.0, 1.0, 0.0), (0.0, 1.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
           (0.0, 1.0, 0.0), (0.0, 0.0, -1.0))
# the full PanelFrame (a, b, n), for L_z
PANEL_FRAME = (((0, 1, 0), (0, 0, 1), (1, 0, 0)), ((-1, 0, 0), (0, 0, 1), (0, 1, 0)),
               ((0, -1, 0), (0, 0, 1), (-1, 0, 0)), ((0, 1, 0), (-1, 0, 0), (0, 0, 1)),
               ((1, 0, 0), (0, 0, 1), (0, -1, 0)), ((0, 1, 0), (1, 0, 0), (0, 0, -1)))


def corner_omega(X, Y):
    return np.arctan(X*Y/np.sqrt(1.0 + X*X + Y*Y))


def panel_cols(k0, k1, j0, j1, n2, n3, panel):
    """exact solid angle (GnomonicSolidAngle), cos_cell and sin^2(colatitude) of the panel
    columns [k0,k1) x [j0,j1) of a panel with n2 x n3 cells; arrays (k1-k0, j1-j0)"""
    xe = np.tan(0.25*np.pi*(-1.0 + 2.0*np.arange(j0, j1 + 1)/n2))
    ye = np.tan(0.25*np.pi*(-1.0 + 2.0*np.arange(k0, k1 + 1)/n3))
    om = (corner_omega(xe[None, 1:], ye[1:, None])
          - corner_omega(xe[None, :-1], ye[1:, None])
          - corner_omega(xe[None, 1:], ye[:-1, None])
          + corner_omega(xe[None, :-1], ye[:-1, None]))
    X = np.tan(0.25*np.pi*(-1.0 + 2.0*(np.arange(j0, j1) + 0.5)/n2))[None, :]
    Y = np.tan(0.25*np.pi*(-1.0 + 2.0*(np.arange(k0, k1) + 0.5)/n3))[:, None]
    C, D = np.sqrt(1 + X*X), np.sqrt(1 + Y*Y)
    cosc = -X*Y/(C*D)
    a2, b2, n2z = PANEL_Z[panel]
    cz = (a2*X + b2*Y + n2z)/np.sqrt(1.0 + X*X + Y*Y)
    sin2 = 1.0 - cz*cz
    return om, cosc + 0*om, sin2 + 0*om


def zrot_cov(k0, k1, j0, j1, n2, n3, panel):
    """covariant components (w1, w2) = ((z x rhat).e1, (z x rhat).e2) at the cell centres
    of panel columns [k0,k1) x [j0,j1): solid rotation about z per unit Omega and r"""
    a, b, n = (np.array(v, float) for v in PANEL_FRAME[panel])
    x = np.tan(0.25*np.pi*(-1.0 + 2.0*(np.arange(j0, j1) + 0.5)/n2))[None, :, None]
    y = np.tan(0.25*np.pi*(-1.0 + 2.0*(np.arange(k0, k1) + 0.5)/n3))[:, None, None]
    dl = np.sqrt(1.0 + x*x + y*y)
    rh = (a*x + b*y + n)/dl
    e1 = (a*dl - rh*x)/np.sqrt(1.0 + y*y)
    e2 = (b*dl - rh*y)/np.sqrt(1.0 + x*x)
    zr = np.stack((-rh[..., 1], rh[..., 0], 0.0*rh[..., 0]), axis=-1)
    return (zr*e1).sum(-1), (zr*e2).sum(-1)


def lz_coef(k0, k1, j0, j1, n2, n3, panel):
    """L_z per unit volume = r (A2 m2 + A3 m3) for the stored (covariant, unit tangent
    basis of PanelTangents) momenta m2, m3 at the cell centres of panel columns
    [k0,k1) x [j0,j1): A2 = (w1 - c w2)/det, A3 = (w2 - c w1)/det, w = (rhat x e)_z,
    c = cos_cell, det = 1 - c^2; arrays (k1-k0, j1-j0)"""
    a, b, n = (np.array(v, float) for v in PANEL_FRAME[panel])
    x = np.tan(0.25*np.pi*(-1.0 + 2.0*(np.arange(j0, j1) + 0.5)/n2))[None, :, None]
    y = np.tan(0.25*np.pi*(-1.0 + 2.0*(np.arange(k0, k1) + 0.5)/n3))[:, None, None]
    dl = np.sqrt(1.0 + x*x + y*y)
    rh = (a*x + b*y + n)/dl
    e1 = (a*dl - rh*x)/np.sqrt(1.0 + y*y)
    e2 = (b*dl - rh*y)/np.sqrt(1.0 + x*x)
    w1 = rh[..., 0]*e1[..., 1] - rh[..., 1]*e1[..., 0]
    w2 = rh[..., 0]*e2[..., 1] - rh[..., 1]*e2[..., 0]
    c = (e1*e2).sum(-1)
    det = 1.0 - c*c
    return (w1 - c*w2)/det, (w2 - c*w1)/det


def eos_temperature(eos, rho, eint):
    """T [K] of (rho, eint) from the dhjcs.EOS table dump, bilinear in (log rho, log T)
    like dhjcs.EOS.invert, but below eos_logd_min LINEARLY EXTRAPOLATED in log rho from
    the first table interval (the code's off-table rule, eos_table.hpp Interpolate)
    instead of clamped to the edge column"""
    rho = np.asarray(rho, float)
    sh = rho.shape
    x = (np.log10(rho.ravel()) - eos.xmin)/eos.dx
    ix = np.clip(np.floor(x).astype(int), 0, eos.nx - 2)
    fx = np.where(x < 0.0, x, np.clip(x - ix, 0.0, 1.0))
    le = eos.le[:, ix]*(1 - fx) + eos.le[:, ix + 1]*fx            # (ny, N)
    tgt = np.log10(np.asarray(eint, float).ravel()/rho.ravel())
    fin = np.isfinite(le)
    lef = np.where(fin, le, -np.inf)
    first = np.argmax(fin, axis=0)
    k = np.clip((lef < tgt[None, :]).sum(axis=0), first + 1, eos.ny - 1)
    n = np.arange(k.size)
    e0, e1 = le[k - 1, n], le[k, n]
    f = np.clip((tgt - e0)/(e1 - e0), 0.0, 1.0)
    return (10.0**(eos.lt[k - 1] + f*eos.dy)).reshape(sh)


def panel_arrays(h, q, n2, n3):
    """per-panel array (6, n3+2, n2+2, O1) of a restart field q (nmb, O3, O2, O1): the
    active cells plus a one-cell face-ghost ring from the restart's own ghosts (at a panel
    seam these are the code's resampled cross-panel ghosts, vectors already in THIS
    panel's basis). Diagonal corners stay NaN: the dimension-by-dimension slopes never
    read them."""
    ng, nx2, nx3 = h['ng'], h['nx2'], h['nx3']
    P = np.full((6, n3 + 2, n2 + 2, q.shape[-1]), np.nan)
    a2, a3 = slice(ng, ng + nx2), slice(ng, ng + nx3)
    for m in range(h['nmb']):
        _, lx2, lx3, _, pn = h['lloc'][m]
        j0, k0 = lx2*nx2, lx3*nx3
        sj, sk = slice(1 + j0, 1 + j0 + nx2), slice(1 + k0, 1 + k0 + nx3)
        P[pn, sk, sj] = q[m, a3, a2]
        if j0 == 0:
            P[pn, sk, 0] = q[m, a3, ng - 1]
        if j0 + nx2 == n2:
            P[pn, sk, n2 + 1] = q[m, a3, ng + nx2]
        if k0 == 0:
            P[pn, 0, sj] = q[m, ng - 1, a2]
        if k0 + nx3 == n3:
            P[pn, n3 + 1, sj] = q[m, ng + nx3, a2]
    return P


def mc(dm, dp):
    """monotonized-central limited slope (per parent cell width)"""
    return np.where(dm*dp > 0.0,
                    np.sign(dm)*np.minimum(np.minimum(2*np.abs(dm), 2*np.abs(dp)),
                                           0.5*np.abs(dm + dp)), 0.0)


def slopes(P, pn, k0, k1, j0, j1):
    """cell value and MC slopes in j (x2) and k (x3) of panel parents [k0,k1) x [j0,j1)"""
    c = P[pn, 1 + k0:1 + k1, 1 + j0:1 + j1]
    sj = mc(c - P[pn, 1 + k0:1 + k1, j0:j1], P[pn, 1 + k0:1 + k1, 2 + j0:2 + j1] - c)
    sk = mc(c - P[pn, k0:k1, 1 + j0:1 + j1], P[pn, 2 + k0:2 + k1, 1 + j0:1 + j1] - c)
    return c, sj, sk


class HProlong:
    """conservative limited piecewise-linear prolongation of one panel patch by an
    integer factor f per direction.  In each parent the linear profile
    q + s_j (d_a - d_j0) + s_k (d_b - d_k0) is evaluated at the child offsets d_a =
    (a + 1/2)/f - 1/2 (parent units, the equiangular chart is uniform in xi, eta), with
    d_j0, d_k0 the EXACT-solid-angle weighted means of the child offsets, so sum_children
    q_c dOmega_c = q dOmega_parent to round-off (the volume is dOmega x the radial shell,
    so this holds per radial cell)."""

    def __init__(self, f, pk0, pk1, pj0, pj1, n2, n3, pn):
        self.f = f
        self.shape = (pk1 - pk0, pj1 - pj0)
        w, cosc, sin2 = panel_cols(pk0*f, pk1*f, pj0*f, pj1*f, n2*f, n3*f, pn)
        kp, jp = self.shape
        w4 = w.reshape(kp, f, jp, f)
        d = (np.arange(f) + 0.5)/f - 0.5
        W = w4.sum(axis=(1, 3))
        dj0 = (w4*d[None, None, None, :]).sum(axis=(1, 3))/W
        dk0 = (w4*d[None, :, None, None]).sum(axis=(1, 3))/W
        self.oj = d[None, None, None, :] - dj0[:, None, :, None]     # (kp, f, jp, f)
        self.ok = d[None, :, None, None] - dk0[:, None, :, None]
        self.w4 = w4
        self.cosc, self.sin2 = cosc, sin2          # of the children (kp f, jp f)

    def children(self, c, sj, sk, positive=False, floor=0.1, qmin=0.0):
        """children (kp f, jp f, O1) of parents c (kp, jp, O1) with slopes sj, sk.
        positive: both slopes scaled down in any parent whose smallest child would drop
        below floor x parent (conservation kept; returns the number of parents scaled)."""
        f = self.f
        kp, jp = self.shape
        oj = self.oj[..., None]
        ok = self.ok[..., None]
        sj4 = sj[:, None, :, None, :]
        sk4 = sk[:, None, :, None, :]
        nsc = 0
        if positive:
            dmin = (sj4*oj + sk4*ok).min(axis=(1, 3))       # (kp, jp, O1)
            lim = np.minimum(np.maximum((floor - 1.0)*c, qmin - c), 0.0)
            bad = dmin < lim - 1e-12*np.abs(c)
            nsc = int(bad.sum())
            if nsc:
                th = np.where(bad, lim/np.where(bad, dmin, -1.0), 1.0)
                sj4 = sj4*th[:, None, :, None, :]
                sk4 = sk4*th[:, None, :, None, :]
        q = c[:, None, :, None, :] + sj4*oj + sk4*ok
        return q.reshape(kp*f, jp*f, -1), nsc

    def parent_sum(self, q):
        """sum_children q dOmega per parent: (kp f, jp f, O1) -> (kp, jp, O1)"""
        kp, jp = self.shape
        f = self.f
        return (q.reshape(kp, f, jp, f, -1)*self.w4[..., None]).sum(axis=(1, 3))


def zorder(lloc_list, level):
    """sort (lx1, lx2, lx3, level, panel) tuples as CreateZOrderedLLList numbers them:
    panel by panel, then the Morton key at the root level (x1 lowest bit)"""
    def key(ll):
        c = 0
        for b in range(level - 1, -1, -1):
            c = (c*8 + (((ll[2] >> b) & 1) << 2) + (((ll[1] >> b) & 1) << 1)
                 + ((ll[0] >> b) & 1))
        return (ll[4], c)
    return sorted(lloc_list, key=key)


def root_level(nb):
    lev = 0
    while (1 << lev) < max(nb):
        lev += 1
    return lev


def selftest():
    """gate (c): a smooth analytic field on one panel, exact cell averages (Gauss
    quadrature in dOmega) on coarse and fine grids; the prolonged coarse averages against
    the exact fine ones, L1 and Linf, for n = 8, 16, 32 per panel and f = 2, 8"""
    xg, wg = np.polynomial.legendre.leggauss(6)

    def fun(X, Y):
        return (2.0 + np.sin(1.3*np.arctan(X) + 0.4)*np.cos(0.9*np.arctan(Y) - 0.2)
                + 0.3*X*Y)

    def avgs(n, lo, hi):
        """exact averages of cells lo..hi-1 (x and y) incl. a ring; returns (hi-lo)^2"""
        e = 0.25*np.pi*(-1.0 + 2.0*np.arange(lo, hi + 1)/n)
        A = np.empty((hi - lo, hi - lo))
        xi = 0.5*(e[1:] + e[:-1])[:, None] + 0.5*np.diff(e)[:, None]*xg[None, :]
        wx = 0.5*np.diff(e)[:, None]*wg[None, :]
        X = np.tan(xi)
        for kk in range(hi - lo):
            Y = X[kk][:, None, None]
            Xx = X[None, :, :]
            dens = (1 + Xx*Xx)*(1 + Y*Y)/(1 + Xx*Xx + Y*Y)**1.5
            ww = wx[kk][:, None, None]*wx[None, :, :]*dens
            A[kk] = (fun(Xx, Y)*ww).sum(axis=(0, 2))/ww.sum(axis=(0, 2))
        return A
    out = []
    for f in (2, 8):
        prev = None
        for n in (8, 16, 32):
            Pc = avgs(n, -1, n + 1)[None, :, :, None]           # ring = extended chart
            hp = HProlong(f, 0, n, 0, n, n, n, 0)
            c, sj, sk = slopes(Pc, 0, 0, n, 0, n)
            q, _ = hp.children(c, sj, sk)
            ex = avgs(n*f, 0, n*f)
            err = np.abs(q[..., 0] - ex)
            cons = np.abs(hp.parent_sum(q)[..., 0]/hp.w4.sum(axis=(1, 3))
                          - c[..., 0]).max()
            l1 = err.mean()
            out.append('  f %d n %3d: L1 %.3e  Linf %.3e  ratio L1 %s  cons %.1e'
                       % (f, n, l1, err.max(),
                          '%.2f' % (prev/l1) if prev else '  - ', cons))
            prev = l1
    print('gate (c) smooth field, prolonged vs exact fine averages (order 2 = ratio 4):')
    print('\n'.join(out))


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument('rst', nargs='?')
    ap_.add_argument('outdir', nargs='?')
    ap_.add_argument('--grid', default='')
    ap_.add_argument('--method', default='avg', choices=('avg', 'point'))
    ap_.add_argument('--interp-identity', action='store_true',
                     help='run the interpolation even when the grid is unchanged')
    ap_.add_argument('--henergy', default='eint', choices=('total', 'eint'),
                     help='horizontal prolongation: conserve e per parent (default; '
                     'E_tot '
                     'then gains the sub-cell KE and the rho Phi quadrature change) or '
                     'the total energy e + KE + rho Phi (e rescaled, |f - 1| <= FTOL)')
    ap_.add_argument('--hthermal', default='specific', choices=('specific', 'density'),
                     help='horizontal prolongation of e and p: through e/rho, p/rho plus '
                     'a per-parent shift (default; T bounded by the neighbours) or as '
                     'densities (the 9896d452 scheme)')
    ap_.add_argument('--hlz', default='shell', choices=('shell', 'parent', 'none'),
                     help='horizontal: restore L_z per radial shell (default; one '
                     'uniform rotation increment per old radial cell), per parent, '
                     'or not')
    ap_.add_argument('--hse-sample', type=int, default=16,
                     help='HSE statistics from every n-th new MeshBlock')
    ap_.add_argument('--eos-table', default='',
                     help='general-EOS table dump (log10 e/rho, log10 p/rho on log rho, '
                     'log T) for T of the old top cell; required when x1max is raised')
    ap_.add_argument('--tfill-min', type=float, default=1000.0,
                     help='raised top: isothermal hydrostatic fill at T = max(T_top, '
                     'this) [K] per column (p continuous at the old top centroid)')
    ap_.add_argument('--check-only', action='store_true')
    ap_.add_argument('--selftest', action='store_true')
    a = ap_.parse_args()
    if a.selftest:
        selftest()
        return

    h = read_rst(a.rst)
    par = h['par']
    ng, nmb = h['ng'], h['nmb']
    nx1o, r0, r1, co = grid_of(par)
    if nx1o != h['nx1'] or int(par['meshblock']['nx1']) != nx1o:
        raise SystemExit('one MeshBlock in r is required')
    if 'cubed_sphere' in par['mesh'] and \
            par['mesh']['cubed_sphere'].lower() not in ('true', '1'):
        raise SystemExit('cubed-sphere restarts only')
    # target grid
    ov = dict(kv.split('=', 1) for kv in a.grid.split())
    parn = {b: dict(v) for b, v in par.items()}
    for k, v in ov.items():
        blk, key = k.split('/')
        parn.setdefault(blk, {})[key] = v
    nx1n, r0n, r1n, cn = grid_of(parn)
    if abs(r0n/r0 - 1) > 1e-12 or r1n < r1*(1.0 - 1e-12):
        raise SystemExit('x1min must not change and x1max must not drop (%g %g -> %g %g)'
                         % (r0, r1, r0n, r1n))
    # RAISED TOP (x1max_new > x1max_old): above the old top-cell centroid every column is
    # filled hydrostatically and isothermally (g_eff of the code's potential,
    # rot_potential incl.) at T = max(T_top, --tfill-min), p continuous at the old top
    # centroid, e/p and u_h of the top cell, u_r = 0.  All budget factors act on the old
    # domain only.
    ext = r1n > r1*(1.0 + 1e-12)
    if ext and not a.eos_table:
        raise SystemExit('raised x1max needs --eos-table (T of the old top cell)')
    # horizontal: per-panel cells and MeshBlock size, old -> new
    n2o, n3o = int(par['mesh']['nx2']), int(par['mesh']['nx3'])
    n2n, n3n = int(parn['mesh']['nx2']), int(parn['mesh']['nx3'])
    b2n, b3n = int(parn['meshblock']['nx2']), int(parn['meshblock']['nx3'])
    if n2n % n2o or n3n % n3o or n2n//n2o != n3n//n3o:
        raise SystemExit('horizontal: nx2, nx3 must grow by the same integer factor')
    if int(parn['meshblock']['nx1']) != nx1n:
        raise SystemExit('one MeshBlock in r is required (meshblock/nx1 = mesh/nx1)')
    if n2n % b2n or n3n % b3n:
        raise SystemExit('meshblock nx2/nx3 must divide mesh nx2/nx3')
    fac = n2n//n2o
    same_layout = (fac == 1 and b2n == h['nx2'] and b3n == h['nx3'])
    identity = (nx1n == nx1o and np.array_equal(cn, co) and same_layout)

    g = float(par['problem']['grav'])
    apl = float(par['problem']['ap'])
    etg = par['hydro'].get('etotgrav', 'false').lower() in ('true', '1')
    rotp = par['problem'].get('rot_potential', 'false').lower() in ('true', '1')
    om2 = float(par['problem'].get('omega', '0'))**2 if rotp else 0.0
    if par['problem'].get('stellar_tide', 'false').lower() in ('true', '1') and etg:
        raise SystemExit('stellar_tide with etotgrav: tidal potential not handled')

    eo = edges(nx1o, r0, r1, co, ng)            # nx1o + 2 ng + 1 edges
    ro = centroids(eo)
    en = edges(nx1n, r0, r1n, cn, 0)
    rn = centroids(en)
    # new grid WITH ghosts: the average <-> centroid conversion on the new grid uses the
    # same stencil context as on the old one (gate a: own grid -> round-off)
    eng = edges(nx1n, r0, r1n, cn, ng)
    rng = centroids(eng)
    dxng = np.diff(eng)
    actn = slice(ng, ng + nx1n)
    dvo = (eo[1:]**3 - eo[:-1]**3)/3.0
    dvn = (en[1:]**3 - en[:-1]**3)/3.0
    if ext:
        # volume of each new cell below the old top r1 (the straddling cell split by
        # volume); budgets "below the old top" and every rescaling use dvb
        enc = np.minimum(en, r1)
        dvb = (enc[1:]**3 - enc[:-1]**3)/3.0
        ia = slice(0, ng + nx1o)                 # bottom ghosts + active (no top ghosts)
        rtop = ro[ng + nx1o - 1]                 # old top-cell centroid
        rsyn = rng[rng > rtop]                    # incl. the new top ghosts
        eos = None
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import dhjcs                              # noqa: E402
        eos = dhjcs.EOS(a.eos_table)
        print('raised top: x1max %.6e -> %.6e; new cells with centroid above the old top '
              'centroid %.6e: %d (synthetic points %d); straddling cell %d (%.3f of its '
              'volume below)' % (r1, r1n, rtop, int((rn > rtop).sum()), rsyn.size,
                                 int(np.nonzero((dvb > 0) & (dvb < dvn))[0][0])
                                 if ((dvb > 0) & (dvb < dvn)).any() else -1,
                                 (dvb/dvn)[(dvb > 0) & (dvb < dvn)][0]
                                 if ((dvb > 0) & (dvb < dvn)).any() else 1.0))
    else:
        dvb = dvn
        ia = slice(None)
    s3 = slice(ng, -ng)
    act = slice(ng, ng + nx1o)
    phig_o = g*apl*(1.0 - apl/ro)
    phig_n = g*apl*(1.0 - apl/rn)
    if h['eint'] is None:
        raise SystemExit('restart without the EINTRST1 slab (not a general-EOS run)')

    # ---------------------------------------------------------- old state checks
    u = h['u'][:, :, s3, s3, :]                 # (nmb, 5, nx3, nx2, O1)
    rho, m1, m2, m3, E = (u[:, n] for n in range(5))
    ei = h['eint'][:, s3, s3, :]
    geo = [panel_cols(ll[2]*h['nx3'], (ll[2] + 1)*h['nx3'], ll[1]*h['nx2'],
                      (ll[1] + 1)*h['nx2'], n2o, n3o, ll[4]) for ll in h['lloc']]
    dom = np.stack([gg[0] for gg in geo])
    lzo = [lz_coef(ll[2]*h['nx3'], (ll[2] + 1)*h['nx3'], ll[1]*h['nx2'],
                   (ll[1] + 1)*h['nx2'], n2o, n3o, ll[4]) for ll in h['lloc']]
    A2o = np.stack([z[0] for z in lzo])
    A3o = np.stack([z[1] for z in lzo])
    cosc = np.stack([gg[1] for gg in geo])
    Bana = -0.5*om2*np.stack([gg[2] for gg in geo])
    ke = kinetic(rho, m1, m2, m3, cosc[..., None])
    if etg:
        phi_c = (E - ke - ei)/rho
        w = ro[act]**2
        B = ((phi_c[..., act] - phig_o[act])*w).sum(-1)/(w*w).sum()
        resid = phi_c[..., act] - phig_o[act] - B[..., None]*w
        print('Phi fit: max |Phi_code - (g ap (1-ap/r) + B r^2)| / |Phi| = %.2e; '
              'max |B_fit - B_analytic| / (Omega^2/2) = %.2e (rot_potential %s)'
              % (np.abs(resid/phi_c[..., act]).max(),
                 np.abs(B - Bana).max()/max(0.5*om2, 1e-300), rotp))
    else:
        B = np.zeros(rho.shape[:-1])
    if not same_layout:
        B = Bana                                  # the analytic B of TotPotAt
    phio = phig_o[None, None, None, :] + B[..., None]*ro**2
    rel_c2p = np.abs(E - ke - ei - (rho*phio if etg else 0.0))[..., act] / ei[..., act]
    print('state check: max |E - KE - e - rho Phi| / e (active) = %.2e' % rel_c2p.max())

    # ---------------------------------------------------------- budgets
    KEYS = ('mass', 'mom1', 'mom2', 'mom3', 'lz', 'eint', 'ekin', 'egrav', 'etot')

    def budget(r_, a1, a2, a3, e_, c_, B_, dv, rc, pg, dom_, A2_, A3_):
        k_ = kinetic(r_, a1, a2, a3, c_[..., None])
        eg = r_*(pg + B_[..., None]*rc**2) if etg else 0.0*r_
        lz_ = rc*(A2_[..., None]*a2 + A3_[..., None]*a3)
        cols = dict(mass=r_, mom1=a1, mom2=a2, mom3=a3, lz=lz_, eint=e_, ekin=k_,
                    egrav=eg, etot=e_ + k_ + eg)
        return {k: ((v*dv).sum(-1)*dom_).sum() for k, v in cols.items()}
    bo = budget(rho[..., act], m1[..., act], m2[..., act], m3[..., act], ei[..., act],
                cosc, B, dvo[act], ro[act], phig_o[act], dom, A2o, A3o)
    bh = dict.fromkeys(KEYS, 0.0)
    bn = dict.fromkeys(KEYS, 0.0)          # new grid, below the old top (dvb)
    ba = dict.fromkeys(KEYS, 0.0)          # new grid, above the old top (dvn - dvb)
    colmax = dict(mass=0.0, eint=0.0, mom1=0.0)

    # ---------------------------------------------------------- new layout
    nb2, nb3 = n2n//b2n, n3n//b3n
    if same_layout:
        lloc_new = [tuple(int(x) for x in ll) for ll in h['lloc']]
    else:
        lev = root_level((1, nb2, nb3))
        lloc_new = zorder([(0, j, k, lev, p) for p in range(6) for k in range(nb3)
                           for j in range(nb2)], lev)
        # the ordering rule reproduces the restart's own list
        lev_o = root_level((1, n2o//h['nx2'], n3o//h['nx3']))
        chk = zorder([tuple(int(x) for x in ll) for ll in h['lloc']], lev_o)
        if lev_o != h['lloc'][0, 3] or \
                chk != [tuple(int(x) for x in ll) for ll in h['lloc']]:
            raise SystemExit('Z-order / root level rule does not reproduce the restart')
    nmbn = len(lloc_new)
    print('horizontal: C%d -> C%d (factor %d), MeshBlocks %d x %d x %d -> %d x %d x %d '
          '(%d -> %d); henergy %s'
          % (n2o, n2n, fac, nx1o, h['nx2'], h['nx3'], nx1n, b2n, b3n, nmb, nmbn,
             a.henergy))

    dfl = float(parn['hydro'].get('dfloor', '0'))
    pfl = float(parn['hydro'].get('pfloor', '0'))
    ldmin = float(parn['hydro'].get('eos_logd_min', '-14'))
    if fac > 1 or not same_layout:
        names = ('rho', 'm1', 'm2', 'm3', 'ei', 'p', 'E')
        full = dict(rho=h['u'][:, 0], m1=h['u'][:, 1], m2=h['u'][:, 2], m3=h['u'][:, 3],
                    ei=h['eint'], p=h['p'], E=h['u'][:, 4])
        PA = {k: panel_arrays(h, full[k], n2o, n3o) for k in names}
        PU = {k: PA[k]/PA['rho'] for k in ('m1', 'm2', 'm3')}
        PS = {k: PA[k]/PA['rho'] for k in ('ei', 'p')}
        hlz_max = [0.0]
        nsc_tot = dict(rho=0, ei=0, p=0)
        fstat = [np.inf, -np.inf, 0, 0.0]
        fhist = np.zeros(nx1o + 2*ng, int)

    method = a.method
    radial_identity = (nx1n == nx1o and np.array_equal(cn, co))
    if radial_identity and not a.interp_identity:
        method = 'copy'
    dxo = np.diff(eo)

    def remap_ln(q, syn=None):
        L = np.log(q[..., ia])
        x = ro[ia]
        if method == 'avg':
            # centroid values Lc with Lc + avg_corr(Lc) = L (fixed point), i.e. the exact
            # inverse of the conversion applied on the new grid below: an own-grid remap
            # then returns the old averages to round-off (gate a)
            La, dxa = L, dxo[ia]
            L = La - avg_corr(La, x, dxa)
            for _ in range(AVG_ITMAX):
                Lk = La - avg_corr(L, x, dxa)
                dmax = np.abs(Lk - L).max()
                L = Lk
                if dmax < 1e-15*max(1.0, np.abs(La).max()):
                    break
            avg_it[0] = max(avg_it[0], _ + 1)
            avg_it[1] = max(avg_it[1], dmax)
        if syn is not None:                     # raised top: centroid values above
            L = np.concatenate((L, np.log(syn)), axis=-1)
            x = np.concatenate((x, rsyn))
        Ln = PchipInterpolator(x, L, axis=-1)(rng)
        if method == 'avg':
            Ln = Ln + avg_corr(Ln, rng, dxng)
        return np.exp(Ln[..., actn])

    avg_it = [0, 0.0]
    fill = dict(ncol=0, nhot=0, tmin=np.inf, tmax=-np.inf, fmax=1.0, ptop=[])

    def top_fill(c_rho, c_ei, c_p, c_m2, c_m3, c_B):
        """hydrostatic isothermal column above the old top-cell centroid at the points
        rsyn: T = max(T_top, tfill_min) at the top cell's mu (p/rho x fT), p continuous
        at rtop, e/p of the top cell, u_r = 0, u_h of the top cell"""
        it = ng + nx1o - 1
        rt, et, pt = c_rho[..., it], c_ei[..., it], c_p[..., it]
        Tt = eos_temperature(eos, rt, et)
        fT = np.maximum(1.0, a.tfill_min/Tt)
        c2 = (pt/rt*fT)[..., None]
        dphi = g*apl*apl*(1.0/rtop - 1.0/rsyn) + c_B[..., None]*(rsyn**2 - rtop**2)
        p_ = pt[..., None]*np.exp(-dphi/c2)
        fill['ncol'] += Tt.size
        fill['nhot'] += int((fT > 1.0).sum())
        fill['tmin'] = min(fill['tmin'], Tt.min())
        fill['tmax'] = max(fill['tmax'], Tt.max())
        fill['fmax'] = max(fill['fmax'], fT.max())
        fill['ptop'].append(p_[..., rsyn.size - ng - 1].ravel())   # last new centroid
        one = np.ones_like(p_)
        return dict(rho=p_/c2, p=p_, ei=p_*(et/pt)[..., None], u0=0.0*one,
                    u1=(c_m2/c_rho)[..., it, None]*one,
                    u2=(c_m3/c_rho)[..., it, None]*one)

    fmr, fer = [np.inf, -np.inf], [np.inf, -np.inf]
    ncl = ncl_b = npf = noff = noff_a = 0
    noff_c = np.zeros(nx1n, int)
    hse_o, hse_n = [], []
    rmin, emin, pmin, pmax = np.inf, np.inf, np.inf, -np.inf

    def hse(rr, pp, r, B_):
        rf = 0.5*(r[1:] + r[:-1])
        geff = g*apl**2/rf**2 + 2.0*B_[..., None]*rf
        dpdr = np.diff(pp, axis=-1)/np.diff(r)
        rhof = 0.5*(rr[..., 1:] + rr[..., :-1])
        return dpdr/(-rhof*geff) - 1.0, np.sqrt(pp[..., 1:]*pp[..., :-1])
    for m in range(nmb):
        if m % a.hse_sample == 0:
            hse_o.append(hse(rho[m][..., act], h['p'][m, s3, s3, act], ro[act], B[m]))

    fo = None
    if not a.check_only:
        os.makedirs(a.outdir, exist_ok=True)
        fn = os.path.join(a.outdir, 'remap.dat')
        fo = open(fn, 'wb')
        fo.write(b'RADREMPE')
        fo.write(struct.pack('<4i', nmbn, nx1n, b2n, b3n))

    # horizontal L_z restoration, 'shell' (default): one pre-pass over the new MeshBlocks
    # collects, per old radial cell, the L_z deficit of the prolonged children and
    # sum rho sin^2 dOmega; the main pass adds the uniform rotation increment that
    # restores L_z of every radial shell (all six panels) to round-off
    passes = (0, 1) if (a.hlz == 'shell' and fac > 1) else (1,)
    lz_def = np.zeros(nx1o + 2*ng)
    lz_den = np.zeros(nx1o + 2*ng)
    dOm_sh = None
    for pass_, mn, ll in [(q_, i_, l_) for q_ in passes
                          for i_, l_ in enumerate(lloc_new)]:
        _, lx2, lx3, _, pn = ll
        if pass_ == 1 and mn == 0 and passes[0] == 0:
            dOm_sh = lz_def/lz_den
            for k in nsc_tot:
                nsc_tot[k] = 0
        # ---- horizontal: this block's columns on the OLD radial grid (O1, ghosts incl.)
        if same_layout and fac == 1:
            c_rho, c_m1, c_m2, c_m3 = rho[mn], m1[mn], m2[mn], m3[mn]
            c_ei, c_p = ei[mn], h['p'][mn, s3, s3, :]
            c_cos, c_B, c_dom = cosc[mn], B[mn], dom[mn]
            c_A2, c_A3 = A2o[mn], A3o[mn]
        else:
            J0, K0 = lx2*b2n, lx3*b3n
            pj0, pj1 = J0//fac, -(-(J0 + b2n)//fac)
            pk0, pk1 = K0//fac, -(-(K0 + b3n)//fac)
            hp = HProlong(fac, pk0, pk1, pj0, pj1, n2o, n3o, pn)
            sl = (slice(K0 - pk0*fac, K0 - pk0*fac + b3n),
                  slice(J0 - pj0*fac, J0 - pj0*fac + b2n))
            ch = {}
            kp, jp = hp.shape
            Wp = hp.w4.sum(axis=(1, 3))[..., None]
            pidx = (pn, slice(1 + pk0, 1 + pk1), slice(1 + pj0, 1 + pj1))
            for k in (('rho', 'ei', 'p') if a.hthermal == 'density' else ('rho',)):
                c, sj, sk = slopes(PA[k], pn, pk0, pk1, pj0, pj1)
                ch[k], nsc = hp.children(c, sj, sk, positive=True,
                                         qmin=(dfl if k == 'rho' else 0.0))
                nsc_tot[k] += nsc
            rW = hp.parent_sum(ch['rho'])
            if a.hthermal == 'specific':
                # e and p through the SPECIFIC values e/rho, p/rho (MC-limited linear,
                # bounded by the neighbours, i.e. T is) plus one shift per parent that
                # restores sum e dOmega (sum p dOmega); e, p and rho prolonged
                # independently put T far below every neighbour into children where rho
                # rises and e falls (night top: T_child 72 K next to 200 K parents)
                for k in ('ei', 'p'):
                    c, sj, sk = slopes(PS[k], pn, pk0, pk1, pj0, pj1)
                    sc, _ = hp.children(c, sj, sk)
                    sc = sc.reshape(kp, fac, jp, fac, -1)
                    ds = (PA[k][pidx]*Wp - hp.parent_sum(
                        ch['rho']*sc.reshape(kp*fac, jp*fac, -1)))/rW
                    q = sc + ds[:, None, :, None, :]
                    bad = (q.min(axis=(1, 3)) < 0.5*c)            # (kp, jp, O1)
                    nsc_tot[k] += int(bad.sum())
                    q = np.where(bad[:, None, :, None, :], c[:, None, :, None, :], q)
                    ch[k] = ch['rho']*q.reshape(kp*fac, jp*fac, -1)
            # momentum through the velocity u = m/rho (MC-limited linear, bounded by the
            # neighbours' u) plus one shift per parent that restores sum m dOmega; m and
            # rho prolonged independently put huge u into low-density children
            for k in ('m1', 'm2', 'm3'):
                c, sj, sk = slopes(PU[k], pn, pk0, pk1, pj0, pj1)
                uc, _ = hp.children(c, sj, sk)
                mp = PA[k][pn, 1 + pk0:1 + pk1, 1 + pj0:1 + pj1]*hp.w4.sum(axis=(1, 3))[
                    ..., None]
                du = (mp - hp.parent_sum(ch['rho']*uc))/rW
                ch[k] = ch['rho']*(uc.reshape(kp, fac, jp, fac, -1) +
                                   du[:, None, :, None, :]).reshape(kp*fac, jp*fac, -1)
            if a.hlz != 'none' and fac > 1:
                # L_z per parent and radial cell: the stored components are conserved
                # above, but the tangent basis turns across a parent, so sum L_z of the
                # children differs from the parent's (-2e-4 of L_z at C32 -> C128).
                # Restore it with a uniform rotation increment dOm about z per parent:
                # dm = rho_c dOm (w1, w2), the covariant components of z x rhat, which
                # add dL_z = r rho_c dOm sin^2(theta_c) (r common to the parent)
                A2p, A3p = lz_coef(pk0, pk1, pj0, pj1, n2o, n3o, pn)
                A2c, A3c = lz_coef(pk0*fac, pk1*fac, pj0*fac, pj1*fac, n2o*fac,
                                   n3o*fac, pn)
                lzp = (A2p[..., None]*PA['m2'][pidx] + A3p[..., None]*PA['m3'][pidx])*Wp
                lzc = hp.parent_sum(A2c[..., None]*ch['m2'] + A3c[..., None]*ch['m3'])
                den = hp.parent_sum(ch['rho']*hp.sin2[..., None])
                if a.hlz == 'parent':
                    dOm = (lzp - lzc)/den
                    dOm = np.repeat(np.repeat(dOm, fac, axis=0), fac, axis=1)
                elif pass_ == 0:
                    if (b2n % fac) or (b3n % fac):
                        raise SystemExit('--hlz shell: meshblock nx2, nx3 % factor != 0')
                    lz_def += (lzp - lzc).sum(axis=(0, 1))
                    lz_den += den.sum(axis=(0, 1))
                    continue
                else:
                    dOm = dOm_sh[None, None, :]
                w1c, w2c = zrot_cov(pk0*fac, pk1*fac, pj0*fac, pj1*fac, n2o*fac,
                                    n3o*fac, pn)
                ch['m2'] = ch['m2'] + ch['rho']*dOm*w1c[..., None]
                ch['m3'] = ch['m3'] + ch['rho']*dOm*w2c[..., None]
                dv = np.abs(dOm*np.sqrt(hp.sin2)[..., None])[..., act]
                hlz_max[0] = max(hlz_max[0], dv.max())
            B_f = -0.5*om2*hp.sin2
            if a.henergy == 'total' and fac > 1:
                # e x f per parent so that sum_children (e + KE + rho Phi) dOmega is the
                # parent's E dOmega (active radial cells; radial ghosts keep f = 1)
                kef = kinetic(ch['rho'], ch['m1'], ch['m2'], ch['m3'], hp.cosc[..., None])
                rpf = ch['rho']*(phig_o + B_f[..., None]*ro**2) if etg else 0.0*kef
                Ep = PA['E'][pn, 1 + pk0:1 + pk1, 1 + pj0:1 + pj1]
                Wp = hp.w4.sum(axis=(1, 3))[..., None]
                ff = (Ep*Wp - hp.parent_sum(kef + rpf))/hp.parent_sum(ch['ei'])
                ff[..., :ng] = 1.0
                ff[..., ng + nx1o:] = 1.0
                fstat[0] = min(fstat[0], ff.min())
                fstat[1] = max(fstat[1], ff.max())
                fstat[2] += int((ff < 0.0).sum())
                off = np.abs(ff - 1.0) > FTOL
                fhist += off.sum(axis=(0, 1))
                if off.any():
                    fstat[3] = max(fstat[3], PA['p'][pn, 1 + pk0:1 + pk1,
                                                     1 + pj0:1 + pj1][off].max())
                ff = np.clip(ff, 1.0 - FTOL, 1.0 + FTOL)
                ch['ei'] = (ch['ei'].reshape(kp, fac, jp, fac, -1) *
                            ff[:, None, :, None, :]).reshape(kp*fac, jp*fac, -1)
            c_rho, c_m1, c_m2, c_m3 = (ch[k][sl] for k in ('rho', 'm1', 'm2', 'm3'))
            c_ei, c_p = ch['ei'][sl], ch['p'][sl]
            c_cos, c_B = hp.cosc[sl], B_f[sl]
            c_dom = panel_cols(K0, K0 + b3n, J0, J0 + b2n, n2n, n3n, pn)[0]
            c_A2, c_A3 = lz_coef(K0, K0 + b3n, J0, J0 + b2n, n2n, n3n, pn)
        bb = budget(c_rho[..., act], c_m1[..., act], c_m2[..., act], c_m3[..., act],
                    c_ei[..., act], c_cos, c_B, dvo[act], ro[act], phig_o[act], c_dom,
                    c_A2, c_A3)
        for k in KEYS:
            bh[k] += bb[k]
        # ---- radial, per column (the dhj-remap method, unchanged)
        if method == 'copy':
            rn_, m1n, m2n, m3n = (q[..., act] for q in (c_rho, c_m1, c_m2, c_m3))
            ein, pn_ = c_ei[..., act], c_p[..., act]
        else:
            syn = None
            if ext:
                syn = top_fill(c_rho, c_ei, c_p, c_m2, c_m3, c_B)
            rr = remap_ln(c_rho, None if syn is None else syn['rho'])
            ee = remap_ln(c_ei, None if syn is None else syn['ei'])
            pp = remap_ln(c_p, None if syn is None else syn['p'])
            inb = dvb > 0.0                       # new cells (partly) below the old top
            f_m = (c_rho[..., act]*dvo[act]).sum(-1)/(rr*dvb).sum(-1)
            f_e = (c_ei[..., act]*dvo[act]).sum(-1)/(ee*dvb).sum(-1)
            fmr = [min(fmr[0], f_m.min()), max(fmr[1], f_m.max())]
            fer = [min(fer[0], f_e.min()), max(fer[1], f_e.max())]
            rn_ = np.where(inb, f_m[..., None], 1.0)*rr
            ein = np.where(inb, f_e[..., None], 1.0)*ee
            pn_ = np.where(inb, f_e[..., None], 1.0)*pp
            mom = []
            for nq, q in enumerate((c_m1, c_m2, c_m3)):
                uo = q/c_rho
                if syn is None:
                    un = PchipInterpolator(ro, uo, axis=-1)(rn)
                else:
                    un = PchipInterpolator(np.concatenate((ro[ia], rsyn)),
                                           np.concatenate((uo[..., ia], syn['u%d' % nq]),
                                                          axis=-1), axis=-1)(rn)
                # m1: column momentum; m2, m3: column r x momentum (L_z per column),
                # both below the old top
                wo = 1.0 if nq == 0 else ro[act]
                wn = 1.0 if nq == 0 else rn
                du = ((q[..., act]*wo*dvo[act]).sum(-1) - (rn_*un*wn*dvb).sum(-1)) / \
                    (rn_*wn*dvb).sum(-1)
                mom.append(rn_*(un + np.where(inb, du[..., None], 0.0)))
            m1n, m2n, m3n = mom
        nc = int((rn_ < dfl).sum())
        if nc:
            # raise to dfloor keeping v and T (e and p scale with rho)
            s_ = np.maximum(dfl/rn_, 1.0)
            rn_, m1n, m2n, m3n = rn_*s_, m1n*s_, m2n*s_, m3n*s_
            ein, pn_ = ein*s_, pn_*s_
            ncl += nc
            ncl_b += int(((s_ > 1.0) & (dvb > 0.0)).sum())
        npf += int((pn_ < pfl).sum())
        noff += int((rn_ < 10.0**ldmin).sum())
        noff_a += int(((rn_ < 10.0**ldmin) & (dvb <= 0.0)).sum())
        noff_c += np.count_nonzero(rn_ < 10.0**ldmin, axis=(0, 1))
        bb = budget(rn_, m1n, m2n, m3n, ein, c_cos, c_B, dvb, rn, phig_n, c_dom, c_A2,
                    c_A3)
        for k in KEYS:
            bn[k] += bb[k]
        if ext:
            bb = budget(rn_, m1n, m2n, m3n, ein, c_cos, c_B, dvn - dvb, rn, phig_n,
                        c_dom, c_A2, c_A3)
            for k in KEYS:
                ba[k] += bb[k]
        for k, q_o, q_n in (('mass', c_rho, rn_), ('eint', c_ei, ein),
                            ('mom1', c_m1, m1n)):
            co_ = (q_o[..., act]*dvo[act]).sum(-1)
            cn_ = (q_n*dvb).sum(-1)
            ref = np.abs(co_) if k != 'mom1' else (c_rho[..., act]*dvo[act]).sum(-1)*1e5
            colmax[k] = max(colmax[k], np.abs((cn_ - co_)/ref).max())
        if mn % a.hse_sample == 0:
            hse_n.append(hse(rn_, pn_, rn, c_B))
        rmin, emin = min(rmin, rn_.min()), min(emin, ein.min())
        pmin, pmax = min(pmin, pn_.min()), max(pmax, pn_.max())
        if not (np.isfinite(rn_).all() and np.isfinite(ein).all() and
                np.isfinite(m1n).all() and np.isfinite(m2n).all() and
                np.isfinite(m3n).all()):
            raise SystemExit('non-finite value in new MeshBlock %d %s' % (mn, ll))
        if fo is not None:
            fo.write(np.array(ll, '<i4').tobytes())
            rec = np.stack([rn_, m1n/rn_, m2n/rn_, m3n/rn_, ein])
            fo.write(np.ascontiguousarray(rec, '<f8').tobytes())
    if fo is not None:
        fo.close()

    # ---------------------------------------------------------- report
    print('rst %s: t = %.6e s, cycle %d, sections %s'
          % (a.rst, h['t'], h['ncycle'], h['sections']))
    print('grid nx1 %d -> %d (%s); method %s; dr min/max new %.3e / %.3e cm'
          % (nx1o, nx1n, 'identity' if identity else 'remap', method,
             np.diff(en).min(), np.diff(en).max()))
    if fac > 1 or not same_layout:
        print('horizontal: parents with slopes scaled for positivity: rho %d, e %d, p %d'
              % (nsc_tot['rho'], nsc_tot['ei'], nsc_tot['p']))
        if a.hlz != 'none' and fac > 1:
            print('  L_z restoration (%s): max added |v_phi| %.3e cm/s (active cells)'
                  % (a.hlz, hlz_max[0]))
            if dOm_sh is not None:
                print('  shell rotation increments dOmega = dv/(r sin): |max| %.3e rad/s '
                      '(Omega %.3e)' % (np.abs(dOm_sh[act]/ro[act]).max(), np.sqrt(om2)))
        if a.henergy == 'total' and fac > 1:
            print('  energy factor f (e x f per parent) %.6f..%.6f (%d < 0)'
                  % (fstat[0], fstat[1], fstat[2]))
        if fhist.any():
            ii = np.nonzero(fhist)[0]
            print('  parents with |f - 1| > %g (f clipped there, E not conserved): %d, '
                  'old radial cells %d..%d, highest parent p %.3e bar'
                  % (FTOL, fhist.sum(), ii[0] - ng, ii[-1] - ng, fstat[3]/BAR))
    if method == 'avg':
        print('average -> centroid inversion: max %d iterations, last change %.1e'
              % (avg_it[0], avg_it[1]))
    if method != 'copy':
        print('column factors: f_m %.6f..%.6f, f_e %.6f..%.6f (1 = the interpolation '
              'alone conserved it)' % (fmr[0], fmr[1], fer[0], fer[1]))
    print('cells raised to dfloor = %g: %d (of them below the old top: %d); cells with '
          'p < pfloor = %g: %d' % (dfl, ncl, ncl_b, pfl, npf))
    print('cells below the EOS table (rho < 1e%g): %d of %d (above the old top: %d); '
          'per new radial cell i: %s' % (ldmin, noff, nmbn*b2n*b3n*nx1n, noff_a,
                                         ', '.join('%d:%d' % (i, n) for i, n in
                                                   enumerate(noff_c) if n)))
    if ext:
        pt_ = np.concatenate(fill['ptop'])
        print('top fill: %d columns, T_top %.0f..%.0f K, %d columns (%.2f %%) raised to '
              '%g K (max factor %.2f); p at the top centroid: max %.3e, median %.3e bar'
              % (fill['ncol'], fill['tmin'], fill['tmax'], fill['nhot'],
                 100.0*fill['nhot']/fill['ncol'], a.tfill_min, fill['fmax'],
                 pt_.max()/BAR, np.median(pt_)/BAR))
    print('%-7s %16s %16s %16s %10s %10s %16s' % ('budget', 'old', 'horizontal',
                                                  'new<oldtop', 'h/old-1', 'new/h-1',
                                                  'new>oldtop'))
    for k in KEYS:
        print('%-7s %16.9e %16.9e %16.9e %+10.2e %+10.2e %16.9e'
              % (k, bo[k], bh[k], bn[k], bh[k]/bo[k] - 1 if bo[k] else 0.0,
                 bn[k]/bh[k] - 1 if bh[k] else 0.0, ba[k]))
    print('  (momenta: sums of the stored covariant components; relative to |mom| they '
          'are near-cancelling sums: vs mass x 1 km/s: h %.1e %.1e %.1e)'
          % tuple(abs(bh[k] - bo[k])/(abs(bo['mass'])*1e5) for k in
                  ('mom1', 'mom2', 'mom3')))
    print('  radial step, max |column change| (mass, eint rel; mom1 vs column mass x 1 '
          'km/s): %.2e %.2e %.2e' % (colmax['mass'], colmax['eint'], colmax['mom1']))
    print('HSE residual |dp/dr / (-rho g_eff) - 1| at faces, median / p90 (old -> new, '
          'every %d-th MeshBlock)' % a.hse_sample)
    ho = np.concatenate([x[0].ravel() for x in hse_o])
    po = np.concatenate([x[1].ravel() for x in hse_o])
    hn = np.concatenate([x[0].ravel() for x in hse_n])
    pnf = np.concatenate([x[1].ravel() for x in hse_n])
    for lo, hi in ((1e-9, 1e-6), (1e-6, 1e-4), (1e-4, 1e-2), (1e-2, 1.0), (1.0, 10.0),
                   (10.0, 100.0), (100.0, 1e4)):
        x = np.abs(ho[(po >= lo*BAR) & (po < hi*BAR)])
        y = np.abs(hn[(pnf >= lo*BAR) & (pnf < hi*BAR)])
        if x.size and y.size:
            print('  %8.0e-%-8.0e bar  %8.2e / %8.2e  ->  %8.2e / %8.2e' % (
                lo, hi, np.median(x), np.percentile(x, 90), np.median(y),
                np.percentile(y, 90)))
    print('rho min %.3e -> %.3e; e min %.3e -> %.3e; p range new %.3e..%.3e bar'
          % (rho[..., act].min(), rmin, ei[..., act].min(), emin, pmin/BAR, pmax/BAR))
    if a.check_only:
        return
    txt = h['text']
    txt = txt[:txt.find('<par_end>')]
    for k, v in ov.items():
        blk, key = k.split('/')
        txt = set_key(txt, blk, key, v, 'remap target grid')
    txt = set_key(txt, 'time', 'start_time', repr(h['t']),
                  'remap: t of %s' % os.path.basename(a.rst))
    txt = set_key(txt, 'problem', 'remap_file', os.path.abspath(fn), 'remap')
    with open(os.path.join(a.outdir, 'remap.athinput'), 'w') as fo:
        fo.write(txt)
    print('wrote', fn, 'and', os.path.join(a.outdir, 'remap.athinput'))


if __name__ == '__main__':
    sys.exit(main())
