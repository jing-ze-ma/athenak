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
The horizontal grid, the MeshBlock layout and x1min/x1max must be those of the restart.

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
    m = par['mesh']
    c = [float(m.get('f_stretch_r_c%d' % k, 0.0)) for k in range(1, NSTR+1)]
    poly = m.get('use_grid_stretch_r_poly', 'false').lower() in ('true', '1')
    if not poly or m.get('use_grid_stretch_r', 'false').lower() in ('true', '1'):
        raise SystemExit('only the polynomial radial stretch is supported')
    return (int(m['nx1']), float(m['x1min']), float(m['x1max']), np.array(c))


def edges(nx1, r0, r1, c, ng=0):
    """cell edges as Coordinates builds them (LeftEdgeX + StretchRPoly), ghosts incl."""
    i = np.arange(-ng, nx1+ng+1)
    x = (i/nx1)*r1 + (1.0 - i/nx1)*r0          # LeftEdgeX
    xi = (x - r0)/(r1 - r0)
    u = xi.copy()
    xik = xi.copy()
    for k in range(NSTR):
        u += c[k]*xik*(1.0 - xi)
        xik = xik*xi
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


