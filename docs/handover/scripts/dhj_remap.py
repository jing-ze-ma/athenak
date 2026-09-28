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


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument('rst')
    ap_.add_argument('outdir')
    ap_.add_argument('--grid', default='')
    ap_.add_argument('--method', default='avg', choices=('avg', 'point'))
    ap_.add_argument('--interp-identity', action='store_true',
                     help='run the interpolation even when the grid is unchanged')
    ap_.add_argument('--check-only', action='store_true')
    a = ap_.parse_args()

    h = read_rst(a.rst)
    par = h['par']
    ng, nmb = h['ng'], h['nmb']
    nx1o, r0, r1, co = grid_of(par)
    if nx1o != h['nx1'] or int(par['meshblock']['nx1']) != nx1o:
        raise SystemExit('one MeshBlock in r is required')
    # target grid
    ov = dict(kv.split('=', 1) for kv in a.grid.split())
    parn = {b: dict(v) for b, v in par.items()}
    for k, v in ov.items():
        blk, key = k.split('/')
        parn.setdefault(blk, {})[key] = v
    nx1n, r0n, r1n, cn = grid_of(parn)
    if abs(r0n/r0 - 1) > 1e-12 or abs(r1n/r1 - 1) > 1e-12:
        raise SystemExit('x1min/x1max must not change (%g %g -> %g %g)'
                         % (r0, r1, r0n, r1n))
    identity = (nx1n == nx1o and np.array_equal(cn, co))

    g = float(par['problem']['grav'])
    apl = float(par['problem']['ap'])
    etg = par['hydro'].get('etotgrav', 'false').lower() in ('true', '1')

    eo = edges(nx1o, r0, r1, co, ng)            # nx1o + 2 ng + 1 edges
    ro = centroids(eo)
    en = edges(nx1n, r0, r1, cn, 0)
    rn = centroids(en)
    dvo = (eo[1:]**3 - eo[:-1]**3)/3.0
    dvn = (en[1:]**3 - en[:-1]**3)/3.0
    s3 = slice(ng, -ng)
    act = slice(ng, ng + nx1o)

    u = h['u'][:, :, s3, s3, :]                 # (nmb, 5, nx3, nx2, O1)
    rho, m1, m2, m3, E = (u[:, n] for n in range(5))
    p = h['p'][:, s3, s3, :]
    cosc, dom = cs_columns(h)
    c3 = cosc[..., None]
    ke = kinetic(rho, m1, m2, m3, c3)
    if h['eint'] is None:
        raise SystemExit('restart without the EINTRST1 slab (not a general-EOS run)')
    ei = h['eint'][:, s3, s3, :]
    # rho*Phi as the code has it, and the per-column centrifugal coefficient B:
    # Phi = g ap (1 - ap/r) + B r^2  (TotPotAt, grav_point_mass; B = -Omega^2 sin^2/2)
    phig = g*apl*(1.0 - apl/ro)
    if etg:
        phi_c = (E - ke - ei)/rho                # the code's Phi (+ c2p round-off)
        w = ro[act]**2
        B = ((phi_c[..., act] - phig[act])*w).sum(-1)/(w*w).sum()
        resid = phi_c[..., act] - phig[act] - B[..., None]*w
        print('Phi fit: max |Phi_code - (g ap (1-ap/r) + B r^2)| / |Phi| = %.2e '
              '(B/(-Omega^2/2) = sin^2 theta in [%.4f, %.4f])'
              % (np.abs(resid/phi_c[..., act]).max(),
                 (B/(-0.5*float(par['problem']['omega'])**2)).min(),
                 (B/(-0.5*float(par['problem']['omega'])**2)).max()))
    else:
        B = np.zeros(rho.shape[:-1])
    phio = phig[None, None, None, :] + B[..., None]*ro**2
    phin = g*apl*(1.0 - apl/rn)[None, None, None, :] + B[..., None]*rn**2
    rel_c2p = np.abs(E - ke - ei - (rho*phio if etg else 0.0))[..., act] / ei[..., act]
    print('state check: max |E - KE - e - rho Phi| / e (active) = %.2e' % rel_c2p.max())

    method = a.method
    if identity and not a.interp_identity:
        method = 'copy'
        rn_, m1n, m2n, m3n = rho[..., act], m1[..., act], m2[..., act], m3[..., act]
        ein = ei[..., act]
        pn = p[..., act]
    else:
        # the restart holds cell AVERAGES; 'avg' turns them into centroid values, puts a
        # monotone cubic (PCHIP, ghosts included so nothing is extrapolated) through the
        # centroid values of ln rho, ln e (and ln p, diagnostics) and u = m/rho, and turns
        # the new centroid values back into averages of the new cells ('point' skips the
        # two conversions: the radab_0924 prototype).  Then per column: rho x f_m so the
        # column mass is exact, e x f_e so the column internal energy is exact (f_m, f_e
        # printed), u + du_k so each column momentum is exact.
        dxo = np.diff(eo)
        dxn = np.diff(en)

        def remap_ln(q):
            L = np.log(q)
            if method == 'avg':
                L = L - avg_corr(L, ro, dxo)
            Ln = PchipInterpolator(ro, L, axis=-1)(rn)
            if method == 'avg':
                Ln = Ln + avg_corr(Ln, rn, dxn)
            return np.exp(Ln)
        rr = remap_ln(rho)
        ee = remap_ln(ei)
        pp = remap_ln(p)
        f_m = (rho[..., act]*dvo[act]).sum(-1)/(rr*dvn).sum(-1)
        f_e = (ei[..., act]*dvo[act]).sum(-1)/(ee*dvn).sum(-1)
        rn_ = f_m[..., None]*rr
        ein = f_e[..., None]*ee
        pn = f_e[..., None]*pp
        mom = []
        for q in (m1, m2, m3):
            un = PchipInterpolator(ro, q/rho, axis=-1)(rn)
            du = ((q[..., act]*dvo[act]).sum(-1) - (rn_*un*dvn).sum(-1)) / \
                (rn_*dvn).sum(-1)
            mom.append(rn_*(un + du[..., None]))
        m1n, m2n, m3n = mom
        print('column factors: f_m %.6f..%.6f, f_e %.6f..%.6f (1 = the interpolation '
              'alone conserved it)' % (f_m.min(), f_m.max(), f_e.min(), f_e.max()))
    dfl = float(par['hydro'].get('dfloor', '0'))
    ncl = int((rn_ < dfl).sum())
    if ncl:
        s_ = np.maximum(dfl/rn_, 1.0)            # keep u: scale the momentum too
        rn_, m1n, m2n, m3n = rn_*s_, m1n*s_, m2n*s_, m3n*s_
    print('cells raised to dfloor = %g: %d' % (dfl, ncl))
    ken = kinetic(rn_, m1n, m2n, m3n, c3)
    if os.environ.get("REMAP_DEBUG"):
        np.savez(os.environ["REMAP_DEBUG"], rn=rn_, ein=ein, pn=pn, rho=rho[..., act],
                 ei=ei[..., act], p=p[..., act], ro=ro[act], rnc=rn, B=B)

    # ------------------------------------------------------------ budgets
    def col(q, dv):
        return (q*dv).sum(-1)
    ob = dict(mass=col(rho[..., act], dvo[act]), mom1=col(m1[..., act], dvo[act]),
              eint=col(ei[..., act], dvo[act]), ekin=col(ke[..., act], dvo[act]),
              egrav=col((rho*phio)[..., act], dvo[act]))
    nb = dict(mass=col(rn_, dvn), mom1=col(m1n, dvn), eint=col(ein, dvn),
              ekin=col(ken, dvn), egrav=col(rn_*phin, dvn))
    for d in (ob, nb):
        d['e+k'] = d['eint'] + d['ekin']
        d['etot'] = d['e+k'] + d['egrav']
    print('rst %s: t = %.6e s, cycle %d, sections %s'
          % (a.rst, h['t'], h['ncycle'], h['sections']))
    print('grid nx1 %d -> %d (%s); method %s; dr min/max new %.3e / %.3e cm'
          % (nx1o, nx1n, 'identity' if identity else 'remap', method,
             np.diff(en).min(), np.diff(en).max()))
    print('%-7s %15s %15s %11s %11s %11s' % ('budget', 'old (global)', 'new (global)',
                                             'glob rel', 'col max|rel|', 'col median'))
    for k in ('mass', 'eint', 'ekin', 'e+k', 'egrav', 'etot', 'mom1'):
        go = (ob[k]*dom).sum()
        gn = (nb[k]*dom).sum()
        ref = np.abs(ob[k]) if k != 'mom1' else np.abs(ob['mass'])*1e5
        rc = (nb[k] - ob[k])/np.where(ref > 0, ref, 1)
        print('%-7s %15.8e %15.8e %+11.3e %11.3e %+11.3e' % (k, go, gn, gn/go - 1,
                                                             np.abs(rc).max(),
                                                             np.median(rc)))
    geg = abs((ob['egrav']*dom).sum())
    print('  (etot change relative to |egrav|: global %+.3e; the mom1 column error is '
          'relative to column mass x 1 km/s)'
          % (((nb['etot'] - ob['etot'])*dom).sum()/geg))

    # HSE residual at interior faces, g_eff = dPhi/dr
    def hse(rr, pp, r):
        rf = 0.5*(r[1:] + r[:-1])
        geff = g*apl**2/rf**2 + 2.0*B[..., None]*rf
        dpdr = np.diff(pp, axis=-1)/np.diff(r)
        rhof = 0.5*(rr[..., 1:] + rr[..., :-1])
        return dpdr/(-rhof*geff) - 1.0, np.sqrt(pp[..., 1:]*pp[..., :-1])
    ho, po = hse(rho[..., act], p[..., act], ro[act])
    hn, pnf = hse(rn_, pn, rn)
    print('HSE residual |dp/dr / (-rho g_eff) - 1| at faces, median / p90 (old -> new)'
          + '  [new p = p interpolated x f_e; the code takes p from (rho, e)]')
    for lo, hi in ((1e-9, 1e-6), (1e-6, 1e-4), (1e-4, 1e-2), (1e-2, 1.0), (1.0, 10.0),
                   (10.0, 100.0), (100.0, 1e4)):
        x = np.abs(ho[(po >= lo*BAR) & (po < hi*BAR)])
        y = np.abs(hn[(pnf >= lo*BAR) & (pnf < hi*BAR)])
        if x.size and y.size:
            print('  %8.0e-%-8.0e bar  %8.2e / %8.2e  ->  %8.2e / %8.2e' % (
                lo, hi, np.median(x), np.percentile(x, 90), np.median(y),
                np.percentile(y, 90)))
    print('rho min %.3e -> %.3e; e min %.3e -> %.3e; p range new %.3e..%.3e bar'
          % (rho[..., act].min(), rn_.min(), ei[..., act].min(), ein.min(),
             pn.min()/BAR, pn.max()/BAR))
    if a.check_only:
        return

    # ------------------------------------------------------------ write
    os.makedirs(a.outdir, exist_ok=True)
    fn = os.path.join(a.outdir, 'remap.dat')
    use_p = False
    th = pn if use_p else ein
    with open(fn, 'wb') as fo:
        fo.write(b'RADREMAP' if use_p else b'RADREMPE')
        fo.write(struct.pack('<4i', nmb, nx1n, h['nx2'], h['nx3']))
        for m in range(nmb):
            fo.write(h['lloc'][m].astype('<i4').tobytes())
            rec = np.stack([rn_[m], m1n[m]/rn_[m], m2n[m]/rn_[m], m3n[m]/rn_[m], th[m]])
            fo.write(np.ascontiguousarray(rec, '<f8').tobytes())
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
