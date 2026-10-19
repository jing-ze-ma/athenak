#!/usr/bin/env python3
"""remap_e1w.py (e1remap-1014): AthenaK restart (he_star_m1 + rad_m1, general EOS, spherical polar) ->
2x-lateral narrow-wedge initial state for problem/he_ic_remap_file (pgen he_star_m1, branch e1remap-1014).

usage: python3 remap_e1w.py RST TEMPLATE_DECK OUTDIR [--x2min A --x2max B] [--dphi D] [--phi0 auto|IDX]
                            [--refine 2] [--blend W] [--name NAME]
  RST            parent restart file (single file, double precision)
  TEMPLATE_DECK  child deck (narrow wedge; its <mesh> x2min/x2max/x3min/x3max/nx2/nx3 define the child grid;
                 radial keys must equal the parent's)
  OUTDIR         writes NAME.remap (binary, see FORMAT), NAME.athinput (deck + start keys), NAME.remap.txt (diagnostics)
  --phi0 auto    parent phi offset (in parent cells) of the sub-wedge with the smallest seam mismatch; or an integer
  --blend W      seam blend half-width in CHILD cells (default 8)

Method (all per radial index i, including the x1 ghosts; the radial grid is copied):
  1. parent sub-wedge = parent cells covering the child's theta/phi range (theta must coincide with parent faces).
  2. conservative prolongation x2 laterally: q_c = q_p + s_th (thbar_c - thbar_p) + s_ph (phbar_c - phbar_p) with volume
     centroids (spherical weights cos a - cos b), MC-limited slopes, then a Barth-Jespersen clip to the parent 3x3
     neighbourhood (monotone).  Variables: rho, m1, m2, m3, eint (restart slab w0(IEN)), E, F1, F2, F3, f0x1.
     rho*Phi per parent cell from the parent (u0(IEN) - KE - eint); ET_c = eint_c + KE_c + rho_c Phi_p, and the KE
     non-linearity of each parent cell is added uniformly to its children's eint (ET exact per parent cell).
  3. periodic seams: pairwise exchange across each seam: the cells d = 0..W-1 from the seam on either side
     (equal volumes: phi uniform; theta range symmetric about pi/2) mix as q_L' = (1-a) q_L + a q_R, a = 0.5 (1 - d/W)
     -> seam jump 0, every pair sum (so every shell sum) exact.  Variables: rho, m1-3, ET, E, F1-3, f0x1.
FORMAT (little endian): 'HSREMAP1' | int32 [ver=1, ng, nx1, nx2, nx3, nhyd=5, nm1=4, 0] | float64 [t, x1min, x1max, x2min,
  x2max, x3min, x3max] | uint64 npgen | npgen bytes (parent PGEN block: Picard window + he_ic_equil a_ref trailer) |
  float64 hydro u0 [5][nx3+2ng][nx2+2ng][nx1+2ng] | m1 u0 [4][..] | f0x1 [nx3+2ng][nx2+2ng][nx1+2ng+1]
  (lateral ghosts = periodic images of the child, x1 ghosts = the parent's, prolonged).
"""
import argparse
import hashlib
import math
import os
import re
import struct
import sys

import numpy as np

MAGICS = [b'PGENST01', b'RTWARM01', b'M1PRED01', b'M1TIME2A', b'M1ONEP01', b'M1MRWIN1', b'EINTRST1', b'CKRST001',
          b'M1CTRLX1', b'M1RSSIG1', b'M1VGDAS1', b'M1CSK1S1', b'M1ROTAV1']


def parse_params(txt):
    blk, out = None, {}
    for line in txt.splitlines():
        s = line.split('#')[0].strip()
        if not s:
            continue
        m = re.match(r'<\s*([^>]+)\s*>', s)
        if m:
            blk = m.group(1).strip()
            continue
        if '=' in s and blk is not None:
            k, v = s.split('=', 1)
            out[blk + '/' + k.strip()] = v.strip()
    return out


def read_rst(fn):
    """returns dict: params, time, dt, ncycle, nmb, lloc, blocks, data_off, data_size, mesh/mb dims, slab offsets."""
    fsz = os.path.getsize(fn)
    with open(fn, 'rb') as f:
        head = f.read(min(fsz, 4 << 20))
    pe = head.find(b'<par_end>')
    if pe < 0:
        raise SystemExit('no <par_end> in ' + fn)
    pos = head.find(b'\n', pe) + 1
    ptxt = head[:pos].decode('latin-1')
    par = parse_params(ptxt)
    with open(fn, 'rb') as f:
        f.seek(pos)
        h = f.read(252)
        nmb, root_level = struct.unpack_from('<ii', h, 0)
        rsz = struct.unpack_from('<9d', h, 8)
        mi = struct.unpack_from('<19i', h, 80)
        bi = struct.unpack_from('<19i', h, 156)
        t, dt = struct.unpack_from('<dd', h, 232)
        ncyc, = struct.unpack_from('<i', h, 248)
        lloc = np.frombuffer(f.read(nmb * 20), dtype='<i4').reshape(nmb, 5)
        f.read(nmb * 4)  # cost
        blocks = {}
        order = []
        while True:
            mg = f.read(8)
            if mg in MAGICS:
                nb, = struct.unpack('<Q', f.read(8))
                blocks[mg] = f.read(nb)
                order.append(mg.decode())
                continue
            data_size, = struct.unpack('<Q', mg)
            data_off = f.tell()
            break
    if data_off + nmb * data_size != fsz:
        raise SystemExit('restart size mismatch: off %d + %d x %d != %d' % (data_off, nmb, data_size, fsz))
    ng = mi[0]
    n1, n2, n3 = bi[1] + 2 * ng, bi[2] + 2 * ng, bi[3] + 2 * ng
    nc = n1 * n2 * n3

    def i32(mg, k):
        return struct.unpack_from('<i', blocks[mg], 4 * k)[0] if mg in blocks else 0
    npred = i32(b'M1PRED01', 1)
    nt2 = i32(b'M1TIME2A', 1)
    neint = i32(b'EINTRST1', 0)
    nctr = i32(b'M1CTRLX1', 1)
    nrss = i32(b'M1RSSIG1', 1)
    nvgd = i32(b'M1VGDAS1', 1)
    ncsk = i32(b'M1CSK1S1', 1)
    nra = struct.unpack_from('<q', blocks[b'M1ROTAV1'], 8)[0] if b'M1ROTAV1' in blocks else 0
    nwarm = i32(b'RTWARM01', 0)
    if b'CKRST001' in blocks:
        raise SystemExit('implicit-ck restart state not supported')
    if neint != 1:
        raise SystemExit('expected one general-EOS internal energy slab, header says %d' % neint)
    tail = 3 + nwarm + npred + nt2 + neint + nctr + nrss + nvgd + ncsk + nra  # wtemp, p, Gamma1 + ...
    sol = None
    for nh in (5, 6, 7, 8):
        for f2 in (0, 1):
            for f3 in (0, f2):
                tot = nc * (nh + 4) + n3 * n2 * (n1 + 1) + f2 * n3 * (n2 + 1) * n1 + f3 * (n3 + 1) * n2 * n1 + tail * nc
                if tot * 8 == data_size:
                    sol = (nh, f2, f3)
    if sol is None:
        raise SystemExit('cannot decode the per-MeshBlock record (data_size %d)' % data_size)
    nh, f2, f3 = sol
    o = 0
    off = {}
    off['hyd'] = o; o += nh * nc
    off['m1'] = o; o += 4 * nc
    off['f0x1'] = o; o += n3 * n2 * (n1 + 1)
    o += f2 * n3 * (n2 + 1) * n1 + f3 * (n3 + 1) * n2 * n1
    off['wtemp'] = o; o += nc
    o += 2 * nc + nwarm * nc + npred * nc + nt2 * nc
    off['eint'] = o; o += nc
    return dict(fn=fn, params=par, ptxt=ptxt, time=t, dt=dt, ncycle=ncyc, nmb=nmb, lloc=lloc, blocks=blocks,
                order=order, data_off=data_off, data_size=data_size, ng=ng, mesh=mi[1:4], mb=bi[1:4],
                rsz=rsz, n=(n1, n2, n3), nh=nh, trans=(f2, f3), off=off,
                counts=dict(npred=npred, nt2=nt2, neint=neint, nctr=nctr, nrss=nrss, nvgd=nvgd, ncsk=ncsk, nra=nra))


def gather(R, key, nvar, face1=False):
    """global array [nvar][NX3][NX2][n1(+1)]: lateral ACTIVE cells, radial incl. ghosts (one MeshBlock radially)."""
    ng = R['ng']
    n1, n2, n3 = R['n']
    NX1, NX2, NX3 = R['mesh']
    b1, b2, b3 = R['mb']
    if b1 != NX1:
        raise SystemExit('parent MeshBlocks must span the full radial range')
    m1 = n1 + 1 if face1 else n1
    out = np.empty((nvar, NX3, NX2, m1))
    mm = np.memmap(R['fn'], dtype='<f8', mode='r', offset=R['data_off'])
    rec = R['data_size'] // 8
    for m in range(R['nmb']):
        lx1, lx2, lx3, lev, _ = R['lloc'][m]
        base = m * rec + R['off'][key]
        a = np.asarray(mm[base:base + nvar * n3 * n2 * m1]).reshape(nvar, n3, n2, m1)
        out[:, lx3 * b3:(lx3 + 1) * b3, lx2 * b2:(lx2 + 1) * b2, :] = a[:, ng:ng + b3, ng:ng + b2, :]
    return out


def th_geom(f):
    w = np.cos(f[:-1]) - np.cos(f[1:])
    cen = ((np.sin(f[1:]) - f[1:] * np.cos(f[1:])) - (np.sin(f[:-1]) - f[:-1] * np.cos(f[:-1]))) / w
    return w, cen


def mc(dl, dr, dc):
    s = np.where(dl * dr > 0, 1.0, 0.0)
    return s * np.sign(dc) * np.minimum(np.minimum(2 * np.abs(dl), 2 * np.abs(dr)), np.abs(dc))


def prolong(q, thc_p, thc_c, phc_p, phc_c):
    """q [K+2, J+2, I] parent with a 1-cell lateral halo -> children [2K, 2J, I].
    thc_p [J+2] parent theta centroids, thc_c [2J] child centroids; phc likewise."""
    qc = q[1:-1, 1:-1]
    # theta slope
    tL = (q[1:-1, 1:-1] - q[1:-1, :-2]) / (thc_p[1:-1] - thc_p[:-2])[None, :, None]
    tR = (q[1:-1, 2:] - q[1:-1, 1:-1]) / (thc_p[2:] - thc_p[1:-1])[None, :, None]
    tC = (q[1:-1, 2:] - q[1:-1, :-2]) / (thc_p[2:] - thc_p[:-2])[None, :, None]
    st = mc(tL, tR, tC)
    pL = (q[1:-1, 1:-1] - q[:-2, 1:-1]) / (phc_p[1:-1] - phc_p[:-2])[:, None, None]
    pR = (q[2:, 1:-1] - q[1:-1, 1:-1]) / (phc_p[2:] - phc_p[1:-1])[:, None, None]
    pC = (q[2:, 1:-1] - q[:-2, 1:-1]) / (phc_p[2:] - phc_p[:-2])[:, None, None]
    sp = mc(pL, pR, pC)
    K, J, I = qc.shape
    dth = (thc_c.reshape(J, 2) - thc_p[1:-1, None])  # [J,2]
    dph = (phc_c.reshape(K, 2) - phc_p[1:-1, None])  # [K,2]
    # Barth-Jespersen clip on the 3x3 neighbourhood
    nb = np.stack([q[a:a + K, b:b + J] for a in range(3) for b in range(3)])
    qmax = nb.max(axis=0)
    qmin = nb.min(axis=0)
    alpha = np.ones_like(qc)
    for a in range(2):
        for b in range(2):
            d = st * dth[None, :, b, None] + sp * dph[:, a, None, None]
            with np.errstate(divide='ignore', invalid='ignore'):
                r = np.where(d > 0, (qmax - qc) / d, np.where(d < 0, (qmin - qc) / d, 1.0))
            alpha = np.minimum(alpha, np.clip(r, 0.0, 1.0))
    st = st * alpha
    sp = sp * alpha
    out = np.empty((2 * K, 2 * J, I))
    for a in range(2):
        for b in range(2):
            out[a::2, b::2] = qc + st * dth[None, :, b, None] + sp * dph[:, a, None, None]
    return out


def blend_axis(q, axis, W):
    """pairwise conservative seam exchange along axis (0 = phi, 1 = theta) of q [K, J, I] (in place)."""
    N = q.shape[axis]
    for d in range(W):
        a = 0.5 * (1.0 - d / W)
        iL, iR = N - 1 - d, d
        sl = [slice(None)] * 3
        sr = [slice(None)] * 3
        sl[axis] = iL
        sr[axis] = iR
        L = q[tuple(sl)].copy()
        Rr = q[tuple(sr)].copy()
        q[tuple(sl)] = (1 - a) * L + a * Rr
        q[tuple(sr)] = (1 - a) * Rr + a * L
    return q


def seam_metric(q, axis, isl):
    """rms seam jump / rms interior adjacent jump along axis, over radial slice isl."""
    qq = q[:, :, isl]
    N = qq.shape[axis]
    seam = np.take(qq, 0, axis=axis) - np.take(qq, N - 1, axis=axis)
    inner = np.diff(qq, axis=axis)
    return float(np.sqrt(np.mean(seam ** 2)) / np.sqrt(np.mean(inner ** 2))), float(np.sqrt(np.mean(seam ** 2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('rst')
    ap.add_argument('deck')
    ap.add_argument('outdir')
    ap.add_argument('--phi0', default='auto')
    ap.add_argument('--refine', type=int, default=2)
    ap.add_argument('--blend', type=int, default=8)
    ap.add_argument('--name', default='e1w_remap')
    ap.add_argument('--kefrac', type=float, default=0.01, help='cap of the sub-cell KE taken from eint (fraction of eint_p)')
    ap.add_argument('--rmin_metric', type=float, default=1.0e13, help='seam/phi0 metric radial range (cm)')
    ap.add_argument('--rmax_metric', type=float, default=3.5e13)
    a = ap.parse_args()
    if a.refine != 2:
        raise SystemExit('refine 2 only')
    os.makedirs(a.outdir, exist_ok=True)
    log = open(os.path.join(a.outdir, a.name + '.remap.txt'), 'w')

    def P(*s):
        line = ' '.join(str(x) for x in s)
        print(line)
        log.write(line + '\n')

    R = read_rst(a.rst)
    par = R['params']
    deck_txt = open(a.deck).read()
    dk = parse_params(deck_txt)
    P('parent rst', a.rst, 'size', os.path.getsize(a.rst))
    P('  t = %.17g  dt = %.17g  ncycle %d  nmb %d  mesh %s  meshblock %s  ng %d' %
      (R['time'], R['dt'], R['ncycle'], R['nmb'], R['mesh'], R['mb'], R['ng']))
    P('  blocks', R['order'], ' nhydro', R['nh'], ' f0x2/f0x3', R['trans'], ' counts', R['counts'])
    if 'problem/he_ic_rad_atm_na' in par:
        raise SystemExit('parent used he_ic_rad_atm_post (a_ref in the input): not supported')
    if par.get('problem/mlt_closure', 'frozen') != 'frozen':
        raise SystemExit('parent mlt_closure must be frozen')
    pg = R['blocks'].get(b'PGENST01', b'')
    if dk.get('problem/he_ic_equil', 'false').lower() in ('true', '1') and pg[-8:] != b'HSDICE01':
        raise SystemExit('child deck has he_ic_equil = true but the parent rst carries no he_ic_equil a_ref trailer')
    for k in ('problem/mlt_ramp_start', 'problem/mlt_ramp_time', 'problem/he_gm', 'problem/he_ic_file',
              'rad_m1/implicit_flux_x1min', 'problem/he_sponge_r0', 'hydro/wb_rmax'):
        if par.get(k) != dk.get(k):
            P('  NOTE key differs parent -> child deck: %s: %s -> %s' % (k, par.get(k), dk.get(k)))
    # radial grid identical
    rkeys = sorted(k for k in set(list(par) + list(dk)) if k.startswith('mesh/') and
                   (k.startswith('mesh/x1') or 'stretch_r' in k or k in ('mesh/nx1', 'mesh/nghost')))
    bad = []
    for k in rkeys:
        pv, cv = par.get(k), dk.get(k)
        dflt = '2' if k == 'mesh/nghost' else '0'   # GridKeyDefault (main.cpp): absent use_* = false, coeffs 0
        pv = {'true': '1', 'false': '0'}.get((pv or dflt).lower(), pv or dflt)
        cv = {'true': '1', 'false': '0'}.get((cv or dflt).lower(), cv or dflt)
        try:
            same = (pv is not None and cv is not None and float(pv) == float(cv))
        except ValueError:
            same = (pv == cv)
        if not same and not (pv is None and cv is None):
            bad.append((k, pv, cv))
    if bad:
        raise SystemExit('radial grid keys differ parent vs child deck: %s' % bad)
    if dk.get('mesh/use_grid_stretch_theta', 'false').lower() in ('true', '1') or \
       par.get('mesh/use_grid_stretch_theta', 'false').lower() in ('true', '1'):
        raise SystemExit('theta stretching not supported')
    NX1, NX2, NX3 = R['mesh']
    ng = R['ng']
    p2a, p2b = float(par['mesh/x2min']), float(par['mesh/x2max'])
    p3a, p3b = float(par['mesh/x3min']), float(par['mesh/x3max'])
    c2a, c2b = float(dk['mesh/x2min']), float(dk['mesh/x2max'])
    c3a, c3b = float(dk['mesh/x3min']), float(dk['mesh/x3max'])
    cn2, cn3 = int(dk['mesh/nx2']), int(dk['mesh/nx3'])
    if int(dk.get('mesh/nghost', '2')) != ng:
        raise SystemExit('nghost differs')
    dth_p = (p2b - p2a) / NX2
    dph_p = (p3b - p3a) / NX3
    j0f = (c2a - p2a) / dth_p
    nj = cn2 // 2
    nk = cn3 // 2
    j0 = int(round(j0f))
    if abs(j0f - j0) > 1e-6 or abs((c2b - c2a) / dth_p - nj) > 1e-6 or abs((c3b - c3a) / dph_p - nk) > 1e-6:
        raise SystemExit('child wedge is not an integer number of parent cells (or nx not 2x)')
    if abs(c2a + c2b - math.pi) > 1e-12:
        raise SystemExit('child theta range must be symmetric about pi/2 (equal-volume seam pairs)')
    per2 = abs((p2b - p2a) - math.pi / 3) < 1e-12 or par.get('mesh/ix2_bc') == 'periodic'
    P('  parent theta [%.15g, %.15g] %d, phi [%.15g, %.15g] %d; child theta [%.15g, %.15g] %d, phi [%.15g, %.15g] %d'
      % (p2a, p2b, NX2, p3a, p3b, NX3, c2a, c2b, cn2, c3a, c3b, cn3))
    P('  sub-wedge: parent theta cells %d..%d (%d), phi %d cells' % (j0, j0 + nj - 1, nj, nk))

    # ---- parent arrays
    hyd = gather(R, 'hyd', R['nh'])[:5]
    m1 = gather(R, 'm1', 4)
    f01 = gather(R, 'f0x1', 1, face1=True)[0]
    eint = gather(R, 'eint', 1)[0]
    rho = hyd[0]
    if R['nmb'] > 1:   # MeshBlock placement check: jumps across MeshBlock borders vs elsewhere
        b2, b3 = R['mb'][1], R['mb'][2]
        isl0 = slice(R['ng'], R['ng'] + R['mesh'][0])
        for nm, q in (('rho', rho), ('eint', None), ('E', None)):
            q = rho if nm == 'rho' else (gather(R, 'eint', 1)[0] if nm == 'eint' else gather(R, 'm1', 1)[0])
            q = q[:, :, isl0]
            out = []
            for ax, b in ((1, b2), (0, b3)):
                dq = np.abs(np.diff(q, axis=ax))
                N = dq.shape[ax]
                bord = [n for n in range(N) if (n + 1) % b == 0]
                oth = [n for n in range(N) if (n + 1) % b != 0]
                out.append(float(np.sqrt(np.mean(np.take(dq, bord, axis=ax) ** 2)) /
                                 np.sqrt(np.mean(np.take(dq, oth, axis=ax) ** 2))) if bord else float('nan'))
            P('  MeshBlock-border check %-4s rms jump across borders / elsewhere: theta %.3f phi %.3f' % (nm, out[0], out[1]))
    ke = 0.5 * (hyd[1] ** 2 + hyd[2] ** 2 + hyd[3] ** 2) / rho
    rphi = hyd[4] - ke - eint
    phi = rphi / rho
    n1t = rho.shape[2]
    isl = slice(ng, ng + NX1)
    ps = phi[:, :, isl]
    spread = np.max(np.abs(ps - ps.mean(axis=(0, 1))[None, None]) / np.maximum(np.abs(ps.mean(axis=(0, 1))), 1e-300))
    P('  rho*Phi from the parent (ET - KE - eint): Phi lateral spread (max rel, active) %.3e; min eint %.4e, min rho %.4e'
      % (spread, eint.min(), rho.min()))

    # parent geometry (with wrap for the halo)
    thf_p = p2a + dth_p * np.arange(NX2 + 1)
    wth_p, thc_p = th_geom(thf_p)
    phc_p = p3a + dph_p * (np.arange(NX3) + 0.5)

    # ---- phi offset
    jsl = slice(j0, j0 + nj)
    # metric radial range: cells whose... use the parent x1 centres are unknown; use index range of the radial window by
    # the IC's cell radii is not available -> metric over all active radial cells weighted equally (rmin/rmax unused)
    def mism(k0):
        kL = (k0 + nk - 1) % NX3
        s = 0.0
        for q in (rho, eint, hyd[1] / rho, m1[0]):
            qs = q[:, jsl, isl]
            d = (qs[k0] - qs[kL])
            typ = np.diff(qs[:, :, :], axis=0)
            s += np.sqrt(np.mean(d ** 2)) / np.sqrt(np.mean(typ ** 2))
        return s / 4.0
    if a.phi0 == 'auto':
        scores = [mism(k) for k in range(NX3)]
        k0 = int(np.argmin(scores))
        P('  phi0 auto: parent phi cell offset %d (seam mismatch / mean adjacent jump %.3f; k0=0: %.3f; median %.3f)'
          % (k0, scores[k0], scores[0], float(np.median(scores))))
    else:
        k0 = int(a.phi0)
        P('  phi0 given: parent phi cell offset %d (mismatch %.3f)' % (k0, mism(k0)))
    phi_shift = p3a + k0 * dph_p - c3a   # child phi = parent phi - phi_shift

    kidx = np.arange(k0 - 1, k0 + nk + 1) % NX3
    jidx = np.arange(j0 - 1, j0 + nj + 1)
    if jidx[0] < 0 or jidx[-1] >= NX2:
        if not per2:
            raise SystemExit('sub-wedge touches the parent theta edge')
        jidx = jidx % NX2

    def sub(q):
        return q[np.ix_(kidx, jidx)] if q.ndim == 3 else None

    # halo centroids (unwrapped coordinates)
    thf_h = p2a + dth_p * np.arange(j0 - 1, j0 + nj + 2)
    _, thc_h = th_geom(thf_h)
    phc_h = p3a + dph_p * (np.arange(k0 - 1, k0 + nk + 1) + 0.5)
    thf_c = c2a + (c2b - c2a) / cn2 * np.arange(cn2 + 1)
    wth_c, thc_c = th_geom(thf_c)
    thc_c_par = thc_c  # child theta coincides with parent theta (no shift)
    dph_c = (c3b - c3a) / cn3
    phc_c = c3a + dph_c * (np.arange(cn3) + 0.5) + phi_shift
    wth_s = wth_p[j0:j0 + nj]

    # ---- prolong
    names = ['rho', 'eint', 'E', 'F1', 'F2', 'F3', 'f0x1']
    src = [rho, eint, m1[0], m1[1], m1[2], m1[3], f01]
    ch = {}
    for nm, q in zip(names, src):
        ch[nm] = prolong(q[np.ix_(kidx, jidx)], thc_h, thc_c_par, phc_h, phc_c)
    rphi_s = rphi[np.ix_(kidx[1:-1], jidx[1:-1])]
    phi_s = rphi_s / rho[np.ix_(kidx[1:-1], jidx[1:-1])]
    ke_s = ke[np.ix_(kidx[1:-1], jidx[1:-1])]
    e_s = eint[np.ix_(kidx[1:-1], jidx[1:-1])]

    def up(x):  # parent [K,J,I] -> children shape
        return np.repeat(np.repeat(x, 2, axis=0), 2, axis=1)
    rc = ch['rho']
    # volume weights (lateral only; radial factor common within a radial index)
    vc = (wth_c[None, :] * dph_c)[:, :, None] * np.ones((cn3, 1, 1))
    vp = (wth_s[None, :] * dph_p)[:, :, None]

    def coarsen(x):  # children [2K,2J,I] (volume-weighted sum) -> parent [K,J,I] sums
        y = x * vc
        return y[0::2, 0::2] + y[1::2, 0::2] + y[0::2, 1::2] + y[1::2, 1::2]
    # velocities: v_c = vbar_p + alpha_p dv_c, vbar_p = m_p/rho_p (mass-weighted mean), dv_c = the limited-linear velocity
    # deviation with its mass-weighted mean removed (momentum exact for any alpha); the KE this adds (alpha^2 dKE1, >= 0)
    # is capped at kefrac x eint_p per parent cell and taken from the children's eint (ET exact per parent cell)
    mass_c = coarsen(rc)
    dke1 = np.zeros_like(rc[0::2, 0::2])
    dv = {}
    for n in (1, 2, 3):
        vb = up(hyd[n][np.ix_(kidx[1:-1], jidx[1:-1])] / rho[np.ix_(kidx[1:-1], jidx[1:-1])])
        vl = prolong((hyd[n] / rho)[np.ix_(kidx, jidx)], thc_h, thc_c_par, phc_h, phc_c) - vb
        vl = vl - up(coarsen(rc * vl) / mass_c)
        dv[n] = (vb, vl)
        dke1 += coarsen(0.5 * rc * vl ** 2)
    dke1 = dke1 / vp
    cap = a.kefrac * e_s
    with np.errstate(divide='ignore', invalid='ignore'):
        alpha = np.where(dke1 > cap, np.sqrt(cap / dke1), 1.0)
    for n in (1, 2, 3):
        vb, vl = dv[n]
        ch['m%d' % n] = rc * (vb + up(alpha) * vl)
    kec = 0.5 * (ch['m1'] ** 2 + ch['m2'] ** 2 + ch['m3'] ** 2) / rc
    dke = (ke_s * vp - coarsen(kec)) / vp  # = -alpha^2 dKE1 (<= 0) up to round-off
    ch['eint'] = ch['eint'] + up(dke)
    P('  velocity slopes: KE added by the sub-cell velocity structure taken from eint, capped at %.3g eint_p: max %.3e of eint_p;'
      ' alpha < 1 in %d of %d parent cells (min alpha %.3f)' %
      (a.kefrac, float(np.max(np.abs(dke) / e_s)), int(np.sum(alpha < 1)), alpha.size, float(alpha.min())))
    for n in (1, 2, 3):
        mp = hyd[n][np.ix_(kidx[1:-1], jidx[1:-1])] * vp
        P('  momentum m%d per parent cell: max |sum children - parent| / max|m| %.3e' %
          (n, float(np.max(np.abs(coarsen(ch['m%d' % n]) - mp)) / np.max(np.abs(mp)))))
    if ch['eint'].min() <= 0 or rc.min() <= 0 or ch['E'].min() <= 0:
        raise SystemExit('non-positive rho/eint/E after prolongation')
    ch['ET'] = ch['eint'] + kec + rc * up(phi_s)

    # conservation check of the prolongation (per parent cell and per shell)
    def relerr_shell(xc, xp):
        sc = coarsen(xc).sum(axis=(0, 1))
        sp_ = (xp * vp).sum(axis=(0, 1))
        return float(np.max(np.abs(sc - sp_) / np.maximum(np.abs(sp_), 1e-300)))
    psub = {'rho': rho, 'ET': hyd[4], 'E': m1[0], 'eint': eint, 'm1': hyd[1], 'F1': m1[1], 'f0x1': f01}
    for nm in ('rho', 'ET', 'E', 'eint', 'f0x1'):
        P('  prolongation conservation per shell (max rel over all radial cells) %-5s %.3e' %
          (nm, relerr_shell(ch[nm], psub[nm][np.ix_(kidx[1:-1], jidx[1:-1])])))
    # seams before
    msl = slice(ng, ng + NX1)
    pre = {}
    for nm in ('rho', 'ET', 'E', 'm1'):
        pre[nm] = (seam_metric(ch[nm], 0, msl)[0], seam_metric(ch[nm], 1, msl)[0])
    vr0 = ch['m1'] / ch['rho']
    pre['vr'] = (seam_metric(vr0, 0, msl)[0], seam_metric(vr0, 1, msl)[0])
    # ---- blend seams (conservative)
    W = a.blend
    shell0 = {nm: (ch[nm] * vc).sum(axis=(0, 1)) for nm in ('rho', 'ET', 'E', 'm1', 'f0x1')}
    for nm in ('rho', 'm1', 'm2', 'm3', 'ET', 'E', 'F1', 'F2', 'F3', 'f0x1'):
        blend_axis(ch[nm], 0, W)
        blend_axis(ch[nm], 1, W)
    for nm in ('rho', 'ET', 'E', 'm1', 'f0x1'):
        s1 = (ch[nm] * vc).sum(axis=(0, 1))
        P('  blend conservation per shell %-5s max rel %.3e' %
          (nm, float(np.max(np.abs(s1 - shell0[nm]) / np.maximum(np.abs(shell0[nm]), 1e-300)))))
    kec = 0.5 * (ch['m1'] ** 2 + ch['m2'] ** 2 + ch['m3'] ** 2) / ch['rho']
    eic = ch['ET'] - kec - ch['rho'] * up(phi_s)
    P('  after blend: min eint %.4e (min eint/parent-min %.3f), min rho %.4e, min E %.4e' %
      (eic.min(), eic.min() / eint.min(), ch['rho'].min(), ch['E'].min()))
    if eic.min() <= 0:
        raise SystemExit('non-positive eint after the blend')
    vr1 = ch['m1'] / ch['rho']
    P('  seam metric = rms(seam jump)/rms(adjacent jump), active radial cells: [phi seam, theta seam] before -> after (W=%d)' % W)
    for nm in ('rho', 'ET', 'E', 'm1'):
        post = (seam_metric(ch[nm], 0, msl)[0], seam_metric(ch[nm], 1, msl)[0])
        P('    %-4s [%.3f, %.3f] -> [%.2e, %.2e]' % (nm, pre[nm][0], pre[nm][1], post[0], post[1]))
    post = (seam_metric(vr1, 0, msl)[0], seam_metric(vr1, 1, msl)[0])
    P('    v_r  [%.3f, %.3f] -> [%.2e, %.2e]' % (pre['vr'][0], pre['vr'][1], post[0], post[1]))
    # blend-zone gradient: max adjacent jump inside the blend band vs interior (rms)
    for nm, q in (('rho', ch['rho']), ('v_r', vr1)):
        qq = q[:, :, msl]
        dphi_ = np.abs(np.diff(qq, axis=0))
        band = np.concatenate([dphi_[:W], dphi_[-W:]])
        P('    %s phi-adjacent jumps: rms in blend band %.3e, rms interior %.3e' %
          (nm, float(np.sqrt(np.mean(band ** 2))), float(np.sqrt(np.mean(dphi_[W:-W] ** 2)))))
    # totals (wedge, active cells), as a parent-sub-wedge check
    P('  child wedge totals (lateral-weighted sums, radial factor omitted): rho %.10e ET %.10e E %.10e' %
      tuple(float((ch[nm][:, :, msl] * vc).sum()) for nm in ('rho', 'ET', 'E')))

    # ---- write (with lateral periodic ghosts)
    def padlat(x):
        return np.pad(x, ((ng, ng), (ng, ng), (0, 0)), mode='wrap')
    t = R['time']
    fn = os.path.join(a.outdir, a.name + '.remap')
    hdr = b'HSREMAP1' + struct.pack('<8i', 1, ng, NX1, cn2, cn3, 5, 4, 0) + \
        struct.pack('<7d', t, float(dk['mesh/x1min']), float(dk['mesh/x1max']), c2a, c2b, c3a, c3b)
    with open(fn, 'wb') as f:
        f.write(hdr)
        f.write(struct.pack('<Q', len(pg)))
        f.write(pg)
        for nm in ('rho', 'm1', 'm2', 'm3', 'ET'):
            f.write(np.ascontiguousarray(padlat(ch[nm]), dtype='<f8').tobytes())
        for nm in ('E', 'F1', 'F2', 'F3'):
            f.write(np.ascontiguousarray(padlat(ch[nm]), dtype='<f8').tobytes())
        f.write(np.ascontiguousarray(padlat(ch['f0x1']), dtype='<f8').tobytes())
    md5 = hashlib.md5(open(fn, 'rb').read()).hexdigest()
    P('  wrote %s (%d bytes) md5 %s; pgen block %d bytes' % (fn, os.path.getsize(fn), md5, len(pg)))

    # ---- deck
    def setkey(txt, blk, key, val, cmt):
        lines = txt.split('\n')
        bi = None
        for n, l in enumerate(lines):
            if re.match(r'\s*<\s*%s\s*>' % re.escape(blk), l):
                bi = n
                continue
            if bi is not None and re.match(r'\s*<', l):
                break
            if bi is not None and re.match(r'\s*%s\s*=' % re.escape(key), l):
                lines[n] = '%-24s = %s   # %s' % (key, val, cmt)
                return '\n'.join(lines)
        if bi is None:
            raise SystemExit('block <%s> not in deck' % blk)
        lines.insert(bi + 1, '%-24s = %s   # %s' % (key, val, cmt))
        return '\n'.join(lines)
    d = deck_txt
    d = setkey(d, 'time', 'start_time', repr(t), 'e1remap-1014: = parent rst time (checked bitwise by the pgen)')
    d = setkey(d, 'problem', 'he_ic_remap_file', os.path.abspath(fn), 'e1remap-1014: md5 ' + md5)
    for blk in sorted(b for b in set(k.split('/')[0] for k in dk) if b.startswith('output')):
        odt = float(dk[blk + '/dt'])
        nlast = math.floor(t / odt + 1e-9)
        g = (nlast + 1) * odt
        d = setkey(d, blk, 'last_time', repr(g - 2 * odt), 'e1remap-1014: next scheduled output at %.6e (parent grid)' % g)
        d = setkey(d, blk, 'file_number', str(nlast), 'e1remap-1014: parent numbering (dump n at ~n x dt)')
    dfn = os.path.join(a.outdir, a.name + '.athinput')
    open(dfn, 'w').write(d)
    P('  wrote deck %s md5 %s' % (dfn, hashlib.md5(d.encode()).hexdigest()))
    log.close()


if __name__ == '__main__':
    main()
