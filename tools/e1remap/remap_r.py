#!/usr/bin/env python3
"""remap_r.py (s3grid-1017): AthenaK restart (he_star_m1 + rad_m1, general EOS, spherical polar, ONE MeshBlock radially)
-> the same state on a NEW RADIAL GRID (same x1min/x1max, same lateral grid) as a problem/he_ic_remap_file start
(pgen he_star_m1, e1remap-1014 path; file format of remap_e1w.py).

usage: python3 remap_r.py PARENT_RST PARENT_X1GRID CHILD_TEMPLATE_RST CHILD_X1GRID CHILD_DECK OUTDIR [--name NAME]
  PARENT_RST           the state to remap (e.g. S3 rst 00046, t = 3.0e7)
  PARENT_X1GRID        its <basename>.x1grid.txt (faces incl. ghosts, written at startup)
  CHILD_TEMPLATE_RST   any restart of a FRESH run of CHILD_DECK (new grid): supplies the child's Phi_eff per cell
                       (u0(IEN) - KE - eint; the new grid's he_ic_equil a_ref and Phi_eff) and its PGEN block (the a_ref
                       trailer of the new grid; the Picard window trailer is stripped by the pgen)
  CHILD_X1GRID         the template run's x1grid.txt
  CHILD_DECK           the child deck (new grid keys); the generated deck = it + start_time + he_ic_remap_file + output
                       schedules aligned to the parent
Method (per lateral column, radial index incl. x1 ghosts):
  cell quantities rho, m1, m2, m3, eint, E, F1, F2, F3: conservative remap of the parent's piecewise-LINEAR profile
  (slopes at the volume centroids x1v, MC-limited: monotone, positive) integrated exactly over each child shell
  (int q r^2 dr); child cells beyond the parent's ghost range take the edge parent value (ghosts, refilled by the BCs).
  ET_child = eint + |m|^2/(2 rho) + rho Phi_child (Phi_child from the template: the code's own Phi_eff of the new grid).
  f0x1 (x1 face fluxes): linear interpolation in r between the parent faces.  Lateral ghosts: periodic images.
Conservation (mass, momentum, eint, E per column) is exact up to round-off where both grids cover the range; the total
energy changes only by rho (Phi_child - Phi_parent) (the two grids' discrete equilibria), reported.
"""
import argparse
import hashlib
import math
import os
import re
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from remap_e1w import read_rst, gather, parse_params  # noqa: E402


def read_x1grid(fn):
    rows = [ln.split() for ln in open(fn) if not ln.startswith('#') and ln.strip()]
    four = [r for r in rows if len(r) == 4]
    two = [r for r in rows if len(r) == 2]
    xf = np.array([float(r[1]) for r in four] + [float(two[-1][1])])
    xv = np.array([float(r[2]) for r in four])
    return xf, xv


def mc_slopes(q, xv):
    """MC-limited slopes along the last axis (cells), one-sided zero at the ends"""
    dl = np.diff(q, axis=-1) / np.diff(xv)
    s = np.zeros_like(q)
    a, b = dl[..., :-1], dl[..., 1:]
    c = (q[..., 2:] - q[..., :-2]) / (xv[2:] - xv[:-2])
    sg = np.where(a * b > 0, np.sign(c), 0.0)
    s[..., 1:-1] = sg * np.minimum(np.minimum(2 * np.abs(a), 2 * np.abs(b)), np.abs(c))
    return s


def remap_cells(q, pf, pv, cf):
    """q [..., Np] parent cell averages on faces pf [Np+1] (centroids pv) -> child averages on faces cf [Nc+1]"""
    s = mc_slopes(q, pv)
    Nc = len(cf) - 1
    out = np.zeros(q.shape[:-1] + (Nc,))
    for n in range(Nc):
        a, b = cf[n], cf[n + 1]
        lo = max(a, pf[0]); hi = min(b, pf[-1])
        tot = np.zeros(q.shape[:-1])
        if hi > lo:
            i0 = int(np.searchsorted(pf, lo, side='right')) - 1
            i1 = int(np.searchsorted(pf, hi, side='left'))
            for i in range(max(i0, 0), min(i1, len(pf) - 1)):
                x0, x1 = max(lo, pf[i]), min(hi, pf[i + 1])
                if x1 <= x0:
                    continue
                v3 = (x1 ** 3 - x0 ** 3) / 3.0
                v4 = (x1 ** 4 - x0 ** 4) / 4.0
                tot += q[..., i] * v3 + s[..., i] * (v4 - pv[i] * v3)
        # outside the parent range: edge values
        if a < pf[0]:
            tot += q[..., 0] * (min(b, pf[0]) ** 3 - a ** 3) / 3.0
        if b > pf[-1]:
            tot += q[..., -1] * (b ** 3 - max(a, pf[-1]) ** 3) / 3.0
        out[..., n] = tot / ((b ** 3 - a ** 3) / 3.0)
    return out


def main():
    ap = argparse.ArgumentParser()
    for k in ('prst', 'pgrid', 'crst', 'cgrid', 'deck', 'outdir'):
        ap.add_argument(k)
    ap.add_argument('--name', default='s3grid_remap')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    log = open(os.path.join(a.outdir, a.name + '.remap.txt'), 'w')

    def P(*s):
        line = ' '.join(str(x) for x in s)
        print(line)
        log.write(line + '\n')

    R = read_rst(a.prst)
    C = read_rst(a.crst)
    dk_txt = open(a.deck).read()
    dk = parse_params(dk_txt)
    P('parent', a.prst, 't = %.17g ncycle %d mesh %s mb %s ng %d nmb %d' % (R['time'], R['ncycle'], R['mesh'], R['mb'], R['ng'],
                                                                       R['nmb']))
    P('child template', a.crst, 't = %.17g ncycle %d mesh %s mb %s ng %d nmb %d' % (C['time'], C['ncycle'], C['mesh'], C['mb'],
                                                                               C['ng'], C['nmb']))
    if R['ng'] != C['ng'] or R['mesh'][1:] != C['mesh'][1:]:
        raise SystemExit('lateral grid / nghost differ')
    for k in ('mesh/x1min', 'mesh/x1max', 'mesh/x2min', 'mesh/x2max', 'mesh/x3min', 'mesh/x3max'):
        if float(R['params'][k]) != float(C['params'][k]) or float(dk[k]) != float(C['params'][k]):
            raise SystemExit('extent %s differs' % k)
    if int(dk['mesh/nx1']) != C['mesh'][0]:
        raise SystemExit('deck nx1 != template nx1')
    pf, pv = read_x1grid(a.pgrid)
    cf, cv = read_x1grid(a.cgrid)
    ng = R['ng']
    n1p, n1c = R['n'][0], C['n'][0]
    if len(pv) != n1p or len(cv) != n1c:
        raise SystemExit('x1grid sizes %d/%d vs rst %d/%d' % (len(pv), len(cv), n1p, n1c))
    if abs(pf[ng] / float(R['params']['mesh/x1min']) - 1) > 1e-12 or abs(cf[ng] / float(dk['mesh/x1min']) - 1) > 1e-12:
        raise SystemExit('x1grid files do not match x1min')
    P('  parent faces incl. ghosts %.6e .. %.6e (%d cells); child %.6e .. %.6e (%d cells)' % (pf[0], pf[-1], n1p, cf[0], cf[-1],
                                                                                            n1c))
    hyd = gather(R, 'hyd', R['nh'])[:5]
    m1 = gather(R, 'm1', 4)
    f01 = gather(R, 'f0x1', 1, face1=True)[0]
    eint = gather(R, 'eint', 1)[0]
    rho = hyd[0]
    ke = 0.5 * (hyd[1] ** 2 + hyd[2] ** 2 + hyd[3] ** 2) / rho
    phi_p = ((hyd[4] - ke - eint) / rho).mean(axis=(0, 1))
    # child Phi_eff from the template
    ch_h = gather(C, 'hyd', C['nh'])[:5]
    ch_e = gather(C, 'eint', 1)[0]
    ch_ke = 0.5 * (ch_h[1] ** 2 + ch_h[2] ** 2 + ch_h[3] ** 2) / ch_h[0]
    phc3 = (ch_h[4] - ch_ke - ch_e) / ch_h[0]
    phi_c = phc3.mean(axis=(0, 1))
    spr = np.max(np.abs(phc3 - phi_c) / np.maximum(np.abs(phi_c), 1e-300))
    P('  child Phi_eff from the template: lateral spread (max rel) %.3e; Phi range %.4e .. %.4e' % (spr, phi_c.min(), phi_c.max()))
    # parent Phi interpolated to child centroids, for the energy report
    phi_pc = np.interp(cv, pv, phi_p)
    out = {}
    names = ['rho', 'm1', 'm2', 'm3', 'eint', 'E', 'F1', 'F2', 'F3']
    src = [rho, hyd[1], hyd[2], hyd[3], eint, m1[0], m1[1], m1[2], m1[3]]
    for nm, q in zip(names, src):
        out[nm] = remap_cells(q, pf, pv, cf)
    out['f0x1'] = np.stack([np.stack([np.interp(cf, pf, f01[k, j]) for j in range(f01.shape[1])]) for k in range(f01.shape[0])])
    # conservation per column over the common active range [x1min, x1max]
    Vp = (pf[1:] ** 3 - pf[:-1] ** 3) / 3.0
    Vc = (cf[1:] ** 3 - cf[:-1] ** 3) / 3.0
    ap_ = slice(ng, ng + R['mesh'][0]); ac_ = slice(ng, ng + C['mesh'][0])
    for nm, q in zip(names, src):
        sp = (q[..., ap_] * Vp[ap_]).sum(axis=-1)
        sc = (out[nm][..., ac_] * Vc[ac_]).sum(axis=-1)
        P('  conservation %-4s per column over [x1min, x1max]: max |child-parent|/max|parent| %.3e' %
          (nm, float(np.max(np.abs(sc - sp)) / max(float(np.max(np.abs(sp))), 1e-300))))
    if out['rho'].min() <= 0 or out['eint'].min() <= 0 or out['E'].min() <= 0:
        raise SystemExit('non-positive rho/eint/E after the remap')
    kec = 0.5 * (out['m1'] ** 2 + out['m2'] ** 2 + out['m3'] ** 2) / out['rho']
    out['ET'] = out['eint'] + kec + out['rho'] * phi_c
    dEpot = float(((out['rho'] * (phi_c - phi_pc))[..., ac_] * Vc[ac_]).sum())
    Etot_p = float((hyd[4][..., ap_] * Vp[ap_]).sum())
    Etot_c = float((out['ET'][..., ac_] * Vc[ac_]).sum())
    P('  total energy (lateral factor omitted): parent %.10e child %.10e (rel %.3e); of which rho (Phi_c - Phi_p) %.3e (rel %.3e)'
      % (Etot_p, Etot_c, Etot_c / Etot_p - 1, dEpot, dEpot / Etot_p))
    P('  min rho %.3e eint %.3e E %.3e; max |v_r| %.3e km/s' % (out['rho'].min(), out['eint'].min(), out['E'].min(),
                                                               np.abs(out['m1'] / out['rho']).max() / 1e5))
    # shell-mean density at a few radii, parent vs child (sanity)
    for x in (0.9, 1.0, 1.05, 1.1, 1.12, 1.15, 1.2):
        r = x * 2.7025e13
        ip = int(np.argmin(abs(pv - r))); ic = int(np.argmin(abs(cv - r)))
        P('    r/Rph %.2f: <rho> parent %.4e (r %.4e) child %.4e (r %.4e)' % (x, rho[..., ip].mean(), pv[ip], out['rho'][..., ic].mean(),
                                                                          cv[ic]))

    def padlat(x):
        return np.pad(x, ((ng, ng), (ng, ng), (0, 0)), mode='wrap')
    pg = C['blocks'].get(b'PGENST01', b'')
    t = R['time']
    NX1, NX2, NX3 = C['mesh']
    fn = os.path.join(a.outdir, a.name + '.remap')
    hdr = b'HSREMAP1' + struct.pack('<8i', 1, ng, NX1, NX2, NX3, 5, 4, 0) + \
        struct.pack('<7d', t, float(dk['mesh/x1min']), float(dk['mesh/x1max']), float(dk['mesh/x2min']),
                    float(dk['mesh/x2max']), float(dk['mesh/x3min']), float(dk['mesh/x3max']))
    with open(fn, 'wb') as f:
        f.write(hdr)
        f.write(struct.pack('<Q', len(pg)))
        f.write(pg)
        for nm in ('rho', 'm1', 'm2', 'm3', 'ET', 'E', 'F1', 'F2', 'F3'):
            f.write(np.ascontiguousarray(padlat(out[nm]), dtype='<f8').tobytes())
        f.write(np.ascontiguousarray(padlat(out['f0x1']), dtype='<f8').tobytes())
    md5 = hashlib.md5(open(fn, 'rb').read()).hexdigest()
    P('  wrote %s (%d bytes) md5 %s; pgen block (child template) %d bytes, tail %s' % (fn, os.path.getsize(fn), md5, len(pg),
                                                                                    pg[-8:]))

    def setkey(txt, blk, key, val, cmt):
        lines = txt.split('\n')
        bi = None
        for n, ln in enumerate(lines):
            if re.match(r'\s*<\s*%s\s*>' % re.escape(blk), ln):
                bi = n
                continue
            if bi is not None and re.match(r'\s*<', ln):
                break
            if bi is not None and re.match(r'\s*%s\s*=' % re.escape(key), ln):
                lines[n] = '%-24s = %s   # %s' % (key, val, cmt)
                return '\n'.join(lines)
        if bi is None:
            raise SystemExit('block <%s> not in deck' % blk)
        lines.insert(bi + 1, '%-24s = %s   # %s' % (key, val, cmt))
        return '\n'.join(lines)
    d = dk_txt
    d = setkey(d, 'time', 'start_time', repr(t), 's3grid-1017: = parent rst time (checked bitwise by the pgen)')
    d = setkey(d, 'problem', 'he_ic_remap_file', os.path.abspath(fn), 's3grid-1017 remap_r.py: md5 ' + md5)
    for blk in sorted(b for b in set(k.split('/')[0] for k in dk) if b.startswith('output')):
        odt = float(dk[blk + '/dt'])
        nlast = math.floor(t / odt + 1e-9)
        g = (nlast + 1) * odt
        d = setkey(d, blk, 'last_time', repr(g - 2 * odt), 's3grid-1017: next scheduled output at %.6e' % g)
        d = setkey(d, blk, 'file_number', str(nlast), 's3grid-1017: parent numbering (dump n at ~n x dt)')
    dfn = os.path.join(a.outdir, a.name + '.athinput')
    open(dfn, 'w').write(d)
    P('  wrote deck %s md5 %s' % (dfn, hashlib.md5(d.encode()).hexdigest()))
    log.close()


if __name__ == '__main__':
    main()
