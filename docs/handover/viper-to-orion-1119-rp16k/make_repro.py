#!/usr/bin/env python3
"""make_repro.py RST PANEL OUTNAME [--n 8] [--vt keep|zero]: reduced RSG repro deck.
A narrow theta/phi-periodic spherical-polar wedge (n x n cells of the cs cell size (pi/2)/128 about the
equator) with the parent's radial grid, filled with the n x n cells at the centre of cs panel PANEL of the
restart (per radial index incl. ghosts): rho, radial momentum (cs IM1 = x1 = r), E, F_r = the comoving
radial flux (f0x1 at the cell centre) + 4/3 E v_r (lab frame), f0x1; tangential momenta and flux = 0 (--vt zero) and ET = u0(IEN) minus the parent's
tangential KE.  Writes OUTNAME.remap (HSREMAP1, pgen he_star_m1 problem/he_ic_remap_file) and
OUTNAME.athinput (parent deck from the restart header with the wedge mesh and start keys)."""
import argparse, hashlib, math, os, re, struct, sys
sys.path.insert(0, '/viper/ptmp2/jinma/ccint_1118/wt/tools/e1remap')
import numpy as np
import remap_e1w as RE
ap = argparse.ArgumentParser(); ap.add_argument('rst'); ap.add_argument('panel', type=int); ap.add_argument('out')
ap.add_argument('--n', type=int, default=8); ap.add_argument('--vt', default='zero'); a = ap.parse_args()
R = RE.read_rst(a.rst); ng = R['ng']; n1, n2, n3 = R['n']; NX1, NX2, NX3 = R['mesh']; b1, b2, b3 = R['mb']
mm = np.memmap(a.rst, dtype='<f8', mode='r', offset=R['data_off']); rec = R['data_size'] // 8
def panel(key, nvar, m1n):
    out = np.empty((nvar, NX3, NX2, m1n))
    for m in range(R['nmb']):
        lx1, lx2, lx3, lev, p = R['lloc'][m]
        if p != a.panel: continue
        base = m * rec + R['off'][key]
        q = np.asarray(mm[base:base + nvar * n3 * n2 * m1n]).reshape(nvar, n3, n2, m1n)
        out[:, lx3 * b3:(lx3 + 1) * b3, lx2 * b2:(lx2 + 1) * b2, :] = q[:, ng:ng + b3, ng:ng + b2, :]
    return out
c0 = NX2 // 2 - a.n // 2; sl = slice(c0, c0 + a.n)
hyd = panel('hyd', 5, n1)[:, sl, sl, :]; m1 = panel('m1', 4, n1)[:, sl, sl, :]
f0 = panel('f0x1', 1, n1 + 1)[0, sl, sl, :]
rho = hyd[0]; ax = int(np.argmax([abs(np.mean(m1[1 + c][:, :, ng:ng + 100])) for c in range(3)]))
sg = math.copysign(1.0, np.mean(m1[1 + ax][:, :, ng:ng + 100]))
ket = (hyd[2] ** 2 + hyd[3] ** 2) / (2 * rho)
fc = 0.5 * (f0[..., :-1] + f0[..., 1:])   # comoving radial flux at the cell centre
keep = (a.vt == 'keep')   # keep: the cs tangential momentum components as the wedge's theta, phi ones
ch = dict(rho=rho, m1=hyd[1], m2=hyd[2] if keep else 0 * rho, m3=hyd[3] if keep else 0 * rho,
          ET=hyd[4] if keep else hyd[4] - ket, E=m1[0],
          F1=fc + 4.0 / 3.0 * m1[0] * hyd[1] / rho,
          F2=0 * rho, F3=0 * rho, f0x1=f0)
print('panel %d radial F axis %d sign %+d; |F_t|/|F| inner %.2e; KE_t/ET inner %.2e' % (a.panel, ax, sg,
      float(np.mean(np.sqrt(sum(m1[1 + c] ** 2 for c in range(3) if c != ax))[:, :, ng:ng + 100]) /
            np.mean(np.abs(m1[1 + ax][:, :, ng:ng + 100]))), float(np.mean(ket[:, :, ng:ng + 100] / hyd[4][:, :, ng:ng + 100]))))
d = math.pi / 2 / NX2; c2a, c2b = math.pi / 2 - a.n / 2 * d, math.pi / 2 + a.n / 2 * d; c3a, c3b = 0.0, a.n * d
pad = lambda x: np.pad(x, ((ng, ng), (ng, ng), (0, 0)), mode='wrap')
txt = open(a.rst, 'rb').read(4 << 20); txt = txt[:txt.find(b'<par_end>')].decode(errors='ignore')
par = RE.parse_params(txt); pg = R['blocks'].get(b'PGENST01', b''); t = R['time']
fn = os.path.abspath(a.out + '.remap')
with open(fn, 'wb') as f:
    f.write(b'HSREMAP1' + struct.pack('<8i', 1, ng, NX1, a.n, a.n, 5, 4, 0) +
            struct.pack('<7d', t, float(par['mesh/x1min']), float(par['mesh/x1max']), c2a, c2b, c3a, c3b))
    f.write(struct.pack('<Q', len(pg))); f.write(pg)
    for nm in ('rho', 'm1', 'm2', 'm3', 'ET', 'E', 'F1', 'F2', 'F3', 'f0x1'):
        f.write(np.ascontiguousarray(pad(ch[nm]), dtype='<f8').tobytes())
md5 = hashlib.md5(open(fn, 'rb').read()).hexdigest()
def setkey(s, blk, key, val):
    L = s.split('\n'); bi = None
    for n, l in enumerate(L):
        if re.match(r'\s*<\s*%s\s*>' % re.escape(blk), l): bi = n; continue
        if bi is not None and re.match(r'\s*<', l): break
        if bi is not None and re.match(r'\s*%s\s*=' % re.escape(key), l):
            L[n] = '%-24s = %s   # make_repro' % (key, val); return '\n'.join(L)
    if bi is None: L += ['<%s>' % blk, '%s = %s' % (key, val)]
    else: L.insert(bi + 1, '%-24s = %s   # make_repro' % (key, val))
    return '\n'.join(L)
s = txt
for k, v in [('use_cubed_sphere', 'false'), ('use_spherical_polar', 'true'), ('nx2', a.n), ('nx3', a.n),
             ('x2min', repr(c2a)), ('x2max', repr(c2b)), ('x3min', repr(c3a)), ('x3max', repr(c3b)),
             ('ix2_bc', 'periodic'), ('ox2_bc', 'periodic'), ('ix3_bc', 'periodic'), ('ox3_bc', 'periodic')]:
    s = setkey(s, 'mesh', k, v)
s = setkey(s, 'meshblock', 'nx2', a.n); s = setkey(s, 'meshblock', 'nx3', a.n)
s = setkey(s, 'time', 'start_time', repr(t)); s = setkey(s, 'problem', 'he_ic_remap_file', fn)
open(a.out + '.athinput', 'w').write(s)
print('wrote %s md5 %s and %s.athinput (wedge %.4f..%.4f x %.4f, %dx%d)' % (fn, md5, a.out, c2a, c2b, c3b, a.n, a.n))
