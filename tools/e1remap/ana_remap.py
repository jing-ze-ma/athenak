#!/usr/bin/env python3
"""ana_remap.py RUNDIR [W] (e1remap-1014): start-up transient of a remapped narrow-wedge run, from bin/*.hydro_w.*.bin
(+ agcar3d.x1grid.txt). Per dump: t, rms v_r [km/s] in 0.4-0.95 R_ph (convective envelope) and 1.0-1.2 R_ph
(surface), separately in the SEAM band (cells within 3 of a periodic seam, theta or phi) and the INTERIOR
(>= W+4 cells from both seams); lateral KE (sum rho (v_th^2 + v_ph^2)/2 dV, wedge) and radial KE; max |v_r|;
rho seam jump metric (rms jump across the phi / theta seam / rms adjacent jump) over 0.4-1.2 R_ph.
R_ph = 2.7022e13 cm (AG Car A IC)."""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bin_convert as bc  # noqa: E402

RPH = 2.7022e13
D = sys.argv[1]
W = int(sys.argv[2]) if len(sys.argv) > 2 else 8
g = np.array([[float(v) for v in l.split()] for l in open(os.path.join(D, 'agcar3d.x1grid.txt'))
              if not l.startswith('#') and len(l.split()) == 4])
hdr = open(os.path.join(D, 'agcar3d.x1grid.txt')).readline().split()
is_ = int(hdr[3].rstrip(";"))
ie_ = int(hdr[6].rstrip(";"))
x1v = g[is_:ie_ + 1, 2]
x1f = np.append(g[is_:ie_ + 1, 1], g[ie_, 1] + g[ie_, 3])


def glob3(x, v):
    n1, n2, n3 = x['Nx1'], x['Nx2'], x['Nx3']
    out = np.empty((n3, n2, n1))
    b1, b2, b3 = x['nx1_out_mb'], x['nx2_out_mb'], x['nx3_out_mb']
    for m, a in enumerate(x['mb_data'][v]):
        lx1, lx2, lx3 = x['mb_logical'][m][:3]
        out[lx3 * b3:(lx3 + 1) * b3, lx2 * b2:(lx2 + 1) * b2, lx1 * b1:(lx1 + 1) * b1] = a
    return out


fs = sorted(glob.glob(os.path.join(D, 'bin', '*.hydro_w.*.bin')))
print('# %s  W=%d  dumps %d' % (D, W, len(fs)))
print('#  t(s)        cyc  vr_env_seam vr_env_int vr_sur_seam vr_sur_int  max|vr|  KE_lat     KE_rad     rho_seam_phi rho_seam_th')
for f in fs:
    x = bc.read_binary(f)
    d = glob3(x, 'dens')
    vr = glob3(x, 'velx')
    vt = glob3(x, 'vely')
    vp = glob3(x, 'velz')
    n3, n2, n1 = d.shape
    th_f = np.linspace(x['x2min'], x['x2max'], n2 + 1)
    wth = np.cos(th_f[:-1]) - np.cos(th_f[1:])
    dph = (x['x3max'] - x['x3min']) / n3
    vol = (wth[None, :, None] * dph) * ((x1f[1:] ** 3 - x1f[:-1] ** 3) / 3.0)[None, None, :]
    k = np.arange(n3)
    j = np.arange(n2)
    dk = np.minimum(k, n3 - 1 - k)
    dj = np.minimum(j, n2 - 1 - j)
    dist = np.minimum(dk[:, None], dj[None, :])
    seam = dist < 3
    inner = dist >= W + 4

    def rms(mask_lat, rlo, rhi):
        ri = (x1v >= rlo * RPH) & (x1v < rhi * RPH)
        q = vr[:, :, ri][mask_lat]
        return np.sqrt(np.mean(q ** 2)) / 1e5
    ri = (x1v >= 0.4 * RPH) & (x1v < 1.2 * RPH)
    ds = d[:, :, ri]

    def sm(axis):
        N = ds.shape[axis]
        s = np.take(ds, 0, axis=axis) - np.take(ds, N - 1, axis=axis)
        a = np.diff(ds, axis=axis)
        return np.sqrt(np.mean(s ** 2)) / np.sqrt(np.mean(a ** 2))
    kel = float(np.sum(0.5 * d * (vt ** 2 + vp ** 2) * vol))
    ker = float(np.sum(0.5 * d * vr ** 2 * vol))
    print('%.6e %6d %10.4f %10.4f %10.4f %10.4f %9.3f %.3e %.3e %8.3f %8.3f' %
          (x['time'], x['cycle'], rms(seam, 0.4, 0.95), rms(inner, 0.4, 0.95), rms(seam, 1.0, 1.2),
           rms(inner, 1.0, 1.2), np.max(np.abs(vr[:, :, ri])) / 1e5, kel, ker, sm(0), sm(1)))
