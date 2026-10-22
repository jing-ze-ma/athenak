#!/usr/bin/env python3
"""photinst-1022: interior cube-pattern check vs time.  For each hydro_w (+m1) dump, at r/Rph0 in RADII:
  v_r: (mean_FC - mean_ED)/rms and rms_FC/rms_ED - 1 (cell-size-dependent convective speed);
  Trad and e/rho: fc/ed - 1 in units of the lateral sd (z = (mean_FC - mean_ED)/sd).
usage: interior.py RUNDIR [every]"""
import sys
import os
import glob
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert as bc  # noqa: E402
A, RPH0 = 7.565733e-15, 5.80019e13
RADII = (0.3, 0.5, 0.7, 0.8, 0.9, 0.95)


def assemble(d, var):
    nmb = len(d['mb_logical'])
    per = nmb // 6
    n3, n2, n1 = d['mb_data'][var][0].shape
    nb = int(round(np.sqrt(per)))
    N = n2 * nb
    P = np.empty((6, N, N, n1), np.float64)
    for m in range(nmb):
        p = m // per
        lx1, lx2, lx3, lev = d['mb_logical'][m]
        P[p, lx3 * n3:(lx3 + 1) * n3, lx2 * n2:(lx2 + 1) * n2, :] = d['mb_data'][var][m]
    return P


run = sys.argv[1]
every = int(sys.argv[2]) if len(sys.argv) > 2 else 1
files = sorted(glob.glob(os.path.join(run, 'bin', 'rsg.hydro_w.*.bin')))[::every]
print('# t_d | per r/Rph0: vr z(FC-ED) rmsFC/rmsED-1 vr_rms[cm/s] | Trad z | e/rho z')
for fh in files:
    h = bc.read_binary(fh)
    t = float(h['time']) / 86400
    r = np.asarray(h['mb_x1v'][0], float)
    rho = assemble(h, 'dens')
    er = assemble(h, 'eint') / rho
    vr = assemble(h, 'velx')
    fm = fh.replace('hydro_w', 'm1')
    Tr = (np.maximum(assemble(bc.read_binary(fm), 'm1_e'), 1e-30) / A) ** 0.25 if os.path.exists(fm) else None
    N = rho.shape[1]
    x2c = -1 + (np.arange(N) + 0.5) * 2.0 / N
    ax2 = np.abs(x2c)[None, :]
    ax3 = np.abs(x2c)[:, None]
    FC = np.broadcast_to((ax2 < 0.3) & (ax3 < 0.3), (6, N, N))
    ED = np.broadcast_to(np.maximum(ax2, ax3) > 0.8, (6, N, N))
    out = []
    for x in RADII:
        i = np.searchsorted(r, x * RPH0)
        v = vr[..., i]
        rms = np.sqrt((v ** 2).mean())
        s = '%.2f: %+.2f %+.2f %.1e' % (x, (v[FC].mean() - v[ED].mean()) / max(rms, 1e-30),
                                        np.sqrt((v[FC] ** 2).mean()) / max(np.sqrt((v[ED] ** 2).mean()), 1e-30) - 1, rms)
        for q in (Tr, er):
            if q is None:
                s += ' | -'
                continue
            a = q[..., i]
            s += ' | %+.2f' % ((a[FC].mean() - a[ED].mean()) / max(a.std(), 1e-300))
        out.append(s)
    print('%6.2f || %s' % (t, ' || '.join(out)), flush=True)
