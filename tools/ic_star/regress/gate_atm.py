#!/usr/bin/env python3
"""gate_atm.py RUNDIR RPH [--bands a-b,c-d,...] [--vmax 1] [--dmax 0.1]: gate (a) atmosphere check.
Per hydro_w dump: per band (r/R_ph) the max over the band's radii of |<v_r>| (mass-weighted shell
mean, km/s) and of |<rho>/<rho>_0 - 1| (shell means vs the first dump).  PASS for the band
0.99-1.1 R_ph: |<v_r>| < vmax and |drho| < dmax at every dump, and no growth (the max over the
last third of the run not above 1.5x the max over the middle third).
Since 2026-10-17 (USER) the ADOPTED atmosphere gate is gate (a'), gate_aprime.py (README sect. 8a);
this table stays as a diagnostic."""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..',
                                'vis', 'python'))
import bin_convert as bc  # noqa: E402


def shell(f):
    h = bc.read_binary_as_athdf(f)
    x2f = np.asarray(h['x2f'], float)
    dom = (np.cos(x2f[:-1]) - np.cos(x2f[1:]))[None, :, None] * \
        np.diff(np.asarray(h['x3f'], float))[:, None, None]
    rho = np.asarray(h['dens'], float)
    w = np.broadcast_to(dom, rho.shape)
    vr = np.asarray(h['velx'], float)
    m = (rho*w).sum(axis=(0, 1))
    return (float(h['Time']), np.asarray(h['x1v'], float), m/w.sum(axis=(0, 1)),
            (rho*vr*w).sum(axis=(0, 1))/m/1e5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run')
    ap.add_argument('rph', type=float)
    ap.add_argument('--bands', default='0.5-0.9,0.9-0.99,0.99-1.01,1.01-1.03,1.03-1.1,0.99-1.1')
    ap.add_argument('--vmax', type=float, default=1.0)
    ap.add_argument('--dmax', type=float, default=0.1)
    a = ap.parse_args()
    fs = sorted(glob.glob(os.path.join(a.run, 'bin', '*.hydro_w.*.bin')))
    bands = [tuple(float(x) for x in b.split('-')) for b in a.bands.split(',')]
    t0, r, rho0, _ = shell(fs[0])
    x = r/a.rph
    rows = []
    print('# %s: %d dumps, R_ph %.6e; per band max|<v_r>| [km/s] / max|<rho>/<rho>_0 - 1|'
          % (a.run, len(fs), a.rph))
    print('# t[d]   ' + '  '.join('%9.3f-%-5.3f' % b for b in bands))
    for f in fs:
        t, _, rho, vr = shell(f)
        row = []
        for lo, hi in bands:
            k = (x >= lo) & (x < hi)
            row.append((np.abs(vr[k]).max() if k.any() else np.nan,
                        np.abs(rho[k]/rho0[k] - 1).max() if k.any() else np.nan))
        rows.append((t, row))
        print('%7.3f  ' % (t/86400) + '  '.join('%6.3f/%-7.4f' % v for v in row))
    ib = len(bands) - 1
    v = np.array([rw[ib][0] for _, rw in rows])
    d = np.array([rw[ib][1] for _, rw in rows])
    n = len(v)
    grow = v[2*n//3:].max() > 1.5*max(v[n//3:2*n//3].max(), 1e-3) if n >= 6 else False
    ok = (v.max() < a.vmax) and (d.max() < a.dmax) and not grow
    print('# band %.2f-%.2f R_ph: max|<v_r>| %.3f km/s (at t %.3f d), max|drho| %.4f, growth %s '
          '-> %s' % (bands[ib][0], bands[ib][1], v.max(), rows[int(np.argmax(v))][0]/86400,
                     d.max(), grow, 'PASS' if ok else 'FAIL'))


if __name__ == '__main__':
    main()
