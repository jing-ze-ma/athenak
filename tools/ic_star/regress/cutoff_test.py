#!/usr/bin/env python3
"""cutoff_test.py RUN RPH T1 T2 [RUN RPH T1 T2 ...]: dominant period of the shell-mean <v_r>
(Hann-windowed, 32x zero-padded FFT, peak per radius) over t in [T1, T2] days, for r/R_ph in
0.99..1.03, and the saturated amplitude (max and rms |<v_r>| in the window, 1.0-1.03 R_ph)."""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..',
                                'vis', 'python'))
import bin_convert as bc  # noqa: E402


def series(run):
    ts, V = [], []
    for f in sorted(glob.glob(os.path.join(run, 'bin', '*.hydro_w.*.bin'))):
        h = bc.read_binary_as_athdf(f)
        d = np.asarray(h['dens'], float)
        ts.append(float(h['Time']))
        V.append((d*np.asarray(h['velx'], float)).mean(axis=(0, 1))/d.mean(axis=(0, 1))/1e5)
        r = np.asarray(h['x1v'], float)
    return np.array(ts), np.array(V), r


def main():
    a = sys.argv[1:]
    for q in range(0, len(a), 4):
        run, R, t1, t2 = a[q], float(a[q+1]), float(a[q+2]), float(a[q+3])
        ts, V, r = series(run)
        m = (ts >= t1*86400) & (ts <= t2*86400)
        t = ts[m]
        per = []
        for x in (0.99, 1.0, 1.01, 1.02, 1.03):
            i = int(np.argmin(abs(r - x*R)))
            v = V[m, i] - V[m, i].mean()
            n = 32*len(v)
            p = np.abs(np.fft.rfft(v*np.hanning(len(v)), n=n))**2
            f = np.fft.rfftfreq(n, d=t[1] - t[0])
            k = 1 + int(np.argmax(p[1:]))
            per.append(1/f[k]/86400)
        band = (r >= R) & (r <= 1.03*R)
        vb = np.abs(V[m][:, band])
        print('%s: t %g-%g d, R_ph %.4e; dominant period at r/R_ph 0.99/1.00/1.01/1.02/1.03: %s d; '
              'median(1.00-1.03) %.3f d; |<v_r>| 1.00-1.03 R_ph max %.1f rms %.1f km/s'
              % (run, t1, t2, R, '/'.join('%.2f' % p_ for p_ in per), np.median(per[1:]),
                 vb.max(), np.sqrt((vb**2).mean())))


if __name__ == '__main__':
    main()
