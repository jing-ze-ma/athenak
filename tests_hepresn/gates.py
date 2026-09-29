#!/usr/bin/env python3
"""Scalar gates of the M1 (frozen-MLT ramp) 3-D he_star_m1 wedge runs.
usage: gates.py RUNDIR [RUNDIR ...] [--tt 4700] [--every 0.5]
Concatenates the user.hst files (restart chain order as given, later runs override
overlapping times) and the run.log files; prints per window of --every turnovers:
M_tot, mass change per turnover (%), Mdot_top (wedge, g/s), L_top/L_in, L_MLT/L_in, w,
mean dt, Picard, sqrt(v1sq_wall), counts of NEWTON-FALLBACK, Picard NON-CONVERGED and
hesdirk2 stage messages (from run.log, by the time of the preceding cycle line)."""
import glob
import re
import sys

import numpy as np


def read_hst(fn):
    lab = None
    for ln in open(fn):
        if ln.startswith('#') and '[1]=' in ln:
            lab = re.findall(r'\[\d+\]=(\S+)', ln)
            break
    d = np.loadtxt(fn, comments='#', ndmin=2)
    return {k: d[:, n] for n, k in enumerate(lab)}


def read_log(fn):
    t, ev = 0.0, []
    for ln in open(fn, errors='replace'):
        m = re.match(r'elapsed=\S+ cycle=(\d+) time=(\S+) dt=(\S+)', ln)
        if m:
            t = float(m.group(2))
            continue
        if 'NEWTON-FALLBACK' in ln:
            ev.append((t, 0))
        elif 'NON-CONVERGED' in ln:
            ev.append((t, 1))
        elif 'hesdirk2 stage' in ln:
            ev.append((t, 2))
    return ev


def main():
    runs = [a for a in sys.argv[1:] if not a.startswith('--') and not a[0].isdigit()]
    tt = float(sys.argv[sys.argv.index('--tt') + 1]) if '--tt' in sys.argv else 4700.0
    evy = float(sys.argv[sys.argv.index('--every') + 1]) if '--every' in sys.argv else 0.5
    H, ev = None, []
    for r in runs:
        h = read_hst(r + '/hepresn.user.hst')
        if H is None:
            H = h
        else:
            keep = H['time'] < h['time'][0]
            H = {k: np.concatenate([H[k][keep], h[k]]) for k in h if k in H}
        for lf in sorted(glob.glob(r + '/run*.log')):
            ev += read_log(lf)
    ev = np.array(ev) if ev else np.zeros((0, 2))
    t = H['time']
    Lin = H['L_in']
    M0 = H['M_tot'][0]
    print('M_tot(0) = %.4e g, L_in = %.4e erg/s (wedge), turnover %.0f s'
          % (M0, Lin[0], tt))
    print('   t/tt      M_tot/M0  dM/M%%/tt  Mdot_top   Ltop/Lin  LMLT/Lin   w    '
          '<dt>   Picard  vwall    NFB   NC  stg')
    edges = np.arange(0.0, t[-1] + evy * tt, evy * tt)
    for a, b in zip(edges[:-1], edges[1:]):
        s = (t >= a) & (t < b)
        if s.sum() < 2:
            continue
        ts = t[s]
        m = H['M_tot'][s]
        dm = (m[-1] - m[0]) / M0 * 100.0 / ((ts[-1] - ts[0]) / tt)
        e = ev[(ev[:, 0] >= a) & (ev[:, 0] < b)] if len(ev) else ev
        cnt = [int((e[:, 1] == q).sum()) if len(e) else 0 for q in (0, 1, 2)]
        pic = H['Picard'][s]
        pic = pic[pic > 0]
        w = H['w_mlt'][s].mean() if 'w_mlt' in H else float('nan')
        lm = (H['L_MLT'][s] / Lin[s]).mean() if 'L_MLT' in H else float('nan')
        print('%6.2f-%5.2f %9.6f %+8.3f %10.3e %8.4f %8.4f %6.3f %6.1f %6.2f %8.2e '
              '%5d %4d %4d'
              % (a / tt, b / tt, m[-1] / M0, dm, H['Mdot_top'][s].mean(),
                 (H['L_top'][s] / Lin[s]).mean(), lm, w, H['dt'][s].mean(),
                 pic.mean() if len(pic) else 0.0,
                 np.sqrt(max(H['v1sq_wall'][s].mean(), 0.0)), *cnt))


if __name__ == '__main__':
    main()
