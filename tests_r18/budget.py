#!/usr/bin/env python3
"""Part 2: global energy budget of the wedge from hst + face-budget log lines."""
import numpy as np, re

D = '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r11/wedge9d/'
TURN = 4705.0
LSTAR = 2.3066e38          # full-sphere L printed at setup
OMEGA = 0.17678            # wedge solid-angle share
LW = OMEGA * LSTAR         # wedge luminosity
RIN = 1.18585e11
FIN = 1.305278e15
print('L_wedge (share x 4pi rin^2 F) = %.4e ; share x L_star = %.4e'
      % (OMEGA * 4 * np.pi * RIN**2 * FIN, LW))

# ---- face budget lines from stdout log
pat = re.compile(r'face budget: gain/L\s+inner=(\S+) outer=(\S+) \| erg in=(\S+) '
                 r'out=(\S+) \| g in=(\S+) out=(\S+) \| L_rad,out/L = (\S+) '
                 r'L_rad,cut/L = (\S+) \(t = (\S+) s\)')
rows = []
for line in open(D + '../wedge9d.log', errors='replace'):
    m = pat.search(line)
    if m:
        g = [float(x) for x in m.groups()]
        rows.append(g)
fb = np.array(rows)
print('face-budget samples: %d, t %.1f..%.1f' % (len(fb), fb[0, 8], fb[-1, 8]))
tfb = fb[:, 8]
lout = fb[:, 6] / OMEGA          # L_out / L_wedge
ergin = fb[:, 2]                 # cumulative erg through inner face (hydro)
ergout = fb[:, 3]

# ---- history
h = np.loadtxt(D + 'he4.hydro.hst')
th, dt, mass, m1 = h[:, 0], h[:, 1], h[:, 2], h[:, 3]
etot = h[:, 6]
ke = h[:, 7] + h[:, 8] + h[:, 9]

# ---- event counters
ev = np.loadtxt(D + 'he4.log')
cyc, dfl, efl, vceil, fofc, ede = ev[:, 0], ev[:, 1], ev[:, 2], ev[:, 4], ev[:, 7], ev[:, 8]


def interp(t, x, y):
    return np.interp(t, x, y)


def seg(t0, t1, label):
    e0, e1 = interp(t0, th, etot), interp(t1, th, etot)
    dE = e1 - e0
    dtw = t1 - t0
    dEdt = dE / dtw
    sel = (tfb >= t0) & (tfb <= t1)
    lo = np.mean(lout[sel])
    Lout = lo * LW
    # hydro (advective) energy through the faces: cumulative, differenced
    hin = interp(t1, tfb, ergin) - interp(t0, tfb, ergin)
    hout = interp(t1, tfb, ergout) - interp(t0, tfb, ergout)
    # energy floor input over the interval: ede is per-interval (between log rows)
    s = (cyc >= 0)
    m = (np.interp(t0, th, np.arange(len(th))), )
    print('\n=== %s : t %.0f -> %.0f s  (%.3f -> %.3f turnovers, dtw=%.0f s)'
          % (label, t0, t1, t0 / TURN, t1 / TURN, dtw))
    print('  tot-E  %.6e -> %.6e   dE = %+.4e erg   dE/dt = %+.4e erg/s = %+.3f L_w'
          % (e0, e1, dE, dEdt, dEdt / LW))
    print('  KE     %.4e -> %.4e' % (interp(t0, th, ke), interp(t1, th, ke)))
    print('  1-mom  %.4e -> %.4e' % (interp(t0, th, m1), interp(t1, th, m1)))
    print('  mass   %.6e -> %.6e  (%+.3e)'
          % (interp(t0, th, mass), interp(t1, th, mass),
             interp(t1, th, mass) - interp(t0, th, mass)))
    print('  <L_out/L_w> = %.4f  ->  L_out = %.4e erg/s' % (lo, Lout))
    print('  L_in - L_out = %+.4e erg/s = %+.3f L_w  (integrated %+.4e erg)'
          % (LW - Lout, (LW - Lout) / LW, (LW - Lout) * dtw))
    print('  hydro advective through faces: in %+.3e  out %+.3e erg (cumulative diff)'
          % (hin, hout))
    print('  residual dE/dt - (L_in - L_out) = %+.4e erg/s = %+.3f L_w'
          % (dEdt - (LW - Lout), (dEdt - (LW - Lout)) / LW))
    return dEdt, Lout


seg(11280.0, 14115.0, '2.40 -> 3.00 turnovers')
seg(14115.0, 16464.0, '3.00 -> 3.50 turnovers')
seg(11280.0, 16464.0, '2.40 -> 3.50 turnovers (whole)')

print('\n# L_out/L_w time series (every ~0.1 turnover)')
for w in [2.4 + 0.05 * k for k in range(23)]:
    tt = w * TURN
    if tt > tfb[-1]:
        break
    j = int(np.argmin(np.abs(tfb - tt)))
    k = int(np.argmin(np.abs(th - tt)))
    print('  turn %5.3f  t=%8.1f  L_out/L_w = %6.4f   tot-E = %.5e  KE = %.3e  '
          '1-mom = %+.3e  dt = %5.2f'
          % (w, tfb[j], lout[j], etot[k], ke[k], m1[k], dt[k]))

print('\n# event counters: cumulative efloor_de and vceil over the run')
print('  rows %d, cycles %d..%d' % (len(cyc), cyc[0], cyc[-1]))
print('  sum(efloor_de) over all logged intervals = %.4e erg  (as printed)' % ede.sum())
print('  sum(vceil events) = %.4e' % vceil.sum())
print('  sum(dfloor) = %.4e   sum(fofc) = %.4e' % (dfl.sum(), fofc.sum()))
# split by time: map cycles to time via hst (hst has no cycle column) -> use fraction
