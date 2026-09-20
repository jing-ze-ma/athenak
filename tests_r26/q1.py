"""Q1: steadiness -- L_out/L_w statistics, energy budget, mass loss, KE/eint trends."""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r26')
from common import (RUNS, TURN, LW, hst, facebudget, rtprof, shellvol, stats)  # noqa

np.set_printoptions(precision=4, suppress=False)

for run in ('lr1t', 'lr_f100t'):
    print("\n" + "=" * 78)
    print("RUN", run)
    tf, Ein, Eout, Lout, gin, gout = facebudget(run)
    tt = tf / TURN
    h = hst(run)
    th, dt_h, mass, p1, etot = h[:, 0], h[:, 1], h[:, 2], h[:, 3], h[:, 6]
    ke = h[:, 7] + h[:, 8] + h[:, 9]
    print("  t_end = %.1f s = %.3f turnovers ; %d face-budget samples, %d hst rows"
          % (tf[-1], tt[-1], len(tf), len(th)))

    # ---- L_out/L_w statistics per window
    print("\n  L_out/L_w statistics (from L_rad,out/L x L_star / L_w):")
    print("   %-12s %6s %8s %8s %8s %8s %8s" %
          ("window[turn]", "n", "mean", "rms", "min", "max", "median"))
    for a, b in [(0, 0.2), (1.9, 2.1), (2, 5), (5, 6), (6, 7), (7, 8), (8, 9),
                 (9, 10), (6, 10), (8, 10)]:
        s = (tt >= a) & (tt < b)
        if s.sum() < 2:
            continue
        st = stats(Lout[s])
        print("   %-12s %6d %8.4f %8.4f %8.4f %8.4f %8.4f"
              % ("%.1f-%.1f" % (a, b), st['n'], st['mean'], st['rms'],
                 st['min'], st['max'], st['med']))

    # spikes
    s = (tt >= 6) & (tt <= 10)
    hi = np.where(Lout[s] > 3.0)[0]
    ts, Ls = tt[s], Lout[s]
    if len(hi):
        print("  samples with L_out/L_w > 3 in 6-10 turn: %d of %d" % (len(hi), s.sum()))
        # group into bursts
        grp = np.split(hi, np.where(np.diff(hi) > 3)[0] + 1)
        for g in grp[:14]:
            print("     burst t = %.4f-%.4f turn (%d samples, %.0f s), peak %.2f L_w"
                  % (ts[g[0]], ts[g[-1]], len(g),
                     (ts[g[-1]] - ts[g[0]]) * TURN, Ls[g].max()))
    else:
        print("  no samples above 3 L_w in 6-10 turnovers")

    # periodogram on 6-10
    y = Ls - Ls.mean()
    dtm = np.median(np.diff(ts[np.isfinite(Ls)])) * TURN
    yy = np.nan_to_num(y)
    F = np.fft.rfft(yy * np.hanning(len(yy)))
    fr = np.fft.rfftfreq(len(yy), dtm)
    P = np.abs(F)**2
    k = np.argsort(P[1:])[-5:][::-1] + 1
    print("  top periods in L_out (6-10 turn, dt_samp %.1f s): " % dtm
          + ", ".join("%.0f s (%.3f turn, pow %.2g)"
                      % (1 / fr[i], 1 / fr[i] / TURN, P[i]) for i in k))

    # ---- energy budget per turnover
    print("\n  ENERGY BUDGET per turnover (units L_w); "
          "dE/dt from hst col7, L_out from face budget, hyd_out = d(erg out)/dt")
    print("   %-9s %8s %8s %8s %8s %9s" %
          ("window", "dE/dt", "L_in", "-L_out", "hyd_out", "residual"))
    for a in np.arange(2.0, 10.0, 1.0):
        b = a + 1.0
        t0, t1 = a * TURN, b * TURN
        dtw = t1 - t0
        dE = (np.interp(t1, th, etot) - np.interp(t0, th, etot)) / dtw
        lo = np.nanmean(Lout[(tt >= a) & (tt < b)]) * LW
        ho = (np.interp(t1, tf, Eout) - np.interp(t0, tf, Eout)) / dtw
        res = dE - (LW - lo + ho)
        print("   %-9s %8.3f %8.3f %8.3f %8.3f %9.3f"
              % ("%.0f-%.0f" % (a, b), dE / LW, 1.0, -lo / LW, ho / LW, res / LW))

    # ---- mass
    print("\n  MASS: M(0) = %.6e g ; M(2 turn) = %.6e ; M(end) = %.6e"
          % (mass[0], np.interp(2 * TURN, th, mass), mass[-1]))
    for a, b in [(0, 2), (2, 5), (5, 8), (6, 10), (8, 10)]:
        m0, m1 = np.interp(a * TURN, th, mass), np.interp(b * TURN, th, mass)
        print("     %.0f-%.0f turn: dM = %+.4e g  (%+.4f %% of M0)  rate %.3e g/s"
              " = %.3e Msun/yr" % (a, b, m1 - m0, 100 * (m1 - m0) / mass[0],
                                   (m1 - m0) / ((b - a) * TURN),
                                   (m1 - m0) / ((b - a) * TURN) * 3.156e7 / 1.989e33))
    print("     hyd mass out (face budget g out, cum) = %.3e g ; g in = %.3e g"
          % (gout[-1], gin[-1]))

    # ---- KE / eint / potential trends
    print("\n  TRENDS (hst; KE1 = radial, KE2+KE3 = horizontal):")
    print("   %-7s %10s %10s %10s %10s %10s %8s" %
          ("t[turn]", "tot-E", "KEr", "KEh", "1-mom", "mass", "dt"))
    for a in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
        i = np.argmin(np.abs(th - a * TURN))
        print("   %-7.1f %10.4e %10.4e %10.4e %+10.3e %10.4e %8.3f"
              % (th[i] / TURN, etot[i], h[i, 7], h[i, 8] + h[i, 9], p1[i],
                 mass[i], dt_h[i]))
    for a, b in [(6, 10), (8, 10)]:
        s = (th >= a * TURN) & (th <= b * TURN)
        for nm, y in [('tot-E', etot), ('KEr', h[:, 7]), ('KEh', h[:, 8] + h[:, 9]),
                      ('1-mom', p1), ('mass', mass)]:
            p = np.polyfit(th[s] / TURN, y[s], 1)
            print("   fit %s on %d-%d turn: slope %+.4e per turnover"
                  " (%+.3f %%/turn of the mean)"
                  % (nm, a, b, p[0], 100 * p[0] / np.mean(np.abs(y[s]))))
        break

    # eint from rt_profile (shell means)
    t, r, A = rtprof(run)
    vol, rf = shellvol(r)
    Eint = (A[:, 6, :] * vol).sum(axis=1)
    Mp = (A[:, 0, :] * vol).sum(axis=1)
    print("\n  eint (rt_profile shell means x vol): "
          + ", ".join("%.1ft: %.4e" % (tx, Eint[np.argmin(np.abs(t - tx * TURN))])
                      for tx in (0, 2, 5, 6, 8, 10)))
    print("  mass  (rt_profile):                  "
          + ", ".join("%.1ft: %.4e" % (tx, Mp[np.argmin(np.abs(t - tx * TURN))])
                      for tx in (0, 2, 5, 6, 8, 10)))
