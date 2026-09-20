"""Q1b/Q4: per-dump structural time series + event-counter / dt health."""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r26')
from three import profile, RS, TURN, LW, OMFRAC  # noqa
from common import events, facebudget, rtprof, shellvol  # noqa

run = sys.argv[1]
print("#" * 96)
print("RUN", run)
print("\n--- structural time series (one line per 0.5-turnover dump) ---")
print("%3s %6s %10s %8s %8s %8s %8s %9s %9s %9s %8s %8s %8s"
      % ("d", "turn", "Mtot", "M<0.7R", "M<0.55", "M>1.0R", "r_ph/R", "rho(0.9R)",
         "T(0.9R)", "v1(0.9R)", "vrms0.9", "nceil", "Tmax>1R"))
rows = []
for i in range(21):
    p = profile(run, i)
    r = p['r']
    k9 = int(np.argmin(np.abs(r - 0.9*RS)))
    kph = int(np.argmin(np.abs(np.log(p['tau']/(2./3.)))))
    out = (i, p['t']/TURN, p['Mtot'],
           p['Msh'][r < 0.70*RS].sum()/p['Mtot'],
           p['Msh'][r < 0.55*RS].sum()/p['Mtot'],
           p['Msh'][r > 1.00*RS].sum()/p['Mtot'],
           r[kph]/RS, p['rho'][k9], p['T'][k9], p['v1'][k9], p['vr_rms'][k9],
           p['nceil'].sum(), p['Tmax'][r > RS].max())
    rows.append(out)
    print("%3d %6.2f %10.4e %8.4f %8.4f %8.2e %8.4f %9.2e %9.3e %+9.2e %8.2e %8.0f %8.2e"
          % out)

ev, tev, cl = events(run)
print("\n--- event counters (per %g s event-log interval) ---" % np.median(np.diff(tev)))
print("%6s %9s %9s %9s %9s %9s %11s %11s"
      % ("turn", "dfloor", "efloor", "vceil", "fofc", "tclamp", "efloor_de", "vceil_de"))
for a in np.arange(0.0, 10.0, 0.5):
    s = (tev >= a*TURN) & (tev < (a+0.5)*TURN)
    if not s.any():
        continue
    print("%6.1f %9.0f %9.0f %9.0f %9.0f %9.0f %11.3e %11.3e"
          % (a, ev[s, 1].mean(), ev[s, 2].mean(), ev[s, 4].mean(), ev[s, 7].mean(),
             ev[s, 9].mean(), ev[s, 8].mean(), ev[s, 11].mean()))
s = (tev >= 6*TURN)
print("  6-10 turn totals: efloor_de = %.4e, vceil_de = %.4e (code = cgs energy"
      " DENSITY summed over cells, per interval)" % (ev[s, 8].sum(), ev[s, 11].sum()))
p = np.polyfit(tev[s]/TURN, ev[s, 11], 1)
print("  vceil_de trend 6-10 turn: %+.3e per turnover on a mean of %.3e"
      % (p[0], ev[s, 11].mean()))

print("\n--- dt history (from the cycle log) ---")
for a in np.arange(0.0, 10.0, 1.0):
    s = (cl[:, 1] >= a*TURN) & (cl[:, 1] < (a+1)*TURN)
    if s.any():
        print("  %4.0f-%.0f turn: dt mean %.4f min %.4f max %.4f s, %d cycles"
              % (a, a+1, cl[s, 2].mean(), cl[s, 2].min(), cl[s, 2].max(), s.sum()))
print("  total cycles %d ; cycles in the last turnover %d"
      % (len(cl), ((cl[:, 1] >= 9*TURN)).sum()))
