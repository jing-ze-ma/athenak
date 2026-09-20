"""Q2/Q3: structure of the settled envelope vs IC and 2.0 turnovers; transport."""
import sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r26')
from three import (profile, avg, RS, TURN, LW, OMFRAC, SIGSB)  # noqa

TARGET = [0.51, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.05, 1.1, 1.2]
run = sys.argv[1]
SET = [int(x) for x in (sys.argv[2:] or ['16', '17', '18', '19', '20'])]

print("#" * 90)
print("RUN", run)
p0 = profile(run, 0)
p2 = profile(run, 4)
ps = [profile(run, i) for i in SET]
pm = avg(ps)
r = p0['r']
idx = [int(np.argmin(np.abs(r/RS - x))) for x in TARGET]

print("\nIC t=%.3f turn | 2.0-turn dump t=%.3f | settled average of dumps %s "
      "(t = %.2f-%.2f turn)"
      % (p0['t']/TURN, p2['t']/TURN, SET, ps[0]['t']/TURN, ps[-1]['t']/TURN))
print("total wedge mass: IC %.5e g | 2.0t %.5e | settled %.5e (%.4f of IC)"
      % (p0['Mtot'], p2['Mtot'], pm['Mtot'], pm['Mtot']/p0['Mtot']))

hdr = ("%-6s %5s | %10s %9s %6s %7s | %10s %9s %6s %7s | %10s %9s %6s %7s "
       "%8s %8s %8s %7s")
print("\n--- SHELL STRUCTURE (rho, T, beta, M(<r)/M ; then v_rms and Mach settled) ---")
print(hdr % ("r/R", "i", "rho_IC", "T_IC", "beta", "Mf", "rho_2t", "T_2t", "beta",
             "Mf", "rho_set", "T_set", "beta", "Mf", "vr_rms", "vh_rms", "Mach", "v1"))
for x, i in zip(TARGET, idx):
    print(("%-6.2f %5d | %10.3e %9.3e %6.3f %7.4f | %10.3e %9.3e %6.3f %7.4f | "
           "%10.3e %9.3e %6.3f %7.4f %8.2e %8.2e %8.4f %+7.1e")
          % (x, i,
             p0['rho'][i], p0['T'][i], p0['beta'][i], p0['menc'][i]/p0['Mtot'],
             p2['rho'][i], p2['T'][i], p2['beta'][i], p2['menc'][i]/p2['Mtot'],
             pm['rho'][i], pm['T'][i], pm['beta'][i], pm['menc'][i]/pm['Mtot'],
             pm['vr_rms'][i], pm['vh_rms'][i], pm['mach'][i], pm['v1'][i]))

print("\n--- MASS DISTRIBUTION ---")
for nm, p in (('IC', p0), ('2.0t', p2), ('settled', pm)):
    fr = lambda a, b: p['Msh'][(r >= a*RS) & (r < b*RS)].sum()/p['Mtot']
    print("  %-8s M(<0.55R) %.4f  M(<0.675R) %.4f  M(<0.70R) %.4f  M(0.7-1.0R) %.4f"
          "  M(>1.0R) %.4f  M(>1.1R) %.3e"
          % (nm, fr(0, .55), fr(0, .675), fr(0, .70), fr(.70, 1.0), fr(1.0, 9),
             fr(1.1, 9)))

print("\n--- PHOTOSPHERE (tau from the shell-mean Rosseland kappa*rho) ---")
for nm, p in (('IC', p0), ('2.0t', p2), ('settled', pm)):
    k = int(np.argmin(np.abs(np.log(p['tau']/(2./3.)))))
    Lo = 4*np.pi*(p['r'][k])**2*OMFRAC*p['Fdiff'][k]
    Teff_l = (LW/(4*np.pi*(p['r'][k])**2*OMFRAC*SIGSB))**0.25
    print("  %-8s tau=2/3 at i=%d r/R=%.4f  rho=%.3e T=%.4e kappa=%.3e | "
          "T_eff(L_w at that r) = %.0f K | tau(top)=%.2e tau(r=R)=%.3e"
          % (nm, k, p['r'][k]/RS, p['rho'][k], p['T'][k], p['kap'][k], Teff_l,
             p['tau'][-1], p['tau'][int(np.argmin(np.abs(r-RS)))]))

print("\n--- ENTROPY / CONVECTION ZONE (settled; ds/dr < 0 = unstable) ---")
print("%5s %7s %11s %11s %11s %9s %9s %9s"
      % ("i", "r/R", "s_IC", "s_2t", "s_set", "dsdr_set", "vr_rms", "vh_rms"))
ds = np.gradient(pm['s'], r)
ds0 = np.gradient(p0['s'], r)
for i in range(0, len(r), 6):
    print("%5d %7.4f %11.4e %11.4e %11.4e %+9.2e %9.2e %9.2e"
          % (i, r[i]/RS, p0['s'][i], p2['s'][i], pm['s'][i], ds[i],
             pm['vr_rms'][i], pm['vh_rms'][i]))
unst = (ds < 0) & (r > 0.505*RS*1.01) & (r < 1.15*RS)
if unst.any():
    print("  settled ds/dr < 0 between r/R = %.4f and %.4f (%d of %d shells)"
          % (r[unst].min()/RS, r[unst].max()/RS, unst.sum(), len(r)))
unst0 = (ds0 < 0) & (r > 0.52*RS) & (r < 1.15*RS)
if unst0.any():
    print("  IC       ds/dr < 0 between r/R = %.4f and %.4f (%d shells)"
          % (r[unst0].min()/RS, r[unst0].max()/RS, unst0.sum()))
drho = np.gradient(pm['rho'], r)
inv = (drho > 0) & (r > 0.52*RS) & (r < 1.2*RS)
print("  settled density inversion (drho/dr > 0) at r/R: "
      + (", ".join("%.3f" % (x/RS) for x in r[inv][:40]) if inv.any() else "none"))
drho0 = np.gradient(p0['rho'], r)
inv0 = (drho0 > 0) & (r > 0.52*RS) & (r < 1.2*RS)
print("  IC       density inversion at r/R: "
      + (", ".join("%.3f" % (x/RS) for x in r[inv0][:40]) if inv0.any() else "none"))

print("\n--- Q3 TRANSPORT (settled average; fluxes in L_w through 4 pi r^2 x 0.17678) ---")
print("%5s %7s %9s %9s %9s %9s %9s | %7s %8s %8s %9s"
      % ("i", "r/R", "L_diff", "L_mean", "L_enth", "L_kin", "L_tot",
         "P1", "rms_rho", "rms_T", "tau"))
for i in range(0, len(r), 6):
    A = pm['Lsh'][i]
    print("%5d %7.4f %9.3f %9.3f %9.3f %9.3f %9.3f | %7.3f %8.4f %8.4f %9.3e"
          % (i, r[i]/RS, A*pm['Fdiff'][i]/LW, A*pm['Fmean'][i]/LW,
             A*pm['Fenth'][i]/LW, A*pm['Fkin'][i]/LW,
             A*(pm['Fdiff'][i]+pm['Fenth'][i]+pm['Fkin'][i])/LW,
             pm['P1'][i], pm['rrms'][i], pm['trms'][i], pm['tau'][i]))

print("\n--- ceiling cells and dt (settled average of the same dumps) ---")
for p, nm in ((p2, '2.0t'), (pm, 'settled')):
    tot = p['nceil'].sum()
    print("  %-8s cells at >0.99 vceil: %.1f of %d (%.2e); >0.5 vceil: %.1f"
          % (nm, tot, p['ncell']*len(r), tot/(p['ncell']*len(r)), p['nceil5'].sum()))
    k = np.argsort(p['nceil'])[-6:][::-1]
    print("           top radii r/R: "
          + ", ".join("%.3f(%.0f)" % (r[i]/RS, p['nceil'][i]) for i in k if p['nceil'][i] > 0))
for p, nm in zip(ps, SET):
    print("  dump %2d t=%.2ft : dt_r %.3f s (i=%d r/R=%.3f) dt_th %.3f s (i=%d) "
          "dt_ph %.3f s (i=%d r/R=%.3f)"
          % (nm, p['t']/TURN, p['dtmin'][0], p['dtargmin'][0],
             r[int(p['dtargmin'][0])]/RS, p['dtmin'][1], p['dtargmin'][1],
             p['dtmin'][2], p['dtargmin'][2], r[int(p['dtargmin'][2])]/RS))
np.save('/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r26/prof_%s.npy' % run,
        np.array([p0, p2, pm], dtype=object), allow_pickle=True)
