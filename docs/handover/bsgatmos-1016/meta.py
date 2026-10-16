"""meta.py RUNDIR GM RPH_RSUN: cells per pressure scale height at R_ph from dump 0
(H_P = P_tot/(rho g), P_tot = (2/3) eint + E/3 [ideal-gas approx], g = GM/r^2)."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bin_convert.py')) else '/raven/u/jinma/ATHENAK/wt-prodfinal/vis/python')
import bin_convert as bc
run, gm, rph = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])*6.957e10
bd = os.path.join(run, 'bin')
h = bc.read_binary_as_athdf(sorted(glob.glob(bd + '/*.hydro_w.*.bin'))[0])
m = bc.read_binary_as_athdf(sorted(glob.glob(bd + '/*.m1.*.bin'))[0])
r = np.asarray(h['x1v'], float); rf = np.asarray(h['x1f'], float)
rho = np.asarray(h['dens'], float).mean(axis=(0, 1)); e = np.asarray(h['eint'], float).mean(axis=(0, 1))
E = np.asarray(m['m1_e'], float).mean(axis=(0, 1))
i = int(np.argmin(abs(r - rph)))
H = (2/3*e[i] + E[i]/3)/(rho[i]*gm/r[i]**2)
dr = rf[i+1] - rf[i]
print('%s: R_ph %.2f Rsun, dr(R_ph) %.4f Rsun, H_P %.4f Rsun, cells per H %.2f, nx1 %d, r %.1f-%.1f Rsun'
      % (os.path.basename(run.rstrip('/')), rph/6.957e10, dr/6.957e10, H/6.957e10, H/dr, len(r), rf[0]/6.957e10, rf[-1]/6.957e10))
