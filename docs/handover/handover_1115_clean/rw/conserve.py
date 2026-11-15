"""conserve.py RUNDIR...: total energy and momentum (gas + radiation) of periodic rw runs over time.
E_tot = sum(eint + rho v^2/2 + E_rad), M_tot = sum(rho v + F/c^2); prints the drift relative to the start."""
import glob, re, sys
import numpy as np

def totals(d, c):
    out = []
    for fh in sorted(glob.glob(d + '/tab/rw.hydro_w.*.tab')):
        fr = fh.replace('hydro_w', 'm1')
        h, r = np.loadtxt(fh), np.loadtxt(fr)
        rho, v, eint = h[:, 3], h[:, 4], h[:, 7]
        erad, frad = r[:, 3], r[:, 4]
        out.append((np.sum(eint + 0.5*rho*v*v + erad), np.sum(rho*v + frad/c**2),
                    np.sum(np.abs(eint) + np.abs(erad)), np.sum(np.abs(rho*v) + np.abs(frad)/c**2)))
    return np.array(out)

for d in sys.argv[1:]:
    c = float(re.search(r'c_light\s*=\s*(\S+)', open(d + '/rw.athinput').read()).group(1))
    t = totals(d, c)
    de = np.max(np.abs(t[:, 0] - t[0, 0]))/t[0, 2]
    dm = np.max(np.abs(t[:, 1] - t[0, 1]))/max(t[:, 3].max(), 1e-300)
    print('%-60s energy drift %.2e  momentum drift %.2e  (%d outputs)' % (d.split('/')[-1][:60], de, dm, len(t)))
