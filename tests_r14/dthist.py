"""tests_r14: dt + floor-event + mass history, and the photosphere/top tracking."""
import sys

import numpy as np

from dtdiag import rd, faces, weight, TURN, RS, ARAD, KB, MU, MH, CFL


def hist(arm):
    h = np.loadtxt(arm + '/he4.hydro.hst')
    e = np.loadtxt(arm + '/he4.log')
    n = min(len(h), len(e))
    print("\n== %s  history (every 47 s); event counters are per-interval" % arm)
    print("  turn   t[s]     dt[s]     mass       tot-E      1-KE      "
          "dfloor efloor tfloor vceil fail fofc tclamp tset")
    for i in range(0, n, 10):
        print("  %5.3f %8.1f %9.3e %10.4e %10.3e %9.3e  %6d %6d %6d %5d %4d %4d %6d %5d"
              % (h[i, 0]/TURN, h[i, 0], h[i, 1], h[i, 2], h[i, 6], h[i, 7],
                 e[i, 1], e[i, 2], e[i, 3], e[i, 4], e[i, 5], e[i, 7],
                 e[i, 9], e[i, 10]))
    print("  last: turn %.3f dt %.3e mass %.4e" % (h[-1, 0]/TURN, h[-1, 1], h[-1, 2]))


def shell(arm):
    R = rd(arm + '/rt_profile.bin')
    r = R[0][1]
    rf = faces(r)
    vol = 4*np.pi/3*(rf[1:]**3 - rf[:-1]**3)
    print("\n== %s  where the dense gas is" % arm)
    print("  turn  rmax(rho>1e-10)/R  rmax(rho>1e-12)/R  rho[top]   v1[top]   "
          "M(>1.0R)[g]  Mdot_top[g/s]  M_tot[g]")
    for t, rr, q in R:
        rho, v1 = q[0], q[1]
        a = np.where(rho > 1e-10)[0]
        b = np.where(rho > 1e-12)[0]
        mtop = (rho*vol)[r > RS].sum()
        mdot = 4*np.pi*rf[-1]**2*rho[-1]*v1[-1]
        print("  %5.3f      %7.4f            %7.4f       %9.3e %9.2e %10.3e %12.3e"
              " %10.4e"
              % (t/TURN, r[a[-1]]/RS if len(a) else 0, r[b[-1]]/RS if len(b) else 0,
                 rho[-1], v1[-1], mtop, mdot, (rho*vol).sum()))


def dmass(arm, t0, t1):
    R = rd(arm + '/rt_profile.bin')
    r = R[0][1]
    rf = faces(r)
    vol = 4*np.pi/3*(rf[1:]**3 - rf[:-1]**3)
    ka = int(np.argmin([abs(x[0]-t0*TURN) for x in R]))
    kb = int(np.argmin([abs(x[0]-t1*TURN) for x in R]))
    ra, rb = R[ka][2], R[kb][2]
    dm = (rb[0]-ra[0])*vol
    print("\n== %s  per-cell mass change %.3f -> %.3f turn  (total %+.4e g)"
          % (arm, R[ka][0]/TURN, R[kb][0]/TURN, dm.sum()))
    print("   i  r/R      dm[g]      rho_a      rho_b      T_a       eint_a    w_a"
          "    rho_min=(10/9)e/cs2  v1_a")
    for i in np.argsort(-np.abs(dm))[:20]:
        rmin = (10.0/9.0)*ra[6][i]/1.0e16
        print("  %3d %7.4f %+10.3e %10.3e %10.3e %9.3e %9.3e %5.2f  %12.3e %10.2e"
              % (i, r[i]/RS, dm[i], ra[0][i], rb[0][i], ra[5][i], ra[6][i],
                 weight(np.array([ra[0][i]]))[0], rmin, ra[1][i]))


if __name__ == '__main__':
    a = sys.argv[1]
    hist(a)
    shell(a)
    if len(sys.argv) > 3:
        dmass(a, float(sys.argv[2]), float(sys.argv[3]))
