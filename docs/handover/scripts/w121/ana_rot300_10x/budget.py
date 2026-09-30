"""budget.py: hst energy/mass budget per rotation window (the ana_rot300 RESULTS section 4 method).
usage: python budget.py RUNDIR.  dE/dt, dM/dt = linear fits of hydro.hst tot-E, mass over the window;
booked = Etot_bot - Etot_top + Efloor (+Efloor_rt shown separately); mass: -Mdot_top + Mdot_bot + Mfloor."""
import sys
import numpy as np
PR = 1.101535e5
D = sys.argv[1]
h = np.loadtxt(D + '/dhj.hydro.hst')
u = np.loadtxt(D + '/dhj.user.hst')
r = h[:, 0]/PR
KE = h[:, 7:10].sum(1)
print('# win  Lsw_abs Lir_top (Lir/Lsw-1) Lrad_bot/Lsw Efloor Efloor_rt Etot_top Etot_bot Mdot_top'
      ' Mdot_bot Mfloor | dE/dt booked resid | dM/dt Mresid | KEmean')
for lo, hi in [(a, a + 10) for a in range(0, 300, 10)] + [(100, 150), (150, 200), (200, 250),
                                                          (250, 300), (200, 300)]:
    s = (r >= lo) & (r <= hi)
    m = u[s].mean(0)
    dE = np.polyfit(h[s, 0], h[s, 6], 1)[0]
    dM = np.polyfit(h[s, 0], h[s, 2], 1)[0]
    bk = m[7] - m[6] + m[11]
    mb = -m[8] + m[9] + m[12]
    print('%3d-%3d %.4e %.4e %+.2f%% %.2f%% %.2e %.2e %+.2e %+.2e %+.2e %+.2e %.2e | %+.2e %+.2e %+.2e'
          ' | %+.2e %+.2e | %.3e' % (lo, hi, m[4], m[2], 100*(m[2]/m[4] - 1), 100*m[5]/m[4],
                                     m[11], m[13], m[6], m[7], m[8], m[9], m[12], dE, bk,
                                     dE - bk, dM, dM - mb, KE[s].mean()))
