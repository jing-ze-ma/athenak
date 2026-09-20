"""tests_r25: energy budget of the vceil_thermalise arms.  Read-only.

For each arm: d(tot-E) over the run against the net energy the radial faces carried
(problem/face_budget cumulative 'erg in/out'), the event-log totals, and the state of
the cells the ceiling clipped.
"""
import re
import sys
import numpy as np

B = '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r25'
PAT = re.compile(r'erg in=(\S+) out=(\S+).*?\(t = (\S+) s\)')


def arm(d):
    h = np.loadtxt('%s/%s/he4.hydro.hst' % (B, d))
    t, etot = h[:, 0], h[:, 6]
    ein = eout = tf = np.nan
    for line in open('%s/%s/full.log' % (B, d), errors='ignore'):
        m = PAT.search(line)
        if m:
            ein, eout, tf = (float(m.group(1)), float(m.group(2)),
                             float(m.group(3)))
    dE = etot[-1] - etot[0]
    net = ein - eout
    # event log: last row is cumulative only within the output interval, so sum
    g = np.loadtxt('%s/%s/he4.log' % (B, d), ndmin=2)
    nv = g[:, 4].sum() if g.size else 0.0
    efde = g[:, 8].sum() if g.size else 0.0
    vcde = g[:, 11].sum() if g.size > 0 and g.shape[1] > 11 else 0.0
    print('%-6s t=%8.1f  dE=%+11.4e  L-net=%+11.4e  resid=%+11.4e (%7.3f%% of |dE|)'
          % (d, t[-1], dE, net, dE - net,
             100.0*abs(dE - net)/max(abs(dE), 1e-99)))
    print('       eos_vceil=%d  efloor_de(sum)=%+.4e  vceil_de(sum)=%+.4e  E0=%.5e'
          % (nv, efde, vcde, etot[0]))
    return dE, net


for d in sys.argv[1:]:
    arm(d)
