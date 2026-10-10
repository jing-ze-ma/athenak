"""order_1009: grey atmosphere (G1v) vs the exact Hopf solution and self-convergence"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import od

def hopf(z):
    tau = 0.128*0.113*(np.exp((1.0 - z)/0.113) - 1.0)
    return 3.0*(tau + 0.710446 - 0.133054*np.exp(-3.4488*tau)), tau

for arm in sys.argv[1:] or ('cen', 'hr', 'hrb', 'hrx0'):
    e1, em, prev = [], [], None
    for n in od.NX:
        q = od.cols(od.rundir('atm', arm, 'be', 'S', n))
        z, E = q['x1v'], q['m1_e']
        Ex, tau = hopf(z)
        r = np.abs(E/Ex - 1)
        # restricted to tau > 0.1 / tau < 0.1
        e1.append(np.abs(E - Ex).sum()/np.abs(Ex).sum())
        print('%s n=%4d  L1 %.3e  max|dE/E| %.3e  max(tau<1) %.3e  max(tau>1) %.3e  top %.3e  tau_cell top %.1e bot %.1e' % (
            arm, n, e1[-1], r.max(), r[tau < 1].max(), r[tau > 1].max(), r[-1],
            0.128/n, 0.128*np.exp(1/0.113)/n))
    e1 = np.array(e1)
    print(arm, 'Hopf L1 orders', ' '.join('%.2f' % x for x in np.log2(e1[:-1]/e1[1:])))
    res, levs = od.evaluate('atm', arm, 'be', 'S')
    for k, e in res.items():
        e = np.array(e)
        print(arm, 'self', k, ' '.join('%.2e' % x for x in e), '| ord', ' '.join('%.2f' % x for x in np.log2(e[:-1]/e[1:])))
