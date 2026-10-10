"""grey atmosphere: per-region SELF-convergence (the Hopf reference has a quadrature/Marshak floor
in the thin top).  e_n = sum|E_n - R E_2n| / sum|E_2n| over the band (tau from the top), R the
2:1 restriction.  Bands by tau, and by the Knudsen number Kn = |dE/dz|/(chi E) of the finest run.
usage: ana_atm_self.py arm ..."""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import od  # noqa: E402

TB = [(0, 0.1), (0.1, 1.0), (1.0, 10.0), (10.0, 1e9)]
KB = [(0, 0.03), (0.03, 0.3), (0.3, 1e9)]


def tau_of(z):
    return 0.128*0.113*(np.exp((1.0 - z)/0.113) - 1.0)


def chi_of(z):
    return 0.128*np.exp((1.0 - z)/0.113)


for arm in sys.argv[1:]:
    qs = []
    for n in od.NX:
        try:
            q = od.cols(od.rundir('atm', arm, 'be', 'S', n))
            qs.append((q['x1v'], q['m1_e']))
        except Exception:
            qs.append(None)
    print(arm)
    # Kn of each level-n coarse cell from the finer run (restricted)
    for name, bands, key in (('tau', TB, 'tau'), ('Kn', KB, 'kn')):
        for lo, hi in bands:
            e = []
            for a in range(len(qs) - 1):
                if qs[a] is None or qs[a+1] is None:
                    e.append(np.nan)
                    continue
                z, Ec = qs[a]
                zf, Ef = qs[a+1]
                Er = 0.5*(Ef[0::2] + Ef[1::2])
                if key == 'tau':
                    v = tau_of(z)
                else:
                    g = np.abs(np.gradient(Er, z))
                    v = g/(chi_of(z)*Er)
                s = (v >= lo) & (v < hi)
                e.append(np.abs(Ec[s] - Er[s]).sum()/max(np.abs(Er[s]).sum(), 1e-300) if s.any()
                         else np.nan)
            o = [np.log2(e[i]/e[i+1]) for i in range(len(e) - 1)]
            print('  %-3s %-10s e %s | ord %s' % (name, '%g-%g' % (lo, hi),
                                                 ' '.join('%.2e' % x for x in e),
                                                 ' '.join('%5.2f' % x for x in o)))
