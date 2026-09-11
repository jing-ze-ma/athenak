import sys
import glob
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/cs_ens/ctl2_s01_twin2/analysis')
from root2_fast import read_full
import numpy as np


def counts(arm, step):
    fs = sorted(glob.glob(arm+'/cyclediag/*.dat'))
    out = []
    for f in fs[::step]+([fs[-1]] if (len(fs)-1) % step else []):
        r = read_full(f, ['rt_T', 'w_dens'])
        s = (slice(r['ks'], r['ke']+1), slice(r['js'],
             r['je']+1), slice(r['is'], r['ie']+1))
        T = r['rt_T'][s]
        rho = r['w_dens'][s]
        out.append((r['cycle'], r['time'], int(np.sum(np.abs(T-200.) < 0.5)),
                   int(np.sum(rho < 1e-9)), float(T.max()), float(rho.min())))
    return out


res = {a: counts(a, 20) for a in ('exp_g3', 'si_g6')}
print('active cells per block = 16*16*128 = %d' % (16*16*128))
for a in ('exp_g3', 'si_g6'):
    print('\n=== arm %s (every 20 cycles) ===' % a)
    print('%8s %13s %12s %12s %12s %12s' %
          ('cycle', 'time', 'n(T=200K)', 'n(rho<1e-9)', 'max T[K]', 'min rho'))
    for c, t, n1, n2, tm, rm in res[a]:
        print('%8d %13.6f %12d %12d %12.4e %12.4e' % (c, t, n1, n2, tm, rm))
