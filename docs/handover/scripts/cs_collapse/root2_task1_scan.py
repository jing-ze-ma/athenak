import sys
import glob
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/cs_ens/ctl2_s01_twin2/analysis')
from root2_fast import read_rows
import numpy as np
arm, k, j, thr = 'si_g6', 17, 2, 5000.
fs = sorted(glob.glob(arm+'/cyclediag/*.dat'))
r0 = read_rows(fs[0], ['rt_T'], k, j)
is_, ie = r0['is'], r0['ie']
hot0 = set(int(i) for i in np.where(r0['rt_T'][is_:ie+1] > thr)[0]+is_)
print('cycle %d: cells >5000 K already: %s' % (r0['cycle'], sorted(hot0)))
first = None
for f in fs:
    r = read_rows(f, ['rt_T'], k, j)
    T = r['rt_T']
    cand = [i for i in range(is_, ie+1) if i not in hot0 and T[i] > thr]
    if cand:
        print('FIRST NEW crossing cycle %d cells %s  T=%s' %
              (r['cycle'], cand, [float('%.5g' % T[i]) for i in cand]))
        first = (r['cycle'], cand)
        break
if first is None:
    print('no new crossing in %d files (last cycle %d)' % (len(fs), r['cycle']))
