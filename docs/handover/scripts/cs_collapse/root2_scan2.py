import sys
import glob
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
from read_cyclediag import load
import numpy as np
arm = sys.argv[1]
k = int(sys.argv[2])
j = int(sys.argv[3])
thr = float(sys.argv[4])
fs = sorted(glob.glob(arm+'/cyclediag/*.dat'))
d0 = load(fs[0])
is_, ie = d0['is'], d0['ie']
T0 = d0['rt_T'][k, j, is_:ie+1]
hot0 = set(np.where(T0 > thr)[0]+is_)
print('cells >%.0f K at first cycle %d: %s' % (thr, d0['cycle'], sorted(hot0)))
for f in fs:
    d = load(f)
    T = d['rt_T'][k, j]
    cand = [i for i in range(is_, ie+1) if i not in hot0 and T[i] > thr]
    if cand:
        print('first NEW crossing: cycle %d  cells %s' % (d['cycle'], cand))
        for i in cand:
            print('  i=%d T=%.5e rho=%.5e' % (i, T[i], d['w_dens'][k, j, i]))
        break
else:
    print('no new crossing')
