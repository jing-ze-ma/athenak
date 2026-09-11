import sys
import glob
import os
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
from read_cyclediag import load
import numpy as np
arm = sys.argv[1]
k = int(sys.argv[2])
j = int(sys.argv[3])
thr = float(sys.argv[4])
fs = sorted(glob.glob(arm+'/cyclediag/*.dat'))
for f in fs:
    d = load(f)
    T = d['rt_T'][k, j, d['is']:d['ie']+1]
    m = np.where(T > thr)[0]
    if m.size:
        print('first crossing cycle %d file %s' % (d['cycle'], os.path.basename(f)))
        for ii in m:
            print('  i=%d T=%.5e rho=%.5e' %
                  (ii+d['is'], T[ii], d['w_dens'][k, j, ii+d['is']]))
        break
else:
    print('none above', thr)
