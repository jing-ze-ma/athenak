import sys
import os
import glob
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/cs_ens/ctl2_s01_twin2/analysis')
from root2_fast import read_full
import numpy as np
arm = 'exp_g3'
fs = sorted(glob.glob(arm+'/cyclediag/*.dat'))
cyc = [int(os.path.basename(f).split('.')[2]) for f in fs]
p = {c: f for c, f in zip(cyc, fs)}
sel = [c for c in cyc if 138700 <= c <= 138780 and c %
       10 == 0]+[c for c in cyc if 138780 < c <= 138867]
print('%8s %13s %11s %11s %11s %11s %11s' %
      ('cycle', 'time', 'n(T=200K)', 'n(rho<1e-9)', 'n(rt_clip)', 'sum|clip|', 'maxT[K]'))
for c in sorted(sel):
    r = read_full(p[c], ['rt_T', 'w_dens', 'rt_clip'])
    s = (slice(r['ks'], r['ke']+1), slice(r['js'], r['je']+1), slice(r['is'], r['ie']+1))
    T = r['rt_T'][s]
    rho = r['w_dens'][s]
    cl = r['rt_clip'][s]
    print('%8d %13.4f %11d %11d %11d %11.4g %11.4e' % (
        c, r['time'], int(np.sum(np.abs(T-200.) < 0.5)), int(np.sum(rho < 1e-9)),
        int(np.sum(cl != 0)), float(np.abs(cl).sum()), float(T.max())))
