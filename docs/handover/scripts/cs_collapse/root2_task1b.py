import sys
import os
import glob
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/cs_ens/ctl2_s01_twin2/analysis')
from root2_fast import read_rows
arm, k, j, I = 'si_g6', 17, 2, 68
fs = sorted(glob.glob(arm+'/cyclediag/*.dat'))
p = {int(os.path.basename(f).split('.')[2]): f for f in fs}
print('--- gap fill 138855-138858 ---')
print('%8s %8s %11s %11s %12s %12s %13s %12s %12s %8s %10s %10s %10s %10s' % (
    'cycle', 'dt', 'T[K]', 'rho', 'v_r', 'e_int', 'delta_u', 'rt_src', 'rt_de',
    'rt_clip', 'T[i-4]', 'T[i-2]', 'T[i+2]', 'T[i+4]'))
for c in range(138855, 138859):
    r = read_rows(p[c], ['u_ener', 'rt_T', 'w_dens', 'w_velx',
                  'w_eint', 'rt_src', 'rt_de', 'rt_clip'], k, j)
    du = r['u_ener'][I]-read_rows(p[c-1], ['u_ener'], k, j)['u_ener'][I]
    T = r['rt_T']
    print('%8d %8.4f %11.4e %11.4e %12.4e %12.4e %13.5e %12.4e'
          ' %12.4e %8.3f %10.4g %10.4g %10.4g %10.4g' % (
              c, r['dt'], T[I], r['w_dens'][I], r['w_velx'][I], r['w_eint'][I], du,
              r['rt_src'][I], r['rt_de'][I], r['rt_clip'][I], T[I-4], T[I-2], T[I+2],
              T[I+4]))
print('--- earlier history at i=68, every 20 cycles 138354-138679 ---')
print('%8s %11s %11s %12s %10s %10s' %
      ('cycle', 'T[K]', 'rho', 'v_r', 'T[i+2]', 'T[i+4]'))
for c in range(138354, 138680, 20):
    r = read_rows(p[c], ['rt_T', 'w_dens', 'w_velx'], k, j)
    print('%8d %11.4e %11.4e %12.4e %10.4g %10.4g' %
          (c, r['rt_T'][I], r['w_dens'][I], r['w_velx'][I], r['rt_T'][I+2],
              r['rt_T'][I+4]))
