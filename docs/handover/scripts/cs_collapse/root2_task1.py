import sys
import os
import glob
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/cs_ens/ctl2_s01_twin2/analysis')
from root2_fast import read_rows
arm, k, j, I = 'si_g6', 17, 2, 68
C0 = 138839
fs = sorted(glob.glob(arm+'/cyclediag/*.dat'))
cyc = [int(os.path.basename(f).split('.')[2]) for f in fs]
pathof = {c: f for c, f in zip(cyc, fs)}
sel = set(c for c in cyc if c >= C0-160 and (c-(C0-160)) % 5 == 0)
sel |= set(c for c in cyc if C0-30 <= c <= C0+15)
sel = sorted(sel)
V = ['u_ener', 'rt_T', 'w_dens', 'w_velx', 'w_eint',
     'rt_src', 'rt_de', 'rt_Fb', 'rt_Ft', 'rt_Em', 'rt_clip']
print('%8s %9s %11s %11s %12s %12s %13s %12s %12s %12s'
      ' %12s %12s %8s %10s %10s %10s %10s' % (
          'cycle', 'dt', 'T[K]', 'rho', 'v_r', 'e_int', 'delta_u', 'rt_src', 'rt_de',
          'rt_Fb',
          'rt_Ft', 'rt_Em', 'rt_clip', 'T[i-4]', 'T[i-2]', 'T[i+2]', 'T[i+4]'))
for c in sel:
    r = read_rows(pathof[c], V, k, j)
    if c-1 in pathof:
        pu = read_rows(pathof[c-1], ['u_ener'], k, j)['u_ener'][I]
        du = r['u_ener'][I]-pu
    else:
        du = float('nan')
    T = r['rt_T']
    print('%8d %9.4f %11.4e %11.4e %12.4e %12.4e %13.5e %12.4e %12.4e'
          ' %12.4e %12.4e %12.4e %8.3f %10.4g %10.4g %10.4g %10.4g' % (
              c, r['dt'], T[I], r['w_dens'][I], r['w_velx'][I], r['w_eint'][I], du,
              r['rt_src'][I], r['rt_de'][I],
              r['rt_Fb'][I], r['rt_Ft'][I], r['rt_Em'][I], r['rt_clip'][I], T[I-4],
              T[I-2],
              T[I+2], T[I+4]))
