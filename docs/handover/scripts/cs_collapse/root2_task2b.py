import sys
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/cs_ens/ctl2_s01_twin2/analysis')
from root2_fast import read_rows
arm, k, j = 'exp_g3', 16, 17
ii = list(range(48, 59))
rows = {}
for c in range(138769, 138801):
    r = read_rows('%s/cyclediag/dhj.cyclediag.%08d.dat' %
                  (arm, c), ['rt_T', 'w_dens', 'w_velx', 'rt_de', 'u_ener'], k, j)
    rows[c] = r
for var, fmt, lab in [('rt_T', '%9.4g', 'T [K]'), ('w_dens', '%9.3e', 'rho'),
                      ('w_velx', '%9.2e', 'v_r'), ('rt_de', '%9.2e', 'rt_de'),
                      ('delta_u', '%9.2e', 'delta_u')]:
    print('\n--- exp_g3 k=16 j=17 : %s vs i (columns i=48..58) ---' % lab)
    print('%8s %7s ' % ('cycle', 'dt')+' '.join('%9d' % i for i in ii))
    for c in range(138770, 138801):
        r = rows[c]
        if var == 'delta_u':
            v = [r['u_ener'][i]-rows[c-1]['u_ener'][i] for i in ii]
        else:
            v = [r[var][i] for i in ii]
        print('%8d %7.3f ' % (c, r['dt'])+' '.join(fmt % x for x in v))
