import sys
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
from read_cyclediag import load
arm = 'exp_g3'
k, j = 16, 17
cyc = range(138769, 138801)
ii = list(range(48, 59))
prev = None
for c in cyc:
    d = load('%s/cyclediag/dhj.cyclediag.%08d.dat' % (arm, c))
    ue = d['u_ener'][k, j]
    if prev is not None and c >= 138770:
        print('=== cycle %d  time %.6f  dt %.5f' % (c, d['time'], d['dt']))
        print('%4s %11s %11s %12s %13s %13s' %
              ('i', 'T[K]', 'rho', 'v_r', 'rt_de', 'delta_u'))
        for i in ii:
            print('%4d %11.4e %11.4e %12.4e %13.5e %13.5e' % (
                i, d['rt_T'][k, j, i], d['w_dens'][k, j, i], d['w_velx'][k, j, i],
                d['rt_de'][k, j, i], ue[i]-prev[i]))
    prev = ue
