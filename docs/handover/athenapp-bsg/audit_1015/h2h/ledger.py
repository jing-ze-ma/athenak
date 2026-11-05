import numpy as np, sys, glob
EU = 1e-9*(6.957e10)**3*(6.957e10/1.899769e4)**2   # code energy unit (rho_u l_u^3 v_u^2) [erg]
SR = 4*np.pi/((np.cos(np.pi/3)-np.cos(2*np.pi/3))*(np.pi/3))   # wedge -> full sphere
NAMES = ['hydflux','gravmass','pointmass','radT_bill','radT_cons','clip','reset','floorE','floorM','radD','gasbill']
def load(d):
    a = np.loadtxt(d + '/out/ledger.txt')
    a = a[np.all(np.isfinite(a), axis=1)]
    return a
for d in sys.argv[1:]:
    a = load(d); t = a[:, 0]/4.547922
    print('== %s  (t %.3f .. %.2f d), energies in 1e45 erg full sphere, change since first row' % (d, t[0], t[-1]))
    for r, rn in enumerate(['r<30 Rsun', 'cells 0-10', 'all']):
        b = 1 + r*12
        Em = (a[:, b] - a[0, b])*EU*SR/1e45
        T = (a[:, b+1:b+12] - a[0, b+1:b+12])*EU*SR/1e45
        clos = T[:, [0, 1, 2, 6, 7, 8, 9, 10]].sum(axis=1)
        for tt in [1, 2, 5, 10, 20, 30]:
            i = np.searchsorted(t, tt)
            if i >= len(t): continue
            x = T[i]
            print('%-10s t %5.2f dE_meas %+9.3f closure_sum %+9.3f (resid %+.2e) | hyd %+8.3f grav(mass+pm) %+8.3f [mass %+8.3f pm %+8.3f] '
                  'radT_bill %+8.3f radT_cons %+8.3f nontelesc %+8.3f clip %+7.3f reset %+7.3f floor %+7.3f'
                  % (rn, t[i], Em[i], clos[i], (Em[i]-clos[i])/max(abs(Em[i]), 1e-30), x[0], x[1]+x[2], x[1], x[2],
                     x[3], x[4], x[3]-x[4], x[5], x[6], x[7]+x[8]))
