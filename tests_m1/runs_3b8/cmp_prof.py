"""Horizontal-mean / fluctuation profiles of two seeded-slab arms, side by side."""
import sys
import glob
import numpy as np
sys.path.insert(0, '../../vis/python')
import bin_convert as bc  # noqa: E402

c = 2.99792458e10
arms = sys.argv[1:]


def load(d, n):
    w = bc.read_binary_as_athdf(sorted(glob.glob(d + '/bin/*hydro_w*.bin'))[n])
    r = bc.read_binary_as_athdf(sorted(glob.glob(d + '/bin/*.m1.*.bin'))[n])
    return w, r


nd = min(len(glob.glob(a + '/bin/*hydro_w*.bin')) for a in arms)
for n in range(1, nd):
    print('==== dump %d' % n)
    for a in arms:
        w, r = load(a, n)
        rho, v1, v2, e = w['dens'][0], w['velx'][0], w['vely'][0], w['eint'][0]
        E, F1, F2 = r['m1_e'][0], r['m1_f1'][0], r['m1_f2'][0]
        m1 = (rho*v1).mean(0)/rho.mean(0)
        ke_mean = (rho.mean(0)*m1**2).sum()
        ke_fluc = (rho*(v1 - m1)**2).mean(0).sum()
        print('%-8s t=%5.1f  KE1 mean-part %.2e  fluct-part %.2e' %
              (a, w['Time'], ke_mean, ke_fluc))
        for i in (83, 80, 76, 72, 68, 64, 58, 50, 40, 25, 5):
            print('   i=%2d <v1> %+9.2e  rms(v1-<v1>) %8.2e  rms v2 %8.2e  '
                  'dE/E %8.2e  de/e %8.2e  <f1> %.4f  max f2 %.2e'
                  % (i, m1[i], np.sqrt(((v1[:, i] - m1[i])**2).mean()),
                     np.sqrt((v2[:, i]**2).mean()), E[:, i].std()/E[:, i].mean(),
                     e[:, i].std()/e[:, i].mean(), (F1[:, i]/(c*E[:, i])).mean(),
                     (np.abs(F2[:, i])/(c*E[:, i])).max()))
