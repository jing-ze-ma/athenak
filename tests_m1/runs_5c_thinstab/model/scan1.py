import numpy as np
from fourier import gmat, gmax, pde_rate


def g1d(f0, tau, nu, **kw):
    return max(np.max(np.abs(np.linalg.eigvals(gmat(al, 0.0, f0, tau, nu, **kw))))
               for al in np.linspace(0.05, np.pi, 60))


print('PDE: max Re(lambda) over k in [0.01,100]^2, K in {0.1,1,10}, f0 in (0.1..0.99)')
w = -1e9
for f0 in (0.1, 0.5, 0.8, 0.95, 0.99):
    for K in (0.1, 1, 10):
        for kx in np.geomspace(1e-2, 1e2, 25):
            for ky in np.geomspace(1e-2, 1e2, 25):
                for cl in ('m1', 'ker'):
                    w = max(w, pde_rate(kx, ky, f0, K, cl) / K)
print('  max Re(lambda)/K = %.3e' % w)
print('f0  tau | 2-D lag full | od none | dchi off | 1-D lag | ker 2-D | edd 2-D')
for f0 in (0.3, 0.5, 0.65, 0.8, 0.95):
    for tau in (0.05, 0.125, 0.25, 0.5, 1.0):
        r = [gmax(f0, tau, 1e6)[0], gmax(f0, tau, 1e6, od=False)[0],
             gmax(f0, tau, 1e6, dchi_on=False)[0], g1d(f0, tau, 1e6),
             gmax(f0, tau, 1e6, closure='ker')[0], gmax(f0, tau, 1e6, closure='edd')[0]]
        print('%.2f %.3f | ' % (f0, tau) + ' | '.join('%7.3f' % x for x in r))
print('dt dependence, f0=0.5 tau=0.125: nu -> gmax')
for nu in (0.1, 1, 3, 10, 100, 1e4, 1e6):
    print('  nu=%g  g=%.3f' % (nu, gmax(0.5, 0.125, nu)[0]))
