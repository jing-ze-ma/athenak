"""model numbers next to the measured runs: nu = 300, r = 1 (dy = dx1), wm = 1."""
import numpy as np
from fourier import gmax, gmat


def gmax_modes(f0, tau, nu, alphas, betas, **kw):
    return max(np.max(np.abs(np.linalg.eigvals(gmat(a, b, f0, tau, nu, **kw))))
               for a in alphas for b in betas if (a, b) != (0, 0))


al = np.linspace(0, np.pi, 33)
print('tau  | f0=0.3  0.4   0.5   (all k)   | f0=0.5 by seed beta (ny=16): k=1,2,4,8')
for tau in (0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0):
    r = [gmax(f, tau, 300.0, n=32)[0] for f in (0.3, 0.4, 0.5)]
    rk = [gmax_modes(0.5, tau, 300.0, al, [2 * np.pi * k / 16]) for k in (1, 2, 4, 8)]
    print('%.4f | ' % tau + ' '.join('%5.2f' % x for x in r) + '   | '
          + ' '.join('%5.2f' % x for x in rk))
print('nu scan f0=0.5 tau=0.125:', ' '.join('%g:%.2f' % (nu, gmax(0.5, 0.125, nu, n=32)[0])
                                         for nu in (1, 3, 30, 300, 3000, 3e4)))
for cl in ('ker', 'edd'):
    print('closure', cl, 'tau=0.125 f0=0.5: %.2f' % gmax(0.5, 0.125, 300.0, n=32, closure=cl)[0])
print('od none: %.2f   wm=0: %.2f' % (gmax(0.5, 0.125, 300.0, n=32, od=False)[0],
                                     gmax(0.5, 0.125, 300.0, n=32, wm=0.0)[0]))
print('implicit (converged Picard = nonlinear BE): %.3f'
      % gmax(0.5, 0.125, 300.0, n=32, mode='imp')[0])
