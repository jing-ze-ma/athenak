"""Exact (transport) J/B for the two-level-atom atmosphere of the nltek decks.

tau(x) = 1e-3 (exp(10 - x) - 1), x in [-10, 10]; S = (1 - eps) J + eps B, B = 1;
no incoming radiation at the top, I = B entering from below at tau_max.
J = Lambda[S] + (1/2) E2(tau_max - tau), Lambda with piecewise-constant S per cell
(exact E2 integrals), solved directly. Compares with the Eddington formula
1 - exp(-sqrt(3 eps) tau)/(1 + sqrt(eps)) used by the scorer and with code output."""
import sys, glob
import numpy as np
from scipy.special import expn

def exact(eps, n):
    tmax = 1e-3*(np.exp(20.0) - 1.0)
    edges = np.concatenate([[0.0], np.logspace(-7, np.log10(tmax), n)])
    tc = 0.5*(edges[1:] + edges[:-1])
    a, b = edges[:-1], edges[1:]
    d1 = tc[:, None] - a[None, :]
    d2 = tc[:, None] - b[None, :]
    # (1/2) int_a^b E1(|t - tau|) dt
    L = 0.5*np.where(d2 >= 0, expn(2, d2) - expn(2, d1),
                     np.where(d1 <= 0, expn(2, -d1) - expn(2, -d2),
                              2.0 - expn(2, d1) - expn(2, -d2)))
    bnd = 0.5*expn(2, tmax - tc)
    A = np.eye(n) - (1.0 - eps)*L
    J = np.linalg.solve(A, eps*L.sum(axis=1) + bnd)
    return tc, J

def edd(eps, tau):
    return 1.0 - np.exp(-np.sqrt(3*eps)*tau)/(1.0 + np.sqrt(eps))

if __name__ == "__main__":
    sys.path.insert(0, '/viper/ptmp2/jinma/clean_1114/ana')
    if __name__ == "__main__":
      pass
    for eps, deck in ((0.1, 'nltek1_FIN3_d'), (1e-4, 'nltek4_FIN3_d')):
        t1, J1 = exact(eps, 1200)
        t2, J2 = exact(eps, 2400)
        print('eps %g  exact J(0)/B = %.5f (sqrt(eps) %.5f)  Eddington(0) = %.5f' %
              (eps, J2[0], np.sqrt(eps), edd(eps, 0.0)))
        ex = lambda t, tt=t2, JJ=J2: np.interp(np.log(t), np.log(tt), JJ)
        chk = np.abs(np.interp(np.log(t1), np.log(t2), J2) - J1).max()
        print('   reference self-check (N 1200 vs 2400) max |dJ| %.1e' % chk)
        for tq in (1e-2, 0.1, 1, 3, 10, 100):
            print('   tau %-6g exact %.5f  eddington %.5f  ratio %.4f' % (tq, ex(tq), edd(eps, tq), edd(eps, tq)/ex(tq)))
        import s2
        for arm in ('L26_x1', 'L26_x4', 'CC5_x1', 'CC5_x2'):
            d = '/viper/ptmp2/jinma/clean_1114/suite/%s/%s' % (deck, arm)
            if not glob.glob(d + '/**/*.tab', recursive=True):
                continue
            m = s2.nlte(d)
            s = (m['tau'] > 1e-3) & (m['tau'] < 1e5)
            e = ex(m['tau'][s])
            print('   %-7s L1 vs Eddington %.4f | L1 vs exact transport %.4f | E/J_exact top %.4f' %
                  (arm, m['L1'], np.abs(m['E'][s] - e).sum()/e.sum(), m['E'][-1]/ex(m['tau'][-1])))
