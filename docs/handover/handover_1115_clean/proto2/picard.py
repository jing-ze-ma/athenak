"""picard.py: does the lagged limiter/closure Picard iteration converge per stage?
(a) per-stage relative change after passes 1,2,3,5,10 and passes to 1e-12 (npic='conv'),
    scattering slab (beam into chi dx 10) and smooth thin wave, c dt/dx 0.5..1000;
(b) slab |F|-cE violations vs npic = 1, 3, 10, conv."""
import sys
import numpy as np
import scheme as S
import tests as T


def stats(H):
    H = [h for h in H if len(h)]
    out = []
    for k in (1, 2, 3, 5, 10):
        v = [h[k - 1] for h in H if len(h) >= k]
        out.append('%.0e/%.0e' % (np.median(v), np.max(v)) if v else '   -   ')
    n = np.array([len(h) for h in H])
    nc = sum(1 for h in H if h[-1] >= 1e-12)
    return ' '.join(out), int(np.median(n)), n.max(), nc, len(H)


def slab_run(opt, cfl, chidx=10.0, n=100):
    x = (np.arange(n) + 0.5)/n
    pb = S.Problem(n, ks=chidx*n*(x > 0.5), bc='inflow')
    worst = [-np.inf]

    def chk(U, Y1):
        for V in (U, Y1):
            if V is not None:
                worst[0] = max(worst[0], (np.abs(V[n:]) - V[:n]).max(), -V[:n].min())
    h = cfl/n
    S.run(pb, np.zeros(n), np.zeros(n), max(0.6, 8*h), h, opt, check=chk)
    return worst[0]


def wave_run(opt, cfl, n=128):
    pb = S.Problem(n, ks=0.5*6.4, ka=0.5*6.4, B=1.0)
    x = pb.x
    S.run(pb, 1 + 0.3*np.sin(2*np.pi*x), 0.2*np.sin(2*np.pi*x + 0.5), max(0.5, 8*cfl/n),
          cfl/n, opt)


if __name__ == '__main__':
    lims = [None, 'minmod', 'mc', 'vl']
    print('(a) per-stage rel. change after pass 1,2,3,5,10 (median/max over stages); passes to'
          ' 1e-12 median, max; stages not converged in 60 passes')
    for case, fn in (('slab chi dx 10', slab_run), ('thin wave chi dx 0.05', wave_run)):
        for lim in lims:
            for cfl in (0.5, 4.0, 100.0, 1000.0):
                S.HIST = []
                fn(dict(cand='A3', recon=lim, Dguard=True, npic='conv'), cfl)
                st = stats(S.HIST)
                print('%-22s %-7s cfl %-6g %s | passes %d max %d | unconverged %d/%d'
                      % (case, lim or '1st', cfl, st[0], st[1], st[2], st[3], st[4]))
                sys.stdout.flush()
    S.HIST = None
    print('(b) slab worst max(|F|-cE, -E) vs npic (A3, Dg, no fallback)')
    for chidx in (1.0, 10.0, 100.0):
        for lim in lims:
            for cfl in (0.5, 4.0, 100.0, 1000.0):
                r = [slab_run(dict(cand='A3', recon=lim, Dguard=True, npic=k), cfl, chidx)
                     for k in (1, 3, 10, 'conv')]
                print('slab chi dx %-4g %-7s cfl %-6g npic 1 %+.2e  3 %+.2e  10 %+.2e  conv %+.2e'
                      % (chidx, lim or '1st', cfl, *r))
                sys.stdout.flush()
