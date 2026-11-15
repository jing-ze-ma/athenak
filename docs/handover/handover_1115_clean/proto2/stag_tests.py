"""stag_tests.py: tests of the staggered scheme (stag.py), same setups as tests.py.
Usage: python3 stag_tests.py pos | thin <chi> | quick | front"""
import sys
import itertools
import numpy as np
from multiprocessing import Pool
import scheme as S
import stag as G
import tests as T

OPT = dict(recon='vl', mood=True, lo_mode='char')          # A3 + van Leer + limiter
OPT1 = dict(recon=None, mood=True, lo_mode='char')         # first order


def setup(test, chidx, kind, n=100):
    """Problem and initial (E, F at faces) for the realizability cases."""
    dx = 1.0/n
    chi = chidx/dx
    x = (np.arange(n) + 0.5)*dx
    ka, ks = (chi, 0.0) if kind == 'abs' else (0.0, chi)
    if test in ('beam', 'slab'):
        if test == 'slab':
            ka, ks = ka*(x > 0.5), ks*(x > 0.5)
        pb = S.Problem(n, ka=ka, ks=ks, bc='inflow')
        return pb, np.zeros(n), np.zeros(n + 1)
    pb = S.Problem(n, ka=ka, ks=ks)
    xf = np.arange(n)*dx
    if test == 'pulse':
        m, mf = (x > 0.2) & (x < 0.4), (xf > 0.2) & (xf < 0.4)
        return pb, 1.0*m, 1.0*mf
    m1, m2 = (x > 0.15) & (x < 0.3), (x > 0.7) & (x < 0.85)
    f1, f2 = (xf > 0.15) & (xf < 0.3), (xf > 0.7) & (xf < 0.85)
    return pb, 1.0*(m1 | m2), 1.0*f1 - 1.0*f2


def pos_case(arg):
    name, (test, chidx, kind, cfl) = arg
    opt = {'vl': OPT, 'vlupw': OPT, '1st': OPT1}[name]
    G.FACECHECK = 'upw' if name == 'vlupw' else 'avg'
    for k in G.STATS:
        G.STATS[k] = 0
    pb, E, F = setup(test, chidx, kind)
    w = [0.0, -np.inf, -np.inf, -np.inf]

    def chk(U, Y1):
        for V in (U, Y1):
            if V is None:
                continue
            e, f, _ = G.split(pb, V)
            if not np.all(np.isfinite(V)):
                w[:] = [-np.inf, np.inf, np.inf, np.inf]
                return
            w[0] = min(w[0], e.min())
            w[1] = max(w[1], (np.abs(G.cellF(pb, f)) - e).max())
            w[2] = max(w[2], (np.abs(f) - G.faceE(pb, e, f)).max())
            w[3] = max(w[3], (np.abs(f) - G.upwindE(pb, e, f)).max())
    h = cfl/100
    G.run(pb, E, F, max(0.6, 8*h), h, opt, check=chk)
    return name, (test, chidx, kind, cfl), w, dict(G.STATS)


def pos():
    cases = [(t, x, k, c) for t in ('beam', 'pulse', 'counter', 'slab')
             for x in (0.0, 0.01, 0.1, 1.0, 10.0, 100.0)
             for k in (('scat',) if x == 0 else ('scat', 'abs'))
             for c in (0.5, 4.0, 100.0, 1000.0)]
    with Pool(32) as p:
        out = p.map(pos_case, list(itertools.product(('vl', 'vlupw'), cases)), chunksize=1)
    for name in ('vl', 'vlupw'):
        rows = [r for r in out if r[0] == name]
        bad = [r for r in rows if r[2][0] < -1e-12 or max(r[2][1], r[2][2]) > 1e-12]
        badu = [r for r in rows if r[2][0] < -1e-12 or r[2][3] > 1e-12]
        print('stag %-5s violations of minE + upwind-face test only: %d, worst face-upw %.2e'
              % (name, len(badu), max(r[2][3] for r in rows)))
        st = {k: sum(r[3][k] for r in rows) for k in rows[0][3]}
        print('stag %-4s cases %d violations %d worst minE %.2e max(|Fcell|-cE) %.2e '
              'max(|F|-cEface) %.2e | steps %d marked %.2f/step re-solves %.2f/stage BE %.1f%% '
              'all-low %d BE still bad %d'
              % (name, len(rows), len(bad), min(r[2][0] for r in rows),
                 max(r[2][1] for r in rows), max(r[2][2] for r in rows), st['steps'],
                 st['marked']/st['steps'], st['resolves']/st['stages'],
                 100.0*st['be']/st['steps'], st['allfaces'], st['befail']))
        for r in bad:
            print('    %-7s chi dx %-6g %-4s cfl %-6g minE %+.2e cell %+.2e face %+.2e upw %+.2e'
                  % (*r[1], *r[2]))
        sys.stdout.flush()


def thin_run(n, chidx64, cfl, opt, tend=0.5):
    chi = chidx64*64.0
    pb = S.Problem(n, ks=0.5*chi, ka=0.5*chi, B=1.0)
    xf = np.arange(n)*pb.dx
    U = G.run(pb, 1.0 + 0.3*np.sin(2*np.pi*pb.x), 0.2*np.sin(2*np.pi*xf + 0.5), tend,
              cfl*pb.dx, opt)
    return U[:n], U[n:]


def thin(chidx, opt=OPT):
    Er, Fr = thin_run(2048, chidx, 0.5, opt)
    errs = []
    for n in (64, 128, 256):
        E, F = thin_run(n, chidx, 0.5, opt)
        errs.append(np.abs(E - T.coarsen(Er, n)).mean() + np.abs(F - Fr[::2048//n]).mean())
    o = np.log2(np.array(errs[:-1])/np.array(errs[1:]))
    print('stag thin wave chi dx(64) %-5g L1(E)+L1(F) %s orders %s | %s'
          % (chidx, ' '.join('%.2e' % e for e in errs), ' '.join('%.2f' % v for v in o), G.STATS))


def thick_run(n, opt, nst_mode='dx', nst=None, chi=6400.0):
    pb = S.Problem(n, ks=chi)
    x, xf = pb.x, np.arange(n)*pb.dx
    kap = 1.0/(3.0*chi)
    s0 = 0.05
    tend = 0.5*s0**2/kap
    if nst is None:
        nst = 16*(n//64) if nst_mode == 'dx' else 16*(n//64)**2
    E = np.exp(-0.5*((x - 0.5)/s0)**2)
    F = kap*(xf - 0.5)/s0**2*np.exp(-0.5*((xf - 0.5)/s0)**2)
    U = G.run(pb, E, F, tend, tend/nst, opt)
    sig2 = s0**2 + 2*kap*tend
    ex = s0/np.sqrt(sig2)*np.exp(-0.5*(x - 0.5)**2/sig2)
    E = U[:n]
    var = np.sum(E*(x - 0.5)**2)/np.sum(E)
    return np.abs(E - ex).mean(), (var - s0**2)/(2*kap*tend) - 1.0, E


def quick():
    for name, opt in (('vl', OPT), ('1st', OPT1)):
        for mode in ('dx', 'dx2'):
            r = [thick_run(n, opt, mode) for n in (64, 128, 256)]
            e = np.array([q[0] for q in r])
            print('stag %-4s thick dt~%-3s L1 %s orders %s kappa_num/kappa-1 %s'
                  % (name, mode, ' '.join('%.2e' % v for v in e),
                     ' '.join('%.2f' % v for v in np.log2(e[:-1]/e[1:])),
                     ' '.join('%+.3f' % q[1] for q in r)))
            sys.stdout.flush()
    # time order: thin wave N 128 chi dx(64) 0.1 and thick Gaussian N 128
    sol = [thin_run(128, 0.1, cfl, OPT, tend=1.0)[0] for cfl in (4, 2, 1, 0.5, 0.25)]
    d = [np.abs(sol[k] - sol[k + 1]).mean() for k in range(4)]
    print('stag vl time order thin N 128 cfl 4..0.25 diffs %s orders %s'
          % (' '.join('%.2e' % v for v in d),
             ' '.join('%.2f' % v for v in np.log2(np.array(d[:-1])/np.array(d[1:])))))
    sol = [thick_run(128, OPT, nst=k)[2] for k in (8, 16, 32, 64, 128)]
    d = [np.abs(sol[k] - sol[k + 1]).mean() for k in range(4)]
    print('stag vl time order thick N 128 c dt/dx 384..24 diffs %s orders %s'
          % (' '.join('%.2e' % v for v in d),
             ' '.join('%.2f' % v for v in np.log2(np.array(d[:-1])/np.array(d[1:])))))
    # odd-even
    for cfl in (0.5, 100.0):
        n = 64
        pb = S.Problem(n, ks=100.0*n)
        sgn = (-1.0)**np.arange(n)
        U = np.r_[1.0 + 0.1*sgn, np.zeros(n)]
        amps = []
        for _ in range(20):
            U, _ = G.step(pb, U, cfl*pb.dx, OPT)
            amps.append(abs(np.mean(U[:n]*sgn)))
        print('stag vl odd-even chi dx 100 cfl %-5g 0.1 -> step1 %.2e step5 %.2e step20 %.2e'
              % (cfl, amps[0], amps[4], amps[-1]))
    # conservation
    for lab, ks in (('vacuum pulse', 0.0), ('scattering pulse chi dx 1', 100.0)):
        pb, E, F = setup('pulse', ks/100, 'scat')
        for cfl in (0.5, 100.0):
            rec = []
            G.run(pb, E, F, max(0.6, 8*cfl/100), cfl/100, OPT,
                  check=lambda U, Y: rec.append((U[:100].sum(), U[100:].sum())))
            dE = max(abs(r[0] - E.sum()) for r in rec)/E.sum()
            dF = max(abs(r[1] - F.sum()) for r in rec)/F.sum()
            print('stag vl %-26s cfl %-5g rel drift sum E %.1e  sum F %.1e%s'
                  % (lab, cfl, dE, dF, '' if ks == 0 else ' (F has the -c chi F sink)'))
    # radiating wave with gas, D fixed 1/3 (analytic dispersion)
    for ka in (1.0, 100.0):
        errs = []
        for n in (32, 64, 128):
            pb, lam, v, _ = T.rwave_setup(n, ka=ka)
            xf = np.arange(n)*pb.dx
            eps = 1e-3
            ph, phf = np.exp(2j*np.pi*pb.x), np.exp(2j*np.pi*xf)
            U = G.run(pb, 1 + eps*np.real(v[0]*ph), eps*np.real(v[1]*phf), 1.0, 0.5*pb.dx,
                      dict(OPT, mood=False, Dfix=1.0/3.0), T=1 + eps*np.real(v[2]*ph))
            Eh = np.sum((U[:n] - 1.0)*np.conj(ph))*2/n/eps
            ln = np.log(Eh/v[0])
            ln += 2j*np.pi*np.round((lam.imag - ln.imag)/(2*np.pi))
            errs.append((ln, abs(ln - lam)))
        e = np.array([q[1] for q in errs])
        print('stag vl rwave (D 1/3) ka %-5g exact %.5f%+.5fi numeric %s |err| %s orders %s'
              % (ka, lam.real, lam.imag, ' '.join('%.5f%+.5fi' % (q[0].real, q[0].imag)
                                                  for q in errs),
                 ' '.join('%.1e' % v for v in e),
                 ' '.join('%.2f' % v for v in np.log2(e[:-1]/e[1:]))))


def front():
    """Beam into vacuum (front width) and beam onto an absorbing slab (transmitted energy)."""
    n = 200
    for cfl in (0.5, 4.0):
        for lab, ka in (('vacuum', 0.0), ('absorber chi dx 10 at x > 0.5', 10.0*n)):
            x = (np.arange(n) + 0.5)/n
            kav = ka*(x > 0.5)
            pb = S.Problem(n, ka=kav, bc='inflow')
            Ug = G.run(pb, np.zeros(n), np.zeros(n + 1), 0.4, cfl/n, OPT)
            Uc = S.run(pb, np.zeros(n), np.zeros(n), 0.4, cfl/n,
                       dict(cand='A3', recon='vl', mood=True, lo_mode='char', mood_all=True))
            for nm, E in (('stag', Ug[:n]), ('cell', Uc[:n])):
                if ka == 0:
                    i1, i9 = np.argmax(E < 0.9), np.argmax(E < 0.1)
                    print('%s cfl %-4g beam into vacuum t 0.4 (front at x 0.4*0.989): 0.9->0.1 '
                          'width %d cells, E(0.2) %.4f' % (nm, cfl, i9 - i1, E[int(0.2*n)]))
                else:
                    print('%s cfl %-4g beam onto absorber: E just inside slab %.3e, E 5 cells in '
                          '%.3e, max E behind %.3e' % (nm, cfl, E[n//2], E[n//2 + 5],
                                                     E[n//2 + 10:].max()))


if __name__ == '__main__':
    if sys.argv[1] == 'pos':
        pos()
    elif sys.argv[1] == 'thin':
        thin(float(sys.argv[2]))
    elif sys.argv[1] == 'quick':
        quick()
    elif sys.argv[1] == 'front':
        front()
