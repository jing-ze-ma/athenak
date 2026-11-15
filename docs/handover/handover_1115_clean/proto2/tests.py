"""tests.py: test drivers for scheme.py.  Usage: python3 tests.py <test> [args]
  pos       realizability sweep (test 1), all variants, parallel; -> pos.out
  thin      smooth thin wave order (test 2)
  torder    time order at fixed grid (SDIRK2 check)
  thick     thick diffusion limit (test 3) + compact-flux check
  oddeven   checkerboard decay (test 4)
  rwave     linear radiating wave with gas (test 5)
  cons      conservation (test 6)"""
import sys
import time
import numpy as np
from multiprocessing import Pool
import scheme as S

VARIANTS = [dict(cand=c, recon=r, rc='oma') for c in ('A1', 'A2', 'A2max', 'A3', 'A3b')
            for r in (None, 'minmod')]
VARIANTS += [dict(cand=c, recon=r, rc=w) for c in ('A1', 'A3') for r in (None, 'minmod')
             for w in ('omt', 'one')]
VARIANTS += [dict(v, Dguard=True) for v in VARIANTS[:10]]
VARIANTS += [dict(v, Dguard=True, fallback=True) for v in VARIANTS[:10]]
VARIANTS += [dict(cand='one', recon=None, Dguard=True, fallback=True, aspeed='c'),
             dict(cand='A1', recon='minmod', Dguard=True, fallback=True, aspeed='c'),
             dict(cand='A3', recon='minmod', Dguard=True, fallback=True, aspeed='c')]


def vname(o):
    return '%s/%s/rc=%s%s%s%s%s' % (o['cand'], o.get('recon') or '1st', o.get('rc', 'oma'),
                                  '/np%s' % o['npic'] if 'npic' in o else '',
                                  '/Dg' if o.get('Dguard') else '',
                                  '/fb' if o.get('fallback') else '',
                                  '/a=c' if o.get('aspeed') == 'c' else '')


# ------------------------------------------------------------------ test 1
def pos_case(test, chidx, kind, cfl, opt, n=100):
    """One realizability case; returns (minE, max(|F|-cE)) over every stage and step,
    both relative to the reference energy 1."""
    c, dx = 1.0, 1.0/n
    chi = chidx/dx
    x = (np.arange(n) + 0.5)*dx
    ka, ks = (chi, 0.0) if kind == 'abs' else (0.0, chi)
    E, F = np.zeros(n), np.zeros(n)
    if test in ('beam', 'slab'):
        if test == 'slab':
            msk = x > 0.5
            ka, ks = ka*msk, ks*msk
        pb = S.Problem(n, ka=ka, ks=ks, bc='inflow')
    else:
        pb = S.Problem(n, ka=ka, ks=ks)
        if test == 'pulse':
            m = (x > 0.2) & (x < 0.4)
            E[m], F[m] = 1.0, c
        else:
            m1, m2 = (x > 0.15) & (x < 0.3), (x > 0.7) & (x < 0.85)
            E[m1], F[m1] = 1.0, c
            E[m2], F[m2] = 1.0, -c
    worst = [0.0, -np.inf]

    def chk(U, Y1):
        for V in (U, Y1):
            if V is None:
                continue
            e, f = V[:n], V[n:2*n]
            if not np.all(np.isfinite(V)):
                worst[0], worst[1] = -np.inf, np.inf
                continue
            worst[0] = min(worst[0], e.min())
            worst[1] = max(worst[1], (np.abs(f) - c*e).max())
    h = cfl*dx/c
    tend = max(0.6, 8*h)
    S.run(pb, E, F, tend, h, opt, check=chk)
    return worst


CFLS = (0.5, 4.0, 100.0)


def pos_variant(opt):
    res = []
    for test in ('beam', 'pulse', 'counter', 'slab'):
        for chidx in (0.0, 0.01, 0.1, 1.0, 10.0, 100.0):
            for kind in (('scat',) if chidx == 0 else ('scat', 'abs')):
                for cfl in CFLS:
                    w = pos_case(test, chidx, kind, cfl, opt)
                    res.append((test, chidx, kind, cfl, w[0], w[1]))
    return opt, res


def test_pos():
    tol = 1e-12
    with Pool(24) as p:
        out = p.map(pos_variant, VARIANTS)
    for opt, res in out:
        bad = [r for r in res if r[4] < -tol or r[5] > tol]
        wE = min(r[4] for r in res)
        wF = max(r[5] for r in res)
        print('%-32s cases %d violations %d  worst minE %.2e  worst max(|F|-cE) %.2e'
              % (vname(opt), len(res), len(bad), wE, wF))
        for r in bad:
            print('    %-7s chi dx %-6g %-4s cfl %-5g minE %+.2e max(|F|-cE) %+.2e' % r)
        sys.stdout.flush()


# ------------------------------------------------------------------ test 2 / time order
def thin_run(n, chidx64, cfl, opt, tend=0.5):
    c = 1.0
    chi = chidx64*64.0
    pb = S.Problem(n, ks=0.5*chi, ka=0.5*chi, B=1.0)
    x = pb.x
    E = 1.0 + 0.3*np.sin(2*np.pi*x)
    F = 0.2*c*np.sin(2*np.pi*x + 0.5)
    U = S.run(pb, E, F, tend, cfl*pb.dx/c, opt)
    return U[:n], U[n:2*n]


def coarsen(u, n):
    return u.reshape(n, -1).mean(axis=1)


def test_thin(opts, cfl=0.5, chis=(0.01, 0.1, 1.0)):
    for opt in opts:
        for chidx in chis:
            t0 = time.time()
            Er, Fr = thin_run(2048, chidx, cfl, opt)
            errs = []
            for n in (64, 128, 256):
                E, F = thin_run(n, chidx, cfl, opt)
                errs.append(np.abs(E - coarsen(Er, n)).mean() + np.abs(F - coarsen(Fr, n)).mean())
            o = np.log2(np.array(errs[:-1])/np.array(errs[1:]))
            print('%-22s cfl %-4g chi dx(64) %-5g L1(E)+L1(F) %s  orders %s  (%.0f s)'
                  % (vname(opt), cfl, chidx, ' '.join('%.2e' % e for e in errs),
                     ' '.join('%.2f' % v for v in o), time.time() - t0))
            sys.stdout.flush()


def test_torder_thick(opts, n=128):
    for opt in opts:
        sol = []
        for nst in (8, 16, 32, 64, 128):
            c, chi = 1.0, 6400.0
            pb = S.Problem(n, ks=chi)
            x = pb.x
            kap = c/(3.0*chi)
            s0 = 0.05
            tend = 0.5*s0**2/kap
            E = np.exp(-0.5*((x - 0.5)/s0)**2)
            F = kap*(x - 0.5)/s0**2*E
            sol.append(S.run(pb, E, F, tend, tend/nst, opt)[:n])
        d = [np.abs(sol[k] - sol[k + 1]).mean() for k in range(4)]
        o = np.log2(np.array(d[:-1])/np.array(d[1:]))
        print('%-22s thick Gaussian N 128 chi dx 50, c dt/dx 384..24 successive diffs %s orders %s'
              % (vname(opt), ' '.join('%.2e' % v for v in d), ' '.join('%.2f' % v for v in o)))
        sys.stdout.flush()


def test_torder(opts):
    for opt in opts:
        for chidx in (0.1, 100.0):
            sol = [thin_run(128, chidx, cfl, opt, tend=1.0)[0] for cfl in (4, 2, 1, 0.5, 0.25)]
            d = [np.abs(sol[k] - sol[k + 1]).mean() for k in range(4)]
            o = np.log2(np.array(d[:-1])/np.array(d[1:]))
            print('%-22s npic %d N 128 chi dx %-5g cfl 4..0.25 successive diffs %s  orders %s'
                  % (vname(opt), opt.get('npic', 3), chidx, ' '.join('%.2e' % v for v in d),
                     ' '.join('%.2f' % v for v in o)))
            sys.stdout.flush()


# ------------------------------------------------------------------ test 3
def thick_run(n, opt, dtmode, chi=6400.0):
    c = 1.0
    pb = S.Problem(n, ks=chi)
    x, dx = pb.x, pb.dx
    kap = c/(3.0*chi)
    s0 = 0.05
    tend = 0.5*s0**2/kap
    nst = 16*(n//64) if dtmode == 'dx' else 16*(n//64)**2
    E = np.exp(-0.5*((x - 0.5)/s0)**2)
    F = kap*(x - 0.5)/s0**2*E
    U = S.run(pb, E, F, tend, tend/nst, opt)
    sig2 = s0**2 + 2*kap*tend
    ex = s0/np.sqrt(sig2)*np.exp(-0.5*(x - 0.5)**2/sig2)
    E = U[:n]
    var = np.sum(E*(x - 0.5)**2)/np.sum(E)
    ratio = (var - s0**2)/(2*kap*tend) - 1.0        # numerical / physical diffusion
    return np.abs(E - ex).mean(), ratio, tend/nst*c/dx


def test_thick(opts):
    for opt in opts:
        for dtmode in ('dx', 'dx2'):
            r = [thick_run(n, opt, dtmode) for n in (64, 128, 256)]
            errs = np.array([q[0] for q in r])
            o = np.log2(errs[:-1]/errs[1:])
            print('%-22s dt~%-3s chi dx 100,50,25 L1 %s orders %s  kappa_num/kappa-1 %s'
                  ' (c dt/dx %s)'
                  % (vname(opt), dtmode, ' '.join('%.2e' % e for e in errs),
                     ' '.join('%.2f' % v for v in o), ' '.join('%+.3f' % q[1] for q in r),
                     ' '.join('%.3g' % q[2] for q in r)))
            sys.stdout.flush()


def compact_diffusion(n, chi=6400.0, dtmode='dx'):
    """Reference: dE/dt = d/dx(kappa dE/dx), compact 3-point, same SDIRK2 and steps."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    c, dx = 1.0, 1.0/n
    x = (np.arange(n) + 0.5)*dx
    kap = c/(3.0*chi)
    s0 = 0.05
    tend = 0.5*s0**2/kap
    nst = 16*(n//64) if dtmode == 'dx' else 16*(n//64)**2
    h = tend/nst
    L = sp.diags([1, -2, 1], [-1, 0, 1], shape=(n, n)).tolil()
    L[0, n - 1] = L[n - 1, 0] = 1
    L = (kap/dx**2)*L.tocsc()
    A = (sp.identity(n) - S.GAM*h*L).tocsc()
    E = np.exp(-0.5*((x - 0.5)/s0)**2)
    for _ in range(nst):
        Y1 = spla.spsolve(A, E)
        k1 = (Y1 - E)/(S.GAM*h)
        E = spla.spsolve(A, E + (1 - S.GAM)*h*k1)
    return E


def test_compact():
    """Thick-limit check: alpha = 0 must reproduce the compact diffusion stencil."""
    for n in (64, 128):
        pb_opt = dict(cand='zero', recon=None, rc='one')
        c, chi = 1.0, 6400.0
        pb = S.Problem(n, ks=chi)
        x = pb.x
        kap = c/(3.0*chi)
        s0 = 0.05
        tend = 0.5*s0**2/kap
        E = np.exp(-0.5*((x - 0.5)/s0)**2)
        F = kap*(x - 0.5)/s0**2*E
        U = S.run(pb, E, F, tend, tend/(16*(n//64)), dict(pb_opt, Dfix=1.0/3.0))
        Ec = compact_diffusion(n)
        print('compact-limit check N %d: alpha=0, D=1/3: max|E - E_compact| = %.2e (max E %.2f)'
              % (n, np.abs(U[:n] - Ec).max(), Ec.max()))


# ------------------------------------------------------------------ test 4
def test_oddeven(opts, nsteps=20):
    for opt in opts:
        for cfl in (0.5, 100.0):
            n = 64
            pb = S.Problem(n, ks=100.0*n)
            sgn = (-1.0)**np.arange(n)
            E = 1.0 + 0.1*sgn
            F = np.zeros(n)
            amps = []
            U = np.r_[E, F]
            for k in range(nsteps):
                U, _ = S.step(pb, U, cfl*pb.dx, opt)
                amps.append(np.abs(np.mean(U[:n]*sgn)))
            print('%-22s odd-even chi dx 100 cfl %-5g amplitude 0.1 -> step1 %.2e step5 %.2e'
                  ' step%d %.2e  minE %.4f' % (vname(opt), cfl, amps[0], amps[4], nsteps,
                                                amps[-1], U[:n].min()))


# ------------------------------------------------------------------ test 5
def rwave_setup(n, ka=10.0, ks=0.0, rho_cv=1.0, D=1.0/3.0):
    c, k = 1.0, 2*np.pi
    gas = dict(rho_cv=rho_cv, T0=1.0, B0=1.0, Bp=4.0)
    chi = ka + ks
    A = np.array([[-c*ka, -1j*k, c*ka*4.0],
                  [-1j*k*c*c*D, -c*chi, 0.0],
                  [c*ka/rho_cv, 0.0, -c*ka*4.0/rho_cv]])
    lam, V = np.linalg.eig(A)
    cplx = np.abs(lam.imag) > 1e-6
    j = np.argmax(np.where(cplx, lam.real, -np.inf)) if cplx.any() else np.argmax(lam.real)
    pb = S.Problem(n, ka=ka, ks=ks, gas=gas)
    return pb, lam[j], V[:, j]/V[0, j], lam


def test_rwave(opts):
    for opt in opts:
      for ka in (1.0, 100.0):
        for Dmode in ('Dfix', 'rays'):
            o = dict(opt, Dfix=1.0/3.0) if Dmode == 'Dfix' else dict(opt)
            errs = []
            for n in (32, 64, 128):
                pb, lam, v, alll = rwave_setup(n, ka=ka)
                x = pb.x
                eps = 1e-3
                ph = np.exp(2j*np.pi*x)
                E = 1.0 + eps*np.real(v[0]*ph)
                F = eps*np.real(v[1]*ph)
                T = 1.0 + eps*np.real(v[2]*ph)
                tend = 1.0
                U = S.run(pb, E, F, tend, 0.5*pb.dx, o, T=T)
                Eh = np.sum((U[:n] - 1.0)*np.conj(ph))*2/n/eps     # complex amplitude
                lam_num = np.log(Eh/v[0])/tend
                lam_num += 2j*np.pi*np.round((lam.imag - lam_num.imag)*tend/(2*np.pi))/tend
                errs.append((n, lam_num, abs(lam_num - lam)))
            e = np.array([q[2] for q in errs])
            print('ka %-5g %-22s %-4s exact lambda %.5f%+.5fi | numeric %s | |err| %s orders %s'
                  % (ka, vname(opt), Dmode, lam.real, lam.imag,
                     ' '.join('%.5f%+.5fi' % (q[1].real, q[1].imag) for q in errs),
                     ' '.join('%.1e' % q for q in e),
                     ' '.join('%.2f' % q for q in np.log2(e[:-1]/e[1:]))))


# ------------------------------------------------------------------ test 6
def test_cons(opts):
    for opt in opts:
        n, c = 100, 1.0
        x = (np.arange(n) + 0.5)/n
        m = (x > 0.2) & (x < 0.4)
        E0, F0 = np.where(m, 1.0, 0.0), np.where(m, c, 0.0)
        for lab, pb, totF in (('vacuum pulse', S.Problem(n), True),
                              ('scattering pulse chi dx 1', S.Problem(n, ks=100.0), False)):
            for cfl in (0.5, 100.0):
                rec = []
                S.run(pb, E0, F0, max(0.6, 8*cfl/n), cfl/n, opt,
                      check=lambda U, Y: rec.append((U[:n].sum(), U[n:].sum())))
                dE = max(abs(r[0] - E0.sum()) for r in rec)/E0.sum()
                dF = max(abs(r[1] - F0.sum()) for r in rec)/F0.sum() if totF else np.nan
                print('%-22s %-26s cfl %-5g rel drift sum E %.1e  sum F %.1e'
                      % (vname(opt), lab, cfl, dE, dF))
        pb, lam, v, _ = rwave_setup(64)
        ph = np.exp(2j*np.pi*pb.x)
        E = 1.0 + 0.1*np.real(v[0]*ph) + 0.05
        F = 0.1*np.real(v[1]*ph)
        T = 1.0 + 0.1*np.real(v[2]*ph)
        rec = []
        S.run(pb, E, F, 1.0, 4.0/64, opt, T=T,
              check=lambda U, Y: rec.append(U[:64].sum() + U[128:].sum()))
        tot0 = E.sum() + T.sum()
        print('%-22s gas+rad (ka=10, rho_cv=1) cfl 4  rel drift sum(E + rho_cv T) %.1e'
              % (vname(opt), max(abs(r - tot0) for r in rec)/tot0))


if __name__ == '__main__':
    which = sys.argv[1]
    base = [dict(cand='A1', recon=None), dict(cand='A1', recon='minmod'),
            dict(cand='A1', recon='mc'), dict(cand='A3', recon='minmod'),
            dict(cand='A2', recon='minmod')]
    base = [dict(b, Dguard=True) for b in base]
    t0 = time.time()
    if which == 'pos':
        test_pos()
    elif which == 'thin':
        test_thin(base[int(sys.argv[2]):int(sys.argv[2]) + 1])
    elif which == 'torder':
        test_torder([dict(cand='A1', recon='minmod', Dguard=True, npic=k) for k in (2, 3, 5)]
                    + [dict(cand='A1', recon=None, Dguard=True)]
                    + [dict(cand='A1', recon='minmod', Dguard=True, ti='be')])
    elif which == 'torder2':
        test_torder_thick([dict(cand='A1', recon='mc', Dguard=True),
                           dict(cand='A1', recon='mc', Dguard=True, ti='be')])
    elif which == 'thick':
        test_compact()
        test_thick(base)
    elif which == 'oddeven':
        test_oddeven(base)
    elif which == 'rwave':
        test_rwave([base[0], base[2]])
    elif which == 'cons':
        test_cons(base[:2])
    print('elapsed %.0f s' % (time.time() - t0))
