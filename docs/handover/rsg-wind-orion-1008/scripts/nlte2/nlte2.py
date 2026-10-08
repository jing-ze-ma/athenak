"""Stage 1: atomic non-LTE heating on the (star, r, P, T) grid for a set of collision /
ionisation / velocity configurations.   python nlte2.py atoms  -> atoms_res.pkl
Stage 2 (combine with molecules, bands, continuum, roots): see thermal.py."""
import os
import sys
import pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
import numpy as np
from multiprocessing import Pool
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
sys.path.insert(0, HERE)
import atoms as A   # noqa: E402
import comp as CP   # noqa: E402

TG = np.geomspace(150., 3500., 64)
PBAR = (1e-11, 1e-10, 1e-9)
XR = (2., 3., 4., 5.)
RSUN = 6.957e10
STARS = {'golden16': dict(Teff=4106., R=669.875*RSUN),
         'm20lgl5.5': dict(Teff=3776., R=None)}
NEUT = ['Fe', 'Ca', 'Ti', 'Cr', 'Mn', 'Ni', 'Mg', 'Na', 'K', 'Al', 'Si']
IONS = ['Fe', 'Ca', 'Ti', 'Mg', 'Si', 'Cr']
EXTRA = ['O', 'C']              # neutral, not photoionised by the star
EMAX = {0: 6.0, 1: 7.5}

# configurations: (SH, qforb, ion mode, v km/s, uvfac)
CONFIGS = [(sh, qf, 'pi', 10., 1.) for sh in (0.01, 0.1, 1.) for qf in (1e-12, 1e-11, 1e-10)]
CONFIGS += [(0.1, 1e-11, 'lte', 10., 1.), (0.1, 1e-11, 'pi', 30., 1.),
            (0.1, 1e-11, 'pi', 10., 0.01), (1., 1e-10, 'lte', 10., 1.)]
QFS = 1e-9


def dilution(x):
    return 0.5*(1 - np.sqrt(max(1 - 1/x**2, 0.)))


def alpha_rr(T):
    """radiative recombination X+ + e -> X, rough common value (x3)."""
    return 3e-12*(T/1000.)**-0.7


_AT = {}


def atom(el, ion):
    k = (el, ion)
    if k not in _AT:
        a = A.Atom(el, ion, 11.0 if el in EXTRA else EMAX[ion])
        if el == 'Fe' and ion == 0:
            a.add_barklem()
        _AT[k] = a
    return _AT[k]


def star_R(name):
    if STARS[name]['R'] is not None:
        return STARS[name]['R']
    sys.path.insert(0, '/orion/ptmp/jinma/rsg_wind_1008/nlte')
    os.environ.setdefault('CK_DATA', '/orion/ptmp/jinma/rsg_wind_1008/lowP/data/')
    import rsglib as L
    return L.star(name)['R']


def species():
    s = ['H', 'H2', 'He', 'e-']
    for el in set(NEUT + IONS + EXTRA):
        s += [el, el + '1+']
    return s


def point(T, P, W, Teff, dvdr, cfg):
    SH, qf, imode, v, uvf = cfg
    c = CP.comp(T, P, species())
    nH, nH2 = c['H'], c['H2']
    Ntot = {el: c[el] + c[el + '1+'] for el in set(NEUT + IONS + EXTRA)}
    ne0 = max(c['e-'], 1e-10)
    out = dict(T=T, nH=nH, nH2=nH2, ne_lte=c['e-'])
    if imode == 'lte':
        ne = ne0
        x0 = {el: c[el]/max(Ntot[el], 1e-300) for el in Ntot}
        gam = {}
    else:
        gam = {}
        for el in NEUT:
            a = atom(el, 0)
            s = a.solve(T, nH, nH2, ne0, max(c[el], 1e-30), W, Teff, dvdr, SH=SH, qforb=qf,
                        qfs=QFS, uvfac=uvf)
            gam[el] = a.gamma_pi(s['n'], W, Teff, uvfac=uvf)
        al = alpha_rr(T)
        tot = sum(Ntot[el] for el in NEUT)

        def f(lne):
            ne = np.exp(lne)
            return ne - sum(Ntot[el]*gam[el]/(gam[el] + ne*al) for el in NEUT) - 0.0
        lo, hi = np.log(1e-12), np.log(max(tot, 1e-10)*1.01)
        if f(lo) > 0:
            ne = np.exp(lo)
        else:
            for _ in range(80):
                mid = 0.5*(lo + hi)
                if f(mid) > 0:
                    hi = mid
                else:
                    lo = mid
            ne = np.exp(0.5*(lo + hi))
        x0 = {el: ne*al/(gam[el] + ne*al) for el in NEUT}
        for el in EXTRA:
            x0[el] = c[el]/max(Ntot[el], 1e-300)
    out['ne'] = ne
    out['x0'] = x0
    out['gam'] = gam
    hv = {}
    diag = {}
    for el in NEUT + EXTRA:
        n0 = Ntot[el]*x0[el]
        if n0 < 1e-20:
            hv[el + ' I'] = 0.
            continue
        a = atom(el, 0)
        s = a.solve(T, nH, nH2, ne, n0, W, Teff, dvdr, SH=SH, qforb=qf, qfs=QFS, uvfac=uvf)
        hv[el + ' I'] = n0*s['heat']
        if el in ('Fe', 'Ca', 'O'):
            diag[el + ' I'] = (s['heat'], a.heat_lte(T, W, Teff, uvfac=uvf),
                               a.heat_lte(T, W, Teff, beta=s['beta'], uvfac=uvf),
                               s['absorb'], n0)
    for el in IONS:
        n1 = Ntot[el]*(1 - x0[el])
        if n1 < 1e-20:
            hv[el + ' II'] = 0.
            continue
        a = atom(el, 1)
        s = a.solve(T, nH, nH2, ne, n1, W, Teff, dvdr, SH=SH, qforb=qf, qfs=QFS, uvfac=uvf)
        hv[el + ' II'] = n1*s['heat']
        if el in ('Fe', 'Ca'):
            diag[el + ' II'] = (s['heat'], a.heat_lte(T, W, Teff, uvfac=uvf),
                                a.heat_lte(T, W, Teff, beta=s['beta'], uvfac=uvf),
                                s['absorb'], n1)
    out['hv'] = hv
    out['diag'] = diag
    return out


def task(args):
    sname, x, P, ic = args
    Teff = STARS[sname]['Teff']
    R = star_R(sname)
    W = dilution(x)
    cfg = CONFIGS[ic]
    dvdr = cfg[3]*1e5/(x*R)
    res = [point(T, P, W, Teff, dvdr, cfg) for T in TG]
    return dict(star=sname, x=x, P=P, ic=ic, cfg=cfg, W=W, Teff=Teff, dvdr=dvdr, R=R,
                pts=res)


def run_atoms():
    jobs = [(s, x, P, ic) for s in STARS for x in XR for P in PBAR
            for ic in range(len(CONFIGS))]
    with Pool(int(os.environ.get('NPROC', 24))) as p:
        res = p.map(task, jobs, chunksize=1)
    with open(HERE + 'atoms_res.pkl', 'wb') as fh:
        pickle.dump(dict(TG=TG, CONFIGS=CONFIGS, res=res), fh)


if __name__ == '__main__':
    if sys.argv[1] == 'atoms':
        run_atoms()
    elif sys.argv[1] == 'test':
        r = task(('golden16', 3., 1e-10, 1))
        for p in r['pts'][::6]:
            print('%.0f ne %.2e (lte %.2e) xFe %.1e  ' % (p['T'], p['ne'], p['ne_lte'],
                                                       p['x0']['Fe']),
                  ' '.join('%s:%.1e' % (k, v) for k, v in p['hv'].items() if abs(v) > 0))
