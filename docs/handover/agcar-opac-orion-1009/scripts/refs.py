"""Reference tables: TOPS (LANL, X0.36 Z0.02 GS98, AthenaK table_rho format) and Ferguson
et al. 2005 GS98 (Wichita f05.g98.pl / f05.gs98, X=0.35 and 0.5 at Z=0.02, linearly
interpolated in X to 0.36).  log10 kappa [cm^2/g] at (log T, log rho)."""
import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI

TOPS = '/orion/ptmp/jinma/agcar_1008/files/'          # read-only
FERG = '/orion/ptmp/jinma/agcar_opac_1009/dl/ferg/'


def read_repo(fn):
    hdr, vals = None, []
    for ln in open(fn):
        if ln.startswith('#'):
            p = ln[1:].split()
            if hdr is None and len(p) == 6:
                try:
                    hdr = (int(p[0]), int(p[1]), *[float(x) for x in p[2:]])
                except ValueError:
                    pass
            continue
        if ln.strip():
            vals.append(float(ln))
    nT, nD, lt0, dlt, ld0, dld = hdr
    return lt0 + dlt*np.arange(nT), ld0 + dld*np.arange(nD), np.array(vals).reshape(nT, nD)


def read_ferg(fn):
    ln = open(fn).read().splitlines()
    i0 = [i for i, s in enumerate(ln) if s.lstrip().startswith('log T')][0]
    lR = np.array(ln[i0].split()[2:], float)
    rows = [s.split() for s in ln[i0+1:] if s.strip()]
    lT = np.array([r[0] for r in rows], float)
    K = np.array([r[1:] for r in rows], float)
    o = np.argsort(lT)
    return lT[o], lR, K[o]


class Tops:
    VALID_T, VALID_D = (3.764, 7.065), (-14.0, 0.0)

    def __init__(self, kind):
        t, d, K = read_repo(TOPS + f'{kind}_tops_gs98_x0.36_z0.02.txt')
        self.f = RGI((t, d), K, bounds_error=False, fill_value=np.nan)

    def __call__(self, lT, lD):
        lT, lD = np.broadcast_arrays(lT, lD)
        v = self.f(np.stack([lT, lD], -1))
        ok = (lT >= self.VALID_T[0]) & (lD >= self.VALID_D[0])
        return np.where(ok, v, np.nan)


class Ferg:
    def __init__(self, kind, X=0.36):
        nm = {'planck': 'g98.pl.%s.02.tpon', 'rosseland': 'g98.%s.02.tron'}[kind]
        t1, r1, K1 = read_ferg(FERG + nm % '35')
        t2, r2, K2 = read_ferg(FERG + nm % '5')
        assert np.allclose(t1, t2) and np.allclose(r1, r2)
        K = K1 + (X - 0.35)/(0.5 - 0.35)*(K2 - K1)
        self.f = RGI((t1, r1), K, bounds_error=False, fill_value=np.nan)
        self.lR = r1

    def __call__(self, lT, lD):
        lT, lD = np.broadcast_arrays(lT, lD)
        lR = lD - 3*lT + 18
        return self.f(np.stack([lT, lR], -1))
