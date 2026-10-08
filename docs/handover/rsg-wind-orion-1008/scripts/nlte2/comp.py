"""Composition from pyfastchem 4.0.3 equilibrium condensation (lowP/chem/fc_eq.npz),
log-interpolated in (T, P)."""
import numpy as np
F = np.load('/orion/ptmp/jinma/rsg_wind_1008/lowP/chem/fc_eq.npz')
NM = list(F['names'])
LT, LP = np.log(F['T']), np.log(F['P'])
LN = np.log(np.maximum(F['n'], 1e-300))
LNG = np.log(F['ng'])


def _interp(arr, T, Pbar):
    lt, lp = np.log(T), np.log(Pbar)
    i = np.clip(np.searchsorted(LT, lt) - 1, 0, len(LT) - 2)
    j = np.clip(np.searchsorted(LP, lp) - 1, 0, len(LP) - 2)
    a = (lt - LT[i])/(LT[i+1] - LT[i])
    b = (lp - LP[j])/(LP[j+1] - LP[j])
    return ((1-a)*(1-b)*arr[i, j] + a*(1-b)*arr[i+1, j] + (1-a)*b*arr[i, j+1]
            + a*b*arr[i+1, j+1])


def comp(T, Pbar, species):
    """number densities (cm^-3) of the listed species at (T, P bar)."""
    lv = _interp(LN, T, Pbar)
    return {s: float(np.exp(lv[NM.index(s)])) for s in species}


def ntot(T, Pbar):
    return Pbar*1e6/(1.380649e-16*T)
