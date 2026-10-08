"""Premix DACE 0.01 cm^-1 cross sections (cm^2/g of species) with FastChem VMRs and
sort to the Exo-FMS 8 g-points.  k per gram of gas = sum_i X_i m_i kappa_i / mu."""
import glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/orion/u/jinma/ATHENAK/athenak/docs/handover/rsg-ck-1008/scripts')
import ck_lib as ck  # noqa
from species import SP  # noqa
R = '/orion/ptmp/jinma/rsg_wind_1008/dace/raw/'
KT = ck.read_ktable('/orion/u/jinma/ATHENAK/athenak/docs/handover/rsg-ck-1008/data/ck/'
                    'Premixed_1x_g8_11_hiT2.txt')
WN = KT['wn']                            # ascending band edges, cm^-1 (= table band order)
NB = len(WN) - 1
DNU = 0.01
NMAX = int(np.ceil(WN[-1]/DNU)) + 2
I0 = np.round(WN/DNU).astype(int)        # grid index of band edges (band = [I0[b], I0[b+1]))
GX, GW = KT['gx'], KT['gw']
GE = np.concatenate([[0.], np.cumsum(GW)])
GE[-1] = 1.
FCN = {'H2O': 'H2O1', 'CH4': 'C1H4', 'C2H2': 'C2H2', 'C2H4': 'C2H4', 'CO': 'C1O1',
       'CO2': 'C1O2', 'NH3': 'H3N1', 'H2S': 'H2S1', 'HCl': 'Cl1H1', 'HCN': 'C1H1N1_hcn',
       'HF': 'F1H1', 'PH3': 'H3P1', 'SiO': 'O1Si1', 'TiO': 'O1Ti1', 'VO': 'O1V1',
       'MgH': 'H1Mg1', 'CaH': 'Ca1H1', 'TiH': 'H1Ti1', 'CrH': 'Cr1H1', 'FeH': 'Fe1H1',
       'OH': 'H1O1', 'SH': 'H1S1', 'SiH': 'H1Si1', 'SiN': 'N1Si1', 'SiS': 'S1Si1',
       'C2': 'C2', 'CH': 'C1H1', 'CN': 'C1N1', 'Fe': 'Fe', 'FeII': 'Fe1+', 'Na': 'Na',
       'K': 'K', 'Ti': 'Ti', 'TiII': 'Ti1+', 'Cr': 'Cr', 'Ca': 'Ca', 'CaII': 'Ca1+',
       'Mg': 'Mg', 'MgII': 'Mg1+', 'Si': 'Si', 'SiII': 'Si1+', 'Al': 'Al', 'V': 'V',
       'Mn': 'Mn', 'Ni': 'Ni'}
SPD = {s[0]: s for s in SP}


def avail_T(name):
    fs = glob.glob(f'{R}{name}/Out_*_n800.bin')
    return np.array(sorted(int(os.path.basename(f).split('_')[3]) for f in fs))


def load(name, T, below='hold'):
    """kappa (cm^2/g) on 0..NMAX at node T; above T_max: held at T_max; below T_min:
    'hold' at T_min or 'zero'. returns (array, T_used)"""
    Ta = avail_T(name)
    if T > Ta[-1]:
        Tu = Ta[-1]
    elif T < Ta[0]:
        if below == 'zero':
            return None, None
        Tu = Ta[0]
    else:
        Tu = Ta[np.argmin(abs(Ta-T))]
        assert Tu == T, (name, T, Tu)
    f = glob.glob(f'{R}{name}/Out_*_{int(Tu):05d}_n800.bin')[0]
    a = np.fromfile(f, '<f4', count=NMAX)
    out = np.zeros(NMAX, np.float32)
    out[:len(a)] = a
    return out, Tu


def sortk(kmix):
    """per band: (value at Gauss g, mean over g sub-interval, band mean)"""
    A = np.zeros((NB, len(GX)))
    B = np.zeros((NB, len(GX)))
    M = np.zeros(NB)
    for b in range(NB):
        s = np.sort(kmix[I0[b]:I0[b+1]].astype(np.float64))
        n = len(s)
        g = (np.arange(n)+0.5)/n
        A[b] = np.interp(GX, g, s)
        c = np.concatenate([[0.], np.cumsum(s)])/n
        e = np.interp(GE, np.arange(n+1)/n, c)
        B[b] = np.diff(e)/GW
        M[b] = s.mean()
    return A, B, M


def vmr(fc, it, ip, names):
    z = fc
    nm = list(z['names'])
    ng = z['ng'][it, ip]
    return {s: z['n'][it, ip, nm.index(FCN[s])]/ng for s in names}, z['mu'][it, ip]


def defloor(a, w=10000):
    """subtract the smooth floor (linear interp of 100 cm^-1 block minima), clip >= 0"""
    n = len(a)//w*w
    mn = a[:n].reshape(-1, w).min(1)
    x = (np.arange(len(mn))+0.5)*w
    fl = np.interp(np.arange(len(a)), x, mn)
    return np.clip(a-fl, 0, None)
