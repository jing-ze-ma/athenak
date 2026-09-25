"""FastChem 3 (pyfastchem, gas phase) VMRs of H2 He H e- H- and mu on a (T,p) grid."""
import numpy as np
import pyfastchem
import os
D = os.environ.get('FASTCHEM_INPUT',
                   '/viper/ptmp2/jinma/ckhitemp_0925/fastchem_input').rstrip('/') + '/'


def run(T, P):
    fc = pyfastchem.FastChem(
        D+'element_abundances/asplund_2009.dat', D+'logK/logK.dat', 0)
    TT, PP = np.meshgrid(T, P, indexing='ij')
    inp = pyfastchem.FastChemInput()
    out = pyfastchem.FastChemOutput()
    inp.temperature = TT.ravel()
    inp.pressure = PP.ravel()
    fc.calcDensities(inp, out)
    n = np.array(out.number_densities)
    ngas = PP.ravel()*1e6/(1.380649e-16*TT.ravel())
    idx = [fc.getGasSpeciesIndex(s) for s in ('H2', 'He', 'H', 'e-', 'H1-')]
    res = np.empty((TT.size, 6))
    res[:, 0] = np.array(out.mean_molecular_weight)
    for k, i in enumerate(idx):
        res[:, k+1] = n[:, i]/ngas
    return res.reshape(len(T), len(P), 6), np.array(out.fastchem_flag)


def read_old(fn=None):
    from ck_lib import DATA
    fn = fn or DATA+'CE_tables/FastChem_ck_1x_int.txt'
    tok = open(fn).read().split()
    nT, nP, nrec, nsp = map(int, tok[:4])
    p = 4+nsp
    T = np.array(tok[p:p+nT], float)
    p += nT
    P = np.array(tok[p:p+nP], float)
    p += nP
    C = np.array(tok[p:p+nT*nP*6], float).reshape(nT, nP, 6)
    return T, P, C


if __name__ == '__main__':
    T, P, C = read_old()
    sel = T >= 2000
    R, flag = run(T[sel], P)
    print('flags nonzero', np.count_nonzero(flag))
    names = ['mu', 'H2', 'He', 'H', 'e-', 'H-']
    for it in [np.argmin(abs(T[sel]-t)) for t in (2000, 3000, 4000, 5000, 6100)]:
        for ip in (0, 18, 24, 30, 33):
            o = C[sel][it, ip]
            r = R[it, ip]
            print(f"T={T[sel][it]:.0f} p={P[ip]:.1e} " +
                  ' '.join(f"{n}:{o[k]:.3e}/{r[k]:.3e}" for k, n in enumerate(names)))
