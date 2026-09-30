"""chem.py: FastChem 3 (pyfastchem 3.1.3) equilibrium tables on a (log T, log p) grid.

  python chem.py            -> chem_1x.npz, chem_10x.npz

Species stored (number densities / n_gas, log10): Fe (Fe I), Na (Na I), e-, H1- (H-), H, H2, He,
plus the mean molecular weight.  Two variants per metallicity: gas-only ('gas') and local
equilibrium condensation ('cond', no rainout).  10x = all elements heavier than He x10
relative to H (He/H unchanged, as in the run: Y/X = 0.3367 in both inputs).
Data: FastChem input files (Stock+2018, Stock+2022, Kitzmann+2024) from the FastChem GitHub.
"""
import numpy as np
import pyfastchem

D = '/viper/ptmp2/jinma/w121prod_0929/synth2_trans/data/'
SP = ['Fe', 'Na', 'e-', 'H1-', 'H', 'H2', 'He', 'Fe1+', 'Na1+']
LT = np.arange(2.30, 3.901, 0.01)       # 200 .. 7943 K
LP = np.arange(-13.0, 3.01, 0.1)        # bar


def table(met, cond):
    fc = pyfastchem.FastChem(D + 'asplund_2009.dat', D + 'logK.dat',
                             D + 'logK_condensates.dat', 0)
    ab = np.array(fc.getElementAbundances())
    for i in range(fc.getElementNumber()):
        if fc.getElementSymbol(i) not in ('H', 'He', 'e-'):
            ab[i] *= met
    fc.setElementAbundances(ab)
    fc.setParameter('accuracyChem', 1e-5)
    idx = [fc.getGasSpeciesIndex(s) for s in SP]
    assert all(i != pyfastchem.FASTCHEM_UNKNOWN_SPECIES for i in idx), idx
    TT, PP = np.meshgrid(10**LT, 10**LP, indexing='ij')
    inp = pyfastchem.FastChemInput()
    out = pyfastchem.FastChemOutput()
    inp.temperature = TT.ravel()
    inp.pressure = PP.ravel()
    if cond:
        inp.equilibrium_condensation = True
        inp.rainout_condensation = False
    st = fc.calcDensities(inp, out)
    nd = np.array(out.number_densities)
    ntot = np.array(out.total_element_density)*0 + PP.ravel()*1e6/(1.380649e-16*TT.ravel())
    res = {s: np.log10(np.maximum(nd[:, i]/ntot, 1e-99)).reshape(TT.shape)
           for s, i in zip(SP, idx)}
    res['mu'] = np.array(out.mean_molecular_weight).reshape(TT.shape)
    res['status'] = np.array(out.fastchem_flag).reshape(TT.shape)
    print('met', met, 'cond', cond, 'status', st, 'nonconv', (res['status'] != 0).sum(), flush=True)
    return res


if __name__ == '__main__':
    for tag, met in (('1x', 1.0), ('10x', 10.0)):
        d = {}
        for cond in (False, True):
            r = table(met, cond)
            for k, v in r.items():
                d[('cond_' if cond else 'gas_') + k] = v
        np.savez_compressed('chem_%s.npz' % tag, LT=LT, LP=LP, **d)
