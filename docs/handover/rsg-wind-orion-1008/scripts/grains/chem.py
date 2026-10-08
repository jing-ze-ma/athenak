"""Gas-phase equilibrium (FastChem 4.0.3, Asplund 2009, JANAF-based logK, no condensation)
on a (T, P) grid -> forsterite activity a0 = K_c(T) p_Mg^2 p_Si p_O^4 (atoms, bar; ln K_c
from logK_condensates.dat Mg2SiO4(s,l), JANAF 1998), Mg/SiO/H2O/H2 fractions.
Cross-check: FastChem equilibrium condensation, first T where Mg-silicates appear.
Writes out/chem.npz. usage: venv/bin/python chem.py"""
import numpy as np
import pyfastchem

D = '/orion/ptmp/jinma/rsg_wind_1008/lowP/FastChem/input/'     # read-only
O = '/orion/ptmp/jinma/rsg_wind_1008/grains/out/'
kB = 1.380649e-16


def lnK_forsterite(T):
    a = [4.6852445155461133e+05, -1.2607653994717314e+01, -4.5339086543552966e+01,
         1.4741562650506954e-02, -2.0198470902129289e-06]   # solid branch, T < 2171 K
    return a[0] / T + a[1] * np.log(T) + a[2] + a[3] * T + a[4] * T**2


def run(T, P, cond=False):
    fc = pyfastchem.FastChem(D + 'element_abundances/asplund_2009.dat', D + 'logK/logK.dat',
                             D + 'logK/logK_condensates.dat', 0)
    inp, out = pyfastchem.FastChemInput(), pyfastchem.FastChemOutput()
    inp.temperature, inp.pressure = T, P
    inp.equilibrium_condensation = cond
    fc.calcDensities(inp, out)
    return fc, out


if __name__ == '__main__':
    Tg = np.arange(300., 3001., 10.)
    lPg = np.arange(-16., -1.9, 0.5)            # bar
    TT, PP = np.meshgrid(Tg, 10**lPg, indexing='ij')
    fc, out = run(TT.ravel(), PP.ravel())
    nd = np.array(out.number_densities)        # [pt, species]
    idx = {s: fc.getGasSpeciesIndex(s) for s in ['Mg', 'Si', 'O', 'H', 'H2', 'O1Si1',
                                                  'H2O1', 'Mg1+', 'e-', 'C1O1']}
    for s, i in idx.items():
        if i == pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
            print('unknown', s)
    nt = PP.ravel() * 1e6 / (kB * TT.ravel())
    p = {s: nd[:, i] * kB * TT.ravel() / 1e6 for s, i in idx.items()}   # bar
    lna = (lnK_forsterite(TT.ravel()) + 2 * np.log(p['Mg']) + np.log(p['Si'])
           + 4 * np.log(p['O']))
    eH = 1.0
    ne = np.array(out.total_element_density)
    # element abundances per H (Asplund 2009 as read by FastChem)
    eps = np.array([fc.getElementAbundance(i) for i in range(fc.getElementNumber())])
    iH = fc.getElementIndex('H')
    nH = ne * eps[iH] / eps.sum()
    sh = TT.shape
    res = dict(T=Tg, lP=lPg, lna=lna.reshape(sh), nH=nH.reshape(sh), ntot=nt.reshape(sh),
               mu=np.array(out.mean_molecular_weight).reshape(sh))
    for s in idx:
        res['x_' + s.replace('+', 'p').replace('-', 'm')] = (nd[:, idx[s]] / nH).reshape(sh)
    names = [fc.getElementSymbol(i) for i in range(fc.getElementNumber())]
    res['eps_Mg'] = eps[names.index('Mg')] / eps[iH]
    res['eps_Si'] = eps[names.index('Si')] / eps[iH]
    res['eps_O'] = eps[names.index('O')] / eps[iH]
    res['eps_C'] = eps[names.index('C')] / eps[iH]
    from_mass = {'H': 1.008, 'He': 4.0026}
    np.savez(O + 'chem.npz', **res)
    print('eps Mg Si O C', res['eps_Mg'], res['eps_Si'], res['eps_O'], res['eps_C'])
    # condensation T of forsterite (a0 = 1) at several P, gas at T
    for lp in [-4, -8, -10, -11, -12, -13]:
        j = np.argmin(abs(lPg - lp))
        la = res['lna'][:, j]
        k = np.where(np.diff(np.sign(la)))[0]
        print('P=1e%d bar  T(a=1) =' % lp, Tg[k] if len(k) else None,
              ' x_Mg, x_SiO, x_H2O, x_H2 at 1000 K:',
              ['%.2e' % res[q][70, j] for q in ['x_Mg', 'x_O1Si1', 'x_H2O1', 'x_H2']])
    # FastChem equilibrium condensation check at 1e-4 bar
    Tc = np.arange(1000., 2000., 5.)
    fc2, out2 = run(Tc, np.full(len(Tc), 1e-4), cond=True)
    ncond = np.array(out2.number_densities_cond)
    cn = [fc2.getCondSpeciesSymbol(i) for i in range(fc2.getCondSpeciesNumber())]
    for s in ['Mg2SiO4(s,l)', 'MgSiO3(s,l)', 'MgAl2O4(s)', 'Al2O3(s,l)']:
        if s in cn:
            v = ncond[:, cn.index(s)]
            nz = np.where(v > 0)[0]
            print('FastChem eq-cond 1e-4 bar: %s first at T <= %s' %
                  (s, Tc[nz.max()] if len(nz) else None))
