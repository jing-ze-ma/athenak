"""chem.py: FastChem 3 equilibrium chemistry (gas + local equilibrium condensation, ions) on the
(T, p) grid of the emission opacity table. Asplund 2009 1x solar (= run eos_xh 0.7381, eos_yhe 0.2485).
-> tables/chem_1x.npz : T, p[bar], mass fractions per pRT species, MMW"""
import numpy as np, pyfastchem, sys
D = '/resnick/groups/carnegie_poc/jingze/obscheck_1102/ext/fastchem/input/'
T = np.concatenate([np.arange(300, 1500, 100), np.arange(1500, 6101, 200)]).astype(float)
P = np.logspace(-8, 3, 34)
SP = {'H2O': 'H2O1', 'CO': 'C1O1', 'CO2': 'C1O2', 'CH4': 'C1H4', 'OH': 'H1O1', 'HCN': 'C1H1N1_hcn',
      'TiO': 'O1Ti1', 'VO': 'O1V1', 'FeH': 'Fe1H1', 'Na': 'Na', 'K': 'K', 'SiO': 'O1Si1',
      'H2': 'H2', 'He': 'He', 'H': 'H', 'H-': 'H1-', 'e-': 'e-', 'Fe': 'Fe', 'Ti': 'Ti'}
MW = {'H2O': 18.015, 'CO': 28.01, 'CO2': 44.01, 'CH4': 16.04, 'OH': 17.007, 'HCN': 27.025, 'TiO': 63.866,
      'VO': 66.94, 'FeH': 56.853, 'Na': 22.99, 'K': 39.098, 'SiO': 44.085, 'H2': 2.016, 'He': 4.0026,
      'H': 1.008, 'H-': 1.008, 'e-': 5.4858e-4, 'Fe': 55.845, 'Ti': 47.867}
cond = len(sys.argv) < 2 or sys.argv[1] != 'nocond'
fc = pyfastchem.FastChem(D + 'element_abundances/asplund_2009.dat', D + 'logK/logK.dat',
                         D + 'logK/logK_condensates.dat', 0)
TT, PP = np.meshgrid(T, P, indexing='ij')
inp = pyfastchem.FastChemInput(); out = pyfastchem.FastChemOutput()
inp.temperature = TT.ravel(); inp.pressure = PP.ravel()
inp.equilibrium_condensation = cond; inp.rainout_condensation = False
flag = fc.calcDensities(inp, out)
print('FastChem:', pyfastchem.FASTCHEM_MSG[flag], 'cond', cond)
n = np.array(out.number_densities); mmw = np.array(out.mean_molecular_weight)
ntot = n.sum(1)
res = {}
for k, s in SP.items():
    i = fc.getGasSpeciesIndex(s)
    if i == pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
        print('missing', k, s); continue
    res[k] = (n[:, i]/ntot*MW[k]/mmw).reshape(TT.shape)
np.savez('/resnick/groups/carnegie_poc/jingze/obscheck_1102/tables/chem_1x%s.npz' % ('' if cond else '_nocond'),
         T=T, p=P, mmw=mmw.reshape(TT.shape), **{'X_' + k: v for k, v in res.items()})
for k in ('H2O', 'CO', 'TiO', 'VO', 'H-', 'e-'):
    print(k, ' '.join('%.1e' % res[k][j, 25] for j in range(0, len(T), 4)))
