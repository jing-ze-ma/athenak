"""FastChem 4.0.3 (gas only, with ions) for the AG Car mixture on the cool part of the
grid: which molecules exist at log rho -21..-12, T 2818-5012 K; and the electron density
vs this work's own Saha.  writes work/chem.npz and prints a summary."""
import sys
import numpy as np
import pyfastchem
sys.path.insert(0, '/orion/ptmp/jinma/agcar_opac_1009/scripts')
import lowrho_lib as L

D = '/orion/ptmp/jinma/rsg_wind_1008/lowP/FastChem/input/'      # read-only
W = L.B + 'work/'
SYM = {1: 'H', 2: 'He', 6: 'C', 7: 'N', 8: 'O', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al',
       14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc',
       22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
       30: 'Zn'}
with open(W + 'agcar_abund.dat', 'w') as f:
    f.write('# AG Car X=0.36 Z=0.02 GS98 (LANL grevsau1), log eps, H=12\ne-  0.00\n')
    for z in L.ZL:
        f.write('%-3s %.4f\n' % (SYM[z], 12 + np.log10(L.NPG[z]/L.NPG[1])))
fc = pyfastchem.FastChem(W + 'agcar_abund.dat', D + 'logK/logK.dat', 0)
at = L.Atoms()
lT = np.round(np.arange(3.45, 3.951, 0.025), 3)
lr = np.round(np.arange(-21.0, -11.99, 0.05), 2)
mols = ['H2', 'C1O1', 'H2O1', 'H1O1', 'O1Ti1', 'O1Si1', 'C1N1', 'N2', 'O1V1', 'C1H1',
        'Fe1H1', 'H1Mg1', 'O2', 'C2', 'H1Si1', 'O1S1', 'S1Si1', 'N1O1', 'H-', 'e-']
res = {m: np.zeros((len(lT), len(lr))) for m in mols}
ne_s = np.zeros((len(lT), len(lr)))
mfrac = np.zeros((len(lT), len(lr)))       # fraction of all nuclei bound in molecules
for i, t in enumerate(lT):
    T = 10**t
    for j, r in enumerate(lr):
        rho = 10**r
        ne, n, _ = L.saha(at, T, rho)
        ne_s[i, j] = ne
        ntot = sum(n.values()) + ne
        P = ntot*L.kB*T/1e6
        for it in range(3):
            inp, out = pyfastchem.FastChemInput(), pyfastchem.FastChemOutput()
            inp.temperature, inp.pressure = [T], [P]
            fc.calcDensities(inp, out)
            mu = out.mean_molecular_weight[0]
            rho_fc = P*1e6/(L.kB*T)*mu*L.amu
            P *= rho/rho_fc
        nd = np.array(out.number_densities[0])
        for m in mols:
            k = fc.getGasSpeciesIndex(m)
            res[m][i, j] = nd[k]/rho if k != pyfastchem.FASTCHEM_UNKNOWN_SPECIES else np.nan
        # nuclei in molecules: total element density minus atoms+atomic ions
        nat = 0.0
        for z in L.ZL:
            for s in ['', '1+', '2+']:
                k = fc.getGasSpeciesIndex(SYM[z] + s)
                if k != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
                    nat += nd[k]
        mfrac[i, j] = 1 - nat/out.total_element_density[0]
np.savez(W + 'chem.npz', lT=lT, lr=lr, ne_saha=ne_s, mfrac=mfrac,
         **{'n_' + m.replace('-', 'm'): v for m, v in res.items()})
np.set_printoptions(precision=2, linewidth=200)
print('log rho', lr[::4])
for m in ['H2', 'C1O1', 'H2O1', 'H1O1', 'O1Ti1', 'O1Si1', 'C1N1']:
    print(m, 'log n/rho [per g] (vs N_el/g)')
    print(np.log10(res[m][:, ::4] + 1e-300))
print('fraction of nuclei in molecules (log)')
print(np.log10(mfrac[:, ::4] + 1e-300))
print('log ne FastChem/Saha')
print(np.log10(res['e-'][:, ::4]*10**lr[::4]/ne_s[:, ::4]))
print('N_C/g %.3e  N_O/g %.3e N_Ti/g %.3e' % (L.NPG[6], L.NPG[8], L.NPG[22]))
