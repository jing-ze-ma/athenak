"""FastChem carrier groups for the line taper: VMR/mu summed over the neutral species
that carry the line opacity of each band."""
import numpy as np
import pyfastchem
import os
D = os.environ.get('FASTCHEM_INPUT',
                   '/viper/ptmp2/jinma/ckhitemp_0925/fastchem_input').rstrip('/') + '/'
IR = ['H2O1', 'H1O1', 'C1O1', 'Fe1H1', 'O1Ti1',
      'O1V1', 'H1Si1', 'Ca1H1', 'H1Mg1', 'O1Si1']
OPT = ['Fe', 'Ti', 'Cr', 'Mg', 'Si', 'Ca', 'Na', 'K', 'O1Ti1', 'O1V1', 'Fe1H1']
_fc = None


def group_vmr(T, P, nb=11, nIR=8):
    global _fc
    if _fc is None:
        _fc = pyfastchem.FastChem(
            D+'element_abundances/asplund_2009.dat', D+'logK/logK.dat', 0)
    T = np.atleast_1d(np.asarray(T, float))
    P = np.broadcast_to(P, T.shape).astype(float)
    inp = pyfastchem.FastChemInput()
    out = pyfastchem.FastChemOutput()
    inp.temperature = T
    inp.pressure = P
    _fc.calcDensities(inp, out)
    n = np.array(out.number_densities)
    ng = P*1e6/(1.380649e-16*T)
    mu = np.array(out.mean_molecular_weight)
    def g(L): return sum(n[:, _fc.getGasSpeciesIndex(s)] for s in L)/ng/mu
    xi, xo = g(IR), g(OPT)
    return np.array([xi if b < nIR else xo for b in range(nb)]).T
