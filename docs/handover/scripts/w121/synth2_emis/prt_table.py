"""prt_table.py: premixed pRT correlated-k (R1000) absorption opacity on a (T, p) grid, with
FastChem equilibrium chemistry (gas + local equilibrium condensation, no rainout; ions incl. e-,
H-; thermal dissociation of H2 / H2O etc. included) at the run's metallicity.

  python prt_table.py 1x|10x eq|nodiss
eq     : chemistry at the cell T
nodiss : chemistry at min(T, 2000 K) (ablation: no thermal dissociation / ionisation above
         2000 K; opacities still at the real T)
Output tables/prt_<met>_<chem>.npz: lT, lP(bar), lam_edges(um), lk (nT, nP, nfreq, ng) float32
log10 cm^2/g (random-overlap combined lines + CIA + H- bf/ff, pRT's own mixing; Rayleigh
scattering kept separately in ks (nT, nP, nfreq), NOT used in the emission RT), gw (ng), mu.
Species line lists (pRT keeper, R1000): H2O POKAZATEL, CO HITEMP, CO2 UCL-4000, CH4 HITEMP,
OH MoLLIST, TiO McKemmish (Toto), VO VOMYT, FeH MoLLIST, Fe Kurucz, Na Allard, K Allard,
SiO SiOUVenIR; CIA H2-H2, H2-He (BoRi); H- (pRT built-in, Gray / John); Rayleigh H2, He.
pRT clamps T and p to each table's range (T max 2995-4000 K, p 1e-6/1e-5 .. 100/1000 bar).
"""
import os
import sys
import numpy as np
os.environ['pRT_emcee_mode'] = 'True'          # skip mpi4py (no libmpi on the login node)
import pyfastchem   # noqa: E402
from petitRADTRANS.radtrans import Radtrans   # noqa: E402
from petitRADTRANS.config import petitradtrans_config_parser as cp   # noqa: E402

SD = '/viper/ptmp2/jinma/w121prod_0929/synth2_emis/'
FC = {'1x': '/viper/ptmp2/jinma/ckhitemp_0925/fastchem_input/',
      '10x': '/viper/ptmp2/jinma/wasp121_0925/met10/fc_in/'}   # Asplund09; 10x = metals +1 dex
# pRT line species name -> FastChem symbol, molar mass
LS = [('1H2-16O__POKAZATEL', 'H2O1', 18.0106), ('12C-16O__HITEMP', 'C1O1', 27.9949),
      ('12C-16O2__UCL-4000', 'C1O2', 43.9898), ('12C-1H4__HITEMP', 'C1H4', 16.0313),
      ('16O-1H__MoLLIST', 'H1O1', 17.0027), ('48Ti-16O__McKemmish', 'O1Ti1', 63.9429),
      ('51V-16O__VOMYT', 'O1V1', 66.9389), ('56Fe-1H__MoLLIST', 'Fe1H1', 56.9427),
      ('56Fe__Kurucz', 'Fe', 55.9349), ('23Na__Allard', 'Na', 22.9898),
      ('39K__Allard', 'K', 38.9637), ('28Si-16O__SiOUVenIR', 'O1Si1', 43.9719)]
CONT = [('H2', 'H2', 2.01565), ('He', 'He', 4.00260), ('H', 'H', 1.00783),
        ('H-', 'H1-', 1.00838), ('e-', 'e-', 5.4858e-4)]
LT = np.linspace(np.log10(150.0), np.log10(6000.0), 48)
LP = np.linspace(-12.0, 3.0, 46)
WL = [0.3, 28.0]


def chem(met, T, P):
    """T, P(bar) 1-D -> dict of mass fractions per pRT key, mean molar mass."""
    d = FC[met]
    fc = pyfastchem.FastChem(d + 'element_abundances/asplund_2009.dat', d + 'logK/logK.dat',
                             d + 'logK/logK_condensates.dat', 0)
    inp = pyfastchem.FastChemInput()
    out = pyfastchem.FastChemOutput()
    inp.temperature = np.asarray(T, float)
    inp.pressure = np.asarray(P, float)
    inp.equilibrium_condensation = True
    inp.rainout_condensation = False
    fc.calcDensities(inp, out)
    n = np.array(out.number_densities)
    mu = np.array(out.mean_molecular_weight)
    ntot = np.asarray(P)*1e6/(1.380649e-16*np.asarray(T))
    mf = {}
    for key, sym, m in LS + CONT:
        mf[key] = n[:, fc.getGasSpeciesIndex(sym)]/ntot*m/mu
    return mf, mu, np.array(out.fastchem_flag)


if __name__ == '__main__':
    met, mode = sys.argv[1], sys.argv[2]
    cp.set_input_data_path('/viper/ptmp2/jinma/prt_data')
    P = 10**LP
    rt = Radtrans(pressures=P, wavelength_boundaries=WL, line_species=[s for s, _, _ in LS],
                  gas_continuum_contributors=['H2--H2', 'H2--He', 'H-'],
                  rayleigh_species=['H2', 'He'], line_opacity_mode='c-k',
                  scattering_in_emission=True)
    lk = None
    for it, T in enumerate(10**LT):
        Tc = min(T, 2000.0) if mode == 'nodiss' else T
        mf, mu, flag = chem(met, np.full(P.size, Tc), P)
        _, _, add = rt.calculate_flux(temperatures=np.full(P.size, T), mass_fractions=mf,
                                      mean_molar_masses=mu, reference_gravity=1000.0,
                                      return_opacities=True)
        op = add['opacities']                 # (g, freq, species, layer); slot 0 = total
        ks = add['continuum_opacities_scattering']
        if lk is None:
            ng, nf = op.shape[0], op.shape[1]
            lk = np.zeros((len(LT), len(LP), nf, ng), np.float32)
            KS = np.zeros((len(LT), len(LP), nf), np.float32)
            MU = np.zeros((len(LT), len(LP)))
        lk[it] = np.log10(np.maximum(op[:, :, 0, :], 1e-60)).transpose(2, 1, 0)
        KS[it] = np.asarray(ks).reshape(nf, -1)[:, :len(LP)].T
        MU[it] = mu
        print('T %.0f  flags %d  mu %.3f..%.3f  k(1 bar, band mean) %.3e' % (
            T, np.count_nonzero(flag), mu.min(), mu.max(),
            np.mean(10**lk[it, np.argmin(abs(LP)), :, :])), flush=True)
    c = 2.99792458e10
    fe = rt._frequency_bins_edges if hasattr(rt, '_frequency_bins_edges') else rt.frequency_bins_edges
    lam_e = c/np.asarray(fe)*1e4
    np.savez(SD + 'tables/prt_%s_%s.npz' % (met, mode), lT=LT, lP=LP, lam_edges=lam_e, lk=lk,
             ks=KS, gw=np.asarray(rt._lines_loaded_opacities['weights_gauss']), mu=MU,
             lam=c/np.asarray(rt._frequencies)*1e4)
    print('wrote', lk.shape)
