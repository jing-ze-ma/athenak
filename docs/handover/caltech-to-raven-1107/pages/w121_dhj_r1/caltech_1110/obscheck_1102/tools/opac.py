"""opac.py: premixed correlated-k emission table from petitRADTRANS 3.4 (its own opacity code:
species k-tables interpolated to (T, p), mass-fraction weighted, combined by random overlap with
resort-rebin (pRT fcore.combine_ck_opacities), plus CIA H2-H2 / H2-He and H- bf/ff continuum) on
the FastChem grid of chem.py; then rebinned R1000 -> R1000/NREB by resort-rebin of NREB adjacent
bins (16 g each -> 16 g).  Absorption only (no Rayleigh, no clouds).
  python opac.py [chem_npz] [out_npz] [NREB]
"""
import sys
import numpy as np
from petitRADTRANS.radtrans import Radtrans

O = '/resnick/groups/carnegie_poc/jingze/obscheck_1102/'
chemf = sys.argv[1] if len(sys.argv) > 1 else O + 'tables/chem_1x.npz'
outf = sys.argv[2] if len(sys.argv) > 2 else O + 'tables/opac_R200.npz'
NREB = int(sys.argv[3]) if len(sys.argv) > 3 else 5
LINES = {'H2O': 'H2O__POKAZATEL', 'CO': 'CO-NatAbund__HITEMP', 'CO2': 'CO2__UCL-4000', 'CH4': 'CH4__HITEMP',
         'OH': 'OH__MoLLIST', 'HCN': 'HCN__Harris', 'TiO': 'TiO-NatAbund__McKemmish', 'VO': 'VO__VOMYT',
         'FeH': 'FeH__MoLLIST', 'Na': 'Na__Allard', 'K': 'K__Allard', 'SiO': 'SiO__SiOUVenIR'}
ch = np.load(chemf)
T, P = ch['T'], ch['p']
rt = Radtrans(pressures=P, line_species=list(LINES.values()), rayleigh_species=[],
              gas_continuum_contributors=['H2--H2', 'H2--He', 'H-'], wavelength_boundaries=[0.5, 5.6],
              line_opacity_mode='c-k')
nu = rt.frequencies
lam = 2.99792458e10/nu*1e4                       # um, pRT orders by frequency
g = rt._lines_loaded_opacities['g_gauss']
w = rt._lines_loaded_opacities['weights_gauss']
print('bins', nu.size, 'g', g.size, 'lam', lam.min(), lam.max(), flush=True)
kap = np.zeros((len(T), len(P), nu.size, g.size), np.float64)
for i, t in enumerate(T):
    mf = {LINES[k]: np.maximum(ch['X_' + k][i], 1e-40) for k in LINES}
    for k in ('H2', 'He', 'H', 'H-', 'e-'):
        mf[k] = np.maximum(ch['X_' + k][i], 1e-40)
    res = rt._calculate_opacities(temperatures=np.full(P.size, t), mass_fractions=mf,
                                  mean_molar_masses=ch['mmw'][i], reference_gravity=1e3)
    op = res[0] if isinstance(res, tuple) else res   # (g, freq, species, layers); total in species 0
    kap[i] = np.transpose(op[:, :, 0, :], (2, 1, 0))
    print('T', t, 'kappa(1 bar) median %.3e' % np.median(kap[i, np.argmin(abs(P - 1))]), flush=True)
# order by wavelength, rebin NREB adjacent bins by resort-rebin
o = np.argsort(lam)
lam = lam[o]
kap = kap[:, :, o]
edges_nu = rt.frequency_bins_edges
le = np.sort(2.99792458e10/edges_nu*1e4)
nb = lam.size//NREB
cum = np.concatenate([[0], np.cumsum(w)])
kr = np.zeros((len(T), len(P), nb, g.size), np.float32)
lamr = np.zeros(nb); ler = np.zeros((nb, 2))
for j in range(nb):
    sl = slice(j*NREB, (j + 1)*NREB)
    # equal spectral weight per sub-bin (they span ~equal d ln lambda)
    kk = kap[:, :, sl].reshape(len(T), len(P), -1)          # NREB*NG samples
    ww = np.tile(w, NREB)/NREB
    srt = np.argsort(kk, axis=-1)
    ks = np.take_along_axis(kk, srt, -1)
    cw = np.cumsum(ww[srt], -1)
    # value at the target g quadrature points: interpolate the sorted k at cumulative weight g
    for n, gg in enumerate(g):
        idx = (cw < gg).sum(-1)
        idx = np.clip(idx, 0, kk.shape[-1] - 1)
        kr[:, :, j, n] = np.take_along_axis(ks, idx[..., None], -1)[..., 0]
    lamr[j] = np.exp(np.log(lam[sl]).mean())
    ler[j] = le[j*NREB], le[min((j + 1)*NREB, le.size - 1)]
np.savez(outf, T=T, p_bar=P, lam_um=lamr, lam_edges_um=ler, wg=w, g=g, lnk=np.log(np.maximum(kr, 1e-60)),
         species=str(LINES), chem=chemf)
print('wrote', outf, kr.shape)
