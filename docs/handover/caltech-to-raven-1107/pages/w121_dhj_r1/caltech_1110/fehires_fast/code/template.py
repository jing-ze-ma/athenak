"""template.py: Fe I CCF template with pRT 3.4 lbl (R 1e6): isothermal 3000 K, solar Fe
(all neutral), H2/He background (Rayleigh + CIA), g = 955 cm/s^2, R_p 1.742 R_J.
Template = transit depth with Fe minus the same without Fe, on the fe_core ln-lambda grid,
broadened by a Gaussian of FWHM 4.28 km/s (Seidel+2026 4UT MR).
-> template.npz (lam grid parameters, raw and broadened template); fe_hires reads only 'raw'
(shipped as data/template_raw.npz) and applies hipass + the planet LSF at run time
(fe_core.template_hp).  The template is planet-independent: the CCF divides by sum(template).
Needs petitRADTRANS 3.4 and its lbl Fe file (README).
  python template.py OUT.npz"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fe_core as F  # noqa: E402
from petitRADTRANS.radtrans import Radtrans  # noqa: E402

RJ = 7.1492e9
p = np.logspace(-10, 2, 121)
T = np.full_like(p, 3000.0)
xfe = F.A_FE_SUN*55.845*F.XH/1.00794
mf = {'H2': np.full_like(p, F.XH), 'He': np.full_like(p, F.YHE), 'Fe': np.full_like(p, xfe)}
out = {}
rt = Radtrans(pressures=p, line_species=['Fe'], rayleigh_species=['H2', 'He'],
              gas_continuum_contributors=['H2--H2', 'H2--He'], line_opacity_mode='lbl',
              wavelength_boundaries=[0.3795, 0.7905])
for k, f in (('fe', 1.0), ('cont', 1e-30)):
    m = dict(mf)
    m['Fe'] = mf['Fe']*f
    w, r, _ = rt.calculate_transit_radii(temperatures=T, mass_fractions=m,
                                         mean_molar_masses=np.full_like(p, 2.33),
                                         reference_gravity=955.0, planet_radius=1.742*RJ,
                                         reference_pressure=0.01)
    out[k] = (w, r)
w = out['fe'][0]
dep = (out['fe'][1]**2 - out['cont'][1]**2)/(1.461*6.957e10)**2
N = int(np.log(F.LAM1/F.LAM0)/F.DLN) + 1
lg = np.log(F.LAM0) + F.DLN*np.arange(N)
o = np.argsort(w)
tmp = np.interp(lg, np.log(w[o]), dep[o])
sig = 4.28e5/2.3548/F.DV
k = np.arange(-int(5*sig), int(5*sig) + 1)
ker = np.exp(-0.5*(k/sig)**2)
ker /= ker.sum()
tb = np.convolve(tmp, ker, mode='same')
np.savez(sys.argv[1], lnl0=np.log(F.LAM0), dln=F.DLN, raw=tmp, tb=tb,
         wprt=w, dep=dep)
print('N', N, 'max depth', tmp.max(), 'prt points', len(w), 'dlnl', np.median(np.diff(np.log(w[o]))))
