import sys
from petitRADTRANS.radtrans import Radtrans
import numpy as np
sp = sys.argv[1:]
rt = Radtrans(pressures=np.logspace(-6,2,10), line_species=sp, rayleigh_species=['H2','He'],
              gas_continuum_contributors=['H2-H2','H2-He','H-'], wavelength_boundaries=[0.5,5.5],
              line_opacity_mode='c-k')
print('OK', sp, rt.frequencies.size, flush=True)
