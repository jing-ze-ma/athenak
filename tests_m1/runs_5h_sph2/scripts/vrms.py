#!/usr/bin/env python3
"""rms and max |v| and max |T/T0-1| of the last hydro_w bin dump (T = p/rho)"""
import glob
import sys
import numpy as np
import binlib
fs = sorted(glob.glob(sys.argv[1] + '/bin/*.bin'))
h0, h = binlib.load(fs[0]), binlib.load(fs[-1])
v2 = h['velx']**2 + h['vely']**2 + h['velz']**2
k = 'eint' if 'eint' in h else 'press'
t = h[k]/h['dens']
t0 = h0[k]/h0['dens']
print(f"cycles {h0['cycle']}..{h['cycle']} rms|v|={np.sqrt(v2.mean()):.3e} "
      f"max|v|={np.sqrt(v2.max()):.3e} max|T/T0-1|={np.max(np.abs(t/t0-1)):.2e}")
