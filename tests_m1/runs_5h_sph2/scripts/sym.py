#!/usr/bin/env python3
"""T-sym: the radial tab slices at several (j,k) of the last dump must agree"""
import glob
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/sph2_0924/scripts')
from ana import tab  # noqa: E402
d = sys.argv[1]
num = sorted(glob.glob(d + '/tab/*.j0k0.*.tab'))[-1].split('.')[-2]
fs = sorted(glob.glob(f'{d}/tab/*.j*k*.{num}.tab'))
es = np.array([tab(f)['m1_e'] for f in fs])
e0 = tab(fs[0].replace(f'.{num}.', '.00000.'))['m1_e']
sp = np.max((es.max(axis=0) - es.min(axis=0))/es.mean(axis=0))
msg = (f"t={tab(fs[0])['time']:.4g} slices={len(fs)} spread(E)/E={sp:.2e} "
       f"evolution={np.max(np.abs(es[0]-e0))/np.max(e0):.2e}")
hs = sorted(glob.glob(f'{d}/tab/*.hj*k*.{num}.tab'))
if hs:
    for v in ('dens', 'velx'):
        a = np.array([tab(f)[v] for f in hs])
        msg += f" {v} spread {np.max(a.max(axis=0)-a.min(axis=0))/np.max(np.abs(a)):.2e}"
print(msg)
