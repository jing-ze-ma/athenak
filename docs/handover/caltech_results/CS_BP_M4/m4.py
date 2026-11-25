#!/usr/bin/env python3
"""zzzzd3: Bp dumps, m4(r) at every radial cell below 100 Rsun, with the merge factor.  m4.py phase_script dumps..."""
import sys, importlib.util, numpy as np
spec = importlib.util.spec_from_file_location('imp', sys.argv[1]); imp = importlib.util.module_from_spec(spec); spec.loader.exec_module(imp)
RS = 6.957e10; RC = (2.0629e12, 3.4307e12, 5.7211e12); F = (8, 4, 2)
for p in sys.argv[2:]:
    cyc, t, r, a4, a8, rs, ra, c4, c8, rp = imp.analyse(p)
    R = r * RS; fac = np.where(R < RC[0], 8, np.where(R < RC[1], 4, np.where(R < RC[2], 2, 1)))
    print('\n# %s cycle %d; merge rc/Rsun: %.2f (f 8), %.2f (f 4), %.2f (f 2), then f 1' % (p.split('Bp/')[-1], cyc, RC[0] / RS, RC[1] / RS, RC[2] / RS))
    print('#  region                 cells  max m4 (r of max)  sym_rms max   m4 at the first 3 cells outside')
    for lab, q in (('f=8 (inside 2.063e12)', R < RC[0]), ('f=4', (R >= RC[0]) & (R < RC[1])), ('f=2 (to 5.72e12)', (R >= RC[1]) & (R < RC[2])), ('f=1 to 100 Rsun', (R >= RC[2]) & (r < 100))):
        if q.any():
            k = np.where(q)[0][np.argmax(a4[q])]; print('  %-22s %4d   %.3e (%.1f Rsun)   %.3e' % (lab, q.sum(), a4[k], r[k], rs[q].max()))
    print('# every cell, r < 100 Rsun: r[Rsun] f m4 phi0[deg] sym_rms')
    for k in np.where(r < 100)[0]:
        print('%7.3f %d %.3e %6.1f %.3e' % (r[k], fac[k], a4[k], np.degrees(-np.angle(c4[k]) / 4) % 90, rs[k]))
