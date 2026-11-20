#!/usr/bin/env python3
"""real_diag.py RUN_DIR: sanity of real.py's |F_r|/cE on one m1 bin (first and last): variables, block layout,
radial range of the 'top 120 shells', percentiles of E, |F1|, |F1|/cE, |F|/cE (all 3 comps), E<=0 count; deep shells too."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.environ.get("BCDIR", os.path.dirname(os.path.abspath(__file__)))); import bin_convert as bc
C = 2.99792458e10; NT = 120
fs = sorted(glob.glob(os.path.join(sys.argv[1], 'a', 'bin', '*.m1.*.bin')) or glob.glob(os.path.join(sys.argv[1], 'bin', '*.m1.*.bin')))
for f in (fs[0], fs[-1]):
    d = bc.read_binary(f); md = d['mb_data']
    print('==', os.path.basename(f), 't', d['time'], 'vars', d['var_names'], 'Nx', d['Nx1'], d['Nx2'], d['Nx3'], 'n_mbs', d['n_mbs'],
          'mb', d['nx1_mb'], d['nx2_mb'], d['nx3_mb'], 'out', d['nx1_out_mb'], d['nx2_out_mb'], d['nx3_out_mb'])
    x1 = np.asarray(d['mb_x1v']); print('  mb x1v range first MB', x1[0][0], x1[0][-1], ' distinct radial ranges', len({(r[0], r[-1]) for r in x1}))
    arr = lambda v: np.stack([np.asarray(md[v][m]) for m in range(d['n_mbs'])])
    E, F1 = arr('m1_e'), arr('m1_f1'); F2, F3 = arr('m1_f2'), arr('m1_f3')
    print('  array shape', E.shape, 'dtype', E.dtype)
    for nm, sl in (('top120', np.s_[..., -NT:]), ('shells300-400', np.s_[..., 300:400])):
        e, f1 = E[sl], F1[sl]; fm = np.sqrt(F1[sl]**2 + F2[sl]**2 + F3[sl]**2)
        q = lambda a: ' '.join('%.3e' % v for v in np.percentile(a, [1, 50, 99, 100]))
        print('  %s: E<=0 %d  E p1/50/99/max %s | |F1| %s | |F1|/cE %s | |F|/cE %s' % (nm, (e <= 0).sum(), q(e), q(np.abs(f1)),
              q(np.abs(f1)/(C*np.abs(e) + 1e-300)), q(fm/(C*np.abs(e) + 1e-300))))
