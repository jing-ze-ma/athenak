#!/usr/bin/env python3
"""real_nb.py RUN_DIR (CONV3 2026-11-20): is the |F|/cE >> 1 of the output cell vector an unrealisable staggered state, or the
cell average of faces next to a near-floor cell?  The staggered x1 face F between cells i-1, i is realisable iff |F| <= c E_donor
(donor = upwind cell, at most c max(E_{i-1}, E_i)).  The output m1_f1 of cell i is ~ the mean of its two x1 faces, so a bound for
the faces is |F1_i| <= c max(E_{i-1}, E_i, E_{i+1}) (and laterally for F2/F3).  Prints, for cells with |F|/cE > 10: the counts and
percentiles of |F1|/(c max_x1-neighbours E), |F|/(c max_3x3x3 E), and how many are at the floor (E <= 1.001e-4)."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.environ.get("BCDIR", os.path.dirname(os.path.abspath(__file__)))); import bin_convert as bc
C = 2.99792458e10
fs = sorted(glob.glob(os.path.join(sys.argv[1], 'a', 'bin', '*.m1.*.bin')) or glob.glob(os.path.join(sys.argv[1], 'bin', '*.m1.*.bin')))
for f in (fs[0], fs[-1]):
    d = bc.read_binary(f); md = d['mb_data']
    print('==', os.path.basename(f), 't', d['time'])
    arr = lambda v: np.stack([np.asarray(md[v][m]) for m in range(d['n_mbs'])]).astype(np.float64)
    # (cubed sphere: m1_f*/m1_phi* are the orthonormal (r, theta, phi) components of the frame of y, CleanStagHydro/
    #  CleanCsFrame, so |F| = sqrt(sum of squares) is physical and *1 is radial; no metric conversion needed)
    E = arr('m1_e')
    for tag in ('f', 'phi'):
      if 'm1_%s1' % tag not in d['var_names']: continue
      F1, F2, F3 = arr('m1_%s1' % tag), arr('m1_%s2' % tag), arr('m1_%s3' % tag)   # (mb, k, j, i)
      fm = np.sqrt(F1**2 + F2**2 + F3**2); r = fm/(C*np.abs(E) + 1e-300)
      print('  -- m1_%s' % tag)
      pad = np.pad(E, ((0, 0), (1, 1), (1, 1), (1, 1)), mode='edge')   # (MeshBlock edges: own edge, a lower bound of the bound)
      ex1 = np.maximum(np.maximum(pad[:, 1:-1, 1:-1, :-2], pad[:, 1:-1, 1:-1, 2:]), E)
      e27 = E.copy()
      for dk in (0, 1, 2):
        for dj in (0, 1, 2):
            for di in (0, 1, 2):
                e27 = np.maximum(e27, pad[:, dk:dk + E.shape[1], dj:dj + E.shape[2], di:di + E.shape[3]])
      q = lambda a: ' '.join('%.3e' % v for v in np.percentile(a, [50, 90, 99, 100])) if a.size else '-'
      for nm, sl in (('all', np.s_[...]), ('top120', np.s_[..., -120:]), ('shells300-400', np.s_[..., 300:400])):
          m = r[sl] > 10.0
          a1 = (np.abs(F1[sl])/(C*ex1[sl]))[m]; a3 = (fm[sl]/(C*e27[sl]))[m]
          print('  %s: cells %d |F|/cE>10 %d (at floor %d) | |F1|/(c max E_x1nb) p50/90/99/max %s (>1: %d) | |F|/(c max E_27) %s (>1: %d)'
                % (nm, r[sl].size, m.sum(), (E[sl][m] <= 1.001e-4).sum(), q(a1), (a1 > 1).sum(), q(a3), (a3 > 1).sum()))
