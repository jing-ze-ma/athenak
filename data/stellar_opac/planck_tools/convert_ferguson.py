#!/usr/bin/env python3
"""Ferguson et al. (2005) low-temperature GS98 tables -> the AthenaK table
format, optionally SPLICED onto a TOPS table at high T.

Ferguson files (https://www.math.wichita.edu/~ferguson/):
  f05.g98.pl.tar.gz  -> g98.pl.<X><Z>.tpon   PLANCK  means, log10 kappa
  f05.gs98.tar.gz    -> g98.<X><Z>.tron      ROSSELAND means, log10 kappa
Both are log10 kappa [cm^2/g] on rows of log T and columns of
log R = log10 rho - 3 log10 T + 18.  Planck rows: log T 2.70..4.50 step 0.05
(55).  Rosseland rows: same range, finer near the top (85).  The X and Z in a
file name are the decimal digits: g98.pl.7.004 = X=0.7, Z=0.004,
g98.pl.0.02 = X=0.0, Z=0.02.  Z grid: 0, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 2e-3,
4e-3, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1 -- there is NO Z=0.014, so
that composition is interpolated linearly in log10 Z between 0.01 and 0.02.

The repo tables are on (log T, log rho), read with rad_kappa_src = table_rho,
so no EOS is needed: log R is converted with log rho = log R + 3 log T - 18.

usage:
  convert_ferguson.py FERG_FILE[,FERG_FILE2,Z1,Z2,Ztarget] OUT.txt [TOPS.txt]
If a third argument is given it is an AthenaK-format table (e.g. the TOPS
Planck table) that supplies log T >= 4.2; the two are blended linearly in
log T over 4.0..4.2, exactly as tools/stellar_opac/merge_rosseland.py blends
AESOPUS onto OPLIB for the Rosseland table.
"""
import sys
import numpy as np
from convert_tops import LT0, DLT, NT, LD0, DLD, ND, bilin, write_table
from compare import read_ferg, read_repo

args = sys.argv[1].split(',')
out = sys.argv[2]
if len(args) == 1:
    lT_s, lR_s, F = read_ferg(args[0])
    src = args[0].split('/')[-1]
else:
    f1, f2, z1, z2, zt = args[0], args[1], *[float(x) for x in args[2:5]]
    lT_s, lR_s, F1 = read_ferg(f1)
    lT2, lR2, F2 = read_ferg(f2)
    assert np.allclose(lT_s, lT2) and np.allclose(lR_s, lR2)
    w = (np.log10(zt) - np.log10(z1))/(np.log10(z2) - np.log10(z1))
    F = (1-w)*F1 + w*F2
    src = '%s and %s, interpolated in log Z to Z=%g (w=%.3f)' % (
        f1.split('/')[-1], f2.split('/')[-1], zt, w)

lT = LT0 + DLT*np.arange(NT)
lD = LD0 + DLD*np.arange(ND)
TT, DD = np.meshgrid(lT, lD, indexing='ij')
RR = DD - 3*TT + 18.0
K = bilin(lT_s, lR_s, F, TT, RR)
notes = ['source: Ferguson et al. (2005) low-temperature table %s' % src,
         'VALID: log T %.2f..%.2f, log R %.1f..%.1f (logR = log rho - 3 log T '
         '+ 18); elsewhere EDGE-FILLED, not data.'
         % (lT_s[0], lT_s[-1], lR_s[0], lR_s[-1])]

if len(sys.argv) > 3:
    hT, hD, hK = read_repo(sys.argv[3])
    KH = bilin(hT, hD, hK, TT, DD)
    w = np.clip((TT - 4.0)/0.2, 0.0, 1.0)
    K = (1-w)*K + w*KH
    notes.append('spliced onto %s for log T >= 4.2, linear blend in log T over '
                 '4.0..4.2 (same recipe as tools/stellar_opac/merge_rosseland.py)'
                 % sys.argv[3].split('/')[-1])

write_table(out, K, ['AthenaK stellar opacity table, log10 kappa [cm^2/g] on '
                     '(log10 T[K], log10 rho[g/cm3])'] + notes +
            ['Read with <hydro>/rad_kappa_src = table_rho.'])
print('wrote', out)
