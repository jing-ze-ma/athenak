#!/usr/bin/env python3
"""Overlap of the low-T sources with TOPS where BOTH are data, per log T band."""
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/scripts')
import build_ext as b
f = np.load(sys.argv[1])
lT, lD = f['lT'], f['lD']
info = f['tops_info']
floor = np.interp(lT, info[:, 0], info[:, 1])
TT, DD = np.meshgrid(lT, lD, indexing='ij')
RR = DD - 3*TT + 18
for nm, kt, kl, rmax in (('Rosseland AESOPUS-TOPS', f['KR_t'], f['KR_l'], 6.0),
                         ('Planck Ferguson-TOPS', f['KP_t'], f['KP_l'], 1.0)):
    print('==', nm, '(dex; both data: log R in [-8,%g], log rho >= TOPS floor, log rho <= 0)' % rmax)
    for t0 in np.arange(3.775, 4.5, 0.1):
        m = (TT >= t0) & (TT < t0 + 0.1) & (RR >= -8) & (RR <= rmax) & (DD >= floor[:, None])
        d = (kl - kt)[m]
        # restrict to the AG Car-relevant density band too
        m2 = m & (DD <= -8)
        d2 = (kl - kt)[m2]
        print('log T %.2f-%.2f n %4d  median %+.3f  |d| 90%% %.3f max %.3f || rho<=1e-8: median %+.3f 90%% %.3f max %.3f' % (
            t0, t0+0.1, m.sum(), np.median(d), np.percentile(abs(d), 90), abs(d).max(),
            np.median(d2), np.percentile(abs(d2), 90), abs(d2).max()))
