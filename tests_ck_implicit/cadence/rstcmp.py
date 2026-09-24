"""rstcmp.py <root>: straight (s8/bin) vs restarted (r4/bin) hydro_w dumps, matched by
time; max relative difference of density and internal energy per dump."""
import glob
import sys

import numpy as np

sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert  # noqa: E402

R = sys.argv[1]
for arm in ['e1', 'e4', 'e4g2']:
    S = {}
    for f in glob.glob(R + '/' + arm + '/s8/bin/*.bin'):
        S[bin_convert.read_binary(f)['time']] = f
    out = []
    for f in sorted(glob.glob(R + '/' + arm + '/r4/bin/*.bin')):
        B = bin_convert.read_binary(f)
        if B['time'] not in S:
            out.append('t=%.4f unmatched' % B['time'])
            continue
        A = bin_convert.read_binary(S[B['time']])
        d = [np.max(np.abs(np.asarray(B['mb_data'][v]) / np.asarray(A['mb_data'][v]) - 1))
             for v in ('dens', 'eint')]
        out.append('cyc %d d %.1e e %.1e' % (B['cycle'], d[0], d[1]))
    print(arm, '; '.join(out))
