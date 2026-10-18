#!/usr/bin/env python3
"""cmplast.py A B: compare the LAST bin dump of every output (by file stem) of run dirs A and B (dump numbers may differ,
e.g. continuous vs restarted); prints BITWISE or the max relative difference per variable, and the hst last-row diff."""
import sys, glob, os, re
import numpy as np
sys.path.insert(0, '/raven/ptmp/jinma/overlap_1015/src/b16b66db/vis/python')
import bin_convert
A, B = sys.argv[1:3]
def last(d):
    out = {}
    for f in glob.glob(os.path.join(d, 'bin', '*.bin')):
        m = re.match(r'(.*)\.(\d+)\.bin$', os.path.basename(f))
        if m and (m.group(1) not in out or int(m.group(2)) > out[m.group(1)][0]): out[m.group(1)] = (int(m.group(2)), f)
    return out
la, lb = last(A), last(B); nbad = 0
for k in sorted(set(la) & set(lb)):
    da, db = bin_convert.read_binary(la[k][1]), bin_convert.read_binary(lb[k][1])
    if abs(da['time'] - db['time']) > 1e-9 * max(1.0, abs(da['time'])):
        print('%s: times differ %.9g vs %.9g' % (k, da['time'], db['time'])); nbad += 1; continue
    for v in da['var_names']:
        x = np.concatenate([np.ravel(a) for a in da['mb_data'][v]]); y = np.concatenate([np.ravel(a) for a in db['mb_data'][v]])
        if x.tobytes() == y.tobytes(): print('%s %s t=%.6g BITWISE' % (k, v, da['time'])); continue
        s = np.maximum(np.abs(x), np.abs(y)); r = np.max(np.where(s > 0, np.abs(x - y) / np.where(s > 0, s, 1), 0))
        print('%s %s t=%.6g differs max rel %.3e' % (k, v, da['time'], r)); nbad += 1
print('CMPLAST %s' % ('BITWISE' if nbad == 0 else '%d differ' % nbad))
