# amlzf-1112: time-order ladder from the final m1 dumps (so.m1.00001.bin) of <dir>_T_{1,2,4,8,16}
# usage: python3 order.py <base dir prefix, e.g. /viper/ptmp2/jinma/pulsecs_fin/stiff/FIN>
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/amcons_1112/wt/vis/python')
import bin_convert as bc
pre = sys.argv[1]
levs = [1, 2, 4, 8, 16]
dat = {}
for l in levs:
    try:
        d = bc.read_binary('%s_T_%d/bin/so.m1.00001.bin' % (pre, l))
        dat[l] = d
    except Exception as e:
        print('lev', l, 'missing', e)
if not dat:
    sys.exit(0)
names = [n for n in dat[min(dat)]['var_names']]
for v in names:
    errs = []
    for a, b in zip(levs[:-1], levs[1:]):
        if a in dat and b in dat:
            ea = np.concatenate([np.ravel(x) for x in dat[a]['mb_data'][v]])
            eb = np.concatenate([np.ravel(x) for x in dat[b]['mb_data'][v]])
            errs.append(np.mean(np.abs(ea - eb)))
    ords = ['%.2f' % np.log2(errs[k]/errs[k+1]) if errs[k+1] > 0 else 'inf' for k in range(len(errs)-1)]
    print('%-6s dL1(T1-2,2-4,4-8,8-16) %s  orders %s' % (v, ' '.join('%.3e' % e for e in errs), ' '.join(ords)))
