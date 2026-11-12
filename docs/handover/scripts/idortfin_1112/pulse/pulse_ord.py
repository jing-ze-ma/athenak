"""idortfin-1112 pulse_cs ladder orders (as st17.evaluate on cs):
T ladder (dt / l, fixed grid): d(l, 2l) = mean|Q_l - Q_2l| / mean|Q_2l - mean Q_2l| per field over every cell; p = log2(d(1,2)/d(2,4)) ...
C ladder (N 32/64/128): printed as L1 of the shell-mean E, F1 differences is NOT done here (st17 restricts in r); only the T orders.
usage: python3 pulse_ord.py <root with B_T_1 .. FIN_T_8 run dirs>   (ATHENAK_VIS = athenak vis/python for bin_convert)"""
import math, os, sys
import numpy as np
sys.path.insert(0, os.environ.get('ATHENAK_VIS', '/orion/ptmp/jinma/Athenak/wt_fin/vis/python'))
import bin_convert


def binflat(d):
    out = {}
    for fn in sorted(os.listdir(os.path.join(d, 'bin'))):
        if fn.endswith('.bin'):
            b = bin_convert.read_binary(os.path.join(d, 'bin', fn))
            for v in b['var_names']:
                out[v] = np.concatenate([np.ravel(x) for x in b['mb_data'][v]])
    return out


R = sys.argv[1]
for arm in ('B', 'FIN'):
    ls = [1, 2, 4, 8]
    Q = {l: binflat(os.path.join(R, '%s_T_%d' % (arm, l))) for l in ls}
    for k in ('m1_e', 'm1_f1', 'm1_f2', 'm1_f3', 'dens', 'eint', 'velx'):
        if k not in Q[1]:
            continue
        d = [np.abs(Q[a][k] - Q[b][k]).mean() / np.abs(Q[b][k] - Q[b][k].mean()).mean() for a, b in zip(ls[:-1], ls[1:])]
        p = [math.log2(d[i]/d[i+1]) for i in range(len(d) - 1)]
        print('%-4s %-6s d %s  p %s' % (arm, k, ' '.join('%.3e' % x for x in d), ' '.join('%.2f' % x for x in p)))
