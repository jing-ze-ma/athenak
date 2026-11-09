"""raycheap2: self-convergence of E (all cells) over runs given finest-last: python3 conv.py RUN_c0.15 RUN_c0.075 ...
d_k = |E_k - E_k+1|_1/|E_finest|_1, p_k = log2(d_k/d_k+1)  (needs vis/python/bin_convert.py)"""
import sys, glob, os, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'vis', 'python'))
import bin_convert as bc
def E(run):
    f = sorted(glob.glob(f'{run}/bin/*.m1.*.bin'))[-1]; d = bc.read_binary(f)
    o = sorted(range(len(d['mb_geometry'])), key=lambda m: tuple(np.round(d['mb_geometry'][m], 12)))
    return d['time'], np.concatenate([np.asarray(d['mb_data']['m1_e'][m]).ravel() for m in o])
U = [E(r) for r in sys.argv[1:]]
n = np.abs(U[-1][1]).sum(); d = [np.abs(U[k][1] - U[k+1][1]).sum()/n for k in range(len(U) - 1)]
print('t_end', [round(u[0], 3) for u in U]); print('d', ['%.3e' % x for x in d])
print('p', ['%.2f' % math.log2(d[k]/d[k+1]) for k in range(len(d) - 1)])
