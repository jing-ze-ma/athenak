# usage: python3 bincmp.py dirA dirB : max |a-b|/max|a| per variable over all bin files (data only)
import sys, glob, numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert as bc
worst = 0.0; nsame = 0; n = 0
for fa in sorted(glob.glob(sys.argv[1] + '/bin/*.bin')):
    fb = fa.replace(sys.argv[1], sys.argv[2], 1)
    A = bc.read_binary(fa); B = bc.read_binary(fb)
    for v in A['var_names']:
        a = np.concatenate([np.ravel(x) for x in A['mb_data'][v]])
        b = np.concatenate([np.ravel(x) for x in B['mb_data'][v]])
        n += 1
        if np.array_equal(a, b): nsame += 1; continue
        d = np.abs(a - b).max() / max(np.abs(a).max(), 1e-300)
        worst = max(worst, d); print(fa.split('/')[-1], v, f'{d:.2e}')
print(f'variables compared {n}, bitwise identical {nsame}, worst rel diff {worst:.2e}')
