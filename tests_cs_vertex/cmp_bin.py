import sys, glob, os
sys.path.insert(0, os.path.abspath("../vis/python"))
import bin_convert
import numpy as np
def cmp(fa, fb):
    a = bin_convert.read_binary(fa); b = bin_convert.read_binary(fb)
    worst = 0.0; names=[]
    for va, vb in zip(a['mb_data'].values(), b['mb_data'].values()):
        d = np.max(np.abs(np.asarray(va)-np.asarray(vb)))
        worst = max(worst, d)
    return worst
for fa in sorted(glob.glob(sys.argv[1])):
    fb = fa.replace(sys.argv[2], sys.argv[3])
    print(os.path.basename(fa), "max|diff| = %.6e" % cmp(fa, fb))
