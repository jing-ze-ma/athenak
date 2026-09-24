# usage: python ts.py run...   max/min E and max |F|/E per dump
import glob, sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/vetdo_0924/src/vis/python')
import bin_convert as bc
for r in sys.argv[1:]:
    print('==', r)
    out = []
    for f in sorted(glob.glob(r + '/bin/*.bin')):
        d = bc.read_binary_as_athdf(f)
        e = d['m1_e'][0]; fx = d['m1_f1'][0]; fy = d['m1_f2'][0]
        ff = np.sqrt(fx**2 + fy**2) / np.maximum(e, 1e-300)
        out.append('t=%5.2f maxE=%.3e minE=%.3e max f=%.3f' % (d['Time'], e.max(), e.min(), ff.max()))
    print('\n'.join(out))
