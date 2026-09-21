import sys
import glob
import numpy as np
sys.path.insert(0, '../../../vis/python')
import bin_convert as bc  # noqa: E402

c = 2.99792458e10
for f in sorted(glob.glob('bin/*.m1.*.bin')):
    d = bc.read_binary_as_athdf(f)
    names = [k for k in d.keys() if k.startswith('m1') or k in ('E', 'F1', 'F2')]
    E = d[names[0]][0]
    F1 = d[names[1]][0]
    F2 = d[names[2]][0]
    x1 = d['x1v']
    print(f, names, 't=%g' % d['Time'])
    f1 = F1/(c*E)
    f2 = np.abs(F2)/(c*E)
    for i in list(range(len(x1)-12, len(x1))) + [len(x1)-20, len(x1)-30, 40, 10]:
        print('  i=%3d z=%.4e  Emin/Emean %.3e  <f1> %.3f  max f2 %.3e  max|f| %.3f'
              % (i, x1[i], E[:, i].min()/E[:, i].mean(), f1[:, i].mean(),
                 f2[:, i].max(), np.sqrt(f1[:, i]**2 + f2[:, i]**2).max()))
