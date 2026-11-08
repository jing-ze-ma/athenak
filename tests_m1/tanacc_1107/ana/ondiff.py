"""max and mean |X_on - X_off| / J0(r) (E) and r^2 |F_on - F_off| over every cell, per pair; + D_rr check."""
import sys, glob, os
import numpy as np
from evalrun import load
R = os.environ.get('TA_RUNS', '/viper/ptmp2/jinma/openitems_1107/tanacc/runs')
for on in sorted(glob.glob(R + '/*_on*')):
    off = on.replace('_on', '_off')
    if not os.path.exists(off + '/bin'):
        continue
    a, b = load(on, 'm1'), load(off, 'm1')
    r = load(on, 'm1_lat')['r']
    J0 = 2*(1 - np.sqrt(1 - 1/r**2))
    dE = np.abs(a['m1_e'] - b['m1_e'])/J0
    dF = np.abs(a['m1_f1'] - b['m1_f1'])*r**2
    print('%-22s dE/J0 max %.2e mean %.2e | r2 dFr max %.2e mean %.2e | bitwise %s' % (
        os.path.basename(on).replace('_on', ''), dE.max(), dE.mean(), dF.max(), dF.mean(),
        bool((a['m1_e'] == b['m1_e']).all())))
