"""order_1009: L1|arm - cen|/L1|cen - <cen>| per level (E for pulses, rho for rw)"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import od

def rel(a, b):
    return np.abs(a - b).sum()/np.abs(b - b.mean()).sum()

for case in ('pulse_k0.128', 'pulse_k12.8', 'pulse_k128', 'pulse_k1280', 'pulse_k12800',
             'rw_t10', 'rw_t1000'):
    var, key = ('hydro_w', 'dens') if case.startswith('rw') else ('m1', 'm1_e')
    arms = ['hr', 'hrb']
    for arm in arms:
        for sch in ('hesdirk2', 'be'):
            for mode in 'XCT':
                levs = od.TL if mode == 'T' else od.NX
                if case == 'rw_t1000':
                    levs = [2, 4, 8, 16] if mode == 'T' else [32, 64, 128, 256]
                out = []
                for n in levs:
                    try:
                        a = od.cols(od.rundir(case, arm, sch, mode, n), var)[key]
                        b = od.cols(od.rundir(case, 'cen', sch, mode, n), var)[key]
                        out.append('%.2e' % rel(a, b))
                    except Exception:
                        out.append('   -    ')
                print('%-13s %-4s %-8s %s |%s-cen| %s' % (case, arm, sch, mode, arm, ' '.join(out)))
