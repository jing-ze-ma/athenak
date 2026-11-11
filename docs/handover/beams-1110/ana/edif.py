#!/usr/bin/env python3
"""max/L1 relative E difference of arm dumps vs arm A. usage: edif.py DECKDIR FILE"""
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/blend_1022/wt/vis/python')
import bin_convert as bc
d, f = sys.argv[1], sys.argv[2]
ref = bc.read_binary_as_athdf('%s/A/rep1/%s' % (d, f))['m1_e'][0]
for a in sys.argv[3:] if len(sys.argv) > 3 else ('A/rep2', 'B/rep1', 'C/rep1', 'D/rep1'):
    try:
        e = bc.read_binary_as_athdf('%s/%s/%s' % (d, a, f))['m1_e'][0]
    except Exception as x:
        print(a, 'ERR', x); continue
    s = ref > 1e-6*ref.max()
    print('%s %s: L1 rel %.3e  max rel (E>1e-6 max) %.3e' % (d.split('/')[-1], a, np.abs(e-ref).sum()/np.abs(ref).sum(),
          (np.abs(e-ref)[s]/ref[s]).max()))
