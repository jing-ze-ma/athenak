"""make_ovr_1x.py: the 1x-solar sensitivity deck overrides = make_ovr.py output (planet keys
from the 1x setup + 1x grid) plus the metallicity keys of the WASP-121b 1x production input
(EOS composition, met, rad_met) and the 1x ck tables.
usage: make_ovr_1x.py SETUP_1x.txt GRIDKEYS_1x.txt OUT.json"""
import json
import os
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
subprocess.run([sys.executable, os.path.join(here, 'make_ovr.py')] + sys.argv[1:4], check=True)
o = json.load(open(sys.argv[3]))
c = '1x SOLAR sensitivity arm (README sec. 3d); value of the WASP-121b 1x input'
o['comment/problem'][0] = o['comment/problem'][0].replace('PURE HYDRO', 'PURE HYDRO, 1x solar')
o['hydro/eos_xh'] = ['0.7381', c]
o['hydro/eos_yhe'] = ['0.2485', c]
o['hydro/rad_met'] = ['0.0', c]
o['problem/met'] = ['0.0', c]
o['problem/ck_table'] = ['ck/CK1X/ck/Premixed_1x_g8_11_hiT2.txt', '1x premixed ck table (repo data/exo_fms_ck)']
o['problem/ck_ce_table'] = ['CE_tables/FastChem_ck_1x_int_hiT2.txt', '1x FastChem CE table']
o['problem/ck_data_dir'] = ['ck/CK1X', '1x table set']
o['problem/ic_profile'] = ['ic/ic_toi1518b_1x_f5.txt', 'face-5 exact RCE IC, 1x solar']
json.dump(o, open(sys.argv[3], 'w'), indent=1)
print('1x keys added:', sys.argv[3], len(o))
