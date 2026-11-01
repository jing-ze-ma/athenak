"""make_ovr.py: the TOI-1518b key overrides (JSON for make_deck.py) from the final
planet_setup.py output (setup.txt) and the grid keys (grid112.py gridkeys_<N>.txt).
usage: make_ovr.py SETUP.txt GRIDKEYS.txt OUT.json"""
import json
import re
import sys

P_ROT = 164385.6          # s, 2 pi / omega, P = 1.90261131 d (Simonnin+2025)
NROT = 200                # tlim in rotations (README section 9)
st = open(sys.argv[1]).read()
m = re.search(r'<problem> keys: ap = (\S+)\s+grav = (\S+)\s+Teq = (\S+)\s+omega = (\S+)', st)
ap, grav, teq, om = m.groups()
kv = dict(t.split('=') for t in open(sys.argv[2]).read().split())
src = 'grid112.py (W121 hfresh cells/H at equal pressure; README sec. 7)'
o = {}
o['comment/problem'] = ['TOI-1518b fresh start, general EOS + metal ionization, PURE HYDRO, '
                        'am_vector (dhj_planets_1101)', 'planet']
for k, v in kv.items():
    if k.startswith('mesh/') or k.startswith('meshblock/'):
        o[k] = [v, src]
o['time/tlim'] = ['%.6e' % (NROT*P_ROT), '%d rot (README sec. 9); wall clock stops it' % NROT]
o['problem/Teq'] = [teq, 'TOI-1518b: Teff sqrt(R*/2a), Simonnin+2025 Table 7 (7299 K, a/R* 4.109; '
                    'zero albedo, full redistribution)']
o['problem/omega'] = [om, 'TOI-1518b: 2 pi / P_orb, P = 1.90261131 d (Simonnin+2025); tidally '
                      'locked, e < 0.01 (Cabot+2021)']
o['problem/grav'] = [grav, 'planet_setup.py: g at ap, M_p 1.83 M_J (Simonnin+2025)']
o['problem/ap'] = [ap, 'planet_setup.py: r(250 bar), R_p 1.8783 R_J at the TESS p_ref']
o['problem/ic_profile'] = ['ic/ic_toi1518b_3x_f5.txt', 'face-5 exact RCE IC of TOI-1518b']
o['problem/ck_swflux'] = ['sw_band_flux_T1518_11.txt', 'ATLAS9 7299 K, log g 4.136, [M/H] -0.10 '
                          '(sed/make_sed.py), under ck_data_dir/sw_flux']
o['problem/albedo'] = ['0.10', 'TOI-1518b Bond albedo: choice, README sec. 4']
o['output1/dt'] = ['%.2f' % (0.01*P_ROT), '0.01 rot']
o['output2/dt'] = ['%.5e' % (0.1*P_ROT), '0.1 rot']
o['output3/dt'] = ['%.1f' % (0.5*P_ROT), '0.5 rot']
o['output4/dt'] = ['%.1f' % P_ROT, '1 rot']
json.dump(o, open(sys.argv[3], 'w'), indent=1)
print('wrote', sys.argv[3], len(o), 'keys')
