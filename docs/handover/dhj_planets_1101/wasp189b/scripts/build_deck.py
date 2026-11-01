"""build_deck.py (dhj_planets_1101/wasp189b): WASP-189b deck from the W121 hydro fresh deck.
usage: build_deck.py MET(1x|3x) SETUP.txt GRIDKEYS.txt OUT.athinput
Reads ap/grav/Teq/omega from planet_setup's '<problem> keys' line and the radial grid from
grid112.py's KEYS file, writes OUT.json (the overrides) and calls make_deck.py."""
import json
import os
import re
import subprocess
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
met, setup, gk, out = sys.argv[1:5]
P = 2.7240330*86400.0                       # s, Lendl+2020
s = open(setup).read()
m = re.search(r'ap = ([0-9.e+-]+)\s+grav = ([0-9.e+-]+)\s+Teq = ([0-9.]+)\s+omega = ([0-9.e+-]+)', s)
ap, grav, teq, om = m.groups()
keys = dict(kv.split('=') for kv in open(gk).read().split())
T = 'W189'
o = {}
for k, v in keys.items():
    o[k] = [v, 'WASP-189b radial grid (grid112.py: W121 cells/H carried over)']
o['comment/problem'] = ['WASP-189b fresh start, general EOS + metal ionization, PURE HYDRO, am_vector, %s solar' % met,
                        'planet']
o['time/tlim'] = ['%.6e' % (300*P), '300 rot of P = 2.7240330 d (Lendl+2020)']
o['problem/Teq'] = [teq, 'Teff sqrt(R*/2a), Lendl+2020 (= Seidel+2026 2638 +- 28 K)']
o['problem/omega'] = [om, '2 pi / P_orb, P = 2.7240330 d (Lendl+2020); tidally locked']
o['problem/grav'] = [grav, 'GM_p/ap^2, M_p 1.99 M_J (Lendl+2020)']
o['problem/ap'] = [ap, 'inner radius from planet_setup.py (R_p 1.619 R_J at the CHEOPS p_ref)']
o['problem/albedo'] = ['0.10', 'Bond albedo adopted (no measurement; Deline+2022 A_g < 0.48); see README']
o['problem/ck_swflux'] = ['sw_band_flux_W189_11.txt', 'ATLAS9 Teff 8000 logg 4.0 [M/H] +0.29 band fluxes (sed/make_sed.py)']
o['problem/ic_profile'] = ['ic/ic_w189_%s.txt' % met, 'face-5 RCE IC for WASP-189b (rce/)']
if met == '1x':
    o['hydro/rad_met'] = ['0.0', '1x solar']
    o['hydro/eos_xh'] = ['0.7381', '1x solar (as the W121 1x input)']
    o['hydro/eos_yhe'] = ['0.2485', '1x solar (as the W121 1x input)']
    o['problem/met'] = ['0.0', '1x solar']
    o['problem/ck_table'] = ['ck/CK1X/ck/Premixed_1x_g8_11_hiT2.txt', '1x solar tables']
    o['problem/ck_ce_table'] = ['CE_tables/FastChem_ck_1x_int_hiT2.txt', '1x solar tables']
    o['problem/ck_data_dir'] = ['ck/CK1X', '1x solar tables (repo data/exo_fms_ck)']
o['output1/dt'] = ['%.4f' % (0.01*P), '0.01 rot']
o['output2/dt'] = ['%.5e' % (0.1*P), '10 per rotation']
o['output3/dt'] = ['%.2f' % (0.5*P), '0.5 rot']
o['output4/dt'] = ['%.1f' % P, '1 rot']
js = out.replace('.athinput', '.json')
json.dump(o, open(js, 'w'), indent=1)
subprocess.run([sys.executable, HERE + '/make_deck.py', js, out, T,
                'WASP-189b (Seidel+2026 hot end) hydro fresh-start deck, %s solar; package '
                '/viper/ptmp2/jinma/dhj_planets_1101/wasp189b (README.md).' % met,
                'Run dir needs ck -> package ck/ and ic -> package ic/ (table paths are relative).'],
               check=True)
