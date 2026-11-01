"""make_ovr.py (dhj_planets_1101): override JSON for make_deck.py from the chain outputs.
usage: make_ovr.py SETUP.txt GRIDKEYS.txt PROT[s] ALBEDO SWFILE ICFILE NAME REFS OUT.json [MET]
MET = 1x | 3x (default, = the template) | 10x: ck/CE tables, ck_data_dir, EOS composition, met keys
SETUP.txt  = planet_setup.py output (the '<problem> keys:' line: ap grav Teq omega)
GRIDKEYS   = grid112.py gridkeys_<N>.txt (mesh/nx1 ... f_stretch_r_*; mesh/x1max)"""
import json
import re
import sys

st, gk, prot, alb, sw, ic, name, refs, out = sys.argv[1:10]
met = sys.argv[10] if len(sys.argv) > 10 else '3x'
prot = float(prot)
line = [s for s in open(st) if s.startswith('<problem> keys:')][0]
v = dict(re.findall(r'(\w+) = ([0-9.eE+-]+)', line))
kv = dict(t.split('=') for t in open(gk).read().split())
o = {}
for k, x in kv.items():
    b, kk = k.split('/')
    if b == 'mesh' and kk == 'x1min':
        x = v['ap']
    o['%s/%s' % (b, kk)] = [x, 'grid112.py (W121 hfresh cells/H carried to %s)' % name]
o['problem/ap'] = [v['ap'], 'planet_setup.py: r(250 bar) on the polar IC column (%s)' % refs]
o['problem/grav'] = [v['grav'], 'planet_setup.py: GM_p/ap^2 (%s)' % refs]
o['problem/Teq'] = [v['Teq'], 'Teff sqrt(R*/2a), zero albedo (%s)' % refs]
o['problem/omega'] = [v['omega'], '2 pi / P_orb, tidally locked (%s)' % refs]
o['problem/albedo'] = [alb, 'Bond albedo, README section 3']
o['problem/ck_swflux'] = [sw, 'ATLAS9 host SED, sed/make_sed.py (under ck/CK3X/sw_flux)']
o['problem/ic_profile'] = [ic, 'face-5 RCE IC of this planet (package-relative)']
o['time/tlim'] = ['%.6e' % (300*prot), '300 rot']
o['output1/dt'] = ['%.6g' % (0.01*prot), '0.01 rot']
o['output2/dt'] = ['%.6g' % (0.1*prot), '0.1 rot']
o['output3/dt'] = ['%.6g' % (0.5*prot), '0.5 rot']
o['output4/dt'] = ['%.6g' % (1.0*prot), '1 rot']
o['comment/problem'] = ['%s fresh start, general EOS + metal ionization, PURE HYDRO, am_vector '
                        '(dhj_planets_1101/%s)' % (name, name.lower().replace('-', '')), name]
# metallicity set (EOS dumps: 1x deepconv_0925, 3x w121prod_0929/met3, 10x wasp121_0925/met10)
MET = {'1x': ('0.0', '0.7381', '0.2485', 'CK1X'), '3x': ('0.4771', '0.7188', '0.2420', 'CK3X'),
       '10x': ('1.0', '0.6588', '0.2218', 'CK10X')}
if met != '3x':
    mh, xh, yhe, d = MET[met]
    c = '%s solar (README section 2)' % met
    o['problem/met'] = [mh, '[M/H] = log10(%s): EOS metal donors (eos_metal_mh default) ' % met + c]
    o['hydro/rad_met'] = [mh, 'unused under route B, kept consistent; ' + c]
    o['hydro/eos_xh'] = [xh, 'X of the %s EOS dump; ' % met + c]
    o['hydro/eos_yhe'] = [yhe, 'Y of the %s EOS dump; ' % met + c]
    o['problem/ck_table'] = ['ck/%s/ck/Premixed_%s_g8_11_hiT2.txt' % (d, met), c]
    o['problem/ck_ce_table'] = ['CE_tables/FastChem_ck_%s_int_hiT2.txt' % met, c]
    o['problem/ck_data_dir'] = ['ck/%s' % d, c]
json.dump(o, open(out, 'w'), indent=1)
print('wrote', out, len(o), 'keys')
