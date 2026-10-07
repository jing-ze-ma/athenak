import re
s = open('/viper/ptmp2/jinma/accretor_1006/bin/plaskett_env9.athinput').read()
g = open('/viper/ptmp2/jinma/accretor_1006/scripts/grid_plaskett_500.txt').read()
keys = g[g.index('use_grid_stretch_r_poly'):g.index('r_acc (face)')].strip()
racc = re.search(r'r_acc \(face\) = (\S+)', g).group(1)
rep = [
    (r'\nnx1       = 246', '\nnx1       = 500', 2),
    (r'\nx1min     = 6.304002334569', '\nx1min     = 6.3', 1),
    (r'\nuse_grid_stretch_r  = true\nf_stretch_r         = .*',
     '\n# option 2 (10-07): polynomial + plateau radial stretch (scripts/grid_design.py 500\n'
     '# 2.65e-3 170 12, GD_R0 6.3 GD_R1 13.5015 GD_RACC 9 GD_ZOFF -0.03): dr 2.5-2.8e-3 Rsun on\n'
     '# 8.80..9.13 (H_p/4 for c_ph 20 km/s), 8.7e-3 at depth 0.3, dr/r 2.2e-3 outside, ratio\n'
     '# <= 1.081\n' + keys, 1),
    (r'\nr_acc       = .*', '\nr_acc       = ' + racc + '         # Rsun, the face nearest the quoted '
     '9 (index 248)', 1),
    (r'\nenv_cs_ph   = .*', '\nenv_cs_ph   = 20.0                   # km/s REAL surface (Teff ~30 kK, '
     'mu 0.62)', 1),
    (r'\nenv_r_spin  = .*', '\nenv_r_spin  = ' + racc + '         # rotation only below R_acc', 1),
]
for a, b, n in rep:
    s, k = re.subn(a, b, s)
    assert k == n, (a, k)
s = s.replace('hot-surface envelope', 'RESOLVED real-surface envelope')
open('/viper/ptmp2/jinma/accretor_1006/bin/plaskett_env10.athinput', 'w').write(s)
print('ok', racc)
