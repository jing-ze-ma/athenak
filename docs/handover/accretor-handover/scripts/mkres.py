import re
s = open('/viper/ptmp2/jinma/accretor_1006/bin/ry_per_accretor_env5d.athinput').read()
g = open('/tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/4ac381a7-66c0-44a1-af49-2c70066c365c/'
         'scratchpad/g900.txt').read()
keys = g[g.index('use_grid_stretch_r_poly'):g.index('r_acc (face)')].strip()
racc = re.search(r'r_acc \(face\) = (\S+)', g).group(1)
old = ('use_grid_stretch_r  = true\n'
       'f_stretch_r         = -1.7453230974697702  # -ln(x1max/x1min): exact log grid\n')
assert old in s
s = s.replace(old, '# option 3 (10-07): polynomial + plateau radial stretch, dr 5.9e-4..9.8e-4 Rsun on\n'
              '# R_acc-0.05..+0.08 (3.4-5.6 cells per H_p), log-like dr/r 2.8e-3 outside, neighbour\n'
              '# ratio <= 1.093 (scripts/grid_design.py 900 8.0e-4 300 14)\n' + keys + '\n')
n0 = s.count('nx1       = 564')
assert n0 == 2
s = s.replace('nx1       = 564', 'nx1       = 900')
old = 'r_acc       = 4.06                   # Rsun (a log-grid x1 face)'
assert old in s
s = s.replace(old, 'r_acc       = ' + racc + '        # Rsun, the x1 face nearest 4.06 (face 259)')
open('/viper/ptmp2/jinma/accretor_1006/bin/ry_per_accretor_env8.athinput', 'w').write(s)
print('ok', racc)
