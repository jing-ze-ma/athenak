#!/usr/bin/env python3
"""mk3d.py OUT_INPUT IC GRID_NPY NX1 [key=value ...]: 3-D wedge input for he_star_m1
(template inputs/radiation/he_presn_m1_wedge.athinput), nx2 = nx3 = 64 in 4 MeshBlocks of
32 x 32, the radial grid (8 poly coefs + optional 2 bumps), he_ic_cols = 5, the FeCZ seed,
mlt_flux_frozen off, outputs: hst 25 s, rst 2350 s, bin (hydro_w, m1) 1175 s.
Extra key=value pairs go to <block>/key (e.g. time/cfl_number=0.9)."""
import re
import sys
import numpy as np
R = 2.3717e11
out, ic, grid, nx1 = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
s = open('/viper/ptmp2/jinma/wt_hepresn/inputs/radiation/'
         'he_presn_m1_wedge.athinput').read()
s = s.replace('R_TOP', '2.425e11').replace('HE_IC_FILE', ic)
s = s.replace('HE_TABLES',
              '/viper/ptmp2/jinma/caltech_handover_0926/athenak_data/he_box')
cs = np.load(grid)
ins = 'use_grid_stretch_r_poly = true\n' + ''.join(
    'f_stretch_r_c%d = %.12e\n' % (k+1, c) for k, c in enumerate(cs[:8]))
for b in range(2):
    if len(cs) > 8 and cs[8+3*b] != 0.0:
        ins += ''.join('f_stretch_r_b%d_%s = %.12e\n' % (b+1, nm, cs[8+3*b+q])
                       for q, nm in enumerate(('amp', 'x', 'w')))
s = s.replace('use_spherical_polar = true\n', 'use_spherical_polar = true\n' + ins, 1)
s = s.replace('nx1       = 640', 'nx1       = %s' % nx1)
s = s.replace('he_ic_cols  = 4 ', 'he_ic_cols  = 5 ')
s = s.replace('he_nfine    = 16384', 'he_nfine    = 262144')
s = s.replace('<rad_m1>\n', '<rad_m1>\nreport_newton_fb = true\n', 1)
seed = ('he_seed     = 3.5e-3\nhe_seed_nk   = 16\nhe_seed_kmin = 1\n'
        'he_seed_kmax = 2\n'
        'he_seed_rad  = true\nhe_seed_rlo  = %.6e\nhe_seed_rhi  = %.6e\n'
        % (0.636*R, 0.933*R))
s = re.sub(r'he_seed     = 0\.0[^\n]*\n', '', s)
for k in ('he_seed_nk', 'he_seed_kmin', 'he_seed_kmax', 'he_seed_rad'):
    s = re.sub(r'\n%s\s*=[^\n]*' % k, '', s)
s = s.replace('<problem>\n', '<problem>\n' + seed, 1)
s = s.replace('<output1>\nfile_type   = hst\ndt          = 5.0',
              '<output1>\nfile_type   = hst\ndt          = 25.0')
s = s.replace('<output2>\nfile_type   = rst\ndt          = 1.0e30',
              '<output2>\nfile_type   = rst\ndt          = 2350.0')
s += ('\n<output3>\nfile_type = bin\nvariable  = hydro_w\ndt = 1175.0\n'
      '\n<output4>\nfile_type = bin\nvariable  = m1\ndt = 1175.0\n')
for kv in sys.argv[5:]:
    k, v = kv.split('=', 1)
    blk, key = k.split('/')
    m = re.search(r'(<%s>.*?\n)(%s\s*=[^\n]*)' % (blk, key), s, re.S)
    if m and '\n<' not in s[m.start(1):m.start(2)][1:]:
        s = s[:m.start(2)] + '%s = %s' % (key, v) + s[m.end(2):]
    else:
        s = s.replace('<%s>\n' % blk, '<%s>\n%s = %s\n' % (blk, key, v), 1)
open(out, 'w').write(s)
