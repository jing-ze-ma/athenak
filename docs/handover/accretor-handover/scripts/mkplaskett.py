import re
s = open('/viper/ptmp2/jinma/accretor_1006/bin/ry_per_accretor_env7.athinput').read()
rep = [
    (r'\nnx1       = 564', '\nnx1       = 246', 2),
    (r'\nx1min     = 2.835492977186198', '\nx1min     = 6.304002334569', 1),
    (r'\nx1max     = 16.241', '\nx1max     = 13.5015                  # 0.85 d_L1 (d_L1 15.884)', 1),
    (r'\nf_stretch_r         = .*', '\nf_stretch_r         = -0.761615688648    # -ln(x1max/x1min): '
     'log grid, R_acc = face 115', 1),
    (r'\nm_acc       = .*', '\nm_acc       = 16.0                   # Msun, Wade+2026: gainer at '
     'the start of accretion', 1),
    (r'\nm_don       = .*', '\nm_don       = 18.0                   # Msun (M1,i = 18.2 less winds; '
     'ESTIMATE)', 1),
    (r'\na_sep       = .*', '\na_sep       = 32.5575                # Rsun, Kepler for 34 Msun at '
     '3.69 d', 1),
    (r'\nperiod      = .*', '\nperiod      = 3.69                   # d, Wade+2026 P_i', 1),
    (r'\nr_acc       = .*', '\nr_acc       = 9.0                    # Rsun, Wade+2026 "R~9Rsun"', 1),
    (r'\nt_don       = .*', '\nt_don       = 33000.0                # K, O-type MS donor (ESTIMATE)', 1),
    (r'\nmu_don      = .*', '\nmu_don      = 0.62', 1),
    (r'\nenv_rho_ph  = .*', '\nenv_rho_ph  = 0.3                    # rho_ph/stream peak at '
     'Mdot 1e-4 (ESTIMATE, PLASKETT.md)', 1),
    (r'\nenv_cs_ph   = .*', '\nenv_cs_ph   = 64.8                   # km/s HOT surface: H_p = 4 dr '
     '(real ~20 km/s)', 1),
    (r'\nenv_cs_stream = .*', '\nenv_cs_stream = 21.0                 # km/s, donor Teff 33 kK, '
     'mu 0.62', 1),
    (r'\nenv_r_spin  = .*', '\nenv_r_spin  = 9.0                    # rotation only below R_acc', 1),
]
for a, b, n in rep:
    s, k = re.subn(a, b, s)
    assert k == n, (a, k)
s = s.replace('problem = RY Per accretor, resolved envelope (ry_per_accretor, inner: envelope)',
              'problem = Plaskett progenitor at the start of accretion, hot-surface envelope '
              '(ry_per_accretor, inner: envelope)')
open('/viper/ptmp2/jinma/accretor_1006/bin/plaskett_env9.athinput', 'w').write(s)
print('ok')
