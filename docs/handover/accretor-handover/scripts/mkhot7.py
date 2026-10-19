s = open('/viper/ptmp2/jinma/accretor_1006/bin/ry_per_accretor_env6.athinput').read()
old = 'env_wall_slip = free '
assert s.count(old) == 1
s = s.replace(old, 'env_r_spin  = 4.06                   # the star rotates (spin term of Phi_wb) only '
                   'below R_acc; atmosphere synchronous\n' + old)
open('/viper/ptmp2/jinma/accretor_1006/bin/ry_per_accretor_env7.athinput', 'w').write(s)
print('ok')
