import re
s = open('/viper/ptmp2/jinma/accretor_1006/bin/ry_per_accretor_env5d.athinput').read()
n0 = len(s)
s, k = re.subn(r'\nenv_cs_ph   = 15.5 .*',
               '\nenv_cs_ph   = 85.3                   # km/s HOT star (H_p = 4 dr at half res; '
               '60.3 at full res)\nenv_cs_stream = 15.5                 # km/s, the stream (window '
               'and dense gas beyond env_r_hot = r_top)', s)
assert k == 1
open('/viper/ptmp2/jinma/accretor_1006/bin/ry_per_accretor_env6.athinput', 'w').write(s)
print('ok', n0, len(s))
