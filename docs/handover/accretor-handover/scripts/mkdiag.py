p = '/viper/ptmp2/jinma/accretor_1006/bin/ry_per_accretor_env5.athinput'
s = open(p).read()
s = s.replace('wb_option       = isothermal\n',
              'wb_option       = isothermal\nwb_rmax         = 0.0                # diagnostic: 0 = WB everywhere\n')
s = s.rstrip('\n') + '''

<output3>
file_type  = rst
dt         = 0.0485
'''
open('/viper/ptmp2/jinma/accretor_1006/bin/ry_per_accretor_env5d.athinput', 'w').write(s)
print('ok')
