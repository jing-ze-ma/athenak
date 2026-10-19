s = open('/viper/ptmp2/jinma/accretor_1006/bin/ry_per_accretor_env5d.athinput').read()
old = 'wb_rmax         = 0.0                # diagnostic: 0 = WB everywhere\n'
assert old in s
s = s.replace(old,
              'wb_rmax         = 4.06               # WB on the envelope only (r <= R_acc)\n'
              'fofc            = true               # first-order flux correction (as BSG)\n'
              'fofc_report     = 282                # fofc firing counts per global x1 cell\n'
              'fofc_report_dt  = 0.002\n')
old = 'nghost    = 2\n'
assert old in s
s = s.replace(old, 'nghost    = 3                        # FOFC with PLM needs 3\n')
open('/viper/ptmp2/jinma/accretor_1006/bin/ry_per_accretor_env5g.athinput', 'w').write(s)
print('ok')
