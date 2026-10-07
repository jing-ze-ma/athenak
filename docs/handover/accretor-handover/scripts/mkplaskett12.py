import re
WB = '8.52'   # R_acc - 1.5 d_pen, d_pen ~0.32 Rsun (ram-pressure balance, PLASKETT.md)
for src, dst in (('plaskett_env10', 'plaskett_env12'), ('plaskett_env11', 'plaskett_env12f')):
    s = open('/viper/ptmp2/jinma/accretor_1006/bin/%s.athinput' % src).read()
    s, k = re.subn(r'\nwb_rmax         = .*',
                   '\nwb_rmax         = ' + WB + '               # WB only below R_acc - 1.5 d_pen', s)
    assert k == 1
    open('/viper/ptmp2/jinma/accretor_1006/bin/%s.athinput' % dst, 'w').write(s)
print('ok')
