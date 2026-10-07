import re
s = open('/viper/ptmp2/jinma/accretor_1006/bin/plaskett_env10.athinput').read()
rep = [
    (r'\nnghost    = 2', '\nnghost    = 3                        # FOFC with PLM needs 3', 1),
    (r'\nwb_rmax         = 0.0 .*', '\nwb_rmax         = 9.00129260997162   # WB on the star only (R_acc)\n'
     'fofc            = true               # first-order flux correction\n'
     'fofc_report     = 500                # firing counts per global x1 cell\n'
     'fofc_report_dt  = 0.0005', 1),
]
for a, b, n in rep:
    s, k = re.subn(a, b, s)
    assert k == n, (a, k)
open('/viper/ptmp2/jinma/accretor_1006/bin/plaskett_env11.athinput', 'w').write(s)
print('ok')
