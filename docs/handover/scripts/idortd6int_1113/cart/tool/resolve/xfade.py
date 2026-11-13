#!/usr/bin/env python3
"""idort-x1-1112 xfade.py <deck>: print 'yes' if vet_gd_rotate_xfade=true is valid on the FINAL deck (vet_gd = true AND
vet_gd_rotate_dt > 0 (number) AND none of the VetGdXfade FATAL combinations: vet_gd_async true, vet_gd_hr_split != 0,
vet_gd_hr_latfit >= 0, vet_rotate_avg_ema_dema true), else 'no <reason>'."""
import sys
sys.path.insert(0, '/viper/ptmp2/jinma/d6int_1113/orion_task/tool')
from blk import read
r = read(sys.argv[1], 'rad_m1')
t = lambda k: r.get(k, 'false').strip().lower() in ('true', '1')
if not t('vet_gd'):
    print('no vet_gd_off'); sys.exit()
dt = r.get('vet_gd_rotate_dt', '')
try:
    on = float(dt) > 0
except ValueError:
    print('no rotate_dt=%s(not_numeric)' % (dt or 'absent')); sys.exit()
if not on:
    print('no rotate_dt=%s' % dt); sys.exit()
if t('vet_gd_async'):
    print('no vet_gd_async'); sys.exit()
if float(r.get('vet_gd_hr_split', '0') or 0) != 0:
    print('no vet_gd_hr_split'); sys.exit()
if float(r.get('vet_gd_hr_latfit', '-1') or -1) >= 0:
    print('no vet_gd_hr_latfit'); sys.exit()
if t('vet_rotate_avg_ema_dema'):
    print('no ema_dema'); sys.exit()
print('yes rotate_dt=%s' % dt)
