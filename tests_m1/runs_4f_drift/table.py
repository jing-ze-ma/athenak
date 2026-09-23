#!/usr/bin/env python3
"""python3 table.py <run dir> ...: L1 rho/T_gas/T_rad at the last dump against the
reference (t7_radshock.py for M5, t7x_radshock.py M0=7 P0=1 for M7), and the shock
drift in cells = shift(last) - shift(t=0), as in runs_3z_validate/t7_table.py."""
import glob
import os
import re
import subprocess
import sys

T7 = '/viper/ptmp2/jinma/wt_drift/tests_m1/t7_radshock.py'
T7X = '/viper/ptmp2/jinma/wt_drift/tests_m1/runs_3z_validate/t7x_radshock.py'
PAT = r'L1 rho=([0-9.e+-]+) T_gas=([0-9.e+-]+) T_rad=([0-9.e+-]+).*shift=([0-9.e+-]+)'


def run(cmd, f):
    o = subprocess.run(cmd + [f, '--hydro', f.replace('.m1.', '.hydro_w.')],
                       capture_output=True, text=True).stdout
    return re.search(PAT, o)


for d in sys.argv[1:]:
    name = os.path.basename(d.rstrip('/'))
    tabs = sorted(glob.glob(os.path.join(d, 'tab', '*.m1.*.tab')))
    if name[1] == '7':
        cmd = ['python3', T7X, '--m0', '7', '--p0', '1.0', '--a-rad',
               '7.5657332502e-11', '--max-shift', '400']
    else:
        cmd = ['python3', T7, '--m0', name[1]]
    a, b = run(cmd, tabs[0]), run(cmd, tabs[-1])
    if not (a and b):
        print('%-26s analysis failed' % name)
        continue
    print('%-26s L1 rho %s Tg %s Tr %s drift %+.2f cells'
          % (name, b.group(1), b.group(2), b.group(3),
             float(b.group(4)) - float(a.group(4))))
