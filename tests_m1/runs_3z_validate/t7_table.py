#!/usr/bin/env python3
"""L1 table of the runs_3z T7 runs: python3 t7_table.py <t7 dir> [idx ...]
For each run dir m{2,5,7}_{h2,be,expl}_n*_c*: L1 rho/T_gas/T_rad vs the reference
(t7_radshock.py; t7x_radshock.py for M0=7 P0=1), shift, spike, precursor, dt,
NON-CONVERGED Picard solves, stage fallbacks and BE steps from run.log."""
import glob
import os
import re
import subprocess
import sys

T7 = '/viper/u2/jinma/ATHENAK/athenak/tests_m1/t7_radshock.py'
T7X = '/viper/ptmp2/jinma/validate_0923/scripts/t7x_radshock.py'
top = sys.argv[1]
idxs = sys.argv[2:] or ['last']


def num(pat, s, default='-'):
    m = re.findall(pat, s)
    return m[-1] if m else default


for d in sorted(glob.glob(os.path.join(top, 'm[257]_*'))):
    name = os.path.basename(d)
    log = open(os.path.join(d, 'run.log')).read() if os.path.exists(
        os.path.join(d, 'run.log')) else ''
    if 'exit=' not in log:
        print('%-22s RUNNING/INCOMPLETE' % name)
        continue
    tabs = sorted(glob.glob(os.path.join(d, 'tab', '*.m1.*.tab')))
    cyc = int(num(r'cycle=(\d+)', log, '0'))
    tend = float(num(r'time=([0-9.eE+-]+) cycle=', log, '0'))
    nc = num(r'NON-CONVERGED=([0-9.eE+-]+)', log)
    fb = num(r'stage fallbacks=([0-9.eE+-]+)', log)
    be = num(r'backward-Euler steps=([0-9.eE+-]+)', log)
    pic = num(r'Picard iterations mean=([0-9.eE+-]+)', log)
    for ix in idxs:
        f = tabs[-1] if ix == 'last' else [t for t in tabs if t.endswith(
            '.%05d.tab' % int(ix))][0]
        h = f.replace('.m1.', '.hydro_w.')
        m0 = name[1]
        if m0 == '7':
            cmd = ['python3', T7X, '--m0', '7', '--p0', '1.0', '--a-rad',
                   '7.5657332502e-11', '--max-shift', '400']
        else:
            cmd = ['python3', T7, '--m0', m0]
        out = subprocess.run(cmd + [f, '--hydro', h], capture_output=True,
                             text=True).stdout
        if ix == idxs[0]:
            o0 = subprocess.run(cmd + [tabs[0], '--hydro', tabs[0].replace(
                '.m1.', '.hydro_w.')], capture_output=True, text=True).stdout
            s0 = re.search(r'shift=([0-9.e+-]+)', o0)
            s0 = float(s0.group(1)) if s0 else 0.0
        l1 = re.search(r'L1 rho=([0-9.e+-]+) T_gas=([0-9.e+-]+) '
                       r'T_rad=([0-9.e+-]+).*shift=([0-9.e+-]+)', out)
        sp = re.search(r'spike T_gas = [0-9.e+]+ \(reference [0-9.e+]+, '
                       r'([+-][0-9.]+) %\)', out)
        pr = re.search(r'precursor length.*reference [0-9.e+-]+, '
                       r'([+-][0-9.]+) %\)', out)
        tt = re.search(r't=([0-9.e+-]+)\s+shock', out)
        if not l1:
            print('%-22s %s analysis failed: %s' % (name, os.path.basename(f),
                                                     out[-300:]))
            continue
        print('%-22s t=%-8s L1 rho %s Tg %s Tr %s dshift %6s | spike %6s%% '
              'prec %6s%% | dt %.3e Picard %s NC %s fallb %s BE %s'
              % (name, tt.group(1) if tt else '?', l1.group(1), l1.group(2),
                 l1.group(3), '%+.2f' % (float(l1.group(4)) - s0), sp.group(1) if sp else '-',
                 pr.group(1) if pr else '-', tend / max(cyc, 1), pic, nc, fb,
                 be))
