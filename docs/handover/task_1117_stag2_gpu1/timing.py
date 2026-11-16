#!/usr/bin/env python3
"""timing.py [runs]: the T table (staggered vs cell layout of the clean scheme).

Per run dir runs/t_*:
  s/cycle   = (elapsed at the last cycle - elapsed at cycle 5) / (cycles after 5), from the
              driver's 'elapsed=... cycle=N' lines (cycles 0-5 = warm-up incl. initialisation)
  kernels   = kokkos-tools simple kernel timer: kernels.txt (kp_reader output of the run's
              *.dat files; made here if $KP_READER points at kp_reader); top 15 by time
  FV build  = time and passes of the FV ray sweep per ray build: cln_fv time / calls of
              cln_gd_chi (one per ray build), cln_fv calls / cln_gd_chi calls
  memory    = lines of run.log / mem*.txt with 'high water' (space-time-stack) or 'vram'
and staggered / cell per deck.  The cycle-5 warm-up is the same (steady start, ray tables).  Writes nothing; print to timing.txt.
Kernel times are summed over the ranks' .dat files."""
import glob
import os
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ana'))
R = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'runs')


def kernels(d):
    f = d + '/kernels.txt'
    dats = sorted(glob.glob(d + '/*.dat'))
    if not os.path.exists(f) and dats and os.environ.get('KP_READER'):
        with open(f, 'w') as fp:
            subprocess.run([os.environ['KP_READER']] + dats, stdout=fp, stderr=subprocess.STDOUT)
    if not os.path.exists(f):
        return None
    out, name = {}, None
    for ln in open(f, errors='ignore'):
        m = re.match(r'^- (\S.*?)\s*$', ln)
        if m:
            name = m.group(1)
            continue
        m = re.match(r'^\s*\((\w+)\)\s+([-+0-9.eE]+)\s+(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+'
                     r'([-+0-9.eE]+)', ln)
        if m and name:
            t, n = float(m.group(2)), int(m.group(3))
            o = out.setdefault(name, [m.group(1), 0.0, 0, float(m.group(6))])
            o[1] += t
            o[2] += n
            name = None
    return out


t = {}
def accuracy(d, tag):
    try:
        if '_thin' in tag:
            import ct
            return 'L1(E) thin %.4e' % ct.score('thin', d)[0]
        fu = glob.glob(d + '/*.user.hst')
        if fu:
            lab = [ln for ln in open(fu[0]) if ln.startswith('#')][-1]
            names = [re.sub(r'\[\d+\]=', '', x) for x in lab[1:].split()]
            a = np.loadtxt(fu[0], comments='#', ndmin=2)
            return 'max v/c_s %.3e (final %.3e)' % (a[1:, names.index('vcs_max')].max(),
                                                  a[-1, names.index('vcs_max')])
    except Exception as e:  # noqa: BLE001
        return 'accuracy error: %s' % str(e)[:80]
    return None


for d in sorted(glob.glob(R + '/t_*')):
    tag = os.path.basename(d)
    if '.old.' in tag:
        continue
    txt = open(d + '/run.log', errors='ignore').read() if os.path.exists(d + '/run.log') else ''
    e = {int(c): float(s) for s, c in re.findall(r'elapsed=([-+0-9.eE]+) cycle=(\d+)', txt)}
    rc = re.findall(r'^rc=(\d+)', txt, re.M)
    print('\n## %s  rc=%s' % (tag, rc[-1] if rc else '?'))
    if 5 in e and max(e) > 5:
        n = max(e)
        t[tag] = (e[n] - e[5])/(n - 5)
        print('s/cycle %.4e over cycles 6-%d; initialisation + cycles 0-5 %.1f s' % (t[tag], n, e[5]))
    else:
        print('not enough cycles: last elapsed lines %s' % sorted(e.items())[-3:])
    k = kernels(d)
    if k:
        tot = sum(v[1] for v in k.values())
        print('| kernel | type | total s | calls | s/call | % of kernels |')
        print('|---|---|---|---|---|---|')
        for name, (ty, s, c, pk) in sorted(k.items(), key=lambda x: -x[1][1])[:15]:
            print('| %s | %s | %.3f | %d | %.3e | %.1f |' % (name, ty, s, c, s/max(c, 1),
                                                          100.0*s/max(tot, 1e-300)))
        if 'cln_fv' in k and 'cln_gd_chi' in k:
            nb = max(k['cln_gd_chi'][2], 1)
            print('FV ray sweep: %d builds, %.3e s per build, %.1f cln_fv launches per build (passes x wavefront levels)' % (
                nb, k['cln_fv'][1]/nb, k['cln_fv'][2]/nb))
    for f in [d + '/run.log'] + glob.glob(d + '/mem*.txt'):
        if not os.path.exists(f):
            continue
        for ln in open(f, errors='ignore'):
            if re.search(r'high water|vram', ln, re.I):
                print('memory: ' + ln.rstrip())
print('\n## staggered / cell (s/cycle)')
for tag in sorted(t):
    if tag.endswith('_stag'):
        p = tag[:-len('_stag')] + '_cell'
        if p in t:
            print('%-14s staggered %.3e  cell %.3e  ratio %.2f' % (tag[:-5], t[tag], t[p], t[tag]/t[p]))
