#!/usr/bin/env python3
"""cmp2.py REFDIR ARMDIR [ARMDIR...] [--rst N]: per arm at restart N (default last common):
deep horizontal v rms by pressure band vs REF, p/rho change, jet band (0.002-0.02 bar) vh rms/max,
plus hst radial KE (col 8) and median dt over the run.  Same proxy vh as oddeven_0930/cmp.py."""
import sys, glob, os
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import dhj_remap as R  # noqa: E402
BANDS = [(100, 1e4), (30, 100), (10, 30), (1, 10), (0.1, 1), (0.02, 0.1), (0.002, 0.02)]
args = [a for a in sys.argv[1:] if not a.startswith('--')]
N = None
if '--rst' in sys.argv:
    N = int(sys.argv[sys.argv.index('--rst') + 1]); args = [a for a in args if a != str(N)]


def prof(fn):
    h = R.read_rst(fn); ng = h['ng']
    s = (slice(None), slice(ng, -ng), slice(ng, -ng), slice(ng, -ng))
    u = h['u']; rho = u[:, 0][s]; p = h['p'][s]
    vh = np.sqrt(u[:, 2][s]**2 + u[:, 3][s]**2)/rho
    return dict(pb=p.mean((0, 1, 2))/1e6, pr=(p/rho).mean((0, 1, 2)),
                vh2=(vh**2).mean((0, 1, 2)), vmax=vh.max((0, 1, 2)), t=h['t'])


def hst(d):
    a = np.loadtxt(os.path.join(d, 'dhj.hydro.hst'))
    return np.median(a[:, 1]), a[-1, 7], a[0, 7]


if N is None:
    N = min(max(int(f[-9:-4]) for f in glob.glob(d + '/rst/dhj.*.rst')) for d in args)
ref = prof('%s/rst/dhj.%05d.rst' % (args[0], N))
print('rst %05d  rot %.3f  ref %s' % (N, ref['t']/1.101535e5, args[0]))
for d in args:
    P = prof('%s/rst/dhj.%05d.rst' % (d, N)); dt, ke1, ke0 = hst(d)
    line = '%-28s dt_med %.3f KE1 %.3e (start %.3e) |' % (d.rstrip('/').split('/')[-1] if 'x1phi' in d else d, dt, ke1, ke0)
    for lo, hi in BANDS:
        k = (ref['pb'] >= lo) & (ref['pb'] < hi)
        v = np.sqrt(P['vh2'][k].mean()); v0 = np.sqrt(ref['vh2'][k].mean())
        dp = np.abs(P['pr'][k]/ref['pr'][k] - 1).max()
        line += ' %g-%g: vh %.1f (%+.1f%%) dpr %.1e' % (lo, hi, v/100, 100*(v/v0 - 1), dp)
        if lo == 0.002:
            line += ' vmax %.0f' % (P['vmax'][k].max()/100)
        line += ' |'
    print(line)
