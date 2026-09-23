"""ck-fast production A/B: every arm's final restart against the T4-every-stage arm.

usage: ab_ana.py <root> [ref=t4]
<root>/<arm>/rst/*.rst (the restart at tlim) and <arm>/run.log (ck_impl_verbose lines).
Prints per arm: T max / rms relative difference to the reference over the active cells
(all, day, night, and three pressure bands), the energy gap ckdesum (mean and max |.|
over calls), passes per call, and the kinked-column counts of scan_dn.py (columns whose
radial osc |T_i - (T_i-1 + T_i+1)/2|/T_i exceeds 0.1 / 0.3, per band, day+night).
"""
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/impl_t4p4_0923')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
from scan4 import decode  # noqa: E402
import bin_convert  # noqa: E402
import dhjcs  # noqa: E402

ROOT = sys.argv[1]
REF = sys.argv[2] if len(sys.argv) > 2 else 't4'
G = np.asarray(bin_convert.read_binary(
    '/viper/u2/jinma/ATHENAK/bench/cs_mhd_prod3/bin/dhj.mhd_w_bcc.00000.bin')['mb_geometry'])
LON = np.empty((24, 16, 16))
for m in range(24):
    xi = np.pi/4*(G[m, 2]+(G[m, 3]-G[m, 2])*(np.arange(16)+0.5)/16)
    eta = np.pi/4*(G[m, 4]+(G[m, 5]-G[m, 4])*(np.arange(16)+0.5)/16)
    XI, ETA = np.meshgrid(xi, eta)
    X = dhjcs.panel_map(m//4, XI, ETA)
    LON[m] = np.arctan2(-X[..., 1], -X[..., 0])
DAY = (np.abs(LON) < np.pi/2)[..., None]
BANDS = [('p>1e-3', 1e30, 1e-3), ('1e-3..1e-5', 1e-3, 1e-5), ('1e-5..1e-7', 1e-5, 1e-7),
         ('<1e-7', 1e-7, 0.0)]    # bar


def final_rst(arm):
    f = sorted(glob.glob(os.path.join(ROOT, arm, 'rst', '*.rst')))
    return f[-1] if f else None


def log_stats(arm):
    txt = open(os.path.join(ROOT, arm, 'run.log')).read()
    gap = [abs(float(g)) for g in re.findall(r'ckdesum=(\S+)', txt)]
    ps = [int(p) for p in re.findall(r'passes=(\d+)', txt)]
    nc = txt.count('NOT-CONVERGED')
    el = re.findall(r'elapsed=(\S+) cycle=(\d+) time=(\S+)', txt)
    return gap, ps, nc, el


def kinks(T, p):
    osc = np.zeros_like(T)
    osc[..., 1:-1] = np.abs(T[..., 1:-1]-0.5*(T[..., :-2]+T[..., 2:]))/T[..., 1:-1]
    out = []
    for nm, a, b in BANDS[1:]:
        sel = (p < a*1e6) & (p >= b*1e6)
        cm = np.where(sel, osc, 0).max(-1)
        out.append('%d/%d' % ((cm > .1).sum(), (cm > .3).sum()))
    return ' '.join(out)


fr = final_rst(REF)
r0 = decode(fr)
T0, p0 = r0['T'], r0['p']
print('reference %s: %s' % (REF, fr))
hdr = ('%-5s %9s %9s %9s %9s %9s' % ('arm', 'max', 'max day', 'rms', 'rms day', 'rms ngt')
       + ''.join(' %17s' % (b[0] + ' d/n') for b in BANDS)
       + ' %9s %9s %6s %4s %s' % ('|gap|mean', '|gap|max', 'pass', 'nc', 'kinks>0.1/>0.3'))
print('T relative difference to %s (active cells; bands: rms day/night by pressure '
      'band [bar])' % REF)
print(hdr)
arms = sorted(os.path.basename(d) for d in glob.glob(os.path.join(ROOT, '*'))
              if os.path.isdir(d))
for a in arms:
    f = final_rst(a)
    if f is None:
        print('%-5s no restart' % a)
        continue
    r = decode(f)
    T, p = r['T'], r['p']
    d = np.abs(T/T0 - 1.0)
    row = '%-5s %9.2e %9.2e %9.2e %9.2e %9.2e' % (
        a, d.max(), d[np.broadcast_to(DAY, d.shape)].max(), np.sqrt((d**2).mean()),
        np.sqrt((d**2)[np.broadcast_to(DAY, d.shape)].mean()),
        np.sqrt((d**2)[~np.broadcast_to(DAY, d.shape)].mean()))
    dayb = np.broadcast_to(DAY, d.shape)
    for nm, hi, lo in BANDS:
        sel = (p0 < hi*1e6) & (p0 >= lo*1e6)
        sd = sel & dayb
        sn = sel & ~dayb
        row += ' %8.1e/%8.1e' % (np.sqrt((d[sd]**2).mean()) if sd.any() else 0.0,
                                 np.sqrt((d[sn]**2).mean()) if sn.any() else 0.0)
    gap, ps, nc, el = log_stats(a)
    row += ' %9.2e %9.2e %6.2f %4d %s' % (
        np.mean(gap) if gap else 0, max(gap) if gap else 0,
        np.mean(ps) if ps else 0, nc, kinks(T, p))
    print(row)
