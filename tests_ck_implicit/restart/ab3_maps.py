"""ab3_maps.py <root> <rot> <out.png> [ref=c2] [arm=e4]: lat-lon maps at <rot> of T and
the zonal wind u at 1e-3 and 1e-6 bar for ref, arm and arm - ref.  T and p come from the
restart (general-EOS caches), u from the hydro_w bin dump of the same output time.
Longitude 0 = substellar; the terminators are at +-90."""
import glob
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import bin_convert  # noqa: E402
import dhjcs  # noqa: E402
from ab3_ana import decode, rsttime, P_ROT  # noqa: E402

ROOT, ROT, OUT = sys.argv[1], float(sys.argv[2]), sys.argv[3]
REF = sys.argv[4] if len(sys.argv) > 4 else 'c2'
ARM = sys.argv[5] if len(sys.argv) > 5 else 'e4'
LEV = [1e-3, 1e-6]   # bar


def load(arm):
    fr = [f for f in glob.glob(os.path.join(ROOT, arm, 'rst', '*.rst'))
          if abs(rsttime(f)/P_ROT - ROT) < 0.01][0]
    fb = None
    for f in glob.glob(os.path.join(ROOT, arm, 'bin', '*.bin')):
        raw = bin_convert.read_binary(f)
        if abs(raw['time']/P_ROT - ROT) < 0.01:
            fb = raw
    r = decode(fr)
    geo = dhjcs.cs_geometry(fb)
    u, v = dhjcs.winds(fb, geo)
    lp = np.log10(r['p'])
    lv = np.log10(np.array(LEV)*1e6)
    T = dhjcs.level_interp(r['T'], lp, lv)
    U = dhjcs.level_interp(u, lp, lv)
    out = {}
    for n, q in (('T', T), ('u', U)):
        la, lo, g = dhjcs.latlon_bin(q, geo['lat'], geo['lon'])
        out[n] = g
    out['lat'], out['lon'] = la, lo
    return out


A, B = load(REF), load(ARM)
fig, ax = plt.subplots(4, 3, figsize=(15, 13), constrained_layout=True)
for il, lev in enumerate(LEV):
    qs = (('T', 'inferno', 'wtemp code units'), ('u', 'RdBu_r', 'code units'))
    for iq, (q, cm, un) in enumerate(qs):
        row = 2*il + iq
        a, b = A[q][il], B[q][il]
        lim = np.nanmax(np.abs(a)) if q == 'u' else None
        for ic, (f, tt) in enumerate(((a, REF), (b, ARM))):
            kw = dict(cmap=cm, shading='auto')
            if q == 'u':
                kw.update(vmin=-lim, vmax=lim)
            else:
                kw.update(vmin=np.nanmin(a), vmax=np.nanmax(a))
            im = ax[row, ic].pcolormesh(A['lon'], A['lat'], f, **kw)
            ax[row, ic].set_title('%s  %s at %.0e bar [%s], rot %.1f'
                                  % (tt, q, lev, un, ROT))
            plt.colorbar(im, ax=ax[row, ic])
        dd = (b/a - 1.0) if q == 'T' else (b - a)
        dl = np.nanmax(np.abs(dd))
        im = ax[row, 2].pcolormesh(A['lon'], A['lat'], dd, cmap='RdBu_r', vmin=-dl,
                                   vmax=dl, shading='auto')
        ax[row, 2].set_title('%s - %s: %s (%s)' % (ARM, REF, q, 'relative' if q == 'T'
                                                   else un))
        plt.colorbar(im, ax=ax[row, 2])
for x in ax.ravel():
    for t in (-90, 90):
        x.axvline(t, color='k', lw=0.6, ls='--')
    x.set_xlabel('longitude from substellar [deg]')
    x.set_ylabel('latitude [deg]')
fig.savefig(OUT, dpi=90)
print('wrote', OUT)
