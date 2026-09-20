"""Q1c: the emergent top-face flux map (rt_surface.bin) -- is the lr1t 7-8 turnover
excess a few bursting columns or the whole surface?"""
import struct, sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r26')
from common import RUNS, TURN, LW, OMFRAC  # noqa

RTOP = 2.964625e11


def recs(run, want):
    """yield (t, rows) for the records whose time is nearest each wanted turnover."""
    f = open(RUNS[run] + '/rt_surface.bin', 'rb')
    out = []
    while True:
        h = f.read(8)
        if len(h) < 8:
            break
        n = struct.unpack('<q', h)[0]
        t = struct.unpack('<d', f.read(8))[0]
        b = f.read(8*4*n)
        if len(b) < 8*4*n:
            break
        out.append((t, f.tell(), n))
    times = np.array([x[0] for x in out])
    res = []
    for w in want:
        k = int(np.argmin(np.abs(times - w*TURN)))
        t, pos, n = out[k]
        f.seek(pos - 8*4*n)
        a = np.frombuffer(f.read(8*4*n), '<f8').reshape(n, 4)
        res.append((t, a.copy()))
    f.close()
    return res


for run in ('lr1t', 'lr_f100t'):
    print("\n" + "=" * 88)
    print("RUN", run, " emergent top-face flux per angular cell (rt_surface)")
    want = [0.5, 2.0, 4.0, 6.0, 6.8, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
    for t, a in recs(run, want):
        th, ph, F = a[:, 1], a[:, 2], a[:, 3]
        n = len(F)
        # spherical: x2v = theta, x3v = phi; area weight ~ sin(theta)
        wgt = np.sin(th)
        wgt = wgt/wgt.sum()
        Fbar = (F*wgt).sum()
        L = Fbar*4*np.pi*RTOP**2*OMFRAC
        o = np.argsort(F*wgt)[::-1]
        c = np.cumsum((F*wgt)[o])/max((F*wgt).sum(), 1e-300)
        f1 = c[max(int(0.01*n)-1, 0)]
        f5 = c[max(int(0.05*n)-1, 0)]
        f10 = c[max(int(0.10*n)-1, 0)]
        print("  t=%6.3f turn  n=%d  L_out=%7.3f L_w | F: mean %.3e med %.3e max %.3e "
              "min %+.2e | share of L from the brightest 1%%/5%%/10%% of columns: "
              "%.3f %.3f %.3f | frac(F<0) %.3f"
              % (t/TURN, n, L/LW, Fbar, np.median(F), F.max(), F.min(),
                 f1, f5, f10, (F < 0).mean()))
