import numpy as np
import sys
import os
import json
import struct


def hst(d, base):
    f = os.path.join(d, base + '.hydro.hst')
    a = np.loadtxt(f)
    return a[:, 0], a[:, 7], a[:, 8]      # t, KE1, KE2


def fit(t, k, t0=1.0, t1=20.0):
    m = (t >= t0) & (t <= t1) & (k > 0)
    if m.sum() < 4:
        return float('nan')
    return float(np.polyfit(t[m], np.log(k[m]), 1)[0])


def prof(d):
    f = os.path.join(d, 'rt_profile.bin')
    recs = []
    if not os.path.exists(f):
        return recs
    with open(f, 'rb') as fh:
        while True:
            h = fh.read(16)
            if len(h) < 16:
                break
            t, = struct.unpack('<d', h[:8])
            nx1, nv = struct.unpack('<ii', h[8:])
            x = np.frombuffer(fh.read(8 * nx1), '<f8')
            a = np.frombuffer(fh.read(8 * nx1 * nv), '<f8').reshape(nv, nx1)
            recs.append((t, x.copy(), a.copy()))
    return recs


out = {}
for d in sys.argv[1:]:
    name = os.path.basename(d.rstrip('/'))
    base = 'm1slab' if os.path.exists(os.path.join(d, 'm1slab.hydro.hst')) else 'hslab'
    t, k1, k2 = hst(d, base)
    g1, g2 = fit(t, k1), fit(t, k2)
    r = prof(d)
    out[name] = {
        't': t.tolist(), 'KE1': k1.tolist(), 'KE2': k2.tolist(), 'g1': g1, 'g2': g2, \
              'prof_t': [
            float(
                x[0]) for x in r], 'x1v': (
            r[0][1].tolist() if r else []), 'v1sq': [
                    x[2][3].tolist() for x in r], 'v2sq': [
                        x[2][4].tolist() for x in r]}
    print("%-12s g1=%7.4f g2=%7.4f  KE1(end)=%.3e KE2(end)=%.3e  tend=%.1f" %
          (name, g1, g2, k1[-1], k2[-1], t[-1]))
json.dump(out, open(sys.argv[0].replace('an.py', '_tmp.json'), 'w'))
