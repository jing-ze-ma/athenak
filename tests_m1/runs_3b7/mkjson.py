import numpy as np
import json
import struct
import os
import glob
import sys
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert as bc  # noqa: E402
ic = np.loadtxt('/viper/u2/jinma/ATHENAK/bench/m1_stage2/ic/ic_m1_V2.txt')
zic, tauic = ic[:, 0], ic[:, 10]


def hst(d):
    b = 'm1slab' if os.path.exists(d + '/m1slab.hydro.hst') else 'hslab'
    return np.loadtxt(d + '/' + b + '.hydro.hst')


def fit(t, k, t0=1., t1=20.):
    m = (t >= t0) & (t <= t1) & (k > 0)
    return float(np.polyfit(t[m], np.log(k[m]), 1)[0]) if m.sum() > 3 else None


def rtprof(f):
    out = []
    with open(f, 'rb') as fh:
        while True:
            h = fh.read(16)
            if len(h) < 16:
                break
            t, = struct.unpack('<d', h[:8])
            nx1, nv = struct.unpack('<ii', h[8:])
            x = np.frombuffer(fh.read(8 * nx1), '<f8')
            a = np.frombuffer(fh.read(8 * nx1 * nv), '<f8').reshape(nv, nx1)
            out.append((float(t), x.copy(), a[3].copy(), a[4].copy()))
    return out


def binprof(d):
    out = []
    for f in sorted(glob.glob(d + '/bin/*.bin')):
        r = bc.read_binary(f)
        A = np.array(r['mb_data']['velx'])[0, 0]
        B = np.array(r['mb_data']['vely'])[0, 0]
        nx1 = r['Nx1']
        x = r['x1min'] + (np.arange(nx1) + 0.5) * (r['x1max'] - r['x1min']) / nx1
        out.append((float(r['time']), x, (A**2).mean(axis=0), (B**2).mean(axis=0)))
    return out


ARMS = [('A0_base', '../runs_3b6/s_base', 'rt'), ('A1_hydro', 'A1_hydro', 'bin'),
        ('A1_static', 'A1_static', None), ('A2_1d', 'A2_1d', None), ('A3_lo', 'A3_lo',
                                                                     None),
        ('A3_hi', 'A3_hi', None), ('A3_cfl', 'A3_cfl', None), ('A4_nophi', 'A4_nophi',
                                                               'bin'),
        ('A4_norad', 'A4_norad', 'bin'), ('A5_cut', 'A5_cut', 'rt'), ('A6_seedcons',
                                                                      'A6_seedcons', 'rt')]
out = {'note': 'milestone 3b phase F: where the seeded 2-D He slab grows',
       'v_MLT_cgs': 1.86e4,
       'conv_rate_per_s': 2.1e-3, 'arms': {}}
for n, d, mode in ARMS:
    if not os.path.exists(d):
        continue
    a = hst(d)
    t, k1, k2 = a[:, 0], a[:, 7], a[:, 8]
    s = max(1, len(t) // 200)
    e = {'g1': fit(t, k1), 'g2': fit(t, k2), 't': t[::s].round(4).tolist(
    ), 'KE1': k1[::s].tolist(), 'KE2': k2[::s].tolist(), 'KE1_end': k1[-1], 'KE2_end': k2[-1]}
    recs = []
    if mode == 'rt' and os.path.exists(d + '/rt_profile.bin'):
        recs = rtprof(d + '/rt_profile.bin')
    elif mode == 'bin':
        recs = binprof(d)
    if recs:
        want = [0., 25., 100., 200.]
        sel = [min(recs, key=lambda r: abs(r[0] - w)) for w in want]
        e['prof_t'] = [r[0] for r in sel]
        e['x1v'] = sel[0][1].tolist()
        e['tau'] = np.interp(sel[0][1], zic, tauic).tolist()
        e['v1sq'] = [r[2].tolist() for r in sel]
        e['v2sq'] = [r[3].tolist() for r in sel]
    out['arms'][n] = e
out['H1'] = json.load(open('h1.json'))
json.dump(
    out,
    open(
        '/viper/u2/jinma/ATHENAK/athenak/tests_m1/plots/impl3b7_growth.json',
        'w'))
print('ok',
      os.path.getsize('/viper/u2/jinma/ATHENAK/athenak/tests_m1/plots/impl3b7_growth.json'))
