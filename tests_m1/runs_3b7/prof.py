import numpy as np
import struct
import os
import glob
import sys
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert as bc  # noqa: E402
ic = np.loadtxt('/viper/u2/jinma/ATHENAK/bench/m1_stage2/ic/ic_m1_V2.txt')
zic, tauic = ic[:, 0], ic[:, 10]


def from_rtprof(f):
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


def from_bin(d):
    out = []
    for f in sorted(glob.glob(os.path.join(d, 'bin', '*.bin'))):
        r = bc.read_binary(f)
        nx1 = r['Nx1']
        x = r['x1min'] + (np.arange(nx1) + 0.5) * (r['x1max'] - r['x1min']) / nx1
        v1 = np.zeros(nx1)
        v2 = np.zeros(nx1)
        for mb in r['mb_data']['velx']:
            pass
        A = np.array(r['mb_data']['velx'])
        B = np.array(r['mb_data']['vely'])
        # single meshblock spanning x1
        v1 = (A[0, 0]**2).mean(axis=0)
        v2 = (B[0, 0]**2).mean(axis=0)
        out.append((float(r['time']), x, v1, v2))
    return out


def report(name, recs, times):
    print("=== %s" % name)
    for t, x, v1, v2 in recs:
        if times is not None and not any(abs(t - q) < 3.0 for q in times):
            continue
        tau = np.interp(x, zic, tauic)
        v = np.sqrt(v1 + v2)
        i = int(np.argmax(v))
        # fraction of column-integrated v^2 above tau=1
        thin = tau < 1.0
        fr = (v1 + v2)[thin].sum() / max((v1 + v2).sum(), 1e-300)
        print(
            "  t=%6.1f  max vrms=%9.3e at z=%10.3e tau=%8.4g   frac(v^2, tau<1)=%.3f" %
            (t, v[i], x[i], tau[i], fr))
    return recs


R = {}
R['A0_s_base'] = report('A0_s_base(M1 full)', from_rtprof(
    '../runs_3b6/s_base/rt_profile.bin'), None)
R['A5_cut'] = report('A5_cut(M1, top at tau=1)', from_rtprof(
    'A5_cut/rt_profile.bin'), [0, 25, 50, 100, 200])
for n in ['A1_hydro', 'A4_nophi', 'A4_norad']:
    R[n] = report(n, from_bin(n), [0, 25, 100, 200])
np.save('prof_cache.npy', R, allow_pickle=True)
