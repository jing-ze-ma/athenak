# add_shell.py X idx [idx ...] : add the native map shell X R_ph(IC) (sh_vr/drho/dT/Tfl_<tag>) to existing red_<idx>.npz
# without a full re-reduce (EOS on that one shell only). USER 10-22: map 0.98 -> 0.96 (0.98 cuts the photospheric T front).
import os, sys, numpy as np
sys.path.insert(0, os.environ['VISPY']); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bin_convert as bc, eos_gen
X = float(sys.argv[1]); TAG = '%03d' % int(round(100 * X)); RPH0 = float(os.environ['RPH0']); TFLOOR = float(os.environ.get('TFLOOR', '2000'))
TBL = eos_gen.Table(); OUT = os.path.join(os.environ['PGOUT'], 'data'); BD = os.environ['BINDIR'].split(':')
def shell(d, var, i):
    nmb = len(d['mb_logical']); per = nmb // 6; n3, n2, n1 = d['mb_data'][var][0].shape; N = n2 * int(round(np.sqrt(per)))
    P = np.empty((6, N, N), np.float64)
    for m in range(nmb):
        p = m // per; lx1, lx2, lx3, lev = d['mb_logical'][m]
        P[p, lx3 * n3:(lx3 + 1) * n3, lx2 * n2:(lx2 + 1) * n2] = d['mb_data'][var][m][..., i]
    return P
for idx in sys.argv[2:]:
    o = os.path.join(OUT, 'red_%s.npz' % idx)
    q = dict(np.load(o))
    if 'sh_vr_' + TAG in q:
        continue
    f = next(os.path.join(b, 'rsg.hydro_w.%s.bin' % idx) for b in BD if os.path.exists(os.path.join(b, 'rsg.hydro_w.%s.bin' % idx)))
    h = bc.read_binary(f); r = np.asarray(h['mb_x1v'][0], float); i = int(np.argmin(abs(r - X * RPH0)))
    rho, vr, ei = shell(h, 'dens', i), shell(h, 'velx', i), shell(h, 'eint', i); del h
    N = rho.shape[1]; x2c = -1 + (np.arange(N) + 0.5) * 2.0 / N
    XI, ETA = np.meshgrid(x2c * np.pi / 4, x2c * np.pi / 4, indexing='xy'); Xt, Yt = np.tan(XI), np.tan(ETA)
    W = np.broadcast_to((1 + Xt**2) * (1 + Yt**2) / np.sqrt(1 + Xt**2 + Yt**2)**3, (6, N, N)); sm = lambda a: (W * a).sum() / W.sum()
    Tg, Tlo = TBL.temperature(rho, ei); Tfl = Tlo | (Tg <= TFLOOR * 1.002)
    q['sh_vr_' + TAG] = (vr / 1e5).astype(np.float32)
    q['sh_drho_' + TAG] = (np.log10(rho) - np.log10(sm(rho))).astype(np.float32)
    q['sh_dT_' + TAG] = (Tg / sm(Tg) - 1).astype(np.float32)
    q['sh_Tfl_' + TAG] = Tfl
    q['map_r'] = np.array(sorted(set(np.atleast_1d(q.get('map_r', [])).tolist()) | {X}))
    np.savez(o + '.tmp.npz', **q); os.replace(o + '.tmp.npz', o)
    print('added', TAG, idx, flush=True)
