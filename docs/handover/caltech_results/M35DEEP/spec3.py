"""spec3.py RUN: per shell from the inner wall to beyond the photosphere: band power and axis power (|kx|<=1 or |ky|<=1) of v_r - <v_r> and of T/<T> - 1. k256 = 256 x cycles/cell. read only."""
import os, re, sys, glob, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d')
import bin_convert as bc
J = '/resnick/groups/carnegie_poc/jingze/'
R = dict(m35zams=('prod_c1_1124/m35zams', 'm35zams', 1.01 * 5.4972e11, 8, 8), t9l62=('prod_c4_1124/t9l62', 't9l62', 3.7e13, 4, 1), t9mlt=('prod_mlt_1125/t9mlt', 't9l62', 5.2e12, 3, 1),
         bsg=('prod_c1_1124/bsg', 'bsg3d', 3.75e12, 4, 1), m35mms=('prod_c1_1124/m35mms', 'm35mms', 1.03 * 9.554e11, 4, 4))
run = sys.argv[1]; D, base, rmax, st, n0 = R[run]; D = J + D + '/'; OUT = J + 'm35deep_1125/stars3/%s/' % run
os.makedirs(OUT, exist_ok=True)
rows = [l.split() for l in open(D + base + '.x1grid.txt') if not l.startswith('#')]; n_ = int(re.search(r'is = 3 ie = (\d+)', open(D + base + '.x1grid.txt').readline()).group(1)) - 2
g = np.array([x[:4] for x in rows if len(x) == 4], float); rc = g[3:3 + n_, 2]; NI = min(len(rc), int(np.searchsorted(rc, rmax)))
BANDS = [(1, 3), (4, 8), (9, 15), (16, 40), (20, 28), (41, 64), (65, 128)]
def asm(d, v, ni):
    lg = np.asarray(d['mb_logical']); n3, n2, n1 = d['mb_data'][v][0].shape
    P = np.empty((n3 * (lg[:, 2].max() + 1), n2 * (lg[:, 1].max() + 1), ni), np.float32)
    for m in range(len(lg)): P[lg[m, 2] * n3:(lg[m, 2] + 1) * n3, lg[m, 1] * n2:(lg[m, 1] + 1) * n2] = np.asarray(d['mb_data'][v][m])[:, :, :ni]
    return P
def one(f):
    n = int(f.split('.')[-2]); o = OUT + 'b_%05d.npz' % n
    if os.path.exists(o): return
    h = bc.read_binary(f); t = float(h['time']); rho = asm(h, 'dens', NI); v = asm(h, 'velx', NI); ei = asm(h, 'eint', NI); del h
    N3, N2 = v.shape[:2]; ky = np.fft.fftfreq(N3, 1. / N3).astype(int); kx = np.fft.fftfreq(N2, 1. / N2).astype(int); KX, KY = np.meshgrid(kx, ky); b = np.rint(256 * np.hypot(KX / N2, KY / N3)).astype(int); ax = (abs(KX) <= 1) | (abs(KY) <= 1)
    out = {}
    nb = len(BANDS); B = [(np.where((b >= lo) & (b <= hi))) for lo, hi in BANDS]
    count = np.array([len(x[0]) for x in B]); cax = np.array([ax[x].sum() for x in B])
    pv = np.zeros((NI, nb)); av = np.zeros((NI, nb)); pT = np.zeros((NI, nb)); aT = np.zeros((NI, nb)); rmsv = np.zeros(NI); rmsT = np.zeros(NI); T0 = np.zeros(NI)
    for i in range(NI):
        x = v[:, :, i].astype(float); x = x - x.mean(); T = ei[:, :, i].astype(float) / rho[:, :, i].astype(float); T0[i] = T.mean(); T = T / T.mean() - 1
        for fld, p, a, rr in ((x, pv, av, 'v'), (T, pT, aT, 'T')):
            F = abs(np.fft.fft2(fld) / (N3 * N2)) ** 2
            for q, (idx) in enumerate(B): p[i, q] = F[idx].sum(); a[i, q] = F[idx][ax[idx]].sum()
        rmsv[i] = x.std(); rmsT[i] = T.std()
    np.savez_compressed(o, t=t, n=n, r=rc[:NI], pv=pv, av=av, pT=pT, aT=aT, rmsv=rmsv, rmsT=rmsT, count=count, cax=cax, N=[N3, N2])
    print('done', run, n, '%.4e' % t, flush=True)
if __name__ == '__main__':
    fs = sorted(glob.glob(D + 'bin/' + base + '.hydro_w.*.bin')); fs = [f for f in fs if time.time() - os.path.getmtime(f) > 180]
    fs = [f for f in fs if (int(f.split('.')[-2]) >= n0 and (int(f.split('.')[-2]) - n0) % st == 0)] + fs[-1:]
    with Pool(int(os.environ.get('NPROC', '4'))) as p: p.map(one, sorted(set(fs)), chunksize=1)
    print('ALL', run, len(fs))
