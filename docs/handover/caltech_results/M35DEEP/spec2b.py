"""spec2.py RUN: lateral spectra of v_r - <v_r> per shell for any wedge run (hydro_w bins). k256 = 256 x cycles per cell (cells per wavelength = 256 / k256). read only."""
import os, re, sys, glob, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d')
import bin_convert as bc
J = '/resnick/groups/carnegie_poc/jingze/'
# name: (dir, basename, rmax [cm], dump stride, first index)
R = dict(m35zams=('prod_c1_1124/m35zams', 'm35zams', 0.985 * 5.4972e11, 8, 8), t9l62=('prod_c4_1124/t9l62', 't9l62', 5.5e12, 5, 1), t9mlt=('prod_mlt_1125/t9mlt', 't9l62', 3.6e12, 5, 1),
         bsg=('prod_c1_1124/bsg', 'bsg3d', 3.4e12, 4, 1), m35mms=('prod_c1_1124/m35mms', 'm35mms', 0.93 * 9.554e11, 12, 10),
         old1121b=('prod_fresh_1121b/m35mms', 'm35mms', 0.93 * 9.554e11, 1, 1), old1121=('prod_fresh_1121/m35mms', 'm35mms', 0.93 * 9.554e11, 1, 1), ia11=('prod_ia11_1123/m35mms', 'm35mms', 0.93 * 9.554e11, 4, 4),
         ia122=('prod_ia122_1123/m35mms', 'm35mms', 0.93 * 9.554e11, 2, 2), rc3=('prod_rc3_1124/m35mms', 'm35mms', 0.93 * 9.554e11, 8, 8), i1=('prod_fresh_1121b/i1_m35mms', 'm35mms', 0.93 * 9.554e11, 1, 0), i1f1k=('prod_fresh_1121b/i1_m35mms_f1k', 'm35mms', 0.93 * 9.554e11, 1, 0), fD0=('m35fast_1123/D0', 'm35mms', 0.93 * 9.554e11, 1, 0), fD1=('m35fast_1123/D1', 'm35mms', 0.93 * 9.554e11, 1, 0), fD2=('m35fast_1123/D2', 'm35mms', 0.93 * 9.554e11, 1, 0), fE0=('m35fast_1123/E0', 'm35mms', 0.93 * 9.554e11, 1, 0), fE4=('m35fast_1123/E4', 'm35mms', 0.93 * 9.554e11, 1, 0), fD1h=('m35fast_1123/D1h', 'm35mms', 0.93 * 9.554e11, 1, 0), rc4v=('prod_rc4v_1124/m35mms', 'm35mms', 0.93 * 9.554e11, 4, 4))
run = sys.argv[1]; D, base, rmax, st, n0 = R[run]; D = J + D + '/'; OUT = J + 'm35deep_1125/stars/%s/' % run
os.makedirs(OUT, exist_ok=True)
XG = D + base + '.x1grid.txt'
if not os.path.exists(XG): XG = J + 'prod_c1_1124/m35mms/m35mms.x1grid.txt'
rows = [l.split() for l in open(XG) if not l.startswith('#')]; n_ = int(re.search(r'is = 3 ie = (\d+)', open(XG).readline()).group(1)) - 2
g = np.array([x[:4] for x in rows if len(x) == 4], float); rc = g[3:3 + n_, 2]; NI = int(np.searchsorted(rc, rmax)); NB = 182
def asm(d, v, ni):
    lg = np.asarray(d['mb_logical']); n3, n2, n1 = d['mb_data'][v][0].shape
    P = np.empty((n3 * (lg[:, 2].max() + 1), n2 * (lg[:, 1].max() + 1), ni))
    for m in range(len(lg)): P[lg[m, 2] * n3:(lg[m, 2] + 1) * n3, lg[m, 1] * n2:(lg[m, 1] + 1) * n2] = np.asarray(d['mb_data'][v][m])[:, :, :ni]
    return P
def one(f):
    n = int(f.split('.')[-2]); o = OUT + 's_%05d.npz' % n
    if os.path.exists(o): return
    h = bc.read_binary(f); t = float(h['time']); ext = [float(h[k]) for k in ('x2min', 'x2max', 'x3min', 'x3max')]; dt_ = str(np.asarray(h['mb_data']['velx'][0]).dtype)
    rho = asm(h, 'dens', NI); v = asm(h, 'velx', NI); ei = asm(h, 'eint', NI); vt = asm(h, 'vely', NI); vp = asm(h, 'velz', NI); del h
    N3, N2 = v.shape[:2]; th = np.linspace(ext[0], ext[1], N2 + 1); w = np.cos(th[:-1]) - np.cos(th[1:]); w = w / w.sum() / N3
    mean = np.einsum('kji,j->i', v, w); rm = np.einsum('kji,j->i', rho, w); em = np.einsum('kji,j->i', ei, w); dv = v - mean; rms = np.sqrt(np.einsum('kji,j->i', dv * dv, w)); vtr = np.sqrt(np.einsum('kji,j->i', vt ** 2 + vp ** 2, w))
    KE = 0.5 * np.einsum('kji,j->i', rho * v * v, w)
    ky = np.fft.fftfreq(N3, 1. / N3).astype(int); kx = np.fft.fftfreq(N2, 1. / N2).astype(int); KX, KY = np.meshgrid(kx, ky)
    b = np.rint(256 * np.hypot(KX / N2, KY / N3)).astype(int); ax = (abs(KX) <= 1) | (abs(KY) <= 1)
    cnt = np.bincount(b.ravel(), minlength=400)[:NB]; cax = np.bincount(b[ax].ravel(), minlength=400)[:NB]
    sp = np.zeros((NB, NI), np.float32); sa = np.zeros((NB, NI), np.float32); top = np.zeros((NI, 6, 4), np.float32)
    for i in range(NI):
        F = abs(np.fft.fft2(dv[:, :, i]) / (N3 * N2)) ** 2
        sp[:, i] = np.bincount(b.ravel(), F.ravel(), minlength=400)[:NB]; sa[:, i] = np.bincount(b[ax], F[ax], minlength=400)[:NB]
        q = np.argpartition(F.ravel(), -6)[-6:]; q = q[np.argsort(F.ravel()[q])[::-1]]
        top[i] = np.stack([KX.ravel()[q], KY.ravel()[q], b.ravel()[q], F.ravel()[q] / F.sum()], 1)
    np.savez_compressed(o, t=t, n=n, r=rc[:NI], mean=mean, rms=rms, vt=vtr, rho=rm, ei=em, KE=KE, sp=sp, sa=sa, cnt=cnt, cax=cax, top=top, ext=ext, N=[N3, N2], dtype=dt_)
    print('done', run, n, '%.4e' % t, flush=True)
if __name__ == '__main__':
    fs = sorted(glob.glob(D + 'bin/' + base + '.hydro_w.*.bin')); fs = [f for f in fs if time.time() - os.path.getmtime(f) > 180]
    fs = [f for f in fs if (int(f.split('.')[-2]) >= n0 and (int(f.split('.')[-2]) - n0) % st == 0)] + fs[-1:]
    with Pool(int(os.environ.get('NPROC', '4'))) as p: p.map(one, sorted(set(fs)), chunksize=1)
    print('ALL', run, len(fs))
