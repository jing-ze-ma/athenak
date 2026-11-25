"""spec_k.py RUN: lateral spectra of v_r - <v_r> per shell in k per wedge period (20 deg), for m35mms at any lateral resolution; read only."""
import os, re, sys, glob, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d'); import bin_convert as bc
J = '/resnick/groups/carnegie_poc/jingze/'; R = dict(m35c='m35c_1126/M35C', c1='prod_c1_1124/m35mms'); run = sys.argv[1]; D = J + R[run] + '/'; OUT = J + 'm35deep_1125/m35c/%s/' % run; os.makedirs(OUT, exist_ok=True)
rows = [l.split() for l in open(J + 'prod_c1_1124/m35mms/m35mms.x1grid.txt') if not l.startswith('#')]; g = np.array([x[:4] for x in rows if len(x) == 4], float); n_ = int(re.search(r'is = 3 ie = (\d+)', open(J + 'prod_c1_1124/m35mms/m35mms.x1grid.txt').readline()).group(1)) - 2; rc = g[3:3 + n_, 2]; NI = int(np.searchsorted(rc, 0.93 * 9.554e11)); NB = 66
def asm(d, v, ni):
    lg = np.asarray(d['mb_logical']); n3, n2, n1 = d['mb_data'][v][0].shape; P = np.empty((n3 * (lg[:, 2].max() + 1), n2 * (lg[:, 1].max() + 1), ni))
    for m in range(len(lg)): P[lg[m, 2] * n3:(lg[m, 2] + 1) * n3, lg[m, 1] * n2:(lg[m, 1] + 1) * n2] = np.asarray(d['mb_data'][v][m])[:, :, :ni]
    return P
def one(f):
    n = int(f.split('.')[-2]); o = OUT + 's_%05d.npz' % n
    if os.path.exists(o): return
    h = bc.read_binary(f); t = float(h['time']); v = asm(h, 'velx', NI); del h; N3, N2 = v.shape[:2]
    ky = np.fft.fftfreq(N3, 1. / N3).astype(int); kx = np.fft.fftfreq(N2, 1. / N2).astype(int); KX, KY = np.meshgrid(kx, ky); b = np.rint(np.hypot(KX, KY)).astype(int); ax = (abs(KX) <= 1) | (abs(KY) <= 1)
    cnt = np.bincount(b.ravel(), minlength=300)[:NB]; cax = np.bincount(b[ax].ravel(), minlength=300)[:NB]; sp = np.zeros((NB, NI), np.float32); sa = np.zeros((NB, NI), np.float32); rms = np.zeros(NI); mean = np.zeros(NI)
    for i in range(NI):
        x = v[:, :, i]; mean[i] = x.mean(); x = x - mean[i]; rms[i] = x.std(); F = abs(np.fft.fft2(x) / (N3 * N2)) ** 2
        sp[:, i] = np.bincount(b.ravel(), F.ravel(), minlength=300)[:NB]; sa[:, i] = np.bincount(b[ax], F[ax], minlength=300)[:NB]
    np.savez_compressed(o, t=t, n=n, r=rc[:NI], sp=sp, sa=sa, cnt=cnt, cax=cax, rms=rms, mean=mean, N=[N3, N2]); print('done', run, n, '%.4e' % t, flush=True)
if __name__ == '__main__':
    fs = sorted(glob.glob(D + 'bin/m35mms.hydro_w.*.bin'))
    fs = fs if run == 'm35c' else [f for f in fs if int(f.split('.')[-2]) % 4 == 0 and int(f.split('.')[-2]) <= 200]
    with Pool(int(os.environ.get('NPROC', '2'))) as p: p.map(one, fs, chunksize=1)
    print('ALL', run, len(fs))
