"""extract.py: m35mms v_r lateral Fourier coefficients per shell (r < 0.90 R_ph,IC) and spectra per dump (hydro_w bins, float32). read only. one output npz per dump."""
import os, re, sys, glob
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d')
import bin_convert as bc
D = '/resnick/groups/carnegie_poc/jingze/prod_c1_1124/m35mms/'; OUT = '/resnick/groups/carnegie_poc/jingze/m35deep_1125/out/'; RP = 9.554e11
rows = [l.split() for l in open(D + 'm35mms.x1grid.txt') if not l.startswith('#')]; n_ = int(re.search(r'is = 3 ie = (\d+)', open(D + 'm35mms.x1grid.txt').readline()).group(1)) - 2
g = np.array([x[:4] for x in rows if len(x) == 4], float); rc = g[3:3 + n_, 2]; NI = int(np.searchsorted(rc, 0.80 * RP)); K = 30
def asm(d, v, ni):
    lg = np.asarray(d['mb_logical']); b = d['mb_data'][v][0]; n3, n2, n1 = b.shape
    P = np.empty((n3 * (lg[:, 2].max() + 1), n2 * (lg[:, 1].max() + 1), ni))
    for m in range(len(lg)): P[lg[m, 2] * n3:(lg[m, 2] + 1) * n3, lg[m, 1] * n2:(lg[m, 1] + 1) * n2] = np.asarray(d['mb_data'][v][m])[:, :, :ni]
    return P
def one(n):
    f = D + 'bin/m35mms.hydro_w.%05d.bin' % n; o = OUT + 'q_%05d.npz' % n
    if os.path.exists(o): return
    h = bc.read_binary(f); t = float(h['time']); v = asm(h, 'velx', NI); rho = asm(h, 'dens', NI); vt = asm(h, 'vely', NI); vp = asm(h, 'velz', NI); ei = asm(h, 'eint', NI); del h
    N3, N2 = v.shape[:2]; th = np.linspace(1.3962634016, 1.745329252, N2 + 1); w = (np.cos(th[:-1]) - np.cos(th[1:])); w = w / w.sum() / N3
    mean = np.einsum('kji,j->i', v, w); rmean = np.einsum('kji,j->i', rho, w); mf = np.einsum('kji,j->i', rho * v, w); rw = np.einsum('kji,j->i', rho * v * v, w)
    dv = v - mean
    rms = np.sqrt(np.einsum('kji,j->i', dv * dv, w))
    vtrms = np.sqrt(np.einsum('kji,j->i', vt * vt + vp * vp, w))
    dr = rho - rmean; drr = np.sqrt(np.einsum('kji,j->i', dr * dr, w)) / rmean
    F = np.fft.fft2(dv, axes=(0, 1)) / (N3 * N2)
    ky = np.fft.fftfreq(N3, 1. / N3).astype(int); kx = np.fft.fftfreq(N2, 1. / N2).astype(int)
    sel = np.ix_(np.where(abs(ky) <= K)[0], np.where(abs(kx) <= K)[0])
    i0 = int(np.searchsorted(rc, 0.705 * RP)); F = F[:, :, i0:]; C = F[sel[0], sel[1]].astype(np.complex64)   # [2K+1, 2K+1, r] (ky: phi index, kx: theta index, in wedge periods)
    n2 = np.sqrt(ky[:, None] ** 2 + kx[None, :] ** 2); b = np.rint(n2).astype(int).ravel(); sp = np.zeros((1, 1))
    dT = (2. / 3.) * ei / rho; np.savez_compressed(o, i0=i0, t=t, n=n, r=rc[:NI], mean=mean, rmean=rmean, mf=mf, rw=rw, rms=rms, vtrms=vtrms, drr=drr, C=C, ky=ky[abs(ky) <= K], kx=kx[abs(kx) <= K], spec=sp.astype(np.float32), Tm=np.einsum('kji,j->i', dT, w) * 0.6215 * 1.66053906660e-24 / 1.380649e-16)
    print('done', n, '%.6e' % t, 'peak rms %.3f km/s at %.3f (r<0.85)' % (rms[rc[:NI] < 0.85 * RP].max() / 1e5, rc[np.argmax(np.where(rc[:NI] < 0.85 * RP, rms, 0))] / RP), flush=True)
if __name__ == '__main__':
    al = sorted(int(f.split('.')[-2]) for f in glob.glob(D + 'bin/m35mms.hydro_w.*.bin'))
    L = [n for n in al if n >= 90]
    L = L[:-1] if len(al) and os.path.getmtime(D + 'bin/m35mms.hydro_w.%05d.bin' % al[-1]) > __import__('time').time() - 120 else L
    with Pool(int(os.environ.get('NPROC', '8'))) as p: p.map(one, L, chunksize=1)
    print('ALL', len(L))
